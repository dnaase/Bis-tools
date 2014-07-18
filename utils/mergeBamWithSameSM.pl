#!/usr/bin/perl -w
##This script merge the bam files, if no option is assigned, they will be assign their own read group name, 
##if no read group they have, they will have their own read group id based on their file name.
## if option is specified. e.g. --sm IMR90.  then all of bam files' read group id will be their own, but SM tag would be changed to be the same SM tag.
## if some of them do not have read group, then they will be assign read group first, read group id is assigned based on their file name. 

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-1-7

#Usage:  perl mergeBamWithSameSM.pl [option] output.bam bam1 bam2 ... 

use strict;
use Getopt::Long;
use File::Basename;

my $MERGE_BAM = "/home/uec-00/shared/production/software/perl_utils_usc/merge_bams.pl";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";


my $sm="";
my $mdups="";
sub usage {
	print "perl mergeBamWithSameSM.pl [option] output.bam bam1 bam2 ... \n";
	exit(1);
}

GetOptions(
	"sm=s" => \$sm,
	"mdups" => \$mdups,
);
usage() if ( scalar(@ARGV) <= 1 );
my $output=shift(@ARGV);
my @files = @ARGV;

my $out_tmp = $output;
$out_tmp =~ s/\.bam$/.tmp.bam/;

##merge bam files
&merge_bam($output, @files);

#system($cmd);

my $header_sam = "";
if($sm ne ""){
	if($mdups ne ""){
		my $outputdups = $output;
		$outputdups =~ s/bam$/mdups\.bam/;
		runcmd("mv $outputdups $out_tmp");
		$header_sam = &creat_header($out_tmp);
		my $cmd = "$SAMTOOLS reheader $header_sam $out_tmp > $outputdups\n";
		$cmd .= "$SAMTOOLS index $outputdups\n";
		$cmd .= "rm $out_tmp\n";
		$cmd .= "rm $output\n";
		$cmd .= "rm $output".".bai\n";
		runcmd("$cmd");
	}
	else{
		runcmd("mv $output $out_tmp");
	$header_sam = &creat_header($out_tmp);
	my $cmd = "$SAMTOOLS reheader $header_sam $out_tmp > $output\n";
	$cmd .= "$SAMTOOLS index $output\n";
	$cmd .= "rm $out_tmp\n";
	runcmd("$cmd");
	
	}
	runcmd("rm $header_sam");
	
	
}

sub merge_bam{
	my $output = shift @_;
	my @files = @_;
	my $cmd = "VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=false MERGE_SEQUENCE_DICTIONARIES=true CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=2000000 OUTPUT='$output' ";


	for my $file (@files)
	{
		die "$file does not exist\n" if !-e $file;
		if(!&hasReadGroup($file))
		{
			#we are catching the case where a non-rg bam is mixed in with a bunch of bams that do have rg. adding the rg to the non-rg bam before merge.
			print STDERR "some files have readgroups and some dont, so adding readgroup to $file since it didnt\n";
			my $RGfile = basename($file) . ".fixed.bam";
			&addReadGroup($file, $RGfile);
			$cmd .= "INPUT='$RGfile' ";
		}
		else
		{
			$cmd .= "INPUT='$file' ";
		}
	}

	runcmd("$JAVA -Xmx12g -jar $PICARD/MergeSamFiles.jar $cmd");
	my $bai = $output;
	$bai =~ s/bam$/bai/;
	runcmd("mv $bai $output.bai");
	if($mdups ne ""){
		my $outputdups = $output;
		$outputdups =~ s/bam$/mdups\.bam/;
		runcmd("$JAVA -Xms7g -Xmx7g -jar $PICARD/MarkDuplicates.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=dupmets.txt READ_NAME_REGEX=null INPUT=$output OUTPUT=$outputdups");
		my $dupbai = $outputdups;
		$dupbai =~ s/bam$/bai/;
		runcmd("mv $dupbai $outputdups\.bai");
	}
	

	
}


sub creat_header{
	my $input = shift @_;
	my $header = $input;
	$header =~ s/\.bam$/.header.tmp.sam/;
	my $out = $input;
	$out =~ s/\.bam$/.header.sam/;
	my $sub_cmd = "$SAMTOOLS view -H $input > $header\n";
	system($sub_cmd);
	open(FH,"<$header") or die "can't open $header\n";
	open(OUT,">$out") or die "can't open $out\n";
	while(<FH>){
		
		my $line = $_;
		if($line =~ /^\@RG/){
			
			$line =~ s/SM:\S+/SM:$sm/;
		}
		
		print OUT "$line";
	}
	close(FH);
	close(OUT);
	runcmd("rm $header");
	#`rm $header`;
	return $out;
}

sub addReadGroup
{
	my $bamIn = shift @_;
	my $bamOut = shift @_;
	my $date = `date`; chomp $date;
	my $flowcell = "ANALYSIS";
	my $lane = "1";
	my $lib= "UNKNOWN_LIB";
	my $sample = "UNKNOWN_SM";
	if($bamIn =~ /^(.+?)_(.+?)_(\d+)_(.+?)\./)
	{
		$sample = $1;
		$flowcell = $2;
		$lane = $3;
		$lib = $4;
	}
	elsif($output =~ /^(.+?)_(.+?)_(\d+)_(.+?)\./)
	{
		$sample = $1;
		$flowcell = $2;
		$lane = $3;
		$lib = $4;
	}
	else
	{

		$sample = "UNKNOWN_SM";
		
		$flowcell = "ANALYSIS";
		$lane = "1";
		$lib = "UNKNOWN_LIB";
	}
	if($sm ne ""){
		$sample = $sm;
	}
		
	runcmd("$JAVA -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1000000 INPUT='$bamIn' OUTPUT='with_rg_$bamOut' RGID='$flowcell\.$lane' RGLB='$lib' RGPL='illumina Hiseq' RGPU='$flowcell\.$lane' RGSM='$sample' RGCN='USC EPIGENOME CENTER' RGDS='from file $bamIn on $date'");
	#overwrite non-readgroups bams
	runcmd("mv with_rg_$bamOut $bamOut");
	my $bai = "with_rg_$bamOut";
	$bai =~ s/bam$/bai/;
	runcmd("mv $bai $bamOut.bai");
}

sub hasReadGroup
{
	print STDERR "checking for readgroup: ";
	my $bam = shift @_;
	my $header = `$SAMTOOLS view -H $bam`; 
	print "YES, $bam has RG\n" if($header =~ /^\@RG/m);
	print "NO, $bam no RG found\n" if($header !~ /^\@RG/m);
	return 1 if($header =~ /^\@RG/m);
	return 0;
}

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}

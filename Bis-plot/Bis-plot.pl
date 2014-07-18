#!/usr/bin/perl -w
##This script will output epistate plot on each individual region provided for RRBS/WGBS/NOMe-seq.
#It accept bam file and a bed file contains all of interested regions.
#it provide mode for WGBS and NOMe-seq

##example cmd:

##WGBS mode: 		perl Bis-plot.pl [option] --plot_mode 1 region.bed input.bam
##NOMe-seq mode: 	perl Bis-plot.pl [option] --plot_mode 2 region.bed input.bam
##RRBS mode: 		perl Bis-plot.pl [option] --plot_mode 3 region.bed input.bam

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2014-7-16

#Usage: perl Bis-plot.pl [option] region.bed input1.bam input2.bam input3.bam ...

use strict;
use Getopt::Long;
use File::Basename;

#epistate_plot4NOMe.R --f BED_file --l sample_list

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "Currently only support NOMe-seq mode:\n";
    print STDERR "Plot mode 1 (NOMe-seq  mode):	perl Bis-plot.pl [option] --plot_mode 1 region.bed input1.bam input2.bam input3.bam ...\n\n";
    print STDERR "Plot mode 2 (WGBS mode): 		perl Bis-plot.pl [option] --plot_mode 2 region.bed input1.bam input2.bam input3.bam ...\n\n";
    print STDERR "Plot mode 3 (RRBS mode): 		perl Bis-plot.pl [option] --plot_mode 3 region.bed input1.bam input2.bam input3.bam ...\n\n";

    print STDERR " Pre-request:\n"; 
    print STDERR "It requires the installation of R. If you can not start R directly, you could also specify your own version of R with --r FILE\n";
    print STDERR " [Options]:\n\n";
    print STDERR " ###############General options:\n\n";
	print STDERR "  --plot_mode NUM : Specify the plot mode.\n";
	print STDERR "                  Mode 1: --plot_mode 1  output epistate plot on each individual region provided for NOMe-seq\n";
	print STDERR "                  Mode 2: --plot_mode 2  output epistate plot on each individual region provided for WGBS\n";
	print STDERR "                  Mode 3: --plot_mode 3  output epistate plot on each individual region provided for RRBS\n";
	print STDERR "                  (Default: 1.  for NOMe-seq)\n\n";
	print STDERR "  --result_dir DIR : Specify the working directory (Default: use current directory).\n";
	print STDERR "  --bistools_path DIR : Specify the Bis-tools root direcotry (Default: not specified. use environment variable \$BISTOOLS).\n";
	print STDERR "  --r FILE : Specify the full path to use R (Default: not specified. use global setting of R).\n";
	print STDERR "  --mem NUM : Specify the number of gigabytes in the memory to use (Default: 8).\n";
	print STDERR "  --genome FILE: reference genome .fasta file that are used when you map reads in the bam file (Default: not specified. use hg19.fa in Bis-tools' resource directory).\n\n";
	print STDERR "  --dbsnp FILE: dbSNP .vcf file that are used for genotyping (Default: not specified. use dbSNP_135.hg19.sort.vcf in Bis-tools' resource directory).\n\n";
			
    exit(1);
}


##default option setting
my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
my $plot_mode=1;
my $r="R";
my $result_dir=`pwd`;
chomp($result_dir);

my $mem=8;

my $genome="$bistools_path/resource/genome/hg19_rCRSchrm.fa";
my $dbsnp="$bistools_path/resource/dbSNP/dbsnp_135.hg19.sort.vcf";



GetOptions( 
			"result_dir=s" => \$result_dir,
			"bistools_path=s" => \$bistools_path,
			"plot_mode=i" => \$plot_mode,
			"r=s" => \$r,
			"mem=i" => \$mem,
			"genome=s" => \$genome,
			"dbsnp=s" => \$dbsnp,

);

	

##STEP 1: check parameters
check_parameter(@ARGV);


##STEP 2: use Bis-SNP output gch/hcg reads.txt file
my @prefix_lists=wrap_bissnp(@ARGV);

##STEP 3: Input gch/hcg reads.txt file into Huy's R script for the plot
wrap_epiplot($ARGV[0],@prefix_lists);

##########################################All Subroutines###############################################################
sub check_parameter {
	usage() if ( scalar(@ARGV) == 0 );
	if ( scalar(@ARGV) < 2 ) {
    	print STDERR "Need one bed file for the location and at least one bam file for the input NOMe-seq/WGBS/RRBS data\n\n";
    	usage();
	}
}


sub wrap_bissnp{
		my $loc=shift @_;
		my @bams=@_;
		my @prefixs=();
		foreach my $bam(@bams){
			my $dirname = dirname($bam);
			my $prefix=basename($bam);
			$prefix=~s/\.bam$//;
			my $vcf1=$result_dir."/cpg.raw.vcf";
			my $vcf2=$result_dir."/snp.raw.vcf";
			my $hcg=$result_dir."/".$prefix.".hcg.reads.txt";
			my $gch=$result_dir."/".$prefix.".gch.reads.txt";
			my $cmd= "java -jar -Xmx${mem}g $bistools_path/Bis-SNP/Bis-SNP.latest.jar -R $genome -T BisulfiteGenotyper -I $bam -D $dbsnp -vfn1 $vcf1 -vfn2 $vcf2 -L $loc -mmq 30 -mbq 5 -stand_call_conf 20 -stand_emit_conf 0 -cpgreads $hcg -gchreads $gch -out_modes NOMESEQ_MODE -sm GM \n";
			print STDERR "Generating HCG/GCH reads files:\n $cmd\n";
			system($cmd)==0 || die "Unexpected stop when generating HCG/GCH reads files for sample $bam: $! \n";
			push(@prefixs,$prefix);
			`rm $vcf1`;
			`rm $vcf2`;
			`rm $vcf1.MethySummarizeList.txt`;
		}
		return @prefixs;
		
}

sub wrap_epiplot{
	my $loc=shift @_;
	my @prefixs=@_;
	my $sample_list=$result_dir."/".$prefixs[0].".sampleList.".scalar(@prefixs)."samples.txt";
	open(OUT,">$sample_list") or die "can't create file $sample_list:$!\n";
	foreach my $prefix(@prefixs){
		print OUT "$prefix\n";
	}
	close(OUT);

	my $cmd="$r --no-restore --no-save --wd $result_dir --f $loc --l $sample_list < $bistools_path/Bis-plot/epistate_plot4NOMe.R \n";
	print STDERR "Generating NOMe-seq epi-states plot:\n $cmd\n";
	system($cmd)==0 || die "Unexpected stop when generating NOMe-seq epi-states plot: $! \n";
	`rm $sample_list`;
}

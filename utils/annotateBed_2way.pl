#!/usr/bin/perl -w
##This script use BEDtools to annotate bed file with multiple genomic features, with percentage overlap, within same strand, reverse strand or both.
##it extend BEDtools to allow count from 5'end, 3'end, both border or center of bed. it also allows different upstream and downstream. 
##it will generate summary files that count how many regions are overlapped with each feature(p value using hypergenomitric test), 
##how many features are in this region(p value using hypergenomitric test)

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2012-7-25

#Usege:  perl annotateBed.pl [Options] bed/gff_file summary_file.txt --anno feature_bed_file_dir

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print "\nUsage:\n";
    print "perl annotateBed.pl [Options] bed/gff_file summary_file.txt --anno feature_bed_file_dir\n\n";

    print " This script use BEDtools to annotate bed file with multiple genomic features, with percentage overlap, within same strand, reverse strand or both.\n";
	print " it extend BEDtools to allow count from 5'end, 3'end, both border or center of bed. it also allows different upstream and downstream. \n";
	print " it will generate summary files that count how many regions are overlapped with each feature(p value using hypergenomitric test), \n";
	print " how many features are in this region(p value using hypergenomitric test)\n\n";
    

	print "  bed/gff_file     Input bed/gff file which requires the annotation.\n";
	print "  summary_file.txt     Summary statistics about annotation result.\n\n";
	
	print "  [Options]:\n\n";
    print "  --anno FILE : genomic features' bed files, which are used to annotate input bed files, allow multiple files\n\n";
    
     print "  --genome FILE : genome fasta index file, used to generate random number in each of the chromosome\n\n";

    print "  --upstream NUM :  how many bp upstream from  align start position (Default: 0)\n\n";
    print "  --downstream NUM :   how many bp downstream from  align start position (Default: 0)\n\n";
    print "  --overlap_only : only genomic feature overlapped with .bed region will be taken into account\n\n";
    print "  --overlap_percentage NUM :  minimum percentage of .bed regions overlapped with genomic features will be taken into account (Default: 0)\n\n";
    print "  --same_strand : Require same strandedness.  That is, only counts overlaps on the _same_ strand (Default: both strand).\n\n";
	print "  --reverse_strand : Require different strandedness.  That is, only count overlaps on the _opposite_ strand (Default: both strand).\n\n";
	print "  --bedtools_path : path to BEDtools/bin/.if not specified in PATH environment variable\n\n";

	print "  --align_start NUM : (Default: --align_start 1)\n";
	print "                 1) extend upstream, downstream from both of 5', 3' border\n";
	print "                 2) extend upstream, downstream from 5' border\n";
	print "                 3) extend upstream, downstream from 3' border\n";
	print "                 4) extend upstream, downstream from center of .bed region\n";

    exit(1);
}

##default option setting
my $anno_dir = ();
my $genome = "";
my $upstream = 0;
my $downstream = 0;
my $overlap_only = "";
my $overlap_percentage = 0.0;
my $same_strand = "";
my $reverse_strand = "";
my $align_start = 1;

my $bedtools_path = "";


GetOptions( "anno_dir=s" => \$anno_dir,
			"genome=s" => \$genome,
			"upstream=i" => \$upstream,
			"downstream=i" => \$downstream,
			"overlap_only" => \$overlap_only,
			"overlap_percentage=f" => \$overlap_percentage,
			"same_strand" => \$same_strand,
			"reverse_strand" => \$reverse_strand,
			"align_start=i" => \$align_start,
			"bedtools_path=s" => \$bedtools_path);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) != 2 ) {
    print "Wrong number of arguments\n";
    usage();
}

my $input_bed = $ARGV[0];
my $summary_file = $ARGV[1];

my $tmp_file = $input_bed;

##create new file, if not overlap_only mode
if($overlap_only eq ""){
	if($tmp_file =~ /\.bed/){
		$tmp_file =~ s/\.bed//;
		$tmp_file .= ".up$upstream.down$downstream.bed";
	}
	elsif($tmp_file =~ /\.gff/){
		$tmp_file =~ s/\.gff//;
		$tmp_file .= ".up$upstream.down$downstream.gff";
	}
	elsif($tmp_file =~ /\.gtf/){
		$tmp_file =~ s/\.gtf//;
		$tmp_file .= ".up$upstream.down$downstream.gtf";
	}
	elsif($tmp_file =~ /\.vcf/){
		$tmp_file =~ s/\.vcf//;
		$tmp_file .= ".up$upstream.down$downstream.vcf";
	}
	else{
		die "not support this file type yet! Only allow .bed, .gff, .vcf files\n";
	}
	&creat_location_file($input_bed,$tmp_file);
	
}
else{
	
}

my $random_file = "";
if($genome ne ""){
	&creat_abs_random_location_file($tmp_file);
}


opendir(DH,"$anno_dir") or die;
my @anno = readdir(DH);

my @data; # head of each column's name
#my @value_1; #count, percentage of each input bed file's overlap with annotate feature
#my @value_2; #count, percentage of each annotate feature file's overlap with input bed files.

&overlap($tmp_file);

open(OUT,">$summary_file") or die;
print OUT @data;
#print OUT @value_1;
#print OUT @value_2;
close(OUT);
closedir(DH);
if($tmp_file ne $input_bed){
	my $cmd = "rm $tmp_file\n";
	system($cmd);
}
if($random_file ne ""){
	my $cmd = "rm $random_file\n";
	system($cmd);
}


sub overlap{
	my $input_bed = shift(@_);
	my $output_bed = $input_bed.".tmp.annotate.result.txt";
	
	foreach my $ann(@anno){
		next if $ann =~  /^\./;
		my $tmp_ann = $anno_dir."$ann";
		my $header = $bedtools_path;
		$header .= "bedtools intersect -wa -u ";
		if($same_strand ne ""){
			$header .= "-s ";
		}
		elsif($reverse_strand ne ""){
			$header .= "-S ";
		}
		my $cmd = $header."-a $input_bed -b $tmp_ann > $output_bed\n";
		system($cmd);
		$cmd = "wc -l $output_bed\n";
		my $numFeature_1 = `$cmd`;
		$numFeature_1 =~ s/(\d+)\s+\S+/$1/;
		chomp($numFeature_1);
		$cmd = $header."-a $tmp_ann -b $input_bed > $output_bed\n";
		system($cmd);
		$cmd = "wc -l $output_bed\n";
		my $numFeature_2 = `$cmd`;
		$numFeature_2 =~ s/(\d+)\s+\S+/$1/;
		chomp($numFeature_2);
		$cmd = "wc -l $input_bed\n";
		my $numFeature_1_total = `$cmd`;
		$numFeature_1_total =~ s/(\d+)\s+\S+/$1/;
		chomp($numFeature_1_total);
		
		$cmd = "wc -l $tmp_ann\n";
		my $numFeature_2_total = `$cmd`;
		$numFeature_2_total =~ s/(\d+)\s+\S+/$1/;
		chomp($numFeature_2_total);
		
		my $numFeature_1_perc = $numFeature_1/$numFeature_1_total;
		my $numFeature_2_perc = $numFeature_2/$numFeature_2_total;
		$ann =~ s/\.hg19\S+sort\S+bed$//;
		if($random_file eq ""){
			push(@data,"$ann\t$numFeature_1\t$numFeature_1_perc\t$numFeature_2\t$numFeature_2_perc\n");
		}
		else{
			$cmd = $header."-a $random_file -b $tmp_ann > $output_bed\n";
			system($cmd);
			$cmd = "wc -l $output_bed\n";
			my $numFeature_1_rand = `$cmd`;
			$numFeature_1_rand =~ s/(\d+)\s+\S+/$1/;
			chomp($numFeature_1_rand);
			$cmd = $header."-a $tmp_ann -b $random_file > $output_bed\n";
			system($cmd);
			$cmd = "wc -l $output_bed\n";
			my $numFeature_2_rand = `$cmd`;
			$numFeature_2_rand =~ s/(\d+)\s+\S+/$1/;
			chomp($numFeature_2_rand);
			my $numFeature_1_perc_rand = $numFeature_1_rand/$numFeature_1_total;
			my $numFeature_2_perc_rand = $numFeature_2_rand/$numFeature_2_total;
			#my $ratio_1 = $numFeature_1/$numFeature_1_rand;
			#my $ratio_2 = $numFeature_2/$numFeature_2_rand;
			push(@data,"$ann\t$numFeature_1\t$numFeature_1_perc\t$numFeature_2\t$numFeature_2_perc\t$numFeature_1_rand\t$numFeature_1_perc_rand\t$numFeature_2_rand\t$numFeature_2_perc_rand\n");
		}
		
		#push(@value_1,"$numFeature_1\t$numFeature_1_perc\t");
		#push(@value_2,"$numFeature_2\t$numFeature_2_perc\t");
		print "$ann\n";
	}
	#my $cmd = "rm $output_bed\n";
	#`$cmd`;
	
	#push(@head,"\n");
	#push(@value_1,"\n");
	#push(@value_2,"\n");

}


sub creat_location_file{
	my $file = shift(@_);
	
	my $location = shift(@_);

	open(FH,"<$file") or die "can not open file:$!";
	open(OUT,">$location") or die "can not open file:$!";
	while(<FH>){
		chomp;
		my $line = $_;
		my @splitin = split "\t",$line;
		if($align_start == 1){
			my $start = $splitin[1] - $upstream;
			my $end = $splitin[2] + $downstream;
			if($#splitin >= 6){
				if($splitin[5] eq "-"){
					$start = $splitin[1] - $downstream;
					$end = $splitin[2] + $upstream;
				}
			}
			$splitin[1] = $start;
			$splitin[2] = $end;
		}
		elsif($align_start == 2 ){
			my $start = $splitin[1] - $upstream;
			my $end = $splitin[1] + $downstream;
			if($#splitin >= 6){
				if($splitin[5] eq "-"){
					$start = $splitin[2] - $downstream;
					$end = $splitin[2] + $upstream;
				}
			}
			$splitin[1] = $start;
			$splitin[2] = $end;
		}
		elsif($align_start == 3 ){
			my $start = $splitin[2] - $upstream;
			my $end = $splitin[2] + $downstream;
			if($#splitin >= 6){
				if($splitin[5] eq "-"){
					$start = $splitin[1] - $downstream;
					$end = $splitin[1] + $upstream;
				}
			}
			$splitin[1] = $start;
			$splitin[2] = $end;
		}
		elsif($align_start == 4 ){
			my $start = int(($splitin[1] + $splitin[2])/2 -1 - $upstream);
			my $end = int(($splitin[1] + $splitin[2])/2 + $downstream);
			if($#splitin >= 6){
				if($splitin[5] eq "-"){
					$start = int(($splitin[1] + $splitin[2])/2 -1 - $downstream);
					$end = int(($splitin[1] + $splitin[2])/2 + $upstream);
				}
			}
			$splitin[1] = $start;
			$splitin[2] = $end;
		}
		else{
			die "Only accept 1,2,3,4 for --align_start option\n";
		}
		my $newLine = join "\t", @splitin;
		print OUT "$newLine\n";
		
	}
	close(FH);
	close(OUT);
}

sub creat_random_location_file{
	my $file = shift(@_);
	my %genome_hash;
	open(GE,"<$genome") or die "can not open file:$!";
	while(<GE>){
		chomp;
		my @splitin = split "\t";
		$genome_hash{$splitin[0]}=$splitin[1];
	}
	close(GE);
	open(FH,"<$file") or die "can not open file:$!";
	$random_file = $file;
	$random_file =~ s/\.bed//;
	$random_file .= ".random.bed";
	open(OUT,">$random_file") or die "can not open file:$!";
	while(<FH>){
		chomp;
		my @splitin = split "\t";
		my $block_range = $splitin[2] - $splitin[1];
		my $range = $genome_hash{$splitin[0]};
		my $random_number_start = int(rand($range));
		my $end = $random_number_start + $block_range;
		$splitin[1] = $random_number_start;
		$splitin[2] = $end;
		my $newLine = join "\t", @splitin;
		print OUT "$newLine\n";
	}
	close(FH);
	close(OUT);
}

sub creat_abs_random_location_file{
	my $file = shift(@_);
	my %genome_hash;
	my $genome_len=0;
	open(GE,"<$genome") or die "can not open file:$!";
	while(<GE>){
		chomp;
		my @splitin = split "\t";
		$genome_hash{$splitin[0]}=$splitin[1];
		$genome_len+=$splitin[1];
	}
	close(GE);
	open(FH,"<$file") or die "can not open file:$!";
	$random_file = $file;
	$random_file =~ s/\.bed//;
	$random_file .= ".random.bed";
	open(OUT,">$random_file") or die "can not open file:$!";
	while(<FH>){
		chomp;
		my @splitin = split "\t";
		my $block_range = $splitin[2] - $splitin[1];
		my $rand_start = int(rand($genome_len));
		my $rand_chr = "chr1";
		foreach my $key (sort(keys %genome_hash)) {
			my $v = $genome_hash{$key};
			$rand_start = $rand_start-$v;
			if($rand_start < 0 ){
				$rand_start = $rand_start + $v;
				$rand_chr = $key;
				last;
			}
			
		}

		my $end = $rand_start + $block_range;
		if($end > $genome_hash{$rand_chr}){
			$end = $genome_hash{$rand_chr};
		}
		$splitin[0] = $rand_chr;
		$splitin[1] = $rand_start;
		$splitin[2] = $end;
		my $newLine = join "\t", @splitin;
		print OUT "$newLine\n";
	}
	close(FH);
	close(OUT);
}




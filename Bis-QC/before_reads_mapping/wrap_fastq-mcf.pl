#!/usr/bin/perl -w
##This script run fastq-mcf to do adapter sequence trimming (also trimmed bad base quality <= 2 as BSMAP did) on fastq files and then do statistics on adapter contaminated number.
## Allow 10% mismatches of adapter sequences, requires at least 5 bases of adapter sequences, when adapter contamination is more than 0.25%, it will begin to do trimming. 
## skewed option is disabled, quality trimming is <=2
## Total reads number is actually read1+read2 for PE.
## Adapter dimmer here is actually the reads number that is too short after clipping, but most of time it is because of adapter dimmer, not due to quality score trimming.


## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-1-9

#Usage:  perl wrap_fastq-mcf.pl [option] 1st_end.fastq [2nd_end.fastq]

use strict;
use Getopt::Long;
use File::Basename;

my $FASTQMCF = "../../External_tools/ea-utils/fastq-mcf";
my $ADAPTER = "../../resource/illumina_adapter_new.fa";
my $nt=1;
my $adapter_trim=0.01;

GetOptions(
	"fastqmcf=s" => \$FASTQMCF,
	"adapter=s" => \$ADAPTER,
	"nt=i" => \$nt,
	"adapter_trim=f" => \$adapter_trim,
);





my $read1 = $ARGV[0] || die "no input reads specified";
my $read2 = $ARGV[1];
my $out1 = $read1;
$out1=~s/\.\w+$/.fastq-mcf.fastq/;
my $out2 = "";
if($read2 ne ""){
	$out2 = $read2;
	$out2=~s/\.\w+$/.fastq-mcf.fastq/;
}
my $error_log = $read1.".fastq-mcf.log";

my $cmd = "$FASTQMCF $ADAPTER -p $nt -t $adapter_trim -k 0 -q 3 $read1 -o $out1 ";
if($read2 ne ""){
	$cmd .= "$read2 -o $out2";	
}
$cmd .= "> $error_log \n";
system($cmd);

&extract_result($error_log);

sub extract_result{
	my $file = shift @_;
	open(FH,"<$file") || die "Fastq-mcf log file $file is missing";
	my $qual_trimmed=0;
	my $adapter_dimmer=0;
	my $adapter_trimmed=0;
	my $total_reads=0;
	while(<FH>){
		chomp;
		my $line=$_;
		if($line=~/^Total reads: (\d+)/){
			$total_reads = $1;
		}
		elsif($line=~/^Too short after clip: (\d+)/){
			$adapter_dimmer += $1;
		}
		elsif($line=~/^Clipped.*Count (\d+)/){
			$adapter_trimmed += $1;
		}
		elsif($line=~/^Trimmed (\d+) reads/){
			$qual_trimmed += $1;
		}
	}
	close(FH);
	if($read2 ne ""){
		$total_reads *= 2;
		$adapter_dimmer *= 2;
	}
	my $qual_trim_perc=sprintf("%.2f",100*$qual_trimmed/$total_reads);
	my $adapter_dimmer_perc=sprintf("%.2f",100*$adapter_dimmer/$total_reads);
	my $adapter_trimmed_perc=sprintf("%.2f",100*$adapter_trimmed/$total_reads);
	my $stat_log= $read1.".fastq-mcf.log.txt";
	open(OUT,">$stat_log") or die "can't create adapter contamination log file: $!\n";
	print OUT "Total reads: $total_reads\n";
	print OUT "Adapter dimmer: $adapter_dimmer (${adapter_dimmer_perc}%)\n";
	print OUT "Adapter trimmed sequences: $adapter_trimmed (${adapter_trimmed_perc}%)\n";
	print OUT "Quality score trimmed sequences(base qual <= 2): $qual_trimmed (${qual_trim_perc}%)\n";
	close(OUT);
	`rm $error_log`;
}

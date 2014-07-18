#!/usr/bin/perl -w
## Methylation bias plot
##This script run Bis-SNP in a quick way to get the WCW (default, in order also to adapt to NOMe-seq) methylation level along the sequencing cycle
##bam file need to be indexed, sorted and have read group.

## author: Yaping Liu  lyping1986@gmail.com 

#Usage:  perl methylation_bias_plot.pl [option] input.bam

use strict;
use Getopt::Long;
use File::Basename;

my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
my $BISSNP = "$bistools_path/Bis-SNP/Bis-SNP.latest.jar";
my $R = "R";
my $pattern = "WCW";
my $genome="$bistools_path/resource/genome/hg19_rCRSchrm.fa";
my $mem="8"; #how many Giga bytes memory need
my $r_script = "$bistools_path/Bis-QC/after_reads_mapping/methyBiasDistPlot.R";

GetOptions(
	"bissnp=s" => \$BISSNP,
	"r=s" => \$R,
	"pattern=s" => \$pattern,
	"genome=s" => \$genome,
	"mem=i" => \$mem,
);

my $file=$ARGV[0];

##generate pattern methylation matrix file:
my $out=$file;
$out =~ s/\.bam$/.${pattern}.methy.cycle.txt/;
my $cmd="java -Xmx${mem}g -jar $BISSNP -T QuickMethylationLevel -R $genome -I $file -pattern $pattern -patternHist $out\n";
print STDERR $cmd;
system($cmd)==0 || die "can't generate methylation matrix file in methylation bias check part:$!\n";

##generate methylation bias plot:
my $pdf=$file;
$pdf=~s/\.bam$/.${pattern}.methy_bias_plot.pdf/;
my $r_cmd="$R --no-restore --no-save --args input=$out output=$pdf < $r_script\n";
print STDERR $r_cmd;
system($r_cmd)==0 || die "can't generate methylation bias plot in methylation bias check part:$!\n";

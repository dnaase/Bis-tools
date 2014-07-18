#!/usr/bin/perl -w
## Bisulfite conversion distribution within the read plot
##This script run Bis-SNP in a quick way to get the WCW (default, in order also to adapt to NOMe-seq) methylation level distribution
##bam file need to be indexed, sorted and have read group.

## author: Yaping Liu  lyping1986@gmail.com 

#Usage:  perl bisulfite_conv_distribution.pl [option] input.bam

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
my $r_script = "$bistools_path/Bis-QC/after_reads_mapping/bisulfiteConvDistPlot.R";

GetOptions(
	"bissnp=s" => \$BISSNP,
	"r=s" => \$R,
	"pattern=s" => \$pattern,
	"genome=s" => \$genome,
	"mem=i" => \$mem,
);

my $file=$ARGV[0];

##generate pattern methylation histgram file:

my $out_hist=$file;
$out_hist=~ s/\.bam$/.hist.txt/;

my $cmd .= "java -Xmx${mem}g -jar $BISSNP -T BisulfiteConversionCheck -R $genome -I $file -pattern $pattern -patternHist $out_hist \n";
print STDERR $cmd;
system($cmd)==0 || die "can't generate methylation histogram file in bisulfite conversion distribution check part:$!\n";

##generate methylation bias plot:
my $pdf=$file;
$pdf=~s/\.bam$/.${pattern}.bisuflite_conv_distribution_plot.pdf/;
my $r_cmd="$R --no-restore --no-save --args input=$out_hist output=$pdf < $r_script\n";
print STDERR $r_cmd;
system($r_cmd)==0 || die "can't generate methylation bias plot in bisulfite conversion distribution check part:$!\n";



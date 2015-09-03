#!/usr/bin/perl -w
## align big wig file in to defined genomic windows
## author: Yaping Liu  lyping1986@gmail.com
##example script:

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl alignWig2Bed.bigWigAverageOverBed.pl [Options] prefix location.bed input.bw\n\n";

    print STDERR " align big wig file in to genome window\n";
   

	
	print STDERR "  [Options]:\n\n";
 	print STDERR " #####General option:\n\n";
 	
  	print STDERR "  --min_data NUM: the miimum data point required in the region.(Default: 1)\n\n";
  	print STDERR "  --useMean0: average over bases with non-covered bases counting as zeroes.(Default: not enabled, average over just covered bases)\n\n";
  	print STDERR "  --NAas STR: what is NA be treated in output.(Default: NA)\n\n";
  	print STDERR "  --format_value:  just keep 3 digit after decimal.(Default: not enabled, keep all digit in the output)\n\n";
 	 exit(1);
}

##default option setting

my $min_data=1;
my $useMean0="";
my $NAas="NA";
my $format_value="";
my $ucsc_script = "bigWigAverageOverBed";

print STDERR "perl alignWig2Bed.bigWigAverageOverBed.pl ";
my $cmd_root=join " ", @ARGV;
print STDERR "$cmd_root\n\n";

GetOptions( 
			"min_data=i" => \$min_data,
			"useMean0" => \$useMean0,
			"NAas=s" => \$NAas,
			"format_value" => \$format_value,
);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) < 3 ) {
    print STDERR "Wrong number of arguments\n";
    usage();
}
my $prefix=$ARGV[0];
my $loc_bed=$ARGV[1];
my $input_wig=$ARGV[2];

my $bw_prefix=basename($input_wig);
$bw_prefix=~s/(\w+)\S+/$1/;
my $out=$prefix.".alignTo.$bw_prefix.min_data-$min_data.txt";


my $tmp_loc=$loc_bed.".tmp_loc.".rand().".bed";
open(G,"<$loc_bed") or die;
open(TMP,">$tmp_loc") or die "can't create file $tmp_loc: $!\n";
while(<G>){
	chomp;
	my @f=split "\t";
	print TMP "$f[0]\t$f[1]\t$f[2]\t$f[0]:$f[1]-$f[2]\n";
}
close(G);
close(TMP);

my $tmp_out=$prefix.".alignTo.$bw_prefix.min_data-$min_data.tmp_out.".rand().".tab";
my $tmp_paste=$prefix.".alignTo.$bw_prefix.min_data-$min_data.tmp_paste.".rand().".tab";

my $cmd_tmp="/home/unix/lyping/compbio/software/UCSC_tools/bigWigAverageOverBed $input_wig $tmp_loc $tmp_out\n";
system($cmd_tmp) == 0 or die "can't execute $cmd_tmp :$!\n";

my $wc_tmp=`wc $loc_bed`;
chomp($wc_tmp);
my @wc_tmps=split " ",$wc_tmp;
my $colNum=0;
if($wc_tmps[0]>0){
	$colNum=int($wc_tmps[1]/$wc_tmps[0]);
}else{
	die "no input line in $loc_bed\n";
}

my $tab_start=$colNum+3;
my $tab_end=$colNum+6;
`paste $loc_bed $tmp_out | cut -f1-$colNum,$tab_start-$tab_end > $tmp_paste`;


open(F,"<$tmp_paste") or die;
open(OUT,">$out") or die;
while(<F>){
	chomp;
	my @f=split "\t";
	my $mean=$NAas;
	if($f[$colNum] >= $min_data){
		if($useMean0 ne ""){
			if($format_value eq ""){
				$mean=$f[$colNum+2];
			}else{
				$mean=sprintf("%.3f", $f[$colNum+2]);
			}
			
		}else{
			if($format_value eq ""){
				$mean=$f[$colNum+3];
			}else{
				$mean=sprintf("%.3f", $f[$colNum+3]);
			}
			
		}
		
	
	}
	for(my $i=0;$i<$colNum;$i++){
		print OUT "$f[$i]\t";
	}	
	print OUT "$mean\n";
	 
}
close(F);
close(OUT);

#`unlink $tmp_out`;
#`unlink $tmp_paste`;
#`unlink $tmp_loc`;
#-type=coverage 2>/dev/null	

			
			

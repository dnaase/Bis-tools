#!/usr/bin/perl -w
## align big wig file in to defined genomic windows
## author: Yaping Liu  lyping1986@gmail.com
##example script:

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl alignWig2BedWindow.pl [Options] prefix genome.chrom.size input.bw\n\n";

    print STDERR " align big wig file in to genome window\n";
   

	
	print STDERR "  [Options]:\n\n";
 	print STDERR " #####General option:\n\n";
 	
 	print STDERR "  --window NUM: window size to calculate mean signal intensity.(Default: 100000)\n\n";
  	print STDERR "  --min_data NUM: the miimum data point required in the region.(Default: 10)\n\n";
  	print STDERR "  --useMean0 : average over bases with non-covered bases counting as zeroes.(Default: not enable, average over just covered bases)\n\n";	
 	 exit(1);
}

##default option setting
my $window=100000;
my $useMean0="";
my $min_data=10;
my $ucsc_script = "bigWigAverageOverBed";

print STDERR "perl alignWig2BedWindow.pl ";
my $cmd_root=join " ", @ARGV;
print STDERR "$cmd_root\n\n";

GetOptions( 
			"window=i" => \$window,
			"min_data=i" => \$min_data,
			"useMean0" => \$useMean0,
);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) < 3 ) {
    print STDERR "Wrong number of arguments\n";
    usage();
}
my $prefix=$ARGV[0];
my $genome_size=$ARGV[1];
my $input_wig=$ARGV[2];

my $out=$prefix.".window-$window.min_data-$min_data.bed";

my $tmp_loc=$prefix.".window-$window.min_data-$min_data.tmp_loc.bed";
my $tmp_out=$prefix.".window-$window.min_data-$min_data.tmp_out.tab";
my $tmp_paste=$prefix.".window-$window.min_data-$min_data.tmp_paste.tab";
open(G,"<$genome_size") or die;
open(TMP,">$tmp_loc") or die "can't create file $tmp_loc: $!\n";
while(<G>){
	chomp;
	my ($chr,$chr_end)=split "\t";
	for(my $start=0;$start<=$chr_end-$window;$start+=$window){
		my $end=$start+$window;
		print TMP "$chr\t$start\t$end\t$chr:$start-$end\n";
	}
}
close(G);
close(TMP);


`bigWigAverageOverBed $input_wig $tmp_loc $tmp_out`;
`paste $tmp_loc $tmp_out | cut -f1-3,7-10 > $tmp_paste`;



open(F,"<$tmp_paste") or die;
open(OUT,">$out") or die;
while(<F>){
	chomp;
	my @f=split "\t";
	my $mean="NA";
		if($f[3] >= $min_data){
			if($useMean0 ne ""){
				$mean=sprintf("%.3f", $f[5]);
			}else{
				$mean=sprintf("%.3f", $f[6]);
			}
			
		}
		
	
	print OUT "$f[0]\t$f[1]\t$f[2]\t$mean\n";
	 
}
close(F);
close(OUT);

`unlink $tmp_loc`;
`unlink $tmp_out`;
`unlink $tmp_paste`;

#-type=coverage 2>/dev/null	

			
			

#!/usr/bin/perl -w
## align big wig file in to defined genomic windows
## author: Yaping Liu  lyping1986@gmail.com
##example script:

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl alignWigToWindow.pl [Options] prefix genome.chrom.size input.bw\n\n";

    print STDERR " align big wig file in to genome window\n";
   

	
	print STDERR "  [Options]:\n\n";
 	print STDERR " #####General option:\n\n";
 	
 	print STDERR "  --window NUM: window size to calculate mean signal intensity.(Default: 100000)\n\n";
 	print STDERR "  --not_include FILE: the region not included in the calculation.(Default: NA)\n\n";
 	
 	 exit(1);
}

##default option setting
my $window=100000;
my $not_include="";
my $min_data=10;
my $ucsc_script = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/bigWigSummary";

print STDERR "perl alignWigToWindow.pl ";
my $cmd_root=join " ", @ARGV;
print STDERR "$cmd_root\n\n";

GetOptions( 
			"window=i" => \$window,
			"not_include=s" => \$not_include,
			"min_data=i" => \$min_data,
);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) < 2 ) {
    print STDERR "Wrong number of arguments\n";
    usage();
}
my $prefix=$ARGV[0];
my $genome_size=$ARGV[1];
my $input_wig=$ARGV[2];

my $out=$prefix.".window-$window.min_data-$min_data.bed";
open(G,"<$genome_size") or die;
open(OUT,">$out") or die "can't create file $out: $!\n";
while(<G>){
	chomp;
	my ($chr,$chr_end)=split "\t";
	for(my $start=0;$start<=$chr_end-$window;$start+=$window){
		my $end=$start+$window;
		my @locations=();
		if($not_include ne ""){
			@locations=`echo "$chr\t$start\t$end\n" | subtractBed -a stdin -b $not_include`;
		}else{
			push(@locations,"$chr\t$start\t$end\n");
		}
		my $sum=0;
		my $cov=0;
		foreach my $loc(@locations){
			chomp($loc);
			my ($c,$s,$e)=split "\t",$loc;
			my $mean=`$ucsc_script $input_wig $c $s $e 1 2>/dev/null `;
			my $cov_perc=`$ucsc_script $input_wig $c $s $e 1 -type=coverage 2>/dev/null `;
			if($mean ne ""){
				$sum += $mean*($e-$s)*$cov_perc;
				$cov += ($e-$s)*$cov_perc;
			}
			
		}
		if($cov >= $min_data){
			my $total_mean=$sum/$cov;
			print OUT "$chr\t$start\t$end\t$total_mean\n";
		}else{
			print OUT "$chr\t$start\t$end\tNA\n";
		}
		

		
	}
}
close(G);
close(OUT);

#-type=coverage 2>/dev/null	

			
			
			
#!/usr/bin/perl -w
## output only unique line when check the columns user defined
## author: Yaping Liu  lyping1986@gmail.com


use strict;
use Getopt::Long;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl uniqLine.pl [Options] input.bed output.bed \n\n";

    print STDERR " output only unique line when check the columns user defined\n\n";
   	
	print STDERR "  [Options]:\n\n";
	print STDERR "  --c NUM: give the column number that are used to check duplication. could be multiple. (e.g. Default: --c 1 --c 2 --c 3 for bed files and --c 1 --c 4 --c 5 for gtf files)\n\n";
	exit(1);
}

##default option setting

my @columns = ();

GetOptions( 
			"c=i" => \@columns,
			);
usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) < 2 ) {
    print STDERR "Wrong number of arguments\n";
    usage();
}

if(scalar(@columns)==0){
	@columns = (1,2,3);
}

my $input = shift(@ARGV);
my $output = shift(@ARGV);
my %hash=();
open(FH,"<$input") or die "can not open input file: $input $!\n";
open(OUT,">$output") or die "can not open output file: $output $!\n";
my $total=0;
my $unique=0;
while(<FH>){
	chomp;
	next if($_=~/^#/);
	$total++;
	my @f=split "\t";
	
	my $key="";
	foreach my $col(@columns){
		my $col_tmp=$col-1;
		#print "$col\t$f[$col_tmp]\n";
		$key=$key.":".$f[$col_tmp];
	}
#	if($_=~/A0A542/){
#		my $s=(! exists $hash{$key});
#		print "$_\n$key\t$s\n";
#	}
	#print "$key\n";
	if(! exists $hash{$key}){
		my $line=join "\t",@f;
		print OUT "$line\n";
		$hash{$key}=1;
		$unique++;
	}
}
close(FH);
close(OUT);
print STDERR "perl uniqLine.pl ";
foreach my $col(@columns){
	print STDERR "--c $col ";
}
print STDERR "$input $output \n";

print STDERR " $total lines in total input\n";
print STDERR " $unique unique line in total output\n";

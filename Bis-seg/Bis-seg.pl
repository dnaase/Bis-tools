#!/usr/bin/perl -w
##This script will segment accessiblity signal provided from NOMe-seq.
#It accept bed file (.6plus2.bed provided by bis-snp or standard bed format) contains methylation level and coverage information.
#it will output NDR, Linker, short protected regions and Mono-Nucleosome regions.


## author: Yaping Liu  lyping1986@gmail.com 
## time: 2014-7-22

#Usage: perl Bis-seg.pl [option] prefix input.bed

use strict;
use Getopt::Long;
use File::Basename;


sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl Bis-seg.pl [option] prefix input.bed\n\n";

    print STDERR " [Options]:\n\n";
    print STDERR " ###############General options:\n\n";
	print STDERR "  --mode NUM : 1. Training + Decoding; 2. Training. 3:Decoding (Default: 1, training+decoding).\n";
	print STDERR "  --result_dir DIR : Specify the working directory (Default: use current directory).\n";
	print STDERR "  --bistools_path DIR : Specify the Bis-tools root direcotry (Default: not specified. use environment variable \$BISTOOLS).\n";
	print STDERR "  --mem NUM : Specify the number of gigabytes in the memory to use (Default: 15).\n";
	print STDERR "  --bed_format NUM : Specify the input bed format. 1: 6plus2 bed format. 2: standard bed format. (Default: 1).\n";
	print STDERR "  --L FILE/STR: Specify the region for the training/decoding (Default: -L chr21 ).\n";
	print STDERR "  --sig_test NUM: Specify the significant test used for p value calculation. 1: One way Binomial test. 2:Fisher exact test (Default: -sig_test 1).\n";
	print STDERR "  --fdr NUM: False Discovery Rate (FDRs) criteria to call significant segment (Default: 0.01).\n\n";
			
    exit(1);
}


##default option setting
my $mode=1;
my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
my $result_dir=`pwd`;
chomp($result_dir);

my $mem=15;
my $region="chr21";
my $bed_format=1;
my $sig_test=1;
my $fdr=0.01;

GetOptions( 
			"mode=i" => \$mode,
			"result_dir=s" => \$result_dir,
			"bistools_path=s" => \$bistools_path,
			"mem=i" => \$mem,
			"L=s" => \$region,
			"bed_format=i" => \$bed_format,
			"sig_test=i" => \$sig_test,
			"fdr=f" => \$fdr,
);

##STEP 1: check parameters
check_parameter(@ARGV);

if($mode == 1 || $mode == 2){
	training_hmm(@ARGV);
}

if($mode == 1 || $mode == 3){
	decoding_hmm(@ARGV);
}


##########################################All Subroutines###############################################################
sub check_parameter{
	usage() if ( scalar(@ARGV) == 0 );
	if ( scalar(@ARGV) < 2 ) {
    	print STDERR "Need prefix for the new file's name and one input bed file\n\n";
    	usage();
	}
	
	if($mode != 1 && $mode != 2 && $mode != 3){
		print STDERR "Wrong mode number!!\n\n";
   		 usage();
	}
}

sub training_hmm{
	my $prefix=shift @_;
	my $input=shift @_;
	my $hmm= $result_dir."/${prefix}.trainedHMM.model.txt";
	$prefix = $result_dir."/$prefix";	
	my $cmd="java -Xmx${mem}g -jar $bistools_path/Bis-seg/NdrHmmHunter.jar $prefix $input -hmmFile $hmm -onlyTrain -sigTestMode $sig_test -L $region -bedFormat $bed_format -fdr $fdr -adjWindow 1000000 -minGCH 1 -beta\n";
	system($cmd)==0 || die "Unexpected stop when training HMM: $! \n";
}

sub decoding_hmm{
	my $prefix=shift @_;
	my $input=shift @_;
	my $hmm= $result_dir."/${prefix}.trainedHMM.model.txt";
	$prefix = $result_dir."/$prefix";
	my $cmd="java -Xmx${mem}g -jar $bistools_path/Bis-seg/NdrHmmHunter.jar $prefix $input -hmmFile $hmm -onlyDecode -sigTestMode $sig_test -L $region -bedFormat $bed_format -fdr $fdr -adjWindow 1000000 -minGCH 1 -beta\n";
	system($cmd)==0 || die "Unexpected stop when decoding HMM: $! \n";
}


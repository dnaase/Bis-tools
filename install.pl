#!/usr/bin/perl -w
##This script is used for quality control on WGBS/NOMe-seq.


## author: Yaping Liu  lyping1986@gmail.com 
## time: 2014-7-15

#Usage: perl install.pl

use strict;
use Getopt::Long;
use File::Basename;


###check perl environment


###check java environment


###configure $BISTOOLS and $PATH environment parameters
my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
if($bistools_path eq ""){
	my $dir = `pwd`;
	chomp($dir);
	if( -e $ENV{"HOME"}."/.bash_profile"){
		my $cmd = "echo \$BISTOOLS\n";
		$cmd = "echo \"export BISTOOLS=$dir\" >> ~/.bash_profile\n";
		print STDERR $cmd;
		system($cmd)==0 || die "can't add BISTOOLS variable into bash profile:$!\n";
		

	}else{
		my $cmd = "echo \"export BISTOOLS=$dir\" > ~/.bash_profile\n";
		print STDERR $cmd;
		system($cmd)==0 || die "can't write BISTOOLS variable into bash profile:$!\n";
	}
	my $cmd = "source $ENV{\"HOME\"}/.bash_profile\n";
	print STDERR $cmd;
	system($cmd)==0 || die "can't excute bash profile:$!\n";
}
$bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
### chmod the script to be excutable


##########install dependent external software
install_fastqmcf();


sub install_fastqmcf{
	#my $cmd = "tar -xvzf $bistools_path/External_tools/ea-utils/ea-utils.1.1.2-484.tar.gz -C $bistools_path/External_tools/ea-utils/  \n";
	#print STDERR $cmd;
	#system($cmd)==0 || die "can't unzip ea-utils:$!\n";
	`tar -xvzf $bistools_path/External_tools/ea-utils/ea-utils.1.1.2-484.tar.gz -C $bistools_path/External_tools/ea-utils/`;
	`cd $bistools_path/External_tools/ea-utils/ea-utils.1.1.2-484 && PREFIX=$bistools_path/External_tools/ea-utils/ make install`;
	#`PREFIX=$bistools_path/External_tools/ea-utils/ make install`;
}


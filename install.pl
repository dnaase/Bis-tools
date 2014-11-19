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
	my $shell=$ENV{SHELL};
	print "$shell\n";
	if($shell=~/bash$/){
		my $bash=$ENV{"HOME"}."/.bash_profile";
		if( -e $bash){
			my $check=check_bistools_bash($bash);
			if( $check == 0){
				my $cmd = "echo \$BISTOOLS\n";
				$cmd = "echo \"export BISTOOLS=$dir\" >> $bash\n";
				run_cmd($cmd);
			}
		}else{
			my $cmd = "echo \"export BISTOOLS=$dir\" > $bash\n";
			run_cmd($cmd);
		}
		print STDERR "please run command:  source ~/.bash_profile\n";
	}elsif($shell=~/tcsh$/){
		my $bash=$ENV{"HOME"}."/.tcshrc";
		if( -e $bash){
			my $check=check_bistools_bash($bash);
			if( $check == 0){
				my $cmd = "echo \$BISTOOLS\n";
				$cmd = "echo \"setenv BISTOOLS \"$dir\"\" >> $bash\n";
				run_cmd($cmd);
			}
		}else{
			my $cmd = "echo \"setenv BISTOOLS \"$dir\"\" > $bash\n";
			run_cmd($cmd);
		}
		print STDERR "please run command:  source ~/.tcshrc\n";
	}
	
	
	#my $cmd = "chmod 755 $ENV{\"HOME\"}/.bash_profile\n";
	#run_cmd($cmd);
	#$cmd = "/bin/bash $ENV{\"HOME\"}/.bash_profile\n";
	#run_cmd($cmd);
	
	
}else{
	print STDERR "Bis-tools has already been installed at $bistools_path\n";
}
#$bistools_path=`echo \$BISTOOLS`;
#chomp($bistools_path);
### chmod the script to be excutable


##########install dependent external software
#install_fastqmcf();


sub install_fastqmcf{
	#my $cmd = "tar -xvzf $bistools_path/External_tools/ea-utils/ea-utils.1.1.2-484.tar.gz -C $bistools_path/External_tools/ea-utils/  \n";
	#print STDERR $cmd;
	#system($cmd)==0 || die "can't unzip ea-utils:$!\n";
	`tar -xvzf $bistools_path/External_tools/ea-utils/ea-utils.1.1.2-484.tar.gz -C $bistools_path/External_tools/ea-utils/`;
	`cd $bistools_path/External_tools/ea-utils/ea-utils.1.1.2-484 && PREFIX=$bistools_path/External_tools/ea-utils/ make install`;
	#`PREFIX=$bistools_path/External_tools/ea-utils/ make install`;
}

sub run_cmd{
	my $cmd=shift @_;
	print STDERR $cmd;
	system($cmd)==0 || die "can't excute command $cmd : $!\n";
}

sub check_bistools_bash{
	my $bash=shift @_;
	open(FH,"< $ENV{\"HOME\"}/.bash_profile") or die "can't open .bash_profile files:$!\n";
	while(<FH>){
		chomp;
		if($_=~/export BISTOOLS/){
			print STDERR "BISTOOLS variable has already been specified!\n";
			print STDERR "$_\n";
			return 1;
		}
	}
	close(FH);
	return 0;
}

my $wd=shift(@ARGV);
my $prefix=shift(@ARGV);
my $meQTL_hic=shift(@ARGV);
my $random_hic=shift(@ARGV);

##10 random pairs
my $r_cmd="/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/bin/R --no-restore --no-save --args wd=$wd prefix=$prefix meqtl_hic=$meQTL_hic random_hics=$random_hic ";

#for(my $i=1;$i<10;$i++){
#	my $random_tmp=$random_hic;
#	$random_tmp=~s/random_\d+/random_${i}/;	
#	$r_cmd.="random_hics=$random_tmp ";
	
#}
$r_cmd.="< permutation_random_plot.R";
run_cmd($r_cmd);
 
 
 sub run_cmd{
 	my $cmd=shift @_;
 	print STDERR "$cmd\n";
 	system($cmd)==0 || die "can't finish $cmd:$!\n";
 }

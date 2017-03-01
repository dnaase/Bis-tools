my $mark_file=$ARGV[0] or die;
my $recomb_map=$ARGV[1] or die;
my $wd=$ARGV[2] || "/Users/yaping/Documents/workspace/project/meQTL/recombRate/1KGV3/select/";
my $window_size=$ARGV[3] || 2000000;
my $bin_size=$ARGV[4] || 10000;  ##need to make it for 1bp bin later

my $chr="chr1";
#my $chr_start=0;
#my $chr_end=249250621;
my $chr_start=66000000;
my $chr_end=67000000;

my %hash=();
open(FH,"<$mark_file") or die;
while(<FH>){
	chomp;
	my @f=split "\t";
	if(($f[0] eq $chr) and ($f[1]>$chr_start) and ($f[2]<=$chr_end)){
		
		$f[1]=int($f[1]/$bin_size);
		$f[2]=int($f[2]/$bin_size);
		
		my $start=$f[1]<$f[2]?$f[1]:$f[2];
		my $end=$f[1]<$f[2]?$f[2]:$f[1];
		my $key="${start}-${end}";
		$hash{$key}=1;
	}
}
close(FH);

#for each window size, generate one pdf figure
for(my $figure_start = $chr_start;$figure_start<=($chr_end-$window_size);$figure_start+=$window_size){
	my $figure_end=$figure_start+$window_size;
	my $prefix="1kGv3Recomb.$chr.window-${figure_start}-${figure_end}.bin_$bin_size";
	my $pdf="$prefix.pdf";
	my $mark_pdf="$prefix.matrix_marker.$mark_file.pdf";
	###for recombination rate map
	my $matrix_file="$prefix.matrix.txt";
	open(OUT,">${wd}/$matrix_file") or die;
	my $mark_matrix_file="$prefix.matrix_marker.$mark_file.txt";
	open(MARK,">${wd}/$mark_matrix_file") or die;
	
	my @figure_matrix=();
	for(my $bin_1=$figure_start;$bin_1<=($figure_end-$bin_size);$bin_1+=$bin_size){  ##each row of matrix
		for(my $bin_2=$figure_start;$bin_2<=($figure_end-$bin_size);$bin_2+=$bin_size){  ###each column of matrix
			my $start=int($bin_1/$bin_size);
			my $end=int($bin_2/$bin_size);
			my $key="${start}-${end}";
			if(exists $hash{$key}){
				print MARK "1\t";
			}else{
				print MARK "NA\t";
			}
			
			if($bin_2<=$bin_1){
				print OUT "NA\t";
			}else{
				my $v=`bigWigSummary $recomb_map $chr $bin_1 $bin_2 1 2>/dev/null`;
				chomp($v);
				if($v eq ""){
					$v="NA";
				}
				print OUT "$v\t";
			}
			
		}
		print OUT "\n";
		print MARK "\n";
		
	}
	close(OUT);
	close(MARK);
	##
	
	my $cmd="R --no-restore --no-save --args wd=$wd markFile=$mark_matrix_file matrix=$matrix_file figure_file=$pdf figure_mark_file=$mark_pdf bin_size=$bin_size < plot_heatmap.local_norm.grey.R";
	print STDERR "$cmd\n";
	system($cmd);
		
}

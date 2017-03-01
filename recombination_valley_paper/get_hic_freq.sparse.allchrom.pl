use File::stat;
my $resolution=$ARGV[0] or die;  ##chromosome number
my $file=$ARGV[1] or die;
my $hic_dir=$ARGV[2] or die;
my $get_pair_in_same_seg=$ARGV[3];
my $sparse_line_mode=$ARGV[4];
$resolution*=1000;

my @chrs=(1..22,"X");

foreach my $chr(@chrs){
	$chr="chr$chr";
	my $hic=$hic_dir."/${chr}.sparseBed.txt";
	
	##extract bed file with this chromsome
	my $file_tmp=$file;
	$file_tmp.=".$chr";
	my $cmd = "grep '^${chr}\\b' $file > $file_tmp";
	#print STDERR "$cmd\n";
	system($cmd);
	my $sb = stat($file_tmp);
	my $size=$sb->size;
	if($size == 0){
		`unlink $file_tmp`;
		next;
	}
	
	
	##make hash
	my %hic_hash;

	if($sparse_line_mode eq ""){
		open(FH,"zcat $hic |") or die;
		my $num_seg=0;
		my $num=0;
		while(<FH>){
			chomp;
			my @f=split " ";
			if($num_seg==0){
				$num_seg=scalar(@f);
			}
		#print "$f[0]\t$f[1]\t$f[2]\n";
			my $in=0;
			foreach my $n(@f){
				if($n>0){
			
					my $key=$num.":".$in;
					$hic_hash{$key}=$n;
					#print STDERR "$key\t$n\t$hic_hash{$key}\n";
				}
				$in++;
			}
		$num++;
		}	
		close(FH);
	}else{
		open(FH,"<$hic") or die;
		while(<FH>){
			chomp;
			my @f=split "\t";
			my $key=$f[0];
			$hic_hash{$key}=$f[1];
		}	
		close(FH);
	}
	
	###
	open(FH,"<$file_tmp") or die;
	my $pairs=0;
	while(<FH>){
		chomp;
		my @f=split "\t";
		my $start_seg=int($f[1]/$resolution);
		my $end_seg=int($f[2]/$resolution);
		$pairs++;
		my $key=$start_seg.":".$end_seg;
		if($get_pair_in_same_seg ne "" and abs($end_seg-$start_seg)<=$get_pair_in_same_seg){
				#next;
				print join("\t",@f)."\tNA\n";
				next;
		}
		if(exists $hic_hash{$key}){
			print join("\t",@f)."\t".$hic_hash{$key}."\n";
		}else{
			print join("\t",@f)."\t0\n";
		}
	}
	close(FH);

	print STDERR "$chr : Pairs: $pairs\n";
	`unlink $file_tmp`;
}



#print STDERR "ok\n";




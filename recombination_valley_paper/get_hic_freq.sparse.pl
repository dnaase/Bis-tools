#my $ref=$ARGV[0] or die;  ##hg19.fa.fai file or chromesize file
my $resolution=$ARGV[0] or die;  ##chromosome number
my $file=$ARGV[1] or die;
my $hic=$ARGV[2] or die;
my $get_pair_in_same_seg=$ARGV[3];
my $sparse_line_mode=$ARGV[4];
$resolution*=1000;

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

#print STDERR "ok\n";
open(FH,"<$file") or die;
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
	#die "resolution number is wrong!\n" if($num_seg<$end_seg);
	

}
close(FH);

print STDERR "Pairs: $pairs\n";



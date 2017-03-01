my $in=$ARGV[0] or die;
my $top=$ARGV[1] || 1;
my $unsort="";
my $chr="";
if($in=~/(chr\w+)_1kb.OE.KRnorm_raw_vs_expect.sparseBed.txt/){
	$chr=$1;
	$unsort="GM12878_1kb_OE_KRnorm_raw_vs_expect.hg19.$chr.distance1M.top${top}perc_tmp.bed";
	open(IN,"<$in");
	open(O,">$unsort");
	while(<IN>){
		chomp;
		@f=split "\t";
		@s=split ":",$f[0];
		if((abs($s[1]-$s[0])>1) && (abs($s[1]-$s[0])<=1000)){
			if($s[1]>$s[0]){
				$start=$s[0];$end=$s[1];
			}else{
				$start=$s[1];$end=$s[0];
			}
			$start*=1000;$end*=1000;
			print O "$chr\t$start\t$end\t.\t$f[1]\t$chr:$start-$end\n";
		}
	}
	close(IN);
	close(O);
	my $num=`wc -l $unsort`;
	chomp($num);
	$num=int($num*$top/100);
	my $sort=$unsort;
	$sort=~s/\.top${top}perc_tmp.bed$/.top${top}perc.bed/;
	`sort -k5g,5g $unsort | tail -n $num > $sort`;
	`unlink $unsort`;
}



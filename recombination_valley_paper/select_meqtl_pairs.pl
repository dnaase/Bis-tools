my $input=$ARGV[0] or die;
my $mode=$ARGV[1] || 1;

my %cg_hash=();
open(F,"<$input") or die;
while(<F>){
	chomp;
	my @f=split "\t";
	my $key=$f[4];
	if(exists $cg_hash{$key}){
		my @lines=split "\t", $cg_hash{$key};
		if($f[7] < $lines[7]){
			$cg_hash{$key}=join("\t",@f);  ##use the best meqtl
		}
		
	}else{
		$cg_hash{$key}=join("\t",@f);
	}
	
}
close(F);

foreach my $key(keys %cg_hash){
	my @f=split "\t",$cg_hash{$key};
	print "$f[0]\t$f[1]\t$f[2]\t.\t$f[7]\t$f[0]:$f[1]-$f[2]\n";
}
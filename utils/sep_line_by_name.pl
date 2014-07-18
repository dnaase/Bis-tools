my $file=$ARGV[0] || die "no input file!\n";
my $col=$ARGV[1] || die "need to give the column number for seperate line!\n";
my $prefix=$file;
$prefix=~s/\.\w+$//;
my %names;
open(FH,"<$file") or die "can't find $file\n";
while(<FH>){
	chomp;
	my $line=$_;
	my @splitin=split "\t",$line;
	my $key = $splitin[$col];
	$key=~s/\W+/_/g;
	if(exists $names{$key}){
		$names{$key}.="$line\n";
	}
	else{
		$names{$key}="$line\n";
	}
}
close(FH);

foreach my $key(keys(%names)){
	my $out=$prefix.".${key}";
	print STDERR "$out\n";
	open(OUT,">$out") or die "can't write $out \n";
	my $content=$names{$key};
	print OUT "$content";
	close(OUT);
}

##perl /home/unix/lyping/compbio/data/public_data/HiC_2014/normalize_hic2014_to_sparseBed.pl /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.KRnorm /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.RAWobserved > /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.KRnorm_raw.sparseBed.txt
##perl /home/unix/lyping/compbio/data/public_data/HiC_2014/normalize_hic2014_to_sparseBed.pl /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.KRnorm /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.RAWobserved 1000 /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.KRexpected > /home/unix/lyping/compbio/data/public_data/HiC_2014/GM12878/GM12878_combined/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.OE.KRnorm_raw_vs_expect.sparseBed.txt
my $norm_file=$ARGV[0] or die "no normalization vector file:$!\n";
my $raw_file=$ARGV[1] or die "no raw hic file:$!\n";
my $resolution=$ARGV[2] || 1000;
my $norm_expect_file=$ARGV[3];

##open norm file
open(FH,"<$norm_file") or die;
my $line=1;
my %norm_vector=();
while(my $in = <FH>){
	chomp($in);
	if($in ne "NaN"){
		$norm_vector{$line}=$in;
	}
	$line++;
}
close(FH);

my %norm_expect_vector=();
if($norm_expect_file ne ""){
	open(FH,"<$norm_expect_file") or die;
	my $line=1;
	while(my $in = <FH>){
		chomp($in);
		$norm_expect_vector{$line}=$in;
		$line++;
	}
	close(FH);
}

open(FH,"<$raw_file") or die;
while(<FH>){
	chomp;
	my @f=split "\t";
	my $s1=$f[0]/$resolution + 1;
	my $s2=$f[1]/$resolution + 1;
	my $expect_key=($f[1]-$f[0])/$resolution + 1;
	if((exists $norm_vector{$s1}) and (exists $norm_vector{$s2})){
		my $v=$f[2]/($norm_vector{$s1}*$norm_vector{$s2});
		$s1--;
		$s2--;
		if($norm_expect_file eq ""){
			print "$s1:$s2\t$v\n";
		}else{
			$v=$v/$norm_expect_vector{$expect_key};
			print "$s1:$s2\t$v\n";
		}
		
	}
}
close(FH);

my $input=$ARGV[0] or die;
my $tads=$ARGV[1] or die;

my $tmp_key="";
open(F,"intersectBed -loj -a $input -b $tads |") or die;
while(my $in=<F>){
	chomp($in);
	my @f=split "\t",$in;
	if($f[10] eq "-1"){
		next;
	}
	if($tmp_key eq $f[5]){next;}
	my $len=$f[2]-$f[1];
	my ($random_number_start, $random_number_end) = &creat_random_location($f[7],$f[8],$len);
	if($random_number_start == -1){next;}
	print "$f[0]\t$random_number_start\t$random_number_end\t.\t.\t$f[5]\n";
	$tmp_key=$f[5];
}
close(F);






sub creat_random_location{
        my $tads_start = shift(@_);
        my $tads_end = shift(@_);
        my $block_length = shift(@_);
        
    	my $tads_range = $tads_end-$tads_start;
       
       	if($block_length>$tads_range){
       		return (-1,-1);
       	}else{
       		my $random_number_start = int(rand($tads_range)) + $tads_start;
       		my $random_number_end = $random_number_start+$block_length;
       		while($random_number_end>$tads_end){
       			$random_number_start = int(rand($tads_range)) + $tads_start;
       			$random_number_end = $random_number_start+$block_length;
       		}
       		return ($random_number_start, $random_number_end);
       	}
		
}
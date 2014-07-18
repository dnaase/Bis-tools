## transfer gtf to bed, with or without upstream/downstream extention, not translate contig with chr_gl001...
## author: Yaping Liu  lyping1986@gmail.com

$gtf_file=$ARGV[0];
my $upstream = $ARGV[1];
my $downstream = $ARGV[2];
my $use_age = "USAGE: perl gtf2bed.pl input_gtf [upstream_extension] [downstream_extension] > output_bed";
if($ARGV[0] eq ""){
	print "$use_age\n";
	exit(1);
}

open(FH,"<$gtf_file") or die;
while(<FH>){
	chomp;
	my $line=$_;
	my @splitin = split "\t",$line;
	next if($splitin[0]=~/\_/ or $splitin[0]!~/^chr/);
	my $start = $splitin[3]-1;
	my $end = $splitin[4];
	if($splitin[6] eq '-'){
		$start -= $downstream;
		$end += $upstream;
		if($start <= 0){
			$start = $splitin[3];
		}
	}
	else{
		$start -= $upstream;
		$end += $downstream;
		if($start <= 0){
			$start = $splitin[3];
		}
	}
	print "$splitin[0]\t$start\t$end\t$splitin[1]\t$splitin[5]\t$splitin[6]\n";
	#	push(@out,"$splitin[0]\t$splitin[3]\t$splitin[4]\texon\t$splitin[5]\t$splitin[6]\n");
	
} 
close(FH);




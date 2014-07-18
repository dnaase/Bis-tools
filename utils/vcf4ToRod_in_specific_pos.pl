## transfer VCF 4.1 to file in ROD format in provided position list
## author: Yaping Liu  lyping1986@gmail.com


my $input_file_name = $ARGV[0];
my $pos_file_name = $ARGV[1];


my $use_age = "USAGE: perl vcf4ToRod.pl input_file_name pos_file_name";
if($input_file_name eq NULL){
	print "$use_age\n";
	exit(1);
}

my $output_file_name = $input_file_name.".rod";
open(OUT,">$output_file_name") or die;

open(FH,"<$pos_file_name") or die;
my @pos=<FH>;
chomp(@pos);
close(FH);

my %pos_hash = ();
foreach my $line(@pos){
	#next if $line=~/^\#/;
	my @splitvcf=split "\t",$line;
	my $chr = "chr".$splitvcf[1];
	my $key = $chr."~".$splitvcf[2];
	#next if $splitvcf[5] <=10;
	$pos_hash{$key} = $line;
#	print "$key\n";
}
@pos=NULL;
print "finish hashing\n";

open(FH,"<$input_file_name") or die;
while(<FH>){
	chomp;
	my $line=$_;
	next if $line =~ /^#/;
	
	my @splitline=split "\t",$line;
	
	my $chr = $splitline[0];
	my $start = $splitline[1]-1;
	my $end = $splitline[1];
	my $rsID = $splitline[2];
	#print "$end\t";
	my $key = $chr."~".$end;
	if(!exists $pos_hash{$key}){
		next;
	}
	#print "$key\n";
	my $ref = $splitline[3];
	next if $splitline[4] eq '.';
	my $oberserve;
	if(($splitline[3] eq 'G' and $splitline[4] eq 'A') or ($splitline[3] eq 'T' and $splitline[4] eq 'C') or ($splitline[3] eq 'C' and $splitline[4] eq 'A') or ($splitline[3] eq 'T' and $splitline[4] eq 'G') or ($splitline[3] eq 'G' and $splitline[4] eq 'C') or ($splitline[3] eq 'T' and $splitline[4] eq 'A')){
		$oberserve = "$splitline[4]\/$splitline[3]";
	}
	else{
		$oberserve = "$splitline[3]\/$splitline[4]";
	}

	
	my $valid = "unknown";
	my $weight;
	if($splitline[7] =~ /VLD/){
		$valid = "by-cluster";
	}
	if($splitline[7] =~ /WGT=(\d+);/){
		$weight = $1;
	}

	my $output_line = "1\t$chr\t$start\t$end\t$rsID\t0\t+\t$ref\t$ref\t$oberserve\tgenomic\tsingle\t$valid\t0\t0\tunknown\texact\t$weight\n";
	print OUT $output_line;
}

close(FH);
close(OUT);

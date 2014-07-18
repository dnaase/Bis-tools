## transfer TCGA VCF 4.1 to file in bed 6+6 format(TCGA level II data): chr + start + end + methylation_value(100%) + C/T reads number + strand + thick_start + thick_end + rgb + A_reads_number_in_non_bisulfite_conversion_strand + G_reads_number_in_non_bisulfite_conversion_strand
## numC+numT==0 sites are not included. need to clarify cContext and SNP call's meaning from ben? This is wrong, since bedDetail in UCSC only allow bedPlus(4-12) + 2 more additional columns..
## author: Yaping Liu  lyping1986@gmail.com


my $input_file_name = $ARGV[0];
my $type = $ARGV[1];
#my $combine = $ARGV[2];
#my $sample_name = $ARGV[2];

my $use_age = "USAGE: perl vcf2bed.pl input_file_name [CG]";
if($ARGV[0] eq ""){
	print "$use_age\n";
	exit(1);
}
my $cpg_name_output = $input_file_name;
$cpg_name_output =~ s/\.vcf//;
open(FH,"<$input_file_name") or die;
my $bissnp_version;
while(<FH>){
	my $line=$_;
	chomp($line);
	if($line =~ /^##(BisSNP-\S+) Program Args=/){
		$bissnp_version = $1;
		last;
	}
	
}
close(FH);
$cpg_name_output = $cpg_name_output.".$bissnp_version.$type.bed";
open(OUT,">$cpg_name_output") or die;
my $descript;
if($type eq "CG"){
	$type = "CG"; 
	$descript = "CG"; ##only output homozygous CpG
}
else{
	if($type eq ""){
		$descript = "Cytosine";
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
	else{
		$descript = $type;
	}
	
}

my $head_line = "track name=$cpg_name_output type=bedDetail description=\"$descript methylation level\" visibility=3";
#variableStep chrom=chr19 span=150"
print OUT "$head_line\n";
my @palette=("0,240,0", "30,210,0","60,180,0","90,150,0","120,120,0","150,90,0","180,60,0","210,0,0");

open(FH,"<$input_file_name") or die;
while(<FH>){
	$line=$_;
	chomp($line);
	next if $line =~ /^#/;
	my @splitin = split "\t", $line;
	next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity");

	if($splitin[9] =~ /\d+,\d+,\d+,(\d+),(\d+),\d+:(\d+):($type):(\d+):/){
			#print "ok";
		my $num_g = $1;
		my $num_a = $2;
		my $num_c = $3;
		my $context = $4;
		my $num_t = $5;
		my $strand=".";
		
		#my $SNPcall = $context;
		
		if($splitin[7] =~ /CS=([+|-])/){
			$strand = $1;
		}
		if($num_c + $num_t != 0){
			$methy = $num_c/($num_c + $num_t);
			my $rgb = $palette[int($methy*8-0.0001)];
			$methy = sprintf("%.2f",100*$methy);
			my $chr = $splitin[0];
			my $start = $splitin[1]-1;
			my $end = $splitin[1];
		 	my $reads_num = $num_c + $num_t;
		 	
			my $out_line = "$chr\t$start\t$end\t$methy\t$reads_num\t$strand\t$start\t$end\t$rgb\t$num_a\t$num_g";
			print OUT "$out_line\n";
		}
		
	}
	
}

close(FH);
close(OUT);
print "finished!\n";


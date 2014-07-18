## transfer TCGA VCF 4.1 to file in bed 6+2 format(TCGA level III data): chr + start + end + . + methylation_value*10 + strand + methylation_value(%) + numCT_reads
## numC+numT==0 sites are not included.Only homozygous CG if pattern specified as CG.  This program gives +/- strand seperately.
## author: Yaping Liu  lyping1986@gmail.com


my $input_file_name = $ARGV[0];
my $type = $ARGV[1];
#my $sample_name = $ARGV[2];

my $use_age = "USAGE: perl vcf2bed6plus2.strand.pl input_file_name [CG]";
if($ARGV[0] eq ""){
	print "$use_age\n";
	exit(1);
}
my $cpg_name_output = $input_file_name;
$cpg_name_output =~ s/\.vcf//;
$cpg_name_output = $cpg_name_output.".$type.strand.6plus2.bed";
open(OUT,">$cpg_name_output") or die;

if($type eq "CG"){
	$type = "CG"; 
}
else{
	if($type eq ""){
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
}

my $head_line = "track name=$cpg_name_output type=bedDetail description=\"methylation level\" visibility=3";
#variableStep chrom=chr19 span=150"
print OUT "$head_line\n";

open(FH,"<$input_file_name") or die;
while(<FH>){
	$line=$_;
	chomp($line);
	next if $line =~ /^#/;
	my @splitin = split "\t", $line;
	next unless ($splitin[6] eq "PASS");

	if($splitin[9] =~ /:(\d+):($type):(\d+):/){
		my $num_c = $1;
		my $context = $2;
		my $num_t = $3;

		my $strand=".";
		
		if($splitin[7] =~ /CS=([+|-])/){
			$strand = $1;
		}
		if($num_c + $num_t != 0){
			$methy = $num_c/($num_c + $num_t);
			my $score = sprintf("%d",1000*$methy);
			$methy = sprintf("%.2f",100*$methy);
			
			my $chr = $splitin[0];
			my $start = $splitin[1]-1;
			my $end = $splitin[1];
		 	my $reads_num = $num_c + $num_t;
			my $out_line = "$chr\t$start\t$end\t.\t$score\t$strand\t$methy\t$reads_num";
			print OUT "$out_line\n";
		}
		
	}
	
}

close(FH);
close(OUT);
print "finished!\n";


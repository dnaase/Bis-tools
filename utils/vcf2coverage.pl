## transfer TCGA VCF 4.1 to file in bed 6+6 format(TCGA level II data): chr + start + end + C/T reads covergaes
##forward and reverse strand are combined
## author: Yaping Liu  lyping1986@gmail.com


my $input_file_name = $ARGV[0];
my $type = $ARGV[1];

my $use_age = "USAGE: perl vcf2coverae.pl input_file_name [CG]\n Output bedGraph format to represent C/T reads covergae\n";
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
$cpg_name_output = $cpg_name_output.".$bissnp_version.$type.coverage.bedgraph";
open(OUT,">$cpg_name_output") or die;
#my $descript;
if($type eq "CG"){
	$type = "CG"; 
#	$descript = "CG"; ##only output homozygous CpG
}
else{
	if($type eq ""){
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
}

#my $head_line = "track name=$cpg_name_output type=bedDetail description=\"$descript methylation level\" visibility=3";
#variableStep chrom=chr19 span=150"
#print OUT "$head_line\n";


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
		
		#my $SNPcall = $context;
		
		if($splitin[7] =~ /CS=([+|-])/){
			$strand = $1;
			
		}
		
			my $reads = $num_c + $num_t;
			my $chr = $splitin[0];
			my $start = $splitin[1]-1;
			my $end = $splitin[1];
			my $out_line = "$chr\t$start\t$end\t$reads";
			print OUT "$out_line\n";

		
	}
	
}

close(FH);
close(OUT);
print "finished!\n";


## transfer TCGA VCF 4.1 to file in bedGraph format: chr + start + end + methylation_value(%) + C/T reads number
##forward and reverse strand are combined
## author: Yaping Liu  lyping1986@gmail.com

use Getopt::Long;

my $sample_order=1;
my $min_ct=1;

sub usage { ### not use these filters yet...
	print "USAGE: perl vcf2bedGraph.pl [options] input_file_name [CG]\n";
	print "  [Options]:\n\n";
	print "  --sample_order INT : sample order choose to output as bed file when there are multiple samples in the same vcf file. 1 means the first sample (default: 1)\n\n";
	print "  --minCT INT : minimum number of CT reads, otherwise, methy column will be '.' (default: 1)\n\n";
	exit(1);
}

GetOptions(

	"sample_order=i" => \$sample_order,
	"minCT=i" => \$min_ct,
);



my $input_file_name = $ARGV[0];
my $type = $ARGV[1];

my $use_age = "USAGE: perl vcf2bedGraph.pl input_file_name [CG]";
if($ARGV[0] eq ""){
	&usage();
	#exit(1);
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
#$cpg_name_output = $cpg_name_output.".$bissnp_version.$type.bedgraph";
if($min_ct>1){
	$cpg_name_output = $cpg_name_output.".$type.minCT$min_ct.bedgraph";
}else{
	$cpg_name_output = $cpg_name_output.".$type.bedgraph";	
}

open(OUT,">$cpg_name_output") or die;
#my $descript;
if($type eq "CG"){
	$type = "CG"; 
#	$descript = "CG"; ##only output homozygous CpG
}
else{
#	$descript = "Cytosine";
	if($type eq ""){
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
}

my $head_line = "track type=bedGraph name=${cpg_name_output}.${bissnp_version}  description=\"$type methylation level\" visibility=3";
#variableStep chrom=chr19 span=150"
print OUT "$head_line\n";


open(FH,"<$input_file_name") or die;
while(<FH>){
	my $line=$_;
	chomp($line);
	next if $line =~ /^#/;
	my @splitin = split "\t", $line;
	next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity");
	my $sample_id = 8 + $sample_order;
	if($sample_id > $#splitin){
		print "--sample_order is larger than the number of samples in VCF files\n\nExit with error!!\n\n ";
		exit(1);
	}
	if($splitin[$sample_id] =~ /:(\d+):($type):(\d+):/){
		my $num_c = $1;
		my $context = $2;
		my $num_t = $3;
		my $strand=".";
		
		#my $SNPcall = $context;
		
		if($splitin[7] =~ /CS=([+|-])/){
			$strand = $1;
			
		}
		if($num_c + $num_t >= $min_ct){
			$methy = $num_c/($num_c + $num_t);
			$methy = sprintf("%.2f",100*$methy);
			my $chr = $splitin[0];
			my $start = $splitin[1]-1;
			my $end = $splitin[1];
			my $ct_reads = $num_c + $num_t;
			#my $out_line = "$chr\t$start\t$end\t$methy\t$ct_reads";
			my $out_line = "$chr\t$start\t$end\t$methy";
			print OUT "$out_line\n";
		}
		
	}
	
}

close(FH);
close(OUT);
print "finished!\n";

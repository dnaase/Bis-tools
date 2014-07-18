## transfer TCGA VCF 4.1 to file in wiggle format: coordinate + CT reads coverage
##forward and reverse strand are combined
## author: Yaping Liu  lyping1986@gmail.com


my $input_file_name = $ARGV[0];
my $type = $ARGV[1];

my $use_age = "USAGE: perl vcf2wig_ct_coverage.pl input_file_name [CG]";
if($ARGV[0] eq ""){
	print "$use_age\n";
	exit(1);
}
my $cpg_name_output = $input_file_name;
$cpg_name_output =~ s/\.vcf//;
$cpg_name_output = $cpg_name_output.".$type.ct_coverage.wig";
open(OUT,">$cpg_name_output") or die;
#my $descript;
if($type eq "CG"){
	$type = "CG"; 
	$descript = $input_file_name; ##only output homozygous CpG
	$descript =~ s/.vcf/.coverage/;
}
else{
	$descript = $input_file_name;
	$descript =~ s/.vcf/.coverage/;
	#$type = "\\w+"; ##by default, output all sites in VCF file
	if($type eq ""){
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
}

#my $head_line = "track name=$cpg_name_output type=bedDetail description=\"$descript methylation level\" visibility=3";
my $head_line = "track type=wiggle_0 name=".$input_file_name." description=\"$descript CT reads coverage\"";
print OUT "$head_line\n";


my $chr_flag = "";
my $pre = "";
open(FH,"<$input_file_name") or die;
while(<FH>){
	$line=$_;
	chomp($line);
	next if $line =~ /^#/;
	my @splitin = split "\t", $line;
	next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity");
	if($splitin[0] ne $chr_flag){
			$chr_flag = $splitin[0];
			my $chr_head_line = "variableStep chrom=$chr_flag";
			print OUT "$chr_head_line\n";
		}

	if($splitin[9] =~ /:(\d+):($type):(\d+):/){
		my $num_c = $1;
		my $context = $2;
		my $num_t = $3;


		if($num_c + $num_t != 0){
			$cov = $num_c + $num_t;
			my $cor = $splitin[1];
			if($cor eq $pre){
				next;
			}
			my $out_line = "$cor\t$cov";
			print OUT "$out_line\n";
			$pre = $cor;
		}
		
	}
	
}

close(FH);
close(OUT);
print "finished!\n";


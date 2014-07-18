#!/usr/bin/perl

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-2-19

#Usege:  bissnp_trinuc_sample.pl $outputfile bam_file 

use strict;
use Getopt::Long;

my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
my $BISSNP = "$bistools_path/Bis-SNP/Bis-SNP.latest.jar";
my $dbsnp = "$bistools_path/resource/dbSNP/dbsnp_135.hg19.sort.vcf";
my $ref="$bistools_path/resource/genome/hg19_rCRSchrm.fa";
my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
my $mem="8"; #how many Giga bytes memory need
my $interval = "chrM";
my $not_use_bad_mates="";

GetOptions(
	"bissnp=s" => \$BISSNP,
	"genome=s" => \$ref,
	"dbsnp=s" => \$dbsnp,
	"nt=i" => \$numcores,
	"mem=i" => \$mem,
	"interval=s" => \$interval,
	"not_use_bad_mates" => \$not_use_bad_mates,
);



my $metric_unsorted_cpg = $ARGV[0] || die "need output log file";
my $input = $ARGV[1] || die "need input bam file";


my $minPatConv=0.8;


my $vcf_unsorted_cpg=$input;
$vcf_unsorted_cpg =~ s/\.bam$/.cytosine.raw.vcf/;
my $vcf_summary=$vcf_unsorted_cpg.".MethySummarizeList.txt";

my @characters=("A","C","G","T");


&bissnp();
&output_tri_file();
`rm $vcf_unsorted_cpg`;
`rm $vcf_summary`;

sub bissnp{
	my $cmd .= "java -Xmx${mem}G -jar $BISSNP -R $ref ";
	$cmd .= "-I $input ";									
	$cmd .= "-T BisulfiteGenotyper -vfn1 $vcf_unsorted_cpg ";
	
	$cmd .= "-L $interval " if $interval;
	$cmd .= "-D $dbsnp " if $dbsnp;
	$cmd .= "-out_modes EMIT_ALL_CYTOSINES -sm BM -minPatConv $minPatConv ";
	
	my $c_str="";
	foreach my $a(@characters){
		#$c_str.="";
		foreach my $b(@characters){
			$c_str.="-C ${a}C${b},2 ";	
		}
	}
	$cmd .= "$c_str";
	$cmd .= "-toCoverage 99999999 " if ($interval eq "chrM" || $interval eq "MT" || $not_use_bad_mates eq "");
	$cmd .= "-stand_call_conf 20 -stand_emit_conf 0 -nt $numcores -minConv 1 \n";

	print STDERR "$cmd\n";
	system($cmd)==0 || die "can't generate trinucleotide methylaiton level in chromesome $interval:$!\n";
	
}


sub output_tri_file{
	open(FH,"<$vcf_summary") or die "can't open $vcf_summary";
	open(OUT,">$metric_unsorted_cpg") or die "can't open $metric_unsorted_cpg";
	my $out_flag=-1;
	while(<FH>){
		chomp;
		print OUT "$_\n" if($out_flag==1);
		$out_flag=1 if($_=~/^##Methylation summary in Read Group:/);
	}
	close(FH);
	close(OUT);
}


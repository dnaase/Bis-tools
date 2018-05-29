#!/usr/bin/perl -w
## make easy usage of BisSNP
## author: Yaping Liu  lyping1986@gmail.com



use strict;
use Getopt::Long;




sub usage {

    print "\nUsage:\n";
    print "perl bissnp_easy_usage.pl [Options] BISSNP_JAR INPUT_BAM REF_GENOME DBSNP_VCF\n\n";

    print " Provide an easy way to operate BisSNP in one command line\n";
    print " It would automately generate: \n";
    print " 1) cpg.vcf (detailed CpGs methylation and genotyping information)\n";
    print " 2) snp.vcf (detailed SNPs information)\n";
    print " 3) cpg.bed (CpG methylation bed file) and cpg.coverage.bed (CpG C/T reads coverage bedGraph file)\n";
    print " 4) realigned and mark duplicated reads BAM file (if local indel realignment option is enabled)\n";
    print " 5) realigned, mark duplicated reads, base quality score recalibrated BAM file (if base quality score recalibration option is enabled)\n\n";

    
    print " All the BisSNP scripts should be put under the same directory. \n";
    print " It include the following steps:\n";
    print " 1). Indel realignment (Default: not enabled since its time consuming)\n";
    print " 2). Mark duplicated reads (Default: not enabled, if it is enabled, need to provide the path of MarkDuplicates.jar in picard tools)\n";
    print " 3). base quality recalibration (Default: not enabled since its time consuming)\n";
    print " 4). genotyping, get cpg.vcf and snp.vcf file\n";
    print " 5). filter out fake SNPs inside cpg.vcf and snp.vcf file\n";
    print " 6). generate cpg.methylation.bed(only homozygous CpGs, CpG adjacent in two strands are combined), \n";
    print " 	cpg.coverage.bedgraph\n\n";

	print "  BISSNP_JAR     Bis-SNP .jar name.\n";
    print "  INPUT_BAM     input BAM files to do genotyping. It should contain ReadGroup tag inside BAM.\n";
    print "  REF_GENOME    .fa/.fasta file, the reference genome file which should be the \n";
    print "              same when you do reads mapping to obtain BAM files.\n";
    print "  DBSNP_VCF     dbSNP.vcf file used for genotyping and base quality recalibration. Its chromosome order should be the same as REF_GENOME.\n\n";
	
	print "  [Options]:\n\n";
    print "  --nomeseq : enable NOMe-seq mode\n\n";
    print "  --rrbs RRBS_READ_GROUPD_ID: enable RRBS mode only to the reads group that specified\n\n";
    print "  --calmd :   enable MD tag calculation, which will benefit bisulfite incomplete conversion filter\n\n";
    print "  --lowCov :  enable low coverage mode (--qual default value will be set to 10)\n\n";
    print "  --indel :   enable indel local realignment\n\n";
    print "  --knownIndel INDEL_FILE : indel.vcf file name, when --indel enabled, this option need to be provided\n\n";
    print "  --recal :   enable base quality recalibration\n\n";
    print "  --duplicate DIR : directory of MarkDuplicates.jar in picard tools\n\n";
	print "  --qual NUM : genotyping quality score (Default: 20)\n\n";
	print "  --interval FILE : interval .bed file name\n\n";
	print "  --minCT NUM : minimum number of C or T reads to call methylation (Default: 1)\n\n";
	print "  --mmq NUM : minimum mapping quality score of reads used for genotyping (Default: 30)\n\n";
	print "  --mbq NUM : minimum base quality score of bases used for genotyping and methylation (Default: 5)\n\n";
	print "  --minConv NUM : minimum number of converted C in the begining of reads required\n";
	print "                  to begin count for genotyping and methylation (Default: 1)\n\n";
	print "  --nt NUM : number of CPU threads to use (Default: 1)\n\n";
	print "  --mem NUM : how much gigabytes memory to use (Default: 10)\n\n";
	print "  --allC : output all of cyosines in filtered VCF file into bed file, otherwise, it will only output CpG into bed file except for NOMe-seq mode\n";
	print "  --cytosines STRING : cytosine pattern to be analyzed (Default: --cytosine CG,1 --cytosine CH,1)\n";
	print "                      you can specify '--cytosine CG,1 --cytosine CHG,1 --cytosine CHH,1' to get CG, CHG and CHH information \n";
	print "                      NOMe-seq mode will change the default value to \n";
	print "                      '--cytosine GCH,2 --cytosine HCH,2 --cytosine HCG,2 --cytosine GCG,2'\n\n";
	print "  --outMode NUM : (Default: --outMode 1)\n";
	print "                 1) output CpG and SNPs\n";
	print "                 2) output Cytosine and SNPs\n";
	print "                 3) output Cytosines\n";
	print "                 4) output CpGs\n";
	print "                 5) output SNPs\n";
	print "                 6) output heterozygous SNPs\n";
	print "                 7) output all confident sites' genotype\n";
	print "                 8) output all callable sites' genotype\n";
	print "                 9) output NOMe-seq related GCH, HCG, GCG\n\n";
	print "                 10) output epistate CpG reads file\n\n";
	print "                 11) output epistate HCG and GCH reads file (for NOMe-seq)\n\n";
    exit(1);
}

##default option setting
my $nomeseq = "";
my @rrbs = ();
my $lowCov = "";
my $indel = "";
my @knownIndels = ();
my $recal = "";
my $duplicate = "";
my $qual = 20;
my $interval = "";
my $minCT = 1;
my $mmq = 30;
my $mbq = 5;
my $minConv = 1;
my $nt = 1;
my $mem = 10;
my @cytosines = ();
my $outMode = 1;
my $allC = "";
my $calmd = "";

GetOptions( "nomeseq" => \$nomeseq,
			"rrbs=s" => \@rrbs,
			"calmd" => \$calmd,
			"lowCov" => \$lowCov,
			"indel" => \$indel,
			"knownIndels=s" => \@knownIndels,
			"recal" => \$recal,
			"duplicate=s" => \$duplicate,
			"qual=i" => \$qual,
			"interval=s" => \$interval,
			"minCT=i" => \$minCT,
			"mmq=i" => \$mmq,
			"mbq=i" => \$mbq,
			"minConv=i" => \$minConv,
			"nt=i" => \$nt,
			"mem=i" => \$mem,
			"allC" => \$allC,
	    "cytosines=s" => \@cytosines,
	    "outMode=i" => \$outMode);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) != 4 ) {
    print "Wrong number of arguments\n";
    usage();
}

my $bissnp_jar = $ARGV[0];
my $input_file = $ARGV[1];
my $ref = $ARGV[2];
my $dbsnp = $ARGV[3];

my $PICARD = $duplicate;
my $BISSNP = $bissnp_jar;
my $JAVA = "java -Xmx$mem" . "G";

my $VCFTOOLS = "/home/uec-00/yapingli/code/mytools/perl/sortByRefAndCor.pl";
my $VCF2BED = "/home/uec-00/shared/production/software/bissnp/vcf2bed6plus2.pl";
my $VCF2COV = "/home/uec-00/yapingli/code/mytools/perl/vcf2coverage.pl";
my $VCF2BEDGRAPH = "/home/uec-00/yapingli/code/mytools/perl/vcf2bedGraph.pl";
my $BEDGRAPH2BIGWIG = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/bedGraphToBigWig";
my $CHROMSIZE = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/hg19.chrom.sizes";
my $VCF2WIG = "/home/uec-00/shared/production/software/bissnp/vcf2wig.pl";
my $VCF2WIGCOV = "/home/uec-00/shared/production/software/bissnp/vcf2wig_ct_coverage.pl";


unless (-e $BISSNP) {
 	print "\n $BISSNP jar file not exists!\n\n";
 	exit(1);
}

my %out_modes_hash = (
	1 => "DEFAULT_FOR_TCGA",
	2 => "EMIT_VARIANT_AND_CYTOSINES",
	3 => "EMIT_ALL_CYTOSINES",
	4 => "EMIT_ALL_CPG",
	5 => "EMIT_VARIANTS_ONLY",
	6 => "EMIT_HET_SNPS_ONLY",
	7 => "EMIT_ALL_CONFIDENT_SITES",
	8 => "EMIT_ALL_SITES",
	9 => "NOMESEQ_MODE",
);

if(scalar(@cytosines) == 0){
	@cytosines = ("CG,1", "CH,1");
}

if($allC ne ""){
		if($outMode == 1){
			$outMode = 2;
		}
		elsif($outMode == 4){
			$outMode = 3;
		}
}


if($nomeseq ne ""){
	@cytosines = ("GCH,2", "HCH,2","HCG,2", "GCG,2");
	$outMode = 9;
}

if($lowCov ne ""){
	if($qual == 20){
		$qual=10;
	}
}

## intermdediate file name 


my $tmp_dir="./";
my $original_input = $input_file;

if($indel ne ""){
	&bam_indel_realign();
}
if($duplicate ne ""){
	&bam_mdups();
}
if($recal ne ""){
	&bam_base_recalibration();
}
if($calmd ne ""){
	&bam_calmd();
}

my $vcf_unsorted_cpg = $input_file;
my $vcf_unsorted_snp = $input_file;
$vcf_unsorted_cpg =~ s/\.bam//;
if($nomeseq ne "" or $outMode == 2 or $outMode == 3){
	$vcf_unsorted_cpg .= ".cytosine.raw.vcf";
}
elsif($outMode == 5 or $outMode == 6){
	$vcf_unsorted_cpg .= ".snp.raw.vcf";
}
elsif($outMode == 7 or $outMode == 8){
	$vcf_unsorted_cpg .= ".genotype.raw.vcf";
}
else{
	$vcf_unsorted_cpg .= ".cpg.raw.vcf";
}

$vcf_unsorted_snp =~ s/\.bam//;
$vcf_unsorted_snp .= ".snp.raw.vcf";

my $vcf_sorted_cpg = $vcf_unsorted_cpg;
my $vcf_sorted_snp = $vcf_unsorted_snp;
$vcf_sorted_cpg =~ s/.vcf//;
$vcf_sorted_cpg .= ".sort.vcf";
$vcf_sorted_snp =~ s/.vcf//;
$vcf_sorted_snp .= ".sort.vcf";

my $vcf_sorted_cpg_filtered = $vcf_sorted_cpg;
my $vcf_sorted_snp_filtered = $vcf_sorted_snp;
$vcf_sorted_cpg_filtered =~ s/.raw.sort.vcf//;
$vcf_sorted_cpg_filtered .= ".filtered.sort.vcf";
$vcf_sorted_snp_filtered =~ s/.raw.sort.vcf//;
$vcf_sorted_snp_filtered .= ".filtered.sort.vcf";

my $wig_cpg = $vcf_sorted_cpg_filtered;
my $wig_cpg_cov = $vcf_sorted_cpg_filtered;
my $wig_hcg = $vcf_sorted_cpg_filtered;
my $wig_hcg_cov = $vcf_sorted_cpg_filtered;
my $wig_hcg_tdf;
my $wig_hcg_tdf_cov;
$wig_cpg =~ s/.vcf//;
if($nomeseq ne ""){
	$wig_cpg .= ".GCH.wig";
	$wig_cpg_cov =~ s/.vcf//;
	$wig_cpg_cov .= ".GCH.ct_coverage.wig";
	
	$wig_hcg =~ s/.vcf//;
	$wig_hcg .= ".HCG.wig";
	$wig_hcg_cov =~ s/.vcf//;
	$wig_hcg_cov .= ".HCG.ct_coverage.wig";
	
	$wig_hcg_tdf = $wig_hcg;
	$wig_hcg_tdf =~ s/.wig//;
	$wig_hcg_tdf_cov = $wig_hcg_cov;
	$wig_hcg_tdf_cov =~ s/.wig//;
	$wig_hcg_tdf .= ".tdf";
	$wig_hcg_tdf_cov .= ".tdf";

}
else{
	$wig_cpg .= ".CG.wig";
	$wig_cpg_cov =~ s/.vcf//;
	$wig_cpg_cov .= ".CG.ct_coverage.wig";
}

my $wig_tdf = $wig_cpg;
$wig_tdf =~ s/.wig//;
my $wig_tdf_cov = $wig_cpg_cov;
$wig_tdf_cov =~ s/.wig//;



my $genome_version;
if($ref =~/hg18/){
	$genome_version = "hg18";
	$CHROMSIZE = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/hg18.chrom.sizes";
}
elsif($ref =~/hg19/){
	$genome_version = "hg19";
}
elsif($ref =~/37/){
	$genome_version = "b37";
	
}
elsif($ref =~/mm9/){
	$genome_version = "mm9";
	$CHROMSIZE = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/mm9.chrom.sizes";
}
$wig_tdf .= ".tdf";
$wig_tdf_cov .= ".tdf";

&bissnp();
&vcf_sort();
my $cmd="rm $vcf_unsorted_cpg";
system($cmd);
$cmd="rm $vcf_unsorted_snp";
system($cmd);
if($outMode != 3 and $outMode != 4 and $outMode != 7 and $outMode != 8 ){
	&vcf_filter();
}
&vcf2bigwig($vcf_sorted_cpg_filtered);
#&vcf2wig();
&vcf2bed();



sub bam_indel_realign{
	if(scalar(@knownIndels) == 0){
		die "  Enabled indel realignemnt steps, but not provide known indels vcf file for the alignment!\n\n";
	}
	my $indel_target_interval = $input_file;
	$indel_target_interval =~ s/\.bam//;
	$indel_target_interval .= ".indels.intervals";
	my $cmd .= "$JAVA -jar $BISSNP -R $ref ";	
	$cmd .= "-I $input_file ";
	$cmd .= "-T BisulfiteRealignerTargetCreator ";
	if($interval ne ""){
		$cmd .= "-L $interval ";
	}
	foreach my $knownIndel(@knownIndels){
		$cmd .= "-known $knownIndel ";
	}
	
	$cmd .= "-o $indel_target_interval -nt $nt\n ";
										
	#realign 
	
	$cmd .= "$JAVA -jar $BISSNP -R $ref ";	
	$cmd .= "-I $input_file ";
	$cmd .= "-T BisulfiteIndelRealigner -targetIntervals $indel_target_interval ";
	foreach my $knownIndel(@knownIndels){
		$cmd .= "-known $knownIndel ";
	}
	if($lowCov ne ""){
		$cmd .= "-LOD 2.0 ";
	}
	$cmd .= "-compress 5 -cigar ";
	$input_file =~ s/\.bam//;
	$input_file .= ".realign.bam";
	$cmd .= "-o $input_file\n ";
	print "$cmd\n";
	system($cmd);
	
}

sub bam_mdups{
	my $remove_bam = $input_file;
	my $md = "$PICARD/MarkDuplicates.jar";
	unless (-e $md) {
 		print "\n MarkDuplicates.jar file not exists!\n\n";
 		exit(1);
	}
	my $cmd .= "$JAVA -jar $PICARD/MarkDuplicates.jar I=$input_file "; 
	$input_file =~ s/\.bam//;
	$input_file .= ".mdups.bam";
	$cmd .= "O=$input_file ";
	$cmd .= "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=3000000 ";
	my $metrics_file = $input_file.".dupsMetrics.txt";
	$cmd .= "METRICS_FILE=$metrics_file ";
	$cmd .= "TMP_DIR=$tmp_dir \n";
	print "$cmd\n";
	system($cmd);
	if($original_input ne $remove_bam){
		system("rm $remove_bam");
		$remove_bam =~ s/\.bam/.bai/;
		system("rm $remove_bam");
	}
	
	
}

sub bam_base_recalibration{


	my $remove_bam = $input_file;
	my $recalFile_before = $input_file;
	$recalFile_before =~ s/\.bam//;
	$recalFile_before .= ".beforeRecal.txt";
	my $recalFile_after=$input_file;
	$recalFile_after =~ s/\.bam//;
	$recalFile_after .= ".afterRecal.txt";
	
##1 countCovariant
	my $cmd .= "$JAVA -jar $BISSNP -R $ref ";	
	$cmd .= "-I $input_file ";
	$cmd .= "-T BisulfiteCountCovariates -nt $nt ";
	$cmd .= "-knownSites $dbsnp ";
	$cmd .= "-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate ";
	$cmd .= "-recalFile $recalFile_before\n ";

##2 TableRecalibration
	$cmd .= "$JAVA -jar $BISSNP -R $ref ";	
	$cmd .= "-I $input_file ";
	$cmd .= "-T BisulfiteTableRecalibration ";
	$input_file =~ s/\.bam//;
	$input_file .= ".recal.bam";
	$cmd .= "-o $input_file ";
	$cmd .= "-recalFile $recalFile_before -maxQ 40\n ";

##3 countCovariantAfterRecalibrate
	$cmd .= "$JAVA -jar $BISSNP -R $ref ";	
	$cmd .= "-I $input_file ";
	$cmd .= "-T BisulfiteCountCovariates -nt $nt ";
	$cmd .= "-knownSites $dbsnp ";
	$cmd .= "-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate ";
	$cmd .= "-recalFile $recalFile_after\n ";
	
	print "$cmd\n";
	system($cmd);
	if($original_input ne $remove_bam){
		system("rm $remove_bam");
		$remove_bam =~ s/\.bam/.bai/;
		system("rm $remove_bam");
	}
	
}

sub bam_calmd{
	my $remove_bam = $input_file;
	my $cmd .= "samtools calmd -b $input_file $ref 2>/dev/null > "; 
	$input_file =~ s/\.bam//;
	$input_file .= ".calmd.bam";
	$cmd .= "$input_file \n";
	print "$cmd\n";
	system($cmd);
	if($original_input ne $remove_bam){
		system("rm $remove_bam");
		$remove_bam =~ s/\.bam/.bai/;
		system("rm $remove_bam");
	}
	$cmd = "samtools index $input_file\n";
	print "$cmd\n";
	system($cmd);

}

sub bissnp{
	my $cmd .= "$JAVA -jar $BISSNP -R $ref ";
	$cmd .= "-I $input_file ";									
	$cmd .= "-D $dbsnp -T BisulfiteGenotyper -vfn1 $vcf_unsorted_cpg -vfn2 $vcf_unsorted_snp ";
	if($nomeseq ne ""){
		$cmd .= "-out_modes $out_modes_hash{ $outMode } -sm GM ";
	}
	else{
		foreach my $cytosine(@cytosines){
			$cmd .= "-C $cytosine ";
		}
		
		$cmd .= "-out_modes $out_modes_hash{ $outMode } ";
	}
	if($interval ne ""){
		$cmd .= "-L $interval ";
	}
	$cmd .= "-stand_call_conf $qual -nt $nt ";
	$cmd .= "-minConv $minConv -vcfCache 1000000 ";
	$cmd .= "-mmq $mmq ";
	$cmd .= "-mbq $mbq\n";
	print "$cmd\n";
	system($cmd);

}
 
sub vcf_sort{
	unless (-e $VCFTOOLS) {
 		print "\n $VCFTOOLS file not exists!\n\n";
 		exit(1);
	}
	my $cmd .= "perl $VCFTOOLS --k 1 --c 2 ";
	$cmd .= "--tmp $tmp_dir ";
	$cmd .= "$vcf_unsorted_cpg ";
	$cmd .= "$ref";
	$cmd .= ".fai ";
	$cmd .= "> $vcf_sorted_cpg \n";
	$cmd .= "perl $VCFTOOLS --k 1 --c 2 ";
	$cmd .= "--tmp $tmp_dir ";
	$cmd .= "$vcf_unsorted_snp ";
	$cmd .= "$ref";
	$cmd .= ".fai ";
	$cmd .= "> $vcf_sorted_snp \n";
	print "$cmd\n";
	system($cmd);
	
}

sub vcf_filter{
	my $cmd .= "$JAVA -jar $BISSNP -R $ref -T VCFpostprocess ";
	$cmd .= "-qual $qual ";
	if($interval ne ""){
		$cmd .= "-L $interval ";
	}
	if($nomeseq ne ""){
		$cmd .= "-C GCH -C HCH -C GCG -C HCG ";
	}
	else{
		foreach my $cytosine(@cytosines){
			$cytosine =~ s/\,\d+//;
			$cmd .= "-C $cytosine ";
		}
	}
	if($lowCov ne ""){
		$cmd .= "-sb 0.0 ";
	}
	$cmd .= "-oldVcf $vcf_sorted_cpg ";
	$cmd .= "-snpVcf $vcf_sorted_snp ";
	$cmd .= "-newVcf $vcf_sorted_cpg_filtered ";
	$cmd .= "-o $vcf_sorted_cpg_filtered.cpgSummary.txt -minCT $minCT\n";
	
	$cmd .= "$JAVA -jar $BISSNP -R $ref -T VCFpostprocess ";
	$cmd .= "-qual $qual ";
	if($interval ne ""){
		$cmd .= "-L $interval ";
	}
	if($nomeseq ne ""){
		$cmd .= "-C GCH -C HCH -C GCG -C HCG ";
	}
	else{
		foreach my $cytosine(@cytosines){
			$cytosine =~ s/\,\d+//;
			$cmd .= "-C $cytosine ";
		}
	}
	if($lowCov ne ""){
		$cmd .= "-sb 0.0 ";
	}
	$cmd .= "-oldVcf $vcf_sorted_snp ";
	$cmd .= "-snpVcf $vcf_sorted_snp ";
	$cmd .= "-newVcf $vcf_sorted_snp_filtered ";
	$cmd .= "-o $vcf_sorted_snp_filtered.cpgSummary.txt \n";
	print "$cmd\n";
	system($cmd);
	
	
}

sub vcf2bed{
	unless (-e $VCF2BED) {
 		print "\n $VCF2BED file not exists! \n\n";
 		exit(1);
	}
	my $cmd .= "perl $VCF2BED --only_good_call $vcf_sorted_cpg_filtered ";
	if($nomeseq ne ""){
		$cmd .= "GCH --seperate_strand \n";
		$cmd .= "perl $VCF2BED --only_good_call --seperate_strand $vcf_sorted_cpg_filtered HCG \n";
	}
	else{
		if($allC eq ""){
			$cmd .= "CG \n";
		}
		else{
			$cmd = "";
			foreach my $cytosine(@cytosines){
				$cytosine =~ s/\,\d+//;
				$cmd .= "perl $VCF2BED --only_good_call --seperate_strand $vcf_sorted_cpg_filtered $cytosine \n";
			}
			
		}
			
	}
	print "$cmd\n";
	system($cmd);
	
	unless (-e $VCF2COV) {
 		print "\n $VCF2COV file not exists! \n\n";
 		exit(1);
	}
	$cmd = "perl $VCF2COV $vcf_sorted_cpg_filtered ";
	if($nomeseq ne ""){
		$cmd .= "GCH \n";
		$cmd .= "perl $VCF2COV $vcf_sorted_cpg_filtered HCG \n";
	}
	else{
		if($allC eq ""){
			$cmd .= "CG \n";
		}
		else{
			$cmd = "";
			foreach my $cytosine(@cytosines){
				$cytosine =~ s/\,\d+//;
				$cmd .= "perl $VCF2COV $vcf_sorted_cpg_filtered $cytosine \n";
			}
			
		}
	}
	
	print "$cmd\n";
	system($cmd);
	
	unless (-e $VCF2BEDGRAPH) {
 		print "\n $VCF2BEDGRAPH file not exists! \n\n";
 		exit(1);
	}
	$cmd = "perl $VCF2BEDGRAPH $vcf_sorted_cpg_filtered ";
	if($nomeseq ne ""){
		$cmd .= "GCH \n";
		$cmd .= "perl $VCF2BEDGRAPH $vcf_sorted_cpg_filtered HCG \n";
	}
	else{
		if($allC eq ""){
			$cmd .= "CG \n";
		}
		else{
			$cmd = "";
			foreach my $cytosine(@cytosines){
				$cytosine =~ s/\,\d+//;
				$cmd .= "perl $VCF2BEDGRAPH $vcf_sorted_cpg_filtered $cytosine \n";
			}
			
		}
			
	}
	print "$cmd\n";
	system($cmd);
	
	my $bedgraph_cpg=$vcf_sorted_cpg_filtered;
	my $bedgraph_hcg=$vcf_sorted_cpg_filtered;
	
	$cmd .= "$BEDGRAPH2BIGWIG ";
	unless (-e $BEDGRAPH2BIGWIG) {
 		print "\n $BEDGRAPH2BIGWIG file not exists! \n\n";
 		exit(1);
	}
	if($nomeseq ne ""){
		$bedgraph_cpg=~s/\.vcf$/.GCH.bedgraph/;
		my $bigwig_cpg=$bedgraph_cpg;
		$bigwig_cpg=~s/\.bedgraph$/.bw/;
		$cmd .= "$bedgraph_cpg $CHROMSIZE $bigwig_cpg \n";
		$bedgraph_hcg=~s/\.vcf$/.HCG.bedgraph/;
		my $bigwig_hcg=$bedgraph_hcg;
		$bigwig_hcg=~s/\.bedgraph$/.bw/;
		$cmd .= "$BEDGRAPH2BIGWIG $bedgraph_hcg $CHROMSIZE $bigwig_hcg \n";
	}
	else{
		$bedgraph_cpg=~s/\.vcf$/.CG.bedgraph/;
		my $bigwig_cpg=$bedgraph_cpg;
		$bigwig_cpg=~s/\.bedgraph$/.bw/;
		$cmd .= "$bedgraph_cpg $CHROMSIZE $bigwig_cpg \n";
	}
	
	system($cmd);
	print "$cmd\n";
	
}


sub vcf2bigwig{
	my $vcf_file=shift @_;
	unless (-e $VCF2BEDGRAPH) {
 		print "\n $VCF2BEDGRAPH file not exists! \n\n";
 		exit(1);
	}
	my $cmd .= "perl $VCF2BEDGRAPH $vcf_file ";
	if($nomeseq ne ""){
		$cmd .= "GCH \n";
		$cmd .= "perl $VCF2BEDGRAPH $vcf_file HCG \n";
	}
	else{
		$cmd .= "CG \n";
	}
	my $bedgraph_cpg=$vcf_file;
	my $bedgraph_hcg=$vcf_file;

	$cmd .= "$BEDGRAPH2BIGWIG ";
	unless (-e $BEDGRAPH2BIGWIG) {
 		print "\n $BEDGRAPH2BIGWIG file not exists! \n\n";
 		exit(1);
	}
	if($nomeseq ne ""){
		$bedgraph_cpg=~s/\.vcf$/.GCH.bedgraph/;
		my $bigwig_cpg=$bedgraph_cpg;
		$bigwig_cpg=~s/\.bedgraph$/.bw/;
		$cmd .= "$bedgraph_cpg $CHROMSIZE $bigwig_cpg \n";
		$bedgraph_hcg=~s/\.vcf$/.HCG.bedgraph/;
		my $bigwig_hcg=$bedgraph_hcg;
		$bigwig_hcg=~s/\.bedgraph$/.bw/;
		$cmd .= "$BEDGRAPH2BIGWIG $bedgraph_hcg $CHROMSIZE $bigwig_hcg \n";
	}
	else{
		$bedgraph_cpg=~s/\.vcf$/.CG.bedgraph/;
		my $bigwig_cpg=$bedgraph_cpg;
		$bigwig_cpg=~s/\.bedgraph$/.bw/;
		$cmd .= "$bedgraph_cpg $CHROMSIZE $bigwig_cpg \n";
	}
	
	system($cmd);
	print "$cmd\n";
	
	
	
}

##-cpgreads /Volumes/HD_2/Documents/workspace/hcg.reads.txt -gchreads /Volumes/HD_2/Documents/workspace/gch.reads.txt -out_modes NOMESEQ_MODE -sm GM
##-cpgreads /Volumes/HD_2/Documents/temp/lei/cpg.reads.txt  -stand_call_conf 20
sub bissnp_epistate{
	my $cmd .= "$JAVA -jar $BISSNP -R $ref ";
	$cmd .= "-I $input_file ";									
	$cmd .= "-D $dbsnp -T BisulfiteGenotyper -vfn1 $vcf_unsorted_cpg -vfn2 $vcf_unsorted_snp ";
	if($nomeseq ne ""){
		$cmd .= "-out_modes $out_modes_hash{ $outMode } -sm GM ";
	}
	else{
		foreach my $cytosine(@cytosines){
			$cmd .= "-C $cytosine ";
		}
		
		$cmd .= "-out_modes $out_modes_hash{ $outMode } ";
	}
	if($interval ne ""){
		$cmd .= "-L $interval ";
	}
	$cmd .= "-stand_call_conf $qual -nt $nt ";
	$cmd .= "-minConv $minConv -vcfCache 1000000 ";
	$cmd .= "-mmq $mmq ";
	$cmd .= "-mbq $mbq\n";
	print "$cmd\n";
	system($cmd);

}



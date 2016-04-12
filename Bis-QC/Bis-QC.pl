#!/usr/bin/perl -w
##This script is used for quality control on WGBS/NOMe-seq.
#It accept .fastq/.fa file for the QC before reads mapping (check adapter contamination and inverted duplication reads)
#It accept .bam file for the QC after reads mapping (check bisulfite conversion rate, coverage distribution, NOMe Enzyme efficiency, NOMe Enzyme CCG leak)


##example cmd:

##mode 0 (before reads mapping): perl Bis-QC.pl [option] --QC_mode 0 input.fastq/input.fq
##mode 1 (after reads mapping): perl Bis-QC.pl [option] --QC_mode 1 input.bam


## author: Yaping Liu  lyping1986@gmail.com 
## time: 2014-7-11

#Usage: perl Bis-QC.pl [option] input.bam/input.fastq

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "QC mode 0 (before reads mapping): perl Bis-QC.pl [option] R1.fastq/R1.fq [R2.fastq/R2.fq]\n";
    print STDERR "                                  It will generate summary statistic report on adapter contamination and inverted dulication reads percentage \n";
    print STDERR "                                  You could also enable the trimming step with --adapter_trim or --invdups_trim.\n\n";
    print STDERR "QC mode 1 (after reads mapping): perl Bis-QC.pl [option] input.bam \n";
    print STDERR "                                 It will generate summary statistic report on bisulfite conversion rate(Methylation bias check or whole reads level), coverage distribution, NOMe Enzyme efficiency, NOMe Enzyme CCG leak \n\n";
    print STDERR " Pre-request:\n"; 
    print STDERR "It requires the installation of R. If you can not start R directly, you could also specify your own version of R with --r FILE\n";
    print STDERR "It requires the installation of fastq-mcf in ea-utils. We provided a pre-compiled version in External_tools directory. You could also specify your own version with --fastqmcf FILE\n";
    print STDERR "We provided a illumina adapter sequence file in External_tools directory. You could also specify your own .fasta file with --adapter FILE.\n\n";
    print STDERR " [Options]:\n\n";
    print STDERR " ###############General options:\n\n";
	print STDERR "  --QC_mode NUM : Specify QC mode. Otherwise, it will automately enable the mode by the file name's suffix (.fastq/.fq file will enable mode 0, .bam file will enable mode 1).\n";
	print STDERR "                  Mode 0: --QC_mode 0  QC before reads mapping. (check adapter contamination and inverted duplication reads)\n";
	print STDERR "                  Mode 1: --QC_mode 1  QC after reads mapping. (check bisulfite conversion rate(Five prime or whole reads level), coverage distribution, NOMe Enzyme efficiency, NOMe Enzyme CCG leak)\n";
	print STDERR "                  (Default: not specified, automately enable the mode by the file name's suffix)\n\n";
	print STDERR "  --bistools_path DIR : Specify the Bis-tools root direcotry (Default: not specified. use environment variable \$BISTOOLS).\n";
	print STDERR "  --r FILE : Specify the full path to use R (Default: not specified. use global setting of R).\n";
	print STDERR "  --nt NUM : Specify the number of cpu to use (Default: not specified. use all of the cpu in the computer).\n";
	print STDERR "  --mem NUM : Specify the number of gigabytes in the memory to use (Default: 8).\n";
		
	
    print STDERR " ############### Mode 0 options:\n\n";
	print STDERR "  --fastqmcf FILE : Specify the full path to use fastq-mcf in ea-utils (Default: not specified. use pre-compiled version in Bis-tools' External_tools directory).\n";
	print STDERR "  --adapter FILE : Specify the full path to the illumina adapter .fasta file (Default: not specified. use default version in Bis-tools' resource directory).\n";
	print STDERR "  --libcomplexity FILE : Specify the full path to use  (Default: not specified. use pre-compiled version in Bis-tools' External_tools directory).\n";
	print STDERR "  --disable_adapter : Disable the adapter contamination check step (Default: not enabled).\n\n";
	print STDERR "  --disable_invdups : Disable the inverted duplication reads check step (Default: not enabled).\n\n";
	print STDERR "  --disable_lib_complexity : Disable the library complexity check step (Default: not enabled).\n\n";
	print STDERR "  --adapter_trim NUM: If adapter contamination are more than this criteria, it will enable the adapter trimming step (Default: 0.01,   it means more than 1% adapter contamination).\n\n";
	print STDERR "  --invdups_trim NUM: If inverted duplication reads contamination are more than this criteria, it will enable the adapter trimming step (Default: 0.1,   it means more than 10% invert-dups contamination).\n\n";
    
    	
    print STDERR " ############### Mode 1 options:\n\n";	
	print STDERR "  --disable_bs_conv_check : Disable the bisulfite conversion check step (Default: not enabled).\n\n";
	print STDERR "  --disable_coverage_check : Disable the coverage distribution in the genome check step (Default: not enabled).\n\n";
	print STDERR "  --disable_enzyme_eff_check : Disable the NOMe enzyme efficiency check step (Default: not enabled).\n\n";
	print STDERR "  --disable_trinuc_check : Disable the step to count all of the trinucleotides' methylaiton level at chr21 and chrM. it could be used to estimate NOMe enzyme CCG off-target (Default: not enabled).\n\n";
	print STDERR "  --pattern STR: Cytosine pattern that used to judge the bisulfite conversion rate (Default: WCW.  W:A/T).\n\n";
	print STDERR "  --genome FILE: reference genome .fasta file that are used when you map reads in the bam file (Default: not specified. use hg19.fa in Bis-tools' resource directory).\n\n";
	print STDERR "  --dbsnp FILE: dbSNP .vcf file that are used for genotyping (Default: not specified. use dbSNP_135.hg19.sort.vcf in Bis-tools' resource directory).\n\n";
	print STDERR "  --enzyme_eff_regions FILES: Region to test enzyme efficiency. (Default: just use CTCF conserved motif and CGI promoters. allow multiple location files).\n\n";
 	
    exit(1);
}


##default option setting
my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);
my $QC_mode="";
my $r="R";
my $nt=`cat /proc/cpuinfo | grep processor -c`;
chomp($nt);
if($nt eq "" || $nt==0){
	$nt=1;
}
my $mem=8;

my $fastqmcf="$bistools_path/External_tools/ea-utils/fastq-mcf";
my $adapter="$bistools_path/resource/illumina_adapter_new.fa";
my $libcomplexity="";
my $disable_adapter="";
my $disable_invdups="";
my $adapter_trim=0.01;
my $invdups_trim=0.1;

my $disable_bs_conv_check="";
my $disable_coverage_check="";
my $disable_enzyme_eff_check="";
my $disable_trinuc_check="";

my $pattern = "WCW";
my $genome="$bistools_path/resource/genome/hg19_rCRSchrm.fa";
my $dbsnp="$bistools_path/resource/dbSNP/dbsnp_135.hg19.sort.vcf";
my @enzyme_eff_regions=();


GetOptions( 
			"bistools_path=s" => \$bistools_path,
			"QC_mode=i" => \$QC_mode,
			"r=s" => \$r,
			"nt=i" => \$nt,
			"mem=i" => \$mem,
			"fastqmcf=s" => \$fastqmcf,
			"adapter=s" => \$adapter,
			"libcomplexity=s" => \$libcomplexity,
			"disable_adapter" => \$disable_adapter,
			"disable_invdups" => \$disable_invdups,
			"adapter_trim=f" => \$adapter_trim,
			"invdups_trim=f" => \$invdups_trim,
			"disable_bs_conv_check" => \$disable_bs_conv_check,
			"disable_coverage_check" => \$disable_coverage_check,
			"disable_enzyme_eff_check" => \$disable_enzyme_eff_check,
			"disable_trinuc_check" => \$disable_trinuc_check,
			"pattern=s" => \$pattern,
			"genome=s" => \$genome,
			"dbsnp=s" => \$dbsnp,
			"enzyme_eff_regions=s" => \@enzyme_eff_regions,
			#"r=s" => \$r,
);

	

##STEP 1: check parameters
check_parameter();

##STEP 2: determine mode
if($QC_mode eq ""){
	if(scalar(@ARGV)==2){
		if($ARGV[0] =~/\.bam$/i or $ARGV[1] =~/\.bam$/i){
			print STDERR "Wrong mode number for BAM files\n\n";
			usage();
		}else{
			mode_0(@ARGV);
		}
	}else{
		if($ARGV[0] =~/\.bam$/i){
			mode_1(@ARGV);
		}elsif($ARGV[0] =~/\.fastq$/i or $ARGV[0] =~/\.fq/i or $ARGV[0] =~/\.txt/i){
			mode_0(@ARGV);
		}else{
			print STDERR "Not recoganized file format!!\n\n";
			usage();
		}
	}
}elsif($QC_mode == 0){
	mode_0(@ARGV);
}elsif($QC_mode == 1){
	mode_1(@ARGV);
}else{
	print STDERR "Wrong mode number\n\n";
	usage();
}



##########################################All Subroutines###############################################################

sub check_parameter {
	usage() if ( scalar(@ARGV) == 0 );
	if ( scalar(@ARGV) > 2 || (scalar(@ARGV)>1 && $QC_mode==1)) {
    	print STDERR "Mode 0 only allow maximum 2 fastq files. Mode 1 only allow maximum 1 bam file\n\n";
    	usage();
	}
}

sub mode_0 {
	print STDERR "Mode 0: QC on reads before mapping \n\n\n";
	##STEP 1 (Mode 0): adapter check & trimming
	if($disable_adapter eq ""){
		check_adapter(@ARGV);
	}
	
	##STEP 2 (Mode 0): inverted duplication reads check & trimming
	if($disable_invdups eq ""){
		check_inv_dups(@ARGV);
	}
	
	##STEP 3 (Mode 0): library complexity check
	##TODO: need to implement later
	#check_lib_complexity(@ARGV);
}

sub mode_1 {
	##STEP 0 (Mode 1): check bam pre-requirement:
	check_bam_preprocess(@ARGV);
	
	
	##STEP 1 (Mode 1): bisulfite conversion check
	if($disable_bs_conv_check eq ""){
		check_bs_conv(@ARGV);
	}
	
	##STEP 2 (Mode 1): coverage distribution check
	if($disable_coverage_check eq ""){
		check_cov_dist(@ARGV);
	}
	
	##STEP 3 (Mode 1): NOMe enzyme efficiency check. Plot average NOMe-seq signal around CTCF and CGI_TSS
	if($disable_enzyme_eff_check eq ""){
		check_enzyme_eff(@ARGV);
	}
	
	##STEP 4 (Mode 1): NOMe enzyme CCG off-target check. check trinucleotide methylation level in chr21 and chrM
	if($disable_trinuc_check eq ""){
		check_trinuc(@ARGV);
	}
	
}



##########################################Mode 0###############################################################

sub check_adapter {
	my @reads=@_;
	my $cmd="perl $bistools_path/Bis-QC/before_reads_mapping/wrap_fastq-mcf.pl --fastqmcf $fastqmcf --adapter $adapter --nt $nt --adapter_trim $adapter_trim $reads[0] ";
	if(scalar(@reads)>1){
		$cmd .= "$reads[1]\n";
	}
	$cmd .= "\n";
	print STDERR "check adapter contamination:\n$cmd\n";
	system($cmd)==0 || die "Unexpected stop at adapter contamination step: $! \n";
}


sub check_inv_dups {
	my @reads=@_;
	my $log=$reads[0].".inv-dups.log.txt";
	my $cmd="perl $bistools_path/Bis-QC/before_reads_mapping/invert_dups_check.pl $log @reads \n";
	print STDERR "check inverted duplication reads:\n$cmd\n";	
	system($cmd)==0 || die "Unexpected stop at inverted duplication reads check step: $! \n";
	
	open(IN,"<$log") or die "can't find inverted duplication reads log file:$!\n";
	my $inv_dups_perc=0;
	while(<IN>){
		if($_=~/inverted Pair Percentage=(\S+)/){
			$inv_dups_perc=$1;
		}
	}
	close(IN);
	if($inv_dups_perc > $invdups_trim*100){
		print STDERR "There are $inv_dups_perc % Inverted duplication reads, which is more than the trimming threshold.\n\n Trimming inverted duplication reads:\n";
		trim_inv_dups(@reads);
	}
}

sub trim_inv_dups {
	my @reads=@_;
	my $log=$reads[0].".inv-dups.log.txt";
	my $cmd="perl $bistools_path/Bis-QC/before_reads_mapping/invert_dups_check.pl --trim_invert_dups --output_clean $log @reads \n";
	print STDERR "$cmd\n";
	system($cmd)==0 || die "Unexpected stop at inverted duplication reads trimming step: $! \n";
}

sub check_lib_complexity {
	
}

##########################################Mode 1###############################################################

sub check_bam_preprocess {
	##check if bam file is sorted or indexed
	
	##check if bam file has read group tag and ask if user would like to add the default read group tag.
	
	##check if reference genome order and bam file header are the same, ask if they would like to reorder the bam file based on the reference genome provided.
	
}

sub check_bs_conv {
	my $bam=shift @_;
	#### Methylation bias check, bisulfite conversion along cycle, check 5' conversion rate
	my $cmd="perl $bistools_path/Bis-QC/after_reads_mapping/methylation_bias_plot.pl --bissnp $bistools_path/Bis-SNP/Bis-SNP.latest.jar --r $r --pattern $pattern --genome $genome $bam\n";
	print STDERR "Methylation bias check:\n $cmd\n";
	system($cmd)==0 || die "Unexpected stop at methylation bias check step: $! \n";
	
	##bisulfite conversion rate for the whole reads
	$cmd="perl $bistools_path/Bis-QC/after_reads_mapping/bisulfite_conv_distribution.pl --bissnp $bistools_path/Bis-SNP/Bis-SNP.latest.jar --r $r --pattern $pattern --genome $genome $bam\n";
	print STDERR "Bisulfite conversion rate distribution of whole reads:\n $cmd\n";
	system($cmd)==0 || die "Unexpected stop at bisulfite conversion rate distribution check step: $! \n";
	
	##Methylation level around different trinucleotide on chrM and chr21
	if($disable_trinuc_check eq ""){
		my $tri_nuc_log=$bam;
		$tri_nuc_log=~s/\.bam$/.trinuc_methy.chrM.txt/;
		$cmd="perl $bistools_path/Bis-QC/after_reads_mapping/bissnp_trinuc_sample.pl --bissnp $bistools_path/Bis-SNP/Bis-SNP.latest.jar --genome $genome --dbsnp $dbsnp --nt $nt --mem $mem --interval chrM $tri_nuc_log $bam\n";
		print STDERR "Bisulfite conversion rate distribution of whole reads at chrM:\n $cmd\n";
		system($cmd)==0 || die "Unexpected stop at methylation level of trinucleotides in chrM check step: $! \n";
	
		$tri_nuc_log=$bam;
		$tri_nuc_log=~s/\.bam$/.trinuc_methy.chr21.txt/;
		$cmd="perl $bistools_path/Bis-QC/after_reads_mapping/bissnp_trinuc_sample.pl --bissnp $bistools_path/Bis-SNP/Bis-SNP.latest.jar --genome $genome --dbsnp $dbsnp --nt $nt --mem $mem --interval chr21 $tri_nuc_log $bam\n";
		print STDERR "Bisulfite conversion rate distribution of whole reads at chr21:\n $cmd\n";
		system($cmd)==0 || die "Unexpected stop at methylation level of trinucleotides in chr21 check step: $! \n";
	}

	
}


sub check_cov_dist {
	##check coverage distribution
	my $bam=shift @_;
	my $prefix=basename($bam);
	my $dirname=dirname($bam);
	$prefix=~s/(\w+)\S+$/Coverage_distribution_$1/;
	my $cmd="perl $bistools_path/Bis-QC/after_reads_mapping/coverageEstimate.pl --covs 1 --covs 3 --covs 5 --covs 7 --covs 10 --covs 15 --covs 20 --covs 30 --covs 60 --mem $mem --cpu $nt --ref $genome --bissnp $bistools_path/Bis-SNP/Bis-SNP.latest.jar --result_dir $dirname  $prefix $bam\n";
	print STDERR "Coverage distribution check:\n $cmd\n";
	system($cmd)==0 || die "Unexpected stop at coverage distribution check step: $! \n";
}

sub check_enzyme_eff {
	##check NOMe Enyzme efficiency methylation/accessibility level around conserved CTCF motif and CGI promoters
	my $bam=shift @_;
	#--r_script /home/uec-00/yapingli/code/mytools/R/MethyPatternFeaturePlotSinglePlotNoSumFeature.R --sort_perl_script /home/uec-00/yapingli/code/mytools/perl/sortByRefAndCor.pl --pbs --result_dir /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/MAR_MPR/
	my $prefix=basename($bam);
	my $dirname=dirname($bam);
	if(scalar(@enzyme_eff_regions)==0){
		$prefix=~s/(\w+)\S+$/Conserved_CTCF_$1/;
		my $cmd="perl $bistools_path/Bis-QC/after_reads_mapping/MethyPatternAlignEasyUsage.pl --mem $mem --cpu $nt --nomeseq $bistools_path/Bis-SNP/Bis-SNP.latest.jar $bam $bistools_path/resource/ctcf/wgEncodeUwTfbsAllCtcf.hg19.uniq.no_ensembl_r75_tss.chr1.sort.bed $genome $prefix $dbsnp ";
		$cmd.="--r_script $bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotForNOMeSeq.R --sort_perl_script $bistools_path/utils/sortByRefAndCor.pl --result_dir $dirname\n";
		print STDERR "Check NOMe Enyzme efficiency methylation/accessibility level around conserved CTCF motif:\n $cmd\n";
		print STDERR "Default refion is on hg19 !!!\n";
		system($cmd)==0 || die "Unexpected stop at NOMe Enyzme efficiency check step: $! \n";
	
		$prefix=basename($bam);
		$prefix=~s/(\w+)\S+$/CGI_TSS_$1/;
		$cmd="perl $bistools_path/Bis-QC/after_reads_mapping/MethyPatternAlignEasyUsage.pl --mem $mem --cpu $nt --nomeseq $bistools_path/Bis-SNP/Bis-SNP.latest.jar $bam $bistools_path/resource/tss/knownGene-tss-ucsc08082013-unique.tj_gg_plus200bp_cgi.hg19.noChrM.sort.bed $genome $prefix $dbsnp ";
		$cmd.="--r_script $bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotForNOMeSeq.R --sort_perl_script $bistools_path/utils/sortByRefAndCor.pl --result_dir $dirname\n";
		print STDERR "Check NOMe Enyzme efficiency methylation/accessibility level around CGI promoters:\n $cmd\n";
		print STDERR "Default refion is on hg19 !!!\n";
		system($cmd)==0 || die "Unexpected stop at NOMe Enyzme efficiency check step: $! \n";
	}else{
		foreach my $enzyme_eff_region(@enzyme_eff_regions){
			$prefix=~s/(\w+)\S+$/$1/;
			my $loc_prefix=basename($enzyme_eff_region);
			$loc_prefix=~s/(\S+)\.\S+$/$1/;
			my $cmd="perl $bistools_path/Bis-QC/after_reads_mapping/MethyPatternAlignEasyUsage.pl --mem $mem --cpu $nt --nomeseq $bistools_path/Bis-SNP/Bis-SNP.latest.jar $bam $enzyme_eff_region $genome $prefix $dbsnp ";
			$cmd.="--r_script $bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotForNOMeSeq.R --sort_perl_script $bistools_path/utils/sortByRefAndCor.pl --result_dir $dirname\n";
			print STDERR "Check NOMe Enyzme efficiency methylation/accessibility level around $enzyme_eff_region:\n $cmd\n";

			system($cmd)==0 || die "Unexpected stop at NOMe Enyzme efficiency check step: $! \n";
			
		}
	}
	
}


sub check_bam_sort{
	
}

sub check_bam_index{
	
}

sub check_bam_read_group{
	
}

sub check_bam_header{
	
}






#!/usr/bin/perl -w
##This script wrap BisSNPCoverageWalker to estimate the coverage in the whole genome or in different genomic regions. 
##It enable to estimate on raw coverage(filter duplicated and non-unique aligned read) or good reads coverage (used by BisSNP)
##also, it draw .pdf for the accumulative distribution of the coverage



## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-10-30

#Usage: perl coverageEstimate.pl [option] prefix input.bam

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl coverageEstimate.pl [option] prefix input1.bam input2.bam ... \n";
    print STDERR "Wrap BisSNPCoverageWalker to estimate the coverage in the whole genome or in different genomic regions.\n";
    print STDERR "It enable to estimate on raw coverage(filter duplicated and non-unique aligned read) or good reads coverage (used by BisSNP).\n";
    print STDERR "Also, it draw .pdf for the accumulative distribution of the coverage.\n\n";
    print STDERR "  [Options]:\n\n";
    print STDERR "  --omit_coverage: omit coverage calculation step. (Default: not enabled)\n\n";
    print STDERR "  --result_dir STR: directory to put result. (Default: current directory)\n\n";
    print STDERR "  --ref STR: reference genome .fasta file. (Default: use hg19 in hpcc)\n\n";
	print STDERR "  --good_reads : only use proper paired, mismatch low reads and filter out incompletly converted reads. (Default: not enabled)\n\n";
	print STDERR "  --intervals STR : The interval regions to summarize coverage distribution. Allow multiple input. (Default: use default set of interval file in hpcc)\n\n";
	print STDERR "  --whole_genome : It will ignore the --interval, only do summarise information in the whole genome. (Default: not enabled, use default set of interval file in hpcc)\n\n";
	print STDERR "  --covs INT: Coverage cut-off. Allow multiple inputs, it will summarize the proportion of bases over these coverage. (Default: --cov 1 --cov 3 --cov 5 --cov 7 --cov 10 --cov 15 --cov 20 --cov 25 --cov 30 --cov 40 --cov 50 --cov 100)\n\n";
    print STDERR "  --min INT: the minimum accumulative coverage in the plot. (Default: --min 1)\n\n";
    print STDERR "  --max INT: the maximum accumulative coverage in the plot. (Default: --max 20)\n\n";
    print STDERR "  --bin INT: the bin size of the accumulative coverage in the plot. (Default: --bin 5)\n\n";
    print STDERR "  --mmq INT: minimum mapping quality score. (Default: --mmq 1)\n\n";
    print STDERR "  --mbq INT: minimum base quality score. (Default: --mbq 0)\n\n";
    print STDERR "  --split_readgroup: split result by read group. (Default: not enabled)\n\n";
    print STDERR "  --cpu INT: number of cpu to use. (Default: --cpu 1)\n\n";
    print STDERR "  --mem INT: number of Gigabyte memory to use. (Default: --mem 4)\n\n";
    print STDERR "  --bissnp STR: the bissnp .jar file location. (Default: use bissnp-default.jar in the hpcc)\n\n";
    print STDERR "  --r_script STR: the r script to plot accumulative coverage distribution. (Default: use r_script in the hpcc)\n\n";
    exit(1);
}
my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);

my $omit_coverage = "";
my $result_dir = `pwd`;
chomp($result_dir);
$result_dir.="/";

my $ref = "$bistools_path/resource/genome/hg19_rCRSchrm.fa";
my $good_reads = "";
my @intervals = ();
my $whole_genome = "";
my @covs = ();
my $min=1;
my $max=20;
my $bin=19;
my $mmq = 1;
my $mbq = 0;
my $split_readgroup="";
my $r_script = "$bistools_path/Bis-QC/after_reads_mapping/cumCovDistPlot.R";
my $bissnp = "$bistools_path/Bis-SNP/Bis-SNP.latest.jar";
my $mem=4;
my $cpu=1;
my $R = "R";

GetOptions( "omit_coverage" => \$omit_coverage,
			"result_dir=s" => \$result_dir,
			"ref=s" => \$ref,
			"good_reads" => \$good_reads,
			"intervals=s" => \@intervals,
			"whole_genome" => \$whole_genome,
			"covs=i" => \@covs,
			"min=i" => \$min,
			"max=i" => \$max,
			"bin=i" => \$bin,
			"mmq=i" => \$mmq,
			"mbq=i" => \$mbq,
			"split_readgroup" => \$split_readgroup,
			"r_script=s" => \$r_script,
			"r=s" => \$R,
			"bissnp=s" => \$bissnp,
			"cpu=i" => \$cpu,
			"mem=i" => \$mem,
);


usage() if ( scalar(@ARGV) < 2 );

my $prefix=shift(@ARGV);
my @bams=@ARGV;

if(scalar(@covs)==0){
	@covs = (1,3,5,7,10,15,20,25,30,40,50,70,100);
}

if(scalar(@intervals)==0 && $whole_genome ne ""){
	@intervals=("$bistools_path/resource/interval/hg19/whole_genome_interval_list.hg19.bed");
}elsif(scalar(@intervals)==0 && $whole_genome eq ""){
	@intervals=("$bistools_path/resource/interval/hg19/knownGene-tss-ucsc08082013-unique.hg19.noChrM.2kbUpDown.sort.bed",
				"$bistools_path/resource/interval/hg19/knownGene-tss-ucsc08082013-unique.tj_gg_plus200bp_cgi.hg19.noChrM.2kbUpDown.sort.bed",
				"$bistools_path/resource/interval/hg19/knownGene-tss-ucsc08082013-unique.non_tj_gg_irizarry2009_cgi_500bp.hg19.noChrM.2kbUpDown.sort.bed",
				"$bistools_path/resource/interval/hg19/knownGene.genePredToGtf.ExonOnly.uniq.hg19.noChrM.noKnownTss2kb.bed",
				"$bistools_path/resource/interval/hg19/knownGene-intron-ucsc08082013.genePredToGtf.hg19.noChrM.noKnownTss2kb.bed",
				"$bistools_path/resource/interval/hg19/knownGene-intergenic-ucsc08082013.2kbUpstreamFromGenes.uniq.hg19.noChrM.bed",
				"$bistools_path/resource/interval/hg19/Takai_Jones_plus_GG.merged.hg19.bed",
				"$bistools_path/resource/interval/hg19/xie_cuddapeh_plus_kim.orientedOnly.noKnownTss4kb.hg19.1kbUpDown.bed",
				"$bistools_path/resource/interval/hg19/UCSC.wgEncodeBroadHmmGm12878HMM.SE.hg19.bed",						
				"$bistools_path/resource/interval/hg19/wgEncodeUwDnaseGm12878HotspotsRep1.hg19.bed");
}

foreach my $interval(@intervals){
	my $java_cmd = "java -Xmx${mem}G -jar $bissnp -T BisSNPCoverage -R $ref -omitIntervals -omitBaseOutput -L $interval -nt $cpu ";
	foreach my $bam(@bams){
		$java_cmd.="-I $bam ";
	}
	if($good_reads eq ""){ ##use all of reads except duplicated and non-unique aligned reads
		$java_cmd.="-mm 1.1 -minPatConv 1.1 -invDups USE_BOTH_END -badMate ";
	}else{
		$java_cmd.="-invDups NOT_TO_USE ";
	}
	foreach my $cov(@covs){
		$java_cmd.="-ct $cov ";
	}
	if($split_readgroup ne ""){
		$java_cmd.="-pt readgroup ";
	}
	
	$interval=basename($interval);
	$interval=~s/\.\w+$//;
	my $prefix_tmp=$prefix.".".$interval;
	$interval="$result_dir".$prefix_tmp;
	my $nbins=$max-$min;
	$java_cmd.="-o $interval --start $min --stop $max --nBins $nbins\n";
	if($omit_coverage eq ""){
		print STDERR $java_cmd;
		system($java_cmd)==0 || die "can't finish coverage esitmation:$!\n";
	}
	
	
	my $input = $interval.".sample_cumulative_coverage_proportions";
	my $r_cmd="$R --no-restore --no-save --args wd=$result_dir input=$input prefix=$prefix_tmp xmin=$min xmax=$max step=$bin < $r_script\n";
	print STDERR $r_cmd;
	system($r_cmd)==0 || die "can't finish coverage esitmation plot:$!\n";
}




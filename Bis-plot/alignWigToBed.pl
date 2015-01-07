#!/usr/bin/perl -w
## align big wig file in to genome feature bed file
## author: Yaping Liu  lyping1986@gmail.com
##example script:
#####average plot (--loc > 1, it will overlap multiple line in the same pdf):


#####density bar plot:
##perl /home/rcf-40/yapingli/storage/code/mytools/perl/alignWigToBed.pl --result_dir /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/ --locs


#####heatmap plot:
##pure NOME-seq
##perl /home/rcf-40/yapingli/storage/code/mytools/perl/alignWigToBed.pl --result_dir /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/ /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/tss_all/Cluster4.heatmap.NONCGI_TSS.knownGene-tss-UCSC062409-unique.non_irizarry2009_cgi.DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.alignedTo.knownGene.1200.GCH.txt.bed testHeatmap /export/uec-gs1/laird/users/yaping/data/Fides/bam/HCT116/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.hg19_rCRSchrm.fa.mdups.realign.recal.clean.cytosine.filtered.sort.BisSNP-0.81.GCH.bw /export/uec-gs1/laird/users/yaping/data/Fides/bam/HCT116/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.hg19_rCRSchrm.fa.mdups.realign.recal.clean.cytosine.filtered.sort.BisSNP-0.81.HCG.bw /export/uec-gs1/laird/users/yaping/data/Fides/bam/HCT116/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/HCT116_merge_C118BACXX_C18FYACXX_D1W2BACXX.fastq-mcf.hg19.mdups.header.clean.realign.recal.cytosine.filtered.sort.BisSNP-0.81.1.GCH.bw /export/uec-gs1/laird/users/yaping/data/Fides/bam/HCT116/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/HCT116_merge_C118BACXX_C18FYACXX_D1W2BACXX.fastq-mcf.hg19.mdups.header.clean.realign.recal.cytosine.filtered.sort.BisSNP-0.81.1.HCG.bw --nomeseq --heatmap --omit_align_step --heatmap_anno /export/uec-gs1/laird/users/yaping/code/NOMeseq/NDR/HMM/HCT116_hg19/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/methy_plot/DKO1_NDR.fdr0.01.MARsMergedIn80bp.sort.bed --heatmap_anno /export/uec-gs1/laird/users/yaping/code/NOMeseq/NDR/HMM/HCT116_hg19/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/methy_plot/HCT116_NDR.fdr0.01.MARsMergedIn80bp.sort.bed --heatmap_anno_col blue --heatmap_anno_col red --heatmap_anno_name Cluster --heatmap_anno_name DKO1_NDR --heatmap_anno_name HCT116_NDR
#NOME-seq as clustering order, then plot chip-seq heatmap
##perl /home/rcf-40/yapingli/storage/code/mytools/perl/alignWigToBed.pl --result_dir /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/ /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/tss_all/Cluster4.heatmap.NONCGI_TSS.knownGene-tss-UCSC062409-unique.non_irizarry2009_cgi.DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.alignedTo.knownGene.1200.GCH.txt.bed testHeatmap /export/uec-gs1/laird/users/yaping/data/Fides/bam/HCT116/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.hg19_rCRSchrm.fa.mdups.realign.recal.clean.cytosine.filtered.sort.BisSNP-0.81.GCH.bw /export/uec-gs1/laird/users/yaping/data/Fides/chip_seq/wiggler_result/DKO1_H2AZ_D1ALAACXX_5_KEL656A210.bw /export/uec-gs1/laird/users/yaping/data/Fides/chip_seq/wiggler_result/DKO1_Input_D1GW4ACXX_5_KEL656A225.bw  --auto_scale --data_matrix_scale 5200 --plot_x_axis_scale 5000 --plot_x_axis_step 2500 --bin_size_align 20 --colors black --colors red --colors orange  --lengends DKO1_GCH --lengends DKO1_KEL656A239 --lengends DKO1_KEL656A210 --heatmap --heatmap_anno /export/uec-gs1/laird/users/yaping/code/NOMeseq/NDR/HMM/HCT116_hg19/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/methy_plot/DKO1_NDR.fdr0.01.MARsMergedIn80bp.sort.bed --heatmap_anno /export/uec-gs1/laird/users/yaping/code/NOMeseq/NDR/HMM/HCT116_hg19/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/methy_plot/HCT116_NDR.fdr0.01.MARsMergedIn80bp.sort.bed --heatmap_anno_col blue --heatmap_anno_col red --heatmap_anno_name Cluster --heatmap_anno_name DKO1_NDR --heatmap_anno_name HCT116_NDR  --omit_align_step --heatmap_capLimit F --heatmap_capLimit T --heatmap_capLimit T --heatmap_keys T --heatmap_keys T --heatmap_keys F --heatmap_ylab Accessibility --heatmap_ylab Tag_density --heatmap_ylab Tag_density --heatmap_col white2darkgreen --heatmap_col white2red --heatmap_col white2red  
#NOME-seq as clustering order, then plot motif frequency heatmap
##perl /home/rcf-40/yapingli/storage/code/mytools/perl/alignWigToBed.pl --result_dir /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/ /export/uec-gs1/laird/users/yaping/code/NOMeseq/MethyPatternFeatureWalker/HCT116_hg19/MethyPatternResult/tss_all/Cluster4.heatmap.NONCGI_TSS.knownGene-tss-UCSC062409-unique.non_irizarry2009_cgi.DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.alignedTo.knownGene.1200.GCH.txt.bed testHeatmap /export/uec-gs1/laird/users/yaping/data/Fides/bam/HCT116/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.hg19_rCRSchrm.fa.mdups.realign.recal.clean.cytosine.filtered.sort.BisSNP-0.81.GCH.bw /export/uec-gs1/laird/users/yaping/data/Fides/chip_seq/wiggler_result/DKO1_H2AZ_D1ALAACXX_5_KEL656A210.bw /export/uec-gs1/laird/users/yaping/code/NOMeseq/motifFreq/CG.frequency.hg19_rCRSchrm.fa.bw  --auto_scale --data_matrix_scale 5200 --plot_x_axis_scale 5000 --plot_x_axis_step 2500 --bin_size_align 20 --colors black --colors red --colors orange  --lengends DKO1_GCH --lengends DKO1_KEL656A239 --lengends CG_freq --heatmap --heatmap_anno /export/uec-gs1/laird/users/yaping/code/NOMeseq/NDR/HMM/HCT116_hg19/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/methy_plot/DKO1_NDR.fdr0.01.MARsMergedIn80bp.sort.bed --heatmap_anno /export/uec-gs1/laird/users/yaping/code/NOMeseq/NDR/HMM/HCT116_hg19/merge_C118BACXX_C18FYACXX_C18G4ACXX_C1F6FACXX/methy_plot/HCT116_NDR.fdr0.01.MARsMergedIn80bp.sort.bed --heatmap_anno_col blue --heatmap_anno_col red --heatmap_anno_name Cluster --heatmap_anno_name DKO1_NDR --heatmap_anno_name HCT116_NDR --heatmap_capLimit F --heatmap_capLimit T --heatmap_capLimit F --heatmap_keys T --heatmap_keys T --heatmap_keys T --heatmap_ylab Accessibility --heatmap_ylab Tag_density --heatmap_ylab CG_density --heatmap_col white2darkgreen --heatmap_col white2red --heatmap_col white2red --motif_freq 0 --motif_freq 0 --motif_freq 1 --heatmap_autoscale F --heatmap_autoscale T --heatmap_autoscale T --omit_align_step
##pure ChIP-seq. clustering/sorting by some of the files


use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl alignWigToBed.pl [Options] prefix input1.bw input2.bw ...\n\n";

    print STDERR " align big wig file in to genome feature bed file\n";
   
    print STDERR " All the scripts should be put under the same directory. It requires R installed in your machine\n";
	
	print STDERR "  [Options]:\n\n";
 	print STDERR " #####General option:\n\n";
 	
 	print STDERR "  --locs FILES: the genomic feature's location files, which the big wig files will be aligned to. It should be bed/gtf/gff/bedgraph format. It allows multiple input.(Required!)\n\n";
 	print STDERR "  --prefixs STR: the prefix name used for each of the --locs files. It allows multiple input.(Required!)\n\n";
 	print STDERR "  --coverages FILES: the coverage file, once it provided for alignment mode 4, it will adjusted value by coverage. use total_num_C/total_num_CT for mean value in that region. (should be the same order as input big wig files)\n";
 	print STDERR "                     the output line will be methy_perc, num_CT, num_cytosine_pattern \n\n";
 	 		
 	print STDERR " #####Alignment options:\n\n";
	print STDERR "  --nomeseq : when specified, it will automately generate legend and use default colors, GCH should be input first, HCG later\n\n";
	print STDERR "  --low_coverage : when specified, it will eable the low coverage mode. input1.bw will be methy level, input2.bw will be number of total reads. then for each bin, it will pool all of C reads & total reads to calculate methy level\n\n";
	print STDERR "  --standard_bed : when specified, it will output standard bed format with name and score column\n\n";
	
	print STDERR "  --include_no_data_line : when specified, it will output line with no data point completely\n\n";
	
	print STDERR "  --adjust_center NUM: when specified, it will automately using 40bp sliding window to detect the local lowest point to realign the value to the new lowest point as a center. NUM is the value of region range to detect lowest point\n\n";
	print STDERR "  --adjust_to_highest: when specified, it will detect the local highest point rather than the lowest point \n\n";
	print STDERR "  --adjust_center_window NUM:  the sliding window size to detect lowest point nearyby. (Default: 100)\n\n";
	print STDERR "  --location_after_adjust_center FILE: when specified, it will output the new location after adjustment\n\n";
	print STDERR "  --motif_freq NUM: when specified to be 1, it will just calculate the motif frequency in the provided .bw file.could be mutliple, but need to be the same order as input.bw. e.g. count all element not NA as 1, NA as 0\n\n";
	print STDERR "  --omit_align_step : when specified, it will not do alignment but just do plot step\n\n";
	print STDERR "  --omit_plot_step : when specified and in align to bam files, it will not do R plot step but will submit two jobs at the same time for align GCH and HCG to accelerate the speed.\n\n";
	print STDERR "  --auto_scale : automately scale the axis by the maximum in this data. (Default: not enabled)\n\n";
	print STDERR "  --mask_matrix FILE :  the matrix file used to mask the value in the alignment result. All value inside mask_matrix will be masked as NA\n\n";
    print STDERR "  --r_script FILE :  R script's path and name, which is used for plot\n\n";
    print STDERR "  --ucsc_script FILE :  UCSC script's path and name, which is used for align wig file (Default: /home/uec-00/shared/production/software/UCSC_Browser_Tools/default/bigWigSummary)\n\n";
    print STDERR "  --sort_perl_script FILE :  sortByRefAndCor.pl script's path and name, which is used for sorting input feature bed file\n\n";
    print STDERR "  --result_dir DIR :   directory to store alignment matrix results and plots\n\n";
    print STDERR "  --smooth : enable smooth around bin_size region, e.g. each data point represent the avergae of methylation value around 20bp\n\n";
 	print STDERR "  --use_COMPARE : when specified, it will use COMPARE to correct the bias of average\n\n";
 
 	print STDERR " #####Average plot options:\n\n";
 	   
    print STDERR "  --data_matrix_scale NUM :   the maximum distance to the feature center/alignment start loci. e.g. 2000 means +/- 2kb around\n";
    print STDERR "                              (Default: plot_x_axis_scale + 10*bin_size)\n\n";
    print STDERR "  --bin_size NUM : bin size(bp) to average for each data point. (Default: 20)\n\n";
    print STDERR "  --bin_size_align NUM : bin size(bp) to average when doing alignment. (Default: 1)\n\n";
    print STDERR "  --plot_x_axis_scale NUM : maximum x axis distance to plot in the summary plot. (Default: 1000)\n\n";
 	print STDERR "  --plot_x_axis_step NUM : step size of x axis to plot in the summary plot (Default: 250)\n\n";
 	print STDERR "  --plot_y_axis_max NUM : maximum y axis value to plot in the summary plot (Default: 1.0)\n\n";
 	print STDERR "  --plot_y_axis_min NUM : maximum y axis value to plot in the summary plot (Default: 0.0)\n\n";
 	print STDERR "  --plot_y_axis_step NUM : step size of x axis to plot in the summary plot (Default: 0.2)\n\n";
 	print STDERR "  --colors STRING : color of line drawing in the summary plot, could be multiple (Default: black)\n\n";
 	print STDERR "  --lengends STRING : legend to put into the plot, could be multiple, but need to be the same as input big wig order (Default: no)\n\n";
 	print STDERR "  --lty NUM : type of line drawing in the summary plot, could be multiple, it is the same definition as in lty in R (Default: 1)\n\n";
 	print STDERR "  --alignment_mode NUM : (Default: --alignment_mode 1)\n";
 	print STDERR "                         1) Align to the feature center\n";
	print STDERR "                         2) Align to the feature 5' boundary\n";
	print STDERR "                         3) Align to the feature 3' boundary\n";
	print STDERR "                         4) give the summary on the bed/gtf/gff file given region\n";
	print STDERR "                         5) give the closest single value in the data_matrix defined region\n";
	print STDERR "                         6) scale the data on the bed/gtf/gff file given region (like avgerage plot across gene body), it will seperate the region into data_matrix_scale/bin_size_align+1 data point and then extend to +/- data_matrix_scale/bin_size_align of this region\n\n";
				
	print STDERR " #####Heatmap options:\n\n";
	print STDERR "  --heatmap : when specified, it will automately generate heatmap, the order will be the same as the first input.bw\n";
	print STDERR "  --heatmap_anno FILE : provide the annotation file for heatmap side bar. file names could be mutliple, but need to be standard bed format\n";
	print STDERR "  --heatmap_anno_col STRING : provide the annotation color for heatmap side bar. color names could be mutliple, but need to be the same order as heatmap_anno. it is binary, when it is not overlap with feature, it will be white\n";
	print STDERR "  --heatmap_anno_name STRING : provide the annotation name for heatmap side bar. names could be mutliple, but need to be the same order as heatmap_anno. \n";
	print STDERR "  --heatmap_col STRING : provide the heatmap color scheme, e.g. --heatmap_col blue2yellow for DNA methylation or --heatmap_col white2green for accessibility, or white2grey2cyan, it could be multiple, but order should be the same as input.bw file order\n";
	print STDERR "  --heatmap_cluster NUM : how many clusters to display and draw average plot under heatmap bottom. could be multiple, but should be the same number as location files number.. (Default: 4)\n\n";
	print STDERR "  --heatmap_png : if specified, it will output png format instead of pdf format (Default: PDF format for heatmap)\n\n";
	print STDERR "  --preClustering : if specified, it will use both of samples to do clustering, if --twoStepClustering is not \n\n";
	print STDERR "  --twoStepClustering : if specified, it will use DNA methylation data and Accessibility data in both of samples to do clustering\n\n";
	print STDERR "  --multiSampleClustering NUM: if specified, it will use multiple samples to do clustering. It allow multiple input. --multiSampleClustering 1 --multiSampleClusterin 3 (Default: 1)\n\n";
	
	print STDERR " #####Heat density bar plot options:\n\n";
	print STDERR "  --density_bar: when specified, it will automately generate heat density bar plotinput1.bw would become sample1.filelist.txt,sample2.filelist.txt.. each file list contain all the experiment in each sample, the location to store aligned matrix, the file data type(percentage, enrichment)\n\n";
	print STDERR "  --category_names STR: names represent each location file, could be multiple\n\n";
	print STDERR "  --sample_names STR: names represent each sample(file list), could be multiple\n\n";
	print STDERR "  --experiment_names STR: names represent each experiment name in the file list, could be multiple\n\n";
	print STDERR "  --rep_num_experiments NUM: number of replicates in each experiment in the file list, could be multiple\n\n";
	print STDERR "  --enrich_max NUM: the maximum number for enrichment scale. (Default: 4)\n\n";
	print STDERR "  --logscale: when specified, it will apply log2 scale for enrichment type experiment (e.g. chip-seq). (Default: Not enabled)\n\n";
	print STDERR "  --capLimit NUM: when specified, it will capLimit the up NUM part of the data to avoid the outlier for enrichment type experiment (e.g. chip-seq). (Default: -1)\n\n";
	print STDERR "  --median : when specified, it will use median rather than mean value for enrichment type experiment (e.g. chip-seq). (Default: not enabled)\n\n";
	
	
	
    exit(1);
}

##default option setting
my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);


my @locs=();
my @prefixs=();
my @coverages=();
my $omit_align_step = "";
my $omit_plot_step = "";
my $nomeseq = "";
my $low_coverage="";
my $include_no_data_line="";
my $adjust_center="";
my $adjust_to_highest="";
my $adjust_center_window=100;
my $location_after_adjust_center = "";
my @motif_freq = ();
my $auto_scale = "";
my $mask_matrix = "";
my $r_script = "$bistools_path/Bis-plot/src/R/MultipleWigPatternOverBedPlot.R";
my $heatmap_r_script = "$bistools_path/Bis-plot/src/R/calcHeatmapForAlignWig2Bed.R";
my $heatmap_with_reps_r_script ="$bistools_path/Bis-plot/src/R/MultipleWigToHeatmapPlusAve.R";
#my $heatmap_with_reps_r_script ="$bistools_path/Bis-plot/src/R/MultipleWigToHeatmapPlusAve.PaperFigure1.R";
my $densitybar_r_script = "$bistools_path/Bis-plot/src/R/MultipleWigToDensityBar.2.R";
my $ucsc_script = "$bistools_path/External_tools/ucsc_tools/bigWigSummary";
my $sort_perl_script = "$bistools_path/utils/sortByRefAndCor.pl";
my $result_dir = `pwd`;
chomp($result_dir);
$result_dir.="/";
my $smooth = 0;

my $average="";
my $bin_size = 20;
my $bin_size_align = 1;
my $plot_x_axis_scale = 1000;
my $plot_x_axis_step = 250;
my $plot_y_axis_max = 100;
my $plot_y_axis_min = 0;
my $plot_y_axis_step = 20;
my @colors = ();
my @lty = ();
my $alignment_mode = 1;
my $data_matrix_scale = $plot_x_axis_scale + 10*$bin_size;
#my $data_matrix_scale = 2000;
my @lengends=();
 

my $heatmap="";
my $heatmap_with_reps="";
my $addAverage="";
my @heatmap_anno=();
my @heatmap_anno_col=();
my @heatmap_anno_name=();
my @heatmap_col=();
my @heatmap_cluster=();
my @heatmap_keys=();
my @heatmap_ylab=();
my @heatmap_ymin=();
my @heatmap_ymax=();
my @heatmap_logscale=();
my @heatmap_capLimit=();
my @heatmap_capUpLimit="";
my @heatmap_capDownLimit="";
my @heatmap_autoscale=();
my @heatmap_file_num_to_cluster=();
my $heatmap_output_subCluster="";
my $heatmap_regionToCluster_low=""; ## need to change to a more easier understood way..
my $heatmap_regionToCluster_high="";
my $fileNumToPrintSideBar=1;
my $heatmap_png="";
#preClustering <- ifelse(twoStepClustering,TRUE,preClustering)
my $preClustering="";
my $twoStepClustering="";
my $heatmap_row_order="";
my $heatmap_normalization_by_mean="";
my $heatmap_max_occ_order=""; ##use owler paper's way to order

my $heatmap_clustering_scale = 1000;
my $heatmap_clustering_bin_size = 20;
my $heatmap_clustering_bin_size_align = 1;
my $heatmap_clustering_matrix = $heatmap_clustering_scale + 10*$heatmap_clustering_bin_size;
my @multiSampleClustering=();
my $not_plot_index_matrix="";
my $heatmap_breaks=-1;
my @order_dendgram=();


my $density_bar="";
my @category_names=();
my @sample_names=();
my @experiment_names=();
my @rep_num_experiments=();
my $enrich_max=4;
my $logscale="";
my $capLimit=-1;

print STDERR "perl alignWigToBed.pl ";
my $cmd_root=join " ", @ARGV;
print STDERR "$cmd_root\n\n";

GetOptions( 
			"locs=s" => \@locs,
			"prefixs=s" => \@prefixs,
			"coverages=s" => \@coverages,
			"omit_align_step" => \$omit_align_step,
			"omit_plot_step" => \$omit_plot_step,
			"nomeseq" => \$nomeseq,
			"low_coverage" => \$low_coverage,
			"include_no_data_line"=> \$include_no_data_line,
			"adjust_center=i" =>\$adjust_center,
			"adjust_to_highest" => \$adjust_to_highest,
			"adjust_center_window=i" => \$adjust_center_window,
			"location_after_adjust_center=s" => \$location_after_adjust_center,
			"motif_freq=i" => \@motif_freq,
			"mask_matrix=s" => \$mask_matrix,
			"r_script=s" => \$r_script,
			"ucsc_script=s" => \$ucsc_script,
			"sort_perl_script=s" => \$sort_perl_script,
			"result_dir=s" => \$result_dir,
			"smooth" => \$smooth,
			"data_matrix_scale=i" => \$data_matrix_scale,
			"bin_size=i" => \$bin_size,
			"bin_size_align=i" => \$bin_size_align,
			
			"average" => \$average,
			"auto_scale" => \$auto_scale,
			"plot_x_axis_scale=i" => \$plot_x_axis_scale,
			"plot_x_axis_step=i" => \$plot_x_axis_step,
			"plot_y_axis_max=f" => \$plot_y_axis_max,
			"plot_y_axis_min=f" => \$plot_y_axis_min,
			"plot_y_axis_step=f" => \$plot_y_axis_step,
			"colors=s" => \@colors,
			"lengends=s" => \@lengends,
			"lty=s" => \@lty,
			"alignment_mode=i" => \$alignment_mode,
			"heatmap" => \$heatmap,
			"heatmap_anno=s" => \@heatmap_anno,
			"heatmap_anno_col=s" => \@heatmap_anno_col,
			"heatmap_anno_name=s" => \@heatmap_anno_name,
			"heatmap_col=s" => \@heatmap_col,
			"heatmap_cluster=i" => \@heatmap_cluster,
			
			"heatmap_with_reps" => \$heatmap_with_reps,
			"addAverage" => \$addAverage,
			"heatmap_keys=s" => \@heatmap_keys,
			"heatmap_ylab=s" => \@heatmap_ylab,
			"heatmap_ymin=f" => \@heatmap_ymin,
			"heatmap_ymax=f" => \@heatmap_ymax,
			"heatmap_logscale=s" => \@heatmap_logscale,
			"heatmap_capLimit=s" => \@heatmap_capLimit,
			"heatmap_capUpLimit=f" => \@heatmap_capUpLimit,
			"heatmap_capDownLimit=f" => \@heatmap_capDownLimit,
			"heatmap_autoscale=s" => \@heatmap_autoscale,
			"heatmap_output_subCluster" => \$heatmap_output_subCluster,
			"heatmap_regionToCluster_low=i" => \$heatmap_regionToCluster_low,
			"heatmap_regionToCluster_high=i" => \$heatmap_regionToCluster_high,
			"fileNumToPrintSideBar=i" => \$fileNumToPrintSideBar,
			"heatmap_png" => \$heatmap_png,
			"preClustering" => \$preClustering,
			"twoStepClustering" => \$twoStepClustering,
			"heatmap_clustering_matrix=i" => \$heatmap_clustering_matrix,
			"heatmap_clustering_scale=i" => \$heatmap_clustering_scale,
			"heatmap_clustering_bin_size=i" => \$heatmap_clustering_bin_size,
			"heatmap_clustering_bin_size_align=i" => \$heatmap_clustering_bin_size_align,
			"heatmap_file_num_to_cluster=i" => \@heatmap_file_num_to_cluster,
			"heatmap_normalization_by_mean" => \$heatmap_normalization_by_mean,
			"heatmap_row_order=s" => \$heatmap_row_order, 
			"heatmap_max_occ_order" => \$heatmap_max_occ_order,
			"multiSampleClustering=i" => \@multiSampleClustering,
			"not_plot_index_matrix" => \$not_plot_index_matrix,
			"heatmap_breaks=i" => \$heatmap_breaks,
			"order_dendgram=i" => \@order_dendgram,
			
			"density_bar" => \$density_bar,
			"category_names=s" => \@category_names,
			"sample_names=s" => \@sample_names,
			"experiment_names=s" => \@experiment_names,
			"rep_num_experiments=i" => \@rep_num_experiments,
			
			"enrich_max=f" => \$enrich_max,
			"logscale" => \$logscale,
			"capLimit=f" => \$capLimit,
			);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) < 1 ) {
    print STDERR "Wrong number of arguments\n";
    usage();
}
if(scalar(@locs) == 0 || (scalar(@prefixs) == 0 && ($density_bar eq "" && $heatmap_with_reps eq ""))){
	print STDERR "Genomic feature files and its associated prefix have not been specified yet. Please specified it by --locs feature.bed!!!\n\n";
    usage();
}
elsif(scalar(@locs) != scalar(@prefixs) && ($density_bar eq "" && $heatmap_with_reps eq "")){
	print STDERR "Number of genomic feature files and its associated prefix names are different!!!\n\n";
    usage();
}

#my $feature = shift(@ARGV);
#my $prefix = shift(@ARGV);
my @input_wigs = @ARGV;

if(scalar(@coverages) > 0 and scalar(@coverages) != scalar(@input_wigs)){
	print STDERR "Number of coverage files and input big wig files are different!!!\n\n";
    usage();
}

my %alignment_modes_hash = (
	1 => "Center",
	2 => "FiveEnd",
	3 => "ThreeEnd",
);
my @output_matrixs=();

for(my $w=1;$w<=(scalar(@locs)-1);$w++){
	push(@motif_freq,@motif_freq) if scalar(@motif_freq)>0;
	push(@lty,@lty) if scalar(@lty)>0;
	push(@lengends,@lengends) if scalar(@lengends)>0;
	push(@colors,@colors) if scalar(@colors)>0;
	push(@heatmap_col,@heatmap_col) if scalar(@heatmap_col)>0;
	push(@heatmap_keys,@heatmap_keys) if scalar(@heatmap_keys)>0;
	push(@heatmap_ylab,@heatmap_ylab) if scalar(@heatmap_ylab)>0;
	push(@heatmap_ymin,@heatmap_ymin) if scalar(@heatmap_ymin)>0;
	push(@heatmap_ymax,@heatmap_ymax) if scalar(@heatmap_ymax)>0;
	push(@heatmap_logscale,@heatmap_logscale) if scalar(@heatmap_logscale)>0;
	push(@heatmap_capLimit,@heatmap_capLimit) if scalar(@heatmap_capLimit)>0;
	push(@heatmap_autoscale,@heatmap_autoscale) if scalar(@heatmap_autoscale)>0;


}


if(scalar(@heatmap_cluster) <= 0){
		for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
					push(@heatmap_cluster, 4);
			}
}
else{
	my @tmp=@heatmap_cluster;
	@heatmap_cluster=();
	for(my $z=0;$z<scalar(@locs);$z++){
		for(my $ord=1;$ord<=scalar(@input_wigs);$ord++){
					push(@heatmap_cluster, $tmp[$z]);
		}	
	}
}

if(scalar(@multiSampleClustering) <= 0){
	push(@multiSampleClustering, 1);
}

##just do it to make density bar plot
if($density_bar ne "" || $heatmap_with_reps ne ""){
	&make_density_bar_or_heatmap();
	exit(0);
}



for(my $z=0;$z<scalar(@locs);$z++){
	my $feature = $locs[$z];
	my $prefix = $prefixs[$z];	
	my @order=();
	if($nomeseq ne ""){
	##default setting for NOMe-seq average plot
	my $step_tmp=4;
	if($low_coverage ne ""){
		$step_tmp=8;	
	}
	if(scalar(@order)<=0){
		for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord+=$step_tmp){
			if($heatmap ne ""){
				if($low_coverage ne ""){
				push(@order, "HCG");
				push(@order, "HCG_cov");
				push(@order, "HCG");
				push(@order, "HCG_cov");
				push(@order, "GCH");
				push(@order, "GCH_cov");
				push(@order, "GCH");
				push(@order, "GCH_cov");
				
				}
				else{
				push(@order, "HCG");
				push(@order, "HCG");
				push(@order, "GCH");
				push(@order, "GCH");
				
				}
			}
			else{
				if($low_coverage ne ""){
				push(@order, "HCG");
				push(@order, "HCG_cov");
				push(@order, "GCH");
				push(@order, "GCH_cov");
				
				}
				else{
					push(@order, "HCG");
					push(@order, "GCH");
					
				}
			}
			
			
			
		}
	}
	if(scalar(@lengends)<=0){
		for(my $ord=0;$ord<scalar(@input_wigs)*scalar(@locs);$ord+=$step_tmp){
			my $pre=basename($input_wigs[($ord % scalar(@locs))]);
			$pre=~s/(\w+)\S+/$1/;
			if($heatmap ne ""){
				if($low_coverage ne ""){
				push(@lengends, "HCG_${pre}");
				push(@lengends, "");
				push(@lengends, "HCG_${pre}");
				push(@lengends, "");
				push(@lengends, "GCH_${pre}");
				push(@lengends, "");
				push(@lengends, "GCH_${pre}");
				push(@lengends, "");
				
				}
				else{
				push(@lengends, "HCG_${pre}");
				push(@lengends, "HCG_${pre}");
				push(@lengends, "GCH_${pre}");
				push(@lengends, "GCH_${pre}");
				
				}
			}
			else{
				if($low_coverage ne ""){
				push(@lengends, "HCG_${pre}");
				push(@lengends, "");
				push(@lengends, "GCH_${pre}");
				push(@lengends, "");
				
				}
				else{
				push(@lengends, "HCG_${pre}");
				push(@lengends, "GCH_${pre}");
				
				}
			}
			
			
		}
	}
	if(scalar(@colors)<=0){
		for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord+=$step_tmp){
			if($heatmap ne ""){
				if($low_coverage ne ""){
				push(@colors, "black");
				push(@colors, "");
				push(@colors, "black");
				push(@colors, "");
				push(@colors, "#00CC99");
				push(@colors, "");
				push(@colors, "#00CC99");
				push(@colors, "");
				
				}
				else{
				push(@colors, "black");
				push(@colors, "black");
				push(@colors, "#00CC99");
				push(@colors, "#00CC99");
				}
			}
			else{
				if($low_coverage ne ""){
				push(@colors, "black");
				push(@colors, "");
				push(@colors, "#00CC99");
				push(@colors, "");
				
				}
				else{
				push(@colors, "black");
				push(@colors, "#00CC99");
				}
			}
			
			

		}
	}
	if(scalar(@lty)<=0){
		for(my $ord=1, my $lt=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord+=$step_tmp,$lt++){
			if($heatmap ne ""){
				if($low_coverage ne ""){
				push(@lty, $lt);
				push(@lty, 0);
				push(@lty, $lt);
				push(@lty, 0);
				push(@lty, $lt);
				push(@lty, 0);
				push(@lty, $lt);
				push(@lty, 0);
				}
				else{
				push(@lty, $lt);
				push(@lty, $lt);
				push(@lty, $lt);
				push(@lty, $lt);
				}
			}
			else{
				if($low_coverage ne ""){
				push(@lty, $lt);
				push(@lty, 0);
				push(@lty, $lt);
				push(@lty, 0);
				}
				else{
				push(@lty, $lt);
				push(@lty, $lt);
				}
			}
			

		}
	}
			

				
	
	#if(scalar(@input_wigs) != scalar(@lengends)){
	#	die "NOMe-seq mode should  have even number of input big wig files\n";
	#}
}

#if($low_coverage ne ""){
#	if(scalar(@motif_freq)<=0){
#		for(my $ord=0;$ord<scalar(@input_wigs);$ord+=2){
#				my $ord2=$ord+1;
#				$motif_freq[$ord]=0;
#				$motif_freq[$ord2]=1;
#		}
#	}
#}

##default setting for heatmap plot
if($heatmap ne ""){
	if($nomeseq ne ""){
		$fileNumToPrintSideBar=scalar(@input_wigs)*scalar(@locs); ##the last one to print side bar..
		if(scalar(@heatmap_col)<=0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord+=4){
				push(@heatmap_col, "blue2yellow");
				push(@heatmap_col, "blue2yellow");
				push(@heatmap_col, "white2darkgreen");
				push(@heatmap_col, "white2darkgreen");
			}
		}
		if(scalar(@heatmap_ylab)<=0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord+=4){
				push(@heatmap_ylab, "Methylation");
				push(@heatmap_ylab, "Methylation");
				push(@heatmap_ylab, "Accessibility");
				push(@heatmap_ylab, "Accessibility");
			}
		}
		if(scalar(@heatmap_keys)<=0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
				if($ord>=scalar(@input_wigs)*scalar(@locs)-1){
					push(@heatmap_keys, "TRUE");
				}
				else{
					push(@heatmap_keys, "FALSE");
				}

			}
		}

			$heatmap_regionToCluster_low=-240 if($heatmap_regionToCluster_low eq "");
			$heatmap_regionToCluster_high=440 if($heatmap_regionToCluster_high eq "");
		
		
	}
	else{
		if(scalar(@heatmap_col)<=0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
				push(@heatmap_col, "white2red");
			}
		}
		if(scalar(@heatmap_ylab)<=0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
				push(@heatmap_ylab, "Density");
			}
		}
		if(scalar(@heatmap_keys)<=0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
					push(@heatmap_keys, "TRUE");
			}
		}
		
			$heatmap_regionToCluster_low=-200 if($heatmap_regionToCluster_low eq "");
			$heatmap_regionToCluster_high=200 if($heatmap_regionToCluster_high eq "");
		
	}
	if(scalar(@heatmap_ymin) <= 0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
					push(@heatmap_ymin, 0);
			}
		}
		if(scalar(@heatmap_ymax) <= 0){
			for(my $ord=1;$ord<=scalar(@input_wigs)*scalar(@locs);$ord++){
					push(@heatmap_ymax, 100);
			}
		}
	
	##add sepereation of color here..
	#for(my $ord=0;$ord<scalar(@heatmap_col);$ord++){
		#$heatmap_col[$ord].="\\(75\\)";
	#}
}


if(scalar(@order)<=0){
	@order=(1..scalar(@input_wigs));
}

if(scalar(@heatmap_logscale)<=0){
	foreach (1..scalar(@input_wigs)*scalar(@locs)){
		push(@heatmap_logscale,"FALSE");
	}
}
if(scalar(@heatmap_capLimit)<=0){
	foreach (1..scalar(@input_wigs)*scalar(@locs)){
		push(@heatmap_capLimit,"FALSE");
	}
}
if(scalar(@heatmap_autoscale)<=0){
	foreach (1..scalar(@input_wigs)*scalar(@locs)){
		push(@heatmap_autoscale,"FALSE");
	}
}

if(scalar(@colors)<=0){
	foreach (1..scalar(@input_wigs)*scalar(@locs)){
		push(@colors,"black");
	}
}

if(scalar(@lty)<=0){
	foreach (1..scalar(@input_wigs)*scalar(@locs)){
		push(@lty,1);
	}
}
if(scalar(@motif_freq)<=0){
	foreach (1..scalar(@input_wigs)*scalar(@locs)){
		push(@motif_freq,0);
	}
}

my $span = int(($data_matrix_scale*2)/$bin_size_align+1);



my $j=0;
my $total_line=`wc -l $feature`;
chomp($total_line);
$total_line =~ s/(\d+)\s+\S+/$1/;
foreach my $input_wig(@input_wigs){
	my $input_prefix = basename($input_wig);
	$input_prefix =~s/(\w+)\S+/$1/;
	my $prefix_tmp = $prefix.".".$input_prefix;
	my $suffix = basename($feature);
	$suffix =~s/(\w+)\S+/$1/; 
	
	my $output_matrix = $result_dir."$prefix_tmp.alignedTo.$suffix.$data_matrix_scale.$order[$j].txt";
	if($j == 0 and $preClustering eq "" and ($heatmap ne "" or $heatmap_with_reps ne "")){
		$span = int(($heatmap_clustering_matrix*2)/$heatmap_clustering_bin_size_align+1);
		$output_matrix = $result_dir."$prefix_tmp.alignedTo.$suffix.$heatmap_clustering_matrix.$order[$j].txt";
	}
	else{
		$span = int(($data_matrix_scale*2)/$bin_size_align+1);
	}

	
	if($motif_freq[$j]==1){
		$output_matrix = $result_dir."$prefix_tmp.alignedTo.$suffix.$data_matrix_scale.$order[$j].motif_freq.txt";
	}
	push(@output_matrixs,$output_matrix);
	
	if($omit_align_step eq ""){
		open(FH,"<$feature") or die "can't open $feature file: $!\n";
		open(OUT, ">$output_matrix") or die "can't open $output_matrix file: $!\n";
		if($location_after_adjust_center ne ""){
			
			open(LOC, ">$location_after_adjust_center") or die "can't open new location file, $location_after_adjust_center: $!\n";
		}
		my $line=1;
		while(<FH>){
			chomp;
			next if($_=~/^track/ or $_ =~/^\#/);
			my @splitin=split "\t";
			my $cor1;
			my $cor2;
			if($feature =~ /\.bed$/ or $feature =~ /\.bedgraph$/ or $feature =~ /\.bedGraph$/){
				$cor1=$splitin[1];
				$cor2=$splitin[2];
			}
			elsif($feature =~ /\.gtf$/ or $feature =~ /\.gff$/){
				$cor1=$splitin[3]-1;
				$cor2=$splitin[4];
			}
			else{
				die "Not support such a format yet. Please use Bed format or GTF/GFF format!!!\n";
			}
			
			my $chr;
			my $start;
			my $end;
			die "Not valid bed format at Line: $line\n" if(scalar(@splitin) < 3);
			$chr=$splitin[0];
			my $strand=".";
			if(scalar(@splitin) >= 6){
					$strand=$splitin[5];
					if($feature =~ /\.gtf$/ or $feature =~ /\.gff$/){
						$strand=$splitin[6];	
					}
			}
			if((!(scalar(@coverages) > 0 and $alignment_mode == 4)) and (scalar(@splitin) >= 6 && $splitin[5] eq '-' || (($feature =~ /\.gtf$/ or $feature =~ /\.gff$/) && $splitin[6] eq '-'))){
					my $tem = $cor1;
					$cor1=$cor2;
					$cor2=$tem;
			}
			
			if($alignment_mode == 1 || $alignment_mode == 5){
				$start=int(($cor1+$cor2)/2)-$data_matrix_scale;
				$end=int(($cor1+$cor2)/2) + 1 + $data_matrix_scale;
				if($j == 0 and $preClustering eq ""  and $heatmap ne ""){
					$start=int(($cor1+$cor2)/2)-$heatmap_clustering_matrix;
					$end=int(($cor1+$cor2)/2) + 1 + $heatmap_clustering_matrix;
				}
				
			}
			elsif($alignment_mode == 2){
				$start=$cor1-$data_matrix_scale;
				$end=$cor1 + 1 + $data_matrix_scale;
				if($j == 0 and $preClustering eq ""  and $heatmap ne ""){
					$start=$cor1-$heatmap_clustering_matrix;
					$end=$cor1 + 1 + $heatmap_clustering_matrix;
				}
			}
			elsif($alignment_mode == 3){
				$start=$cor2 -1 - $data_matrix_scale;
				$end=$cor2 + $data_matrix_scale;
				if($j == 0 and $preClustering eq ""  and $heatmap ne ""){
					$start=$cor2 -1 - $heatmap_clustering_matrix;
					$end=$cor2 + $heatmap_clustering_matrix;
				}
			}
			elsif($alignment_mode == 4 || $alignment_mode == 6){
				$start=$cor1;
				$end=$cor2;
				if($alignment_mode == 4){
					$span = 1;
				}else{
					$span = int($data_matrix_scale/$bin_size_align+1);
				}
				
				if(scalar(@coverages) > 0){
					$span = $end - $start;
				}
			}
			else{
				print STDERR "unsupport alignment mode!!\n";
				exit(1);
			}


			
		
			my $align_cmd="$ucsc_script $input_wig $chr $start $end $span 2>/dev/null |";
			if($motif_freq[$j]==1){
				$align_cmd="$ucsc_script $input_wig $chr $start $end $span -type=coverage 2>/dev/null |";	
			}
			my $tmp="";
			open(P, $align_cmd) or die;
			while (<P>) {
  				chomp;
  				$tmp = $_;
			}
			close(P);
			
			my $mask="";
			if($mask_matrix ne ""){
				my $mask_cmd="$ucsc_script $mask_matrix $chr $start $end $span 2>/dev/null |";
				if($motif_freq[$j]==1){
					$mask_cmd="$ucsc_script $mask_matrix $chr $start $end $span -type=coverage 2>/dev/null |";
				}
				open(X, $mask_cmd) or die;
				while (<X>) {
  					chomp;
  					$mask = $_;
				}
				close(X);
			}
			
			if($tmp ne ""){
				my @data_array=split "\t",$tmp;
				my @mask_array = ();
				
				##adjust the mean value by coverage in align mode 4
				if(scalar(@coverages) > 0 and $alignment_mode == 4){
					@data_array = &adjust_value_by_coverage($coverages[$j],$chr,$start,$end, $span, $tmp);
				}
				
				##recenter to the lowest point, it can't work with mask matrix right now..
				if($adjust_center ne ""){				
					my @new_data_array = &recenter($input_wig,$chr,$start,$end, $strand, $span, $tmp);
					if(scalar(@new_data_array) == scalar(@data_array)){
						@data_array = @new_data_array;
					}
				}
				
				
				##if mask the data
				if($mask_matrix ne "" && $mask ne ""){
					@mask_array = split "\t",$mask;
					if(scalar(@mask_array) != scalar(@data_array)){
						die "Line $line: mask matrix do not have the same length as data matrix ($input_wig)\n";
					}
					for(my $num=0;$num <= $#data_array; $num++){
						if($mask_array[$num] ne "n/a"){
							$data_array[$num] = "NA";
						}
					}
				}
				
				##if extend data for alignment_mode 6
				if($alignment_mode == 6){
					my @new_data_array = &alignment_extend($input_wig,$chr, $start,$end,int($data_matrix_scale/2),int($data_matrix_scale/($bin_size_align*2)), @data_array);
					@data_array = @new_data_array;
				}
				
				if((!(scalar(@coverages) > 0 and $alignment_mode == 4)) and (scalar(@splitin) >= 6 && $splitin[5] eq '-' || (($feature =~ /\.gtf$/ or $feature =~ /\.gff$/) && $splitin[6] eq '-'))){
					@data_array = reverse(@data_array);
				}
				
				if($heatmap_normalization_by_mean ne ""){
					my @new_data_array = &normalize($input_wig, @data_array);
					@data_array = @new_data_array;
				}
				
				$tmp=join "\t",@data_array;
				if($alignment_mode == 5){
					my $center=int(scalar(@data_array)/2);
					my $closest="";
					#print "$tmp\n$center\n";
					for(my $dist=0;$dist<$center;$dist++){
						#print "$data_array[$center+$dist]\t$data_array[$center-$dist]\t";
						if($data_array[$center+$dist] ne "n/a"){
							$closest=$data_array[$center+$dist];
							last;
						}elsif($data_array[$center-$dist] ne "n/a"){
							$closest=$data_array[$center-$dist];
							last;
						}
					}
					#print "\n";
					if($closest ne ""){
						$tmp=$closest;
					}else{
						next;
					}
					
				}
				
				if($motif_freq[$j]==1){
					#$tmp =~ s/[-]\d+[.\d+]/1/g;
					$tmp =~ s/n\/a/0/g;
					#$tmp =~ s/n\/a/NA/g;
				}
				else{
					$tmp =~ s/n\/a/NA/g;
				}
				
				$tmp = "$chr\t$cor1\t$cor2\t$strand\t$tmp\n";
				print OUT $tmp;
				
			}
			elsif($include_no_data_line ne ""){
				
				if($motif_freq[$j]==1){
					$tmp = join "\t", (0) x $span;
				}
				else{
					$tmp = join "\t", ("NA") x $span;
				}
				$tmp = "$chr\t$cor1\t$cor2\t$strand\t$tmp\n";
				print OUT $tmp;	
			}
			$line++;
			if($line % 1000 == 0){
				print STDERR "Aligning $input_wig to line $line ... ...\t $total_line genomic features in total\n";
			}
		}
		
		
		if($location_after_adjust_center ne ""){
			close(LOC);	
		}
		close(FH);
		close(OUT);
		
	}
	$j++;	
}

}
if($omit_plot_step eq ""){
	##make average plot
	if($average ne ""){
		my $r_cmd = "R --no-restore --no-save --args wd=$result_dir prefix=$prefixs[0] step=$bin_size bin_size_align=$bin_size_align scale=$plot_x_axis_scale axistep=$plot_x_axis_step smooth=$smooth ";	
	if($auto_scale ne ""){
		$r_cmd .= "autoScale=TRUE ";
	}
	else{
		$r_cmd .= "autoScale=FALSE y_scale_max=$plot_y_axis_max y_scale_min=$plot_y_axis_min y_step=$plot_y_axis_step ";
	}
	if($low_coverage ne ""){
		$r_cmd .= "lowCov=TRUE ";
	}
	foreach my $output_matrix(@output_matrixs){
		my $fn=basename($output_matrix);
		$r_cmd .= "fn=$fn ";
	}
	foreach my $color(@colors){
		$r_cmd .= "color=$color ";
	}
	foreach my $legend_name(@lengends){
		$r_cmd .= "legendName=$legend_name ";
	}	
	foreach my $line_type(@lty){
		$r_cmd .= "line_types=$line_type ";
	}

	$r_cmd .= "< $r_script \n";
	print STDERR $r_cmd;
	system($r_cmd);
	}
	
	
	##do heatmap
	if($heatmap ne ""){
		
		my $heatmap_cmd = "R --no-restore --no-save --args wd=$result_dir move_step=$bin_size bin_size_align=$bin_size_align scale=$plot_x_axis_scale  ";
		#$heatmap_cmd .= "y_max=$plot_y_axis_max y_min=$plot_y_axis_min ";
		my $breaks=scalar(@input_wigs);
		$heatmap_cmd .= "regionToClusterLow=$heatmap_regionToCluster_low regionToClusterHigh=$heatmap_regionToCluster_high fileNumToPrintSideBar=$fileNumToPrintSideBar breaks=$breaks ";

		if($heatmap_png ne ""){
			$heatmap_cmd .= "pdfOut=FALSE ";
		}
		if($preClustering ne ""){
			$heatmap_cmd .= "preClustering=TRUE ";
		}
		if($twoStepClustering ne ""){
			$heatmap_cmd .= "twoStepClustering=TRUE ";
		}
		if($heatmap_row_order ne ""){
			$heatmap_cmd .= "heatmap_row_order=$heatmap_row_order ";
		}
		if($heatmap_max_occ_order ne ""){
			$heatmap_cmd .= "orderByMaxOccDisToTss=TRUE ";
		}
		foreach my $sampleNum(@multiSampleClustering){
			$heatmap_cmd .= "multiSampleClustering=$sampleNum ";
		}
		
		foreach my $output_matrix(@output_matrixs){
			my $fn=basename($output_matrix);
			$heatmap_cmd .= "inputFn=$fn ";
		}	
		foreach my $prefix(@prefixs){
			$prefix="heatmap.$prefix";
			foreach (1..scalar(@input_wigs)){
				$heatmap_cmd .= "prefix=$prefix ";
			}
			
		}	
		
		foreach my $clust(@heatmap_cluster){
			$heatmap_cmd .= "clusterNum=$clust ";
		}
		foreach my $col(@heatmap_col){
			$heatmap_cmd .= "heatMapCols=$col ";
		}
		foreach my $key(@heatmap_keys){
			$heatmap_cmd .= "keys=$key ";
		}
		foreach my $ylab(@heatmap_ylab){
			$heatmap_cmd .= "ylabForAvePlot=$ylab ";
		}
		foreach my $ymin(@heatmap_ymin){
			$heatmap_cmd .= "y_min=$ymin ";
		}
		foreach my $ymax(@heatmap_ymax){
			$heatmap_cmd .= "y_max=$ymax ";
		}
		if(scalar(@heatmap_autoscale) > 0){
			foreach my $autoscal(@heatmap_autoscale){
				$heatmap_cmd .= "autoScale=$autoscal ";
			}
		}
		if(scalar(@heatmap_logscale) > 0){
			foreach my $logscal(@heatmap_logscale){
				$heatmap_cmd .= "logScale=$logscal ";
			}
		}
		if(scalar(@heatmap_capLimit) > 0){
			foreach my $cap(@heatmap_capLimit){
				$heatmap_cmd .= "capLimit=$cap ";
			}
		}
		if(scalar(@heatmap_capUpLimit) > 0){
			foreach my $cap(@heatmap_capUpLimit){
				$heatmap_cmd .= "capUpLimit=$cap ";
			}
		}
		if(scalar(@heatmap_capDownLimit) > 0){
			foreach my $cap(@heatmap_capDownLimit){
				$heatmap_cmd .= "capDownLimit=$cap ";
			}
		}
		
		
		if($heatmap_output_subCluster ne ""){
			$heatmap_cmd .= "outputSubClusterCordinate=TRUE ";
		}
	


		##if provide side bar..
		if(scalar(@heatmap_anno) > 0){
			foreach my $anno(@heatmap_anno){
				$heatmap_cmd .= "rowSideFiles=$anno ";
			}
		}
		if(scalar(@heatmap_anno_col) > 0){
			foreach my $anno_col(@heatmap_anno_col){
				$heatmap_cmd .= "colToUses=$anno_col ";
			}
		}
		if(scalar(@heatmap_anno_name) > 0){
			foreach my $anno_name(@heatmap_anno_name){
				$heatmap_cmd .= "RowSideColorsName=$anno_name ";
			}
		}

		if($preClustering eq ""){
			$heatmap_cmd .= "heatmap_clustering_scale=$heatmap_clustering_scale heatmap_clustering_bin_size=$heatmap_clustering_bin_size heatmap_clustering_bin_size_align=$heatmap_clustering_bin_size_align ";
			
			
		}
	
		$heatmap_cmd .= "< $heatmap_r_script \n";
		print STDERR $heatmap_cmd;
		system($heatmap_cmd);
	}
	
	##do heatmap_with_reps
	if($heatmap_with_reps ne ""){
		my $heatmap_cmd = "R --no-restore --no-save --args wd=$result_dir move_step=$bin_size bin_size_align=$bin_size_align scale=$plot_x_axis_scale ";
		#$heatmap_cmd .= "y_max=$plot_y_axis_max y_min=$plot_y_axis_min ";
		my $breaks=scalar(@input_wigs);
		$heatmap_cmd .= "regionToClusterLow=$heatmap_regionToCluster_low regionToClusterHigh=$heatmap_regionToCluster_high fileNumToPrintSideBar=$fileNumToPrintSideBar breaks=$breaks ";
		
		foreach my $output_matrix(@output_matrixs){
			my $fn=basename($output_matrix);
			$heatmap_cmd .= "inputFn=$fn ";
		}	
		foreach my $prefix(@prefixs){
			$prefix="heatmap.$prefix";
			foreach (1..scalar(@input_wigs)){
				$heatmap_cmd .= "prefix=$prefix ";
			}
			
		}	
		
		foreach my $clust(@heatmap_cluster){
			$heatmap_cmd .= "clusterNum=$clust ";
		}
		foreach my $col(@heatmap_col){
			$heatmap_cmd .= "heatMapCols=$col ";
		}


		if(scalar(@heatmap_logscale) > 0){
			foreach my $logscal(@heatmap_logscale){
				$heatmap_cmd .= "logScale=$logscal ";
			}
		}
		if(scalar(@heatmap_capLimit) > 0){
			foreach my $cap(@heatmap_capLimit){
				$heatmap_cmd .= "capLimit=$cap ";
			}
		}
		for(my $i=0;$i<scalar(@category_names);$i++){

			$heatmap_cmd.= "categoryNames=$category_names[$i] ";
		}
		
		foreach (1..scalar(@category_names)){
			foreach my $experiment_name(@experiment_names){			
				$heatmap_cmd.= "sampleNames=$experiment_name ";

			}
		}
		

		$heatmap_cmd .= "heatmap_clustering_scale=$heatmap_clustering_scale heatmap_clustering_bin_size=$heatmap_clustering_bin_size heatmap_clustering_bin_size_align=$heatmap_clustering_bin_size_align ";
		
		if($addAverage ne ""){
			$heatmap_cmd .= "addAverage=TRUE ";
		}

		$heatmap_cmd .= "< $heatmap_with_reps_r_script \n";
		print STDERR $heatmap_cmd;
		system($heatmap_cmd);
	}
}


sub recenter{
		my $input=shift @_;
		my $chr=shift @_;
		my $start=shift @_;
		my $end=shift @_;
		my $strand=shift @_;
		my $span=shift @_; 
		my $data=shift @_;
		my @old_array=split "\t",$data;
		my $old_loc=$chr.":".($start+$data_matrix_scale)."-".($end-$data_matrix_scale).":".$strand;
		my $lowest="";
		my $lowest_point="";
		for(my $i=0+int(($data_matrix_scale-$adjust_center)/$bin_size_align);$i<=$#old_array-int(($data_matrix_scale-$adjust_center)/$bin_size_align);$i++){			
			my $sum=0;
			my $count=0;
			for(my $j=$i;$j<=$i+int($adjust_center_window/$bin_size_align);$j++){
				if($old_array[$j] !~ /n\/a/){
					$sum+=$old_array[$j];
					$count++;
				}	
			}
			if($count>0){
				$sum/=$count;
				if($lowest eq ""){
					$lowest=$sum;
					$lowest_point=$start+$i*$bin_size_align+$adjust_center_window/2;
				}
				elsif($sum != 0){
					if($adjust_to_highest eq "" && $sum < $lowest){
						$lowest=$sum;
						$lowest_point=$start+$i*$bin_size_align+$adjust_center_window/2;
					}elsif($adjust_to_highest ne "" && $sum > $lowest){
						$lowest=$sum;
						$lowest_point=$start+$i*$bin_size_align+$adjust_center_window/2;
					}
					
				}
			}
			
		}
		if($lowest_point ne ""){
			$start=$lowest_point-$data_matrix_scale;
			$end=$lowest_point + 1 + $data_matrix_scale;
			my $align_cmd="$ucsc_script $input $chr $start $end $span |";
			my $tmp="";
			open(P, $align_cmd) or die;
			while (<P>) {
  				chomp;
  				$tmp = $_;
			}
			close(P);
			if($location_after_adjust_center ne ""){
				my $new_start=$lowest_point;
				my $new_end = $lowest_point+1;
			
				print LOC "$chr\t$new_start\t$new_end\t$old_loc\t.\t$strand\n";
			}
			my @new_array=split "\t",$tmp;
			return @new_array;
		}
		else{
			return -1;
		}
		
}


sub recenter_by_2_highest_center{
		my $input=shift @_;
		my $chr=shift @_;
		my $start=shift @_;
		my $end=shift @_;
		my $strand=shift @_;
		my $span=shift @_; 
		my $data=shift @_;
		my @old_array=split "\t",$data;
		my $old_loc=$chr.":".($start+$data_matrix_scale)."-".($end-$data_matrix_scale).":".$strand;
		my $lowest="";
		my $lowest_point="";
		for(my $i=0+int(($data_matrix_scale-$adjust_center)/$bin_size_align);$i<=$#old_array-int(($data_matrix_scale-$adjust_center)/$bin_size_align);$i++){			
			my $sum=0;
			my $count=0;
			for(my $j=$i;$j<=$i+int($adjust_center_window/$bin_size_align);$j++){
				if($old_array[$j] !~ /n\/a/){
					$sum+=$old_array[$j];
					$count++;
				}	
			}
			if($count>0){
				$sum/=$count;
				if($lowest eq ""){
					$lowest=$sum;
					$lowest_point=$start+$i*$bin_size_align+$adjust_center_window/2;
				}
				elsif($sum != 0){
					if($adjust_to_highest eq "" && $sum < $lowest){
						$lowest=$sum;
						$lowest_point=$start+$i*$bin_size_align+$adjust_center_window/2;
					}elsif($adjust_to_highest ne "" && $sum > $lowest){
						$lowest=$sum;
						$lowest_point=$start+$i*$bin_size_align+$adjust_center_window/2;
					}
					
				}
			}
			
		}
		if($lowest_point ne ""){
			$start=$lowest_point-$data_matrix_scale;
			$end=$lowest_point + 1 + $data_matrix_scale;
			my $align_cmd="$ucsc_script $input $chr $start $end $span |";
			my $tmp="";
			open(P, $align_cmd) or die;
			while (<P>) {
  				chomp;
  				$tmp = $_;
			}
			close(P);
			if($location_after_adjust_center ne ""){
				my $new_start=$lowest_point;
				my $new_end = $lowest_point+1;
			
				print LOC "$chr\t$new_start\t$new_end\t$old_loc\t.\t$strand\n";
			}
			my @new_array=split "\t",$tmp;
			return @new_array;
		}
		else{
			return -1;
		}
		
}

sub adjust_value_by_coverage{
		my $input=shift @_;
		my $chr=shift @_;
		my $start=shift @_;
		my $end=shift @_;
		my $span=shift @_; 
		my $data=shift @_;
		my @old_array=split "\t",$data;
		
		my $align_cmd="$ucsc_script $input $chr $start $end $span |";
		my $tmp="";
		open(P, $align_cmd) or die;
		while (<P>) {
  				chomp;
  				$tmp = $_;
		}
		close(P);
		
		$tmp =~ s/n\/a/0/g;
		my @coverag_array=split "\t",$tmp;
		if(scalar(@old_array) != scalar(@coverag_array)){
			die "the coverage track produced array's length is different from original input big wig data\n!";
		}
		my $num_c_pattern=0;
		my $numCT=0;
		my $numC=0;
		my $total_methy=0;
		for(my $i=0;$i<scalar(@old_array); $i++){
			if($coverag_array[$i] > 0){
				$num_c_pattern++;
				$numCT+=$coverag_array[$i];
				$numC+=int($coverag_array[$i]*$old_array[$i]/100);
				$total_methy+=$old_array[$i];
			}
		}
		my $methy_adjust=100*$numC/$numCT;
		#print OUT "\n$data\n$tmp\n";
		my @new=($methy_adjust,$numCT,$total_methy,$num_c_pattern);
		return @new;
}


sub make_density_bar_or_heatmap{
	
	my $sample_order=0;
	my @file_name_lists=();
	#$AoA[$i] = [ @tmp ];
	foreach my $input(@ARGV){ ##each sample's experiment list
		open(FL,"<$input") or die "can't open file list: $input. $!\n";
		my @ins=<FL>;
		chomp(@ins);
		close(FL);
		my $category_order=0;
		foreach my $loc(@locs){
			my $experiment_order=0;
			my $out_list=$result_dir."/";
			foreach my $prefix(@prefixs){
				$out_list.="$prefix.";
			}
			if($heatmap_with_reps ne ""){
				$out_list.=$category_names[$category_order].".filelist.txt";
				$file_name_lists[$category_order]=$out_list;
			}elsif($density_bar ne ""){
				$out_list.=$sample_names[$sample_order].".".$category_names[$category_order].".filelist.txt";
				#push(@file_name_lists,$out_list);
				$file_name_lists[$category_order][$sample_order]=$out_list;
			}
			
			
			open(OL,">$out_list") or die "cant write on file list: $out_list $!\n";

			foreach my $in(@ins){
				my @splitin=split "\t",$in;
				my $file=$splitin[0];
				my $dir=$splitin[1];
				my $file_prefix=basename($file);;
				$file_prefix=~s/([\w|-]+)\S+/$1/;
				my $data_type=$splitin[2];
				if($file ne "" && $file ne "NA"){
					my $output_matrix = $dir."/$prefixs[0].$file_prefix.alignedTo.$category_names[$category_order].$data_matrix_scale.txt";
					if($preClustering eq "" and ($heatmap ne "" or $heatmap_with_reps ne "")){
						foreach my $c(@multiSampleClustering){
							if($experiment_order == ($c-1)){
								$output_matrix = $dir."/$prefixs[0].$file_prefix.alignedTo.$category_names[$category_order].$heatmap_clustering_matrix.txt";
								last;
							}
						}
					}
					
					print OL "$output_matrix\t$data_type\n";
					if($omit_align_step eq ""){						
						&align_wig2loc($output_matrix,$loc, $file,$experiment_order);			
					}
				}else{
					print OL "NA\t$data_type\n";
				}
				$experiment_order++;
			}
			close(OL);
			$category_order++;
		}
		$sample_order++;
	}
	if($omit_plot_step eq ""){
		if($density_bar ne "" ){
			&plot_density_bar(@file_name_lists);
		}
		if($heatmap_with_reps ne ""){
			&plot_heatmap_reps(@file_name_lists);
		}
		
	}
		
}

sub plot_density_bar{
	my @file_name_lists = @_;
	my $plot_prefix="DensityBarPlot";
		foreach my $prefix(@prefixs){
				$plot_prefix.=".$prefix";
			}
			
		#$plot_prefix.=join ".",@sample_names;
		#$plot_prefix.=join ".",@category_names;
		my $r_cmd = "R --no-restore --no-save --args wd=$result_dir prefix=$plot_prefix step=$bin_size bin_size_align=$bin_size_align scale=$plot_x_axis_scale enrichScoreMax=$enrich_max capLimitPerc=$capLimit ";	
		if($logscale eq ""){
			 $r_cmd.= "logscale=F ";
		}else{
			$r_cmd.= "logscale=$logscale ";
		}
		
		#foreach my $file_name_list(@file_name_lists){
		#	$r_cmd.= "file_name_lists=$file_name_list ";
		#}
		for(my $i=0;$i<scalar(@category_names);$i++){
			for(my $j=0;$j<scalar(@sample_names);$j++){
				$r_cmd.= "file_name_lists=$file_name_lists[$i][$j] ";
			}
			$r_cmd.= "category_names=$category_names[$i] ";
		}
		foreach my $sample_name(@sample_names){
			$r_cmd.= "sample_names=$sample_name ";
		}
		foreach my $experiment_name(@experiment_names){
			$r_cmd.= "experiment_names=$experiment_name ";
		}
		
		foreach my $rep_num_experiment(@rep_num_experiments){
			$r_cmd.= "rep_num_experiment=$rep_num_experiment ";
		}

		$r_cmd .= "< $densitybar_r_script \n";
		print STDERR $r_cmd;
		system($r_cmd) == 0 || die "Unexpected termination of R script!!\n";
}

sub plot_heatmap_reps{
	my @file_name_lists = @_;
	my $plot_prefix="Heatmap";
	foreach my $prefix(@prefixs){
			$plot_prefix.=".$prefix";
	}
	
	my $heatmap_cmd = "R --no-restore --no-save --args wd=$result_dir move_step=$bin_size bin_size_align=$bin_size_align scale=$plot_x_axis_scale ";
		my $breaks=1;
		$heatmap_cmd .= "regionToClusterLow=$heatmap_regionToCluster_low regionToClusterHigh=$heatmap_regionToCluster_high fileNumToPrintSideBar=$fileNumToPrintSideBar ";
		if(scalar(@heatmap_capUpLimit) > 0){
			foreach my $cap(@heatmap_capUpLimit){
				$heatmap_cmd .= "capUpLimit=$cap ";
			}
		}
		if(scalar(@heatmap_capDownLimit) > 0){
			foreach my $cap(@heatmap_capDownLimit){
				$heatmap_cmd .= "capDownLimit=$cap ";
			}
		}
		
		if(scalar(@multiSampleClustering) > 0){
			foreach my $c(@multiSampleClustering){
				$heatmap_cmd .= "multiSampleClustering=$c ";
			}
		}
		
		my $category_order=0;
		
		foreach my $file_name_list(@file_name_lists){
			open(F,"<$file_name_list") or die "cant open file name list : $file_name_list: $!\n ";
			my @fin=<F>;
			chomp(@fin);
			close(F);
			my $experiment_order=0;
			foreach my $line(@fin){
				my @splitin=split "\t",$line;
				my $fn = $splitin[0];
				if($fn ne "" && $fn ne "NA"){
					$heatmap_cmd .= "inputFn=$fn ";
					my $prefix.="$plot_prefix.$category_names[$category_order]";
					$heatmap_cmd .= "prefix=$prefix ";
					$heatmap_cmd .= "clusterNum=$heatmap_cluster[$category_order] ";
					$heatmap_cmd .= "heatMapCols=$heatmap_col[$experiment_order] ";
					if(scalar(@heatmap_logscale) > 0){
							$heatmap_cmd .= "logScale=$heatmap_logscale[$experiment_order] ";
					}else{
						$heatmap_cmd .= "logScale=F ";
					}
					if(scalar(@heatmap_capLimit) > 0){
							$heatmap_cmd .= "capLimit=$heatmap_capLimit[$experiment_order] ";
					}else{
						$heatmap_cmd .= "capLimit=F ";
					}
					$heatmap_cmd.= "sampleNames=$experiment_names[$experiment_order] ";
				}
				$experiment_order++;
			}
			$breaks=$experiment_order;
			$category_order++;
		}
		
		if($heatmap_breaks == -1){
			$heatmap_breaks = $breaks;
		}
		
		
		$heatmap_cmd .= "breaks=$heatmap_breaks heatmap_clustering_scale=$heatmap_clustering_scale heatmap_clustering_bin_size=$heatmap_clustering_bin_size heatmap_clustering_bin_size_align=$heatmap_clustering_bin_size_align ";
		
		if($addAverage ne ""){
			$heatmap_cmd .= "addAverage=TRUE ";
		}
		if($not_plot_index_matrix ne ""){
			$heatmap_cmd .= "notPlotIndexMatrix=TRUE ";
		}
		
		
		if(scalar(@heatmap_anno)>0){
			foreach my $anno(@heatmap_anno){
				$heatmap_cmd .= "heatmap_anno=$anno ";
			}
		}
		if(scalar(@heatmap_anno_col)>0){
			foreach my $anno(@heatmap_anno_col){
				$heatmap_cmd .= "heatmap_anno_col=$anno ";
			}
		}
		if(scalar(@heatmap_anno_name)>0){
			foreach my $anno(@heatmap_anno_name){
				$heatmap_cmd .= "heatmap_anno_name=$anno ";
			}
		}
		if(scalar(@order_dendgram)>0){
			foreach my $order(@order_dendgram){
				$heatmap_cmd .= "order_dendgram=$order ";
			}
		}
	
	
		$heatmap_cmd .= "< $heatmap_with_reps_r_script \n";
		print STDERR $heatmap_cmd;
		system($heatmap_cmd);
		

		
}


sub align_wig2loc{
	my $output_matrix=shift @_;
	my $feature=shift @_;
	my $wig=shift @_;
	my $j=shift @_;
	
	my $total_line=`wc -l $feature`;
	chomp($total_line);
	$total_line =~ s/(\d+)\s+\S+/$1/;
	my $span = int(($data_matrix_scale*2)/$bin_size_align+1);
	my $index="";
	foreach my $c(@multiSampleClustering){
			if($j == ($c-1)){
				$index=$j;
				last;
			}
	}
	if($index ne "" and $preClustering eq "" and ($heatmap ne "" or $heatmap_with_reps ne "")){
			$span = int(($heatmap_clustering_matrix*2)/$heatmap_clustering_bin_size_align+1);
	}
	
	open(FH,"<$feature") or die "can't open $feature file: $!\n";
		open(OUT, ">$output_matrix") or die "can't open $output_matrix file: $!\n";
		if($location_after_adjust_center ne ""){
			
			open(LOC, ">$location_after_adjust_center") or die "can't open new location file, $location_after_adjust_center: $!\n";
		}
		my $line=1;
		while(<FH>){
			chomp;
			next if($_=~/^track/ or $_ =~/^\#/);
			my @splitin=split "\t";
			my $cor1;
			my $cor2;
			if($feature =~ /\.bed$/ or $feature =~ /\.bedgraph$/ or $feature =~ /\.bedGraph$/){
				$cor1=$splitin[1];
				$cor2=$splitin[2];
			}
			elsif($feature =~ /\.gtf$/ or $feature =~ /\.gff$/){
				$cor1=$splitin[3]-1;
				$cor2=$splitin[4];
			}
			else{
				die "Not support such a format yet. Please use Bed format or GTF/GFF format!!!\n";
			}
			
			my $chr;
			my $start;
			my $end;
			die "Not valid bed format at Line: $line\n" if(scalar(@splitin) < 3);
			$chr=$splitin[0];
			my $strand=".";
			if(scalar(@splitin) >= 6){
					$strand=$splitin[5];
					if($feature =~ /\.gtf$/ or $feature =~ /\.gff$/){
						$strand=$splitin[6];	
					}
			}
			
			if((!(scalar(@coverages) > 0 and $alignment_mode == 4)) and (scalar(@splitin) >= 6 && $splitin[5] eq '-' || (($feature =~ /\.gtf$/ or $feature =~ /\.gff$/) && $splitin[6] eq '-'))){
					my $tem = $cor1;
					$cor1=$cor2;
					$cor2=$tem;
				}
			
			if($alignment_mode == 1 || $alignment_mode == 5){
				$start=int(($cor1+$cor2)/2)-$data_matrix_scale;
				$end=int(($cor1+$cor2)/2) + 1 + $data_matrix_scale;
				if($index ne "" and $preClustering eq ""  and ($heatmap ne "" or $heatmap_with_reps ne "")){
					$start=int(($cor1+$cor2)/2)-$heatmap_clustering_matrix;
					$end=int(($cor1+$cor2)/2) + 1 + $heatmap_clustering_matrix;
				}
				
			}
			elsif($alignment_mode == 2){
				$start=$cor1-$data_matrix_scale;
				$end=$cor1 + 1 + $data_matrix_scale;
				if($index ne "" and $preClustering eq ""  and ($heatmap ne "" or $heatmap_with_reps ne "")){
					$start=$cor1-$heatmap_clustering_matrix;
					$end=$cor1 + 1 + $heatmap_clustering_matrix;
				}
			}
			elsif($alignment_mode == 3){
				$start=$cor2 -1 - $data_matrix_scale;
				$end=$cor2 + $data_matrix_scale;
				if($index ne "" and $preClustering eq ""  and ($heatmap ne "" or $heatmap_with_reps ne "")){
					$start=$cor2 -1 - $heatmap_clustering_matrix;
					$end=$cor2 + $heatmap_clustering_matrix;
				}
			}
			elsif($alignment_mode == 4 || $alignment_mode == 6){
				$start=$cor1;
				$end=$cor2;
				if($alignment_mode == 4){
					$span = 1;
				}
				
				if(scalar(@coverages) > 0){
					$span = $end - $start;
				}
			}
			else{
				print STDERR "unsupport alignment mode!!\n";
				exit(1);
			}

		
			my $align_cmd="$ucsc_script $wig $chr $start $end $span |";
			#if($motif_freq[$j]==1){
			#	$align_cmd="$ucsc_script $wig $chr $start $end $span -type=coverage |";	
			#}
			my $tmp="";
			open(P, $align_cmd) or die;
			while (<P>) {
  				chomp;
  				$tmp = $_;
			}
			close(P);
			
			my $mask="";
			if($mask_matrix ne ""){
				my $mask_cmd="$ucsc_script $mask_matrix $chr $start $end $span |";
				#if($motif_freq[$j]==1){
				#	$mask_cmd="$ucsc_script $mask_matrix $chr $start $end $span -type=coverage |";
				#}
				open(X, $mask_cmd) or die;
				while (<X>) {
  					chomp;
  					$mask = $_;
				}
				close(X);
			}
			
			if($tmp ne ""){
				my @data_array=split "\t",$tmp;
				my @mask_array = ();
				
				##adjust the mean value by coverage in align mode 4
				if(scalar(@coverages) > 0 and $alignment_mode == 4){
					@data_array = &adjust_value_by_coverage($coverages[$j],$chr,$start,$end, $span, $tmp);
				}
				
				##recenter to the lowest point, it can't work with mask matrix right now..
				if($adjust_center ne ""){				
					my @new_data_array = &recenter($wig,$chr,$start,$end, $strand, $span, $tmp);
					if(scalar(@new_data_array) == scalar(@data_array)){
						@data_array = @new_data_array;
					}
				}
				
				
				##if mask the data
				if($mask_matrix ne "" && $mask ne ""){
					@mask_array = split "\t",$mask;
					if(scalar(@mask_array) != scalar(@data_array)){
						die "Line $line: mask matrix do not have the same length as data matrix ($wig)\n";
					}
					for(my $num=0;$num <= $#data_array; $num++){
						if($mask_array[$num] ne "n/a"){
							$data_array[$num] = "NA";
						}
					}
				}
				
				if((!(scalar(@coverages) > 0 and $alignment_mode == 4)) and (scalar(@splitin) >= 6 && $splitin[5] eq '-' || (($feature =~ /\.gtf$/ or $feature =~ /\.gff$/) && $splitin[6] eq '-'))){
					@data_array = reverse(@data_array);
				}
				
				if($heatmap_normalization_by_mean ne ""){
					my @new_data_array = &normalize($wig, @data_array);
					@data_array = @new_data_array;
				}
				
				$tmp=join "\t",@data_array;
				if($alignment_mode == 5){
					my $center=int(scalar(@data_array)/2);
					my $closest="";
					#print "$tmp\n$center\n";
					for(my $dist=0;$dist<$center;$dist++){
						#print "$data_array[$center+$dist]\t$data_array[$center-$dist]\t";
						if($data_array[$center+$dist] ne "n/a"){
							$closest=$data_array[$center+$dist];
							last;
						}elsif($data_array[$center-$dist] ne "n/a"){
							$closest=$data_array[$center-$dist];
							last;
						}
					}
					#print "\n";
					if($closest ne ""){
						$tmp=$closest;
					}else{
						next;
					}
					
				}
				
				#if($motif_freq[$j]==1){
					#$tmp =~ s/[-]\d+[.\d+]/1/g;
				#	$tmp =~ s/n\/a/0/g;
					#$tmp =~ s/n\/a/NA/g;
				#}
				#else{
					$tmp =~ s/n\/a/NA/g;
				#}
				
				$tmp = "$chr\t$cor1\t$cor2\t$strand\t$tmp\n";
				print OUT $tmp;
				
			}
			elsif($include_no_data_line ne ""){
				
				#if($motif_freq[$j]==1){
				#	$tmp = 0;
				#}
				#else{
					$tmp = "NA";
				#}
				$tmp = "$chr\t$cor1\t$cor2\t$strand\t$tmp\n";
				print OUT $tmp;	
			}
			$line++;
			if($line % 1000 == 0){
				print STDERR "Aligning $wig to line $line ... ...\t $total_line genomic features in total\n";
			}
		}
		
		
		if($location_after_adjust_center ne ""){
			close(LOC);	
		}
		close(FH);
		close(OUT);
		
	
	
}

sub normalize{
	my $wig = shift @_;
	my @data_array = @_;
	my $cmd="bigWigInfo $wig |";
	my $mean="";
	open(P, $cmd) or die;
	while (<P>) {
  			chomp;
  			if($_=~/mean:\s*(\S+)$/){
  				$mean = $1;
  			}
	}
	close(P);
	my @out=();
	if($mean ne ""){
		foreach my $input(@data_array){
			if($input ne "n/a" and $mean != 0){
				push(@out, $input/$mean);
			}else{
				push(@out, "n/a");
			}
			
		}
		return @out;
	}else{
		return @data_array;
	}
}

sub alignment_extend{
	my $wig = shift @_;
	my $chr = shift @_;
	my $start = shift @_;
	my $end = shift @_;
	my $len = shift @_;
	my $span = shift @_;
	my @data_array = @_;
	my @new_array=();
	##left_extend
	my $left=$start-$len;
	my $right=$start+1;
	my $align_cmd="$ucsc_script $wig $chr $left $right $span |";
	my $tmp_left="";
	open(P, $align_cmd) or die;
	while (<P>) {
  		chomp;
  		$tmp_left = $_;
	}
	close(P);
	if($tmp_left ne ""){
			@new_array=(split("\t",$tmp_left),@data_array);
	}else{
		@new_array=@data_array;
		for ((1..$span)){
			unshift(@new_array,"n/a");
		}
	}
	
	##right_extend
	$left=$end-1;
	$right=$end+$len;
	$align_cmd="$ucsc_script $wig $chr $left $right $span |";
	my $tmp_right="";
	open(P, $align_cmd) or die;
	while (<P>) {
  		chomp;
  		$tmp_right = $_;
	}
	close(P);
	if($tmp_right ne ""){
			@new_array=(@new_array, split("\t",$tmp_right));
	}else{
		for ((1..$span)){
			push(@new_array,"n/a");
		}
	}	
	return @new_array;
}



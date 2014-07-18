#!/usr/bin/perl -w
## pipeline script for MethyPattern align in to genome feature
## author: Yaping Liu  lyping1986@gmail.com

use strict;
use Getopt::Long;




sub usage {

    print "\nUsage:\n";
    print "perl MethyPatternAlignEasyUsage.pl [Options] BISSNP_JAR DATA_FILE FEATURE_BED_FILE REF_GENOME SAMPLE_NAME DBSNP_VCF\n\n";

    print " Provide an easy way to operate MethyPatternAlignment in one command line\n";
    print " It would automately generate: \n";
    print " 1) Data matrix around feature, which could be used for heatmap making next step(need to enable --smooth for clustering)\n";
    print " 2) Summary plot of data pattern around feature, .eps file\n";
    print " 3) Data point density around feature, which could bed used for chip-seq or other tag density heatmap plot\n\n";
    

    
    print " All the scripts should be put under the same directory. It requires R installed in your machine\n";

	print "  BISSNP_JAR     Bis-SNP .jar location and name.\n";
    print "  DATA_FILE     input BAM/Bed files to do alignment. For Bed file, it should contain Methylation value, otherwise it will only plot the tag density.\n";
    print "  FEATURE_BED_FILE    Feature .bed file which you want to align to. e.g. Chip-seq peaks, TSS.\n";
    print "  REF_GENOME    .fa/.fasta file, the reference genome file which should be the \n";
    print "              same when you do reads mapping to obtain BAM files.\n";
    print "  SAMPLE_NAME     Sample's name you want to appear in the prefix of data matrix and plot header.\n";
    print "  DBSNP_VCF     dbSNP.vcf file used for genotyping when BAM file rather than Bed file is provided. Its chromosome order should be the same as REF_GENOME.\n\n";
	
	print "  [Options]:\n\n";
	print "  --cytosines STRING : cytosine pattern to be aligned (Default: --cytosine CG,1)\n";
	print "                      The default option will only align CpG methylation pattern around genomic feature provided \n";
	print "                      NOMe-seq mode will change the default value to '--cytosine GCH,2 --cytosine HCG,2'\n\n";
	print "  --pbs : when specified, the job will be submitted to hpcc pbs job\n\n";
	print "  --data_point_density : when specified, it will align the data point density(or chip-seq tag density) instead of methylation value, only support to align to bed yet\n\n";
	print "  --avoid_align_itself : when specified and also --data_point_density specified, it will not align the the feature which is the same as themself\n\n";
	print "  --shrink_dp_size : when specified, and in data_point_density mode, it will shrink the input data block size to its center point\n\n";
	print "  --omit_align_step : when specified, it will not do alignment but just do plot step\n\n";
	print "  --omit_plot_step : when specified and in align to bam files, it will not do R plot step but will submit two jobs at the same time for align GCH and HCG to accelerate the speed.\n\n";
	print "  --mem NUM : How much gigabytes memory required. (Default: 8)\n\n";
	print "  --cpu NUM : How many cpus required. (Default: 6)\n\n";
	print "  --r :specify the path to use R (Default: R)\n\n";
    print "  --nomeseq : enable NOMe-seq mode, which will draw two lines (GCH/HCG methylation) in the same plot\n\n";
    print "  --heatmap : generate matrix file for Cluster3 plotting heatmap\n\n";
    print "  --r_script FILE :  R script's path and name, which is used for plot\n\n";
    print "  --sort_perl_script FILE :  sortByRefAndCor.pl script's path and name, which is used for sorting input feature bed file\n\n";
    print "  --result_dir DIR :   directory to store alignment matrix results and plots\n\n";
    print "  --qual NUM : genotyping quality score (Default: 20)\n\n";
    print "  --smooth : enable smooth around bin_size region, e.g. each data point represent the avergae of methylation value around 20bp\n\n";
    print "  --use_bad_mate : enable the use of unproper paired reads\n\n";
    print "  --invert_dups_usage : when met inverted dups reads (Default: --invert_dups_usage 1)\n";
    print "                         1) USE_ONLY_1ST_END: only use the 1st end of reads  when it was a inverted dups\n";
	print "                         2) USE_BOTH_END: use both ends of reads  when it was a inverted dups\n";
	print "                         3) NOT_TO_USE: don't use inverted dups at all\n\n";
    
    print "  --data_matrix_scale NUM :   the maximum distance to the feature center/alignment start loci. e.g. 2000 means +/- 2kb around\n";
    print "                              (Default: plot_x_axis_scale + 10*bin_size)\n\n";
    print "  --bin_size NUM : bin size(bp) to average for each data point. (Default: 20)\n\n";
    print "  --auto_scale : automately scale the axis by the maximum in this data. (Default: not enabled)\n\n";
    print "  --plot_x_axis_scale NUM : maximum x axis distance to plot in the summary plot. (Default: 1000)\n\n";
 	print "  --plot_x_axis_step NUM : step size of x axis to plot in the summary plot (Default: 250)\n\n";
 	print "  --plot_y_axis_max NUM : maximum y axis value to plot in the summary plot (Default: 1.0)\n\n";
 	print "  --plot_y_axis_min NUM : maximum y axis value to plot in the summary plot (Default: 0.0)\n\n";
 	print "  --plot_y_axis_step NUM : step size of x axis to plot in the summary plot (Default: 0.2)\n\n";
 	print "  --minCT NUM : minimum number of C or T reads to call methylation (Default: 1)\n\n";
 	print "  --mmq NUM : minimum mapping quality score (Default: 30)\n\n";
 	print "  --mbq NUM : minimum base quality score (Default: 0)\n\n";
 	print "  --trim5 NUM : trim 5 end of sequence (Default: 0)\n\n";
 	print "  --trim3 NUM : trim 3 end of sequence (Default: 0)\n\n";
 	print "  --color STRING : color of line drawing in the summary plot (Default: black)\n\n";
 	print "  --minConv NUM : minConv filter to filter out first unconverted cytosines. (Default: 1)\n\n";
 	print "  --alignment_mode NUM : (Default: --alignment_mode 1)\n";
 	print "                         1) Align to the feature center\n";
	print "                         2) Align to the feature 5' boundary\n";
	print "                         3) Align to the feature 3' boundary\n\n";
 	print "  --bed_format_style NUM : (Default: --bed_format_style 2)\n";
 	print "                         1) BED6PLUS2: TCGA bed 6+2 format, whose score column is 0-1000 to represent methylation value\n";
	print "                         2) BED3PLUS2: bed 3+2 format, whose score column is 0-100 to represent methylation value\n";
	print "                         3) ENCODE: Encode bed 9 format, whose score column is 0-1000 to represent methylation value\n\n";
    exit(1);
}

my $bistools_path=`echo \$BISTOOLS`;
chomp($bistools_path);

##default option setting
my $pbs = "";
my $r="R";
my $mem = 8;
my $cpu = 6;
my $data_point_density = "";
my $avoid_align_itself = "";
my $shrink_dp_size = "";
my $omit_align_step = "";
my $omit_plot_step = "";
my $nomeseq = "";
my $heatmap = "";
my $r_script = "$bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotSinglePlotNoSumFeature.R";
my $sort_perl_script = "$bistools_path/Bis-QC/after_reads_mapping/sortByRefAndCor.pl";
my $result_dir = "./";
my $smooth = "";
my $invert_dups_usage = 1;
my $use_bad_mate = "";
my $minConv=1;

my $qual = 20;
my $bin_size = 20;
my $plot_x_axis_scale = 1000;
my $auto_scale = "";
my $plot_x_axis_step = 250;
my $plot_y_axis_max = 100;
my $plot_y_axis_min = 0.0;
my $plot_y_axis_step = 20;
my $minCT = 1;
my $mmq = 30;
my $mbq = 0;
my $trim5 = 0;
my $trim3 = 0;
my $color = "black";
my $alignment_mode = 1;
my $bed_format_style = 2;
my @cytosines = ();
my $data_matrix_scale = $plot_x_axis_scale + 10*$bin_size;
#my $data_matrix_scale = 2000;

GetOptions( "pbs" => \$pbs,
			"r=s" => \$r,
			"mem=i" => \$mem,
			"cpu=i" => \$cpu,
			"data_point_density" => \$data_point_density,
			"avoid_align_itself" => \$avoid_align_itself,
			"shrink_dp_size" => \$shrink_dp_size,
			"omit_align_step" => \$omit_align_step,
			"omit_plot_step" => \$omit_plot_step,
			"nomeseq" => \$nomeseq,
			"heatmap" => \$heatmap,
			"r_script=s" => \$r_script,
			"sort_perl_script=s" => \$sort_perl_script,
			"result_dir=s" => \$result_dir,
			"smooth" => \$smooth,
			"data_matrix_scale=i" => \$data_matrix_scale,
			"bin_size=i" => \$bin_size,
			"qual=i" => \$qual,
			"plot_x_axis_scale=i" => \$plot_x_axis_scale,
			"minCT=i" => \$minCT,
			"use_bad_mate" => \$use_bad_mate,
			"invert_dups_usage" =>\$invert_dups_usage,
			"auto_scale" => \$auto_scale,
			"plot_x_axis_step=i" => \$plot_x_axis_step,
			"plot_y_axis_max=f" => \$plot_y_axis_max,
			"plot_y_axis_min=f" => \$plot_y_axis_min,
			"plot_y_axis_step=f" => \$plot_y_axis_step,
			"color=s" => \$color,
			"minConv=i" => \$minConv,
			"mmq=i" => \$mmq,
			"mbq=i" => \$mbq,
			"trim5=i" => \$trim5,
			"trim3=i" => \$trim3,
			"alignment_mode=i" => \$alignment_mode,
			"cytosines=s" => \@cytosines,
			"bed_format_style=i" => \$bed_format_style);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) < 5 ) {
    print "Wrong number of arguments\n";
    usage();
}

my $bissnp_jar = $ARGV[0];
my $input_file = $ARGV[1];
my $feature = $ARGV[2];
my $ref = $ARGV[3];
my $sample = $ARGV[4];
my $dbsnp = $ARGV[5];

my %alignment_modes_hash = (
	1 => "Center",
	2 => "FiveEnd",
	3 => "ThreeEnd",
);

my %bed_formats_hash = (
	1 => "BED6PLUS2",
	2 => "BED3PLUS2",
	3 => "ENCODE",
);

my %invert_dups_hash = (
	1 => "USE_ONLY_1ST_END",
	2 => "USE_BOTH_END",
	3 => "NOT_TO_USE",
);

if(scalar(@cytosines) == 0){
	@cytosines = ("CG,1");
}
if($nomeseq ne ""){
	@cytosines = ("GCH,2", "HCG,2");
#	if($r_script eq "$bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotSinglePlot.R" || $r_script eq "$bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotSinglePlot.R"){
		$r_script = "$bistools_path/Bis-QC/after_reads_mapping/MethyPatternFeaturePlotForNOMeSeq.R";
#	}
}

$mem .= "G";

my $legend_name="Methylation";

my $sort_feature_file = &sort_bed_file($feature);
my $location_file = &generate_location_file($sort_feature_file); ## generate location file with 1000bp more distance than $data_matrix_scale specified

if ( $input_file =~ /\.bam/) {
    if(scalar(@ARGV) < 6 ){
    	print "Wrong number of arguments\n";
    	usage();
    }

    &align_bam();
    
}
elsif( $input_file =~ /\.bed/){
	my $sorted_input_file = $input_file;
	if($input_file !~ /sort/){
		$sorted_input_file = &sort_bed_file($input_file);
	}

	&align_bed($sorted_input_file);
}
elsif( $input_file =~ /\.gtf$/){
	$input_file = &gtf_to_bed($input_file);
	my $sorted_input_file = $input_file;
	#if($input_file !~ /sort/){
		$sorted_input_file = &sort_bed_file($input_file);
	#}

	&align_bed($sorted_input_file);
}
else{
		print "Only allow bed/gtf format or BAM format as data file currently!!\n\n";
    	exit(1);
}



sub align_bam{
	my $java_cmd;
	my $prefix = $input_file;
	my $java_cmd_hcg;
	$prefix =~ s/.+\///g;
	$prefix =~ s/^(\w+)\.\S+.bam/$1/;
	$prefix = $sample.".$prefix";
	
	my $R_cmd_plot = "";
	my $gch_file;
	my $cpg_file;
	my $hcg_file;
	if($smooth eq ""){
		$smooth = 0;
	}
	else{
		$smooth = 1;
	}
	my $sort_feature_file_notpath = $sort_feature_file;
	$sort_feature_file_notpath =~ s/.+\///g;
	$sort_feature_file_notpath =~ s/minusUpstream\S+.plusDownstream\S+\.bed//;
	if($nomeseq ne ""){
		 $java_cmd = "java -Xmx$mem -jar $bissnp_jar -T MethyPatternFeatureByBam ";
		$java_cmd .= "-R $ref -I $input_file -D $dbsnp -stand_call_conf $qual -stand_emit_conf 0 -orientated -mmq $mmq -mbq $mbq -trim5 $trim5 -trim3 $trim3 -minConv $minConv ";
		$java_cmd .= "-feature $sort_feature_file -distance $data_matrix_scale -minCTdepth $minCT -alignmentType $alignment_modes_hash{$alignment_mode} ";
		$java_cmd .= "-L $location_file -invDups $invert_dups_hash{$invert_dups_usage} ";
		if($use_bad_mate ne ""){
			$java_cmd .= "-badMate ";
		}
		
		$gch_file = $result_dir."$prefix.alignedTo.$sort_feature_file_notpath.gch.$data_matrix_scale.txt";
		$java_cmd = $java_cmd."-methyFile $gch_file ";
		my $gch_dataPoint_file = $result_dir."$prefix.alignedTo.$sort_feature_file_notpath.gchDataPoint.$data_matrix_scale.txt";
		$java_cmd = $java_cmd."-dataPoint $gch_dataPoint_file -C GCH,2 \n";
		if($heatmap ne ""){
			$java_cmd = $java_cmd."perl /home/uec-00/yapingli/code/mytools/perl/generate_heatmap_matrix.pl --smooth $smooth --data_matrix_scale $data_matrix_scale --bin_size $bin_size $gch_file\n";
		}
		
		if($omit_plot_step eq ""){
			$java_cmd .= "java -Xmx$mem -jar $bissnp_jar -T MethyPatternFeatureByBam ";
			$java_cmd .= "-R $ref -I $input_file -D $dbsnp -stand_call_conf $qual -stand_emit_conf 0 -orientated -mmq $mmq -mbq $mbq -trim5 $trim5 -trim3 $trim3 -minConv $minConv ";
			$java_cmd .= "-feature $sort_feature_file -distance $data_matrix_scale -minCTdepth $minCT -alignmentType $alignment_modes_hash{$alignment_mode} ";
			$java_cmd .= "-L $location_file -invDups $invert_dups_hash{$invert_dups_usage} ";
			if($use_bad_mate ne ""){
				$java_cmd .= "-badMate ";
			}
			$hcg_file = $result_dir."$prefix.alignedTo.$sort_feature_file_notpath.hcg.$data_matrix_scale.txt";
			$java_cmd = $java_cmd."-methyFile $hcg_file -C HCG,2 \n";
			if($heatmap ne ""){
				$java_cmd = $java_cmd."perl /home/uec-00/yapingli/code/mytools/perl/generate_heatmap_matrix.pl --smooth $smooth --data_matrix_scale $data_matrix_scale --bin_size $bin_size $hcg_file\n";
			}
			#my $wcg_file = $result_dir."$prefix.alignedTo.$sort_feature_file.wcg.$data_matrix_scale.txt";
			#$java_cmd = $java_cmd."-wcgFile $wcg_file ";
		}
		else{
			$java_cmd_hcg = "java -Xmx$mem -jar $bissnp_jar -T MethyPatternFeatureByBam ";
			$java_cmd_hcg .= "-R $ref -I $input_file -D $dbsnp -stand_call_conf $qual -stand_emit_conf 0 -orientated -mmq $mmq -mbq $mbq -minConv $minConv ";
			$java_cmd_hcg .= "-feature $sort_feature_file -distance $data_matrix_scale -minCTdepth $minCT -alignmentType $alignment_modes_hash{$alignment_mode} ";
			$java_cmd_hcg .= "-L $location_file -invDups $invert_dups_hash{$invert_dups_usage} ";
			if($use_bad_mate ne ""){
				$java_cmd_hcg .= "-badMate ";
			}
			$hcg_file = $result_dir."$prefix.alignedTo.$sort_feature_file_notpath.hcg.$data_matrix_scale.txt";
			$java_cmd_hcg = $java_cmd_hcg."-methyFile $hcg_file -C HCG,2 \n";
			if($heatmap ne ""){
				$java_cmd_hcg = $java_cmd_hcg."perl /home/uec-00/yapingli/code/mytools/perl/generate_heatmap_matrix.pl --smooth $smooth --data_matrix_scale $data_matrix_scale --bin_size $bin_size $hcg_file\n";
			}
		}
		
		
		$R_cmd_plot = "$r --no-restore --no-save --args wd=$result_dir prefix=$prefix gchfn=$gch_file hcgfn=$hcg_file step=$bin_size scale=$plot_x_axis_scale axistep=$plot_x_axis_step smooth=$smooth y_scale_max=$plot_y_axis_max y_scale_min=$plot_y_axis_min y_step=$plot_y_axis_step < $r_script \n";
	
	}
	else{
		foreach my $cytosine(@cytosines){
			my $output_format = $cytosine;
			$output_format =~ s/\,\d+//;
			$java_cmd .= "java -Xmx$mem -jar $bissnp_jar -T MethyPatternFeatureByBam ";
			$java_cmd .= "-R $ref -I $input_file -D $dbsnp -stand_call_conf $qual -stand_emit_conf 0 -orientated -mmq $mmq -mbq $mbq -trim5 $trim5 -trim3 $trim3 -minConv $minConv ";
			$java_cmd .= "-feature $sort_feature_file -distance $data_matrix_scale -minCTdepth $minCT -alignmentType $alignment_modes_hash{$alignment_mode} ";
			$java_cmd .= "-C $cytosine ";
			$java_cmd .= "-L $location_file -invDups $invert_dups_hash{$invert_dups_usage} ";
			if($use_bad_mate ne ""){
				$java_cmd .= "-badMate ";
			}
			$cpg_file = $result_dir."$prefix.alignedTo.$sort_feature_file_notpath.$output_format.$data_matrix_scale.txt";
			$java_cmd = $java_cmd."-methyFile $cpg_file ";
			my $data_point_file = $result_dir."$prefix.alignedTo.$sort_feature_file_notpath.$output_format".".DataPoint.$data_matrix_scale.txt";
			$java_cmd = $java_cmd."-dataPoint $data_point_file \n";
			if($heatmap ne ""){
				$java_cmd = $java_cmd."perl /home/uec-00/yapingli/code/mytools/perl/generate_heatmap_matrix.pl --smooth $smooth --data_matrix_scale $data_matrix_scale --bin_size $bin_size $cpg_file\n";
			}
			$R_cmd_plot .= "$r --no-restore --no-save --args wd=$result_dir prefix=$prefix gchfn=$cpg_file step=$bin_size scale=$plot_x_axis_scale axistep=$plot_x_axis_step smooth=$smooth y_scale_max=$plot_y_axis_max y_scale_min=$plot_y_axis_min y_step=$plot_y_axis_step color=$color < $r_script \n";
			
		}
		
	}
			
	

	if($pbs ne ""){
		my $header = "#PBS -q laird\n";
		$header .= "#PBS -l walltime=168:00:00,mem=$mem,nodes=1:ppn=$cpu\n";
		#$header .= "#PBS -W depend=afterany:1996998.hpc-pbs\n";
		$header .= "#PBS -N $prefix\n";
		$header .= "#PBS -j oe\n";
		$header .= "cd \$PBS_O_WORKDIR\n";
		my $header_2 = $header;
		if($omit_plot_step ne ""){
				$header_2 .= "$java_cmd_hcg\n";
		}

		if($omit_align_step eq ""){
			
			$header .= "$java_cmd\n";
		}
		
		$header .= "$R_cmd_plot \n";	
		my $count = int(rand(10000));
		my $outfileName = "submitPBS" . $count;
			open(OUT, ">$outfileName") or die "can not open file:$!";
			print OUT $header;
			close OUT;
			system("qsub $outfileName");
			$count++;
			print "$prefix is submitted\n";
			if($omit_plot_step ne ""){
				open(OUT, ">$outfileName") or die "can not open file:$!";
				print OUT $header_2;
				close OUT;
				system("qsub $outfileName");
				$count++;
				print "$prefix is submitted\n";	
			}
			unlink $outfileName;
			#print "$header\n";
	}
	else{
		if($omit_align_step eq ""){
			system("$java_cmd");
		}
		system("$R_cmd_plot");
	}

	
}


sub align_bed{
	my $sorted_input_file = shift(@_);
	my $count_cmd="wc -l $sorted_input_file\n";
	my $java_cmd = "java -Xmx$mem -jar $bissnp_jar -T MethyPatternFeatureByBed ";
	$java_cmd .= "-R $ref -orientated ";
	$java_cmd .= "-feature $sort_feature_file -distance $data_matrix_scale -minCTdepth $minCT -alignmentType $alignment_modes_hash{$alignment_mode} ";
	$java_cmd .= "-L $location_file -bedFormat $bed_formats_hash{$bed_format_style} ";	
	my $prefix = $sorted_input_file;
	$prefix =~ s/.+\///g;
	$prefix =~ s/\.sort\.bed//;
	$prefix = $sample.".$prefix";
	my $gch_file;
	$sort_feature_file =~ s/.+\///g;
	$sort_feature_file =~ s/minusUpstream\S+.plusDownstream\S+\.bed//;
	my $shrink_input_file;
	if($data_point_density ne ""){
		$legend_name = "Density";
		if($shrink_dp_size ne ""){
			$shrink_input_file = &shrink_file($sorted_input_file);
			$java_cmd .= "-values $shrink_input_file ";
			$prefix = $shrink_input_file;
			$prefix =~ s/.+\///g;
			$prefix = $sample.".$prefix";
		}
		else{
			$java_cmd .= "-values $sorted_input_file ";
		}
		$gch_file = $result_dir."$prefix.alignedTo.$sort_feature_file.dataPoint.$data_matrix_scale.txt";
		if($avoid_align_itself ne ""){
			$java_cmd = $java_cmd."-noItself ";
		}
		
		$java_cmd = $java_cmd."-dataPoint $gch_file \n";
	}
	else{
		$java_cmd .= "-values $sorted_input_file ";
		$gch_file = $result_dir."$prefix.alignedTo.$sort_feature_file.methy.$data_matrix_scale.txt";
		$java_cmd = $java_cmd."-outFile $gch_file \n";
	}
	
	if($data_point_density ne "" and $shrink_dp_size ne ""){
		$java_cmd .= "rm $shrink_input_file\n";  
	}
	
	if($smooth eq ""){
		$smooth = 0;
	}
	else{
		$smooth = 1;
	}
	
	my $numFeature = `$count_cmd`;
	$numFeature =~ s/(\d+)\s+\S+/$1/;
	chomp($numFeature);
	my $sumFeatureFn = $gch_file.".sumFeature.txt";
	my $perl_cmd = "perl /home/uec-00/yapingli/code/mytools/perl/sumFeature.pl $gch_file $plot_x_axis_scale \n";
	if($auto_scale ne ""){
		$auto_scale="TRUE";
	}
	else{
		$auto_scale="FALSE";
	}
	my $R_cmd_plot = "$r --no-restore --no-save --args wd=$result_dir prefix=$prefix gchfn=$gch_file step=$bin_size scale=$plot_x_axis_scale axistep=$plot_x_axis_step smooth=$smooth y_scale_max=$plot_y_axis_max y_scale_min=$plot_y_axis_min y_step=$plot_y_axis_step color=$color sumFeatureFn=$sumFeatureFn numFeature=$numFeature legendName=$legend_name autoScale=$auto_scale < $r_script \n";
	
			
	if($pbs ne ""){
		my $header = "#PBS -q laird\n";
		$header .= "#PBS -l walltime=168:00:00,mem=$mem,nodes=1:ppn=$cpu\n";
		#$header .= "#PBS -W depend=afterany:1996998.hpc-pbs\n";
		$header .= "#PBS -N $prefix\n";
		$header .= "#PBS -j oe\n";
		$header .= "cd \$PBS_O_WORKDIR\n";
		if($omit_align_step eq ""){
			$header .= "$java_cmd\n";
		}
		$header .= "$perl_cmd\n";
		$header .= "$R_cmd_plot \n";	
		my $count = int(rand(10000));
		my $outfileName = "submitPBS" . $count;
			open(OUT, ">$outfileName") or die "can not open file:$!";
			print OUT $header;
			close OUT;
			system("qsub $outfileName");
			$count++;
			print "$prefix is submitted\n";
			unlink $outfileName;
			#print "$header\n";
	}
	else{
		if($omit_align_step eq ""){
			system("$java_cmd");
		}
		system("$perl_cmd");
		system("$R_cmd_plot");
	}			
	#&generate_heatmap_matrix($gch_file);		
}

sub sort_bed_file{
	my $file = shift(@_);
	if($file =~ /sort/){
		return $file;
	}
	my $sorted_file = $file;
	$sorted_file =~ s/\.bed//;
	$sorted_file .= ".sort.bed";
	unless (-e $sort_perl_script) {
 		print "\n $sort_perl_script file not exists!\n\n";
 		exit(1);
	}
	my $header .= "perl $sort_perl_script --k 1 --c 2 ";
	$header .= "--tmp ./ ";
	$header .= "$file ";
	$header .= "$ref";
	$header .= ".fai ";
	$header .= "> $sorted_file \n";
	print "$header\n";
	system($header);
	return $sorted_file;
}


sub generate_location_file{
	my $file = shift(@_);
	my $ref_index=$ref.".fai";
	open(RI,"<$ref_index") or die "cant find reference genome index file!!\n";
	my @refi=<RI>;
	chomp(@refi);
	close(RI);
	my %genome_size=();
	foreach my $c(@refi){
		my @splitin=split "\t",$c;
		my $key=$splitin[0];
		my $value=$splitin[1];
		$genome_size{$key}=$value;
	}
	
	my $extend_dist = $data_matrix_scale + 2*$bin_size;
	my $location = $file;
	$location =~ s/\.bed//;
	$location .= ".minusUpstream$extend_dist".".plusDownstream$extend_dist".".sort.bed";
	open(FH,"<$file") or die "can not open file:$!";
	open(OUT,">$location") or die "can not open file:$!";
	while(<FH>){
		chomp;
		my $line = $_;
		my @splitin = split "\t",$line;
		my $start = int(($splitin[1] + $splitin[2])/2 -1 - $extend_dist);
		$start = $start < 0 ? 0 : $start;
		my $end = int(($splitin[1] + $splitin[2])/2 + $extend_dist);
		$end = $end > $genome_size{$splitin[0]} ? $genome_size{$splitin[0]} : $end;
		print OUT "$splitin[0]\t$start\t$end\n";
	}
	close(FH);
	close(OUT);
	return $location;
}

sub generate_heatmap_matrix{
	my $file = shift(@_);
	open(FH,"<$file") or die "can not open file:$!";
	my $location = $file;
	$location =~ s/\.txt//;
	$location .= ".heatmap.matrixForCluster3.txt";
	open(OUT,">$location") or die "can not open file:$!";
	print OUT "Genes\t";
	for(my $i=0-$data_matrix_scale;$i<=$data_matrix_scale; $i++){
		if($i==$data_matrix_scale){
			print OUT "$i";
		}
		else{
			print OUT "$i\t";
		}
		
	}
	print OUT "\n";
	while(<FH>){
		chomp;
		my $line = $_;
		my @splitin = split "\t",$line;
		#print "$#splitin\t$data_matrix_scale\n";
		if($#splitin == $data_matrix_scale * 2 + 4){
			my $header = join ":",$splitin[0],$splitin[1],$splitin[2],$splitin[3];
			print OUT "$header\t";
			for(my $i=4;$i<$data_matrix_scale * 2+4; $i++){
				if($smooth == 1){
					my $sum=0;
					my $num=0;
					for(my $j=$i-int($bin_size/2);$j<=$i+int($bin_size/2); $j++){
						if($j>=4 and $j < ($data_matrix_scale * 2+4)){
							if($splitin[$j] ne "NaN"){
								$sum += $splitin[$j];
								$num++;
							}
						
						}
					
					}
					if($num == 0){
						print OUT "NA\t";
					}
					else{
						my $methy = $sum/$num;
						print OUT "$methy\t";
					}
					
				}
				else{
					if($splitin[$i] eq "NaN"){
						print OUT "NA\t";
					}
					else{
						print OUT "$splitin[$i]\t";
					}
					
				}
				
			}
			print OUT "\n";

		}
		
	}
	close(FH);
	close(OUT);
}

sub shrink_file{
	my $input = shift @_;
	my $output = $input;
	$output =~ s/\.bed//;
	$output .= ".center.bed";
	open(FH,"<$input") or die;
	open(OUT,">$output") or die;
	while(<FH>){
		chomp;
		my @splitin=split "\t";
		my $start = int(($splitin[1] + $splitin[2])/2);
		my $end = int(($splitin[1] + $splitin[2])/2)+1;
		$splitin[1] = $start;
		$splitin[2] = $end;
		my $line = join "\t",@splitin;
		print OUT "$line\n";
	}
	close(FH);
	close(OUT);
	return $output;
}


sub gtf_to_bed{
	my $input = shift @_;
	my $output = $input;
	$output =~ s/\.gtf/.bed/;
	open(FH,"<$input") or die;
	open(OUT,">$output") or die;
	while(<FH>){
		chomp;
		my $line=$_;
		my @splitin = split "\t",$line;
		next if($splitin[0]=~/\_/ or $splitin[0]!~/^chr/);
		my $start = $splitin[3]-1;
		my $end = $splitin[4];
		print OUT "$splitin[0]\t$start\t$end\t$splitin[1]\t$splitin[5]\t$splitin[6]\n";
	}

	close(FH);
	close(OUT);
	return $output;
}


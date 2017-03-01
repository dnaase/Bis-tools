#!/usr/bin/perl -w
##This script convert bam file into visualized wig file, and then to binary wiggler (output_prefix.bw) file.
#it is an improved version of wrap_wiggler.pl, since it could use wiggler for normalization, or just normalize count by user provided.
#if do read extesion, it need to have fragment length estimation file from HOMER

##bam to bed
##homer estimate fragment length
##mode0: no normlization step, [reads extension,  bed to bedgraph, bedgraph get equal binned ]
##mode1: use wiggler, do read extension, normlization
##mode2: reads extension,  bed to bedgraph, bedgraph get equal binned, calc RPKM, [mius Input.RPKM] [transform score to Z-score]
##mode3: reads extension,  bed to bedgraph, bedgraph get equal binned, cal ratio to adjacent 10kb region's reads mean number (not implement it yet....)
##mode4: reads extension,  bed to bedgraph, bedgraph get equal binned, quantile normalization
##bedgraph to bigwig


##example cmd:

##mode0:
##mode1:
##mode2:
##mode3:
##mode4:perl ~/storage/code/mytools/perl/bam2normalizedwig.pl test1.extend.sort.bed test2.extend.sort.bed --no_bam2bed --no_read_extension --no_bed_sort --keep_bed --homer ./ --keep_bedgraph --normalize 5 --bin 10


## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-5-10

#Usage: perl bam2normalizedwig.pl [option] input.bam/input.bed

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl bam2normalizedwig.pl [option] input.bam/input.bed \n\n";
    print STDERR "Convert input bam/bed/bedGraph/wig file into bw.\n\n";
    print STDERR "Fully function requirs installing Wiggler, HOMER, UCSC's bedGraphToBigWig, genome's index file.\n\n";
    print STDERR "  [Options]:\n\n";
	print STDERR "  --frag_len NUM : fragment length to extend the reads, allow multiple number, but the order should be the same as input bam or bed files. (Default: no specified --frag_len, no extension)\n\n";
	print STDERR "  --homer DIR: when specified, used homer to do fragment length estimation automately, and HOMER tag directory will be in DIR\n\n";
	print STDERR "  --paired : enable it for paired-end reads, only use proper paired reads\n\n";
	print STDERR "  --rnaseq : enable it for RNA-seq reads, use -split option in bedtools to avoid false coverage caused by CIGAR\n\n";
	print STDERR "  --normalize NUM : the way to do normalization, require input file to be bed/bam file, bedgraph file could not be normalized correctly:\n";
	print STDERR "  					1) use total mapped reads count for normalization, just calculate RPKM out.. \n";
	print STDERR "  					2) no normalization, just raw coverage\n";
	print STDERR "  					3) use wiggler to normalize. good for ChIP-seq, MNase-seq, DNase-seq\n";
 	print STDERR "  					4) use 10kb region reads coverage in the upstream/downstream 20kb(as Ben did..)\n\n";
  	print STDERR "  					5) use quantile normalization, need to provide multiple files at the same time, it will normalize them all together..\n\n";
	print STDERR "  					6) use total mapped reads count for normalization, just calculate RPM/bp out..(Default)\n";
  	
  			
	print STDERR "  --transform NUM   	1) use z-score transformation(as Xie 2013 Cell for ChIP-seq) (Default)\n";  
 	print STDERR "  					2) use Variance Stabilization transformation(as DEseq did for RNA-seq)\n"; 
  	print STDERR "  					3) no transform \n";
  	print STDERR "  					4) use local z-score transformation (local 1Mb region's mean and variance, instead of global mean and variance)\n\n";
   	print STDERR "  --input_bedgraph FILE : RPKM in the bed file will minus the value in input .bedgraph file. this input bedgraph should be normalized RPKM \n\n";
 	print STDERR "  --bin NUM: split bedgraph into NUM bp non-overlapped windows (Default: not binned)\n\n";
 	print STDERR "  --no_bam2bed : when specified, don't do bam to bed step\n\n";
 	print STDERR "  --no_read_extension : when specified, don't do read_extension_step\n\n";
 	print STDERR "  --no_bed_sort : when specified, don't do bed_sort step\n\n";
 	print STDERR "  --no_bed2bedgraph : when specified, don't do bed2bedgraph step\n\n";
 	print STDERR "  --keep_bedgraph : when specified, don't delete the bedgraph file before z-score transform\n\n";
  	print STDERR "  --keep_bed : when specified, don't delete the original input bed file\n\n";
  	print STDERR "  --no_bedgraph2wig : when specified, don't do bedgraph2wig step\n\n";
  	print STDERR "  --useLimma : when specified, use limma package to do quantile normalization\n\n";
  	print STDERR "  --pseudo_count NUM: the number used as a phsedo count for RNAseq's RPKM calculation. (Default: 0)\n\n";
    print STDERR "  --max_read_len NUM: the maximum reads length allowed for chip-seq. specifically designed since wiggler's mappability track is 25-54, while very few reads(100 out of 4M reads) have some gapped aligned, so it is liks more than 54bp long. (Default: 53)\n\n";
  	
 	
    exit(1);
}

my $WIGGLER = "/export/uec-gs1/laird/users/yaping/software/external/align2rawsignal/bin/align2rawsignal";
my $MAPPABILITY = "/export/uec-gs1/laird/users/yaping/data/genome_data/mappability/hg19/globalmap_k20tok54";
my $REFDIR = "/export/uec-gs1/laird/users/yaping/data/genome_data/genome/hg19/ucsc/maleByChr";
my $WIG2BW = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/bedGraphToBigWig";
my $CHROMSIZE="/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/hg19.chrom.sizes";
my $SORTBED = "perl /home/rcf-40/yapingli/storage/code/mytools/perl/sortByRefAndCor.pl";
my $REFINDEX="/export/uec-gs1/laird/users/yaping/data/genome_data/genome/hg19/hg19_rCRSchrm.fa.fai";
my $BEDMINUS = "perl /home/rcf-40/yapingli/storage/code/mytools/perl/bedGraphMinusAnother.pl";
my $BIGWIGSUM = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/bigWigSummary";
my $BEDQUANTNORM = "perl /home/rcf-40/yapingli/storage/code/mytools/perl/bedGraphQuantileNorm.pl";

my @fragment_lens = ();
my $normalize = 6;
my $paired = "";
my $rnaseq = "";
my $homer="";
my $transform=3;
my $input_bedgraph="";
my $no_bam2bed="";
my $no_read_extension="";
my $no_bed_sort="";
my $no_bed2bedgraph="";
my $no_bedgraph2wig="";
my $keep_bedgraph="";
my $keep_bed="";
my $bin="";
my $useLimma="";
my $pseudo_count=0;
my $max_read_len=53;

GetOptions( "frag_len=i" => \@fragment_lens,
			"normalize=i" => \$normalize,
			"transform=i" => \$transform,
			"paired" => \$paired,
			"rnaseq" => \$rnaseq,
			"homer=s" => \$homer,
			"input_bedgraph=s" => \$input_bedgraph,
			"no_bam2bed" => \$no_bam2bed,
			"no_read_extension" => \$no_read_extension,
			"no_bed_sort" => \$no_bed_sort,
			"no_bed2bedgraph" => \$no_bed2bedgraph,
			"no_bedgraph2wig" => \$no_bedgraph2wig,
			"keep_bedgraph" => \$keep_bedgraph,
			"keep_bed" => \$keep_bed,
			"bin=i" => \$bin,
			"useLimma" => \$useLimma,
			"pseudo_count=f" => \$pseudo_count,
			"max_read_len=i" => \$max_read_len,
);


usage() if ( scalar(@ARGV) == 0 );

my @files=@ARGV;
my @bed_file_tmps=@files;



##genome size
my %range=();
open(FH,"<$CHROMSIZE") or die;
while(<FH>){
        chomp;
        my @f=split "\t";
        my $key=$f[0];
        my $value=$f[1];
        $range{$key}=$value;
}
close(FH);

##bam to bed files
for(my $file_num=0;$file_num< scalar(@files);$file_num++){
	
	if($files[$file_num]=~/\.bam$/ and $no_bam2bed eq ""){
		$bed_file_tmps[$file_num]=~s/\.bam$/.bed/;
	
		my $cmd="samtools view -uF 1796 ";
		if( $paired ne ""){
			$cmd="-f 2 "; ##only use proper paired reads..
		}
		$cmd.="$files[$file_num] | bamToBed ";
		if($rnaseq ne ""){
			$cmd.="-split ";
		}
		$cmd.="-i stdin > $bed_file_tmps[$file_num] \n";
		print STDERR $cmd;
		system($cmd)==0 || die;

	}


	
##estimate fragment length by homer
	if($homer ne ""){
		my $homer_out = basename($bed_file_tmps[$file_num]);
		$homer_out =~ s/(\w+)\S+$/$1/;
		$homer_out = $homer."/".$homer_out;
		if(!(-d "$homer_out")){
			my $cmd = "makeTagDirectory $homer_out -genome hg19 $bed_file_tmps[$file_num]\n";
			print STDERR $cmd;
			system($cmd) ==0 || die;
		}
		open(FH,"<$homer_out/tagInfo.txt") or die "no tag info file found\n";

		while(<FH>){
			if($_=~/fragmentLengthEstimate=(\d+)/){
				my $fragment_len=$1;
				push(@fragment_lens,$fragment_len);
			}

		}
		close(FH);
	
	}
	
	##check bed reads length satisfy the maximum read length or not
	if($bed_file_tmps[$file_num]=~/\.bed$/ and $normalize == 3 and $no_bed2bedgraph eq ""){
		my $tmp=$bed_file_tmps[$file_num];
		$tmp=~s/\.bed$/.tmp.bed/;
		open(FH,"<$bed_file_tmps[$file_num]") or die "no bed file found: $bed_file_tmps[$file_num]\n";
		open(OUT,">$tmp") or die "cant write to $tmp\n";
		my $count=0;
		while(<FH>){
			chomp;
			my @splitin=split "\t";
			if($splitin[2]-$splitin[1] > $max_read_len){
				if($splitin[5] eq '-'){
					$splitin[1]=$splitin[2]-$max_read_len;
				}else{
					$splitin[2]=$splitin[1]+$max_read_len;
				}
				
				$count++;
			}
			my $tmpLine=join "\t",@splitin;
			print OUT "$tmpLine\n";
		}
		close(FH);
		close(OUT);
		`mv $tmp $bed_file_tmps[$file_num]`;
		print STDERR "mv $tmp $bed_file_tmps[$file_num]\n";
		print STDERR "File $bed_file_tmps[$file_num] contain $count reads more than $max_read_len bp\n";
	}
		
		
}


##deal with wiggler
my $wig_output="";
if($normalize == 3){
	print STDERR "Normalize reads by Wiggler...\n";
	$ENV{'MCRROOT'} = "/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714/runtime/glnxa64";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714/bin/glnxa644";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714/sys/os/glnxa64";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_thread";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/sys/java/jre/glnxa64/jre/lib/amd64/server";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/sys/java/jre/glnxa64/jre/lib/amd64";
	$ENV{'LD_LIBRARY_PATH'} .= ":/export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714/X11/app-defaults";
	
#`setenv MCRROOT /export/uec-gs1/laird/users/yaping/software/external/MCR2010_ROOT/v714`;
#`setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${MCRROOT}/runtime/glnxa64`;
#`setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${MCRROOT}/bin/glnxa64`;
#`setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${MCRROOT}/sys/os/glnxa64`;
#`setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64/native_threads`;
#`setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64/server`;
#`setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:\${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64`;
#`setenv XAPPLRESDIR \${MCRROOT}/X11/app-defaults`;
 	if(scalar(@fragment_lens)==0){
 		foreach (@bed_file_tmps){
 			push(@fragment_lens,1);
 		}
 	}
	my $cmd = "$WIGGLER -of=bg -n=5 -mm=7 -s=$REFDIR -u=$MAPPABILITY ";
	for(my $i=0;$i<scalar(@bed_file_tmps);$i++){
		my $input = $bed_file_tmps[$i];
		$input =~ s/\.\w+$/.tagalign/;
    	system("ln -s $bed_file_tmps[$i] $input");
		$cmd .= "-i=$input -l=$fragment_lens[$i] ";
	}
	$wig_output=$bed_file_tmps[0];
	if(scalar(@bed_file_tmps)>1){
		$wig_output=~s/\.\w+$/.wiggler_norm.merged.bedgraph/;
	}
	else{
		$wig_output=~s/\.\w+$/.wiggler_norm.bedgraph/;
	}
	 
	$cmd .= "-o=$wig_output\n";
	print STDERR "$cmd\n";

	if(system($cmd)==0 && -e $wig_output){
		for(my $i=0;$i<scalar(@bed_file_tmps);$i++){
			my $input = $bed_file_tmps[$i];
			$input =~ s/\.\w+$/.tagalign/;
    		unlink($input);
		}
		@bed_file_tmps=($wig_output);
		$no_read_extension=1; ##don't do read extension again
		$no_bed_sort=""; ##don;t do bed sort again
		$no_bed2bedgraph=1;
		#$normalize=2; ##don't do normalization again
		#$bin=""; ##dont do bedgraph binned
	}
	@bed_file_tmps=($wig_output);
		$no_read_extension=1; ##don't do read extension again
		$no_bed_sort=""; ##don;t do bed sort again
		$no_bed2bedgraph=1;
}




##read extention


my @bed_files=@bed_file_tmps;
my @bed_file_sorts=@bed_files;
my @bedgraph_files=@bed_file_sorts;
my @bedgraph_file_binneds=@bedgraph_files;


for(my $file_num=0;$file_num< scalar(@bed_file_tmps);$file_num++){
	if(scalar(@fragment_lens)==scalar(@bed_file_tmps) && $no_read_extension eq ""){
		my $fragment_len=$fragment_lens[$file_num];
		print STDERR "read extened to length $fragment_len bp for file: $bed_files[$file_num] ... \n";
		$bed_files[$file_num]=~s/\.\w+$/.extend.bed/;
		open(FH,"<$bed_file_tmps[$file_num]") or die;
		open(OUT,">$bed_files[$file_num]") or die;
		while(<FH>){
        	chomp;
      		 my @f=split "\t";
        	my $chr=$f[0];
        	if($f[5] eq '-'){
                if($f[2]-$fragment_len > 0){
                        $f[1]=$f[2]-$fragment_len;
                }
                else{
                        $f[1]=0;
                }

        	}
        	else{
                if($f[1]+$fragment_len<$range{$chr}){
                        $f[2]=$f[1]+$fragment_len;
                }
                else{
                        $f[2]=$range{$chr};
                }

        	}
        	my $line=join "\t",@f;
        	print OUT "$line\n";
		}
		close(FH);
		close(OUT);
		if( -e $bed_files[$file_num] && $keep_bed eq ""){
			unlink($bed_file_tmps[$file_num]);
		}
	}	
#sort bed file
		$bed_file_sorts[$file_num]=$bed_files[$file_num];	
		if($no_bed_sort eq ""){	
			$bed_file_sorts[$file_num]=~s/\.(\w+)$/.sort.$1/;	
			my $cmd="$SORTBED --k 1 --c 2 --tmp ./ $bed_files[$file_num] $REFINDEX > $bed_file_sorts[$file_num]\n";
			print STDERR $cmd;
			if(system($cmd)==0 && -e $bed_file_sorts[$file_num]){
				unlink($bed_files[$file_num]);
			}
		}
		
##coverage to bedGraph
		$bedgraph_files[$file_num]=$bed_file_sorts[$file_num];
		if($no_bed2bedgraph eq ""){			
			$bedgraph_files[$file_num]=~s/\.\w+$/.bedgraph/;
			my $cmd="genomeCoverageBed -bg -i $bed_file_sorts[$file_num] -g $CHROMSIZE > $bedgraph_files[$file_num]\n";
			print STDERR $cmd;
			if(system($cmd)==0 && -e $bedgraph_files[$file_num]){
			#unlink($bed_file_sorts[$file_num]);
			}
		}
		
##split bedGraph into X bp non-overlap bins, give each bin's reads number
		$bedgraph_file_binneds[$file_num] = $bedgraph_files[$file_num];
		if($bin ne ""){		
			$bedgraph_file_binneds[$file_num] =~ s/\.\w+$/.binned_$bin.bedgraph/;
			my $tmp_bw = $bedgraph_file_binneds[$file_num];
			$tmp_bw =~ s/\.\w+$/.tmp.bw/;
			my $cmd="$WIG2BW $bedgraph_files[$file_num] $CHROMSIZE $tmp_bw\n";
			print STDERR $cmd;
			system($cmd)==0 or die;
			open(OUT,">$bedgraph_file_binneds[$file_num]") or die "cant open file: $bedgraph_file_binneds[$file_num] : $!\n";
			my @chrs=(1..22,"X","Y"); ##currently, only consider autosome+chrX+chrY
			#my @chrs=(1); 
			foreach my $chr(@chrs){
				$chr="chr".$chr;
				my $range_chr = $range{$chr}-1;
				
				
				for(my $bw_start=0; $bw_start<$range_chr; $bw_start+=$bin*100000){ ##when more then 1M data point bigWigSummary program would meet some problem..
						my $bw_end=$bw_start+$bin*100000;
						$bw_end=$bw_end<$range_chr ? $bw_end : $range_chr;
						my $point=int(($bw_end-$bw_start)/10);
						my $tmp =`$BIGWIGSUM $tmp_bw $chr $bw_start $bw_end $point `;
						#print STDERR "$BIGWIGSUM $tmp_bw $chr $bw_start $bw_end $point\n";
						chomp($tmp);

						if($tmp ne ""){
							$tmp =~ s/n\/a/0/g;
								my @data_array=split "\t",$tmp;			
							for(my $start=$bw_start,my $index=0; $start<$bw_end-$bin; $start+=$bin,$index++){
								my $end=$start+$bin;
								print OUT "$chr\t$start\t$end\t$data_array[$index]\n";
							}		
						}
						else{
							for(my $start=$bw_start,my $index=0; $start<$bw_end-$bin; $start+=$bin,$index++){
								my $end=$start+$bin;
								print OUT "$chr\t$start\t$end\t0\n";
							}	
						}
				}
				
			}
			close(OUT);
			#unlink($tmp_bw);
	
		}		

}



##normalized reads
my @bedgraph_file_norms=@bedgraph_file_binneds;
if($normalize == 1 || $normalize == 6){
	for(my $file_num=0;$file_num< scalar(@bedgraph_file_binneds);$file_num++){
		
		#my $total_reads = `awk '{sum+=\$4} END {print sum}' $bedgraph_file`;
		my $total_reads = `wc -l $bed_file_sorts[$file_num]`;
		$total_reads =~ s/(\d+)\s+\S+/$1/;
		chomp($total_reads);
		if($normalize == 6){
			print STDERR "Normalize reads by RPM/bp, the total number of good reads is: $total_reads ...\n";
			$bedgraph_file_norms[$file_num]=~s/\.\w+$/.RPM.bedgraph/;
		}
		else{
			print STDERR "Normalize reads by RPKM, the total number of good reads is: $total_reads ...\n";
			$bedgraph_file_norms[$file_num]=~s/\.\w+$/.RPKM.bedgraph/;
		}
		
		#if( $paired ne ""){
		#	$total_reads*=2;
		#}
	
	
		open(FH,"<$bedgraph_file_binneds[$file_num]") or die "cant open file: $bedgraph_file_binneds[$file_num] : $!\n";
		open(OUT,">$bedgraph_file_norms[$file_num]") or die "cant open file: $bedgraph_file_norms[$file_num] : $!\n";
		while(<FH>){
      	  chomp;
      	  my @f=split "\t";
      	  if($normalize == 6){
      	  	$f[3] = ($f[3]+$pseudo_count)*1000000/$total_reads; # this is euqal to RPM: Reads/($total_reads/1000000)
      	  }
      	  else{
      	  	$f[3] = ($f[3]+$pseudo_count)/(($f[2]-$f[1])/1000); # this is euqal to RPK $f[3]/(length/1000)
      	  	$f[3] = $f[3]*1000000/$total_reads; # this is euqal to RPKM: RPK/($total_reads/1000000)
      	  }
       	 
       	 my $line=join "\t",@f;
      	  print OUT "$line\n";
		}
		close(FH);
		close(OUT);
		if( -e $bedgraph_file_norms[$file_num]){
			#unlink($bedgraph_files[$file_num]);
		}
		
		
##minus RPKM value in input.bedgraph file
	my $bedgraph_file_minus_input=$bedgraph_file_norms[$file_num];
	if($input_bedgraph ne ""){
		print STDERR "minus RPKM in input: $input_bedgraph ...\n";
		$bedgraph_file_minus_input=~s/\.\w+$/.minusInput.bedgraph/;
		my $cmd = "$BEDMINUS --overlap $bedgraph_file_minus_input $bedgraph_file_norms[$file_num] $input_bedgraph \n";
		print STDERR $cmd;
		if(system($cmd)==0 && -e $bedgraph_file_minus_input){
			#unlink($bedgraph_file_norm);
		}
	}


##transform reads
	my $bedgraph_file_transform=$bedgraph_file_minus_input;
	if($transform == 1){
		print STDERR "transform RPKM to Z-score ...\n";	
		my $std_reads = `awk '{delta = \$4 - avg; avg += delta / NR; mean2 += delta * (\$4 - avg); } END { print sqrt(mean2 / NR); }' $bedgraph_file_minus_input`;
		my $mean_reads = `awk '{sum+=\$4} END {print sum/NR}' $bedgraph_file_minus_input`;
		chomp($mean_reads);
		chomp($std_reads);
		print STDERR "RPKM mean is: $mean_reads\nRPKM standard variation is: $std_reads \n ... ...\n";	
		$bedgraph_file_transform=~s/\.\w+$/.Z_score.bedgraph/;
		open(FH,"<$bedgraph_file_minus_input") or die;
		open(OUT,">$bedgraph_file_transform") or die;
		while(<FH>){
     	   chomp;
      	  my @f=split "\t";
      	  $f[3] = ($f[3]-$mean_reads)/$std_reads; 
        	my $line=join "\t",@f;
        	print OUT "$line\n";
		}
		close(FH);
		close(OUT);
		if( -e $bedgraph_file_transform && $keep_bedgraph eq ""){
			unlink($bedgraph_file_minus_input);
		}	
	}
	elsif($transform == 4){
		my $bw=$bedgraph_file_minus_input;
		$bw =~ s/\.\w+$/.meanStdTmp.bw/;
		my $cmd="$WIG2BW $bedgraph_file_minus_input $CHROMSIZE $bw\n";
		print STDERR "Get adjacent 1M region's mean and standard deviation value ... \n";
		print STDERR $cmd;
		system($cmd)==0 or die;
		$bedgraph_file_transform=~s/\.\w+$/.local1M_Z_score.bedgraph/;
		open(FH,"<$bedgraph_file_minus_input") or die;
		open(OUT,">$bedgraph_file_transform") or die;
		while(<FH>){
     	   chomp;
      	  my @f=split "\t";
      	  my ($mean_reads, $std_reads)=&get_mean_std_from_bedgraph($bw, $f[0], int(($f[1]+$f[2])/2));
      	  #print STDERR "Mean: $mean_reads,  STD: $std_reads \n";
      	  $f[3] = ($f[3]-$mean_reads)/$std_reads; 
        	my $line=join "\t",@f;
        	print OUT "$line\n";
		}
		close(FH);
		close(OUT);
		unlink($bw);
		if( -e $bedgraph_file_transform && $keep_bedgraph eq ""){
			#unlink($bedgraph_file_minus_input);
		}	
	}
	elsif($transform == 3){
		#die "not support such a transformation method yet!\n";
	}
	else{
		die "not support such a transformation method yet!\n";
	}
		


##bedGraph to bigwig
	my $bigwig_file=$bedgraph_file_transform;
	if($no_bedgraph2wig eq ""){
		$bigwig_file=~s/\.\w+$/.bw/;
		my $cmd="$WIG2BW $bedgraph_file_transform $CHROMSIZE $bigwig_file\n";
		print STDERR $cmd;

		if( system($cmd)==0 && -e $bigwig_file){
			#unlink($bedgraph_file_transform);
		}
	}

}	
}
elsif($normalize == 2){
##bedGraph to bigwig	
	foreach my $bedgraph_file_norm(@bedgraph_file_norms){
		##minus RPKM value in input.bedgraph file
	my $bedgraph_file_minus_input=$bedgraph_file_norm;
	if($input_bedgraph ne ""){
		print STDERR "minus RPKM in input: $input_bedgraph ...\n";
		$bedgraph_file_minus_input=~s/\.\w+$/.minusInput.bedgraph/;
		my $cmd = "$BEDMINUS --overlap $bedgraph_file_minus_input $bedgraph_file_norm $input_bedgraph \n";
		print STDERR $cmd;
		if(system($cmd)==0 && -e $bedgraph_file_minus_input){
			#unlink($bedgraph_file_norm);
		}
	}


##transform reads
	my $bedgraph_file_transform=$bedgraph_file_minus_input;
	if($transform == 1){
		print STDERR "transform RPKM to Z-score ...\n";	
		my $std_reads = `awk '{delta = \$4 - avg; avg += delta / NR; mean2 += delta * (\$4 - avg); } END { print sqrt(mean2 / NR); }' $bedgraph_file_minus_input`;
		my $mean_reads = `awk '{sum+=\$4} END {print sum/NR}' $bedgraph_file_minus_input`;
		chomp($mean_reads);
		chomp($std_reads);
		print STDERR "RPKM mean is: $mean_reads\nRPKM standard variation is: $std_reads \n ... ...\n";	
		$bedgraph_file_transform=~s/\.\w+$/.Z_score.bedgraph/;
		open(FH,"<$bedgraph_file_minus_input") or die;
		open(OUT,">$bedgraph_file_transform") or die;
		while(<FH>){
     	   chomp;
      	  my @f=split "\t";
      	  $f[3] = ($f[3]-$mean_reads)/$std_reads; 
        	my $line=join "\t",@f;
        	print OUT "$line\n";
		}
		close(FH);
		close(OUT);
		if( -e $bedgraph_file_transform && $keep_bedgraph eq ""){
			#unlink($bedgraph_file_minus_input);
		}	
	}
	elsif($transform == 4){
		my $bw=$bedgraph_file_minus_input;
		$bw =~ s/\.\w+$/.meanStdTmp.bw/;
		my $cmd="$WIG2BW $bedgraph_file_minus_input $CHROMSIZE $bw\n";
		print STDERR "Get adjacent 1M region's mean and standard deviation value ... \n";
		print STDERR $cmd;
		system($cmd)==0 or die;
		$bedgraph_file_transform=~s/\.\w+$/.local1M_Z_score.bedgraph/;
		open(FH,"<$bedgraph_file_minus_input") or die;
		open(OUT,">$bedgraph_file_transform") or die;
		while(<FH>){
     	   chomp;
      	  my @f=split "\t";
      	  my ($mean_reads, $std_reads)=&get_mean_std_from_bedgraph($bw, $f[0], int(($f[1]+$f[2])/2));
      	  #print STDERR "Mean: $mean_reads,  STD: $std_reads \n";
      	  $f[3] = ($f[3]-$mean_reads)/$std_reads; 
        	my $line=join "\t",@f;
        	print OUT "$line\n";
		}
		close(FH);
		close(OUT);
		unlink($bw);
		if( -e $bedgraph_file_transform && $keep_bedgraph eq ""){
			#unlink($bedgraph_file_minus_input);
		}	
	}
	elsif($transform == 3){
		#die "not support such a transformation method yet!\n";
	}
	else{
		die "not support such a transformation method yet!\n";
	}
		my $bigwig_file=$bedgraph_file_transform;
		if($no_bedgraph2wig eq ""){
			$bigwig_file=~s/\.\w+$/.bw/;
			my $cmd="$WIG2BW $bedgraph_file_transform $CHROMSIZE $bigwig_file\n";
			print STDERR $cmd;

			if( system($cmd)==0 && -e $bigwig_file){
				#unlink($bedgraph_file_transform);
			}
		}
	}
	
}
elsif($normalize == 5){

	my $cmd = "$BEDQUANTNORM ";
	foreach my $bedgraph_file_binned(@bedgraph_file_binneds){
		$cmd .= "$bedgraph_file_binned ";
	}
	if($useLimma ne ""){
		$cmd .= "--useR ";
	}
	$cmd .= "\n";
	print STDERR "$cmd";
	system($cmd)==0 or die;
	
##bedGraph to bigwig
	foreach my $bedgraph_file_norm(@bedgraph_file_norms){
		$bedgraph_file_norm=~s/\.(\w+)$/.quantNorm.$1/;
		my $bigwig_file=$bedgraph_file_norm;
		if($no_bedgraph2wig eq ""){
			$bigwig_file=~s/\.\w+$/.bw/;
			my $cmd="$WIG2BW $bedgraph_file_norm $CHROMSIZE $bigwig_file\n";
			print STDERR $cmd;

			if( system($cmd)==0 && -e $bigwig_file){
				#unlink($bedgraph_file_transform);
			}
		}
	}
	
}
elsif($normalize ==3 ){
	##bedGraph to bigwig
	$wig_output=~s/\.(\w+)$/.sort.$1/;
	my $bigwig_file=$wig_output;
	if($no_bedgraph2wig eq ""){
		$bigwig_file=~s/\.\w+$/.bw/;
		my $cmd="$WIG2BW $wig_output $CHROMSIZE $bigwig_file\n";
		print STDERR $cmd;

		if( system($cmd)==0 && -e $bigwig_file){
			#unlink($bedgraph_file_transform);
		}
	}
}
else{
	die "not support such a normalization method yet!\n";
}

#get mean and std from the adjacent 1M region. if no reads there, it will extend 1M more until it get some value there.. 
sub get_mean_std_from_bedgraph{
	my $bw=shift @_;
	my $chr=shift @_;
	my $center=shift @_;
	my $start=$center;
	my $end=$center;
	my $mean="";
	my $std="";
	while($mean eq "" || $std eq "" || $std eq "n/a"|| $std eq "0"){
		$start=$start-500000;
		$start=$start >= 0 ? $start : 0;
		$end=$end+500000;
		$end=$end<($range{$chr}-1) ? $end : ($range{$chr}-1);
		$mean =`$BIGWIGSUM $bw $chr $start $end 1 `;
		chomp($mean);
		$std =`$BIGWIGSUM $bw $chr $start $end 1 -type=std`;
		chomp($std);
	}
	return ($mean, $std);

}



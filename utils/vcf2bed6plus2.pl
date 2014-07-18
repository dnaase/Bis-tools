## transfer TCGA VCF 4.1 to file in bed 6+2 format(TCGA level III data): chr + start + end + . + methylation_value*10 + strand + methylation_value(%) + numCT_reads
## numC+numT==0 sites are not included.Only homozygous CG if pattern specified as CG.  This combined +/- strand together, all of CpG positions are using '+' strand CpG location
## author: Yaping Liu  lyping1986@gmail.com

use strict;
use Getopt::Long;

##copy from Bio::Tools::IUPAC, to avoid the condition that no BioPerl in cluster machine..
my %IUB = ('A' => 'A',
	     'C' => 'C',
	     'G' => 'G',
	     'T' => 'T',
	     'U' => 'U',
	     'M' => 'A C M',
	     'R' => 'A G R',
	     'W' => 'A T W',
	     'S' => 'C G S',
	     'Y' => 'C T Y',
	     'K' => 'G T K',
	     'V' => 'A C G M R S V',
	     'H' => 'A C T H M W Y',
	     'D' => 'A G T D K R W',
	     'B' => 'C G T B K S Y',
	     'X' => 'G A T C B D H K M N R S V W X Y',
	     'N' => 'G A T C B D H K M N R S V W X Y',
);


my $NO_REF_FLAG="no_ref";
my $SNP="snp";
my $LOW_QUAL="low_qual";

my $qual=20;
my $min_ct=1;
my $max_cov=250;
my $sb=-0.02;
my $sample_order=1;
my $only_good_call="";

sub usage { ### not use these filters yet...
	print "USAGE: perl vcf2bed6plus2.pl [options] input_file_name [CG]\n";
	print "  [Options]:\n\n";
	print "  --only_good_call : only output confident call of cytosine pattern. (default: not enabled)\n\n";
	print "  --qual INT : quality score for confidently called cytosine pattern (default: 20)\n\n";
	print "  --minCT INT : minimum number of CT reads, otherwise, methy column will be 'NA' (default: 1)\n\n";
	print "  --maxCov INT : maximum coverage allowed(default: 250)\n\n";
	print "  --sample_order INT : sample order choose to output as bed file when there are multiple samples in the same vcf file. 1 means the first sample (default: 1)\n\n";
	exit(1);
}

GetOptions(
	"only_good_call" => \$only_good_call,
	"qual=i" => \$qual,
	"maxCov=i" => \$max_cov,
	"minCT=i" => \$min_ct,
	"sample_order=i" => \$sample_order,
);
my $input_file_name = $ARGV[0];
my $type = $ARGV[1];
#my $sample_name = $ARGV[2];

my $use_age = "USAGE: perl vcf2bed6plus2.pl input_file_name [CG]";
if($ARGV[0] eq ""){
	&use_age();

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


$cpg_name_output = $cpg_name_output.".$type.6plus2.bed";
open(OUT,">$cpg_name_output") or die;

if($type eq "CG"){
	$type = "CG"; 
}
else{
	if($type eq ""){
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
	
}

my $head_line = "track name=$cpg_name_output type=bedDetail description=\"$type methylation level-$bissnp_version\" visibility=3";
#variableStep chrom=chr19 span=150"
print OUT "$head_line\n";

#only previous + strand CpG is recorded here
my $num_c_pre=-1;
my $num_t_pre=-1;
my $methy_pre="NA";
my $chr_pre="";
my $start_pre="";
my $end_pre="";
my $filter=".";

open(FH,"<$input_file_name") or die;
while(<FH>){
	my $line=$_;
	chomp($line);
	next if $line =~ /^#/;
	my @splitin = split "\t", $line;
	
	
	
	my $sample_id = 8 + $sample_order;
	if($splitin[$sample_id] =~ /:(\d+):($type):(\d+):/){
		my $num_c = $1;
		my $context = $2;
		my $num_t = $3;
		my $strand=".";
		if($splitin[7] =~ /CS=([+|-])/){
			$strand = $1;
		}
		if($type eq "\\w+"){ ##if mode is to output all of cytosine in VCF, then strand combination is not applied
			my $methy = 'NA';
		 	my $score = sprintf("%d",0);
			if($num_c + $num_t >= $min_ct){
				$methy = $num_c/($num_c + $num_t);
				$score = sprintf("%d",1000*$methy);
				$methy = sprintf("%.2f",100*$methy);
			}
				$filter=".";
				my $chr = $splitin[0];
				my $start = $splitin[1]-1;
				my $end = $splitin[1];
		 		my $reads_num = $num_c + $num_t;
				my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
				print OUT "$out_line\n";
				
			#}
			next;
		}
		
		if($strand eq '+'){
			
			if($num_c_pre == -1){ 
				if($chr_pre eq ""){## no previous record, just remember this CpG location
					
				}
				else{## some previous bad called record, just output and then remember this CpG location
						my $chr = $chr_pre;
						my $start = $start_pre;
						my $end = $end_pre;
		 				my $reads_num = 'NA';
		 				my $methy = 'NA';
		 				my $score = sprintf("%d",0);
						my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
						$filter=".";
				}
			}
			else{		## there are previous record, output previous record and remember this CpG location
				my $methy = 'NA';
		 		my $score = sprintf("%d",0);
				if($num_c_pre + $num_t_pre  >= $min_ct and $methy_pre ne 'NA'){
					$methy = $num_c_pre/($num_c_pre + $num_t_pre);
					$score = sprintf("%d",1000*$methy);
					$methy = sprintf("%.2f",100*$methy);
				}
					my $chr = $chr_pre;
					my $start = $start_pre;
					my $end = $end_pre;
		 			my $reads_num = $num_c_pre + $num_t_pre;
					my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
					print OUT "$out_line\n";
					$filter=".";
				#}
			}
			$num_c_pre = $num_c;
			$num_t_pre = $num_t;
			$methy_pre = 1;
			$chr_pre=$splitin[0];
			$start_pre=$splitin[1]-1;
			$end_pre=$splitin[1];
			&pass_filter(@splitin);
		}
		else{ ## - strand
			if($num_c_pre == -1){  ## no previous record, just output this CpG location
				my $methy = 'NA';
		 		my $score = sprintf("%d",0);
				if($num_c + $num_t  >= $min_ct){
					$methy = $num_c/($num_c + $num_t);
					$score = sprintf("%d",1000*$methy);
					$methy = sprintf("%.2f",100*$methy);
				}
				$filter=~s/$LOW_QUAL/./;
					my $chr = $splitin[0];
					my $start = $splitin[1]-2;
					my $end = $splitin[1]-1;
		 			my $reads_num = $num_c + $num_t;
					my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
					print OUT "$out_line\n";
					$filter=".";
				#}
			}
			else{ ## there are previous record, then check continuoucy with current location
				if(($splitin[0] eq $chr_pre) and ($splitin[1] == ($end_pre+1) )){ #same chromosome, and previous record is one base before
					my $methy = 'NA';
		 			my $score = sprintf("%d",0);
		 			my $reads_num = $num_c + $num_t;
		 			if($methy_pre ne 'NA'){
		 				if($num_c_pre + $num_t_pre + $num_c + $num_t  >= $min_ct){
							$methy = ($num_c + $num_c_pre)/($num_c_pre + $num_t_pre + $num_c + $num_t);
							$reads_num = $reads_num + $num_c_pre + $num_t_pre;
							$score = sprintf("%d",1000*$methy);
							$methy = sprintf("%.2f",100*$methy);
						}
		 			}
		 			else{
		 				if($num_c + $num_t  >= $min_ct){
							$methy = ($num_c )/( $num_c + $num_t);
							$score = sprintf("%d",1000*$methy);
							$methy = sprintf("%.2f",100*$methy);
						}
		 			}
					
						my $chr = $chr_pre;
						my $start = $start_pre;
						my $end = $end_pre;
		 				
		 				$filter=~s/$LOW_QUAL/./;
						my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
						$filter=".";
					#}
					
				}
				else{ ##not continuous, then output previous and current location, reset $num_c_pre.
					my $methy = 'NA';
		 			my $score = sprintf("%d",0);
					if($num_c_pre + $num_t_pre  >= $min_ct and $methy_pre ne 'NA'){
						$methy = $num_c_pre/($num_c_pre + $num_t_pre);
						$score = sprintf("%d",1000*$methy);
						$methy = sprintf("%.2f",100*$methy);
					}
						my $chr = $chr_pre;
						my $start = $start_pre;
						my $end = $end_pre;
		 				my $reads_num = $num_c_pre + $num_t_pre;
		 				$filter=~s/$LOW_QUAL/./;
						my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
						$filter=".";
					#}
					$methy = 'NA';
		 			$score = sprintf("%d",0);
					if($num_c + $num_t  >= $min_ct){
						$methy = $num_c/($num_c + $num_t);
						my $score = sprintf("%d",1000*$methy);
						$methy = sprintf("%.2f",100*$methy);
					}
						$chr = $splitin[0];
						$start = $splitin[1]-2;
						$end = $splitin[1]-1;
		 				$reads_num = $num_c + $num_t;
						$out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
						$filter=".";
					#}
				}
			}
			
			$num_c_pre=-1;
			$num_t_pre=-1;
			$methy_pre = 'NA';
			$chr_pre="";
			$start_pre="";
			$end_pre="";
			&pass_filter(@splitin);
		}
		
		
		
	}elsif($only_good_call ne ""){ ## if only output confident call of CG, then do nothing here.
		$filter=".";
	}
	else{
		my $num_c=-1;
		my $num_t=-1;
		if($splitin[$sample_id]=~ /\d\/\d:\S+:(\d+),(\d+),\d+,\d+,\d+,\d+/){
			$num_c=$1;
			$num_t=$2;
		}
		if($splitin[7] =~ /REF=.*${type}.*/ or ($splitin[7] =~ /REF=0/ and $splitin[7] =~ /Context=.*${type}.*/)){	
			
			#if(($splitin[3] eq 'C' and $splitin[$sample_id] =~ /0\//) or  ($splitin[4] =~ /C/ and ($splitin[$sample_id] =~ /1\// or $splitin[$sample_id] =~ /\d\/1/)) or ($splitin[4] =~ /,C/ and ($splitin[$sample_id] =~ /2\// or $splitin[$sample_id] =~ /\d\/2/))){ ##if not confidently called as CG, and it is in forward strand (C strand)
			if(($splitin[3] eq 'C')){	
				if($methy_pre eq 'NA'){ 
					if($chr_pre eq ""){## no previous record, just remember this CpG location
					
					}
					else{## some previous bad called record, just output and then remember this CpG location
						my $chr = $chr_pre;
						my $start = $start_pre;
						my $end = $end_pre;
		 				my $reads_num = 'NA';
		 				if($num_c_pre != -1){
		 					$reads_num = $num_c + $num_t;
		 				}
		 				my $methy = 'NA';
		 				my $score = sprintf("%d",0);
						my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
						$filter=".";
					}
				}
				else{		## there are previous record, output previous record and remember this CpG location
					my $methy = 'NA';
		 			my $score = sprintf("%d",0);
					if($num_c_pre + $num_t_pre  >= $min_ct){
						$methy = $num_c_pre/($num_c_pre + $num_t_pre);
						$score = sprintf("%d",1000*$methy);
						$methy = sprintf("%.2f",100*$methy);
					}
						my $chr = $chr_pre;
						my $start = $start_pre;
						my $end = $end_pre;
		 				my $reads_num = $num_c_pre + $num_t_pre;
						my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
						$filter=".";
				#}
				}
				
				$num_c_pre=$num_c;
				$num_t_pre=$num_t;
				$methy_pre = 'NA';
				$chr_pre=$splitin[0];
				$start_pre=$splitin[1]-1;
				$end_pre=$splitin[1];
				&pass_filter(@splitin);
			}
			#elsif(($splitin[3] eq 'G' and $splitin[$sample_id] =~ /0\//) or ($splitin[4] =~ /G/ and ($splitin[$sample_id] =~ /1\// or $splitin[$sample_id] =~ /\d\/1/)) or ($splitin[4] =~ /,G/ and ($splitin[$sample_id] =~ /2\// or $splitin[$sample_id] =~ /\d\/2/))){##if not confidently called as CG, and it is in reverse strand (G strand), juts output 'NA' in methy column
			elsif(($splitin[3] eq 'G')){	
					
				if(($splitin[0] eq $chr_pre) and ($splitin[1] == ($end_pre+1) )){ #same chromosome, and previous record is one base before
							my $chr = $chr_pre;
							my $start = $start_pre;
							my $end = $end_pre;
		 					my $reads_num = 'NA';
		 					my $methy = 'NA';
		 					my $score = sprintf("%d",0);
						if($num_c_pre != -1){
							if($methy_pre ne 'NA'){
								$reads_num = $num_c_pre + $num_t_pre;
							}
							else{
								$reads_num = $num_c_pre + $num_t_pre +  $num_c + $num_t;
							}
							
							if($num_c_pre + $num_c_pre  >= $min_ct and $methy_pre ne 'NA'){
								
		 						$methy = $num_c_pre/($num_c_pre + $num_t_pre);
								$score = sprintf("%d",1000*$methy);
								$methy = sprintf("%.2f",100*$methy);
								$filter=".";
							}
							
							
						}
						
		 				
						my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
						print OUT "$out_line\n";
							
						$filter=".";
						
				}
				else{##not continuous, then output previous and current location, reset $num_c_pre.
						
		 					my $reads_num = 0;
		 					my $methy = 'NA';
		 					my $score = sprintf("%d",0);
		 					
		 					if($chr_pre ne ""){
		 						my $chr = $chr_pre;
								my $start = $start_pre;
								my $end = $end_pre;
								if($num_c_pre != -1){
									$reads_num = $num_c_pre + $num_t_pre;
									if($num_c_pre + $num_c_pre  >= $min_ct and $methy_pre ne 'NA'){
										
		 								$methy = $num_c_pre/($num_c_pre + $num_t_pre);
										$score = sprintf("%d",1000*$methy);
										$methy = sprintf("%.2f",100*$methy);
									}
							
							
								}
								my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
								print OUT "$out_line\n";
								$filter=".";
		 					}
		 					$reads_num = $num_c + $num_t;
		 					my $chr = $splitin[0];
							my $start = $splitin[1]-2;
							my $end = $splitin[1]-1;
							&pass_filter(@splitin);
							my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
							print OUT "$out_line\n";
		 					$filter=".";
		 					
						
				}
				
				$num_c_pre=-1;
				$num_t_pre=-1;
				$methy_pre = 'NA';
				$chr_pre="";
				$start_pre="";
				$end_pre="";
			}
		}
		elsif($splitin[7] =~ /Context=(.*?);/){
			#if($splitin[1] == 7004290){
			#	print "$methy_pre\t$num_c\t$num_t\t$num_c_pre\t$num_t_pre\n";
			#}
			#print "$splitin[1]\n";
			my @contexts=split ",",$1;
			foreach my $context(@contexts){
				#my $flag=&iupac_match($context,$type);
				
				#print "$splitin[1]\t$flag\t$context\t$type\n";
				if(&iupac_match($context,$type) == 1){
						my $strand=".";
						if($splitin[7] =~ /CS=([+|-])/){
								$strand = $1;
						}
						
						  if($strand eq '+'){
						  	
								if($methy_pre eq 'NA'){ 
									if($chr_pre eq ""){## no previous record, just remember this CpG location
					
									}
									else{## some previous bad called record, just output and then remember this CpG location
										my $chr = $chr_pre;
										my $start = $start_pre;
										my $end = $end_pre;
		 								my $reads_num = 'NA';
		 								if($num_c_pre != -1){
		 										$reads_num = $num_c_pre + $num_t_pre;
		 									
		 								}
		 								my $methy = 'NA';
		 								my $score = sprintf("%d",0);
										my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
										print OUT "$out_line\n";
										$filter=".";
									}
								}
								else{		## there are previous record, output previous record and remember this CpG location
									my $methy = 'NA';
		 							my $score = sprintf("%d",0);
									if($num_c_pre + $num_t_pre  >= $min_ct){
										$methy = $num_c_pre/($num_c_pre + $num_t_pre);
										$score = sprintf("%d",1000*$methy);
										$methy = sprintf("%.2f",100*$methy);
									}
									my $chr = $chr_pre;
									my $start = $start_pre;
									my $end = $end_pre;
		 							my $reads_num = $num_c_pre + $num_t_pre;
									my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
									print OUT "$out_line\n";
									$filter=".";
				#}
								}
								$num_c_pre = $num_c;
								$num_t_pre = $num_t;
								$methy_pre = 1;
								$chr_pre=$splitin[0];
								$start_pre=$splitin[1]-1;
								$end_pre=$splitin[1];
								&pass_filter(@splitin);
						}
						else{
							if($methy_pre eq 'NA'){  ## no previous record, just output this CpG location
									if($chr_pre eq ""){## no previous record, just output this CpG location
					
									}
									else{## some previous bad called record, just output previous and this CpG location
										if(($splitin[0] eq $chr_pre) and ($splitin[1] == ($end_pre+1) )){
											
										}
										else{
											my $chr = $chr_pre;
											my $start = $start_pre;
											my $end = $end_pre;
		 									my $reads_num = 'NA';
		 									if($num_c_pre != -1){
		 										$reads_num = $num_c_pre + $num_t_pre;
		 									}
		 									my $methy = 'NA';
		 									my $score = sprintf("%d",0);
											my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
											print OUT "$out_line\n";
											$filter=".";
											$num_c_pre=-1;
											$num_t_pre=-1;
										}
										
									}
									
								my $methy = 'NA';
		 						my $score = sprintf("%d",0);
								
								my $chr = $splitin[0];
								my $start = $splitin[1]-2;
								my $end = $splitin[1]-1;
		 						my $reads_num = 'NA';
		 						if($num_c != -1){
		 							$reads_num = $num_c + $num_t;	
		 						}
		 						if($num_c_pre != -1){
		 							$reads_num = $num_c_pre + $num_t_pre + $reads_num;

		 						}
		 						$filter=".";
		 						&pass_filter(@splitin);
								my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
								print OUT "$out_line\n";
								$filter=".";
				#}
							}
							else{ ## there are previous record, then check continuoucy with current location
								if(($splitin[0] eq $chr_pre) and ($splitin[1] == ($end_pre+1) )){ #same chromosome, and previous record is one base before
									my $methy = 'NA';
		 							my $score = sprintf("%d",0);
									if($num_c_pre + $num_t_pre >= $min_ct){
										$methy = $num_c_pre/($num_c_pre + $num_t_pre);
										$score = sprintf("%d",1000*$methy);
										$methy = sprintf("%.2f",100*$methy);
									}
									my $chr = $chr_pre;
									my $start = $start_pre;
									my $end = $end_pre;
		 							my $reads_num = $num_c_pre + $num_t_pre;
		 							&pass_filter(@splitin);
									my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
									print OUT "$out_line\n";
									$filter=".";
					#}
					
								}
								else{ ##not continuous, then output previous and current location, reset $num_c_pre.
									my $methy = 'NA';
		 							my $score = sprintf("%d",0);
									if($num_c_pre + $num_t_pre  >= $min_ct and $methy_pre ne 'NA'){
										$methy = $num_c_pre/($num_c_pre + $num_t_pre);
										$score = sprintf("%d",1000*$methy);
										$methy = sprintf("%.2f",100*$methy);
									}
									my $chr = $chr_pre;
									my $start = $start_pre;
									my $end = $end_pre;
		 							my $reads_num = $num_c_pre + $num_t_pre;
									my $out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
									print OUT "$out_line\n";
									$filter=".";
					#}
									$methy = 'NA';
		 							$score = sprintf("%d",0);
									
									$chr = $splitin[0];
									$start = $splitin[1]-2;
									$end = $splitin[1]-1;
		 							$reads_num = 'NA';
		 							if($num_c != -1){
		 								$reads_num = $num_c + $num_t;	
		 							}
		 							&pass_filter(@splitin);
									$out_line = "$chr\t$start\t$end\t$filter\t$score\t\.\t$methy\t$reads_num";
									print OUT "$out_line\n";
									$filter=".";
					#}
								}
							}
			
							$num_c_pre=-1;
							$num_t_pre=-1;
							$methy_pre='NA';
							$chr_pre="";
							$start_pre="";
							$end_pre="";
							$filter=".";
						}
				}
			}
		}
			
	}
	
}

close(FH);
close(OUT);
print "finished!\n";


sub pass_filter{
	my @ins=@_;
	my $flag1 = &low_qual(@ins);
	#my $flag2 = &snps($flag1, @ins);	
	if($flag1==0){
		&not_in_ref($flag1, @ins);
		return 1;
	}
	return 0;
}

## 1 means low qual, 0 means high qual
sub low_qual{
	my @ins=@_;
	my $sample_id = 8 + $sample_order;
	
	if($ins[$sample_id] =~ /\d\/\d:\S+:\S+:\S+:(\S+):\S+:(\d+):\S+:\S+:(\S+):\d+/){
		my $word = $1;
		my $dp=$2;
		my $gq=$3;
		
		#print "$dp\t$max_cov\t$gq\t$qual\n";
		if($dp >=$max_cov || $gq < $qual){
			if($filter eq '.'){
				$filter = $LOW_QUAL;
			}
			else{
				$filter .= ";$LOW_QUAL";
			}
			return 1;
		}
		else{
			if($type ne $word){
				if($filter eq '.'){
					$filter = $SNP;
				}
				else{
					$filter .= ";$SNP";
				}
				return 1;
			}
		}
			
	}
	
	return 0;
}

sub not_in_ref{
	my $flag=shift @_;
	my @ins=@_;
	if($flag == 0 ){
		if($ins[7] =~ /REF=0/ || $ins[7] =~ /REF=.*${type}.*/){
			
		}
		else{
			if($filter eq '.'){
				$filter = $NO_REF_FLAG;
			}
			else{
				$filter .= ";$NO_REF_FLAG";
			}
		}	
	}
}

sub iupac_match{
	my $query=shift @_;
	my $pattern=shift @_;
	if(length($query) != length($pattern)){
		return 0;
	}
	
	my @querys=split "",$query;
	my @patterns=split "",$pattern;
	for(my $i=0; $i<$#patterns; $i++){ 
		
		if($IUB{$querys[$i]} !~ /$patterns[$i]/){
			return 0;
		}
		
	}
	return 1;
	
}






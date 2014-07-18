#!/usr/bin/perl -w

#
# invert_dups_check.pl: This script are used to identify inverted duplication reads for paired end reads 
# on Illumina platform Copyright (C) <2013>
# <Yaping Liu: lyping1986@gmail.com>
# 
# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program. If
# not, see <http://www.gnu.org/licenses/>.
#



## take care, 1st end, T should be also considered as C, while in 2nd end, A should be also considered as G

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-1-8

#Usage:  perl invert_dups_check.pl [option] output_log.txt sample.1st_end.fastq sample.2nd_end.fastq 

##example1: trim_mode_1: ls *_F.*.fastq | perl -ne 'chomp;$r2=$_;$r2=~s/_F\./_R./;$log=$_;$log=~s/\.fastq$/.log.txt/;$start=$_;$start=~s/\.fastq$/.trim_start_pos.txt/;`perl ~/qsubexc4g.pl "perl /home/rcf-40/yapingli/storage/code/mytools/perl/invert_dups_check.pl $log $_ $r2 --trim_invert_dups --output_invert_dups --output_clean --output_dirty --output_trim_start $start" inv_dups_hunter_$_`'

##example2: trim_mode_2: ls *_F.*.fastq | perl -ne 'chomp;$r2=$_;$r2=~s/_F\./_R./;$log=$_;$log=~s/\.fastq$/.log.txt/;$start=$_;$start=~s/\.fastq$/.trim_start_pos.txt/;`perl ~/qsubexc4g.pl "perl /home/rcf-40/yapingli/storage/code/mytools/perl/invert_dups_check.pl $log $_ $r2 --trim_mode_2 --trim_invert_dups --output_clean --output_trim_start $start" inv_dups_hunter_$_`'


use strict;
use Getopt::Long;



my $mismatch=3;
my $first_N=20;
my $non_bisulfite="";
my $output_invert_dups="";
my $output_clean="";
my $output_dirty="";
my $trim_invert_dups="";
my $output_trim_start = "";
my $trim_mode_2 = "";
my $min_len=30;


sub usage {
	
    print STDERR "\nUsage:\n";
    print STDERR "perl invert_dups_check.pl [Options]  output_log.txt sample.1st_end.fastq sample.2nd_end.fastq\n\n";

    print STDERR "This script are used to identify inverted duplication reads for paired end reads. reads in the two fastq file should be in the same order as each other.\n";
    print STDERR "If the first 20bp of the reads have 0 mismatches with their mates, they are considered to be contaminated by inverted dups.If you specify --trim_invert_dups, \n";
    print STDERR "then reads will be trimmed after the mismatches are more than 3. It provides a way before reads mapping to identify inverted dups rather than counting the same starting coordinate for mapped reads.\n\n";
	print STDERR "[Options]:\n\n";
	
	print STDERR "  --mismatch INT : maximum mismathces allowed for trim the inverted dups. (Default: 3)\n\n";
	print STDERR "  --first_N INT : First N bases to count for the matches. (Default: 20)\n\n";
	print STDERR "  --trim_invert_dups : Output reads that has trimmed inverted dups already, seperately on two different ends. (Default: disabled)\n\n";
	print STDERR "  --output_invert_dups : Output the whole reads that contained inverted dups, seperately on two different ends. (Default: disabled)\n\n";
	print STDERR "  --output_clean : Output the whole clean reads not contain any inverted dups. (Default: disabled)\n\n";
	print STDERR "  --output_dirty : Output the part of reads that has been trimmed off. (Default: disabled)\n\n";
	print STDERR "  --output_trim_start FILE : Output the trimming starting position. (Default: disabled)\n\n";
	print STDERR "  --non_bisulfite : Non Bisulfite mode. (Default: disabled)\n\n";
	print STDERR "  --trim_mode_2 : Seperate 1st end into two ends, the beginning of the 1st end as 1st end, while the last part of 1st end reads as 2nd end (Default: disabled)\n\n";
	print STDERR "  --min_len INT: minimum read1+read2 length after trimming. (Default: 30)\n\n";
	
	print STDERR "For any question, please contact us by email: lyping1986\@gmail.com\n\n";
	
	exit(1);
}

GetOptions(
	"mismatch=i" => \$mismatch,
	"first_N=i" => \$first_N,
	"non_bisulfite" => \$non_bisulfite,
	"output_invert_dups" => \$output_invert_dups,
	"output_clean" => \$output_clean,
	"output_dirty" => \$output_dirty,
	"trim_invert_dups" => \$trim_invert_dups,
	"output_trim_start=s" => \$output_trim_start,
	"trim_mode_2" => \$trim_mode_2,
	"min_len=i" => \$min_len,
);

my $output_log = $ARGV[0];
my $first_end = $ARGV[1];
my $second_end = $ARGV[2];

usage() if ( scalar(@ARGV) < 3 );


my $linecount_within = 1;
my $linecount_global = 1;
my $enough_long_pairs=0;
my $num_inv_dups = 0;
my $seq_so_far_1st = "";
my $seq_so_far_2nd = "";
my $inv_dups_flag = 0;
my $trim_start_pos=-1;
my $trim_start_pos_2=-1; #to record start position at the other end, in mode 2

my $seq_within_1st_line1 = "";
my $seq_within_1st_line3 = "";

my $seq_within_2nd_line1 = "";
my $seq_within_2nd_line3 = "";



open(FH1,"<$first_end") or die "can't open read file $first_end: $!\n";
open(FH2,"<$second_end") or die "can't open read file $second_end: $!\n";

if($output_trim_start ne ""){
	open(TRIMSTART,">$output_trim_start") or die "can't open trimming start position file $output_trim_start: $!\n";
	
}

if($output_invert_dups ne ""){
	my $invert_dups_1st = $first_end;
	$invert_dups_1st =~ s/\.\w+$/.invert_dups.1st.fastq/;
	open(INVDUPS1ST, ">$invert_dups_1st") or die "can't open invert dups file $invert_dups_1st: $!\n";
	my $invert_dups_2nd = $second_end;
	$invert_dups_2nd =~ s/\.\w+$/.invert_dups.2nd.fastq/;
	open(INVDUPS2ND, ">$invert_dups_2nd") or die "can't open invert dups file $invert_dups_2nd: $!\n";
}

if($output_clean ne ""){
	my $invert_dups_1st = $first_end;
	$invert_dups_1st =~ s/\.\w+$/.clean.1st.fastq/;
	open(NONCONTAM1ST, ">$invert_dups_1st") or die "can't open invert dups file $invert_dups_1st: $!\n";
	my $invert_dups_2nd = $second_end;
	$invert_dups_2nd =~ s/\.\w+$/.clean.2nd.fastq/;
	open(NONCONTAM2ND, ">$invert_dups_2nd") or die "can't open invert dups file $invert_dups_2nd: $!\n";
}

if($output_dirty ne ""){
	my $invert_dups_1st = $first_end;
	$invert_dups_1st =~ s/\.\w+$/.dirty.1st.fastq/;
	open(CONTAM1ST, ">$invert_dups_1st") or die "can't open invert dups file $invert_dups_1st: $!\n";
	my $invert_dups_2nd = $second_end;
	$invert_dups_2nd =~ s/\.\w+$/.dirty.2nd.fastq/;
	open(CONTAM2ND, ">$invert_dups_2nd") or die "can't open invert dups file $invert_dups_2nd: $!\n";
}


if($trim_invert_dups ne ""){
	my $clean_1st = $first_end;
	$clean_1st =~ s/\.\w+$/.invert_dups_trimmed.fastq/;
	open(AFTERTRIM1ST, ">$clean_1st") or die "can't open invert dups file $clean_1st: $!\n";
	#my $clean_2nd = $second_end;
	#$clean_2nd =~ s/\.\w+$/.invert_dups_trimmed.2nd.fastq/;
	#open(AFTERTRIM2ND, ">$clean_2nd") or die "can't open invert dups file $clean_2nd: $!\n";
	
}

while(my $seq1 = <FH1>, my $seq2 = <FH2>){
	$seq_so_far_1st .= $seq1;
    $seq_so_far_2nd .= $seq2;
    
	chomp($seq1);
	chomp($seq2);
		if ($linecount_within == 1)
        {
            $seq_within_1st_line1 = $seq1;
        	$seq_within_2nd_line1 = $seq2;
            # Double check the format
            if($seq1 =~ /^\@/ and $seq2 =~ /^\@/){
            	my $tmp1 = $seq1;
            	$tmp1=~s/^(\S+)\s+.*$/$1/;
            	my $tmp2 = $seq2;
            	$tmp2=~s/^(\S+)\s+.*$/$1/;
            	$tmp1 =~ s/\/1$/\/2/;

            	if($tmp1 ne $tmp2){
            		print STDERR "Incorrect FASTQ file \nLine ${linecount_global}: reads $seq1 are not paired and the same order in both of ends fastq files : $tmp1 and $tmp2\n";
            	}
            	else{
            		
            	}
            }
            else{
            	print STDERR "Incorrect FASTQ file \nLine ${linecount_global}: $seq1\nMod4 lines should start with \@\n";

            }
           
  
        }
        elsif ($linecount_within == 2)
        {
			if(length($seq1) >= $first_N and length($seq2) >= $first_N){
					if(&match($seq1, $seq2)){
						$num_inv_dups++;
						$inv_dups_flag = 1;
						my $original_seq1 = $seq1;
						my $original_seq2 = $seq2;
						$seq1 = &trim($seq1, $seq2, 1); ##return 1st sequence' content
						$trim_start_pos = length($seq1);
						if($trim_mode_2 ne ""){
							my $tmp1=reverse($original_seq1);
							my $tmp2=reverse($original_seq2);
							$seq2 = &trim($tmp1, $tmp2, 2);##return 2nd sequence' content
							#print "start:".$original_seq1."\n".reverse($original_seq1)."\n".$original_seq2."\n". reverse($original_seq2)."\n$tmp2\n";
							$trim_start_pos_2 = length($seq2);
							$seq2 = reverse($seq2);
							
						}else{
							$seq2 = &trim_reads($trim_start_pos, $seq2);
							##test:
							#$seq2 = &reverse_complement($seq2);
						}
						
						if($output_trim_start ne ""){
							if($trim_mode_2 ne ""){
								print TRIMSTART "$trim_start_pos\t$trim_start_pos_2\n";
							}else{
								print TRIMSTART "$trim_start_pos\n";
							}
							
						}
						
						 if($trim_invert_dups ne "" and (length($seq1)+length($seq2)) >= $min_len and length($seq2)>0){
        					print AFTERTRIM1ST "$seq_within_1st_line1\n$seq1\n";
 		        			#print AFTERTRIM2ND "$seq_within_2nd_line1\n$seq2\n";      		
    					 }
						
						if($output_dirty ne ""){
								my $seq_contam_1 = &reads_trimmed_out($trim_start_pos,$original_seq1);
								my $seq_contam_2;
								if($trim_mode_2 ne ""){
									my $tmp=reverse($original_seq2);
									$seq_contam_2= &reads_trimmed_out($trim_start_pos_2,$tmp);
									$seq_contam_2 = reverse($seq_contam_2);
								}else{
									$seq_contam_2= &reads_trimmed_out($trim_start_pos,$original_seq2);
								}
								 
								if($seq_contam_1 ne '0'){
									print CONTAM1ST "$seq_within_1st_line1\n$seq_contam_1\n";
								}
								if($seq_contam_2 ne '0'){
									print CONTAM2ND "$seq_within_2nd_line1\n$seq_contam_2\n";
								}
						}
						
					}
					$enough_long_pairs++;
			}

        }
        elsif ($linecount_within == 3)
        {
        	$seq_within_1st_line3 = $seq1;
        	$seq_within_2nd_line3 = $seq2;
        }
        elsif ($linecount_within == 4)
        {
        	if($inv_dups_flag != 0 ){
        		my $original_seq1 = $seq1;
				my $original_seq2 = $seq2;
        		$seq1 = &trim_reads($trim_start_pos, $seq1); ##go back to see if we need the 2nd end reads also trimmed here?
        		if($trim_mode_2 ne ""){
					my $tmp2=reverse($seq2);
        			$seq2 = &trim_reads($trim_start_pos_2, $tmp2);
        			$seq2 = reverse($seq2);
        		}else{
        			$seq2 = &trim_reads($trim_start_pos, $seq2);
        			##test:
        			#$seq2 = reverse($seq2)
        		}
				
				
				if($trim_invert_dups ne "" and (length($seq1)+length($seq2)) >= $min_len and length($seq2)>0){
        					print AFTERTRIM1ST "$seq_within_1st_line3\n$seq1\n";
 		        			#print AFTERTRIM2ND "$seq_within_2nd_line3\n$seq2\n";      		
    			}
    					 
        		if($output_invert_dups ne ""){
        			print INVDUPS1ST $seq_so_far_1st;
        			print INVDUPS2ND $seq_so_far_2nd;
        		}
        		
        		if($output_dirty ne ""){
						my $seq_contam_1 = &reads_trimmed_out($trim_start_pos,$original_seq1);
						my $seq_contam_2;
						if($trim_mode_2 ne ""){
									my $tmp=reverse($original_seq2);
									$seq_contam_2= &reads_trimmed_out($trim_start_pos_2,$tmp);
									$seq_contam_2 = reverse($seq_contam_2);
						}else{
									$seq_contam_2= &reads_trimmed_out($trim_start_pos,$original_seq2);
						}
						
						if($seq_contam_1 ne '0'){
							print CONTAM1ST "$seq_within_1st_line3\n$seq_contam_1\n";
						}
						if($seq_contam_2 ne '0'){
							print CONTAM2ND "$seq_within_2nd_line3\n$seq_contam_2\n";
						}
				}
        		
        	}
        	else{
        		if($output_clean ne ""){
        			print NONCONTAM1ST $seq_so_far_1st;
        			print NONCONTAM2ND $seq_so_far_2nd;
        		}
        		
        	}
        	
        	$linecount_within = 0;
        	$seq_so_far_1st = "";
            $seq_so_far_2nd = "";
            $inv_dups_flag = 0;
            $trim_start_pos = -1;
            $trim_start_pos_2 = -1;
            
            $seq_within_1st_line1 = "";
			$seq_within_1st_line3 = "";
			$seq_within_2nd_line1 = "";
			$seq_within_2nd_line3 = "";
        }
        
   if($linecount_global/4 % 1000000 == 0){
   	print STDERR "$linecount_global/4 reads processed ... \n\n";
   }
     
	$linecount_within++;
	$linecount_global++;
}
close(FH1);
close(FH2);
if($output_trim_start ne ""){
	close(TRIMSTART);
}

if($output_invert_dups ne ""){
	close(INVDUPS1ST);
	close(INVDUPS2ND);
}
if($trim_invert_dups ne ""){
	close(AFTERTRIM1ST);
	#close(AFTERTRIM2ND);
}
if($output_clean ne ""){
	close(NONCONTAM1ST);
	close(NONCONTAM2ND);
}
if($output_dirty ne ""){
	close(CONTAM1ST);
	close(CONTAM2ND);
}

#my $percentage = sprintf("%.2f",100*$num_inv_dups/$enough_long_pairs);
#print "There are $enough_long_pairs enough long (>=$first_N) pair of reads in the fastq file\n";
#print "There are $num_inv_dups inverted dups reads (${percentage}%) in the fastq file\n";
#$enough_long_pairs *= 2;
## take care!! Here is not the real mapped reads!! Here is the reads that in the fastq file. This is just for USCEC pipeline convenient. 
$linecount_global--;
$linecount_global /=2 ;
my $perc=100*$num_inv_dups/($linecount_global/2);
if($perc<0.000001){
	$perc=0.000001
}
my $percentage = sprintf("%.6f",$perc);
open(OUT,">$output_log") or die;
#print OUT "mapped reads=$linecount_global\n";
$linecount_global/=2;
print OUT "Total pair of reads=$linecount_global\n";
print OUT "Inverted Read Pairs=$num_inv_dups\n";
print OUT "inverted Pair Percentage=$percentage\n";
close(OUT);

sub match{
	my $seq= shift @_;
	my $query = shift @_;
	my @seqs=split "",$seq;
	my @queries = split "",$query;
	my $mismatch_count = 0;
	for(my $i=0;$i<$first_N;$i++){
		if($non_bisulfite ne ""){
			if($seqs[$i] ne $queries[$i]){

				$mismatch_count++;
			}
		}
		else{
			if(&bisulfite_match($seqs[$i],$queries[$i])){
			
			}
			else{
				$mismatch_count++;
			}
		}
		
		if($mismatch_count > 0){  ## or be 1?
			return 0;
		}
	}
	return 1;
	
}

sub bisulfite_match{
	my $first = shift @_;
	my $second = shift @_;
	if($first eq $second){
		return 1;
	}
	else{
		if(($first eq 'T' and $second eq 'C') or ($first eq 'G' and $second eq 'A')){
			return 1;
		}
	}
	return 0;
}

sub trim_reads{
	my $num= shift @_;
	my $read = shift @_;
	my $trimmed = "";
	if(length($read) < $num){
		return $read;
	}
	
	my @chars = split "",$read;
	
	for(my $i = 0; $i < $num; $i++){
		$trimmed .= $chars[$i];
	}
	
	return $trimmed;
}

sub trim{
	my $seq= shift @_;
	my $query = shift @_;
	my $mode = shift @_;
	my @seqs=split "",$seq;
	my @queries = split "",$query;
	my $mismatch_count = 0;
	my $len = length($seq) >= length($query) ? length($query) : length($seq);
	my $trimmed = "";
	my @chars;
	if($mode == 1){
		 @chars = split "",$seq;
		  ##in mode 2, return 2nd end read sequences..
	}elsif($mode == 2){
		 @chars = split "",$query;
	}else{
		die "not such a trim mode: $mode \n";
	}
	for(my $i=0;$i<$len;$i++){
		$trimmed .= $chars[$i];
		if($non_bisulfite ne ""){
			if($seqs[$i] ne $queries[$i]){

				$mismatch_count++;
			}
		}
		else{
			if(&bisulfite_match($seqs[$i],$queries[$i])){
			
			}
			else{
				$mismatch_count++;
			}
		}
		
		if($mismatch_count > $mismatch){  
			return $trimmed;
		}
	}
	return $trimmed;
	
}

sub reads_trimmed_out{
	my $num= shift @_;
	my $read = shift @_;
	my $trimmed = "";
	if(length($read) <= $num){
		return 0;
	}
	
	my @chars = split "",$read;

	for(my $i = $num-1; $i < $#chars; $i++){
		$trimmed .= $chars[$i];
	}
	
	return $trimmed;
}

sub complement {
        my $dna = shift;


	# complement the DNA sequence
        $dna =~ tr/ACGTacgt/TGCAtgca/;
        return $dna;
}

sub reverse_complement {
        my $dna = shift;
		$dna =  reverse($dna);

	# complement the DNA sequence
        $dna =~ tr/ACGTacgt/TGCAtgca/;
        return $dna;
}


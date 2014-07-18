#!/usr/bin/perl -w
## value in one begraph minus the value in another bedgraph file.
## author: Yaping Liu  lyping1986@gmail.com


use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl bedGraphMinusAnother.pl [Options] output.bedgraph file1.bedgraph file2.bedgraph \n\n";

    print STDERR " value in one begraph minus the value in another bedgraph file\n";
	
	print STDERR "  [Options]:\n\n";
	print STDERR "  --keepA : when specified, it will keep all file1's line, even there is no line overlapped in file2, bu minus 0 to the score field.\n\n";
	print STDERR "  --overlap : when specified, it will only minus file2's score field in the overlapped regions of file1. not overlapped region will be kept the same, and single region in file1 will be split into overlapped/not-overlapped regions...\n\n";

	exit(1);
}

my $keepA="";
my $overlap="";
GetOptions( 
	"keepA" => \$keepA,
	"overlap" => \$overlap,
);

if(scalar(@ARGV) < 3){
	&usage();
}

my $out=$ARGV[0];
my $file1=$ARGV[1];
my $file2=$ARGV[2];
open(OUT,">$out") or die "can not open file $out:$!\n";
if($overlap ne ""){
	&overlap_mode($file1,$file2);
	close(OUT);
	exit(0);
}

if($keepA ne ""){
	my $tmp=$file1;
	$file1=$file2;
	$file2=$tmp;	
}
#print "file1:$file1\n";
#print "file2:$file2\n";



my @chrs=(1..22,"X","Y");
for(my $i=0; $i<=$#chrs;$i++){
	my %f1_hash=();
	my $chr="chr".$chrs[$i];
	#print "$chr\n";
	my $cmd="grep \"^$chr\\b\" $file1 | ";
	#print STDERR "$cmd\n";
	open(FH1,$cmd) or die "can not open file $file1:$!\n";
	while(<FH1>){
		chomp;
		next if($_=~/^#/ or $_=~/^track/);
		my @f=split "\t";
		my $key=$f[0]."-".$f[1]."-".$f[2];
		my $value=$f[3];
		$f1_hash{$key}=$value;
		#print "$key:$value\n";
		
	}
	close(FH1);
	$cmd="grep \"^$chr\\b\" $file2 | ";
	#print STDERR "$cmd\n";
	open(FH2,$cmd) or die "can not open file $file2:$!\n";
	while(<FH2>){
		chomp;
		next if($_=~/^#/ or $_=~/^track/);
		my @f=split "\t";
		my $key=$f[0]."-".$f[1]."-".$f[2];
		my $value=$f[3];
		#print "$key:$value\n";
		if(exists($f1_hash{$key})){
			my $minus=$f1_hash{$key}-$value;
			if($keepA ne ""){
				$minus=$value-$f1_hash{$key};
			}
			#print STDERR "$f[0]\t$f[1]\t$f[2]\t$minus\n";
			print OUT "$f[0]\t$f[1]\t$f[2]\t$minus\n";
		}
		elsif($keepA ne ""){
			my $minus=$value;
			print OUT "$f[0]\t$f[1]\t$f[2]\t$minus\n";
		}
		
		
		
	}
	close(FH2);
}
close(OUT);

#fileB's region should not have any overlap, and should be sorted!!
sub overlap_mode{
	my $fileA=shift @_;
	my $fileB=shift @_;
	my $cmd="intersectBed -loj -sorted -a $fileA -b $fileB | ";
	my $chr_pre="";
	my $start_pre=-1;
	my $end_pre=0;
	my @same_lines_in_A=();
	open(FH,$cmd) or die "can not open file $fileA or $fileB:$!\n";
	while(<FH>){
		chomp;
		next if($_=~/^#/ or $_=~/^track/);
		my @f=split "\t";
		if($f[0] eq $chr_pre && $f[1] eq $start_pre && $f[2] eq $end_pre){
			
		}
		else{
			
			##sort array first, bedtools have some bug here, not completed sort for the output result..
			if(scalar(@same_lines_in_A)>1){
				@same_lines_in_A = sort { (split "\t", $a)[5] <=> (split "\t", $b)[5] } @same_lines_in_A;
			}
			
			
			##extract overlap & no_overlap region
			
			my $preSegEnd=-1;
			foreach my $same_line_in_A(@same_lines_in_A){
				#print OUT "line: $same_line_in_A\t $chr_pre \t $start_pre \t $end_pre\n";
				my @inlines=split "\t",$same_line_in_A;
				if($inlines[5] eq "-1"){##if no overlap at all..
					print OUT "$inlines[0]\t$inlines[1]\t$inlines[2]\t$inlines[3]\n";
					last; 
				}
				
				my @overlap_interval=&overlap_region($inlines[1],$inlines[2],$inlines[5],$inlines[6]);
				my $minus;
				if($preSegEnd != -1){
					if($preSegEnd < $overlap_interval[0]){
						$minus=$inlines[3];
						print OUT "$inlines[0]\t$preSegEnd\t$overlap_interval[0]\t$minus\n";
					}
				}
				$minus=$inlines[3]-$inlines[7];
				print OUT "$inlines[0]\t$overlap_interval[0]\t$overlap_interval[1]\t$minus\n";
				$preSegEnd=$overlap_interval[1];
			}
			
			##clear array tmp
			$chr_pre = $f[0];
			$start_pre = $f[1];
			$end_pre = $f[2];
			@same_lines_in_A=();
			
		}
		push(@same_lines_in_A,$_);
		
	}
	close(FH);
	
	
}

sub overlap_region{
	my $startA=shift @_;
	my $endA=shift @_;
	my $startB=shift @_;
	my $endB=shift @_;
	my $overlap_start=$startA > $startB ? $startA: $startB;
	my $overlap_end=$endA > $endB ? $endB: $endA;
	
	return ($overlap_start, $overlap_end);
	
}


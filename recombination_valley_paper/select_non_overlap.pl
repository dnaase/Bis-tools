my $input=$ARGV[0];
my $random=int(rand(10000000));
my $tmp = "tmp.$random";
#`intersectBed -loj -a $input -b $input > $tmp`;
my %duplicate_hash=();
my %unique_tmp_hash=();
my $tmp_key="";
open(F,"intersectBed -loj -a $input -b $input |") or die;
while(my $in=<F>){
	chomp($in);
	my @f=split "\t",$in;
	my $key1=$f[3];
	my $key2=$f[8];
	
	if($tmp_key eq ""){
		$tmp_key = $key1;
		$unique_tmp_hash{$key2}=$in;
		
	}elsif($tmp_key eq $key1){
			if(! exists $duplicate_hash{$key2}){
				$unique_tmp_hash{$key2}=$in;
			}

	}else{
		if (%unique_tmp_hash) {
			my $randomkey = choose( keys %unique_tmp_hash );
			#print STDERR "$key1\t$key2\t$randomkey\n";
			##put all other keys into duplicated hash
			foreach my $dup_key(keys %unique_tmp_hash){
				if($dup_key ne $randomkey){
					$duplicate_hash{$dup_key}=1;
					#print STDERR "$dup_key\t";
				}
			}
			#print STDERR "$dup_key\n";
			##print out the non-overlapped key
			my @r = split "\t",$unique_tmp_hash{$randomkey};
			print "$r[5]\t$r[6]\t$r[7]\t$r[8]\t$r[9]\n";
			##put this key also into duplicated hash
			$duplicate_hash{$randomkey}=1;
		}
		%unique_tmp_hash=();
		$tmp_key = $key1;
		if(! exists $duplicate_hash{$key2}){
			$unique_tmp_hash{$key2}=$in;
		}
		
		
	}
	
}
close(F);

#`unlink $tmp`;

sub choose {
  my $count = 0;
  my $result = undef;
  for( @_ ) {
    $result = $_ if rand() < 1 / ++$count;
  }
  $result;
}


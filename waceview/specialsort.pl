use strict vars ;

my $n1 = 0 ;
my @aa ;
my @cc ;
my %gg ;

while (<>)
  {
    chomp () ;
    my $line = $_ ;
    next unless ($line =~ /REMOTE_ADDR/ ||
		$line =~  /Server normal exit/ ||
		$line =~  /^#### Server starts /  ||
		$line =~  /^....-..-.._..:..:..New Client / ||
		$line =~  /^....-..-.._..:..:.. New Client / ||
		$line =~  /^....-..-.._..:..:.. Closing Client / ||
		$line =~  /PROXIED_IP/  ||
		$line =~  /HTTP_HOST/
		) ;
    next unless ($line =~ /^([0-9][0-9][0-9][0-9])-([0-9][0-9])-([0-9][0-9])_([0-9][0-9]):([0-9][0-9]):([0-9][0-9])/) ;
    my ($y,$m,$d,$h,$mn,$sec) = ($1,$2,$3,$4,$5,$6) ;
    my $nn = (400*($y-2000) + $m * 32 + $d)*100000 + 3600*$h + 60 * $mn + $sec ;
    
    next if ($nn < 0) ;
    $gg{$n1} = $nn;
    $aa[$n1] = $line ;
    $cc[$n1] = $n1 ;
    $n1++ ;
    #print "$nn $y $m $d $h $mn $sec $line\n" ;
  }
#die "ok\n" ;
my @bb = sort { $gg{$a} <=> $gg{$b} } @cc ;

for (my $i = 0 ; $i < $n1 ; $i++)
  { my $j = $bb[$i] ; print "$aa[$j]\n" ; }


#!/bin/env perl

# @(#)client.pl	1.2 5/24/97 
#from Programming Perl by Larry Wall and Randall Schwartz
#O'Reilly and Associates, Inc 1991
#pp 342-344
# and modified

$| = 1; #turns off output buffering for STDOUT

chop ($uname = `uname`);
chop ($version = `uname -r`); $version =~ s/v//ig; # just want a number

#########################################
# edit this section if client.pl doesn't work on your computer
#
### first try uncommenting the next 3 lines if you're using perl 5
#use Socket;
#$AF_INET = &AF_INET;
#$SOCK_STREAM = &SOCK_STREAM;
###
#
### otherwise, comment out the previous section, and add 
#   values for your computer in the section below.
#   To get these values, check sys/socket.h
#   probably in /usr/include

unless ($AF_INET) { # bypassed if the 'use Socket' section works
    if (($uname =~ /IRIX/) || 
	     (($uname eq 'SunOS') && ($version >= 5))) {
	# works for IRIX64 and SOLARIS (= SunOS 5.x)
	$SOCK_STREAM = 2;
	$AF_INET = 2;
    }  else {
	# default
	# works for SunOS 4, AIX, HP, OSF1
	$SOCK_STREAM = 1;
	$AF_INET = 2;
    }
}

### don't edit past here
#################################################

print STDERR "client.pl running...\n";
($them, $port) = @ARGV;
$port = 2345 unless $port;
$them = 'localhost' unless $them;

$SIG{'INT'} = 'dokill';
sub dokill {
    kill 9,$child if $child;
}

$sockaddr = 'S n a4 x8';

chop($hostname = `hostname`);

($name,$aliases,$proto) = getprotobyname('tcp');
($name,$aliases,$port) = getservbyname($port, 'tcp')
    unless $port =~ /^\d+$/;;
($name,$aliases,$type,$len,$thisaddr) =
    gethostbyname($hostname);
($name,$aliases,$type,$len,$thataddr) = gethostbyname($them);

$this = pack($sockaddr, $AF_INET, 0, $thisaddr);
$that = pack($sockaddr, $AF_INET, $port, $thataddr);

# Make the socket filehandle.

if (socket(S, $AF_INET, $SOCK_STREAM, $proto)) {
    print STDERR "socket ok\n";
}
else {
    die $! . "\nPerhaps \$AF_INET and \$SOCK_STREAM are not set correctly for your system;\nEdit wcripts/client.pl\n";
}

# Give the socket an address.

if (bind (S, $this)) {
   print STDERR "bind ok\n";
}
else {
    die $!;
}

# Call up the server.

if (connect(S, $that)) {
    print STDERR "connect ok\n";
}
else {
    die $!;
}

# Set socket to be command buffered.
# unbuffer!
select(S); $| = 1; select(STDOUT);
#select(S); undef $|; select(STDOUT);

# Avoid deadlock by forking.
# try switching parent and child

if($child = fork) {
    while(<S>){
	print;
#	exit if ($_ eq "201 A bientot\n");
    }
    print STDERR "client.pl done.\n";
}
else
{
    while(<STDIN>)
    {
	#print STDERR "xxx".$_ ;
	print S ;
    }
    #print STDERR "xxxfin" ;
    print S "__22__\n" ; 
}
# if($child = fork) {
#     while(<STDIN>){
# 	print S;
#     }
# #    sleep 100;
# #    do dokill();
# }
# else {
#     while(<S>){
# 	print;
#     }
# }


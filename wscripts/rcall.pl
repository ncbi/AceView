#!/bin/env perl

# @(#)rcall.pl	1.1 5/24/97
# There's a program called expect that you can use to wrap telnet into a
# script.  It's part of the TCL distribution.  However, I'd do the whole
# thing in perl.
# 
# Here's a small perl client that you can use to communicate with a
# server on a port somewhere.  You could open a pipe to it to handle
# all the socket communications, and capture its output to a file.
# 
# If you try the script and it complains that it can't find
# sys/socket.ph, you need to run the perl h2ph program to perl-ize the
# system's /include/sys/socket.h file and install the .ph file in
# /usr/local/lib/perl5/sys/.
# 
# Good luck,
# 
# Lincoln Stein, 
 
$| = 1; #turns off output buffering for STDOUT.

require "sys/socket.ph";

$usage = 'Usage: socketload host port [input files]\n';

$remote_host = shift || die $usage;
$port = shift || die $usage;

# make sure we die when we're supposed to
$SIG{'KILL'} = 'dokill';
$SIG{'TERM'} = 'dokill';

sub dokill {
        kill 'TERM',$child if $child;
}

sub dochild {
    exit 0;
}

sub dopipe {
    warn "Connection closed by foreign host.\n";
    exit 0;
}

# open up a socket to the port
$sockaddr = 'S n a4 x8';
chop($hostname = `hostname`);

($name,$aliases,$proto) = getprotobyname('tcp');
($name,$aliases,$type,$len,$thisaddr) = gethostbyname($hostname);

if ($remote_host =~ /^(\d+)\.(\d+)\.(\d+)\.(\d+)/) {
    $thataddr = pack(CCCC,$1,$2,$3,$4);
} else {
    ($name,$aliases,$type,$len,$thataddr) = gethostbyname($remote_host);
    $name || die "Can't find address of $remote_host";
}


$this = pack($sockaddr, &AF_INET, 0, $thisaddr);
$that = pack($sockaddr, &AF_INET, $port, $thataddr);

# Make the socket filehandle.

socket(S, &AF_INET, &SOCK_STREAM, $proto) || die $!;

# give the socket an address

bind(S, $this) || die $!;

# Call up the server.

connect(S,$that) || die "$remote_host refused connection for port $port.\n";

# set socket to be command buffered.

select(S); $| = 1; select(STDOUT);


# avoid deadlock by forking

if ($child = fork) {
    $SIG{'CHLD'} = 'dochild';
    $SIG{'PIPE'} = 'dopipe';
    while (<>) {
        print S "$_";
    }
    print S "\cD\n";
    wait;

} else {
    $SIG{'PIPE'} = 'dopipe';
    while (<S>) {
        print "$_";
    }
}

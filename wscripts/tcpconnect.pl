#!/bin/env perl
#
# $Id: tcpconnect.pl,v 1.1.1.1 2002/07/19 20:23:33 sienkiew Exp $
#
# tcpconnect: make a connection to remote port by tcp
#
# Copyright (c) 1991,1992 Kazumasa Utashiro <utashiro@sra.co.jp>
# Software Research Associates, Inc., Japan
#
#
# With this command, you can forward some specific protocol
# on gateway machine which doesn't forward IP packet.
#
# Syntax:
#      tcpconnect server[:service]
#      tcpconnect server[:service] client
#      tcpconnect server[:service] client server[:service] client ...
#      tcpconnect -f config-file
#
# + If client is not specified or is '*', any client is allowed
#   to connect.
#
# + Server is choosed on first-hit policy.  So, if there is '*'
#   client at the first, remained servers will not be checked.
#
# + Same port number will be used when no server services is
#   specified.
#
# + Client and server name can be symbolic hostname or IP
#   address in dotted decimal notation.
#
# + Service name can be specified by symbolic name or port-number
#
# + Configuration file contains client and server pair on each line.
#
# EXAMPLES:
#
# All SMTP connection is to be forwarded to mail-server:
#
#      smtp stream tcp nowait root /etc/tcpconnect tcpconnect mail-server
#
# NNTP connectioned is exchanged between A and B:
#
#      nntp stream tcp nowait root /etc/tcpconnect tcpconnect A B B A
#
# Simply connect a terminal to nntp server:
#
#      % tcpconnect nntp-server:nntp
#

# require 'sys/socket.ph';
unless (do 'sys/socket.ph') {
    #print "File sys/socket.ph is not found. Using default...\n";
    eval 'sub SOCK_STREAM {1;} sub AF_INET {2;} sub PF_INET {2;}';
}

while ($_ = $ARGV[0], /^-/) {
    shift;
    if (/-f$/)          {$configfile = shift || &usage; next;}
    if (/-d(\d*)$/)     {$debug = $1||1;                next;}
    if (/-h/)           {&usage;                        next;}
    &usage;
}

sub usage {
    $0 =~ s|.*/||;
    $* = 1;
    ($usage = <<"    EOF") =~ s/^    //g;
    Usage: $0 server[:service] [client] [server[:service] client ...]
           $0 -f config-file

    ($rcsid)
    EOF
    for (@_) { print "ERROR: $_", /\n$/ ? "" : "\n"; }
    print $usage;
    exit(1);
}

unless ($configfile) {
    @forwardlist = @ARGV;
} else {
    open(CONF, $configfile) || die("$configfile: $!\n");
    while(<CONF>) {
        chop;
        s/#.*$//;
        s/^\s*//;
        next if /^$/;
        push(@forwardlist, (split($_))[0,1]);
    }
    close(CONF);
}

$sockaddr='S n a4 x8';
($name, $aliases, $TCP) = getprotobyname('tcp');

chop($localname = `hostname`);

unless ($hersockaddr = getpeername(STDIN)) {
    $server = shift || &usage;
} else {
    open(STDERR, ">/dev/console");
    select(STDERR); $| = 1; select(STDOUT);

    ($family, $herport, $heraddr) = unpack($sockaddr, $hersockaddr);

    $mysockaddr = getsockname(STDIN);
    ($family, $myport, $myaddr) = unpack($sockaddr, $mysockaddr);

    if ($debug) {
        printf STDERR ("$0: Connection from %s(%d)\n",
                       &dotted($heraddr), $herport);
    }
    CHECKADDR: {
        while (($server, $client) = splice(@forwardlist, 0, 2)) {
            if (!defined($client) || $client eq '*') {
                $clientaddr = "\0\0\0\0";
            } else {
                ($clientaddr = &getaddr($client)) || next;
            }
            if ($clientaddr eq "\0\0\0\0" || $clientaddr eq $heraddr) {
                $server = "$server:$myport" unless ($server =~ /:/);
                last CHECKADDR;
            }
        }
        printf STDERR ("Connection from %s is not allowed!\n",
                        &dotted($heraddr));
        exit(1);
    }
}

($servername, $serverport) = split(/:/, $server);
$serverport || &usage("No server port");

($serveraddr = &getaddr($servername)) || die "Unknown server $servername.\n";
$serverport = (getservbyname($serverport, 'tcp'))[2]
    unless $serverport =~ /^\d+$/;

$that = pack($sockaddr, &AF_INET, $serverport, $serveraddr);
socket(S, &PF_INET, &SOCK_STREAM, $TCP) || die "socket: $!";
connect(S, $that) || die "connect: $!";
select(S); $| = 1; select(stdout);

if ($child = fork) {
    &forward(S, STDOUT);
} else {
    &forward(STDIN, S);
}

print STDERR "$0($$): exiting\n" if $debug;
exit(0);

sub forward {
    local($from, $to) = @_;
    select($from); $| = 1;
    select($to); $| = 1;
    if (-t $from) {
        eval "print $to \$_ while(<$from>);";
    } else {
        eval "print $to \$_ while(read($from, \$_, 4096));";
    }
    shutdown($from, 1); shutdown($to, 0);
}

sub getaddr {
    local($_) = @_;
    /^[0-9\.]+$/ ? pack("C4", split(/\./)) : (gethostbyname($_))[4];
}

sub dotted { join('.', unpack('C4', shift)); }



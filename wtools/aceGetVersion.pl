#!/usr/local/bin/perl
#
# Perl script to read the acedb kernel version from the source file that
# holds the version number:   aceversion.c
#
# The idea is that only the version.c source file will be altered when changing
# the ACEDB version number, then the scripts that produce the public tar files
# will use this script to get that version number.
#
# The script returns the version, release and update read from aceversion.c
# in the form:    <version>_<release><update>     e.g.  4_6d
#
# it must be given the name of the directory where aceversion.c can be found.
#

# Initialise...
$debug = 0 ;						    # Debug flag, may be reset from command line...

($cmd = $0) =~ s#.*/## ;				    # Get this commands basename (got this
							    # magic from 'Learning Perl', p. 214.).

$msghdr = $cmd . ': ' ;					    # Set up a prefix for error messages.

$source_file = 'aceversion.c' ;				    # The ACEDB source file to be searched.


# Parse command line options:
#                             -d  turns on debugging
#
require "getopts.pl" ;
&Getopts('d') ;
if ($opt_d == 1) { $debug = $opt_d ;}


# We don't want any arguments to the script.
#
($debug) && print $msghdr . "program arguments - @ARGV\n" ;
($#ARGV < 0 || $#ARGV > 0) && die $msghdr . "Wrong number of parameters supplied, specify JUST the directory where aceversion.c can be found.\n" ;


# This is the location of the aceversion.c file we want to read.
#
$input_file = "@ARGV[0]/$source_file" ;

# Now check out the input file to see if it's OK...
(-f $input_file && -r $input_file) || die $msghdr . "file specified is not a file, or is not readable by this program.\n" ;

# Open the input file...
open(INFILE, $input_file) || die $msghdr . "could not open source code file\n" ;

# Now start parsing the input file...
while (<INFILE>)
  {
  ($debug) && print $msghdr . "Data just read in $&\n" ;

  # Once we find the start tag we know that the version, release & update are on the next
  # three lines.
  if ( /ACE_VERSION_START/ )
    {
    $version_line = <INFILE> ;
    ($junk, $junk, $version, $junk) = split(/\s+/, $version_line) ;
    $release_line = <INFILE> ;
    ($junk, $junk, $release, $junk) = split(/\s+/, $release_line) ;
    $update_line = <INFILE> ;
    ($junk, $update, $junk) = split(/"/, $update_line) ;
    }

  }

# Don't need the file anymore so close it...
close(INFILE) || die $msghdr . "could not close source file\n" ;

($debug) && print $msghdr . "File successfully processed...\n" ;

# Return the result as a single string sent to STDOUT in the form "4_5e"
print "$version\_$release$update\n" ;

exit 0 ;


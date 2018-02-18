#  File: ace_perl.pl
#   Author: Jean Thierry-Mieg, Steve Rosen, Gabor Marth
# -------------------------------------------------------------------
#  This file is part of the client/server interface to
#       the ACEDB genome database package, written by
#  	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
#	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
#   Copyright (C) J Thierry-Mieg and R Durbin, 1989-
#
#  Description:
#  Exported functions:
#   ace_in, ace_new, ace_count, ace_table, ace_script
# -------------------------------------------------------------------

#  $Id: ace_perl.pl,v 1.1.1.1 2002/07/19 20:23:38 sienkiew Exp $ 

#################################################################################
####################### Perl acedb client tool box ##############################
## 
## Wash U, Jan 12, 1995
#################################################################################

$user_name = (getpwuid($<))[0];
$who_and_when = "who \"$user_name $0\"\nwhen now";

#################################################################################

#usage ace_in (ace file commands)
#example which parses data into acedb
#ace_query(<< END
#        Subclone $key
#        Estimated_length $sequence_length
#        $who_and_when
#        State Waiting_for_enough_picking
#END ) ;
#        , "-ace_out -ace_in" ) ;

sub ace_in {
    open(ACECLIENT, 
      "echo '$_[0]' | /home1/crick/mieg/ace/acl crick -ace_out -ace_in|") ||
       die $@;
    local(@r) = <ACECLIENT>;
    die $@ if $@;
    return @r;
}

#################################################################################

#usage ace_new (Class, format_with_%d_inside)
#example $new_name = ace_new ("T_subclone", "toto%d") 
sub ace_new {
    local ($class, $format) = @_ ;
    open(ACECLIENT, "echo New $class $format | /home1/crick/mieg/ace/acl crick |") || die $@;
    local(@r) = <ACECLIENT>;
    die $@ if $@;
    for $rr (@r)
     { 
	if ($rr =~ /$class/) {
          local ($uu, $xx) = split (' ', $rr) ;
	  return $xx ;
	  }
      }
}

#################################################################################

#usage ace_count (Class)
#example ace_count ("In_Line")  
#   ... where In_Line is a subclass declared in wspec/subclasses.wrm
sub ace_count {
    open(ACECLIENT, "echo Find '$_[0]' | /home1/crick/mieg/ace/acl crick |") || die $@;
    local(@r) = <ACECLIENT>;
    die $@ if $@;
    for $rr (@r)
	{ $rr =~ /Found (.*) objects in this class/;
	  if ($1 > 0)  { return $1 ;}
        }
}

#################################################################################

#usage @table = ace_table ("Table.def")
sub ace_table {
    local ($table) = @_ ;
    open(ACECLIENT, "echo Table $table | /home1/crick/mieg/ace/acl crick -ace_out |") || die $@;
    local(@r) = <ACECLIENT>;
    die $@ if $@;

    local (@reply) ;
    for $rr (@r)
     { chop ($rr) ;
       if (! ($rr eq ""))
         {  $rr = $rr . "\n" ;
	     push (@reply, ($rr)) ;
         }
      }
  @reply ;
}

#################################################################################

#usage ace_query (ace_script, command line parameters)
# example of a script
#ace_query(<<END
#        Find Subclone
#        Follow T_subclone
#        Edit State Waiting_for_enough_Mapping
#        Show
#END
#        , "-ace_out" ) ;

sub ace_script {
    open(ACECLIENT, "echo '$_[0]' | /home1/crick/mieg/ace/acl crick $_[1] |") || die $@;
    local(@r) = <ACECLIENT>;
    die $@ if $@;
    return @r;
}

#################################################################################
################## end of perl acedb client tool box ############################
#################################################################################


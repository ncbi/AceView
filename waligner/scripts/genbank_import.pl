#!/usr/local/bin/perl
# $Header: /home/mieg/aa/CVS_ACEVIEW/ace/waligner/scripts/genbank_import.pl,v 1.3 2014/04/10 20:51:15 mieg Exp $

# this was once GenEmbl2ace.pl, but not any more - it is revised
# to more correctly generate the objects that we use at NCBI
#  Mark S. 10/2002

############################# CONFIGURABLES ############################

# These may have been configurable to the original author, but don't
# count on being able to change them any more.  Any other configuration
# is not tested with the current import system.  If you break it, you
# can keep both pieces.
#  Mark S. 10/2002

# use  either database accession number or identifier (LOCUS in GENBANK)
# as basis for ?Sequence object name

# include database record as LongText

$LongText = 'YES';
$LongText = 'NO';

######################### END OF CONFIGURABLES #########################

#
# typical perl-ism:  This seems to mean that the "line seperator" is
# a newline slash slash newline.  This matches a blank comment that
# appears between record in the genbank output.  After we do this,
# all the "operate on line" operations that perl has will actually
# operate on records.

$/="\n//\n";

open(LTFILE, ">longtext.ace");

#
# chitchat_count is so we can print a count of records every thousand
# records processed.  This is just because it normally takes so long
# that you would like some indication that it isn't hung.

$chitchat_count = 0;

#
# Here are some parameters from the bigger system.  They are 
# environment variables.

#
# want_import is either 
#	"NM" if you want NM 
#	"NNM" if want things that are not NM
#

$want=$ENV{want_import} ;

#
# These 3 lines are data for the Imported tag in the sequence object
#

$date_marker=$ENV{date_marker} ;
$import_marker=$ENV{import_marker} ;
$query_marker=$ENV{query_marker} ;
$active_species==$ENV{active_species} ;

#
# tag_insert_[1-3] are tags that are placed on each Sequence object
# created, no matter what
#

$tag_insert_1=$ENV{tag_insert_1} ;
$tag_insert_2=$ENV{tag_insert_2} ;
$tag_insert_3=$ENV{tag_insert_3} ;

$accept_join=$ENV{accept_join} ;
#
# name_prefix is a string stuck on the front of every Sequence object name
# Initially it is either "GB:" or ""
#

$name_prefix=$ENV{name_prefix} ;

#
# objectName can be one of 
#	LOCUS
#	VERSION
#	ACCESSION

$objectName=$ENV{objectName} ;

#
# want_remark is "y" if we want to collect long_remark fields
#
$want_remark=$ENV{want_remark} ;

#
#
#
if ( $want eq "NM" ) { }
elsif ($want eq "NNM") { }
else	 {
	print "setenv want_import { NM | NNM }\n";
	exit 1;
	}


#
# exit with this code at the end.

$exitcode=0;

#
# This is a typical perl main program, but remember that <> is
# now operating on records instead of lines.

RECORDLOOP:
while (<>){
    /^LOCUS/ || die "Entry does not start with LOCUS line: $_" ;

    /(^LOCUS)\s+(\S+)/ || die "Entry does not have an identifier: $_";
    $locus = $2;


    #
    # remind the user that we are not hung

    $chitchat_count = $chitchat_count + 1;
    if (($chitchat_count % 1000) == 0)
	{
	printf(STDERR "%d records\n",$chitchat_count);
	}

    #
    # we have not yet identified this sequence as being bad

    $bad_sequence = 0;

    #
    # Here, it puts SPLITHERE in between various fields
    # of the record.  Later, it will break the record into
    # fields using SPLITHERE.  I think the whole point is
    # to try to use the same split() call for GB or EM.
    # We don't use EMBL any more, so this may go away.

    s/\n(\S+)/\nSPLITHERE\n$1/g;

    # print "vvvv $_\n";

    # $text to store full text of record
    $text = '';

    # $papers to concatenate paper objects
    $papers = '';

    # blank out some variables that should not be inherited
    # from record to record

    $direction = '';

    #
    # iterate over the fields of the record.  Be sure
    # to say  next;  at the end of the code for your field.

    for (split(/SPLITHERE\n/)){

	if (/^ORIGIN/){
            s/ //g ;
            s/\d//g ;
            s/\/\/\n// ;
            s/U/T/g ;
            s/^.*\n// ;
            $seq = $_ ;
	    next;
        }

        if ( $LongText == 'YES' )
	    { $text .= $_ ; }

	if (/(^ACCESSION)\s+(.*)/){
            @ac = split(/\;| /,$2);
	    $ac = shift(@ac);
            #
            # I would explain this better if I really knew what it was about.
            # There are NM records and "not NM" records.  An old version of 
            # the import would just grab everything and then purge them from 
            # the database later.  I just have you tell this script which you
            # want, and don't process the other ones.
         if (1)
	   {
	     if ($ac =~ /^[NX][MR]_*/ )
	       {
		 if ( $want ne "NM") { next RECORDLOOP; }
	       }
            else
	      {
	        if ( $want eq "NM") { next RECORDLOOP; }
	      }
	     next;
	   }
	}

	if (/(^COMMENT)\s+(.*)/){
		if ($want_remark eq "y")  {
			$a = $_;
			$a =~ s/COMMENT/       / ;
			printf("Long_remark sequence_%s\n", $locus);
			printf( LTFILE "LongText sequence_%s\n", $locus);
			printf( LTFILE "%s\n",$a);
			printf( LTFILE "***LongTextEnd***\n\n");
			next;
		}
	}

	if (/(^VERSION)\s+([^\s]*)\s+(.*)/){
	    $ncbiGI = $3 ;
            @versiona = split(/\;| /,$2);
	    $ncbiversion = shift(@versiona);
	    next;
	}

	if (/(^NID)\s+(\w+)/){
	    $ni = $2;
	    next;
	}

        if (/^DEFINITION/){
	  s/DEFINITION  //;
	  s/            //g;
	  chomp ;
	  s/\n/ /g ;
	  $title = $_;

	  # direction is only used for EST - other records
	  # have titles that say 5' or 3' but are not intended
	  # to indicate direction

	  $direction = "Forward" if ($title =~ /cDNA clone .* 5\'/) ;
	  $direction = "Reverse" if ($title =~ /cDNA clone .* 3\'/) ;

	  $isCDNA = 1 if ($title =~ /mRNA/) ;
          $isComplete = 1 if ($title =~ /complete cds/) ;
	  next;
        }

	if (/^KEYWORDS/) {
	    s/KEYWORDS    //;
	    s/\n\s+/ /g;
	    s/^\s+//;
	    s/\n$//;
	    s/\.$//;
	
	    # object name locus or accession ?
            # default
	    $nam = $ac ;
	    $nam = $name_prefix . $nam;

            if  ($locus && ($objectName eq 'LOCUS'))
	      {
		$nam = $locus ;
		$nam = $name_prefix . $nam;
	      }
	    elsif  ($objectName eq 'ACCESSION')
	      {
		$nam = $ac ;
		$nam = $name_prefix . $nam;
	      }
	    elsif  ($ncbiversion && ($objectName eq 'VERSION'))
	      {
		$nam = $ncbiversion ;
	      }

	    $dbnam = 'GENBANK';
	    $nPaper = 0 ;
	    # start output
	    print  "\nSequence : $nam\n" ;
	    undef $geneidNam ;
            if ($date_marker) { print  "Imported $date_marker \"$import_marker\" \"$query_marker\"\n" ;}
	    else { print  "Imported 2006-09-07\n" ; }
	    print  "Not_curated\n";
	    print  "Title \"$title\"\n" ;

	    if ($title =~ /mutant/i ) 
		{ print "Is_mutant\n"; }
	    elsif ($title =~ /complete/i)
		{ print "Is_mrna\n"; }
	    elsif ($title =~ /partial/i)
		{ print "Is_partial\n"; }

	   if ( $tag_insert_1 )
		{ printf("%s\n",$tag_insert_1); }
	   if ( $tag_insert_2 )
		{ printf("%s\n",$tag_insert_2); }
	   if ( $tag_insert_3 )
		{ printf("%s\n",$tag_insert_3); }

	    print "cDNA\nIs_read\n" if ($isCDNA) ; undef $isCDNA ;
	    print "Complete_CDS\n" if ($isComplete) ; undef $isComplete ;
	    print  "DB_annotation $dbnam $nam\n" if ($LongText eq 'YES');
            print  "Database $dbnam \"$ncbiversion\" \"$ncbiGI\"\n";
	    # mieg print  "Identifier $locus\n";
	    # print  "Locus $locus\n"  if ($locus) ; pas correct ?

	    if ( $direction eq "Forward" ) {
		print "Forward\n";
		print "aForward\n";
	    }

	    if ( $direction eq "Reverse" ) {
		print "Reverse\n";
		print "aReverse\n";
	    }

	    for(@ac){
		print "!Accession $_\n";
		print  "Database Genbank \"$ncbiversion\" \"$ncbiGI\" \n";
	    }

	    print  "!Nucleotide_id $ni\n";

	    $have_est = 0;

	    if ( $keyWords eq 'YES' ) {

	        for(split(/\; /)){
	 	    if ( $_ eq "EST" ) 
			{
			$have_est = 1;
			}
		    else
			{
		        print "Keyword \"$_\"\n";
			}
	        }

	        if ( $have_est == 0 ) {
			print "IS_est\n";
		}
 	    }

	    next;
	}

	if (/^SOURCE/) {
	    /ORGANISM  (\w+ \w+)/;

	    print  "Species \"$1\"\n" ;

	    next;
        }

        if (/^REFERENCE/) {
	    &Paper($nam,$_);
	    next;
	}

	if (/^FEATURES/) {

	    # insert dividers for splitting into features
            s/\n\s{5}(\S+)/ZZZZ$1/g;

	    #print "uuuuu $_\n" ;

	    undef $haveExons ;
	    $nsubs = 0 ;

	    # split into features
	    for (split(/ZZZZ/)) {

		#
		# the actual FEATURES line is iterated across first because 
		# of the way we split the features section.  Skip this
		# line
		/^FEATURES/ && next ;

		#
		#
		chomp ;

		s/\n\s{21}\//ZZZZ/g; # insert dividers for splitting into qualifiers
		s/\n\s{21}//g;       # substitute excess space

		(s/(\S+)\s+// && ($key = $1)) || die "No FT key in $locus\n" ;


		@quals = split (/ZZZZ/) ; # split into qualifiers
		($loc = shift (@quals)) || die "No loc in $id - $key\n" ; # pull out location


		if ($key eq "source")
		  {
		    #print "// key=$key\n" ;
		    foreach $qual (@quals)
		      {
			#print "// qual $qual\n" ;
			if ($qual =~ /^clone=/)
			  {
			    $qual =~ s/clone=// ;
			    $qual =~ s/\n/ / ;
			    # $qual =~ s/ 5\'// ;
			    $qual =~ s/Kohara clone// ;
			    $qual =~ s/from Y\. Kohara// ;
			    $qual =~ s/obtained from Yuji Kohara// ;
                            ($clone) = split(/\;/,$qual) ;
			    $clone =~ s/\"//g ;
			    $clone =~ s/\n/ /g ;
			    $clone =~ s/^\s*//g ;
			    chomp ($clone) ;

			    $good_clone = 0;

			    if (index($clone,"IMAGE:") == 0)
				{ $good_clone=1; }

			    if ($clone eq "IMAGE:")
				{ $good_clone=0; }

			    $title_subset = "cDNA clone " . $clone . " 3'";
			    if (index($title, $title_subset) >= 0)
				{ $good_clone=1; }

			    $title_subset = "cDNA clone " . $clone . " 5'";
			    if (index($title, $title_subset) >= 0)
				{ $good_clone=1; }

			    if ($good_clone)
				{ printf("cDNA_clone \"%s\"\n",$clone); }
			    else
			   	{
				# printf( STDERR "bad clone name %s in %s\n",$clone, $nam);
				}
			    printf("cDNA_clone_GB \"%s\"\n",$clone); 

			    next ;
			  }
			if ($qual =~ /^dev_stage=/)
			  {
			    $qual =~ s/dev_stage=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Dev_stage \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^clone_lib=/)
			  {
			    $qual =~ s/clone_lib=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Library \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^tissue_type=/)
			  {
			    $qual =~ s/tissue_type=// ;
			    $qual =~ s/\"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Tissue \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^dev_stage=/)
			  {
			    #$qual =~ s/dev_stage=// ;
			    $qual =~ s/\"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Remark \"dev_stage = $qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^lab_host=/)
			  {
			    $qual =~ s/lab_host=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "!Lab_host \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^note=/)
			  {
			    $qual =~ s/note=// ;
			    $qual =~ s/\n/ / ;
			    $qual =~ s/\"//g ;
			    @notesa = split (/\;\s*/,$qual) ;
			    foreach $note (@notesa)
			      {
				#print "// note=$note\n" ;
				if ($note =~ /^Organ: /)
				  {
				    $note =~ s/^Organ: // ;
				    print "Tissue \"$note\"\n" ;
				  }
				elsif ($note =~ /^Vector: /)
				  {
				    $note =~ s/^Vector: // ;
				    print "Seq_vec \"$note\"\n" ;
				  }
				elsif ($note =~ /^Site_1:/)
				  {
				    $note =~ s/^Site_1:// ;
				    print "!Site_1 \"$note\"\n" ;
				  }
				elsif ($note =~ /^Site_2:/)
				  {
				    $note =~ s/^Site_2:// ;
				    print "!Site_2 \"$note\"\n" ;
				  }
				else
				  {
				    print "Note_gb \"$note\"\n" ;
				  }
			      }
			    next ;
			  }
		      }
		    next ;
		  }

		if ($key eq "gene")
		  {
		    #print "// key=$key\n" ;
		    foreach $qual (@quals)
		      {
			#print "// qual $qual\n" ;
			if ($qual =~ /^gene=/)
			  {
			    $qual =~ s/gene=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    # in worm case, check for structure of the locus name
			    # next if if ($active_species =~ /worm/ && !($qual =~ /^[A-Za-z][A-Za-z][A-Za-z]-[0-9]/ ))
			    print "Locus_GB \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^note=/)
			  {
			    $qual =~ s/note=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Remark \"gene note = $qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^db_xref=/)
			  {
			    $qual =~ s/db_xref=// ;
			    $qual =~ s/\n/ / ;
			    $qual =~ s/\"//g ;
			    @notesa = split (/\;\s*/,$qual) ;
			    foreach $note (@notesa)
			      {
				#print "// note=$note\n" ;
				if ($note =~ /^GeneID:/)
				  {
				    $note =~ s/^GeneID:// ;
				    $geneidNam = $note ;
				    print "GeneID \"$geneidNam\"\n" ;
				  }
				if ($note =~ /^MIM:/ &&  $geneidNam)
				  {
				    $note =~ s/^MIM:// ;
				    # print "\nGeneID \"$geneidNam\"\n" ;
				    print "Extern \"OMIM_$note\"\n" ;
				    # print "\nSequence : $nam\n" ;
				  }
			      }
			    next ;
			  }
		      }
		    next ;
		  }

		if ($key eq "CDS")
		  {
                    if ( $loc =~ /([0-9]+)\.\.([0-9]+)/)
                        {
                           print "Genbank_CDS $1 $2\n" ;
                        }

		    if ( $loc =~ /^join/ )
			{

			# The gene data we want is really mrna, which is
			# usually what is there, but sometimes the whole
			# gene is there with the introns/exons marked.
			# If that is the case, we don't want the sequence
			# at all.  It it too late to avoid printing it to
			# the output file, so we set bad_sequence=1 to
			# cause it to be deleted later.

			# printf(STDERR "CDS JOIN PRESENT IN %s\n",$nam);
			# $bad_sequence=1;
			}

		    foreach $qual (@quals)
		      {
			# print "// qual $qual\n" ;
			if ($qual =~ /^function=/)
			  {
			    $qual =~ s/function=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Function \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^gene=/)
			  {
			    $qual =~ s/gene=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    if ($qual =~ /^[A-Za-z][A-Za-z][A-Za-z]-[0-9]/ )
				{
			    	print "Locus_GB \"$qual\"\n" ;
				}
			    next ;
			  }
			if ($qual =~ /^note=/)
			  {
			    $qual =~ s/note=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Remark \"CDS note = $qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^product=/)
			  {
			    $qual =~ s/product=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Remark \"CDS product = $qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^db_xref=/)
			  {
			    $qual =~ s/db_xref=// ;
			    $qual =~ s/\n/ / ;
			    $qual =~ s/\"//g ;
			    @notesa = split (/\;\s*/,$qual) ;
			    foreach $note (@notesa)
			      {
				#print "// note=$note\n" ;
				if ($note =~ /^GeneID:/)
				  {
				    $note =~ s/^GeneID:// ;
				    $geneidNam = $note ;
				    print "GeneID \"$geneidNam\"\n" ;
				  }
				if ($note =~ /^MIM:/ &&  $geneidNam)
				  {
				    $note =~ s/^MIM:// ;
				    # print "\nGeneID \"$geneidNam\"\n" ;
				    print "Extern \"OMIM_$note\"\n" ;
				    # print "\nSequence : $nam\n" ;
				  }

			      }
			    next ;
			  }
			if ($qual =~ /^clone_lib=/)
			  {
			    $qual =~ s/clone_lib=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Library \"$qual\"\n" ;
			    next ;
			  }
			if ($qual =~ /^cell_type=/)
			  {
			    $qual =~ s/cell_type=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Celltype \"$qual\"\n" ;
			    next ;
			  }
		      }
		    next ;
		  }

		if ($key eq "PolyA_site")
		  {
		    foreach $qual (@quals)
		      {
			#print "// qual $qual\n" ;
			if ($qual =~ /^function=/)
			  {
			    $qual =~ s/function=// ;
			    $qual =~ s/"//g ;
			    $qual =~ s/\n/ /g ;
			    print "Function \"$qual\"\n" ;
			    next ;
			  }
		      }
		   }

		if ($key eq "mRNA" && $accept_join != 1)
		  {
		    if ( $loc =~ /^join/ )
			{

			# The gene data we want is really mrna, which is
			# usually what is there, but sometimes the whole
			# gene is there with the introns/exons marked.
			# If that is the case, we don't want the sequence
			# at all.  It it too late to avoid printing it to
			# the output file, so we set bad_sequence=1 to
			# cause it to be deleted later.

			printf(STDERR "mRNA JOIN PRESENT IN %s\n",$nam);
			$bad_sequence=1;
			}
		  next;
		  }

		if ($key eq "PolyA_signal")
		  {
		  next;
		  }

		if ($key eq "exon")
		  {
		  next;
		  }

		if ($key eq "misc_feature")
		  {
		  next;
		  }

		if ($key eq "intron")
		  {
		  next;

		  }
		if ($key eq "3'UTR")
		  {
		  next;
		  }

		if ($key eq "5'UTR")
		  {
		  next;
		  }

		if ($key eq "-")
		  {
		  next;
		  }

	    }
	next;
	}
    }

    print "\n";

    if ($bad_sequence) {

	# it was a bad sequence that was identified as bad after we already
	# started printing it to the file.  follow it with a delete command.

	print "\n-D Sequence $nam\n\n";

    } else {

	# it was a good sequence, so output the rest of the related records

        print "$papers";

        print  "\nDNA $nam\n" ;
        print  $seq ;
    
        print "\n";

        if ($LongText eq 'YES'){
            $text =~ s|\n//||;
	    print  "\nLongText $nam\n" ;
	    print  $text ;
	    print  "***LongTextEnd***\n" ;
        }
    }
}

exit($exitcode);



sub Paper {

  my($seqNam,$paper) = @_;
  my($medId,$pmId,$authors,$title,$jnl,@authors,@jnl,$jname) = '';
  my($issue,$page1,$pageEnd,$year,$zoo,$foo,$goo,$part) = '';
  my($bin,$authorFirst,$paperFirst) = '';

  $paper =~ s/\n            / /g;

  for(split(/\n/,$paper)) {
	$authors = $1 if /  AUTHORS   (.*)/;
	$title = $1 if /  TITLE     (.*)/;
	$medId = $1 if /  MEDLINE   (\d+)/;
	$pmId = $1 if /  PUBMED   (\d+)/;
	$jnl = $1 if /  JOURNAL   (.*)/;
  }

  # process title
  if ($title =~ /Direct Submission/)
	{ $title = $jnl ; $jnl = '' ; }

  # process authors
  $authors =~ s/ and /, /;
  $authors =~ s/,(\S)/ $1/g;
  @authors = split(/, /,$authors);

  # process journal, date and pages
  $jnl =~ s/ (\d)/ZZZZ$1/;
  ($jname,$zoo) = split(/ZZZZ/,$jnl);
  $zoo =~ s/, /ZZZZ/;
  ($foo,$goo) = split(/ZZZZ/,$zoo);
  $issue = $1 if $foo =~ /^(\d+)/;
  $part = $1 if $foo =~ /\((\d+)\)/;
  if ($goo =~ /(\d+)?-?(\d+)?/) { $page1 = $1; $pageEnd = $2 };
  $year = $1 if $goo =~ /\((\d+)\)/;
  $jnl = '' if ($jnl =~ /Unpublished/) ;

  # $nPaper is a global, yuuk
  $nPaper++ ;
  if ($pmId ne '') { $paperName = "pm$pmId"; }
  else { $paperName = $seqNam.'_'.$nPaper; }
# print "XXXXXXXXXXXXX $nPaper  $paperName\n";

  print "Reference $paperName\n";

  # $papers is a global
  $papers .= "\nPaper $paperName\n"; # ! so i will rather import the paper from medline directly
  if  ($pmId ne '')
    {  $papers .= "PMID $pmId\n" ; }
  else # give details
    {
      $papers .= "Title \"$title\"\n" unless $title eq '' ;
      foreach(@authors)
	{ $papers .= "Author \"$_\"\n";  }
      $papers .= "Journal \"$jname\"\n" unless $jname eq '' ;
      $papers .= "Year $year\n" unless $year eq '' ;
      #$papers .= "Volume $issue\n";
      $papers .= "Page $page1 $pageEnd\n" unless $page eq '' ;
      $papers .= "Medline_acc $medId\n\n" unless $medId eq '';
    }

}



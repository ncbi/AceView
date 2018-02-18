#!/bin/env perl -w

# @(#)embl2ace.pl	1.1 5/24/97
# SCCSID = @(#) Module embl2ace.perl File SCCS/s.embl2ace.perl ()
#          @(#) Version 1.32 created on 05/30/96 at 10:36:47

#========================================================================#
# Program name: embl2ace.perl                                            #
# Author:       Zhuoan Jiao                                              #
# Date:         20 December 1994 - May 1996                              #
# Note:         This program converts flat EMBL sequence files into ace  #
#		format.							 #
#========================================================================#

#------------------------------------------------------------------------
# Use Entrez data provided by Dr. Jaime Prilusky at Weizmann to
# locate the medline unique identifier of a citation. Program still works
# when no Entrez data available.
#------------------------------------------------------------------------

$PATH_Entrez_by_authors = "$ENV{'HOME'}/igd_paths/Entrez/by-authors";
$PATH_Entrez_by_journal = "$ENV{'HOME'}/igd_paths/Entrez/by-journal";

$| = 1; #turns off output buffering for STDOUT.

# initialise object ID
undef $VAR_object_id;

undef %VAR_citation_IDs;
undef %VAR_citations_of_entry;
undef %VAR_citation_Authors;
undef %VAR_citation_Title;
undef %VAR_citation_Journal;
undef %VAR_citation_Volume;
undef %VAR_citation_Issue;
undef %VAR_citation_Pages;
undef %VAR_citation_Year;
undef %VAR_citation_Summary;
undef %VAR_citation_Type;

$VAR_subseq{'CDS'} = 1;
$VAR_subseq{'mRNA'} = 1;
$VAR_subseq{'tRNA'} = 1;
$VAR_subseq{'rRNA'} = 1;
$VAR_subseq{'snRNA'} = 1;

$VAR_DISPLAY_feature{'TSL_site'} = 1;
$VAR_DISPLAY_feature{'TSL'} = 1;
$VAR_DISPLAY_feature{'Possible_exon'} = 1;
$VAR_DISPLAY_feature{'promoter'} = 1;
$VAR_DISPLAY_feature{'polyA_site'} = 1;
$VAR_DISPLAY_feature{'polyA_signal'} = 1;
$VAR_DISPLAY_feature{'misc_signal'} = 1;
$VAR_DISPLAY_feature{'misc_feature'} = 1;
$VAR_DISPLAY_feature{'repeat_region'} = 1;
$VAR_DISPLAY_feature{'repeat_unit'} = 1;
$VAR_DISPLAY_feature{'mutation'} = 1;
$VAR_DISPLAY_feature{'sig_peptide'} = 1;
$VAR_DISPLAY_feature{'mat_peptide'} = 1;
$VAR_DISPLAY_feature{'old_sequence'} = 1;
$VAR_DISPLAY_feature{'modified_base'} = 1;
$VAR_DISPLAY_feature{'TATA_signal'} = 1;

$VAR_sequence_reference{'SWISS-PROT'} = 'sw';


($#ARGV>=0) || die
   "usage: embl2ace [-D] [-o <output_filename>] <filename>
    if -D is presented, the ace file created will have \"-D <Class>: id\" line
    proceeded to each object \"<Class>: id\", as the \"-D\" flag of ACeDB;
    use - for standard input;
    the default output is on the screen;
    options must be given in the order shown above\n";

if ($ARGV[0] =~ /^-D$/) 
{	$VAR_D = 1;
	shift 
}
else {  $VAR_D = 0 };

if ($ARGV[0] =~ /^-o$/)
{       shift;
        # set the filehandle for output as $ARGV[0]
        open(OUTPUT, ">$ARGV[0]");     
        select(OUTPUT);       
        shift
};

SCAN:
while (<>)
{    
	next SCAN if /^\s*$/;			# skip blank lines.
	next SCAN if /^XX/;			# skip XX lines
        next SCAN if /^KW\s+\.\s*$/;            # skip empty KW line
        next SCAN if /^OC\s+\.\s*$/;            # skip empty OC line

        if (/^\/\//)				# // lines
	{		
	    # output (last) entry details.
	    &output_entry_details if ($VAR_object_id);

	    next SCAN
	};

	if (/^ID/) 
	{ 
		# start of an entry, initialize variables.
		&variables_initialisation;

		# split;	# does not work in Perl5
		@ID_elements = split(' ', $_, 3);

		$VAR_entry_id = $ID_elements[1];
		$VAR_cDNA = "cDNA" if ($_ =~ /RNA;/);  

		# $1 refers to the part matching [0-9]*
		$VAR_length = $1 if ($_ =~ /([0-9]+) BP/); 

		next SCAN
	};

	if (/^AC\s+(.+)/) 
	{       
		$VAR_accessions = $1;
 
            # CONTINUE_AC:
                while (<>)
                {       if (/^AC\s+(.+)/) {
                           $VAR_accessions .= $1 }
                        else { last };
                };
 
                @VAR_accession = split(/;\s*/, $VAR_accessions, 2);
                $VAR_primary_accession = $VAR_accession[0];
                $VAR_object_id = "em:$VAR_primary_accession";
		$VAR_accessions =~ s/;/ /g;    # to be the same as GenBank data.
		$VAR_accessions =~ s/\s{2,}/ /g;
		$VAR_accessions =~ s/\s$//g;    # no ending space

		redo SCAN;
	};

	if (/^DT\s+(.+)/)
	{	
		if ($VAR_date) {
			$VAR_date .= "; $1" }
		else {	$VAR_date = $1 };
 
                next SCAN
	};

	if (/^DE\s+(.+)/)
	{
		if ($VAR_desciption) {
			if ($1 =~ /\-$/) { # last character is "-" 
				$VAR_desciption .= $1}
			else {  $VAR_desciption .= " $1"} }
		else { 	$VAR_desciption = $1 };

		next SCAN
	};

	if (/^KW\s+(.+)/)
	{	
		if ($VAR_keywords) {
			$VAR_keywords .= " $1" }
		else {  $VAR_keywords  = $1};

		$VAR_keywords =~ s/\.\s*$//g;      # remove trailing "." if any

		next SCAN
	};
		
        if (/^OS\s+(.+)/)
	{
		next SCAN if ($1 =~ /Unknown/i || /None/i);

		if ($VAR_species) { # OS may scatter around.
		      	$VAR_species .= ";;;$1" }
		else {  $VAR_species = $1 };

	   	next SCAN
	};

	if (/^OC\s+(.+)/)
	{
		$VAR_classifications .= $1;

		$VAR_classifications =~ s/\.\s*$/;/;   # use ";" is useful!

		next SCAN

	};

        if (/^DR\s+(.+)/)
	{  
		&store_cross_reference($1);  
                next SCAN
	};

	if (/^CC\s+(.+)/)
	{
		if ($VAR_comment) {
			if ($1 =~ /\-$/) { # last character is "-"
				$VAR_comment .= $1 }
			else {  $VAR_comment .= " $1"}}
		else {  $VAR_comment = $1};

		next SCAN
	};

        if (/^FT\s{3}\S+/ )
        {   
		# begin of a feature ... variable initialisation
            	$VAR_FT_qualifier = "";
 
                chop;
                $VAR_FT_key = substr($_, 5, 15);
                $VAR_FT_key =~ s/\s+$//g;         # remove trailing spaces
 
                $VAR_FT_location = substr($_, 21);
 
            # continue FT location ...
                while (<>)
                {
                    if (/^FT\s{19}\S+/ && substr($_,21,1) ne '/')
                    {	# not begins with a "/", location line continues ...
                        chop;
                        $VAR_FT_location .= substr($_, 21);
                    }
		    else { last }
		};

                # tidy up location line ... create value for $VAR_FT_id
                $VAR_FT_location =~ s/\s{2,}/ /g;        # single space
                $VAR_FT_location =~ s/\s+$//g;           # remove end space
                $VAR_FT_location =~ s/\"/'/g;
                $VAR_FT_location =~ s/\s*,\s*/,/g;       # no space around ","
                $VAR_FT_location =~ s/\s*\.\s*\.\s*/\.\./g; # remove spaces around ".."
 
                # store information ...
                $VAR_FT_count++;
                $VAR_key_count{$VAR_FT_key}++;
                $VAR_FT_id = $VAR_FT_key.'.'.$VAR_key_count{$VAR_FT_key};
                $VAR_FT_list[$VAR_FT_count] = $VAR_FT_id;
                $VAR_FT_key{$VAR_FT_id} = $VAR_FT_key;
                $VAR_FT_location{$VAR_FT_id} = $VAR_FT_location;
 
             FIND_qualifiers: 
             { 
		if (/^FT\s{19}\//)
		{	# qualification line begins ...
			# initialise variable.
                        $VAR_FT_qualifier = "";
                        $VAR_AA_translation = "";
		     #  undef $VAR_translation_id;  (future version needs this)

                        chop;
			$VAR_FT_qualifier = substr($_, 21);
			$VAR_AA_translation = 1 if ($VAR_FT_qualifier =~ /translation=/);
 
			# continue FT qualifiers ...
                        while (<>)
                        {
			    if (/^FT\s{19}\S+/ && substr($_,21,1) ne '/')
			    {	chop;

                                if ($VAR_AA_translation) {
                                        $VAR_FT_qualifier .= substr($_, 21) }
                                else { 
                                        $VAR_FT_qualifier .= ' '.substr($_, 21);
                                        $VAR_FT_qualifier =~ s/\s{2,}/ /g }; 
                                next
                            };
 
                            last
                        };
 
                        # end of current qualifier line,
                        # tidy up qualifiers
                        $VAR_FT_qualifier =~ s/\"/\'/g;
 
			if ($VAR_FT_qualifier{$VAR_FT_id}) {
				$VAR_FT_qualifier{$VAR_FT_id} .= " $VAR_FT_qualifier"}
			else {	$VAR_FT_qualifier{$VAR_FT_id} = $VAR_FT_qualifier }
 
                        redo FIND_qualifiers      # find next qualifier(s) of the current FT key

		}; # end of current FT key-location-qualifier pair ...
	     };   # end of FIND_qualifiers block
 
	     redo SCAN  # process next feature if any
        };
			
	if (/^SQ\s+Sequence\s+(.+)/)
	{
		$VAR_basepairs = $1;

		# these replacements are optional.
		$VAR_basepairs =~ s/;$/\./;	# replace last ";" with "."
		$VAR_basepairs =~ s/BP/base pairs/;

	   # GET_sequence:
		while (<>)
		{
			last if (/^\/\//);

			$_ =~ s/[0-9]+//g;
			$_ =~ s/\s+//g;

			$VAR_sequence .= $_;
			$VAR_sequence .= "\n";
		};

		redo SCAN
	};

        if (/^RN\s+/)
        { 
                $VAR_citation_id = "";
                $VAR_first_author = "";;
                $VAR_citation_authors = "";
                $VAR_citation_title = "";
                $VAR_location_line = "";
                $VAR_citation_db = "";
                $VAR_citation_db_id = "";
		$VAR_medline_id = "";

	   REFERENCE_BLOCK:
		while(<>)
		{
			next REFERENCE_BLOCK if (/^XX\s+/);
			next REFERENCE_BLOCK if (/^RP\s+/);
			next REFERENCE_BLOCK if (/^RC\s+/);
			next REFERENCE_BLOCK if (/^RA\s+;$/);  # no 'author'
			next REFERENCE_BLOCK if (/^RT\s+;$/);  # no 'title'

			if (/^RX\s+(.+);\s*(.+)/)             # new feature
			{   
			   	$VAR_citation_db = $1;
				$VAR_citation_db_id = $2;
				$VAR_citation_db_id =~ s/\.$//g;

				if ($VAR_citation_db =~ /Medline/i) {
					$VAR_medline_id = "med:$VAR_citation_db_id";
					$VAR_citation_id = $VAR_medline_id;

				#	$VAR_citation_IDs{$VAR_citation_id}++;

					if ($VAR_citations_of_entry{$VAR_object_id}) {
						$VAR_citations_of_entry{$VAR_object_id} .= ";;; $VAR_citation_id" }
					else {  $VAR_citations_of_entry{$VAR_object_id} = $VAR_citation_id};

					last REFERENCE_BLOCK
				}
				else {	next REFERENCE_BLOCK }
			};

			if (/^RA\s+(.+)/)
			{  	
				# authors
				$VAR_citation_authors = $1;
	
		  	    # continue RA line ...
				while(<>)
				{	if (/^RA\s+(.+)/) {
					       $VAR_citation_authors .= " $1" }
					else { last }
				};

			    #	$VAR_citation_authors =~ s/\./ /g;      # replace . with a space.
				$VAR_citation_authors =~ s/\.//g;	# as in Entrez database.
				$VAR_citation_authors =~ s/\s+,/,/g;    # remove space(s) before a ','
				$VAR_citation_authors =~ s/\s{2,}/ /g;  # single space
				$VAR_citation_authors =~ s/[\s;]*$//g;	# remove last ';' and ending spaces.

				redo REFERENCE_BLOCK
			};

			if (/^RT\s+(.+)/)
			{	
				next REFERENCE_BLOCK if ($1 =~ /^;\s*$/);

				# Title
                        	$VAR_citation_title = $1;

   
                   	   # continue RT line ...
                        	while(<>)
                        	{   
				    if (/^RT\s+(.+)/) {
                                           $VAR_citation_title .= " $1" }
                                    else { last }
                        	};

				$VAR_citation_title =~ s/\"//g;		  # remove "
				$VAR_citation_title =~ s/\s{2,}/ /g;      # single space
				$VAR_citation_title =~ s/;\s*$//g;        # remove last ";".
				$VAR_citation_title =~ s/^The //;	  # remove "The" 
									  # ("WashU-Merck EST Project" = "The WashU-Merck EST Project")
                
				redo REFERENCE_BLOCK
			};

			if (/^RL\s+(.+)/)
			{
                	        $VAR_location_line = $1;

			   # continue RL line ...
			        while(<>)
        			{
        	        	   if (/^RL\s+(.+)/) {
        	        	            $VAR_location_line .= " $1" }
                		    else {  last }
        			};

				$VAR_location_line =~ s/\"//g;
				$VAR_location_line =~ s/\./ /g;		# replace "." with " ".
				$VAR_location_line =~ s/\s{2,}/ /g;
				$VAR_location_line =~ s/\s,/,/g;
			
                        	redo REFERENCE_BLOCK
                	};

	                # call subroutine "get_citation_details" to get values for:
           	        # $VAR_citation_id, and
	                # %VAR_citation_{Authors,Title,Journal,Volume,Pages,Year,Summary,Type}
			# and change the "$VAR_use_medline_id' value when necessary.
			&get_citation_details($VAR_citation_authors, $VAR_citation_title, $VAR_location_line);

			if ($VAR_citation_id)
			{   
			    # When "$VAR_medline_id" is not empty, EntrezToAcedb will
			    # produce details for object "$VAR_citation_id".
			    $VAR_citation_IDs{$VAR_citation_id}++ if ($VAR_medline_id eq "");

  			    if ($VAR_citations_of_entry{$VAR_object_id}) {
				    $VAR_citations_of_entry{$VAR_object_id} .= ";;; $VAR_citation_id" }
			    else {  $VAR_citations_of_entry{$VAR_object_id} = $VAR_citation_id}
			};

			last REFERENCE_BLOCK

		};   # end of REFERENCE_BLOCK

		redo SCAN 

	};  # end of RN blocks.
 };


#-------------------------------------------------------------------
# SUBROUTINE: VARIABLES_INITIALISATION
#
# Before processing an EMBL entry, initialise variables.
#-------------------------------------------------------------------
sub variables_initialisation
{
        $VAR_object_id = "";
        $VAR_entry_id = "";
        $VAR_primary_accession = "";    # used to construct $VAR_object_id.
        $VAR_accessions = "";
        $VAR_keywords = "";
        $VAR_date = "";
        $VAR_cDNA = "";
        $VAR_length = "";
        $VAR_desciption = "";
        $VAR_species = "";
        $VAR_classifications = "";
        $VAR_comment = "";
        $VAR_cross_references = "";
        $VAR_sequence = "";
        $VAR_basepairs = "";
 
        $VAR_FT_count = 0;
        undef %VAR_key_count;
        undef %VAR_FT_list;
        undef %VAR_FT_key;
        undef %VAR_FT_location;
        undef %VAR_FT_qualifier;
 
        $VAR_subseq_count = 0;
        $VAR_subseq_offset = 0;
        undef %VAR_subseq_list;
        undef %VAR_subseq_source;
        undef %VAR_subseq_location;
        undef %VAR_subseq_offset;
        undef %VAR_subseq_comment;
 
        undef %VAR_corresponding_peptide;
        undef %VAR_foreign_reference;
};

#-------------------------------------------------------------------
# SUBROUTINE: OUTPUT_ENTRY_DETAILS
# Note: Output detail information on current sequence entry.
#       Variables referred are global variables.
#-------------------------------------------------------------------
sub output_entry_details {

	print "\n-D Sequence : \"$VAR_object_id\"\n" if ($VAR_D);

	print "\nSequence : \"$VAR_object_id\"\n";

	print "generic_properties identification local_id \"$VAR_entry_id\"\n";

	if ($VAR_keywords)
	{   $VAR_keywords =~ s/\"/'/g;
            @VAR_keyword = split(/;\s*/, $VAR_keywords);
            for ($i = 0; $i <= $#VAR_keyword; $i++)
            { print "generic_properties keyword \"$VAR_keyword[$i]\"\n" };
        };

	if ($VAR_date)
	{   @VAR_date = split(/;\s*/, $VAR_date);
            for ($i = 0; $i <= $#VAR_date; $i++)
            {  	if ($VAR_date[$i] =~ /([0-9A-Z\-]+).+Created/i) {
	       		print "generic_properties add_date \"$1\"\n" }
	     	elsif ($VAR_date[$i] =~ /([0-9A-Z\-]+).+updated/i) {
			print "generic_properties mod_date \"$1\"\n" }
		else {  print "Date \"$VAR_date[$i]\"\n"}
	    }
	};

	if ($VAR_desciption)
	{
	    $VAR_desciption =~ s/\"/'/g;
	    $VAR_desciption =~ s/\s{2,}/ /g;   # only single space.
	    print "Title \"$VAR_desciption\"\n"
	};

        if ($VAR_length) {
            print "DNA \"$VAR_object_id\" \"$VAR_length\"\n";
            print "Length \"$VAR_length\"\n" }
	else {  print "DNA \"$VAR_object_id\"\n" };

	#print "Data_Source \"EMBL\" \"$VAR_entry_id\" \"$VAR_primary_accession\"\n";  # use 'Data_Source' to replace 'Library'!
	print "Library \"EMBL\" \"$VAR_entry_id\" \"$VAR_primary_accession\"\n";

	if ($VAR_species) 
	{
		$VAR_species =~ s/\"/'/g;
		$VAR_species =~ s/\s{2,}/ /g;

		@VAR_species = split(/;;;/, $VAR_species);
		for ($i = 0; $i <= $#VAR_species; $i++)
		{ 
			if ($VAR_species[$i] =~ /([^\(\)]+)\((.+)\)/) {
				# e.g. "OS   EQUUS ASINUS (DONKEY) (ANACYSTIS NIDULANS R2)."
				$VAR_species_id = $1;
				$VAR_species_id =~ s/\s+$//g;

			#	print "Species \"$VAR_species_id\" -C \"\($2\)\"\n"} 
				print "Species \"$VAR_species_id\"\n"}
			else {  print "Species \"$VAR_species[$i]\"\n" }
		}
	};

  #             for ($i = 0; $i <= $#VAR_accession; $i++)
  #              {  print "Accessions \"$VAR_accession[$i]\"\n"};
	print "Accessions \"$VAR_accessions\"\n";

	print "$VAR_cDNA\n" if ($VAR_cDNA);		 

	if ($VAR_classifications)
        { 
	    $VAR_classifications =~ s/\"/'/g;
            @VAR_classification =  split(/;\s*/, $VAR_classifications);
            for ($i = 0; $i <= $#VAR_classification; $i++)
            {  print "Classification \"$VAR_classification[$i]\"\n" }
	};
		
	if ($VAR_cross_references)
	{	$VAR_cross_references =~ s/\"/'/g;
		@VAR_cross_reference = split(/;\s*/, $VAR_cross_references);
		for ($i = 0; $i <= $#VAR_cross_reference; $i++) 
		{	print "Cross_reference \"$VAR_cross_reference[$i]\"\n" }
	};	

	if ($VAR_comment)
	{	$VAR_comment =~ s/\"/'/g;       # replace " with '
               	$VAR_comment =~ s/\s{2,}/ /g;  # single space.
               	print "generic_properties remark \"$VAR_comment\"\n"
	};

	if ( $VAR_citations_of_entry{$VAR_object_id})
	{ 	@VAR_citations_in_entry = split(/;;;\s*/, $VAR_citations_of_entry{$VAR_object_id});
		for ($i = 0; $i <= $#VAR_citations_in_entry; $i++)
		{ print "Reference \"$VAR_citations_in_entry[$i]\"\n"}
	};

        if (%VAR_corresponding_peptide)
        {   foreach $VAR_Ref_ID (keys(%VAR_corresponding_peptide)) {
                 print "Corresponding_peptide \"$VAR_Ref_ID\"\n"}
        };
 
        if (%VAR_foreign_reference)
        {   foreach $VAR_Ref_ID (keys(%VAR_foreign_reference)) {
                 print "Foreign_Reference \"$VAR_Ref_ID\"\n"}
        };

	print "BasePairs \"$VAR_basepairs\"\n" if $VAR_basepairs;

	&display_feature_table;

	print "\n";
	&display_subsequence;

        print "\n";
        print "DNA : \"$VAR_object_id\"\n";
	print "$VAR_sequence"
};


#-------------------------------------------------------------------
# Output Citation objects. This is done after all sequence entries
# have been output.
#-------------------------------------------------------------------

print "\n";

foreach $VAR_citation_id (keys(%VAR_citation_IDs))
{
 #   	print "\n-D Citation : \"$VAR_citation_id\"\n" if $VAR_D;  # <== Wrong!
	print "\nCitation : \"$VAR_citation_id\"\n";
	print "type \"$VAR_citation_Type{$VAR_citation_id}\"\n" if $VAR_citation_Type{$VAR_citation_id};
	print "title \"$VAR_citation_Title{$VAR_citation_id}\"\n" if $VAR_citation_Title{$VAR_citation_id};

	if (@VAR_authors) {
		@VAR_authors = split(/,\s*/, $VAR_citation_Authors{$VAR_citation_id});
		for ($i = 0; $i <= $#VAR_authors; $i++)
		{ print "Author \"$VAR_authors[$i]\"\n"};
	};

	print "Journal \"$VAR_citation_Journal{$VAR_citation_id}\"\n" if $VAR_citation_Journal{$VAR_citation_id};
	print "volume \"$VAR_citation_Volume{$VAR_citation_id}\"\n" if $VAR_citation_Volume{$VAR_citation_id};
	print "issue \"$VAR_citation_Issue{$VAR_citation_id}\"\n" if $VAR_citation_Issue{$VAR_citation_id};
 	print "pages \"$VAR_citation_Pages{$VAR_citation_id}\"\n" if $VAR_citation_Pages{$VAR_citation_id};
	print "year \"$VAR_citation_Year{$VAR_citation_id}\"\n" if $VAR_citation_Year{$VAR_citation_id};
	print "generic_properties summary \"$VAR_citation_Summary{$VAR_citation_id}\"\n" if $VAR_citation_Summary{$VAR_citation_id};
}


#-------------------------------------------------------------------------            
# Subroutine:   GET_CITATION_DETAILS
#-------------------------------------------------------------------------
sub get_citation_details {

	local($VAR_citation_authors, $VAR_citation_title, $VAR_location_line) = @_;

	local(@VAR_authors, $VAR_book_details, $VAR_length, $i);
	local($VAR_citation_journal, $VAR_citation_volume, $VAR_citation_page, $VAR_citation_year);

   #---------------------------
   # Find the 1st author ...
   #---------------------------
	if ($VAR_citation_authors) { # if $VAR_citation_authors is not empty.
		@VAR_authors = split(/,\s*/, $VAR_citation_authors);
		$VAR_first_author = $VAR_authors[0]
	};

   #----------------------------
   # process $VAR_location_line 
   #----------------------------

	#---------------------------
	# process "Patent number"
	#---------------------------

	if ($VAR_location_line =~ /^Patent number ([a-zA-Z0-9\-\/]+)/i)
	{	
	   if ($VAR_first_author) {
		  $VAR_citation_id = "$VAR_first_author, patent:$1";
		  $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors}
	   else { $VAR_citation_id = "patent:$1" };

	   $VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title if $VAR_citation_title;
	   $VAR_citation_Summary{$VAR_citation_id} = $VAR_location_line;

	   return
	};

        #---------------------------
        # process "Submitted" data.
        #---------------------------   
	if ($VAR_location_line =~ /^Submitted (\([0-9a-z\-]+\))/i)
	{
	   if ($VAR_first_author) {
                  $VAR_citation_id = "$VAR_first_author, submitted:$1";
		  $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors}
           else { $VAR_citation_id = "submitted:$1" };
                
	   $VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title if $VAR_citation_title;
           $VAR_citation_Summary{$VAR_citation_id} = $VAR_location_line;

           return
        };

        #---------------------------
        # process "Thesis"
        #---------------------------   
	if ($VAR_location_line =~ /^Thesis (\([0-9a-z\-]+\))/i)
	{
           if ($VAR_first_author) {
                  $VAR_citation_id = "$VAR_first_author, thesis:$1";
	   	  $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors}
           else { $VAR_citation_id = "thesis:$1" };

	   $VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title if $VAR_citation_title;
           $VAR_citation_Summary{$VAR_citation_id} = $VAR_location_line;
	   $VAR_citation_Type{$VAR_citation_id} = 'Thesis';        

           return
        };

        #---------------------------
        # process books ("in")
        #---------------------------   
	if ($VAR_location_line =~ /^\(in\)(.+)/i) 
        {
	   $VAR_book_details = $1;
	   $VAR_book_details =~ s/\s{2,}/ /g;

	   # take the first 3 words and the last 2 words from $VAR_book_details to 
	   # construct the citation ID.
	   # @VAR_book_details_partial = split(/\s+/, $VAR_book_details);
	   @VAR_book_details_partial = split(' ', $VAR_book_details);

	   $VAR_length  = @VAR_book_details_partial;

	   if ($VAR_length <= 5)  # used to be <=3
	   {	$VAR_book_details_parts = $VAR_book_details }
	   else 
	   {
	   	$VAR_book_details_parts_1 = $VAR_book_details_partial[0];
	   	for ($i = 1; $i <= 2; $i++) { $VAR_book_details_parts_1 .= " $VAR_book_details_partial[$i]"};
	   	$VAR_book_details_parts = "$VAR_book_details_parts_1 ... $VAR_book_details_partial[$VAR_length-2] $VAR_book_details_partial[$VAR_length-1]";
		$VAR_book_details_parts =~ s/\s+$//g
	   };
	   
           if ($VAR_first_author) {
                  $VAR_citation_id = "$VAR_first_author, in: '$VAR_book_details_parts'";
	   	  $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors }
           else {  if ($VAR_citation_title) {
		 	 $VAR_citation_id = "$VAR_citation_title, in: '$VAR_book_details_parts'" }
		   else{ $VAR_citation_id = "in: '$VAR_book_details_parts'"}
		};

	   $VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title if $VAR_citation_title;           
           $VAR_citation_Type{$VAR_citation_id} = 'Book_section';
	   $VAR_citation_Summary{$VAR_citation_id} = $VAR_book_details;		

           return
        };
 
	#---------------------------
        # process "Unpublished"
        #---------------------------   
	if ($VAR_location_line =~ /^Unpublished/i && $VAR_citation_title) 
        {  
	   # take the first 3 words and the last 2 words from $VAR_book_details to
	   # construct the citation ID.
	   @VAR_title_partial = split(/\s+/, $VAR_citation_title);
	   $VAR_length  = @VAR_title_partial;

	   if ($VAR_length <= 5) {
		   $VAR_title_parts = $VAR_citation_title}
	   else {
	   	   $VAR_title_parts_1 = $VAR_title_partial[0];
	   	   for ($i = 1; $i <= 2; $i++) {
	   		$VAR_title_parts_1 .= " $VAR_title_partial[$i]"};
	   	   $VAR_title_parts = "$VAR_title_parts_1 ... $VAR_title_partial[$VAR_length-2] $VAR_title_partial[$VAR_length-1]";
		   $VAR_title_parts =~ s/\s+$//g};
	
           if ($VAR_first_author) {
                  $VAR_citation_id = "$VAR_first_author, unpublished:'$VAR_title_parts'";
		  $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors }
           else { $VAR_citation_id = "unpublished:'$VAR_title_parts'" };
           
	   $VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title;
           $VAR_citation_Type{$VAR_citation_id} = 'Unpublished';

           return
        };

	if ($VAR_location_line =~ /^Unpublished/i && ! $VAR_citation_title)
        {  
	   # $VAR_citation_title is unknown.
	
           if ($VAR_first_author) {
                  $VAR_citation_id = "$VAR_first_author, unpublished";
		  $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors }
           else { $VAR_citation_id = "unpublished" };
           
           $VAR_citation_Type{$VAR_citation_id} = 'Unpublished';

           return
        };

        #-------------------------------------------------
        # process "journal"  (has 'issue')
	#
	# (1) J. Biol. Chem. 266(c):23053-23059 (1991)
        #-------------------------------------------------
	if ( $VAR_location_line =~ /(.+)\s([0-9]+)\s?\((.+)\)\s?:\s?([0-9\-]+)\s?\(([0-9]+)\)/)
	{ 
	   $VAR_citation_journal = $1;		
	   $VAR_citation_volume  = $2;
	   $VAR_citation_page    = $4;
	   $VAR_citation_year    = $5;
		
	   &find_medline_id($VAR_first_author, $VAR_citation_journal, $VAR_citation_volume, $VAR_citation_page, $VAR_citation_year);

	   if ($VAR_medline_id)
	   {    $VAR_citation_id = $VAR_medline_id;
	   }
	   else
	   {
		if ($VAR_first_author) {
			$VAR_citation_id = "$VAR_first_author, $1, $2:$4 \($5\)";
			$VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors	}
		else {	$VAR_citation_id = "$1, $2:$4 \($5\)"};

       		$VAR_citation_Journal{$VAR_citation_id} = $1;
		$VAR_citation_Volume{$VAR_citation_id}  = $2;
		$VAR_citation_Issue{$VAR_citation_id}   = $3;
		$VAR_citation_Pages{$VAR_citation_id}   = $4;
		$VAR_citation_Year{$VAR_citation_id}    = $5;
		$VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title;
		$VAR_citation_Type{$VAR_citation_id}    = 'Article';
	   };

	   return
	};

        #-------------------------------------------------
        # process "journal"    (no 'issue')
        #
        # (2) J. Biol. Chem. 266:23053-23059(1991)
        #-------------------------------------------------
	if ( $VAR_location_line =~ /(.+)\s+([0-9a-zA-Z]+)\s?:\s?([0-9\-]+)\s?\(([0-9]+)\)/)
	{ 
	   $VAR_citation_journal = $1;		
	   $VAR_citation_volume  = $2;
	   $VAR_citation_page    = $3;
	   $VAR_citation_year    = $4;

	   &find_medline_id($VAR_first_author, $VAR_citation_journal, $VAR_citation_volume, $VAR_citation_page, $VAR_citation_year);

	   if ($VAR_medline_id)
	   {    $VAR_citation_id = $VAR_medline_id;
	   }
	   else
	   {
		if ($VAR_first_author) {
			$VAR_citation_id = "$VAR_first_author, $1, $2:$3 \($4\)";
			$VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors	}
		else {	$VAR_citation_id = "$1, $2:$3 \($4\)"};

       		$VAR_citation_Journal{$VAR_citation_id} = $1;
		$VAR_citation_Volume{$VAR_citation_id}  = $2;
		$VAR_citation_Pages{$VAR_citation_id}   = $3;
		$VAR_citation_Year{$VAR_citation_id}    = $4;
		$VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title;
		$VAR_citation_Type{$VAR_citation_id}    = 'Article';
	   };

	   return
	};

        #-------------------------------------------------
        # process "journal"    (no 'issue' and 'volume')
        #
        # (3) J. CHEM. SOC. CHEM. COMMUN. 679-680(1967)
        #-------------------------------------------------
        if ( $VAR_location_line =~ /(.+)\s+([0-9\-]+)\s?\(([0-9]+)\)/)
        {
           $VAR_citation_journal = $1;
           $VAR_citation_page    = $2;
           $VAR_citation_year    = $3;
           $VAR_citation_volume  = "";
 
           &find_medline_id($VAR_first_author, $VAR_citation_journal, $VAR_citation_volume, $VAR_citation_page, $VAR_citation_year);
 
           if ($VAR_medline_id)
           {    $VAR_citation_id = $VAR_medline_id;
           }
           else
           {
                if ($VAR_first_author) {
                        $VAR_citation_id = "$VAR_first_author, $1, $2 \($3\)";
                        $VAR_citation_Authors{$VAR_citation_id} = $VAR_citation_authors }
                else {  $VAR_citation_id = "$1, $2 \($3\)"};
 
                $VAR_citation_Journal{$VAR_citation_id} = $1;
                $VAR_citation_Pages{$VAR_citation_id}   = $2;
                $VAR_citation_Year{$VAR_citation_id}    = $3;
                $VAR_citation_Title{$VAR_citation_id}   = $VAR_citation_title;
                $VAR_citation_Type{$VAR_citation_id}    = 'Article';
           };
 
           return
        }

};

#-------------------------------------------------------------------------            
# Subroutine:   STORE_CROSS_REFERENCE
#-------------------------------------------------------------------------
sub store_cross_reference {

        local($VAR_input) = @_;
        local(@VAR_cross_ref, $VAR_Ref_ID);
 
        @VAR_cross_ref = split(/\s*;\s*/, $VAR_input);
 
        if ($VAR_sequence_reference{$VAR_cross_ref[0]})
        {
                $VAR_Ref_ID = "$VAR_sequence_reference{$VAR_cross_ref[0]}:$VAR_cross_ref[1]";
 
                $VAR_corresponding_peptide{$VAR_Ref_ID}++;
        }
	else 
	{
	        $VAR_Ref_ID = "$VAR_cross_ref[0]:$VAR_cross_ref[1]";
	        $VAR_foreign_reference{$VAR_Ref_ID}++
	}
};

#--------------------------------------------------------
# Subroutine: DISPLAY_FEATURE_TABLE
#--------------------------------------------------------
sub display_feature_table 
{
 #	for ($i = 1; $i <= $#VAR_FT_list; $i++)  # cause error!
	for ($i = 1; $i <= $VAR_FT_count; $i++)	
	{   
	    $VAR_FT_id = $VAR_FT_list[$i];

	    if ($VAR_FT_location{$VAR_FT_id} =~ /^[\>\<]?([0-9]+)[\>\<]?\.\.[\>\<]?([0-9]+)[\>\<]?$/) 
	    {
		if ($VAR_subseq{$VAR_FT_key{$VAR_FT_id}}) 
		{	
			$VAR_subseq_count++;
			$VAR_subseq_id = "$VAR_object_id:$VAR_FT_id";
			print "Subsequence  \"$VAR_subseq_id\" \"$1\" \"$2\"\n";

			$VAR_subseq_offset  = $1 - 1;
			$VAR_location_list = "$1..$2";  #">" & "<" are removed!

			# Record subsequence details
		        $VAR_subseq_list[$VAR_subseq_count]  = $VAR_subseq_id;
        		$VAR_subseq_source{$VAR_subseq_id}   = $VAR_object_id;
        		$VAR_subseq_key{$VAR_subseq_id}      = $VAR_FT_key{$VAR_FT_id};
			$VAR_subseq_location{$VAR_subseq_id} = $VAR_location_list;
			$VAR_subseq_offset{$VAR_subseq_id}   = $VAR_subseq_offset;
        		$VAR_subseq_comment{$VAR_subseq_id}  = $VAR_FT_qualifier{$VAR_FT_id};
			
			next
		};

		if ($VAR_DISPLAY_feature{$VAR_FT_key{$VAR_FT_list[$i]}})
		{  
		    if ($VAR_FT_qualifier{$VAR_FT_id}) {
			   print "$VAR_FT_key{$VAR_FT_id} \"$1\" \"$2\" \"$VAR_FT_qualifier{$VAR_FT_id}\"\n"}
		    else { print "$VAR_FT_key{$VAR_FT_id} \"$1\" \"$2\"\n"};

		    next
		};

		if ($VAR_FT_qualifier{$VAR_FT_id}) {
		       print "other_features \"$VAR_FT_key{$VAR_FT_id}\" \"$VAR_FT_location{$VAR_FT_id}\" \"$VAR_FT_qualifier{$VAR_FT_id}\"\n" }
		else { print "other_features \"$VAR_FT_key{$VAR_FT_id}\" \"$VAR_FT_location{$VAR_FT_id}\"\n"};

		next
	    };

	    if ($VAR_FT_location{$VAR_FT_id} =~ /^join\(([\>\<]?[0-9]+[\>\<]?\.\.[\>\<]?[0-9]+[\>\<]?,*)+\)$/)
	    {  
                if ($VAR_subseq{$VAR_FT_key{$VAR_FT_id}})
                {
			$VAR_FT_location{$VAR_FT_id} =~ s/[\>\<]+//g;		# remove ">" or "<"

			$VAR_subseq_count++;
                        $VAR_subseq_id = "$VAR_object_id:$VAR_FT_id";

			$VAR_FT_location{$VAR_FT_id} =~ /^join\((.+)\)$/;
			$VAR_location_list = $1;
			@VAR_location_list = split(/,/, $1);

			$VAR_first_location = shift(@VAR_location_list);
			$VAR_last_location  = pop(@VAR_location_list);

			$VAR_first_location =~ /^([0-9]+)/;
			$VAR_subseq_start = $1;
        		$VAR_subseq_offset  = $VAR_subseq_start - 1;

			if ($VAR_last_location) {
				$VAR_last_location  =~ /.+\.\.([0-9]+)$/;
				$VAR_subseq_end = $1;
	                        print "Subsequence  \"$VAR_subseq_id\" \"$VAR_subseq_start\" \"$VAR_subseq_end\"\n"}
		        else {  print "Subsequence  \"$VAR_subseq_id\" \"$VAR_subseq_start\"\n"};

                        # Record subsequence details
                        $VAR_subseq_list[$VAR_subseq_count]  = $VAR_subseq_id;
                        $VAR_subseq_source{$VAR_subseq_id}   = $VAR_object_id;
                        $VAR_subseq_key{$VAR_subseq_id}      = $VAR_FT_key{$VAR_FT_id};
                        $VAR_subseq_location{$VAR_subseq_id} = $VAR_location_list;
                        $VAR_subseq_offset{$VAR_subseq_id}   = $VAR_subseq_offset;
                        $VAR_subseq_comment{$VAR_subseq_id}  = $VAR_FT_qualifier{$VAR_FT_id};

			next
                };

                if ($VAR_DISPLAY_feature{$VAR_FT_key{$VAR_FT_id}})
		{
			$VAR_FT_location{$VAR_FT_id} =~ s/[\<\>]+//g;     # remove ">" & "<".

			$VAR_FT_location{$VAR_FT_id} =~ /^join\((.+)\)$/;
			@VAR_location_list = split(/,/, $1);
			$VAR_list_lenth = @VAR_location_list;

			if ($VAR_FT_qualifier{$VAR_FT_id}) 
			{
			       for ($j = 0; $j <= ($VAR_list_lenth -1); $j++)
			       {
			    	 $VAR_location_list[$j] =~ /([0-9]+)\.\.([0-9]+)/;
                             	 print "$VAR_FT_key{$VAR_FT_id} \"$1\" \"$2\" \"$VAR_FT_qualifier{$VAR_FT_id}\"\n"
			       }
			}
		        else 
			{
			       for ($j = 0; $j <= ($VAR_list_lenth -1); $j++)
                               {
                                 $VAR_location_list[$j] =~ /([0-9]+)\.\.([0-9]+)/;
                                 print "$VAR_FT_key{$VAR_FT_id} \"$1\" \"$2\"\n"
                               }
		        };

			next
		};

		if ($VAR_FT_qualifier{$VAR_FT_id}) {
			print "other_features \"$VAR_FT_key{$VAR_FT_id}\" \"$VAR_FT_location{$VAR_FT_id}\"  \"$VAR_FT_qualifier{$VAR_FT_id}\"\n"}
		else { 	print "other_features \"$VAR_FT_key{$VAR_FT_id}\" \"$VAR_FT_location{$VAR_FT_id}\"\n"};
                
                next
	    };

	    # location is in complex format.
	    if ($VAR_FT_qualifier{$VAR_FT_id}) {
		   print "other_features \"$VAR_FT_key{$VAR_FT_id}\" \"$VAR_FT_location{$VAR_FT_id}\" \"$VAR_FT_qualifier{$VAR_FT_id}\"\n"}
	    else { print "other_features \"$VAR_FT_key{$VAR_FT_id}\" \"$VAR_FT_location{$VAR_FT_id}\"\n"}

	}  # end of "for $i = 0 to ..."
};

#--------------------------------------------------------
# Subroutine: DISPLAY_SUBSEQUENCE
#--------------------------------------------------------
sub display_subsequence
{
	local($j);

	for ($j = 1; $j <= $VAR_subseq_count; $j++)
	{	
		$VAR_subseq_id = $VAR_subseq_list[$j];

	   	print "\n-D Sequence : \"$VAR_subseq_id\"\n" if ($VAR_D); 

		print "\nSequence : \"$VAR_subseq_id\"\n";
		print "Source \"$VAR_subseq_source{$VAR_subseq_id}\"\n";
		&display_exons($VAR_subseq_location{$VAR_subseq_id}, $VAR_subseq_offset{$VAR_subseq_id});
		
		if ($VAR_subseq_key{$VAR_subseq_id} eq 'CDS') 
		{
		    print "CDS\n";
		    if ($VAR_subseq_comment{$VAR_subseq_id}) {
			print "generic_properties description remark \"$VAR_subseq_comment{$VAR_subseq_id}\"\n"};
		    next 
		};

		if ($VAR_subseq_key{$VAR_subseq_id} eq 'mRNA') 
		{
		    print "mRNA\n";
		    if ($VAR_subseq_comment{$VAR_subseq_id}) {
			print "generic_properties description remark \"$VAR_subseq_comment{$VAR_subseq_id}\"\n"};
		    next
		};

		if ($VAR_subseq_comment{$VAR_subseq_id})
		{	print "$VAR_subseq_key{$VAR_subseq_id} \"$VAR_subseq_comment{$VAR_subseq_id}\"\n";
			next
		};

		print "$VAR_subseq_key{$VAR_subseq_id}\n";
	}
}


#--------------------------------------------------------
# Subroutine: DISPLAY_EXONS
#--------------------------------------------------------
sub display_exons
{
	local($VAR_location_list, $VAR_offset) = @_;
	local($VAR_list_length, $VAR_exon_start, $VAR_exon_end, $VAR_exon_length, $i);

	@VAR_location_list = split(/,/, $VAR_location_list);
	$VAR_list_length = @VAR_location_list;

	for ($i = 0; $i <= ($VAR_list_length - 1) ; $i++)
	{ 
		$VAR_location_list[$i] =~ /^([0-9]+)\.\.([0-9]+)$/;
		$VAR_exon_start  = $1 - $VAR_offset;
		$VAR_exon_length = $2 - $1 + 1;
		$VAR_exon_end    = $VAR_exon_start + $VAR_exon_length - 1;
		print "Source_Exons \"$VAR_exon_start\" \"$VAR_exon_end\"\n"
	}
};

#--------------------------------------------------------
# Subroutine: FIND_MEDLINE_ID
#--------------------------------------------------------
sub find_medline_id {

        # Entrez data are not available.
        return if ((! -d "$PATH_Entrez_by_authors") ||
                   (! -d "$PATH_Entrez_by_journal"));
        
	local($VAR_author, $VAR_journal, $VAR_volume, $VAR_page, $VAR_year) = @_;
	local($VAR_matched_line);

	$VAR_journal =~ s#\(#\\\(#g;
	$VAR_journal =~ s#\)#\\\)#g;

	$VAR_year    =~ s#\(#\\\(#g;
	$VAR_year    =~ s#\)#\\\)#g;

	if ($VAR_author)
	{ 
	   $VAR_citation_filename = substr($VAR_author,0,1);
	 # $VAR_citation_filename =~ tr/A-Z/a-z/;   # commented out on 13/07/95.
	   $PATH_Entrez = $PATH_Entrez_by_authors;

           if ($VAR_volume) {
                  $VAR_matched_line = `egrep -i "$VAR_author" "$PATH_Entrez/$VAR_citation_filename" | egrep -i "$VAR_journal" | egrep "$VAR_year" | egrep "$VAR_page" | egrep "$VAR_volume"`}
           else { $VAR_matched_line = `egrep -i "$VAR_author" "$PATH_Entrez/$VAR_citation_filename" | egrep -i "$VAR_journal" | egrep "$VAR_year" | egrep "$VAR_page"`};

	   if ($VAR_matched_line =~ /(med:[0-9]+)/)
	   {
	   	$VAR_medline_id = $1
	   };

	   return
	};

	# If author name is unknown, but journal name is known ...
	if ($VAR_journal)
	{
	   $VAR_citation_filename = substr($VAR_journal,0,1);
	   $PATH_Entrez = $PATH_Entrez_by_journal;;

           if ($VAR_volume) {
                  $VAR_matched_line = `egrep -i "$VAR_journal" "$PATH_Entrez/$VAR_citation_filename" | egrep "$VAR_year" | egrep "$VAR_page" | egrep "$VAR_volume"`}
           else { $VAR_matched_line = `egrep -i "$VAR_journal" "$PATH_Entrez/$VAR_citation_filename" | egrep "$VAR_year" | egrep "$VAR_page"`}

	   if ($VAR_matched_line =~ /(med:[0-9]+)/)
	   { 
	   	$VAR_medline_id       = $1
	   };

	   return
	};
};



	

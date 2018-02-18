
{ gsub (/\"/,"\\\"",$0); }
{ gsub (/http:\//,"http:\\/",$0) ;}

{
    gsub ("&gt;"," __greater_than__ ", $0) ;
    gsub ("&lt;"," __lower_than__ ", $0) ;
    gsub ("&amp;","and", $0) ;
}

{if (biosample == "") next;}

{  # grep the title
    
    i = index ($0, "<h2 class=\\\"title\\\">") ;
    if (i>0) 
    { 
	z = substr ($0, i+20) ;
	j = index (z, "</h2>") ;  # we are garanteed to discover it
	z = substr (z, 1,j-1) ;
	gsub (/Sample from /,"", z) ;
	gsub (/Generic sample from /,"", z) ;
	    if (1)
	    {
		if (z == "na") z = "" ;
		if (z == "NA") z = "" ;
		if (z == "N/A") z = "" ;
		if (z == "n/a") z = "" ;
		if (z == "missing") z = "" ;
		if (z == "not applicable") z = "" ;
		if (z == "non provided") z = "" ;
		if (z == "no description") z = "" ;
		if (z == "unspecified") z = "" ;
		gsub (" __greater_than__ ", ">", z) ;
		gsub (" __lower_than__ ", "<", z) ;
	    }
	    
	if (length (z) > 1)
	    title = "Title \"" z "\"" ;
    }
}
{ # Identifier
    i = index ($0, ">Identifiers</dt><dd>") ;
    if (i > 0)
    {
	z = substr ($0, i+21) ;
	i = index (z, "Sample name") ;
	if (i > 0)
	{
	    z = substr (z, i+13) ;
	    n = split (z, aa, ";") ;
            z = aa[1] ;
	    while (substr(z,1,1) == "<")
	    {
		i = index (z, ">")  ;
		if (i > 0)
		    z = substr (z, i+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0)
		z = substr (z, 1, j-1) ;

	    if (1)
	    {
		if (z == "na") z = "" ;
		if (z == "NA") z = "" ;
		if (z == "N/A") z = "" ;
		if (z == "n/a") z = "" ;
		if (z == "missing") z = "" ;
		if (z == "not applicable") z = "" ;
		if (z == "non provided") z = "" ;
		if (z == "no description") z = "" ;
		if (z == "unspecified") z = "" ;
		gsub (" __greater_than__ ", ">", z) ;
		gsub (" __lower_than__ ", "<", z) ;
	    }
	    if (substr(z,1,6) == "E-MTAB")
	    {
		j = index (z, ":") ;
		if (j > 1) z = substr(z, j+1) ;
	    }
	    if (substr (z, 1, 3) == "GSM")
		z = "" ;
	    if (substr (z, 1, 4) == "SRA:")
		z = "" ;
            gsub ("&gt;"," ", z) ;
	    if (length (z) > 0)
		identifier = "Identifier \"Sample name\" \"" z "\"" ;
	}
    }
}
{ # Attributes

    i = index ($0, ">Attributes<") ;
    if (i > 0)
    {
	Biosample_attribute = "Biosample_attribute" ;

	z = substr ($0, i) ;
	i = index (z, "<table") ;
	z = substr (z, i) ;
	j = index (z, "</table>") ;
	zz = substr (z,1,j-1) ;
	 # print z ;
	n = split (zz, aa, "</tr>") ;
	for (k = 1 ; k <= n ; k++) 
	{
	    split (aa[k], bb, "<td>") ;
	    z = bb[1] ;
	    while (substr(z,1,1) == "<")
	    {
		i = index (z, ">")  ;
		if (i > 0)
		    z = substr (z, i+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0)
		z = substr (z, 1, j-1) ;
	    cc1 = z ;

	    z = bb[2] ;
	    while (substr(z,1,1) == "<")
	    {
		i = index (z, ">")  ;
		if (i > 0)
		    z = substr (z, i+1) ;
	    }
	    j = index (z, "<") ;
	    if (j > 0)
		z = substr (z, 1, j-1) ;

	    for (m = 1 ; m <= 5 ; m++)
	    {
		j = length(z) ;
		if (substr(z,j,1) == " ")
		    z = substr (z, 1, j-1) ; 
	    }
	    if (1)
	    {
		if (z == "na") z = "" ;
		if (z == "NA") z = "" ;
		if (z == "N/A") z = "" ;
		if (z == "n/a") z = "" ;
		if (z == "missing") z = "" ;
		if (z == "not applicable") z = "" ;
		if (z == "non provided") z = "" ;
		if (z == "no description") z = "" ;
		if (z == "unspecified") z = "" ;
		gsub (" __greater_than__ ", ">", z) ;
		gsub (" __lower_than__ ", "<", z) ;
		gsub (" __greater_than__", ">", z) ;
		gsub (" __lower_than__", "<", z) ;
		gsub ("__greater_than__", ">", z) ;
		gsub ("__lower_than__", "<", z) ;
		gsub ("\"","",z) ;
		gsub ("'","",z) ;
		gsub ("\\","",z) ;
	    }
	    cc2 = z ;

	    if (length(cc1) > 0)
	    {
		if (cc1 == "ENA-CHECKLIST")
		    continue ;
                if (0 && cc2 == "Drosophila melanogaster")
		    continue ;
		if (cc2 == cc2 + 0)
		    cc2 = cc1  " " cc2 ;
		if (cc1 == "sample comment")
		    description = "Description \""    cc2    "\"" ;
		else
		{
		    Biosample_attribute = Biosample_attribute "\nBiosample_attribute \"" cc1 ;
		    if (length(cc2) > 0) 
			Biosample_attribute =  Biosample_attribute "\" \"" cc2 "\"" ;
		}
	    }
	}

    }
}

{ # decription

    i = index ($0, ">Description<") ;
    if (i > 0)
    {
	z = substr ($0, i+12) ;
	j = index (z, "</p>") ;
	z = substr (z,1,j-1) ;
	while (substr(z,1,1) == "<")
	{
	    i = index (z, ">")  ;
	    if (i > 0)
		z = substr (z, i+1) ;
	}
	j = index (z, "<") ;
	if (j > 0)
	    z = substr (z, 1, j-1) ;
	for (m = 1 ; m <= 5 ; m++)
	{
	    j = length(z) ;
	    if (substr(z,j,1) == " ")
		z = substr (z, 1, j-1) ; 
	}
	    if (1)
	    {
		if (z == "na") z = "" ;
		if (z == "NA") z = "" ;
		if (z == "N/A") z = "" ;
		if (z == "n/a") z = "" ;
		if (z == "missing") z = "" ;
		if (z == "not applicable") z = "" ;
		if (z == "non provided") z = "" ;
		if (z == "no description") z = "" ;
		if (z == "unspecified") z = "" ;
		gsub (" __greater_than__ ", ">", z) ;
		gsub (" __lower_than__ ", "<", z) ;
	    }
	if (length (z) > 1)
	{
	    description = "Description\nDescription \"" z "\"" ;
	}
    }
}
{ # submission

    i = index ($0, ">Submission<") ;
    if (i > 0)
    {
	z = substr ($0, i+11) ;
	j = index (z, "</dd>") ;
	z = substr (z,1,j-1) ;
	while (substr(z,1,1) == "<")
	{
	    i = index (z, ">")  ;
	    if (i > 0)
		z = substr (z, i+1) ;
	}
	j = index (z, "<") ;
	if (j > 0)
	    z = substr (z, 1, j-1) ;
	    if (1)
	    {
		if (z == "na") z = "" ;
		if (z == "NA") z = "" ;
		if (z == "N/A") z = "" ;
		if (z == "n/a") z = "" ;
		if (z == "missing") z = "" ;
		if (z == "not applicable") z = "" ;
		if (z == "non provided") z = "" ;
		if (z == "no description") z = "" ;
		if (z == "unspecified") z = "" ;
		gsub (" __greater_than__ ", ">", z) ;
		gsub (" __lower_than__ ", "<", z) ;
	    }
	if (length (z) > 1)
	{
	    gsub (/EBI; /, "", z) ;
	    submission = "Submission\nSubmission \"" z "\"" ;
	}
    }
}

END {
    if (biosample) printf ("Biosample %s\n", biosample) ;
    printf ("-D Description\n-D Title\n-D Submission\n-D Population_Description\n-D Biosample_attribute\n-D Identifier\n") ;

    if (length (title) > 1)
	print title ;
    if (length (identifier) > 1)
	print identifier ;
    if (index (Biosample_attribute, "\n") > 1)
	print Biosample_attribute ;
    if (index (description, "\n") > 1)
	print description ;
    if (index (submission, "\n") > 1)
	print submission ;
    printf ("\n") ;
}

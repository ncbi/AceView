
{ gsub (/\"/,"\\\"",$0); }
{ gsub (/http:\//,"http:\\/",$0) ;}

{if (srx == "") next;}

{  # is this $&!@& document there is NO way to discover the title, except SRX is in a special font
    i = index ($0, ">"srx"<") ;
    if (i < 1) next ;
    z = substr ($0, i) ;
    j = index (z, "<") ;  # we are garanteed to discover it
    z = substr (z, j) ;
     while (substr(z,1,1) == "<") # a style or an anchor
    {
	i = index (z, ">") ;
	z = substr (z, i+1) ;
	if (i == 0) break ;
    }
    while (n<25)
    {
	c = substr (z,1,1) ;
	if ((c >= "A" && c <= "Z") || (c >= "a" && c <= "z") || (c >= "0" && c <= "9")) 
	    break 
	z = substr (z, 2) ;
	n++ ;
    }
  
    j = index (z, "<") ;
    z = substr (z,1,j-1) ;
    while (substr (z,1,1) == " ") z = substr (z,2) ;
    j = length (z)+1 ;
    while (j>1 && substr (z,j-1,1) == " ") j-- ;
    if (j > 1)
    {
	z = substr (z,1,j-1) ;
	
	i = index (z, "sequencing of ") ;
	if (i > 0)
	{
	    w = substr (z, i+15,4) ;
	    if (( substr(w,1,3) == "DRS" ||  substr(w,1,3) == "SRS" ||  substr(w,1,3) == "GSM")  &&  0 + substr(w,4) > 1)
		z = "" ;
	    if (( substr(w,1,4) == "SAMD" ||  substr(w,1,4) == "SAME" ||  substr(w,1,4) == "SAMN")  &&  0 + substr(w,5) > 1)
		z = "" ;
	}

	if (length (z) > 1)
	title = "Title \"" substr (z,1,j-1)  "\"" ;
    }
}
{   # >Construction protocol: <span>TruSeq DNA LT Sample Prep Kit</
    i = index ($0, ">Construction protocol: <") ;
    if (i > 1)
    {
	z = substr ($0,i+24) ;
	j = index(z, "</span>") ;
	z = substr (z,1,j-1) ;
	
	while (substr(z,1,1) == "<") # a style or an anchor
	{
	    i = index (z, ">") ;
	    z = substr (z, i+1) ;
	    if (i == 0) break ;
	}   
	j = index (z, "<") ;
	if (j > 1) 
	    z = substr (z,1,j-1) ;
	while (substr (z,1,1) == " ") z = substr (z,2) ;
	j = length (z)+1 ;
	while (j>1 && substr (z,j-1,1) == " ") j-- ;
	z = substr (z,1,j-1) ;
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
	}
	if (length (z) > 1)
	    cons_prot = "Construction_protocol \""z"\"" ;
    }
}

{   # >Construction protocol: <span>TruSeq DNA LT Sample Prep Kit</
    i = index ($0, ">Design: <") ;
    if (i > 1)
    {
	z = substr ($0,i+9) ;

	j = index(z, "</span>") ;
	z = substr (z,1,j-1) ;
	while (substr(z,1,1) == "<") # a style or an anchor
	{
	    i = index (z, ">") ;
	    z = substr (z, i+1) ;
	    if (i == 0) break ;
	}   
	j = index (z, "<") ;
	if (j > 1) 
	    z = substr (z,1,j-1) ;
	while (substr (z,1,1) == " ") z = substr (z,2) ;
	j = length (z)+1 ;
	while (j>1 && substr (z,j-1,1) == " ") j-- ;
	z = substr (z,1,j-1) ;	 
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
	}
	if (length (z) > 1)
	    design = "Design \""z"\"" ;
    }
}
{   # >Construction protocol: <span>TruSeq DNA LT Sample Prep Kit</
    i = index ($0, ">Submitted by: <") ;
    if (i > 1)
    {
	z = substr ($0,i+15) ;
	j = index(z, "</span>") ;
	z = substr (z,1,j-1) ;
	
	while (substr(z,1,1) == "<") # a style or an anchor
	{
	    i = index (z, ">") ;
	    z = substr (z, i+1) ;
	    if (i == 0) break ;
	}   
	j = index (z, "<") ;
	if (j > 1) 
	    z = substr (z,1,j-1) ;
	while (substr (z,1,1) == " ") z = substr (z,2) ;
	j = length (z)+1 ;
	while (j>1 && substr (z,j-1,1) == " ") j-- ;
	z = substr (z,1,j-1) ;
	if (substr (z,1,1) == "(" && substr(z, j-1,1) == ")")
	    z = substr (z,2,j-3) ;

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
	    if (z == "Gene Expression Omnibus (GEO)") z = "" 
	    if (z == "UC SANTA CRUZ") z = "" ;
	}
	
	if (length (z) > 1)
	    author = "Submitted_by \""z"\"" ;
    }
}

END {
    if (srx) printf ("SRX %s\n", srx) ;
    printf ("-D Design\n-D Submitted_by\n-D Construction_protocol\n-D Title\n") ;

    if (length (title) > 1)
	print title ;
    if (length (design) > 1)
	print design ;
    if (length (cons_prot) > 1)
	print cons_prot ;
    if (length (dscr) > 1)
	print dscr ;
    if (length (author) > 1)
	print author ;
    if (length (submission) > 1)
	print submission ;
    printf ("\n") ;
}

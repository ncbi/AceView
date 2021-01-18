function cleanSpace(_s) {
    _s2 = "" ;
    _iMax = length(_s) ;
    for (_i = 1 ; _i <= _iMax ; _i++)
    {
	_z = substr(_s,_i) ;
	_j = index(_z,"  ") ;
	if (_j == 0)
	{
	    _s2 = _s2 _z ;
	    break ;
	}
	_s2 = _s2 substr(_z,1,_j) ;
	_j++ ;
	while (substr(_s2,_j,1) == " ")
	    _j++ ;
	_i = _i + _j - 2 ;
    } 
    while (substr (_s2,1,1) == " ") _s2 = substr (_s2, 2) ;
    while (substr(_s2,1,1) == "<") # jump the decorations 
    {  
	i = index (_s2, ">") ;  
	if (i < 1) i = 1 ; # avoid loops in case of misconfigured html
	_s2 = substr (_s2, i + 1) ;
	while (substr (_s2,1,1) == " ") _s2 = substr (_s2, 2) ;
    }

    return _s2  ;
}
{ gsub (/\"/,"\\\"",$0); }
{ gsub (/http:\//,"http:\\/",$0) ;}

/^<h1>/ {
    i = index($0,"</h1>") ;
    t1 = cleanSpace(substr($0,5,i-5)) ;
}
/<STUDY_ABSTRACT/ {
    i = index($0,"<STUDY_ABSTRACT>") ;  

    j = index($0,"</STUDY_ABSTRACT>") ;
    z = substr($0,i+16,j-i-16) ;

    z = substr(z, 1, j-1) ;
   j = index (z,"<") ; 
 
    z  = cleanSpace(z) ;

    while (length(z) > 1)
    {   
	c = substr (z,1,1) ;
	if (c == "<" || (c >= "A" && c <= "Z") || (c >= "a" && c <= "z") || (c >= "0" && c <= "9")) 
	    break ;
	z = substr (z, 2) ;
    }
    ab = cleanSpace(z) ;
    if (0)
	while (j > 0 && j < 5)
	{   # this is a hack, NCBI abstract have often several non ascii char before a <class> construct
	    ab = cleanSpace(substr(ab, j)) ;
	    j = index (ab,"<") ;
	}
}
/<STUDY_DESCRIPTION/ {
    i = index($0,"<STUDY_DESCRIPTION>") ;
    j = index($0,"</STUDY_DESCRIPTION>") ;
    z = substr($0,i+19,j-i-19) ;
    z  = cleanSpace(z) ; 
    if (z != "none provided") 
	dcr = z ;
}
/Show experiment information in SRA Entrez/ {
    i = index($0,"Show experiment information in SRA Entrez") ;
    z = substr($0,i+20) ;
    i = index(z,"\">") ;
    z = substr(z,i+2) ;
    i = index(z,"</a>") ;
    k = 0+substr(z,1,i-1) ;
    if(k > 0) 
    { 
	if (substr(t1,length(t1)) ==  ".") 
	    t1 = substr(t1,1, length(t1) - 1) ;
	t1  =  t1 " (" k " exp)." ;
    }
}
/<span>Identifiers:/ {
    ok = 0 ; z = $0 ;
    while (ok == 0)
    {
	i = index(z,"<em>") ; if (i == 0) next ;
	z = substr (z, i+4) ;
	j = index(z,":") ; if (j == 0) next ; 
	idx = substr(z, 1, j) ;
	idx =  cleanSpace(idx) ; 
	id0 = index (idx, "GEO") ;
	id1 = index (idx, "SRA") ;
	id2 =  index(idx, "BioProject") ;
	id3 = index (idx, "NCBI") ;
	id4 = index (idx, "DDBJ") ;
	id5 = index (idx, "dbGaP") ;
	id6 = index (idx, "UWGS-RW") ;
	if (id0 > 0)
	{ 
	    z = substr (z, j) ;
	    i = index(z, "acc=") ; if (i > 0)  z = substr (z, i+4) ;
	    j =  index(z, "\"") ;
	    if (j > 1) 
	    { geo = substr (z,1,j-2) ; z = substr (z, j-1) ; }
	}
	else if (id1 + id2 + id6 == 0)
	{ ok = 2 ; id = idx ; } # Good identifier found
	else if (id3 + id4 + id5 > 0)
	{ 
	    i = index (z, ">") ; z = substr(z, i+1) ;#  print "##"z ;
	    j = index (z, "<") ; z = substr(z, 1, j-1) ;
	    if (j > 1) {  ok = 2 ; id = z ; } # Good identifier found
	}
	else
	{ z = substr (z, j) ; }
    }
}
/<span>External Link.*pubmed/{
    nn = split ($0,aa,"pubmed:") ;
    for (ii = 2 ; ii <= nn ; ii++)
    {
	z = aa[ii] ;
	gsub (/ /,"",z) ;
	j = index (z, "<") ;
	if (j > 1) 
	   pm = pm "Reference pm" substr(z,1,j-1) "\n" ;
    }
}
END {
    if (ab && dcr && index(ab,dcr) > 0)
	dcr  =  ""  ;
    printf("SRP %s\n-D Title\n-D Abstract\n-D Description\n-D Identifier\n-D Reference\n-D GEO\n",srp) ;
    
    if (length(t1) > 1)
	printf("Title \"%s\"\n",t1) ;
    if (length(ab) > 1)
	printf("Abstract \"%s\"\n", srp) ;
    if (length(dcr) > 1)
	printf("Description \"%s\"\n",dcr) ;
    if (length(id) > 1)
	printf("Identifier \"%s\"\n",id) ;
    if (length(geo) > 1)
	printf("GEO \"%s\"\n",geo) ;
    if (length(pm) > 1)
	print pm ;
    printf("\n") ;
    if (length(ab) > 1)
    { 
	printf("LongText \"%s\"\n",srp) ;
	print ab ;
        printf("***LongTextEnd***\n\n") ;
    }
}

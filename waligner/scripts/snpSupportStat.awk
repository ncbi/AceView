BEGIN {      
  }
/^#/ { next }
/^ZZZZZ/{zz++ ; next ;}
{
  if (irMax == 0)
    {
      run = "Sum" ; irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; 
      run = "Union" ; irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; 
    }

  if (zz < 1)
    { 
      r = $1 ; ir=r2ir[r];if(ir<1){irMax++;r2ir[r]=irMax;ir2r[irMax]=r;ir=irMax;}
      next ;
    }
}
{ 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
}
{
  if (zz == 1) 
    { 
      r = $1 ; ir=r2ir[r];if(ir<1)next ;
      r = $1 ; ir=r2ir[r];if(ir<1){irMax++;r2ir[r]=irMax;ir2r[irMax]=r;ir=irMax;}
      gsub (/System:/,"",$0);
      gsub(/NULL/,"NA",$2) ; gsub(/ /,"_",$2) ; machine[ir]=$2; if ($2 != "NA") nmachine = 1 ;
      gsub(/NULL/,"NA",$3) ; gsub(/ /,"_",$3) ; if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ;if ($3 != "NA") nsample = 1 ;
      gsub(/NULL/,"NA",$4) ; gsub(/ /,"_",$4) ; if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; if ($4 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$7) ; gsub(/ /,"_",$7) ; if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ;if ($7 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$5) ; gsub(/ /,"_",$5) ; if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ;if ($5 != "NA") ntissue = 1 ;
      gsub(/NULL/,"NA",$6) ; gsub(/ /,"_",$6) ; title[ir]=$6;if ($6 != "NA") ntitle = 1 ;
      gsub(/NULL/,"NA",$8) ; gsub(/ /,"_",$8) ; runid[ir]=$8; if ($8 != "NA") nrunid = 1 ;
      gsub(/NULL/,"NA",$11) ; gsub(/ /,"_",$10) ; stitle[ir]=$11; if ($11 != "NA") nstitle = 1 ; # Sorting_title
      gsub(/NULL/,"NA",$12) ; gsub(/ /,"_",$10) ; stitle2[ir]=$12; if ($12 != "NA") nstitle2 = 1 ; # Sorting_title_2
      gsub(/NULL/,"NA",$13) ; gsub(/ /,"_",$10) ; otitle[ir]=$13; if ($13 != "NA") notitle = 1 ; # Other_title

      next;
    }
}


{
    mult = 1 ; # one snp per line
    if (format == "snp")
    {
      chrom = $1 ; pos = $2;
      run = $6 ;
      ir=r2ir[run];if(ir<1)next ;
      ir = r2ir[run] ; if (ir<1) { irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; }
      a = $5 ;
      b = $6 ;
      if (a == ".")
	z = "Ins " b ;
      else if  (b == ".")
	z = "Del " a  ;
      else if (length (a)+length (b) == 2)
	z =  a ">" b  ;
      else 
	z =  "Other" ;
      mut = $10 ; wild = $11 ; cover = mut + wild ;
      if (cover < 1) next ;
      if (mut < minMut || cover < minCover || 100 * mut < minF * cover) next ;
    }
  else if (format == "txt" || format == "snpq")
    {
      chrom = $1 ; pos = $2;
      run = $6 ;
      ir=r2ir[run];if(ir<1)next ;
      ir = r2ir[run] ; if (ir<1) { irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; }
      z = $3 ;
      gsub (/Ins/, "Ins ", z) ;
      gsub (/Del/, "Del ", z) ;
      mut = 0 + $10 ; wild = 0 + $11 ; 
      cover = mut + wild ;
      # print "XXXXXX", cover,mut,wild ;
      if (cover < 1) next ;
      if (mut < minMut || cover < minCover || 100 * mut < minF * cover) next ;
    }
  else if (format == "Error_type")
    {
      gsub (/\"/,"",$0) ;
      run = $1 ;
      ir=r2ir[run];if(ir<1)next ;
      chrom = 1 ; pos = 1;      z = $3 ; if (substr(z,1,1)=="*") z = substr(z,2) ; 
      gsub (/a/,"A",z) ; gsub (/t/,"T",z) ; gsub (/g/,"G",z) ; gsub (/c/,"C",z) ; 
      gsub(/---/,"--",z) ;      gsub(/--/,"-",z) ;      gsub(/-/,"Del ",z) ;
      gsub(/+++/,"++",z) ;      gsub(/++/,"+",z) ;      gsub(/+/,"Ins ",z) ;
      gsub(/2/,">",z) ;
      ir = r2ir[run] ; if (ir<1) { irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; }
      mut = $4 ; wild = 0 ;
      cover = mut + wild ;
      if (cover < 1) next ;
      if (mut < minMut || cover < minCover || 100 * mut < minF * cover) next ;
    }
   else if (format == "stats")
    {
      gsub (/\"/,"",$0) ;
      run = $1 ;
      ir=r2ir[run];if(ir<1)next ;
      chrom = 1 ; pos = 1;      z = $2 ; if (substr(z,1,1)=="*") z = substr(z,2) ; 
      if (substr(z,1,3) == "Del") z = "Del " substr(z,4) ;
      if (substr(z,1,3) == "Ins") z = "Ins " substr(z,4) ;
      mult = $5 ;  # so many SNPs of this type were counted
      mut = $6 ; wild = 7 ; # covered by so many mutant/wild type reads
      f = $3 ; cover = $4 ;
       if (cover < minCover) next ;
      if (f < minF) next ;
    }
  else if (format == "countEdited")
    {
      chrom = $1 ; pos +=  1 ;
      run = $8 ;
      ir=r2ir[run];if(ir<1)next ;
      ir = r2ir[run] ; if (ir<1) { irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; }
      genotype = $9 ; mut = $12 ; wild = $13 ; cover = mut + wild ;
      if (cover < 1) next ;
      if (mut < minMut || cover < minCover || 100 * mut < minF * cover) next ;


      if (index($6, "n") > 0) next ;
      if (index($7, "n") > 0) next ;

      if (substr($5,2,1) == ">")
	z = $5 ;
      else if (substr($5,1,3) == "mIS")
	z = "Ins " $7 ;
      else if (substr($5,1,3) == "mDL")
	z = "Del " $6 ;
      else 
	z =  "Other" ;
    }
  else
    {
      print "Bad format parameter in snpSupportStat.awk format =#" format "#" ;
      exit (1) ;
    } 
  m[z,ir] += mut ;
  w[z,ir] += wild ;
  c[z,ir]+=mult ; 
  anyC[ir]+=mult ;

  if (1)  # sum and union 
    {
      m[z,1] += mut ;
      w[z,1] += wild ;
      
      m[z,2] += mut ;
      w[z,2] += wild ;
    }
  c[z,1]+=mult ; 
  anyC[1]+=mult ;
  
  crr[z, chrom pos]+=mult ; if (crr[z, chrom pos] == mult) { c[z,2]++ ; cr[z]++ ; anyC[2]++ ; }
}
END  {
   Types = "Any,Substitution,Transition,Transversion,Insertion,Deletion,Double insertion,Double deletion,Triple insertion,Triple deletion,Other,A>G,T>C,G>A,C>T,A>T,T>A,G>C,C>G,A>C,T>G,G>T,C>A,Ins A,Ins T,Ins G,Ins C,Del A,Del T,Del G,Del C,Ins AA,Ins TT,Ins GG,Ins CC,Ins AG,Ins CT,Ins AC,Ins GT,Ins TG,Ins CA,Ins TC,Ins GA,Ins AT,Ins TA,Ins GC,Ins CG,Del AA,Del TT,Del GG,Del CC,Del AG,Del CT,Del AC,Del GT,Del TG,Del CA,Del TC,Del GA,Del AT,Del TA,Del GC,Del CG,Ins AAA,Ins TTT,Ins GGG,Ins CCC,Ins AAT,Ins ATT,Ins AAG,Ins CTT,Ins AAC,Ins GTT,Ins TTA,Ins TAA,Ins TTG,Ins CAA,Ins TTC,Ins GAA,Ins GGA,Ins TCC,Ins GGT,Ins ACC,Ins GGC,Ins GCC,Ins CCA,Ins TGG,Ins CCT,Ins AGG,Ins CCG,Ins CGG,Ins ATA,Ins TAT,Ins ATG,Ins CAT,Ins ATC,Ins GAT,Ins AGA,Ins TCT,Ins AGT,Ins ACT,Ins AGC,Ins GCT,Ins ACA,Ins TGT,Ins ACG,Ins CGT,Ins TAG,Ins CTA,Ins TAC,Ins GTA,Ins TGA,Ins TCA,Ins TGC,Ins GCA,Ins TCG,Ins CGA,Ins GAG,Ins CTC,Ins GAC,Ins GTC,Ins GTG,Ins CAC,Ins GCG,Ins CGC,Ins CAG,Ins CTG,Del AAA,Del TTT,Del GGG,Del CCC,Del AAT,Del ATT,Del AAG,Del CTT,Del AAC,Del GTT,Del TTA,Del TAA,Del TTG,Del CAA,Del TTC,Del GAA,Del GGA,Del TCC,Del GGT,Del ACC,Del GGC,Del GCC,Del CCA,Del TGG,Del CCT,Del AGG,Del CCG,Del CGG,Del ATA,Del TAT,Del ATG,Del CAT,Del ATC,Del GAT,Del AGA,Del TCT,Del AGT,Del ACT,Del AGC,Del GCT,Del ACA,Del TGT,Del ACG,Del CGT,Del TAG,Del CTA,Del TAC,Del GTA,Del TGA,Del TCA,Del TGC,Del GCA,Del TCG,Del CGA,Del GAG,Del CTC,Del GAC,Del GTC,Del GTG,Del CAC,Del GCG,Del CGC,Del CAG,Del CTG,Ambiguous" ;         
  nTypes = split  (Types, aaa, ",") ;
  for (uu = 1 ; uu<= 6 ; uu++)
    {
      if (format == "Error_type")
	{
	  if (0)
	    continue ;
	}
      if (uu == 1)
	printf ("Number of variants") ;
      if (uu == 2)
	printf ("\t\tNumber of reads supporting the variants") ;
      if (uu == 3)
	printf("\t\tNumber of reads supporting the reference") ;
      if (uu == 4)
	printf("\t\tPrevalence of the variant type") ;
      if (uu == 5)
	printf("\t\tAverage  number of reads supporting the variant") ;
      if (uu == 6)
	printf("\t\tAverage  number of reads supporting the reference") ;

      if (nrunid) printf("\tRunId") ;
      if (nmachine == 1) printf("\tMachine") ;
      if (nsample == 1)  printf("\tSample") ;
      if (nsystm == 1)  printf("\tSystem") ;
      if (ntissue == 1) printf("\tTissue") ;
      if (ntitle == 1)  printf("\tTitle") ;
      if (nstitle == 1)  printf("\tSorting title") ;
      if (nstitle2 == 1)  printf("\tSorting title 2") ;
      if (notitle == 1)  printf("\tOther title") ;

      for( t = 1; t <= nTypes ; t++)
	printf("\t%s", aaa[t]);
    }
  printf("\n");


for (ir = 1 ; ir <= irMax ; ir++)
 {
   if (format == "Error_type")
     anyC[ir] = 100 ;
   for (uu = 1 ; uu<= 6 ; uu++)
     {
       if (format == "Error_type")
	 {
	   if (0)
	     if (uu == 1 || uu == 3 || uu == 5 || uu == 6)
	       continue ;
	 }
       for (t =  1 ; t <=  10  ; t++)
	 n[aaa[t]] = 0 ;
       if (uu == 1)
	 {              # number of variants of each type
	   for (z in cr)
	     n[z] = c[z,ir] ;
	 }
       if (uu == 2)
	 {             # number of reads supporting the mutant
	   for (z in cr)
	     {
	       n[z] = m[z,ir] ; 
	       if (c[z,ir] == 0)
		 c[z,ir] = 1 ;
	     }
	 }
       if (uu == 3)
	{             # number of reads supporting the wild type
	  for (z in cr)
	    n[z] = w[z,ir] ;
	}
       if (uu == 4)
	 {             # fraction of mutant of each type
	   if (format == "Error_type")
	     {
	       for (z in cr)
		 {
		   n[z] = mAny[z] ;
		 }
	     }
	   else
	     {
	       for (z in cr)
		 n[z] = c[z,ir] ;
	     }
	 }
       if (uu == 5)
	 {            # average number of reads supporting the mutant
	   for (z in cr)
	     n[z] = m[z,ir] ;
	 }
       if (uu == 6)
	 {      # average number of reads supporting the wild type
	   for (z in cr)
	     n[z] = w[z,ir] ;
	 }
       
       n[aaa[1]] = 0 ;
       for (t =  11 ; t <=  nTypes  ; t++)
	 {
	   n[aaa[1]] += n[aaa[t]] ;
	 }  
       for (t =  12 ; t <=  23  ; t++)
	 {
	   n[aaa[2]] += n[aaa[t]] ;
	 }
       for (t =  12 ; t <=  15  ; t++)
	 {
	   n[aaa[3]] += n[aaa[t]] ;
	 }
       for (t =  16 ; t <=  23  ; t++)
	 {
	   n[aaa[4]] += n[aaa[t]] ;
	 }
       for (t =  24 ; t <=  27  ; t++)
	 {
	   n[aaa[5]] += n[aaa[t]] ;
	 }
       for (t =  28 ; t <=  31  ; t++)
	 {
	   n[aaa[6]] += n[aaa[t]] ;
	 }
       for (t =  32 ; t <=  47  ; t++)
	 {
	   n[aaa[7]] += n[aaa[t]] ;
	 }
       for (t =  48 ; t <=  63  ; t++)
	 {
	   n[aaa[8]] += n[aaa[t]] ;
	 }
       for (t =  64 ; t <=  127  ; t++)
	 {
	   n[aaa[9]] += n[aaa[t]] ;
	 }
       for (t =  128 ; t <=  191  ; t++)
	 {
	   n[aaa[10]] += n[aaa[t]] ;
	 }
       
       if (uu>1) 
	 printf ("\t\t") ;
       printf ("%s",ir2r[ir]) ; 
       if (nrunid) printf("\t%s", runid[ir]) ;
       if (nmachine == 1) printf("\t%s", machine[ir]) ;
       if (nsample == 1)  printf("\t%s", sample[ir]) ;
       if (nsystm == 1)  printf("\t%s", systm[ir]) ;
       if (ntissue == 1) printf("\t%s", tissue[ir]) ;
       if (ntitle == 1)  printf("\t%s", title[ir]) ;  
       if (nstitle == 1)  printf("\t%s", stitle[ir]) ;  
       if (nstitle2 == 1)  printf("\t%s", stitle2[ir]) ;  
       if (notitle == 1)  printf("\t%s", otitle[ir]) ;  

       for (t =  1 ; t <=  nTypes  ; t++)
	 {
	   if (uu == 1)
	     nTR[aaa[t],ir] = n[aaa[t]] ;
	   if (uu == 2 && t >= 11)
	     {
	       u = n[aaa[1]] ; if (u == 0) u = 1 ;
	       mAny[aaa[t]] = 100*n[aaa[t]]/u ;
	     }
	   if (uu < 4)
	     printf ("\t%d", n[aaa[t]]) ;
	   else if (uu == 4)
	     {
	       u = anyC[ir] ; if (u == 0) u = 1 ;
	       printf ("\t%.4f", 100 * n[aaa[t]]/u) ;
	     }
	   else if (uu >= 5)
	     { 
	       u = nTR[aaa[t],ir] ; if (u == 0) u = 1 ;
	       printf ("\t%.2f", n[aaa[t]]/u) ;
	     }
	 }
     }
   printf ("\n") ;
 }
}

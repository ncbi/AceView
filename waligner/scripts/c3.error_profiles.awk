{gsub (/\"/,"",$0);}
/^Error_profile/{
  r = 0 ; t=$2; 
  if ($2 == "Any") 
    {
      nAli += $4 ; # in Mb 
      nn[t] += $3;
      next ;
    }
  if (index (t, "n") > 0) next ;
  t1 = t ; gsub(/[at*g>c+-]/,"",t1); if (length(t1)>0) next;
  if (substr(t,1,1) == "*")
    {
      t = substr(t,2) ; r = 1;
    }

  isI = 10000 * index (t,"+++")  + 100 * index (t,"++")  + index (t,"+") ;
  isD = 10000 * index (t,"---")  + 100 * index (t,"--")  + index (t,"-") ;

  gsub (/a/,"A",t) ;  gsub (/t/,"T",t) ;  gsub (/g/,"G",t) ;  gsub (/c/,"C",t) ;
  gsub (/+++/, "Ins ", t) ; gsub (/++/, "Ins ", t) ; gsub (/+/, "Ins ", t) ;
  gsub (/---/, "Del ", t) ; gsub (/--/, "Del ", t) ; gsub (/-/, "Del ", t) ;
  if (r == 1)
    {
      nr[t] += $3;
      nr["Any"] += $3 ;

      if (isD >= 10000) nr["Triple deletion"] += $3 ;
      else if (isD >= 100) nr["Double deletion"] += $3 ;
      else if (isD >= 1) nr["Deletion"] += $3 ;

      if (isI >= 10000) nr["Triple insertion"] += $3 ;
      else if (isI >= 100) nr["Double insertion"] += $3 ;
      else if (isI >= 1) nr["Insertion"] += $3 ;

      if (isI + isD == 0) 
	{
	  nr["Substitution"] += $3 ;
          if (t == "A>G" || t == "G>A" || t == "T>C" || t == "C>T")
	    nr["Transition"] += $3 ;
	  else
	    nr["Transversion"] += $3 ;
	}
    }
  nn[t] += $3;

      if (isD >= 10000) nn["Triple deletion"] += $3 ;
      else if (isD >= 100) nn["Double deletion"] += $3 ;
      else if (isD >= 1) nn["Deletion"] += $3 ;

      if (isI >= 10000) nn["Triple insertion"] += $3 ;
      else if (isI >= 100) nn["Double insertion"] += $3 ;
      else if (isI >= 1) nn["Insertion"] += $3 ;

      if (isI + isD == 0)
	{
	  nn["Substitution"] += $3 ;
          if (t == "A>G" || t == "G>A" || t == "T>C" || t == "C>T")
	    nn["Transition"] += $3 ;
	  else
	    nn["Transversion"] += $3 ;
	}
  nnn += $3 ;

}
END{
  if (nAli == 0) nAli = 1 ;
  if (nnn == 0) nnn = 1 ;
  printf (" evaluated on %d aligned Mb", nAli) ;
  for (i=1;i<5;i++)printf("\t%s [%dMb]",title,nAli);
  printf ("\nType\tMismatches per million aligned bases\tFrequency\t%% indel identical to a neighbor\tNumber observed\n") ;
  t = "Any" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
  t = "Substitution" ; 
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;

  t = "Transversion" ; 
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;

  t = "Transition" ; 
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;

  t = "Insertion" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
  t = "Deletion" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;

  t = "Double insertion" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
  t = "Double deletion" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;

  t = "Triple insertion" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
  t = "Triple deletion" ;
  n = nn[t] ; if (n == 0) n = 1 ;
  printf ("%s\t%.1f\t%.2f\t%.2f\t%d\n", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;

  if (nnn == 0) nnn = 1 ;
  transits = "A>T T>A G>C C>G A>C T>G G>T C>A A>G T>C G>A C>T"

  for (i = 1 ; i < 48 ; i += 4)
    { 
      t = substr (transits, i, 3) ;  
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }
  for (j = 1 ; j <= 4 ; j++)
    { 
      if (i == j) continue ;
      b = substr ("ATGC", j, 1) ;
      t = "Ins " b ;
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }
  for (j = 1 ; j <= 4 ; j++)
    { 
      if (i == j) continue ;
      b = substr ("ATGC", j, 1) ;
      t = "Del " b ;
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }


  doublets = "AA TT GG CC AT TA GC CG AG CT AC GT TG CA TC GA "
  for (i = 1 ; i <= 46 ; i+=3)
    { 
      ab = substr (doublets, i, 2) ;
      t = "Ins " ab ;
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }
  for (i = 1 ; i <= 46 ; i+=3)
    { 
      ab = substr (doublets, i, 2) ;
      t = "Del " ab ;
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }

  triplets = "AAA TTT GGG CCC AAT AAG AAC TTA TTG TTC GGA GGT GGC CCA CCT CCG ATA ATT ATG ATC AGA AGT AGG AGC ACA ACT ACG ACC TAA TAT TAG TAC TGA TGT TGG TGC TCA TCT TCG TCC GAA GAT GAG GAC GTA GTT GTG GTC GCA GCT GCG GCC CAA CAT CAG CAC CTA CTT CTG CTC CGA CGT CGG CGC" ;
  for (i = 1 ; i <= 4*64-4 ; i+=4)
    { 
      abc = substr (triplets, i, 3) ;
      t = "Ins " abc ;
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }
  for (i = 1 ; i <= 4*64-4 ; i+=4)
    { 
      abc = substr (triplets, i, 3) ;
      t = "Del " abc ;
      n = nn[t] ; if (n == 0) n = 1 ;
      printf ("%s\t%.1f\t%.2f\t%.2f\t%d", t, nn[t]/nAli, 100 * nn[t]/nnn, 100 * nr[t]/n, nn[t]) ;
      printf ("\n") ;
    }



}

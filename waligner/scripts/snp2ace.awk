/^#/ { next }
{
  chrom = $1 ; pos = $2; 
  run = $6 ;
  ir = r2ir[run] ; if (ir<1) { irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; }
  genotype = $7 ; f = $8 ; cover = $9 ; m = $10 ; w = $11 ; mp = $12 ; wp = $13 ; mm = $14 ; wm = $15 ; bad = $16
  
  if (index($4, "n") > 0) next ;
  if (index($5, "n") > 0) next ;
  if (index($4, "N") > 0) next ;
  if (index($5, "N") > 0) next ;
  
  if (substr($3,2,1) == ">")
    { typ = substr($3,1,1) "2" substr($3,3,1)  ; typ2 = typ ; typ3 = typ ; }
  else if (substr($3,1,3) == "Ins")
    { typ = $3 ; typ2 = typ ; typ3 = typ  ; if (length($3) > 4) typ3 = "Multi_insertion " typ ; } 
  else if (substr($3,1,3) == "Del")
    { typ = $3 ; typ2 = typ ; typ3 = typ  ; if (length($3) > 4) typ3 = "Multi_deletion " typ ; }
  else if (substr($3,1,3) == "Dup")
    { ln = pos1 - pos2 ; typ = $3 ln ;typ2 = typ ; typ3 = "Multi_insertion Duplication_" ln "_bp" ; } 
  else 
    { typ = $3 ; typ2 = type ; typ3 = 0 ; }
  
  if (genotype == "+/+") geno2 = "ww" ;
  else if (genotype == "m/m") geno2 = "mm" ;
  else if (genotype == "m/+") geno2 = "wm" ;
  else if (genotype == "-") next ;
  else if (genotype == "+/?") geno2 = "wx" ;
  else if (genotype == "intermediate") geno2 = "Intermediate" ;
  else if (genotype == "m/?") geno2 = "mx" ;


  if (genotype == 0) geno2 = "ww" ;
  else if (genotype == 2) geno2 = "mm" ;
  else if (genotype == 0.5) geno2 = "wx" ;
  else if (genotype == "intermediate") geno2 = "Intermediate" ;
  else if (genotype == 1) geno2 = "wm" ;
  else if (genotype == 1.5) geno2 = "mx" ;
  else if (genotype == 1.3) geno2 = "Intermediate" ;
  else if (genotype == 0.3) geno2 = "Intermediate" ;


  v = chrom ":" pos ":" typ ; pos1 = pos+1 ;
  txt = "-D " targetType "\n" targetType  " " chrom " " pos " 1\ntyp " typ2  ;
  if (typ3 != 0) txt = txt  "\n" typ3  ; 

  if (format == "showCounts") 
  {
      txt2 = geno2 " " run " " f " " cover " " m " " w " " mp " " wp " " mm " " wm " " bad ;
      if ($19 != "-")  { txt2 = txt2  "\n" $19  ; if (index ($19, "100")> 0) txt2 = txt2 " " run ; }
      # if ($20 != "-")  txt2 = txt2  "\n" $20 " " run ;
      counts[ir,v] = txt2 ;
  }
  typV[v] = txt ;
}
END  {
  for (v in typV)
    {
      printf ("Variant \"%s\"\n%s\n", v, typV[v]) ;
      if (targetType == "IntMap") printf ("Found_in_genome\n") ;
      if (targetType == "mRNA") printf ("Found_in_mRNA\n") ;

       if (format == "showCounts") 
	   for (ir = 1 ; ir <= irMax ; ir++)
	   {
	       if (counts[ir,v])
		   printf ("FQ %s\n", counts[ir,v]) ;
	   }
       printf ("\n") ;
    }
}

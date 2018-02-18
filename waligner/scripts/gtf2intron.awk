{ gsub (/\"/,"",$0) ;
  chr = $1 ; gsub (/chr/,"",chr) ; a1 = $4 ; a2 = $5 ; str = $7 ;
  txt = $9 ; split (txt, aa, ";") ;
  oId = 0 ; i = index (txt, "gene_name") ; if (i > 0) { s = substr (txt, i+ 10) ; split (s, aa, ";") ; oId = aa[1] ; }
  if (oId == old)
    {
      if (str == "+") { u1 = b2 + 1 ; u2 = a1 - 1 ; }
      else { u1 = a1 - 1 ; u2 = b2 + 1 ; }
      printf ("Intron %s__%d_%d // %s\n", chr, u1, u2, oId) ;
    }
  old = oId ; b1 = a1 ; b2 = a2 ;
}

  

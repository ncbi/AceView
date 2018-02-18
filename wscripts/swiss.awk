# File wscripts/swiss.awk, postprocessing commands used by blast_search script
# @(#)swiss.awk	1.1 3/28/97

BEGIN { state = 0 ;}
/^AC/ {seq = "SP:"$2 ; gsub(/;/,"",seq) ; printf("\nProtein %s\nDatabase SWISSPROT %s\n",seq,seq); next;}
/^DE/ {xx = substr($0,6) ; printf("Description \"%s\n", xx) ; next ;}
/^RX/ {xx = substr($0,6) ; printf("Reference \"%s\n", xx) ; next ;}
/^RL/ {xx = substr($0,6) ; printf("Reference \"%s\n", xx) ; next ;}
/^GN/ {xx = substr($0,6) ;  gsub(/\./,"",xx) ; gsub (/ OR /,"\nLocus \"",xx) ; printf("Locus \"%s\n", xx) ; next ;}
/^\/\// {state = 0 ; printf("\n") ; next; }

/^SQ/ {printf ("Peptide %s\n\n",seq) ; state = 2 ; next ;}
{ if (state == 2) { printf("Peptide %s\n", seq) ; state = 3 ;  }
  if (state == 3) 
    { pp = $0 ; gsub(/[0-9]/,"",pp) ; gsub(/ /,"",pp) ;
      printf("%s\n", pp) ;
    }
}


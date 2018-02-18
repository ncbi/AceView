# File wscripts/gp.awk, postprocessing commands used by blast_search script
# $Id: gp.awk,v 1.1.1.1 2002/07/19 20:23:33 sienkiew Exp $

BEGIN { state = 0 ;}
/^ACCESSION/ {seq = "GP:"$2 ; printf("\nProtein %s\nDatabase GP %s\n",seq,$2); next;}
/^TITLE/ {printf("Title \"%s\n", substr($0,6)) ; next ;}
/^\/\/\// {state = 0 ; printf("\n") ; next; }

/^SEQUENCE/ {printf ("Peptide %s\n\n",seq) ; state = 2 ; next ;}
{ if (state == 2) { printf("Peptide %s\n", seq) ; state = 3 ; next ; }
  if (state == 3) 
    { pp = $0 ; gsub(/[0-9]/,"",pp) ; gsub(/ /,"",pp) ;
      printf("%s\n", pp) ;
    }
}


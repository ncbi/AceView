# File wscripts/gb.awk, postprocessing commands used by blast_search script
# $Id: gb.awk,v 1.1.1.1 2002/07/19 20:23:33 sienkiew Exp $

BEGIN { state = 0 ;}
/^ACCESSION/ {seq = "GB:"$2 ; printf("\nSequence %s\nDatabase GenBank %s\n",seq,$2); next;}
/^TITLE/ {printf("Title \"%s\n", substr($0,6)) ; next ;}
/^\/\// {state = 0 ; printf("\n") ; next; }

/^ORIGIN/ {printf ("DNA %s\n\n",seq) ; state = 2 ; next ;}
{ if (state == 2) { printf("DNA %s\n", seq) ; state = 3 ; next ; }
  if (state == 3) 
    { pp = $0 ; gsub(/[0-9]/,"",pp) ; gsub(/ /,"",pp) ; gsub(/\//,"",pp) ;
      printf("%s\n", pp) ;
    }
}


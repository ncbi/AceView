# File wscripts/emb.awk, postprocessing commands used by blast_search script
# %W% %G%

BEGIN { state = 0 ;}
/^ACCESSION/ {seq = "EMBL:"$2 ; printf("\nSequence %s\nDatabase EMBL %s\n",seq,$2); next;}
/^TITLE/ {printf("Title \"%s\n", substr($0,6)) ; next ;}
/^\/\/\// {state = 0 ; printf("\n") ; next; }

/^SEQUENCE_JUNK/ {printf ("Sequence %s\n\n",seq) ; state = 2 ; next ;}
{ if (state == 2) { printf("DNA %s\n", seq) ; state = 3 ; next ; }
  if (state == 3) 
    { pp = $0 ; gsub(/[0-9]/,"",pp) ; gsub(/ /,"",pp) ; gsub(/\//,"",pp) ;
      printf("%s\n", pp) ;
    }
}


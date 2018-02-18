# File wscripts/swissbig.awk, postprocessing commands used by blast_search script
# @(#)swiss.awk	1.1 3/28/97

BEGIN { state = 0 ; ce = 0 ; seq="" ; ac = "" ;}
/^ID/ {seq = "SW:"$2 ; gsub(/;/,"",seq) ; next;}
/^AC/ {ac = $2 ; gsub(/;/,"",ac) ; next;}
/^DE/ {de = substr ($0,6) ;  next;}
/^RA/ {ra = substr ($0,6) ;  gsub(/\./,"",ra) ; gsub(/, /,"\nAuthor \"",ra) gsub(/;/,"",ra) ;  next;}
/^RX   MEDLINE/ {rx = substr ($0,15) ; gsub(/\./,"",rx) ;  next;}
/^\/\// {
if (ce)
  {
    printf("\nProtein %s\n",seq);
    printf("Database SWISSPROT %s\n",ac);
    if (de != "") printf("Title \"%s\n",de);
    if (ra != "") printf("Author \"%s\n",ra);
    if (rx != "") printf("Medline_acc %s\n",rx);
    if (state>1)  printf("\n%s\n", aa) ;
  }
state = 0 ; seq="" ; ac = "" ; ra = "" ; rx = "" ;
ce = 0 ;
next; }

/^SQ/ {state = 2 ; next ;}
{ if (state == 2) { aa="Peptide " seq  "\n" ; state = 3 ;  }
  if (state == 3) 
    { pp = $0 ; gsub(/[0-9]/,"",pp) ; gsub(/ /,"",pp) ;
      aa = aa  pp  "\n" ;
    }
}

/^OS   CAENORHABDITIS ELEGANS/ { ce=1 ;}

#include "regular.h"
#include "aceio.h"
#include "mytime.h"

#define PREFIX(b,s) ( strncmp(b,s,sizeof(s)-1) == 0 )


static char *fAceConfig (const char *tag, AC_HANDLE h)
{
  char *val = 0, *cp ;

  cp = sessionFilName ("wspec/fiche", "wrm", 0) ;
  
  if (cp && tag && *tag)
    {
      ACEIN fi = aceInCreateFromFile (cp, "r", 0, h) ;

      while (aceInCard (fi))
	{
	  cp = aceInWord (fi) ;
	  if (cp && ! strcmp (cp, tag))
	    {
	      cp = aceInWord (fi) ;
	      if (cp)
		val = strnew (cp, h) ;
	      break ;
	    }
	}
      aceInDestroy (fi) ;
    }

  return val ;
}

static char *getPrefix (ACEIN fi, char *prefix)
{
  char *cp ;
  
  while (aceInCard (fi))
    {
      cp = aceInPos (fi) ;
      if (!cp || strncmp (cp, prefix, strlen(prefix)))
	continue ;
      return cp + strlen(prefix) ; /* the rest of the line */
    }
  return 0 ;
}

static int aceKog2aceParse (ACEIN fi)
{
  int nKantor = 0, subjectTitle, gene, link, product ;
  int a1, a2, x1, x2,  b1, y1, eValue, nHit ;
  float bitScore = 0, bestScore = 0 ;
  char *cp, *cq, species ;
  Stack s = stackCreate (1000) ;
  int spec[256], tag[256] ;

  pushText (s, "-") ;
  for (a1 = 0 ; a1 < 256 ; a1++)
    spec[a1] = 0 ;
  spec['w'] = stackMark (s) ;
  pushText (s, "Caenorhabditis elegans") ;
  spec['a'] = stackMark (s) ;
  pushText (s, "Arabidopsis thaliana") ;
  spec['h'] = stackMark (s) ;
  pushText (s, "Homo sapiens") ;
  spec['d'] = stackMark (s) ;
  pushText (s, "Drosophila melanogaster") ;
  tag['w'] = stackMark (s) ;
  pushText (s, "AceKogWorm") ;
  tag['a'] = stackMark (s) ;
  pushText (s, "AceKogAra") ;
  tag['h'] = stackMark (s) ;
  pushText (s, "AceKogHuman") ;
  tag['d'] = stackMark (s) ;
  pushText (s, "AceKogDroso") ;
 lao: 

  bestScore = 0 ;
  cp = getPrefix (fi, "<b>Query=</b>") ;
  if (!cp)
    {
      stackDestroy (s) ;
      return nKantor ;
    }

  if (!cp || !*cp++ || !*cp)
    goto lao ;
  nKantor++ ;
  printf ("\nKantor \"%s\"\n", cp) ;
  printf ("-D AKG\n") ;
  printf("AceKog_Date %s\n", timeShowNow ());
  printf("Kantor_Date %s\n", timeShowNow ());
  nHit = 0 ;

 nextHit:
  cp = getPrefix (fi, "><a name =") ;
  if (!cp) goto lao ;

 /* one letter species code */
  cp = strstr (cp, "</a>") ;
  if (!cp) goto lao ;
  cp += 4 ;
  species = *cp++ ; 

 /* first word of title is gene name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  gene = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

 /* next word of title is locuslink/newname name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  link = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

  /* next word of title is product name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  product = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

  /* rest if full title on several lines */
  subjectTitle = stackMark (s) ; 
  pushText (s, aceInUnprotect (fi, cp)) ;
  if (aceInCard (fi) &&
      (cp = aceInPos (fi)) &&
      (cp = strstr (cp, "SET\\/")))
    {
      cp += 4 ; *cp = ' ' ;
      cq = strstr (cp, "query") ;
      if (cq ) *cq = 0 ;
      catText (s, cp) ;
    }

  /* Score =  214 bits (545), Expect = 1e-55 */
  cp = getPrefix (fi, "Score") ;
  if (! aceInWord (fi)  || /* Score */
      ! aceInWord (fi)  || /* = */
      ! aceInFloat (fi, &bitScore)  ||  /* 214 */
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! (cp = aceInWord (fi)))
    {
      fprintf (stderr, "// bad score line %d: %s\n", aceInStreamLine (fi), cp) ;
      goto lao ;
    }
  eValue = stackMark (s) ;
  pushText (s, cp) ; /* very often out of range of a float */

  aceInCard (fi) ;  aceInCard (fi) ;  aceInCard (fi) ; 
  b1 = y1 = -1 ;
  while (1) 
    {
      /* Query: 135 IESSRDLLHRIKDEVGAPGIVVGVSVDGKEVWSEGLGYADVENRVPCKPETVMRIASISK 194 */
      if (!aceInCard (fi)) break ;
      cp = aceInWord (fi) ;
      if (!cp || strcmp (cp, "Query:"))
	break ;
      if (!aceInInt (fi, &a1) || !aceInWord (fi) || !aceInInt (fi, &a2))
	break ;
      if (b1 == -1) b1 = a1 ;
      /* I+ +++L+       G PG+ + VS+DGK VW  G GYA++E+   C  ++VMRIASISK */
      aceInCard (fi) ;  /* jump that line */
      
      /* Sbjct: 115 IKKAKELVETTMAIQGIPGLSIAVSLDGKMVWKSGFGYANLESFAKCTGDSVMRIASISK 294 */
      if (!aceInCard (fi)) break ;
      cp = aceInWord (fi) ;
      if (!cp || strcmp (cp, "Sbjct:"))
	break ;
      if (!aceInInt (fi, &x1) || !aceInWord (fi) || !aceInInt (fi, &x2))
	break ;
      if (y1 == -1) y1 = x1 ; 
      if (!aceInCard (fi)) break ;
    }

  if (y1 != -1 && 100 * bitScore > 88 * bestScore)      /* success export one hit */
    {
      printf ("AceKog %s \"%s\" %g %d %d %d %d %s\n"
	      , stackText (s, product)
	      , stackText (s, spec[(int)species])
	      , bitScore, b1, a2, y1, x2
	      , aceInProtect (fi, messprintf ("%s [%s] eValue = %s",
					      stackText (s, subjectTitle) 
					      , stackText (s, spec[(int)species])
					      , stackText (s, eValue)
					      )
			      )
	      ) ;
      if (1 || !nHit) printf ("%s %s %s\n"
			 , stackText (s, tag[(int)species])
			 , stackText (s,gene)
			 , stackText (s,link)
			 ) ;
      if (! nHit || bitScore > bestScore)
	bestScore = bitScore ;
      nHit++ ;
    }
  if (nHit >= 36)
    goto lao ;
  goto nextHit ;
} /* aceKog2aceParse */

/*
*****************************
* PSORT
*****************************
*/
static char *psortSearchDottedSequence (char * src, char *regexp, int direction, int caseIgn)
{
  int sc, i, len = strlen (src) ;
  char *rexp;
  
  if (direction == -1) 
    i = len - 1 ;
  else 
    i = 0 ;
  
  for (; i>=0 && i<len ; i += direction)
    {
      for (rexp = regexp ; *rexp ; rexp += strlen(rexp) + 1)
	{
	  for (sc=0 ; rexp[sc] && src[i+sc] ; sc++)
	    {  /* compare current position with regexp */
	      if (rexp[sc]=='.') 
		continue;
	      
	      else if (rexp[sc] != src[i+sc])
		break;
	    }
	  if (regexp[sc] == 0)
	    return src + i ; /* correspondence */
	}
    }
  return 0;
}


static char *kantorParsePSORTDomain (char *src, char *strLook, int seq, int st, int ed, int dsc, char *pseq, int * pSt, int *pEd, char *pdsc, int * cnt)
{
  int i;
  char *ptr;
  
  *cnt = 0;
  if (strLook)
    {ptr = strstr (src, strLook) ;if (!ptr)return 0;else ptr += strlen (strLook) ;}
  else ptr = src;
  
  for (i = 0;ptr && (i <= seq || i <= st || i <= ed) ;i++)
    {
      if (i==seq && pseq)
	{sscanf(ptr,"%s ",pseq) ;(*cnt)++;if (strstr(pseq,"none"))
	  {ptr=vstrSkipWords(ptr,1,0) ; break;}}
      if (i==st && pSt)
	{sscanf(ptr,"%d ",pSt) ;(*cnt)++;}
      if (i==ed && pEd)
	{sscanf(ptr,"%d ",pEd) ;(*cnt)++;}
      if (i==dsc && pdsc)
	{sscanf(ptr,"%s ",pdsc) ;(*cnt)++;}
      ptr=vstrSkipWords(ptr,1,0) ;
      
    }
  if (seq != -1 && pseq && ed == -1 && pEd && st != -1 && pSt && i >= st && i >= seq)
    *pEd = *pSt+strlen (pseq)-1;
  return ptr;
}

#define fixPsortCoord(cord)		 (cord >= xPos ? (cord+xLen) : cord)

static int kantorParsePSORTList (vSTR * blk, char *src, char *look, char *section, int iseq, int ist, char *seqBuf)
{
  char   seq[vTEXTMAXFILESIZE], *ptr, *fnd;
  int    qrySt, qryEd, iNum, cnt, xPos, xLen = 0 ;
  
  /* find XXX sequence position and length */
  if ((ptr = strchr (seqBuf, 'X')))
    {
      xPos=ptr-seqBuf;
      while(*ptr=='X')ptr++;
      xLen=ptr-seqBuf-xPos;
    }else {xPos=0;xLen=0;}
  
  /* 2nd peroxisomal targeting signal*/
  for (ptr = src, iNum = 0;ptr && *ptr;iNum++)
    {
      /* scan variables about domain section */
      ptr=kantorParsePSORTDomain(ptr,look,iseq,ist,-1,-1,seq,&qrySt,&qryEd,0,&cnt) ;if (!cnt)break;
      if (strstr(seq,"none"))continue;
      if (strstr(seq,"found"))continue;
      
      if (ist!=-1)
	{qrySt=fixPsortCoord(qrySt) ;qryEd=fixPsortCoord(qryEd) ;}
      /* find the sequence from C-terminus if position is not defined */
      if (ist==-1)
	{vstrCleanEnds(seq,seq,0," \n\t",1) ;
	if (!strstr(seq,"none"))
	  {sscanf("%s ",seq) ;
	  seq[strlen(seq)+1]=0;
	  if ((fnd=psortSearchDottedSequence(seqBuf,seq,-1,1))!=0)
	    {qrySt=fnd-seqBuf;qryEd=qrySt+strlen(seq) ;qrySt++;} 
	  else break;
	  }
        }
      if (! (xLen && qrySt < xPos && qryEd > xPos))
	  vstrPrintf(blk,"Domain \"%s_domain\" \"Psort\" %d. %d %d %d %d \"%s\" \n"
		     ,section,1,qrySt,qryEd,1,1,seq) ;
    }
  return iNum;
}

struct 
{char *st, * ed, *fnd, *nm;int seq, qst;}lkmot[] = 
{
  
  {"nuclear localization signals</A>\n\0\0","content of basic\0<A HREF=\0\0",0,"Nuclear_localization",1,3},
  {"ER retention motif in the C-terminus:\0\0","\n\0<A HREF=\0\0",0,"ER_retention",0,-1},
  {"ER Membrane Retention Signals:</A>\0\0","<A HREF=\0\0","motif in the C-terminus:","ER_membrane",0,-1},
  /* Nakai himself says it is unreliable and it hits 90% of the big human genes, so forget it, jan 2004
  {"Gavel: prediction of cleavage sites for mitochondrial preseq</A>\0\0","<A HREF=\0\0",0,"Mitochondrial_cleavage_site",4,3},
  */
  {"peroxisomal targeting signal in the C-terminus:\0\0","\n\0<A HREF=\0\0",0,"peroxisomal",0,-1},
  {"2nd peroxisomal targeting signal:  found\0\0","<A HREF=\0\0",0,"2nd_peroximal",0,2},
  {"possible vacuolar targeting motif: found\0\0","<A HREF=\0\0",0,"vacuolar",0,2},
  {"RNA-binding motif:</A> found\0\0","<A HREF=\0\0",0,"RNA_binding",0,2},
  {"Actinin-type actin-binding motif:</A>\0\0","type 2\0<A HREF=\0\0","type 1: found\n","actin_binding_1",0,2},
  {"Actinin-type actin-binding motif:</A>\0\0","<A HREF=\0\0","type 2: found\n","actin_binding_2",0,2},
  {"N-myristoylation pattern :\0\0","\n\0<A HREF=\0\0",0,"N_myristoylation",0,-1},
  {"Prenylation motif:</A>\0\0","\n\0<A HREF=\0\0","motif in the C-terminus:","prenylation",0,-1},
  {"transport motif from cell surface to Golgi: found\0\0","<A HREF=\0\0",0,"Golgi_transport",0,2},
  {"Dileucine motif in the tail:</A> found\0\0","<A HREF=\0\0",0,"Dileucine",0,2},
  {"Zinc finger, C2H2 type, domain (PS00028):  *** found ***\0\0","\n\n\0<A HREF=\0\0",0,"Zinc_finger",0,2},
  {"Leucine zipper pattern (PS00029):  *** found ***\0\0","\n\n\0<A HREF=\0\0",0,"Leucine_zipper",0,2},
    {"Ribosomal protein S4e signature (PS00528):  *** found ***\0\0","<A HREF=\0\0",0,"Ribosomal_protein",0,2},
  {"Sigma-54 interaction domain ATP-binding region A signature (PS00675):  *** found ***\0\0","<A HREF=\0\0",0,"ATP_binding",0,2},
  {0,}
};

static int psort2aceParse (char *kantorName, char *seqBuf, ACEIN fi, AC_HANDLE h)
{
  int     is, ln, qrySt, qryEd = 0, sbjSt = 1, sbjEd = 1, iNum = 0, prvS, prvE, cnt = 0, i, clen, lnBuf;
  char *anB, *sctP, *ptr, seq[vTEXTMAXFILESIZE];
  int		xPos, xLen = 0, iScore = 1;
  float fScore = 0 ;
  char  pep = 0;
  char *cp ;

  if (!psortBuf)return 0;
  ln = strlen (psortBuf) ;if (!ln)return 0;
  anB = messalloc (ln+1) ;if (!anB)return 0;
  sctP = psortBuf;
  
  /* find XXX sequence position and length */
  lnBuf = strlen (seqBuf) ;
  xLen = 0 ; xPos = 0  ;
  cp = seqBuf ;
  while (*cp && *cp != 'X')
    { cp++ ; xPos++ ; }
  if (*cp == 'X')
    {
      while (*cp == 'X')
	xLen++ ; cp++ ;
    }
  else
    xPos = 0 ;
  
  /* cleaveable singal peptide */
  sctP=vstrClnMarkup(psortBuf,anB,0,">>> Seems to have a cleavable signal peptide (\0\0",")\0\n\0\0","",1,1,0) ;
  if (anB[0])
    {
        sscanf(anB,"%d to %d",&qrySt,&qryEd) ;
	qrySt=fixPsortCoord(qrySt) ;qryEd=fixPsortCoord(qryEd) ;
	ln=qryEd-qrySt+1;
	if (qryEd>lnBuf)strcpy(seq,"recommended to reKantorize") ;
	else {memcpy(seq,seqBuf+qrySt-1,ln) ;seq[ln]=0;}
	if (! (xLen && qrySt < xPos && qryEd > xPos))
	    vstrPrintf(blk,"Domain \"%s\" \"Psort\" %d. %d %d %d %d \"%s\" \n","Cleavable_signal_peptide",iScore,qrySt,qryEd,sbjSt,sbjEd,seq) ;
        iNum++;
    }
  
  /* transmembrane domains */
  sctP=vstrClnMarkup(psortBuf,anB,0,"Tentative number of TMS(s)\0\0","ALOM score\0<A HREF=\0\0","",1,1,0) ;
  vstrFindReplaceSymbols(anB,anB,0,"-"," ",0,1,0) ;
    for(prvS=-1,prvE=-1,ptr=anB;ptr && *ptr;iNum++)
      {
        ptr=kantorParsePSORTDomain(ptr+1,"Transmembrane",-1,0,1,-1,0,&qrySt,&qryEd,0,&cnt) ;
	if (!ptr)break;
	qrySt=fixPsortCoord(qrySt) ;qryEd=fixPsortCoord(qryEd) ;
	if (ptr && prvE>=qrySt)qrySt=prvS; /* glue if overlap */
        else if (prvS!=-1 && prvE!=-1)
	  {
	    if (prvE>prvS && prvE<prvS+sizeof(seq))ln=prvE-prvS+1; 
	    else ln=0;
	    if (prvE>lnBuf)strcpy(seq,"recommended to reKantorize") ;
	    else { memcpy(seq,seqBuf+prvS-1,ln) ;seq[ln]=0;}
	    if (! (xLen && qrySt < xPos && qryEd > xPos))
		vstrPrintf(blk,"Domain \"%s_domain\" \"Psort\" %d. %d %d %d %d \"%s\" \n","Transmembrane",iScore,prvS,prvE,sbjSt,sbjEd,seq) ;
	  }
        if (!cnt)break;
        prvS=qrySt;prvE=qryEd;
      }
    
    /* other domains */
    for(cnt=0,i=0;lkmot[i].st!=0;i++,iNum++)
      {
        sctP=vstrClnMarkup(psortBuf,anB,0,lkmot[i].st,lkmot[i].ed,"",1,1,0) ;
        vstrClnMarkup(anB,anB,0,"(\0\0",")\0\0","",0,0,0) ;
        if (anB[0])cnt+=kantorParsePSORTList(blk,anB,lkmot[i].fnd,lkmot[i].nm,lkmot[i].seq,lkmot[i].qst,seqBuf) ;
      }
    /* coil coiled regions */
    sctP=vstrClnMarkup(psortBuf,anB,0,"detect coiled-coil regions\0\0","total\0\0","",1,1,0) ;
    vstrClnMarkup(anB,anB,0,"</A>\0\0","\0\0","",1,1,0) ;
    clnVar(anB) ;
    if (anB[0])for(clen=0,is=0,ln=2,ptr=anB,qrySt=-2,i=0;ln==2;i++,is++)
      {
        
        if (ptr)ln=sscanf(ptr,"%d %c",&qryEd,(char *)&pep) ;else ln=0;
        if (!i)
	  {qrySt=qryEd;prvE=qryEd;}
        
        if ((qrySt>=0 && qryEd-prvE>1) || ln!=2)
	  {seq[is]=0;
	  if (! (xLen && qrySt < xPos && qryEd > xPos))
	    vstrPrintf(blk,"Domain \"%s\" \"Psort\" %d. %d %d %d %d \"%s\" \n","Coiled_coil_region",iScore,fixPsortCoord(qrySt),fixPsortCoord(prvE),sbjSt,sbjEd,seq) ;
	    qrySt=qryEd;cnt++;clen+=is;is=0;
	  }
        if (ln!=2)break;
        ptr=vstrSkipWords(ptr,3,0) ;
        prvE=qryEd;seq[is]=pep;
      }
    
    /* localization */
    sctP=vstrClnMarkup(psortBuf,anB,0,"-NN Prediction\0\0",">>\0\0","",1,1,0) ;
    vstrClnMarkup(anB,anB,0,"<PRE>\0\0","\0","",1,1,0) ;
    clnVar(anB) ;
    if (anB[0])for(ptr=anB,i=0;ptr && *ptr;i++)
      {
        if (sscanf(ptr,"%f",&fScore)==0)break;
        ptr=vstrClnMarkup(ptr,seq,0,"%:","\n\0\0","",1,1,0) ;
	clnVar(seq) ;
		vstrFindReplaceStrings(seq,seq,0,"extracellular, including cell wall"_00, "Secreted"_00,0,1) ;
		vstrFindReplaceSymbols(seq,seq,0," ","_",0,1,0) ;
		if (seq[0])vstrPrintf(blk,"Psort Localization \"%s\" %4.1f \n",seq,fScore) ;
      }
    messfree (anB) ;
    return i;
}

/******************************************************************/

static int psort2aceParse (char *kantorName, char *seqBuf, ACEIN fi, AC_HANDLE h)
{
  vSTR * blk=(vSTR * )param;
  
  if(!strstr(content,">> prediction for")){fprintf(stderr,"BROKEN Psort output for query %s\n",nm);return 1;}
  
  printf ("Kantor %s\n", kantorName);
  printf ("-D PSR\n") ;
  psort2aceDoParse (kantorName, seqBuf, fi, h) ;
  vstrPrintf(blk,"Psort_Date %s\n", timeShowNow()) ;
  vstrPrintf(blk,"Kantor_Date %s\n\n", timeShowNow()) ;
  
  return 1;
}

/******************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage psort2ace < x.psort >! x.psort.ace \n") ;
  exit (1) ;
} /* usage */

/******************************************************************/
/* split the stdin file on Query= */
static int aceKog2aceSplit (ACEIN fi)
{
  Stack s = stackCreate (10000) ;
  ACEIN splitFi = 0 ;
  int nn = 0, line = 0, state = 0 ;
  char *cp, *prefix, *kantorName, *seqBuf ; 
  AC_HANDLE h = 0 ;
 
  stackTextOnly (s) ;
  while (TRUE) 
    {
      if (aceInCard (fi))
	cp = aceInPos (fi) ;
      else
	cp = 0 ;
      line++ ;
      switch (state)
	{
	case 0: /* outside */
	  prefix = "<H1>Input Sequence</H1>" ;
	  if (!cp || !strncmp (cp, prefix, strlen(prefix)))
	    continue ;
	  state = 1 ;
	  aceInCard (fi) ; line++ ; /* jump <PRE> line */
	  break ;
	case 4: /* inside result part */
	  prefix = "<H1>Input Sequence</H1>" ;
	  if (!cp || !strncmp (cp, prefix, strlen(prefix)))
	    {
	      if (cp) 
		catText (s, cp) ; 
	      catText (s, "\n") ;
	      continue ;
	    }
	  if (stackMark (s))
	    {
	      splitFi = aceInCreateFromText (stackText (s, 0), 0, 0) ;
	      nn += psort2aceParse (kantorName, seqBuf, splitFi, h) ;
	      aceInDestroy (splitFi) ;
	    }
	  stackClear (s) ;
	  messfree (h) ;
	  state = 1 ;
	  aceInCard (fi) ; line++ ; /* jump <PRE> line */
	  break ;
	case 1: /* wainting for kantor name */
	  if (!cp)
	    {
	      fprintf (stderr, "missing kantor name (i.e. N1359571160 (79 aa)) at line %d\n", line) ;
	      state = 0 ;
	    }
	  h = handleCreate () ;
	  kantorName = strnew (cp, h) ;
	  state = 2 ;
	  break ;
	case 2:  /* waiting for actual peptide sequence */
	  if (!cp) continue ;
	  seqBuf = strnew (cp, h) ;
	  state = 3 ;
	  break ;
	case 3:  /* waiting for start of results */
	  prefix = "<H1>Results of Subprograms</H1>" ;
	  if (!cp || !strncmp (cp, prefix, strlen(prefix)))
	    continue ;
	  state = 4 ;
	  stackClear (s) ;
	  break ; 
	}
    }

  messfree (h) ;
  stackDestroy (s) ;
  return nn ;
}

/******************************************************************/

int main (int argc, char *argv[])
{
  if (argc == 1)
    {
      ACEIN fi = aceInCreateFromStdin (0, 0, 0) ;

      int nn = aceKog2aceSplit (fi) ;
      fprintf (stderr, "// Parse %d blast outputs\n", nn) ;
      aceInDestroy (fi) ;
    }
  else
    usage () ;
    
  return 0 ;
}

/******************************************************************/
/******************************************************************/

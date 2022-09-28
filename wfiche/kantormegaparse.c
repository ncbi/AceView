/*

using the linux cflow, 

cflow -Dcflow -r main -Iwh -a -d 10 wfiche/kantormegaparse.c wfiche/facelib.c > file

*/


#include "ac.h"
#include    <vstd.h>
#include "../wfiche/taxtree.h"
#include "vtxt_.h"
#include "mytime.h"
#include "bitset.h"

static char *selfSpecies = 0 ;

static BOOL DELETE_OLD_DATA = TRUE ; /* is TRUE -D previous results */

#define clnVar(v_var)   vstrCleanEnds((v_var),(v_var),0," \t\n\r",1)

typedef struct  
{
  int     num;
  char *dscr;
  int qrySt, qryEd, sbjSt, sbjEd;
  double iEvalue, iScore;
  int	Identities, Pozitives, Gaps;
  char species[1024];
} pepHIT;

struct stat stb;
/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  PARSING FUNCTIONS                       _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/*
*****************************
* BLAST
*****************************
*/
static int kantorParseBlastOrder (const void *va, const void *vb)
{
  const char *a = *(const char **)va, *b = *(const char **) vb ; 
  char *zz[] = {"ref|NP","sp","gb", "emb", "pir","ref|XP", ""} ;
  int ia, ib ;

  ia = ib = 999 ;
  if (a)
    for (ia = 0 ; *zz[ia]; ia++)
      if (strstr (a, zz[ia])) break ;
  if (b)
    for (ib = 0 ; *zz[ib]; ib++)
      if (strstr (b, zz[ib])) break ;

  return ia - ib ;
}

static char * kantorParseBestName (char *nam, char *untreat)
{
  Array aa = arrayCreate (12, char *) ;
  int n = 0 ;
  char *cp, *cq, *cr, buf[vTEXTMAXFILESIZE] ;
  BOOL debug = FALSE ;
  
  strcpy (buf, nam) ;
  strcat (buf, untreat) ;
  cp = buf ;
  while (cp && *cp)
    {
      cq = strstr (cp, "<>") ;
      if (cp == cq) { cp+= 2 ; continue ; }
      cr = strstr (cp, "|") ;
      if (cr && cq && cr < cq)
	{ *cq = ' ' ; *(cq+1) = ' ' ;  cq = strstr (cp, "<>") ; }
      if (cq) *cq = 0 ;
      if (*cp)
	array (aa, n++, char *) = cp ;
      if (!strncmp (cp, "sp", 2))
	{ while (*cp) { *cp = ace_lower(*cp) ; cp++ ; } }
      cp = cq ? cq+2 : "" ;
    }
  if (debug)
    for (n = 0 ; n < arrayMax (aa) ; n++)
      {
	cp = arr (aa, n, char *) ;
	printf ("VV %s\n", cp) ;
      }
  arraySort (aa, kantorParseBlastOrder) ;
  if (debug)
    for (n = 0 ; n < arrayMax (aa) ; n++)
      {
	cp = arr (aa, n, char *) ;
	printf ("VV2 %s\n", cp) ;
      }
  if (n)
    strcpy (untreat, arr (aa, 0, char*)) ;
  else
    untreat [0] = 0 ;

  arrayDestroy (aa) ;
  
  return cp ;
}

/* analyses the BLAST results and prepares the Resulting file for FICHE */
static int kantorParseBLAST (vSTR * blkp, vSTR * hintBlkp, char *blastBuf, char* separ)
{
  char    *cp, *cq, myname[vTEXTMAXFILESIZE], info[vTEXTMAXFILESIZE], bfr[2*vTEXTMAXFILESIZE], untreat[vTEXTMAXFILESIZE];
  int     ln, lets, iNum = 0, icnt, ih, isSkip, segmentedHit, iIndepHit, absErr, isSameSpecies;
  char *anB, *ptr, *tblPtr, *dscPtr, *nxtDsc ;
  vMEX    ht;
  pepHIT  ph, *pp;
  int isPir, isEmb, isGb ;

  if (!blastBuf)return 0;
  if (strstr (blastBuf, "No hits found"))return -1;
  ln = strlen (blastBuf) ;if (!ln)return 0;
  anB = messalloc (ln+1) ;if (!anB)return 0;
  vSet0 (ht) ;
  
  /* get the total number of letters */
  ptr=vstrSearchSubstring(blastBuf,"Query=\0\0",1,"",0) ;if (!ptr)goto final_kb;
  vstrClnMarkup(ptr,anB,0,"(\0\0","letters\0\0","",1,1,0) ;
  lets=0;sscanf(anB,"%i ",&lets) ;if (!lets)goto final_kb;
  
  /* change all html markup into <> and position to the table */
  vstrClnMarkup(blastBuf,anB,0,"<\0\0",">\0\0","<>",0,0,0) ;
  tblPtr=vstrSearchSubstring(anB,"Sequences producing significant\0\0",1,"",0) ;if (!tblPtr)goto final_kb;
  tblPtr=vstrSearchSubstring(tblPtr,"<>\0\0",1,"",0) ;if (!tblPtr)goto final_kb;
  
  dscPtr = tblPtr+3;
  
  for (iIndepHit = 0, iNum = 0;dscPtr && tblPtr;iNum++)
    {
      memset (&ph, 0, sizeof (ph)) ;
      tblPtr += 2;

      /* find the name */
      tblPtr=vstrSeparateSubstringFromString(tblPtr,myname,0,0,"<>\0\0",0) ;if (!tblPtr)break;
      vstrCleanEnds(myname,myname,0,"\r\n \t",1) ;if (!myname[0])break;
      

      /* find score and evalue */
      tblPtr=vstrSeparateSubstringFromString(tblPtr+1,bfr,0,1,"<>\0\0",0) ;if (!tblPtr)break;
      sscanf(bfr,"%lf",&ph.iScore) ;
      vstrSeparateSubstringFromString(tblPtr+1,bfr,0,0,"<>\0\n\0\0",0) ;
      sscanf(bfr,"%lf",&ph.iEvalue) ;
      
      /* look for information block by the name */
      dscPtr=vstrSearchSubstring(dscPtr,myname,1,"",0) ;if (!dscPtr)break;
      /* for analyzer */
      
      {
	vstrClnMarkup(dscPtr,untreat,0,myname,"Length =\0\0","",1,1,0) ;
	if (0) strcpy (untreat,"aaaaa<hello world>bbbbb<junk>cc  \n  cc</a>dddd") ;
	vstrFindReplaceStrings(untreat,untreat,0,"</a>\0\0","",0,0) ;
	if (0) vstrClnMarkup(untreat,untreat,0,"<\0\0",">\0\0","<>",0,0,0) ;
      /* here we replace any consecutive string of 'spaces' by a single blank */
	vstrFindReplaceSymbols(untreat,untreat,0,"\n \t"," \0\0",0,1,1) ;

      /* origin of homol
       * here we have the full description of all the entries
       * considered equivalent by BLAST
       */
      isGb =  strstr(untreat,"gb|")   ? TRUE : FALSE ;
      isEmb = strstr(untreat,"emb|")  ? TRUE : FALSE ;
      isPir = strstr(untreat,"pir|")  ? TRUE : FALSE ;
      /* isRef = strstr(untreat,"ref|")  ? TRUE : FALSE ; */

      /* species */ 
      ptr = untreat;
      while (ptr && *ptr)
	
	{
	  ptr=vstrClnMarkup(ptr,ph.species,0,"["_00,"]\0Length ="_00,"",1,1,0) ; 
	  vstrFindReplaceSymbols(ph.species,ph.species,0,"\r\n \t"," ",0,1,1) ;
	  vstrCleanEnds(ph.species,ph.species,0,"\n \t",1) ;
	  if (	!strstr(ph.species,"imported") && !strstr(ph.species,"vector"))
	    break;
	  else ph.species[0] = 0;
	}
      if (0) printf ("SPECIES = %s\n", ph.species) ;

      kantorParseBestName (myname, untreat) ;
      }
      /* prepare the name and the information string */
      cp = untreat ;
      while (*cp && (cq = strstr (cp, "|"))) cp = cq + 1 ; 
      while (*cp == ' ') cp++ ;
      strcpy (info, cp ? cp : "") ;

      cp = untreat ;
      if ((cq = strstr (cp, "|")) &&
	  (cp = cq + 1) &&
	  strstr (cp, "|"))
	{
	  while (*cp == ' ') cp++ ;
	  strcpy (myname, cp) ;
	  if ((cq = strstr (myname, "|")))
	    *cq = 0 ;
	}

      /* analyse the alignment paragraphs */
      segmentedHit = 0;
      /* find first and the last aligned aminoacid 	query lines */
      dscPtr=strstr(dscPtr,"Length =") ;if (!dscPtr)break;
      nxtDsc=strstr(dscPtr+10,"><>") ;if (nxtDsc)*nxtDsc=0;
      *dscPtr=0;dscPtr+=1;

      for(ptr=dscPtr;ptr;)
	{
	  
	  ptr=vstrSearchSubstring(dscPtr,"Score ="_00,1,"<>"_00,0) ;if (!ptr)break;ptr+=7;
	  dscPtr=ptr;
	  sscanf(ptr,"%lf",&ph.iScore) ;
	  ptr=vstrSearchSubstring(dscPtr,"Identities ="_00,1,"<>"_00,0) ;if (!ptr)break;ptr+=12;
	  sscanf(ptr,"%i",&ph.Identities) ;
	  ptr=vstrSearchSubstring(dscPtr,"Positives ="_00,1,"<>"_00,0) ;if (!ptr)break;ptr+=11;
	  sscanf(ptr,"%i",&ph.Pozitives) ;
	  ptr=vstrSearchSubstring(dscPtr,"Gaps ="_00,1,"<>"_00,0) ;
	  if (ptr)
	    {ptr+=6;sscanf(ptr,"%i",&ph.Gaps) ;}
	  
	  ptr=vstrSearchSubstring(dscPtr,"Query:"_00,1,"<>"_00,0) ;if (!ptr)break;ptr+=6;
	  sscanf(ptr,"%i ",&ph.qrySt) ;
	  ptr=vstrSearchSubstring(dscPtr,"Query:"_00,-1,"<>"_00,0) ;if (!ptr)break;ptr+=6;
	  sscanf(ptr,"%s %s %i ",bfr,bfr,&ph.qryEd) ;
	  ptr=vstrSearchSubstring(dscPtr,"Sbjct:"_00,1,"<>"_00,0) ;if (!ptr)break;ptr+=6;
	  sscanf(ptr,"%i ",&ph.sbjSt) ;
	  ptr=vstrSearchSubstring(dscPtr,"Sbjct:"_00,-1,"<>"_00,0) ;if (!ptr)break;ptr+=6;
	  sscanf(ptr,"%s %s %i ",bfr,bfr,&ph.sbjEd) ; 
	  
	  /* position to the next section of the same hit */
	  ptr=vstrSearchSubstring(dscPtr,"Score ="_00,1,"Length ="_00,0) ;
	  if (ptr)dscPtr=ptr;
	  
	  /* skipping sequences */
	  isSkip = 0;
				/* define the same species */
	  isSameSpecies = 0;
	  if (!strcasecmp (selfSpecies, "none"))
	    isSameSpecies = 1;					
	  if (!strcasecmp (selfSpecies, "worm") && !strcasecmp (ph.species, "Caenorhabditis elegans"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "human") && !strcasecmp (ph.species, "Homo sapiens"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "hs") && !strcasecmp (ph.species, "Homo sapiens"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "ara") && !strcasecmp (ph.species, "Arabidopsis thaliana"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "droso") && !strcasecmp (ph.species, "Drosophila melanogaster"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "mouse") && !strncasecmp (ph.species, "Mus musculus", 12)) /* protect against mm domesticus */
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "mm") && !strncasecmp (ph.species, "Mus musculus", 12))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "rat") && !strcasecmp (ph.species, "Rattus norvegicus"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "rn") && !strcasecmp (ph.species, "Rattus norvegicus"))
	    isSameSpecies = 1;
	  else if (!strcasecmp (selfSpecies, "banana") && !strcasecmp (ph.species, "Musa acuminata"))
	    isSameSpecies = 1;

	  if (0 && strstr (myname, "XP_" ))
	    isSkip = 1;

	  if (isSameSpecies && 
	      !segmentedHit)	    
	    { /* for segmented hits we do not skip continuations */
	      pp=(pepHIT *)ht.buf; /* count hits overlapping with this */
	      for(icnt=0,ih=0;ih<(ht.cnt) ;ih++)
		{
		  if (pp[ih].qrySt <= ph.qrySt && pp[ih].qryEd >= ph.qryEd)
		    icnt++;
		  if ( abs((pp[ih].Identities-ph.Identities))<2 ||
		       abs((pp[ih].Pozitives-ph.Pozitives))<2 ) 
		    {
		      isSkip = 1;
		      break;
		    }
		}	      
	      
	      /* hit to self if positives and identities are the same - unprobable for different genes */
	      absErr = (ph.Identities-ph.Pozitives)*100/ ((ph.Identities+ph.Pozitives)/2) ;
	      
	      /*	      if (!ph.species[0])		isSkip = 1; */
	      if (!strcasecmp (selfSpecies, "worm") && /* we keep the humyyan, they have nice titles */
		  strstr (myname, "NP_")) /* we do align the NP internally anyway */
		{
		  if ( ((vAbs (absErr)) < 2 && isSameSpecies))
		    isSkip = 1;
		}
	      
	      if (!strcasecmp (selfSpecies, "worm") && 
		  isSameSpecies &&
		  isPir && !isEmb && !isGb)
		isSkip = 1;
	    }
	  
	  /* if there are no other hits with higher score in this area add this hit*/
	  if ( !isSkip )	    
	    {
	      vmexAppend (&ht, &ph, sizeof (ph)) ;
	      
	      /*				vstrPrintf (blkp, "%s%d %d  %s  ", segmentedHit ? "\t" : "" , ph.Identities, ph.Pozitives, ph.species) ;*/
	      sprintf (bfr, "%s eValue=%.2g", info, ph.iEvalue) ;
	      if (iIndepHit<30 && *myname)
		vstrPrintf (blkp, "BlastP \"%s\" \"BlastP\" %.2f %d %d %d %d %s \n"
			    , myname, ph.iScore, ph.qrySt, ph.qryEd, ph.sbjSt, ph.sbjEd, 
			    segmentedHit ? "" : freeprotect (bfr)) ;
	      if (!segmentedHit)		
		{
		  iIndepHit++;
		  vstrPrintf (hintBlkp, "%.2f %s %s\n", ph.iScore, untreat, separ ? separ : "") ;
		  if (0) printf ("UNTREAT: %s", vstrPtr(hintBlkp)) ;
		}
	      
	      segmentedHit = 1; /* we are inside of the segmented hit */
	    }
	  else
	    ptr = 0;
	  
	}
      if (nxtDsc)*nxtDsc='L';
      tblPtr=vstrSearchSubstring(tblPtr+1,"\n"_00,1,"",0) ;
      if (tblPtr && *tblPtr=='\n')tblPtr++;
    }
  
 final_kb:
  messfree (anB) ;
  vmexEmpty (&ht) ;
   if (0) printf ("FINISHUNTREAT: %s", vstrPtr(hintBlkp)) ;
  vstrPrintf (hintBlkp, "\n \n") ;
  
  return iNum;
}													  

/*
*****************************
* TAXBLAST
*****************************
*/

typedef struct 
{
  char * name;
  int		hits, clnHits;
  int 	level;
  char *	anc;
}CRITT;

static int readTaxTree (char *anB, CRITT * critter, int * pMaxLev)
{
  int     iNum, maxLev;
  char    *tblPtr, *ptrSpc;
  
  /* read information about critters */
  for (maxLev = 0, iNum = 0, tblPtr = anB; iNum < 1022 && tblPtr && *tblPtr;)
    {
      
      /* align to the current level */
      for(critter[iNum].level=0;*tblPtr=='.' || *tblPtr==' ';critter[iNum].level++,tblPtr++) ;
      critter[iNum].level/=2;
      if (maxLev<critter[iNum].level)maxLev=critter[iNum].level;
      
      /* the critter name */
      critter[iNum].name = tblPtr;
      for (;*tblPtr && ! (*tblPtr == '.' && (tblPtr[1] == '.' || (tblPtr[1] == ' ' && tblPtr[2] == ' '))) ;tblPtr++)
	{} ;
      *tblPtr = 0;
      clnVar (critter[iNum].name) ;
      /*vstrFindReplaceSymbols (critter[iNum].name, critter[iNum].name, 0, ".", "", 0, 1, 0) ; */
      
      /* get the hits */
      for(tblPtr++; *tblPtr=='.' || *tblPtr==' ';tblPtr++) ;
      sscanf(tblPtr,"%d",&critter[iNum].hits) ;
      critter[iNum].clnHits=critter[iNum].hits;
      
      for(ptrSpc=tblPtr+1;*ptrSpc!='[' && *ptrSpc!='\n';ptrSpc++) ;
      if (*ptrSpc=='[')
	{
	  critter[iNum].anc=ptrSpc+1;
	  for(;*ptrSpc!='\n' && *ptrSpc!=']' ;ptrSpc++) ;
	  if (*ptrSpc==']')
	    {*ptrSpc=0;tblPtr=ptrSpc+1;}
	}
      else critter[iNum].anc = 0;
      
      iNum++;
      if (*tblPtr == '\n' || *tblPtr == '\r')tblPtr++;
      else tblPtr = vstrSkipWords (tblPtr, 1, "\r\n") ;
    }
  
  if (pMaxLev)*pMaxLev = maxLev;
  
  return iNum;
}

static int critterSort (void * adr, int i, int j) 
{
  CRITT * crt = (CRITT * )adr;
  return crt[i].level>crt[j].level ? -1 : (crt[i].level<crt[j].level ? 1 : 0) ;
}



/* remove Http staff */
static char * vonlineCleanHTML(char * srcHtml,char * dstHtml,int sizeDest)
{
    return vstrClnMarkup(srcHtml,dstHtml,sizeDest,"<"_00,">"_00," ",0,0,0);
}

static int kantorParseTAXBLASTTREE (vSTR * blkp, char *taxblastBuf)
{
  int     i, j, cur, ln, iNum, myNum, maxLev, myMaxLev, isAny, par;
  CRITT	critter[1024], myCritter[256], *crt;
  int		crtInd[1024], crtStk[1024];
  char    * anB, *myBuf, *smyBuf, tmpBuf[1024], * ptr, protBuf[1024], *ancPtr;
  
  /* clean the text */
  if (!taxblastBuf || !(ln=strlen(taxblastBuf)) || !(anB=messalloc(ln+1)) )return 0;
  vstrClnMarkup(taxblastBuf,anB,0,"NAME=taxonomy_report>Taxonomy Report</A></B>"_00,"</FONT>"_00,"",1,1,0) ;
  vonlineCleanHTML(anB,anB,0) ;clnVar(anB) ;
  if ( !anB[0] || !myTaxTree || !*myTaxTree || !(myBuf=strnew(myTaxTree, 0)))return 0;
  

  /* read information about critters in this buffer */
  iNum = readTaxTree (anB, critter, &maxLev) ;
  myNum = readTaxTree (myBuf, myCritter, &myMaxLev) ;
  
  vsortQuickCallbackIndex (0, critter, iNum, crtStk, crtInd, (vSORTFunction)critterSort) ;
  
  /* remove hits for critters with dotters */
  /*	for (i = 0;i<iNum-1;i++)if (critter[i+1].level>critter[i].level)critter[i].hits = 0; */
  
  for (i = 0;i<iNum && iNum < 1023;i++)
    {
      cur = crtInd[i];
      crt = &critter[cur];
      for (par = cur-1;par>0 && critter[par].level >= crt->level;par--) ; /* get the parent in critter[cur] */
      isAny = 0;
      
		/* if this critter is found in myCritter list it keeps its hits, otherwise find the parent to keep its hits */
      for (j = 0;j<myNum;j++)
	
	{
	  if (!strcasecmp (myCritter[j].name, crt->name))
	    
	    {
	      isAny++;
	      if (par >= 0)critter[par].clnHits -= crt->hits;
	      if (myCritter[j].hits == -1)
		myCritter[j].hits = -1-crt->clnHits;
	      else 
		myCritter[j].hits = crt->clnHits;
	      break;
	    }
	}
      
      if (!isAny && crt->anc)
	
	{/* couldn't find the critter name in our tree look in ancesstrity line */
	  for (ancPtr = crt->anc+strlen (crt->anc) ;ancPtr != crt->anc;)
	    
	    { /* scan ancesstors in the reverse order */
	      if (! (ancPtr = strrchr (crt->anc, ';')))ancPtr = crt->anc;else 
		{*ancPtr = 0;ancPtr += 2;}
	      
	      for (j = 0;j<myNum;j++)
		
		{
		  if (strstr (ancPtr, myCritter[j].name))
		    
		    {
		      isAny++;
		      if (par >= 0)critter[par].clnHits -= crt->hits;
		      if (myCritter[j].hits == -1)
			myCritter[j].hits = -1-crt->clnHits;
		      else 
			myCritter[j].hits = crt->clnHits;
		      ancPtr = crt->anc;
		      break;
		    }
		}
	    }
	}
      
      if (j == myNum)
	
	{
	  if (par >= 0)
	    critter[par].clnHits -= (crt->hits-crt->clnHits) ;
	}
    }
  
  for (i = 0;i<myNum  && i < 1023;i++)
    {
      crt = &myCritter[i];
      j = crt->hits;
      
      if (j==0 || j==-1)continue;
      if (j<-1)j=-(crt->hits+1) ;
      /*		for(j=0;j<crt->level;j++)printf("..") ;*/
      vstrFindReplaceSymbols(crt->name,tmpBuf,0," /"_00,"_"_00,0,1,1) ;
      vstrFindReplaceSymbols(tmpBuf,tmpBuf,0,".","",0,1,0) ;
      vstrPrintf(blkp,"%s\n",tmpBuf) ; /* the name on the tree */
      
      smyBuf = 0;iNum = 0;
      if (crt->hits<0)
	
	{ /* take the best hit for the current creature  from the next line of the found species specific hit paragraph */
	  smyBuf=protBuf+sprintf(protBuf,">%s</a></b>",crt->name) ;smyBuf[1]=0; /* double zero */
	  if ((ptr=vstrSearchSubstring(taxblastBuf,protBuf,1,_00,0)) && (ptr=vstrSkipWords(ptr,1,"\n")) )
	    {
	      vstrSeparateSubstringFromString(ptr+1,protBuf,0,0,"\n"_00,0) ;
	      vstrClnMarkup(protBuf,protBuf,0,"<"_00,">"_00,0,0,0,0) ;
	      if ((smyBuf=vstrNextList(protBuf,1)) && (ptr=vstrNextList(smyBuf,2)) )
		{
		  sscanf (ptr, "%d", &iNum) ;
		  clnVar (smyBuf) ;
		  freeprotect (smyBuf) ;
		}
	    }
	}
      if (smyBuf && iNum)vstrPrintf(blkp,"N%s %d BlastP %s %d\n",tmpBuf,j,smyBuf,iNum) ;
      else vstrPrintf(blkp,"N%s %d\n",tmpBuf,j) ; 
    }
  
  /* free the buffers */
  messfree (anB) ;
  messfree (myBuf) ;
  
  return iNum;
}
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

static int kantorParsePSORTList (vSTR * blkp, char *src, char *look, char *domain, int iseq, int ist, char *seqBuf)
{
  char   seq[vTEXTMAXFILESIZE], *ptr, *fnd ;
  int    qrySt, qryEd = 0, iNum, cnt, xPos, xLen = 0 ;
  
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
	{
	  vstrCleanEnds(seq,seq,0," \n\t",1) ;
	  if (!strstr(seq,"none"))
	    {
	      seq[strlen(seq)+1]=0;
	      if ((fnd=psortSearchDottedSequence(seqBuf,seq,-1,1))!=0)
		{qrySt=fnd-seqBuf;qryEd=qrySt+strlen(seq) ;qrySt++;} 
	      else break;
	    }
        }
      if (! (xLen && qrySt < xPos && qryEd > xPos))
	  vstrPrintf(blkp,"Domain \"%s_domain\" \"Psort\" %d. %d %d %d %d \"%s\" \n"
		     ,domain,1,qrySt,qryEd,1,1,seq) ;
    }
  return iNum;
}

struct 
{char *st, * ed, *fnd, *nm;int seq, qst;}lkmot[] = 
{
  
  /* Nakai himself says it is unreliable and it hits 90% of the big human genes, so forget it, jan 2004
  {"Gavel: prediction of cleavage sites for mitochondrial preseq</A>\0\0","<A HREF=\0\0",0,"Mitochondrial_cleavage_site",4,3},
 
  2011_-2_15, dan ne veut plus chercher les motifs suivants which are in PFAM or are useless like LL dileucine
  {"Zinc finger, C2H2 type, domain (PS00028):  *** found ***\0\0","\n\n\0<A HREF=\0\0",0,"Zinc_finger",0,2},
  {"Dileucine motif in the tail:</A> found\0\0","<A HREF=\0\0",0,"Dileucine",0,2},
  {"Ribosomal protein S4e signature (PS00528):  *** found ***\0\0","<A HREF=\0\0",0,"Ribosomal_protein",0,2},
  {"Leucine zipper pattern (PS00029):  *** found ***\0\0","\n\n\0<A HREF=\0\0",0,"Leucine_zipper",0,2},
  {"Sigma-54 interaction domain ATP-binding region A signature (PS00675):  *** found ***\0\0","<A HREF=\0\0",0,"ATP_binding",0,2},
  {"Actinin-type actin-binding motif:</A>\0\0","type 2\0<A HREF=\0\0","type 1: found\n","actin_binding_1",0,2},
  {"Actinin-type actin-binding motif:</A>\0\0","<A HREF=\0\0","type 2: found\n","actin_binding_2",0,2},
  {"RNA-binding motif:</A> found\0\0","<A HREF=\0\0",0,"RNA_binding",0,2},
  {"ER Membrane Retention Signals:</A>\0\0","<A HREF=\0\0","motif in the N-terminus:","ER_membrane",0,-1}, // was mistakenly written C-terminus
  */

  {"ER retention motif in the C-terminus:\0\0","\n\0<A HREF=\0\0",0,"ER_retention",0,-1},


  {"peroxisomal targeting signal in the C-terminus:\0\0","\n\0<A HREF=\0\0",0,"peroxisomal",0,-1},
  {"2nd peroxisomal targeting signal:  found\0\0","<A HREF=\0\0",0,"2nd_peroximal",0,2},
  {"possible vacuolar targeting motif: found\0\0","<A HREF=\0\0",0,"vacuolar",0,2},

  {"N-myristoylation pattern :\0\0","\n\0<A HREF=\0\0",0,"N_myristoylation",0,-1},
  {"Prenylation motif:</A>\0\0","\n\0<A HREF=\0\0","motif in the C-terminus:","prenylation",0,-1},
  {"transport motif from cell surface to Golgi: found\0\0","<A HREF=\0\0",0,"Golgi_transport",0,2},


  {0,}
};

static int kantorParsePSORT (vSTR * blkp, char *psortBuf, char *seqBuf)
{
  int     is, ln, qrySt, qryEd = 0, sbjSt = 1, sbjEd = 1, iNum = 0,  prvE = 0, cnt = 0, i, clen, lnBuf;
  char *anB, *ptr, seq[vTEXTMAXFILESIZE];
  int		xPos, xLen = 0, iScore = 1;
  float fScore = 0 ;
  char  pep = 0;
  
  if (!psortBuf)return 0;
  ln = strlen (psortBuf) ;if (!ln)return 0;
  anB = messalloc (ln+1) ;if (!anB)return 0;
  lnBuf = strlen (seqBuf) ;
  
  /* find XXX sequence position and length */
  if ((ptr = strchr (seqBuf, 'X')))
    {
      xPos = ptr-seqBuf;
      while (*ptr == 'X')ptr++;
      xLen = ptr-seqBuf-xPos;
    }else 
      {xLen = 0;xPos = 0;}
  
  /* cleaveable signal peptide */
  vstrClnMarkup(psortBuf,anB,0,">>> Seems to have a cleavable signal peptide (\0\0",")\0\n\0\0","",1,1,0) ;
  if (anB[0])
    {
        sscanf(anB,"%d to %d",&qrySt,&qryEd) ;
	qrySt=fixPsortCoord(qrySt) ;qryEd=fixPsortCoord(qryEd) ;
	ln=qryEd-qrySt+1;
	if (qryEd>lnBuf)strcpy(seq,"recommended to reKantorize") ;
	else {memcpy(seq,seqBuf+qrySt-1,ln) ;seq[ln]=0;}
	if (! (xLen && qrySt < xPos && qryEd > xPos))
	    vstrPrintf(blkp,"Domain \"%s\" \"Psort\" %d. %d %d %d %d \"%s\" \n","Cleavable_signal_peptide",iScore,qrySt,qryEd,sbjSt,sbjEd,seq) ;
        iNum++;
    }
  
  /* transmembrane domains, mieg 2011_02_15 */
  {
    char *cend = 0, ccend = 0 ; /* block and restore */
    char *cp, *cq, *cr, cc ;
    int n, x, y ;
    
    cp = strstr (psortBuf, "Tentative number of TMS(s)") ;
    if (cp)
      {
	cend = strstr (cp, "<A") ;
	if (cend) { ccend = *cend ; *cend = 0 ; }
	while ((cp = strstr(cp, "Transmembrane")))
	  {
	    n = x = y = 0 ;
	    n = sscanf (cp, "Transmembrane  %d - %d\n", &x, &y) ;
	    if (n == 2 && x >= 1 && y >= 1 && x <= lnBuf && y <= lnBuf) /* success */
	      {
		cq = seqBuf + x - 1 ;
		cr = seqBuf + y  ;
		cc = *cr ; *cr = 0 ;
		vstrPrintf (blkp, "Domain Transmembrane_domain Psort 1. %d %d 1 1 \"%s\"\n"
			    , x, y, cq
			    ) ;
		*cr = cc ;
	      }
	    cp += 10 ; /* arbitrarily move a little forward */
	  }
      }
    if (cend) *cend = ccend ;
  }

  /* nuclear localization */
  {
    char *cend = 0, ccend = 0 ; /* block and restore */
    char *cp, *cq, *cr, cc ;
    int n, i, x, y, xmax = 0 ;
    BitSet bb = bitSetCreate (1000, 0) ;    
    /*
     <A HREF="/psort/helpwww2.html#nuc">NUCDISC: discrimination of nuclear localization signals</A>
      pat4: RPRR (4) at    4
      pat4: PRRR (4) at   69
      pat4: RRRK (5) at   70
      pat4: RRKR (5) at   71
      pat4: RKRR (5) at   72
      pat4: KRRR (5) at   73
      pat4: RRRK (5) at   74
      pat4: RRKR (5) at   75
      pat7: PRRRKRR (5) at   69
      pat7: PVATRRR (3) at   91
      bipartite: RRTQISGPRRRKRRRKR at   62
      bipartite: RRPTSTWRWSPVATRRR at   81
      content of basic residues:  24.8%
      NLS Score:  4.86
    */
    
    cp = strstr (psortBuf, "NUCDISC: discrimination of nuclear localization signals") ;
    if (cp)
      {
	cend = strstr (cp, "<A") ;
	if (cend) { ccend = *cend ; *cend = 0 ; }
	while ((cq = strstr(cp, "pat4:")))
	  {
	    cp = cq + 5 ; cq = strstr (cp, "at") ;
	    cr = strstr (cp, "\n") ; 
	    if (cq && cq < cr)
	      {
		cq += 2 ; while (*cq== ' ') cq++ ;
		n = sscanf (cq, "%d", &x) ;
		if (n == 1)
		  {
		    for (i = 0 ; i < 4 ; i++)
		      bitSet (bb, x + i) ;
		    if (x + 1 > xmax) xmax = x + i ;
		  }
	      }
	    cp += 10 ;
	  }
	while ((cq = strstr(cp, "pat7:")))
	  {
	    cp = cq + 4 ; cq = strstr (cp, "at") ;
	    cr = strstr (cp, "\n") ; 
	    if (cq && cq < cr)
	      {
		cq += 2 ; while (*cq== ' ') cq++ ;
		n = sscanf (cq, "%d", &x) ;
		if (n == 1)
		  {
		    for (i = 0 ; i < 7 ; i++)
		      bitSet (bb, x + i) ;
		    if (x + 1 > xmax) xmax = x + i ;
		  }
	      }
	    cp += 10 ;
	  }
    
	while ((cq = strstr(cp, "bipartite:")))
	  {
	    cp = cq + 12;   cq = strstr (cp, "at") ;
	    cr = strstr (cp, "\n") ; 
	    if (cq && cq < cr)
	      {
		y = cq - cp ;
		cq += 2 ; while (*cq== ' ') cq++ ;
		n = sscanf (cq, "%d", &x) ;
		if (n == 1)
		  {
		    for (i = 0 ; i < y ; i++)
		      bitSet (bb, x + i) ;
		    if (x + 1 > xmax) xmax = x + i ;
		  }
	      } 
	  }
      }
    /* export the merged segments */
    bitUnSet (bb, xmax) ;
    for (i = x = y = 0 ; i <= xmax ; i++)
      {
	if (bit (bb, i))
	  {
	    if (x) y++ ;
	    else x = y = i ;
	  }
	else
	  {
	    if (x)
	      {  /* export */
		cq = "" ; cr = 0 ;
		if (x >= 1 && y >= 1 && x <= lnBuf && y <= lnBuf) /* success */
		  {
		    cq = seqBuf + x - 1 ;
		    cr = seqBuf + y ;
		    cc = *cr ; *cr = 0 ;
		  }
		vstrPrintf (blkp,"Domain Nuclear_localization_domain Psort 1. %d %d 1 1 \"%s\"\n"
			    , x, y
			    , cq
			    ) ;
		if (cr) *cr = cc ;
		x = y = 0 ;
	      }
	  }
      }
    bitSetDestroy (bb) ;
    if (cend) *cend = ccend ;
  }
  
  
  /* other domains */
  for(cnt=0,i=0;lkmot[i].st!=0;i++,iNum++)
    {
      char *cp = lkmot[i].nm + strlen(lkmot[i].nm) - 1 ;
      while  (*cp == ' ') *cp-- = 0 ;
      vstrClnMarkup(psortBuf,anB,0,lkmot[i].st,lkmot[i].ed,"",1,1,0) ;
      vstrClnMarkup(anB,anB,0,"(\0\0",")\0\0","",0,0,0) ;
      if (anB[0])cnt+=kantorParsePSORTList(blkp,anB,lkmot[i].fnd,lkmot[i].nm,lkmot[i].seq,lkmot[i].qst,seqBuf) ;
    }
    /* coil coiled regions */
    vstrClnMarkup(psortBuf,anB,0,"detect coiled-coil regions\0\0","total\0\0","",1,1,0) ;
    vstrClnMarkup(anB,anB,0,"</A>\0\0","\0\0","",1,1,0) ;
    clnVar(anB) ;
    if (anB[0])for(clen=0,is=0,ln=2,ptr=anB,qrySt=-2,i=0;ln==2;i++,is++)
      {
        
        if (ptr)ln=sscanf(ptr,"%d %c",&qryEd,(char *)&pep) ;else ln=0;
        if (!i)
	  { qrySt=qryEd ; prvE=qryEd; }
        
        if ((qrySt>=0 && qryEd-prvE>1) || ln!=2)
	  {
	    seq[is]=0;
	    if (! (xLen && qrySt < xPos && qryEd > xPos))
	      vstrPrintf(blkp,"Domain \"%s\" \"Psort\" %d. %d %d %d %d \"%s\" \n","Coiled_coil_region",iScore,fixPsortCoord(qrySt),fixPsortCoord(prvE),sbjSt,sbjEd,seq) ;
	    qrySt=qryEd;cnt++;clen+=is;is=0;
	  }
        if (ln!=2)break;
        ptr=vstrSkipWords(ptr,3,0) ;
        prvE=qryEd;seq[is]=pep;
      }
    
    /* localization */
    vstrClnMarkup(psortBuf,anB,0,"-NN Prediction\0\0",">>\0\0","",1,1,0) ;
    vstrClnMarkup(anB,anB,0,"<PRE>\0\0","\0","",1,1,0) ;
    clnVar(anB) ;
    if (anB[0])for(ptr=anB,i=0;ptr && *ptr;i++)
      {
        if (sscanf(ptr,"%f",&fScore)==0)break;
        ptr=vstrClnMarkup(ptr,seq,0,"%:","\n\0\0","",1,1,0) ;
	clnVar(seq) ;
		vstrFindReplaceStrings(seq,seq,0,"extracellular, including cell wall"_00, "Secreted"_00,0,1) ;
		vstrFindReplaceSymbols(seq,seq,0," ","_",0,1,0) ;
		if (seq[0])vstrPrintf(blkp,"Psort Localization \"%s\" %4.1f \n",seq,fScore) ;
      }
    messfree (anB) ;
    return i;
}


static int parseBlast(char * nm,char * content,char * seq,vSTR * blkp)
{			
  vSTR hints;
  
  if(!vstrSearchSubstring(content,"BLOSUM62\0<PRE><b>BLAST\0</HTML>\0\0" ,1,"",1)){fprintf(stderr,"BROKEN Blastp output for query %s\n",nm);return 1;}
  /*	if(!strstr(content,"<PRE><b>BLASTP")){fprintf(stderr,"BROKEN Blastp output for query %s\n",nm);return 1;}*/
  vSet0(hints);
  
  vstrPrintf(blkp,"\nKantor %s\n",nm);
  if (DELETE_OLD_DATA) vstrPrintf(blkp,"-D BLP\n");	  
  
  kantorParseBLAST(blkp,&hints,content,0);
  vstrPrintf(blkp,"Blastp_Date %s\n",timeShowNow()) ;
  
  vstrPrintf(blkp,"Kantor_Date %s\n",timeShowNow());
  
  vstrPrintf(blkp,"\n");
  
  vstrEmpty(&hints);
  return 1;														   
}
	

/*
* Apparently: 
* The input is a big ASN.1 file.  fAcePieceMakerNamed() breaks it
* apart at some syntactically recognizeable bits, and calls this
* function on discrete subsets of the data.
*
* Here, we run taxblast on one of them, then convert the taxblast
* output into .ace format
*/

char * blastPieceMakerNamed (char * fileName, char * sectionMarker,int * cnt,char * startSearchName,char * endSearchName,char * startSearchSeq,char * endSearchSeq, int (*func)(char * nm,char * cont,char * seq, vSTR *blkp), vSTR *blkp)
{													  
    int     i,ofs=0,st,fn,hFile=-1,pos=0;
    char    * startPiece,*endPiece=0,* ptr ;
    char    nm[vTEXTMAXFILESIZE],seq[vTEXTMAXFILESIZE];
    AC_HANDLE h = 0 ;
							 
    if(!fileName)
      return 0 ;

    /* get the start position in the fasta buffer and in the raw buffer */
   
      {
	if( ! (hFile = vtxtFileOpen (fileName, 0, O_RDONLY)))
	  return 0 ;
	h = handleCreate () ;
	ptr = vtxtFileReadUntil (hFile,&pos,sectionMarker, 0, h) ; 
	messfree (h) ;
	h = handleCreate () ;
	startPiece = vtxtFileReadUntil (hFile,&pos,sectionMarker,1, h);
      }

    for(i=0;(cnt==0 || *cnt==0 || i<(*cnt)) && startPiece;i++){

        /* get the end of the section - and cut it there */
      if(!fileName && (endPiece=strstr(startPiece+ofs,sectionMarker)))*endPiece=0;
      
      ptr=vstrStructuredSearch(startPiece,startSearchName,endSearchName,&st,&fn);
      strncpy(nm,startPiece+st,fn-st);nm[fn-st]=0;
      vstrCleanEnds(nm,nm,0," \t\n",1);
      
      seq[0]=0;
      if(startSearchSeq && endSearchSeq && 
	 (ptr=vstrStructuredSearch(startPiece,startSearchSeq,endSearchSeq,&st,&fn))){
	strncpy(seq,startPiece+st,fn-st);seq[fn-st]=0;
	vstrFindReplaceSymbols(seq,seq,0," \t\n","",0,1,1);
      }
      
      
      if(nm[0] && (func && !func(nm,startPiece,seq, blkp)))break;
      
      /* move to the next section */
      if(fileName){
	messfree(h) ;
	h = handleCreate () ;
	startPiece = vtxtFileReadUntil (hFile,&pos,sectionMarker,1, h);
      }
      else {if(endPiece)*endPiece=sectionMarker[0];startPiece=endPiece;}
    }
    messfree(h) ;
    
    if(fileName)
      {
	vtxtFileClose (hFile) ;
      }
    if(cnt)*cnt=i;
    
    return endPiece;
}



static int parseAsn(char * nm,char * content,char * seq, vSTR *blkp)
{
	char	cmd[2*MAXPATHLEN],
		*resultBufr,flnm[MAXPATHLEN],
		asnf[MAXPATHLEN + 30];

	/*
	* construct a temp file name - I'm using the local host name
	* plus the current pid so that we have some expectation as
	* uniqueness even when running multiple jobs on the compute  
	* farm with the same working directory
	*
	* flnm is a thing that is used as part of the temp file name.
	*/

	gethostname(flnm,sizeof(flnm));
	sprintf(flnm+strlen(flnm),".%d.blastasn.tmp", (int)getpid());

	/*
	* We store some of our data in asnf, then use that as an
	* input file to the command
	*/

	sprintf(asnf,"%s.asnf",flnm);

	/*
	* here is the taxblast command that does that
	*/

	sprintf(cmd,"taxblast -i %s",asnf);

	/*
	* create asnf, containing the bits of asn.1 that we need
	* to parse.
	*/

	vtxtFileSetContent(asnf,content);

	/*
	* run the command; the return string is what the stdout of
	* the command was.
	*/

	resultBufr = vtxtFileExecution (flnm, cmd, "", 0, 0);

	/*
	* we created asnf, so delete it.  (There were also other
	* temp files inside vtxtFileExecution, but they were deleted
	* before it returned.)
	*/

	if (1) remove(asnf);

	/*
	*
	*/

	vstrPrintf(blkp,"\nKantor %s\n",nm);
	if (DELETE_OLD_DATA) 
		vstrPrintf(blkp,"-D TXB\n");
	if(resultBufr && kantorParseTAXBLASTTREE(blkp,resultBufr)){
	  vstrPrintf(blkp,"Taxblast_Date %s\n",timeShowNow());
	  vstrPrintf(blkp,"Kantor_Date %s\n\n",timeShowNow());
	}
	
	messfree(resultBufr);
	
	return 1;
}


static int parsePsort(char * nm,char * content,char * seq, vSTR *blkp)
{  
  char *cp ;
  if(!strstr(content,">> prediction for")){fprintf(stderr,"BROKEN Psort output for query %s\n",nm);return 1;}
  
  vstrPrintf(blkp,"\nKantor %s\n",nm);
  if (DELETE_OLD_DATA) vstrPrintf(blkp,"-D PSR\n");
  
  cp = strstr (content, "<PRE>\n") ;
  if (cp) { cp += 6 ; cp = strstr(cp, "\n") ;}
  if (cp) { cp += 1 ; cp = strstr(cp, "\n") ;}
  if (cp) seq=cp+1 ;
  kantorParsePSORT(blkp,content,seq);
  vstrPrintf(blkp,"Psort_Date %s\n",timeShowNow()) ;
  vstrPrintf(blkp,"Kantor_Date %s\n",timeShowNow());
  
  return 1;
}


int main(int argc, char * argv[],char * envp[])
{
  char    *dstFlnm;
  vSTR	log;
  int		cnt;
  char *operation;

  if (argc != 4)
	{
	fprintf(stderr,"\nError: wrong number of arguments %d (should be 4) on the command line\n", argc);
	fprintf(stderr,"\nThis is the new kantormegaparse\n\n");
	fprintf(stderr,"kantormegaparse operation species inputfile > outputfile\n");
	fprintf(stderr,"operation is one of:\n");
	fprintf(stderr,"\tblast\n");
	fprintf(stderr,"\tblastasn\n");
	fprintf(stderr,"\tpsort\n");
	fprintf(stderr,"species is one of:\n\tworm\tdroso\tara\tbanana\thuman\tmouse\trat\n");
	exit(1);
	}
  
  vSet0(log);
  log.mode|=vSTRECHO;
  log.mode|=vSTRDISCARD;

  operation = argv[1];
  
  selfSpecies = argv[2];
  
  dstFlnm=argv[3];

  if (stat(dstFlnm,&stb))
    {
    perror(dstFlnm);
    exit(1);
    }
    
  if (strcmp(operation,"blastasn") == 0) 
	{
	fprintf(stderr,"\n%s: parsing ASNBLAST    %s\n",dstFlnm,timeShowNow());
	blastPieceMakerNamed(dstFlnm,
		"Seq-annot ::= {",	/* } */
		&cnt,	
		"local\0\"\0\0",
		",\0\"\0\0",
		0,
		0,
		parseAsn,
		&log);
	}
  else if (strcmp(operation, "blast") == 0)
	{
	fprintf(stderr,"\n%s: parsing BLAST    %s\n",dstFlnm,timeShowNow());
	blastPieceMakerNamed(dstFlnm,
		"<b>Query=</b>",
		&cnt,
		"<b>Query=</b>"_00,
		"\n"_00,
		0,
		0,
		parseBlast,
		&log);
	}
  else if (strcmp(operation, "psort") == 0)
	{
	fprintf(stderr,"\n%s: parsing PSORT    %s\n",dstFlnm,timeShowNow());
	blastPieceMakerNamed(dstFlnm,
			     "<H1>Input Sequence</H1>",
			     &cnt,
			     "<PRE>\n"_00,
			     " ("_00,
			     ")"_00,
			     "</PRE>",
			     parsePsort,
			     &log);
	}
  else
	{
	fprintf(stderr,"\nunknown operation: %s\n", operation);
	exit(1);
	}

  vstrPrintf(&log,"\n\n");
  vstrEmpty(&log);

return 0;
}


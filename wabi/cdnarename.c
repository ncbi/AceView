#include <stdio.h>
#include <ctype.h>

#include "acedb.h"
#include "regular.h"
#include "dict.h"
#include "systags.h"
#include "query.h"
#include "parse.h"
#include "cdna.h"
#include "lex.h"
#include "../whooks/tags.h"

static DICT *gene_names_used = 0 ;

/*
* This is a hack.  We do not want to re-use any names that may
* have been assigned by some external entity.  To prevent that,
* we get a list of LocusLink names from <I forget where it came
* from> and a list of Pfam names from some working database.
*
* To make life easier here (so this code can be ready today), I
* put both lists together in a file and strip off all quote marks
* and any ACE file syntax.  The resulting file has one name per
* line.  We should not use any of those names.
*
* The name of this file is provided as a parameter
*/

/***************************************************************************/

static void cdnaRenameInit (char *fName)
{
  int nn = 0 ;
  FILE *f ;
  char b[1000] ;
  gene_names_used = dictCreate (1000) ;

  _LocusLink = str2tag ("LocusLink") ;
  _Pfam = str2tag ("Pfam") ;

  f = filopen (fName,"","r") ;
  if (f)
    while (fgets (b,sizeof (b), f))
      {
	char *s ;
	for (s=b ; *s ; s++)
	  *s = tolower (*s) ;
	s = strchr (b,'\n') ;
	if (s) *s=0 ;
	nn++ ;
	dictAdd (gene_names_used, b, NULL) ;
      }
  filclose (f) ;  /* balance filopen filclose, or fopen fclose */
  messout ("read list of %d reserved names from %s", nn, fName) ;
} /* cdnaRenameInit */

/***************************************************************************/
/* all name reservation are done in full lower case
   although the public names in the data base are mixed upper/lower */
static BOOL myDictAdd (DICT *dict, char *text, int count)
{
  char buf[1000] ;
  char *cp = buf, *cq ;

  strncpy (buf, text, 999) ; buf[999] = 0 ;
  while (*cp) { *cp = tolower (*cp) ; cp++ ; }
     

  switch (count)
    {
    case 0: /* we search name and name.0 */
      if (dictFind (gene_names_used, buf, 0))
	return FALSE ; /* name already known */
      cq = cp ; *cp++ = '.' ; *cp++ = '0' ;
      if (dictFind (gene_names_used, buf, 0))
	return FALSE ; /* name.0 already known */
      *cq = 0 ;  /* we shall insert name, not name.0 in the dict */
      break ;
    case -1:
      cq = cp ; *cp++ = '.' ; *cp++ = '0' ;
      if (dictFind (gene_names_used, buf, 0))
	return FALSE ; /* name.0 already known */
      *cq = 0 ; /* we shall insert name, not name.0 in the dict */
      if (cq > buf + 2 && *(cq-1) == '0' && *(cq-2) == '.') /* name = xx.0 */
	{
	  *(cq-2) = 0 ; /* try to remove the .0 */
	  if (dictFind (gene_names_used, buf, 0))
	    return FALSE ; /* xx already known */
	  *(cq-2) = '.' ; /* reestablish */
	}
      break ;
    default:
      sprintf(cp,".%d", count) ;
      break ;
    }

  return dictAdd (gene_names_used, buf, 0) ; /* FALSE if name already known */
}

/***************************************************************************/

static void cdnaRenameFinish (char *fName)
{
  FILE *f ;
  int x ;
  
  /*
   * save all the newly generated names so that we can know not
   * to re-use them later.
   */
  f=filopen (fName,"","w") ;
  if (!f)
    { 
      messcrash ("cannot save generated names\n") ;
      return ;
    }
  for (x = 1 ; x <= dictMax (gene_names_used) ; x++)
    fprintf (f,"%s\n",dictName (gene_names_used, x)) ;
  filclose (f) ;
  messout ("wrote list of %d reserved names from %s", x, fName) ;
  
  dictDestroy (gene_names_used) ; /* use the macro in the dict.h interface */
} /* cdnaRenameFinish */

/***************************************************************************/
/* in case of several geneid, transfer the geneid into the relevant mRNAs
 * we keep the global acedb name, but on the web, we show the correct sub-name
 */

static void cDnaRenameGeneId2Mrna (KEY gene, KEYSET ks)
{
  Array pgs = arrayCreate (12, HIT) ; /* name and coords of the models */
  Array mrnas = arrayCreate (12, HIT) ; /* name and coords of the models */
  int ii, jj, kk, a1, a2, c1, c2, d1, d2 ;
  OBJ Pg = 0, Mrna = 0 ;
  KEY imap, mrna, pg, geneId, ll ;
  HIT *up, *vp, *bestvp ;
  double z, bestz ;
  KEYSET ks1, ks2, ks3, ks4, ks5 = 0 ;
  
  ks1 = queryKey (gene, ">geneid") ; /* these are relevant */
  ks2 = query (ks, ">geneId") ;    /* these come from the good locuslink names */
  ks3 = keySetAND (ks1, ks2) ;
  ks4 = query (ks3, ">Predicted_gene") ; /* those are the models giving their names */
  if (! keySetMax (ks4))
    goto done ;

  /* collect the coords of the models giving their names */
  for (ii = jj = 0 ; ii < keySetMax (ks4) ; ii++)
    {
      pg = keySet (ks4, ii) ; /* a predicted_gene */
      if ((Pg = bsCreate (pg)))
	{
	  if (bsGetKey (Pg, _IntMap, &imap) &&
	      bsGetData (Pg, _bsRight, _Int, &a1) &&
	      bsGetData (Pg, _bsRight, _Int, &a2))
	    {
	      up = arrayp (pgs, jj++, HIT) ;
	      up->gene = pg ; /* a predicted_gene */
	      if (a1 > a2) { a1 = -a1 ; a2 = -a2 ; }
	      up->a1 = a1 ;
	      up->a2 = a2 ;
	      if (bsGetKey (Pg, str2tag ("GeneId_pg"), &geneId))
		up->est = geneId ;
	    }
	  bsDestroy (Pg) ;
	}
    }
  arraySort (pgs, cDNAOrderGloballyByA1) ;
  
  /* collect the coords of the mrnas */
  ks5 = queryKey (gene, ">transcribed_gene ; > mrna") ;
  for (ii = jj = 0 ; ii < keySetMax (ks5) ; ii++)
    {
      mrna = keySet (ks5, ii) ; /* a mrna */
      if ((Mrna = bsCreate (mrna)))
	{
	  if (bsGetKey (Mrna, _IntMap, &imap) &&
	      bsGetData (Mrna, _bsRight, _Int, &a1) &&
	      bsGetData (Mrna, _bsRight, _Int, &a2))
	    {
	      vp = arrayp (mrnas, jj++, HIT) ;
	      vp->gene = mrna ; /* a mrna */
	      if (a1 > a2) { a1 = -a1 ; a2 = -a2 ; }
	      vp->a1 = a1 ;
	      vp->a2 = a2 ;
	    }
	  bsDestroy (Mrna) ;
	}
    }
  arraySort (mrnas, cDNAOrderGloballyByA1) ;
  /* choose best overlap ratio intersection/union */
  for (ii = 0 ; ii < arrayMax (mrnas) ; ii++)
    {
      up = arrayp (mrnas, ii, HIT) ; 
      mrna = up->gene ;
      if (!(Mrna = bsUpdate (mrna)))
	continue ;
      bestz = 0 ; bestvp = 0 ;
      for (jj = 0 ; jj < arrayMax (pgs) ; jj++)
	{
	  vp = arrayp (pgs, jj, HIT) ;
	  if (!vp->est) continue ;
	  c1 = up->a1 > vp->a1 ? up->a1 : vp->a1 ;
	  c2 = up->a2 < vp->a2 ? up->a2 : vp->a2 ;
	  d1 = up->a1 < vp->a1 ? up->a1 : vp->a1 ;
	  d2 = up->a2 > vp->a2 ? up->a2 : vp->a2 ;
	  z = c2 > c1 ? ((c2 - c1)*1.0)/(d2 - d1 + 1) : 0 ;
	  if (z > bestz) { bestz = z ; bestvp = vp ; }	  
	  if (z > .9)
	    bsAddKey (Mrna, _GeneId, bestvp->est) ;
	}
      if (bestvp && bestvp->est)
	bsAddKey (Mrna, _GeneId, bestvp->est) ;
      bsSave (Mrna) ;	  
    }
  /* reorder the list of name giving locuslink */
  keySetMax (ks) = 0 ;
  for (ii = jj = 0 ; ii < arrayMax (pgs) ; ii++)
    {
      up = arrp (pgs, ii, HIT) ;
      ll = keyGetKey (up->est, _LocusLink) ;
      for (jj = 0 ; ll && jj < keySetMax (ks) ; jj++)
	if (keySet (ks, jj) == ll)
	  ll = 0 ;
      if (ll)
	keySet (ks, jj++) = ll ;
    }
  /* reorder the list inside each mrna */
  for (ii = 0 ; ii < keySetMax (ks5) ; ii++)
    {
      mrna = keySet (ks5, ii) ; /* a mrna */
      keySetDestroy (ks1) ;
      ks1 = queryKey (mrna, ">GeneId") ;
      if ((Mrna = bsUpdate (mrna)))
	{
	  if (bsFindTag (Mrna, _GeneId))
	    bsRemove (Mrna) ;
	  for (jj = 0 ; jj  < arrayMax (pgs) ; jj++)
	    {  
	      up = arrp (pgs, jj, HIT) ;
	      ll = up->est ;
	      for (kk = 0 ; kk < keySetMax (ks1) ; kk++)
		if (keySet(ks1, kk) == ll)
		  bsAddKey (Mrna, _GeneId, ll) ;
	    }
	  bsSave (Mrna) ;
	}
    }
 done:
  arrayDestroy (pgs) ;
  arrayDestroy (mrnas) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
  keySetDestroy (ks4) ;
  keySetDestroy (ks5) ;
} /* cDnaRenameGeneId2Mrna */

/***************************************************************************/
/*
* Given the key of a Gene object, choose a possible new name for it.
*
* If it has a LocusLink field, the name of the locuslink object is
* the new name to use.
*
* If it has Pfam products, find the pfam name with the highest score
* and use that. 
* Otherwise, choose a pseudoword for the new name.  I am just handing
* out the pseudowords in order, skipping anything that is already in
* the dictionary.  Note that if you call this function again for the
* same gene, it will get a different pseudoword.
*
* NOTE: This function is dependent on initialization in cDnaRenameGenes ().
* 
*/

static char *collectLocusLinkNames (KEY gene, int justFromPg)
{
  static char newname[1000] ; /* used to return the result, hence static */ 
  char *cp, *okName = 0 ;
  KEYSET ks = 0 ;
  int ii, jj, tries ;
  KEY locus = 0, locus2 = 0 ;
  BOOL isShed = keyFindTag (gene, str2tag ("Shedded_from")) ;
  /*
   * if the gene has a locuslink name, it uses the locuslink name
   * even if the locuslink name is LOC* (which is usually followed
   * by a number) 
   */
  switch (justFromPg)
    {
    case 1: 
      /* try first a rigorous match */
      ks = queryKey (gene, "> GeneId ; > LocusLink ; ! IS \"* *\" && ! IS \"*:*\"  && ! IS LOC* && ! IS FLJ* && ! IS MGC* ") ;
      if (!keySetMax(ks))
	{
	  keySetDestroy (ks) ;
	  ks = queryKey (gene, ">GeneId ; > LocusLink ; ! IS \"* *\" && ! IS \"*:*\"  ") ;
	}
      break ;
    case 2:
      ks = queryKey (gene, ">locuslink") ; /* obtained from reads aligned in this gene */
      break ;
    case 3:
      /* Construct double names AAandBB; 
	 query used to be before 2006_11_20:
	 
	 >Transcribed_gene ! shedded_from ; >Matching_genefinder_gene ; 
	 COUNT {>Matching_Transcribed_gene ; ! shedded_from } == 1 && NOT COUNT {>model_of; 
	 > gene_model} > 1 ; > GeneId_pg ; > LocusLink ;
	 ! IS \"* *\" && ! IS \"*:*\" && ! IS LOC* && ! IS FLJ* && ! IS MGC*
	 
      */
      ks = queryKey (gene, "> GeneId ; > LocusLink ; ! IS \"* *\" && ! IS \"*:*\" && ! IS LOC* && ! IS FLJ* && ! IS MGC* ") ;
      if (keySetMax(ks) < 2)
	keySetDestroy (ks) ;
      break ;
    default:
      break ;
    }

  if (ks)
    switch (justFromPg)
      {
      case 1: 
      case 2:
	for (ii = 0 ; !okName && ii < keySetMax (ks) ; ii++) 
	  {
	    locus = keySet (ks, ii) ;
	    
	    if (locus && strlen (name (locus)) > 900) /* before arrayCreate */
	      {
		messerror ("Locus name longer than 900 characters\n") ;
		continue ;
	      }
	    
	    if (locus && (strchr (name (locus), ':') || strchr (name (locus), ' ')))
	      continue ;
	    for (tries = isShed ? 1 : 0 ; !okName ; tries++)
	      {
		sprintf (newname, "%s.%d", name (locus), tries) ;
		if (myDictAdd (gene_names_used, name (locus), tries))
		  {
		    okName = strnew (tries ? newname : name (locus), 0) ; 
		  }
	      }
	  } 
	break ;
      case 3: /* composite names */
	/* ncbi now has composite geens named like a    b     a-b
	 * if one of the names is a substring of the other, we forget it
	 */
	for (ii = 0 ; ii < keySetMax (ks) ; ii++) 
	  {
	    if (!keySet (ks, ii))
	      continue ;
	    locus = keySet (ks, ii) ;
	    for (jj = 0 ; jj < keySetMax (ks) ; jj++) 
	      {
		if (ii == jj || !keySet (ks, jj))
		  continue ;
		locus2 = keySet (ks, jj) ;
		if (strstr (name(locus), name(locus2)))
		  keySet (ks, jj) = 0 ;
		else if (strstr (name(locus2), name(locus)))
		  { keySet (ks, ii) = 0 ; break ; }
	      }
	  }
	/* keep happy few */
	for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++) 
	  if (keySet (ks, ii))
	    {
	      if (ii > jj) keySet (ks, jj) = keySet (ks, ii) ;
	      jj++ ;
	    }
	keySetMax (ks) = jj ;
	cDnaRenameGeneId2Mrna (gene, ks) ;
	  
	newname[0] = 0 ;
	for (ii = jj = 0 ; !okName && ii < keySetMax (ks) ; ii++) 
	  {
	    if (!keySet (ks, ii))
	      continue ;
	    locus = keySet (ks, ii) ;
	    
	    if (strlen(newname) + 2 + strlen (name (locus)) > 200) /* before arrayCreate */
	      {
		messerror ("Locus name longer than 200 characters\n") ;
		continue ;
	      }
	    if (jj++) strcat (newname, "--") ;
	    strcat (newname, name (locus)) ;
	  }
	    
	cp = newname + strlen (newname) ;
	for (tries = 0 ; !okName ; tries++)
	  {
	    if (myDictAdd (gene_names_used, newname, tries))
	      {
		if (tries) sprintf (cp, ".%d", tries) ;
		okName = strnew (newname, 0) ; 
	      }
	  } 
	break ;
      }
  
  keySetDestroy (ks) ;
  return okName ;
} /* collectLocusLinkNames  */

/***************************************************************************/

static char *collectPfamNames (KEY gene, BOOL isOld)
{
  KEY high_pfam = 0 ;
  KEYSET ks = 0 ;
  static Array units = 0 ; 
  char newname[1000] ;
  BSunit *uu ;
  char *okName = 0 ;
  int  i, ii ;
  float score = 0, high_score = 0 ;   

  /*
   * if the gene is attached to a PFAM product, find the PFAM
   * with the highest score and use that for the name.
   */
  {
    OBJ Pk ;
    KEY pk ;
  
    ks = queryKey (gene, "follow product ") ;
    
    units = arrayReCreate (units, 12, BSunit) ;
    for (ii = 0 ; ii < keySetMax (ks) ; ii++)
      {
	pk = keySet (ks,ii) ;
	Pk = bsCreate (pk) ;
	if (!Pk)
	  continue ;
	
	if (bsGetArray (Pk, _Pfam, units, 3))
	  for (i = 0 ; i < arrayMax (units) ; i += 3)
	    {
	      uu = arrp (units, i, BSunit) ;
	      score = uu[2].f ;
	      
	      if (score > high_score)
		{
		  high_pfam = uu[0].k ;
		  high_score = score ;
		}
	    }
	bsDestroy (Pk) ;
      }
    keySetDestroy (ks) ;
  }
  
  if (!high_pfam)
    goto done ;
  
  if (isOld) /* try to find this name inside previous names */
    {
      KEY est, chrom, _pnn = str2tag ("previous_newname") ;
      OBJ Est = 0 ;
      char buf[80], *cp, *oldName, *chromnam ;
      int nn = strlen (name(high_pfam)) ;
      
      chrom = keyGetKey (gene, _IntMap) ;
      strncpy (buf, name(chrom), 79) ;
      cp = buf ;
      while (*cp && *cp != '|') cp++ ;
      if (*cp == '|') *cp = 0 ;
      ks = queryKey (gene, ">transcribed_gene ! shedded_from ; >mrna ;{ Mrna_covered_by ; >Mrna_covered_by} SETOR { !Mrna_covered_by ; >cdna_clone ; > read} ; previous_newname") ; 
      for (ii = 0 ; !okName && ii < keySetMax (ks) ; ii++)
	{
	  est  = keySet (ks,ii) ;
	  Est = bsCreate (est) ;
	  if (!Est)
	    continue ;
	  
	  units = arrayReCreate (units, 12, BSunit) ;
	  if (bsGetArray (Est, _pnn, units, 3))
	    for (i = 0 ; !okName &&  i < arrayMax (units) ; i += 3)
	      {
		uu = arrp (units, i, BSunit) ;
		oldName = uu[0].s ;
		/* char *type = uu[1].s ; */
                chromnam = uu[2].s ; 
		/* wrong chromosome */
		cp = chromnam ;
		while (cp && *cp && *cp != '|') cp++ ;
		if (cp && *cp == '|') *cp = 0 ;
		if (!chromnam || strcasecmp (chromnam, buf))
		  continue ;
		/* this pfam was already used last time */
		if (strncasecmp (name(high_pfam), oldName, nn))
		  continue ;
		if (!dictFind (gene_names_used, oldName, 0)) /* not yet used in this build */
		  {
		    char *cp ;
		    myDictAdd (gene_names_used, oldName, -1) ; /* use it */
		    cp = oldName + nn ;
		    if (*cp == '.')
		      okName = strnew (oldName, 0) ; 
		    else  /* dec 2003, always add a .0 at the end of pfam names */
		      okName = strnew (messprintf("%s.0", oldName), 0) ;
		  }
	      }
	  bsDestroy (Est) ;
	}
    }
  else
    { /* loop until we get a non used pfam in this or previous release */
      int c ;
      if (strlen (name (high_pfam)) > 900)
	{
	  /* I don't ever expect to see this error.  */
	  messerror ("pfam name longer than 900 characters\n") ;
	  return NULL ;
	}
      
      for (c = 0 ; !okName  ; c++) 
	{
	  sprintf (newname,"%s.%d",name (high_pfam),c) ;
	  if (myDictAdd (gene_names_used, name (high_pfam), c))
	    okName = strnew (c >= 0 ? newname : name (high_pfam), 0) ;
	}
    }

 done:
  return okName ;
} /* collectPfamNames */

/***************************************************************************/
  /* reserve all previous pfam names, run that on PreviousName Info */
static int reservePreviousNames (char *type)
{ 
  Array units = arrayCreate (12, BSunit) ;
  BSunit *uu ;
  KEYSET ks = query (0, "Find sequence ;  previous_newname") ;
  int i, ii, count = 0 ;
  KEY est, _pnn = str2tag ("previous_newname") ;
  OBJ Est = 0 ;
  char *oldName, *oldType ;
  
  for (ii = 0 ; ii < keySetMax(ks) ; ii++)
    {
      est  = keySet (ks,ii) ;
      Est = bsCreate (est) ;
      if (!Est)
	continue ;
      
      if (bsGetArray (Est, _pnn, units, 2))
	for (i = 0 ; i < arrayMax (units) ; i += 2)
	  {
	    uu = arrp (units, i, BSunit) ;
	    
	    oldName = uu[0].s ;
	    oldType = uu[1].s ;
	    
	    if (type && oldType && !strcasecmp (type, oldType) &&
		myDictAdd (gene_names_used, oldName, -1))  /* reserve it */
	      
	      count++; 
	  }
      bsDestroy (Est) ;
    }
  keySetDestroy (ks) ;
  arrayDestroy (units) ;

  return count ;
} /* reservePreviousNames */

/***************************************************************************/
typedef struct nnStruct { char newname[120] ; KEY est ; int nn ; } NNM ;
static int nnmOrder (const void *va, const void *vb)
{
  const NNM *a = (const NNM *)va, *b = (const NNM *)vb ;
  return strcasecmp (a->newname, b->newname) ;
}
static int nnmOrder2  (const void *va, const void *vb)
{
  const NNM *a = (const NNM *)va, *b = (const NNM *)vb ;
  return b->nn - a->nn ; /* largest first */
}

static char *reusePseudoNames (KEY gene, BOOL is36a)
{
  KEYSET ks = 0 ;
  int i, ii, j, jj, n ;
  Array units = 0 ; 
  KEY est, chrom, _pnn ;
  OBJ Est = 0 ;
  char *cp, buf[80], *type, *oldName, *okName = 0, *chromnam ;
  BSunit *uu ;
  NNM *nnm ;
  Array nnms = arrayCreate (8, NNM) ;

  _pnn = is36a ? str2tag ("previous_newname_36a") : str2tag ("previous_newname") ;
  chrom = keyGetKey (gene, _IntMap) ;
  strncpy (buf, name(chrom), 79) ;
  cp = buf ;
  while (*cp && *cp != '|') cp++ ;
  if (*cp == '|') *cp = 0 ;
#ifdef JUNK
  if (is36a)
    ks = queryKey (gene, ">transcribed_gene  ! shedded_from ; >mrna ;{ Mrna_covered_by ; >Mrna_covered_by} SETOR { !Mrna_covered_by ; >cdna_clone ; >read } ; previous_newname_36a") ; 
  else
    ks = queryKey (gene, ">transcribed_gene  ! shedded_from ; >mrna ;{ Mrna_covered_by ; >Mrna_covered_by} SETOR { !Mrna_covered_by ; >cdna_clone ; >read } ; previous_newname") ; 
#endif
  ks = queryKey (gene, ">transcribed_gene  ! shedded_from ; >read") ;
  
  for (ii = jj = 0 ; !okName && ii < keySetMax (ks) ; ii++)
    {
      est  = keySet (ks,ii) ;
      Est = bsCreate (est) ;
      if (!Est)
	continue ;
      
      units = arrayReCreate (units, 12, BSunit) ;
      if (bsGetArray (Est, _pnn, units, 3))
	for (i = 0 ; !okName &&  i < arrayMax (units) ; i += 3)
	  {
	    uu = arrp (units, i, BSunit) ;
	    oldName = uu[0].s ;
	    type = uu[1].s ;
	    chromnam = uu[2].s ; 
	    /* wrong chromosome */
	    cp = chromnam ;
	    while (cp && *cp && *cp != '|') cp++ ;
	    if (cp && *cp == '|') *cp = 0 ;
	    if (!chromnam || strcasecmp (chromnam, buf))
	      continue ;
	    /* this pfam was already used last time */
	    if (strcasecmp (type, "pseudo"))
	      continue ;
	    if (!dictFind (gene_names_used, oldName, 0)) /* not yet used in this build */
	      {
		nnm = arrayp (nnms, jj++, NNM) ;
		strncpy (nnm->newname, oldName, 119) ;
		nnm->est = est ;
		nnm->nn = 1 ;
	      }
	  }
      bsDestroy (Est) ;
    }

  if (arrayMax (nnms))
    {
      arraySort (nnms, nnmOrder) ;
      for (jj = 0, nnm = arrayp (nnms, jj, NNM) ; jj < arrayMax (nnms) ; nnm++, jj++)
	{
	  oldName = nnm->newname ;
	  for (j = n = 1 ; jj + j <  arrayMax (nnms) ; j++)
	    if (! strcmp (oldName, (nnm+j)->newname))
	      n++ ;
	    else
	      break ;
	  while (j--)
	    (nnm+j)->nn = n ;
	}
      arraySort (nnms, nnmOrder2) ;
      nnm = arrayp (nnms, 0, NNM) ; 
      oldName = nnm->newname ;		   
      myDictAdd (gene_names_used, oldName, -1) ; /* use it */
      okName = strnew (oldName, 0) ; 
    }
  arrayDestroy (nnms) ;
  return okName ;
} /* reusePseudoName  */

/***************************************************************************/
/*   we didn't find a locus or pfam, so make a new name for it. */
static char *createPseudoNames (KEY gene)
{
  int tries ;
  char buffer [80] ;
  BOOL doRandom = FALSE ; /* successive or random names along chromosome */
  static int random_number = -1, e_count = 0, j_count = 0 ;
  enum pseudoword_phoneme_set lang ;
 
  /*
   * Default is an English-like name.  The clone can have a tag
   * "japanese" indicating that it should contribute a 
   * japanese-sounding name to genes is is involved with.  If most 
   * of the clones say "japanese", we will use that.
   *
   * In principle, you could add hints for other language styles
   * to the clones.
   */
  
  /*
   * pick the language
   */
  lang = pseudoword_english ;
  {
    KEYSET ks2 = queryKey (gene, "!NewName ;>Transcribed_gene  ! shedded_from   ; > read ; japanese") ;

    if ( keySetMax (ks2))
      lang = pseudoword_japanese ;
    keySetDestroy (ks2) ;
  }

  /*
   * Pick a name. If it is in the dictionary already, try
   * it again until we get a new one.  It could be in the dictionary
   * because the names were already assigned externally, or because
   * this name was generated by generate_pseudoword () with some
   * other language selected.
   *
   */
  
  tries=1 ;
  if (doRandom || random_number == -1)
    {
      random_number = 1 ;
      /*
	int i = gene % 51 ;
	while (i-- >= 0)
	random_number = randfloat () * 100000 ;
      */
      e_count = j_count = random_number ;
    }
  for ( ; ;)
    {
      char *s ;
      if (doRandom)
	{
	  random_number += 3*7*11*13*17 ; /*  prime with 100000 */
	  e_count = j_count = random_number ;
	}

      if (lang == pseudoword_japanese)
	{
	  if (!doRandom)
	    j_count++ ;
	  if (j_count < pseudoword_japanese_number_base)
	    continue ;
	  random_number = j_count ;
	}
      else
	{
	  if (!doRandom)
	    e_count++ ;
	  if (e_count < pseudoword_english_number_base)
	    continue ;
	  random_number = e_count ;
	}
      
      while (random_number > 500000) random_number -= 500000 ;
      s = generate_pseudoword (buffer, random_number, lang) ;
      if (myDictAdd (gene_names_used, s, 0))
	{
	  /* dictAdd returns true when new added, so we're done */
	  if (tries > 40) printf ("\ntries: %d\n\n",tries) ;
	  return strnew (s, 0) ;
	}
      tries++ ;
      if (tries > 500000)
	{
	  messcrash ("unable to find unused name after 500 000 tries\n") ;
	}
    }
  return 0 ;
} /* createPseudoName */
  
/***************************************************************************/
/*
* The command implementation.  For each gene in the keyset, 
* select a new name for it.  If we actually get a new name (always
* true currently, I think), we make a 2 line "ace file" in memory
* and parse it to add a field called Newname to the object.
*/
static int cDnaRenameSelectNewName (KEYSET genes,  char *fName, int justFromPg)
{
  int ii, count = 0, hasIntron ;
  char *newname = 0, *type = 0 ;
  KEY tg, gene ;
  Stack s = stackCreate (100000) ;
  
  cdnaRenameInit (fName) ;
  
  switch (justFromPg)
    {
    case 4: /* reserve all pfam names from previous build */
      count = reservePreviousNames ("pfam") ;
      break ;
    case 8: /* reserve all pseudo names from previous build */
      count = reservePreviousNames ("pseudo") ;
      break ;
    default:
      /* aug 20, 2005
	 organize the search to give the .0 name to intron gene
      */
      for (hasIntron = 0 ; hasIntron < 2 ; hasIntron++)
	{
	  for (ii = 0 ; ii < keySetMax (genes) ; ii++)
	    {
	      gene = keySet (genes, ii) ; newname = 0 ;

	      if (justFromPg != 10 && keyFindTag (gene, _NewName))
		continue ;

	      tg = keyGetKey (gene, _Transcribed_gene) ;
	      if (! tg)
		continue ;
	      if (keyFindTag (tg, str2tag ("Shedded_from")))
		continue ;
	      switch (hasIntron)
		{
		case 0: /* principal genes with introns */
		  if (!keyFindTag (tg, _gt_ag) && !keyFindTag (tg, _gc_ag))
		    continue ;
		  break ;
		case 1: /* principal genes with out introns */
		  if (keyFindTag (tg, _gt_ag) || keyFindTag (tg, _gc_ag))
		    continue ;
		  break ;
		}
	      
	      switch (justFromPg)
		{
		case 1: /* collect predicted gene name (XM_) */
		  newname = collectLocusLinkNames (gene, 1) ; type = "LocusLink" ;
		  break ;
		case 2: /* collect locuslink from aligned cdna_clone */
		  newname = collectLocusLinkNames (gene, 2) ; type = "LocusLink" ;
		  break ;
		case 3: /* collect justified pfam names */
		  newname = collectPfamNames (gene, TRUE) ; type = "pfam" ;
		  break ;
		case 5: /* collect new pfam names */
		  newname = collectPfamNames (gene, FALSE) ; type = "pfam" ;
		  break ;
		case 6: /* collect justified pseudo names */
		  newname = reusePseudoNames (gene, FALSE) ; type = "pseudo" ;
		  break ;
		case 7: /* collect justified pseudo names of build 36a */
		  newname = reusePseudoNames (gene, TRUE) ; type = "pseudo" ;
		  break ;
		case 9: /* collect new pseudo names */
		  newname = createPseudoNames (gene) ; type = "pseudo" ;
		  break ;
		case 10: /* collect double locuslink name from XM */
		  newname = collectLocusLinkNames (gene, 3) ; type = "DoubleLocusLink" ;
		  break ;
		default:
		  goto abort ;
		}
	      
	      if (newname && *newname)
		{
		  if (0) printf ("%d %s %s\n", ii, name (gene), newname) ;
		  catText (s, messprintf ("Gene %s\nNewname \"%s\" \"%s\"\n\n", name (gene), newname, type)) ;
		  count++ ;
		  messfree (newname) ;
		}
	    }
	  parseBuffer (stackText (s, 0), 0) ;
	}
      break ;
    }
 abort:
  cdnaRenameFinish (fName) ;
  
  stackDestroy (s) ;
  return count ;
}

/***************************************************************************/
/* 2006_11_26
 * new strategy
 *  name as a pool all the mrna of the genebox sorted by CDS length of good product
 *  name all product from thei parent mrna
 *  name the transcribed_gene as their best ranking mrna
 */
typedef struct renameRenameStruct
{ KEY tg, mrna, product ; BOOL isGood, hasIntron ; int ln ; char suffix[32], oldGeneName[1000] ;} RRN ;
static int rrnOrder (const void *va, const void *vb)
{
  const RRN *a = (const RRN *)va, *b = (const RRN *)vb ;
  if (a->isGood != b->isGood)
    return a->isGood ? -1 : 1 ; /* good first */
  return a->ln > b->ln ? -1 : 1 ; /* long first */
  return a->mrna < b->mrna ? -1 : 1 ;
}

/***************************************************************************/

static void cDnaRenameDoRename ( RRN *rrn, int type, const char *newGeneName, char *date, int ntgs)
{
  int ii, jj ;
  KEYSET ks = 0 ;
  KEY key, key2 ;
  char buff1[1000] ;
  char buff2[1060] ;
  char postFix[100] ;
    
  switch (type)
    {
    case 1: 
      {
	OBJ Mrna = bsCreate (rrn->mrna) ;
	ii = 0 ;
	ks = keySetCreate () ;
	for (jj = 0; jj < 4 ; jj++)
	  if (bsGetKey (Mrna, _Product, &key))
	    do
	      {
		switch (jj)
		  {
		  case 0:
		    if (!keyFindTag (key, _Best_product))
		      continue ;
		    break ;
		  case 1:
		    if (keyFindTag (key, _Best_product) ||
			! keyFindTag (key, _Very_good_product)
			)
		      continue ;
		    break ;
		  case 2:
		    if (keyFindTag (key, _Best_product) ||
			keyFindTag (key, _Very_good_product) ||
			! keyFindTag (key, _Good_product)
			)
		      continue ;
		    break ;
		  case 3:
		    if (keyFindTag (key, _Best_product) ||
			keyFindTag (key, _Very_good_product) ||
			keyFindTag (key, _Good_product)
			)
		      continue ;
		    break ;
		  }
		keySet (ks, ii++) = key ;
	      } while (bsGetKey (Mrna, _bsDown, &key)) ;
	bsDestroy (Mrna) ;
	break ;
      }
    case 2: ks = queryKey (rrn->mrna, ">Refseqmaker") ; break ;
    case 3: ks = queryKey (rrn->tg, "IS *") ; break ;
    case 4: ks = queryKey (rrn->mrna, "IS *") ; break ;
    }
      
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      if (ii)
	sprintf (postFix, "%d", ii) ;
      else
	postFix[0] = 0 ;
    
      /*  Construct:  rori.aAug05  7A5.bAug05-unspliced slynar.e2Aug05
       * date : "Dec06" "Aug05" is now inherited from
       *  clone->Aceview_release_suffix 
       */
      if (type == 3 && ntgs == 1)
	sprintf (buff1, "%s"
	       , newGeneName
	       ) ;
      else
	sprintf (buff1, "%s%s.%s%s%s%s"
		 , type == 2 ? "AM_" : ""
		 , newGeneName
		 , rrn->suffix
		 , postFix
		 , date
		 , (type == 1 || type == 4) && !rrn->hasIntron ? "-unspliced" : ""
		 ) ;

      lexAlias (&key, buff1, 0, 0) ;
      switch (type)
	{
	case 1:  /* product */
	  key2 = keyGetKey (key, _Peptide) ;
	  sprintf (buff2, "Product:%s", buff1) ;
	  lexAlias (&key2, buff2,0,0) ;	  
	  break ;
	case 2:  /*.am */
	  key2 = keyGetKey (key, _DNA) ;
	  sprintf (buff2, "%s", buff1) ;
	  lexAlias (&key2, buff2,0,0) ;	  
	  break ;
	case 3: /* tg */
          break ;
	case 4:  /* mrna */
	  key2 = keyGetKey (key, _DNA) ;
	  sprintf (buff2, "mRNA:%s", buff1) ;
	  lexAlias (&key2, buff2,0,0) ;	  
	  break ;
	}
    }
} /* cDnaRenameDoRename */

/***************************************************************************/

static int cDnaRenameRenameOneGene (KEY gene, KEY newName, char *chromName, char *date)
{
  KEYSET tgs = 0, tgsDone = 0, mrnas = 0, products = 0 ;
  Array mrnaArray = arrayCreate (12, RRN) ;
  char newGeneName[1000] ;
  int ii, jj, iTg, iMrna, ln, ntgs = 0 ;
  KEY tg, mrna, product, dummy ;
  RRN *rrn ;
  OBJ Product = 0, Mrna = 0 ;

  /* buffer the new gene name, it is used as the root of all other names */
  if (*chromName)
    sprintf (newGeneName, "ZORGLUB%s.%s", chromName, name (newName)) ; 
  else
    sprintf (newGeneName, "%s", name (newName)) ; 
 
  lexAlias (&gene, newGeneName, 0, 0) ;

  /* collect the properties of all mrnas */
  tgs =  queryKey (gene, ">transcribed_gene") ;
  ntgs = keySetMax (tgs) ;

  for (iTg = iMrna = 0 ; iTg < keySetMax(tgs) ; iTg++)
    {
      tg = keySet (tgs, iTg) ;
      mrnas = queryKey (tg, ">Mrna") ;
      for (ii = 0 ; ii < keySetMax (mrnas) ; ii++)
	{
	  mrna = keySet (mrnas, ii) ;
	  rrn = arrayp (mrnaArray, iMrna++, RRN) ;
	  rrn->mrna = mrna ;
	  rrn->tg = tg ;
	  strncpy (rrn->oldGeneName, name (tg), 999) ;
	  if (keyFindTag (mrna, _gt_ag) || keyFindTag (mrna, _gc_ag))
	    rrn->hasIntron = TRUE ;
	  products = queryKey (mrna, ">Product good_product && best_product") ;
	  if (keySetMax (products)) /* use the best-good product length */
	    {
	      rrn->isGood = TRUE ;
	      for (jj = 0 ; jj < keySetMax (products) ; jj++)
		{
		  product = keySet (products, jj) ;
		  if ((Product = bsCreate (product)))
		    {
		      if (bsGetData (Product, _Coding_length, _Int, &ln) &&
			  ln > rrn->ln)
			rrn->ln = ln ;
		      bsDestroy (Product) ;
		    }
		}
	    }
	  else  /* use the dna length */
	    { 
	      rrn->isGood = FALSE ;
	      if ((Mrna = bsCreate (mrna)))
		{
		  if (bsGetKey (Mrna, _DNA, &dummy) &&
		      bsGetData (Mrna, _bsRight, _Int, &ln) &&
		      ln > rrn->ln)
		    rrn->ln = ln ;
		  bsDestroy (Mrna) ;
		}	  
	    }
	  keySetDestroy (products) ;
	}
    }
  arraySort (mrnaArray, rrnOrder) ;
      
  /* rename the mrnas in order */
  tgsDone = keySetCreate () ;
  for (ii = 0 ; ii < arrayMax (mrnaArray) ; ii++)
    {
      rrn = arrp (mrnaArray, ii, RRN) ;	

      if (ii < 21) /* 21 letters a - u */
	sprintf (rrn->suffix, "%c", 'a' + ii) ;
      else if (ii < 21 + 26) /* va - vz */
	sprintf (rrn->suffix, "v%c", 'a' + ii - 21) ;
      else if (ii < 21 + 2*26) /* wa - wz */
	sprintf (rrn->suffix, "w%c", 'a' + ii - 21 - 26) ;
      else if (ii < 21 + 3*26) /* xa - xz */
	sprintf (rrn->suffix, "x%c", 'a' + ii - 21 - 2*26) ;
      else if (ii < 21 + 4*26) /* ya - yz */
	sprintf (rrn->suffix, "y%c", 'a' + ii - 21 - 3*26) ;
      else if (ii < 21 + 5*26) /* za - zz */
	sprintf (rrn->suffix, "z%c", 'a' + ii - 21 - 4*26) ;
      else  /* need z behind number because prods have a number */
	sprintf (rrn->suffix, "zz_%dz", ii - 21 - 5*26 + 1) ;
  
      /* rename the products */
      cDnaRenameDoRename (rrn, 1, newGeneName, date, ntgs) ; /* products */
      cDnaRenameDoRename (rrn, 2, newGeneName, date, ntgs) ; /* refseqmaker .am */
      /* rename the Tg, but just once */
      if (!keySetFind (tgsDone, rrn->tg, 0))
	{
	  keySetInsert (tgsDone, rrn->tg) ; 
	  cDnaRenameDoRename (rrn, 3, newGeneName, date, ntgs) ; /* tg */
	}
      cDnaRenameDoRename (rrn, 4, newGeneName, date, ntgs) ; /* mrna, should come last */
    }

  keySetDestroy (tgs) ;
  keySetDestroy (tgsDone) ;
  keySetDestroy (mrnas) ;
  keySetDestroy (products) ;
  arrayDestroy (mrnaArray) ;
  return 1 ;
} /* cDnaRenameRenameOneGene */

/***************************************************************************/

static int cDnaRenameRenameGenes (KEYSET ks, char *chromName)
{
  int nn, ii ; 
  KEYSET ks1 ;
  KEY gene, newName ;
  char date[15] ;

  date[0] = 0 ;
  ks1 = query (0, "Find clone Aceview_release_suffix:1") ;
  if (keySetMax(ks1))
    {
      OBJ Clone = bsCreate (keySet (ks1, 0)) ;
      char *cp = 0 ;

      if (Clone)
	{
	  bsGetData (Clone, str2tag ("Aceview_release_suffix"), _Text, &cp) ;
	  if (cp && *cp)
	    strncpy (date, cp, 12) ;
	  bsDestroy (Clone) ;
	}
    }
  nn = 0 ;
  if (!date[0])
    messout ("Missing tag ReleaseDateSuffix in the Strategy clone") ;
  else
    for (ii = 0 ; ii < keySetMax (ks) ; ii++)
      {
	gene = keySet (ks, ii) ;
	newName = keyGetKey (gene, _NewName) ;
	if (newName)
	  {
	    cDnaRenameRenameOneGene (gene, newName, chromName, date) ;
	    nn++ ;
	  }
      }
  keySetDestroy (ks1) ;
  return nn ;
}

/***************************************************************************/
/****************************** public functions ***************************/

int cDnaRenameGenes (KEYSET genes,  char *fName, 
		     char * chromName, BOOL selectNewName, BOOL renameGene)
{
  KEYSET ks = 0 ;
  int nn = 0 ;
  int justFromPg = 0 ;
  
  if (chromName)
    { freeforcecard (chromName) ; freeint (&justFromPg) ; }

  if (selectNewName)
    {
      switch (justFromPg)
	{
	case 10:
	  ks = query (genes, "CLASS Gene") ;
	  break ;
	default:
	  ks = query (genes, "CLASS Gene ; ! NewName") ;
	  break ;
	}
       nn = cDnaRenameSelectNewName (ks, fName, justFromPg) ;
       keySetDestroy (ks) ;
    }

  if (renameGene)
    {
       ks = query (genes, "CLASS Gene ; NewName") ;
       nn = cDnaRenameRenameGenes (ks, chromName) ;
       keySetDestroy (ks) ;
    }

  return nn ;
}

/***************************************************************************/
/***************************************************************************/

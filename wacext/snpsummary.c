/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * 15 Sept 2014
 * snpsummary.c
 *   Input: MetaDB->runs and ali
 *   Output: a wholesale SNP table with all kinds of results
 *      all : all columns
 *        r  readFate
 *        p  pairFate
 *        s  sponge
 */

#include "../wac/ac.h"

typedef struct tsnpStruct {
  const char *project ;
  const char *dbName ;
  const char *export ;
  int snpType ; /* 0: any, 1=sublib; 2=runs ; 3=group */
  int minSnpFrequency, minSnpCover ;
  int allC, allM ;
  ACEOUT ao ;
  KEYSET runs ;
  KEYSET groups ;
  DICT *bloomDict ;
  DICT *runDict ;
  vTXT runsHeader, groupsHeader ;
  AC_DB db ;
  AC_HANDLE h ;
  AC_TABLE snps ;
  const char *Etargets[4] ;
  const char *EtargetsBeau[4] ;
  const char *orderBy ;
  const char *outFileName ;
  BOOL doub ;
  BOOL gzi, gzo ;
  Array rrs ;
  KEYSET covers, mutant ;  /* associated to the current snp */
  KEYSET doubleDetect ;
} TSNP ;

typedef enum { T_ZERO = 0, T_Raw, T_RawA, T_RawKb,  T_Length, T_Seq, T_Read, T_kb, T_Rejected, T_Unaligned, T_aliFrag, T_cFrag, T_uFrag, T_A, T_T, T_G, T_C, T_MAX } T_TYPE ;
typedef enum { F_ZERO = 0, F_Hide, F_p10, F_perCentRead, T_FMAX } T_FORMAT ;
typedef struct snpStruct {
  AC_HANDLE h ;
  AC_OBJ Snp ;
  AC_OBJ Run ;
  float var[T_MAX] ;
} SNP ;

typedef struct runStruct {
  BOOL isGroup ;
  KEYSET r2g ; /* list of groups g of which r is a member */
} RR ;

typedef void (*TSNPFunc)(TSNP *tsnp, SNP *snp) ;
typedef struct mmStruct { char cc ; TSNPFunc f ;} MM ;
typedef struct tagtitleStruct { const char *tag, *title ; int col ; T_TYPE setVar ; T_FORMAT format ;} TT ;

#define EMPTY ""

/*************************************************************************************/
/*************************** actual work *********************************************/
/*************************************************************************************/

static void snpChapterCaption (TSNP *tsnp, TT *tts, const char *caption) 
{
  TT *ti ;
  int n = 0 ;
  char buf[8128] ;
  
  memset (buf, 0, sizeof (buf)) ;

  if (caption)
    aceOutf (tsnp->ao, "\t\t%s", caption) ;

  if (tts)
    {
      for (ti = tts ; ti->tag && n < 8127 ; ti++)
	{
	  const char *cp = ti->title ;
	  
	  buf[n++] = '\t' ; /* count the paragraphs */
	  if (cp)
	    {
	      cp-- ;
	      while ((cp = strchr (cp + 1, '\t')))
		buf[n++] = '\t' ;  /* check for multi columns */
	    }
	}  
      if (n < 2) n = 2 ;
      buf[n-2] = 0 ; /* reserve 2 columns */
      aceOutf (tsnp->ao, "%s", buf) ;
    }
  return ;
} /*  snpChapterCaption */

/*************************************************************************************/

static int snpShowMultiTag (TSNP *tsnp, AC_TABLE tt)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, nn = 0, ir, jr ;
  const char *ccp, *ccq ;
  vTXT txt = vtxtCreate () ;
  Stack s = stackHandleCreate (50, h) ;


  vtxtPrintf (txt, "toto") ;
  /* enter the data backwards, do not duplicate up except for column 1 which we always keep */
  for (ir = tt->rows - 1 ; ir >= 0 ; ir--)
    for (jr = tt->cols - 1 ; jr >= 0 ; jr--)
      {
	ccp = ac_table_printable (tt, ir, jr, 0) ;
	if (ccp)
	  {
	    ccq = hprintf(h, "#%s#", ccp) ;
	    if (
		(jr == 0 && !  strstr (vtxtPtr (txt), ccq) && ( ir == tt->rows - 1 ||  strcmp (ccp, ac_table_printable (tt, ir+1, 0, 0)))) ||
		(jr > 0 && ! strstr (vtxtPtr (txt), ccp))
		)
	      {
		vtxtPrintf (txt, "%s", ccq) ;
		pushText (s, ccp) ;
		nn++ ;
	      }
	  }
      }
  /* export */
  n = nn ;
  while (n--)
    {
      ccp = popText(s) ;
      if (! strcmp (ccp, "polyA")) ccp = "polyA selected" ;
      aceOutf (tsnp->ao, "%s%s"
	       , n == nn - 1 ? "" : ", "
	       , ccp
	       ) ;
      ac_free (txt) ;
    }
  return nn ;
} /* snpShowMultiTag */

/*************************************************************************************/

static int snpShowTable (TSNP *tsnp, AC_TABLE tt, int nw)
{
  int ir, jr ;
  
  if (! tt || ! tt->rows)
    return 0 ;
  
  for (ir = 0 ; ir < tt->rows ; ir++)
    for (jr = 0 ; jr < tt->cols ; jr++)
      {
	const char *ccp = ac_table_printable (tt, ir, jr, 0) ;
	if (ccp)
	  {
	    if (!strncasecmp (ccp, "sra", 3)) ccp += 3 ;
	    if (!strcasecmp (ccp, "sraUnspecified_RNA"))
	      continue ;
	    aceOutf (tsnp->ao, "%s%s", nw++ ? ", " : "",  ccp) ;			  
	  }
      }
  return nw ;
} /* snpShowTable */

/*************************************************************************************/

static void snpShowMillions (TSNP *tsnp, float z)
{
  if (z == 0)
    aceOutf (tsnp->ao, "\t0") ;
  else if (z > 100000 || z < -100000)
    aceOutf (tsnp->ao, "\t%.3f",  z/1000000) ;
  else
    aceOutf (tsnp->ao, "\t%.6f",  z/1000000) ;
} /* snpShowMillions */

/*************************************************************************************/

static void snpShowPercent (TSNP *tsnp, SNP *snp, long int z, BOOL p)
{
  if (p)
    {
      float z1 =  snp->var[T_Read] ? 100 * z / snp->var[T_Read] : 0 ;
      
      if (snp->var[T_Read] > 0) 
	{
	  aceOutf (tsnp->ao, "\t") ;
	  if (0)
	    aceOutf (tsnp->ao, "%f", z1) ;
	  else
	    aceOutPercent (tsnp->ao, z1) ;
	}
      else 
	aceOutf (tsnp->ao, "\t-") ;
    }
  else
    {
      if (z == 0)
	aceOutf (tsnp->ao, "\t0") ;
      else
	aceOutf (tsnp->ao, "\t%ld", z) ;
    }
} /*  snpShowPercent */

/*************************************************************************************/
static void snpShowTag (TSNP *tsnp, SNP *snp, TT *ti)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp = EMPTY ;
  float z = 0 ;

  if (strcmp (ti->tag, "Spacer"))
    {
      ccp = EMPTY ; tt = ac_tag_table (snp->Snp, ti->tag, h) ;
      if (tt && tt->rows && tt->cols >= ti->col)
	{
	  ccp = ac_table_printable (tt, 0, ti->col, EMPTY) ;
	  z = ac_table_float (tt, 0, ti->col, 0) ;
	  if (ti->setVar)
	    snp->var[ti->setVar] = z ;
	}
    }


  switch (ti->format)
    {
    default:
      aceOutf (tsnp->ao, "\t%s", ccp) ;
      break ;
    case F_p10:
      aceOutf (tsnp->ao, "\t%.1f", z/10) ;
      break ;
    case F_perCentRead:
      aceOutf (tsnp->ao, "\t%.4f", snp->var[T_Read] ? 100.0 * z/snp->var[T_Read] : 0) ;
      break ;
    case F_Hide:
      break ;
    }
} /* snpShowTag */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void snpMismatchTypes (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir ;
  float zAny = 0, denominator = 0 ;
  float zTransition = 0, zTransversion = 0 ;
  float zSlidingInsertion = 0, zSlidingDeletion = 0 ;
  float zInsertion = 0, zDeletion = 0 ;
  float zInsertion1 = 0, zInsertion2 = 0, zInsertion3 = 0 ;
  float zInsertionA = 0, zInsertionT = 0, zInsertionG = 0, zInsertionC = 0 ;
  float zDeletion1 = 0, zDeletion2 = 0, zDeletion3 = 0 ; 
  float zDeletionA = 0, zDeletionT = 0, zDeletionG = 0, zDeletionC= 0 ; 
  float zSub[12] ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
 
    { "Compute", "Megabases uniquely aligned and used for mismatch counts", 1, 0, 0} ,
    { "Compute", "Total number of mismatches", 2, 0, 0} ,
    { "Compute", "Mismatches per kb aligned", 21, 0, 0} ,

    { "Compute", "Transitions", 21, 0, 0} ,
    { "Compute", "Transversions", 21, 0, 0} ,

    { "Compute", "1, 2 or 3 bp insertions not in polymers", 12, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions not in polymers", 13, 0, 0} ,
    { "Compute", "1, 2 or 3 bp insertions in polymers", 14, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions in polymers", 15, 0, 0} ,

    { "Compute", "A>G", 21, 0, 0} , /* Transitions */
    { "Compute", "T>C", 21, 0, 0} ,
    { "Compute", "G>A", 21, 0, 0} ,
    { "Compute", "C>T", 21, 0, 0} ,

    { "Compute", "A>T", 21, 0, 0} ,  /* Transversions */
    { "Compute", "T>A", 21, 0, 0} ,
    { "Compute", "G>C", 21, 0, 0} ,
    { "Compute", "C>G", 21, 0, 0} ,
    { "Compute", "A>C", 21, 0, 0} ,
    { "Compute", "T>G", 21, 0, 0} ,
    { "Compute", "G>T", 21, 0, 0} ,
    { "Compute", "C>A", 21, 0, 0} ,

    { "Compute", "Insert A", 21, 0, 0} ,
    { "Compute", "Insert T", 21, 0, 0} ,
    { "Compute", "Insert G", 21, 0, 0} ,
    { "Compute", "Insert C", 21, 0, 0} ,
    { "Compute", "Delete A", 21, 0, 0} ,
    { "Compute", "Delete T", 21, 0, 0} ,
    { "Compute", "Delete G", 21, 0, 0} ,
    { "Compute", "Delete C", 21, 0, 0} ,

    { "Compute", "Single Insertions", 20, 0, 0} ,
    { "Compute", "Single Deletions", 21, 0, 0} ,
    { "Compute", "Double insertions", 22, 0, 0} ,
    { "Compute", "Double deletions", 23, 0, 0} ,
    { "Compute", "Triple insertions", 22, 0, 0} ,
    { "Compute", "Triple deletions", 25, 0, 0} ,
 
    { "Compute", "Transitions per kb", 21, 0, 0} ,
    { "Compute", "Transversions per kb", 21, 0, 0} ,
    { "Compute", "Insertions per kb", 21, 0, 0} ,
    { "Compute", "Deletions per kb", 21, 0, 0} ,

    { "Compute", "% transitions", 21, 0, 0} ,
    { "Compute", "% transversions", 21, 0, 0} ,
    { "Compute", "% indels", 21, 0, 0} ,

    { "Compute", "% insertions not in polymers", 21, 0, 0} ,
    { "Compute", "% deletions not in polymers", 71, 0, 0} ,
    { "Compute", "% insertions in polymers", 21, 0, 0} ,
    { "Compute", "% deletions in polymers", 21, 0, 0} ,

    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "Number of mismatches, per type" ;
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;


  memset (zSub, 0, sizeof (zSub)) ;
  /*
Distribution of mismatches in best unique alignments, absolute, observed, counts
*/

  for (ti = tts ; ti->tag ; ti++)
    {
      char buf[256] ;

      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  aceOutf (tsnp->ao, "\t%s", ac_name(snp->Snp)) ;
	  continue ;
	}

      tt = ac_tag_table (snp->Snp, "Error_profile", h) ;

      if (! strcmp (ti->tag, "Compute"))
	{
	  float zz = 0 ;
	  switch (ti->col)
	    {
	    case 1:   /* Megabases uniquely aligned and used for mismatch counts */
	      for (ir = 0 ; tt &&  ir < tt->rows ; ir++)
		{
		  if (ir == 0)
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz = ac_table_float (tt, ir, 3, 0) ;
		    }
		  else if (strcmp (buf, ac_table_printable (tt, ir, 0, "xxx")))
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz +=  ac_table_float (tt, ir, 3, 0) ;
		    }
		}
	      aceOutf (tsnp->ao, "\t%.0f", zz) ;
	      denominator = zz ;
	      break ;
	      
	    case 2:   /* Total number of mismatches */
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  const char *ccp = ac_table_printable (tt, ir, 1, "xxx") ;
		  zz =  ac_table_float (tt, ir, 2, 0) ;
		  if (! strcasecmp (ccp, "Any"))
		    continue ;
		  zAny += zz ;
		  if (ccp[1] == '>' && ccp[3] == 0) /* substitution */
		    {
		      int i = 0 ;
		      const char *cq, *subTypes = "a>g t>c g>a c>t a>t t>a g>c c>g a>c t>g g>t c>a" ;

		      cq = strstr (subTypes, ccp) ;
		      i = cq ? (cq - subTypes) / 4 : -1 ;
		      if (i >= 4) zTransversion += zz ;
		      else if (i >= 0) zTransition += zz ;

		      if (i >= 0 && i < 12) zSub[i] += zz ;
		    }
		  else if (! strncmp (ccp, "+", 1))
		    {
		      zInsertion += zz ;
		      if (! strncmp (ccp, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "-", 1))
		    {
		      zDeletion += zz ;
		      if (! strncmp (ccp, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*+", 2))
		    {
		      zSlidingInsertion += zz ;
		      if (! strncmp (ccp + 1, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp + 1, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp + 1, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp + 1, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp + 1, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp + 1, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*-", 2))
		    {
		      zSlidingDeletion += zz ;   
		      if (! strncmp (ccp + 1, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp + 1, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp + 1, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp + 1, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp + 1, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp + 1, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		}

	      /* report all at once */
	      aceOutf (tsnp->ao, "\t%.0f", zAny) ;
	      aceOutf (tsnp->ao, "\t%.5f", denominator > 0 ? zAny / (1000 * denominator) : 0) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zTransition, zTransversion) ;

	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion, zDeletion ) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zSlidingInsertion, zSlidingDeletion) ;

	      { 
		int i ; 
		for (i = 0 ; i < 12 ; i++)
		  aceOutf (tsnp->ao, "\t%.0f", zSub[i]) ;
	      }

	      aceOutf (tsnp->ao, "\t%.0f", zInsertionA) ;
	      aceOutf (tsnp->ao, "\t%.0f", zInsertionT) ;
	      aceOutf (tsnp->ao, "\t%.0f", zInsertionG) ;
	      aceOutf (tsnp->ao, "\t%.0f", zInsertionC) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionA) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionT) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionG) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionC) ;

	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion1, zDeletion1) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion2, zDeletion2) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion3, zDeletion3) ;

	      aceOutf (tsnp->ao, "\t%.5f", denominator > 0 ? zTransition / (1000 * denominator) : 0) ;
	      aceOutf (tsnp->ao, "\t%.5f", denominator > 0 ? zTransversion / (1000 * denominator) : 0) ;
	      aceOutf (tsnp->ao, "\t%.5f", denominator > 0 ? (zInsertion + zSlidingInsertion)/ (1000 * denominator) : 0) ;
	      aceOutf (tsnp->ao, "\t%.5f", denominator > 0 ? (zDeletion + zSlidingDeletion)/ (1000 * denominator) : 0) ;

	      if (zAny)
		{
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zTransition / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zTransversion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * (zInsertion + zDeletion + zSlidingInsertion + zSlidingDeletion) / zAny) ;

		  aceOutf (tsnp->ao, "\t%.2f", 100 * zInsertion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zDeletion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zSlidingInsertion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zSlidingDeletion / zAny) ; 
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t\t\t\t") ;
	      break ;
	    default: /* alread  treated in case 2 */
	      break ;

	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* snpMismatchTypes */

/*************************************************************************************/

static void snpSnpTypesDo (TSNP *tsnp, SNP *snp, BOOL isRejected)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir ;
  float zAny = 0, denominator = 0 ;
  float zTransition = 0, zTransversion = 0 ;
  float zSlidingInsertion = 0, zSlidingDeletion = 0 ;
  float zInsertion = 0, zDeletion = 0 ;
  float zInsertion1 = 0, zInsertion2 = 0, zInsertion3 = 0 ;
  float zInsertionA = 0, zInsertionT = 0, zInsertionG = 0, zInsertionC = 0 ;
  float zDeletion1 = 0, zDeletion2 = 0, zDeletion3 = 0 ; 
  float zDeletionA = 0, zDeletionT = 0, zDeletionG = 0, zDeletionC= 0 ; 
  float zSub[12] ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,

    { "Compute", "Total number of variant alleles", 2, 0, 0} ,
    { "Compute", "TSNP support per kb aligned", 21, 0, 0} ,

    { "Compute", "Transitions", 21, 0, 0} ,
    { "Compute", "Transversions", 21, 0, 0} ,

    { "Compute", "1, 2 or 3 bp insertions not in polymers", 12, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions not in polymers", 13, 0, 0} ,
    { "Compute", "1, 2 or 3 bp insertions in polymers", 14, 0, 0} ,
    { "Compute", "1, 2 or 3 bp deletions in polymers", 15, 0, 0} ,

    { "Compute", "A>G", 21, 0, 0} , /* Transitions */
    { "Compute", "T>C", 21, 0, 0} ,
    { "Compute", "G>A", 21, 0, 0} ,
    { "Compute", "C>T", 21, 0, 0} ,

    { "Compute", "A>T", 21, 0, 0} ,  /* Transversions */
    { "Compute", "T>A", 21, 0, 0} ,
    { "Compute", "G>C", 21, 0, 0} ,
    { "Compute", "C>G", 21, 0, 0} ,
    { "Compute", "A>C", 21, 0, 0} ,
    { "Compute", "T>G", 21, 0, 0} ,
    { "Compute", "G>T", 21, 0, 0} ,
    { "Compute", "C>A", 21, 0, 0} ,

    { "Compute", "Insert A", 21, 0, 0} ,
    { "Compute", "Insert T", 21, 0, 0} ,
    { "Compute", "Insert G", 21, 0, 0} ,
    { "Compute", "Insert C", 21, 0, 0} ,
    { "Compute", "Delete A", 21, 0, 0} ,
    { "Compute", "Delete T", 21, 0, 0} ,
    { "Compute", "Delete G", 21, 0, 0} ,
    { "Compute", "Delete C", 21, 0, 0} ,

    { "Compute", "Single Insertions", 20, 0, 0} ,
    { "Compute", "Single Deletions", 21, 0, 0} ,
    { "Compute", "Double insertions", 22, 0, 0} ,
    { "Compute", "Double deletions", 23, 0, 0} ,
    { "Compute", "Triple insertions", 22, 0, 0} ,
    { "Compute", "Triple deletions", 25, 0, 0} ,
 
    { "Compute", "% transitions", 21, 0, 0} ,
    { "Compute", "% transversions", 21, 0, 0} ,
    { "Compute", "% indel", 21, 0, 0} ,

    { "Compute", "% insertions not in polymers", 21, 0, 0} ,
    { "Compute", "% deletions not in polymers", 71, 0, 0} ,
    { "Compute", "% insertions in polymers", 21, 0, 0} ,
    { "Compute", "% deletions in polymers", 21, 0, 0} ,
    { "Compute", "RNA edition: % excess of A>G relative to T>C", 21, 0, 0} ,

    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "Number of rejected variant alleles, per type. Designed for clinical cohorts with sequenced diploid tissues. SNPs which are monomodal across all samples, usually with MAF in the 1-30% range, come most probably from systematic noise and are rejected." ;
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;


  memset (zSub, 0, sizeof (zSub)) ;
  /*
Distribution of mismatches in best unique alignments, abslute, oberved, counts
*/

  for (ti = tts ; ti->tag ; ti++)
    {
      char buf[256] ;

      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s %s", isRejected ? "Rejected " : "",  ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  aceOutf (tsnp->ao, "\t%s", ac_name(snp->Snp)) ;
	  continue ;
	}

      tt = ac_tag_table (snp->Snp, "SNP_profile", h) ;

      if (! strcmp (ti->tag, "Compute"))
	{
	  float zz = 0 ;
	  switch (ti->col)
	    {
	    case 1:   /* Megabases uniquely aligned and used for mismatch counts */
	      for (ir = 0 ; tt &&  ir < tt->rows ; ir++)
		{
		  if (ir == 0)
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz = ac_table_float (tt, ir, 3, 0) ;
		    }
		  else if (strcmp (buf, ac_table_printable (tt, ir, 0, "xxx")))
		    {
		      strcpy (buf, ac_table_printable (tt, ir, 0, "xxx")) ;
		      zz +=  ac_table_float (tt, ir, 3, 0) ;
		    }
		}
	      aceOutf (tsnp->ao, "\t%.0f", zz) ;
	      denominator = zz ;
	      break ;
	      
	    case 2:   /* Total number of mismatches */
	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
		{
		  const char *ccp = ac_table_printable (tt, ir, 1, "xxx") ;
		  zz =  ac_table_float (tt, ir, isRejected ? 3 : 2, 0) ;
		  if (! strcasecmp (ccp, "Any"))
		    continue ;
		  zAny += zz ;
		  if (ccp[1] == '>' && ccp[3] == 0) /* substitution */
		    {
		      int i = 0 ;
		      const char *cq, *subTypes = "a>g t>c g>a c>t a>t t>a g>c c>g a>c t>g g>t c>a" ;

		      cq = strstr (subTypes, ccp) ;
		      i = cq ? (cq - subTypes) / 4 : -1 ;
		      if (i >= 4) zTransversion += zz ;
		      else if (i >= 0) zTransition += zz ;

		      if (i >= 0 && i < 12) zSub[i] += zz ;
		    }
		  else if (! strncmp (ccp, "+", 1))
		    {
		      zInsertion += zz ;
		      if (! strncmp (ccp, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "-", 1))
		    {
		      zDeletion += zz ;
		      if (! strncmp (ccp, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*+", 2))
		    {
		      zSlidingInsertion += zz ;
    		      if (! strncmp (ccp + 1, "+a", 2))
			zInsertionA += zz ;
		      if (! strncmp (ccp + 1, "+t", 2))
			zInsertionT += zz ;
		      if (! strncmp (ccp + 1, "+g", 2))
			zInsertionG += zz ;
		      if (! strncmp (ccp + 1, "+c", 2))
			zInsertionC += zz ;
		      if (! strncmp (ccp + 1, "+++", 3))
			zInsertion3 += zz ;
		      else if (! strncmp (ccp + 1, "++", 2))
			zInsertion2 += zz ;
		      else 
			zInsertion1 += zz ;
		    }
		  else if (! strncmp (ccp, "*-", 2))
		    {
		      zSlidingDeletion += zz ;   
		      if (! strncmp (ccp + 1, "-a", 2))
			zDeletionA += zz ;
		      if (! strncmp (ccp + 1, "-t", 2))
			zDeletionT += zz ;
		      if (! strncmp (ccp + 1, "-g", 2))
			zDeletionG += zz ;
		      if (! strncmp (ccp + 1, "-c", 2))
			zDeletionC += zz ;
		      if (! strncmp (ccp + 1, "---", 3))
			zDeletion3 += zz ;
		      else if (! strncmp (ccp + 1, "--", 2))
			zDeletion2 += zz ;
		      else 
			zDeletion1 += zz ;
		    }
		}

	      /* report all at once */
	      aceOutf (tsnp->ao, "\t%.0f", zAny) ;
	      aceOutf (tsnp->ao, "\t%.5f", denominator > 0 ? zAny / (1000 * denominator) : 0) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zTransition, zTransversion) ;

	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion, zDeletion ) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zSlidingInsertion, zSlidingDeletion) ;

	      { 
		int i ; 
		for (i = 0 ; i < 12 ; i++)
		  aceOutf (tsnp->ao, "\t%.0f", zSub[i]) ;
	      }

	      aceOutf (tsnp->ao, "\t%.0f", zInsertionA) ;
	      aceOutf (tsnp->ao, "\t%.0f", zInsertionT) ;
	      aceOutf (tsnp->ao, "\t%.0f", zInsertionG) ;
	      aceOutf (tsnp->ao, "\t%.0f", zInsertionC) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionA) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionT) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionG) ;
	      aceOutf (tsnp->ao, "\t%.0f", zDeletionC) ;

	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion1, zDeletion1) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion2, zDeletion2) ;
	      aceOutf (tsnp->ao, "\t%.0f\t%.0f", zInsertion3, zDeletion3) ;


	      if (zAny)
		{
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zTransition / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zTransversion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * (zInsertion + zDeletion + zSlidingInsertion + zSlidingDeletion) / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zInsertion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zDeletion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zSlidingInsertion / zAny) ;
		  aceOutf (tsnp->ao, "\t%.2f", 100 * zSlidingDeletion / zAny) ; 
		  if ( zSub[1] > 100)
		    aceOutf (tsnp->ao, "\t%.2f", zSub[0]/zSub[1] -1.0) ;
		  else
		    aceOut (tsnp->ao, "\t") ;
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t\t\t\t\t") ;
	      break ;
	    default: /* alread  treated in case 2 */
	      break ;

	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* snpSnpTypesDo */

/***********/

static void snpSnpTypes (TSNP *tsnp, SNP *snp)
{ 
  snpSnpTypesDo (tsnp, snp, FALSE) ;
}  /* snpSnpTypes */

/***********/

static void snpSnpRejectedTypes (TSNP *tsnp, SNP *snp)
{ 
  snpSnpTypesDo (tsnp, snp, TRUE) ;
}  /* snpSnpRejectedTypes */

/*************************************************************************************/

static void snpSnpCoding (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt, ttS, ttG, ttP ;

  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Number of SNV sites tested\t% of tested SNV sites covered at least 10 times\tMeasured SNV sites\t% rejected monomodals, likely mapping or sequencing errors or RNA-edited sites", 1, 0, 0 } ,
    { "Compute", "Exonic SNV sites\tPure reference SNV (< 5% variant)\tLow frequency SNV (5-20%)\tMid frequency SNV (20-80%)\tHigh frequency SNV (80-95%)\tPure variant SNV (95-100%)", 2, 0, 0} ,
    { "Compute", "Protein changing SNV sites\tPure reference, no change in protein (< 5% variant)\tProtein changing SNV, intermediate (5-95%)\tProtein changing SNV, pure variant (95-100%)", 3, 0, 0} ,
    { "Compute", "Heterozygosity index: heterozygous/homozygous ratio variants", 4, 0, 0 } ,
    /*
    { "Compute", "% Pure reference SNV (< 5% variant)\t% Low frequency SNV (5-20%)\t% Mid frequency SNV (20-80%)\t% High frequency SNV (80-95%)\t% Pure variant SNV (95-100%)", 5, 0, 0} ,
    { "Compute", "% Protein changing SNV sites\t% Pure reference, no change in protein (< 5% variant)\t% Protein changing SNV, intermediate frequency (5-95%)\t% Protein changing SNV, pure variant (95-100%)", 6, 0, 0} ,
    */
    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    "SNV"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  if (snp)
    {
      ttS = ac_tag_table (snp->Snp, "SNP", h) ;
      ttG = ac_tag_table (snp->Snp, "Genomic", h) ;
      ttP = ac_tag_table (snp->Snp, "Protein_changing", h) ;
    }

  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  aceOutf (tsnp->ao, "\t%s", ac_name(snp->Snp)) ;
	  continue ;
	}
            
      tt = ttS ;
      
      if (! strcmp (ti->tag, "Compute"))
	{
	  float z1 = 0, z2 = 0 ;
	  int nT, nC, nR, nnM ;
	  switch (ti->col)
	    {
	    case 1:
	      nT = ac_tag_int (snp->Snp, "Tested_sites", 0) ;
	      nnM = ac_tag_int (snp->Snp, "Not_measurable_sites", 0) ; 
	      nR = ac_tag_int (snp->Snp, "Rejected_sites", 0) ;
	      if (nT > 0)
		{
		  nC = nT - nnM ;
		  z1 = nT ; z1 = 100.0 * nC / z1 ;
		  z2 = nC ; z2 = nC > 0 ? (100.0 * nR / z2) : 0 ;
		  aceOutf (tsnp->ao, "\t%d\t%.2f\t%d\t%.2f", nT, z1, nC, z2) ; 
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t") ;
	      break ;
	    case 2:
	      tt = ttG ;
	      if (tt)
		{
		  aceOutf (tsnp->ao, "\t%d\t%d\t%d\t%d\t%d\t%d"
			   , ac_table_int (tt, 0, 0, 0)
			   , ac_table_int (tt, 0, 2, 0)
			   , ac_table_int (tt, 0, 4, 0)
			   , ac_table_int (tt, 0, 6, 0) 
			   , ac_table_int (tt, 0, 8, 0)
			   , ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t\t\t") ;
	     
	      break ;
	    case 3: 
	      tt = ttP ;
	      if (tt)
		{
		  aceOutf (tsnp->ao, "\t%d\t%d\t%d\t%d"
			   , ac_table_int (tt, 0, 0, 0)
			   , ac_table_int (tt, 0, 2, 0)
			   , ac_table_int (tt, 0, 4, 0) +  ac_table_int (tt, 0, 6, 0) + ac_table_int (tt, 0, 8, 0)
			   , ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t") ;
	      break ;

	    case 4:  
	      nC = nT = 0 ;
	      tt = ac_tag_table (snp->Snp, "Genomic", h) ;
	      if (tt)
		{
		  nC = ac_table_int (tt, 0, 6, 0) +  ac_table_int (tt, 0, 8, 0) ;
		  nT = ac_table_int (tt, 0, 10, 0) ;
		}
	      aceOut (tsnp->ao, "\t") ;  
	      if (nC + nT >= 10)
		{
		  z1 = (nC < 100 * nT ?  nC/(1.0 * nT) : 1000) ;
		  aceOutf (tsnp->ao, "%.2f", z1) ;
		}
		
	      break ;
#ifdef USELESS
	    case 5:
	      tt = ttG ;
	      nT = tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
	      if (nT > 10)
		{
		  float z = 100.0/nT ;
		  aceOutf (tsnp->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			   , z * ac_table_int (tt, 0, 2, 0)
			   , z * ac_table_int (tt, 0, 4, 0)
			   , z * ac_table_int (tt, 0, 6, 0) 
			   , z * ac_table_int (tt, 0, 8, 0)
			   , z * ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t\t") ;
	      break ;
	    case 6:
	      tt = ttG ;
	      nT = tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
	      tt = ttP ;
	      if (nT > 10)
		{
		  float z = 100.0/nT ;
		  aceOutf (tsnp->ao, "\t%.2f\t%.2f\t%.2f\t%.2f"
			   , z * ac_table_int (tt, 0, 0, 0)
			   , z * ac_table_int (tt, 0, 2, 0)
			   , z * (ac_table_int (tt, 0, 4, 0) +  ac_table_int (tt, 0, 6, 0) + ac_table_int (tt, 0, 8, 0))
			   , z * ac_table_int (tt, 0, 10, 0)
			   ) ;
		}
	      else
		aceOut (tsnp->ao, "\t\t\t\t") ;
	      break ;
#endif
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  
  ac_free (h) ;
  return;
} /*  snpSnpCoding */

/*************************************************************************************/

static void snpVCF (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    /* { "Line", "Line Number", 0, 0, 0} , */
    { "VCF", "VCF Chromosome\tVCF position\tVCF ID\tVCF Reference\tVCF Variant", 1, 0, 0} ,
    { "Magic", "Magic Identifier", 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  static int line = 0 ;
  const char *caption =
    "VCF Identifiers"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "Line"))
	{
	  aceOutf (tsnp->ao, "\t%d", ++line) ;
	}
      else if (! strcmp (ti->tag, "VCF"))
	{
	  AC_TABLE tt = ac_tag_table (snp->Snp, "VCF", h) ;
	  if (tt && tt->cols > 2)
	    aceOutf (tsnp->ao
		     , "\t%s\t%d\t.\t%s\t%s"
		     , ac_table_printable (tt, 0, 0, "-")
		     , ac_table_int (tt, 0, 1, 0)
		     , ac_table_printable (tt, 0, 2, "-")
		     , ac_table_printable (tt, 0, 3, "-")
		     ) ;
	  else
	    aceOut (tsnp->ao, "\t\t\t\t\t") ;
	}
      else if (! strcmp (ti->tag, "Magic"))
	aceOutf (tsnp->ao, "\t%s", ac_name (snp->Snp)) ;
      else  /* gName rName pName Dan_Li Typ */
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpVCF */

/*************************************************************************************/

static void snpIdentifiers (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    {  "gName", "genome Name", 0, 0, 0} ,
    {  "rName", "RNA Name", 0, 0, 0} ,  
    {  "pName", "Protein Name", 0, 0, 0} ,  
    {  "Dan_Li", "Dan Li Name", 0, 0, 0} ,  
    {  "Typ", "Type", 0, 0, 0} ,  
    { "Coding", "Coding\tProtein type", 1, 0, 0} ,
    { "Seq_Var", "RNA variation", 10, 0, 0} ,
    { "Seq_Var", "Protein variation", 20, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  static int line = 0 ;
  const char *caption =
    "Variant Identifiers"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "Line"))
	{
	  aceOutf (tsnp->ao, "\t%d", ++line) ;
	}
      else if (! strcmp (ti->tag, "VCF"))
	{
	  AC_TABLE tt = ac_tag_table (snp->Snp, "VCF", h) ;
	  if (tt && tt->cols > 2)
	    aceOutf (tsnp->ao
		     , "\t%s\t%d\t.\t%s\t%s"
		     , ac_table_printable (tt, 0, 0, "-")
		     , ac_table_int (tt, 0, 1, 0)
		     , ac_table_printable (tt, 0, 2, "-")
		     , ac_table_printable (tt, 0, 3, "-")
		     ) ;
	  else
	    aceOut (tsnp->ao, "\t\t\t\t\t") ;
	}
      else if (! strcmp (ti->tag, "Coding"))
	{
	  aceOut (tsnp->ao, "\t") ;
	  if (ac_has_tag (snp->Snp, "Intergenic"))
	    aceOut (tsnp->ao, "Intergenic\t") ;
	  else if (ac_has_tag (snp->Snp, "Intronic"))
	    aceOut (tsnp->ao, "Intronic\t") ;
	  else if (ac_has_tag (snp->Snp, "Non_coding_transcript"))
	    aceOut (tsnp->ao, "Non_coding_transcript\t") ;
	  else if (ac_has_tag (snp->Snp, "UTR_5prime"))
	    aceOut (tsnp->ao, "UTR_5prime\t") ;
	  else if (ac_has_tag (snp->Snp, "UTR_3prime"))
	    aceOut (tsnp->ao, "UTR_3prime\t") ;
	  else if (ac_has_tag (snp->Snp, "Synonymous"))
	    aceOutf (tsnp->ao, "Coding Synonymous\t%s", ac_tag_text (snp->Snp, "Synonymous", "")) ;
	  else if (ac_has_tag (snp->Snp, "AA_substitution"))
	    aceOutf (tsnp->ao, "Coding substitution\t%s", ac_tag_text (snp->Snp, "AA_substitution", "")) ;
	  else if (ac_has_tag (snp->Snp, "Length_variation"))
	    {
	      AC_TABLE tt = ac_tag_table (snp->Snp, "Length_variation", h) ;
	      aceOutf (tsnp->ao, "%s\t%s %s"
		       , ac_table_printable (tt, 0, 0, "")
		       , ac_table_printable (tt, 0, 1, "")
		       , ac_table_printable (tt, 0, 2, "")
		       ) ;
	    }
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else if (! strcmp (ti->tag, "Seq_Var"))
	{
	  int ic ;
	  const char *tag1, *tag2 ;
	  AC_TABLE tt1 = 0, tt2 = 0 ;

	  switch (ti->col)
	    {
	    case 10: 
	      tag1 = "Reference_RNAexon_sequence" ;
	      tag2 = "Observed__RNAexon_sequence" ;
	      ic = 0 ;
	      break ;
	    case 20:
	      tag1 = "Reference_protein_sequence" ;
	      tag2 = "Observed__protein_sequence" ;
	      ic = 1 ;
	    }
	  tt1 = ac_tag_table (snp->Snp, tag1, h) ;
	  tt2 = ac_tag_table (snp->Snp, tag2, h) ;
	  if (tt1 && tt2 && tt1->cols > ic && tt2->cols > ic)
	    aceOutf (tsnp->ao, "\t%s > %s"
		     , ac_table_printable (tt1, 0, ic, "")
		     , ac_table_printable (tt2, 0, ic, "")
		     ) ;
	  else
	    aceOut (tsnp->ao, "\t") ;
	}
      else  /* gName rName pName Dan_Li Typ */
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpIdentifiers */

/*************************************************************************************/

static void snpDanLiCounts (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "AGLR1 m:m+w", 10, 0, 0} ,
    { "Compute", "AGLR2 m:m+w", 11, 0, 0} ,
    { "Compute", "ROCR1 m:m+w", 12, 0, 0} ,
    { "Compute", "ROCR2 m:m+w", 13, 0, 0} ,
    { "Sum", "Any m:m+w", 14, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 
  int mm = 0, cc = 0 ;

  const char *caption =
    "Dan Li counts"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  const char *tag ;
	  char run[6] ;
	  AC_TABLE tt = 0 ;

	  tag = "DanLi_counts" ;
	  strncpy (run, ti->title, 5) ; run[5] = 0 ;
	  tt = ac_tag_table (snp->Snp, tag, h) ;

	  aceOut (tsnp->ao, "\t") ;
	  if (tt)
	    {
	      int ir ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		if (!strcmp (ac_table_printable (tt, ir, 0, "toto"), run))
		  {
		    int m = ac_table_int (tt, ir, 1, 0) ;
		    int c = ac_table_int (tt, ir, 3, 0) ;
		    aceOutf (tsnp->ao, "%d:%d", m, c) ;
		    mm += m ; cc += c ;
		  }
	    }
	}
      else if (! strcmp (ti->tag, "Sum"))
	{
	  aceOutf (tsnp->ao, "\t%d:%d", mm, cc) ;
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpDanLiCounts */

/*************************************************************************************/

static void snpDanLiFrequency (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "AGLR1 m/m+w", 1, 0, 0} ,
    { "Compute", "AGLR2 m/m+w", 2, 0, 0} ,
    { "Compute", "ROCR1 m/m+w", 3, 0, 0} ,
    { "Compute", "ROCR2 m/m+w", 4, 0, 0} ,
    { "Sum", "Any m/m+w", 14, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 
  int mm = 0, cc = 0 ;

  const char *caption =
    "Dan Li counts"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = 0, ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\tDanLi:%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  const char *tag ;
	  char run[6] ;
	  AC_TABLE tt = 0 ;

	  tag = "DanLi_counts" ;
	  strncpy (run, ti->title, 5) ; run[5] = 0 ;
	  tt = ac_tag_table (snp->Snp, tag, h) ;

	  aceOut (tsnp->ao, "\t") ;
	  if (tt)
	    {
	      int ir ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		if (!strcmp (ac_table_printable (tt, ir, 0, "toto"), run))
		  {
		    int m = ac_table_int (tt, ir, 1, 0) ;
		    int c = ac_table_int (tt, ir, 3, 0) ;
		    if  (c == 0) c = 1 ;

		    aceOutf (tsnp->ao, "%.2f", 100.0 * m/c) ;
		    mm += m ; cc += c ;
		  }
	    }
	}
      else if (! strcmp (ti->tag, "Sum"))
	{
	  if (cc == 0) cc = 1 ;
	  aceOutf (tsnp->ao, "\t%.2f", 100.0 * mm / cc) ;
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpDanLiFrequency */

/*************************************************************************************/

static void snpBrsParseCounts (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = 0 ;
  const char *tag = "BRS_counts" ;
  int ir, r, nn = 0 ;
  KEYSET ddd = tsnp->doub ? keySetHandleCreate (h) : 0 ;
  BOOL is20_8 = FALSE ;
  
  tt = ac_tag_table (snp->Snp, tag, h) ;
  for (ir = 0 ; tt && ir < tt->rows ; ir++)
    if (dictFind (tsnp->runDict, ac_table_printable (tt, ir, 0, "toto"), &r))
      {
	RR *rr = arrayp (tsnp->rrs, r, RR) ;
	KEYSET r2g = rr->r2g ;
	int c = ac_table_int (tt, ir, 1, 0) ;
	int m = ac_table_int (tt, ir, 2, 0) ;
	
	if (rr->isGroup) /* avoid double counting */
	  continue ;
	keySet (tsnp->covers, r) = c ;
	keySet (tsnp->mutant, r) = m ;
	tsnp->allC += c ;
	tsnp->allM += m ;
	if (tsnp->doub && 100 * m >= 2 * c)
	  {
	    if (c >= 20 && m >= 4) is20_8 = TRUE ;
	    if (c >= 20 && m >= 4)
	      { keySet (ddd, r) = 1 ; nn++ ; }
	  }
	if (r2g && keySetMax (r2g))
	  {
	    for (int i = 0 ; i < keySetMax (r2g) ; i++)
	      {
		int g = keySet (r2g, i) ;
		keySet (tsnp->covers, g) += c ;
		keySet (tsnp->mutant, g) += m ;
	      }
	  }	      
      }
  if (is20_8)
    {
      int N = dictMax (tsnp->runDict) + 1 ;
      for (int r1 = 1 ; r1 <= N ; r1++)
	for (int r2 = 1 ; r2 <= N ; r2++)
	  if (nn >= 2 || (keySet (ddd, r1) > 0 && keySet (ddd, r2) > 0))
	    keySet (tsnp->doubleDetect, N *r1 + r2) += 1 ;
    }
  ac_free (h) ;
  return;
} /*  snpBrsParseCounts */

/*************************************************************************************/

static void snpBrsFrequencyCounts (TSNP *tsnp, SNP *snp, BOOL isFrequency, BOOL isGroup)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { isGroup ? "Groups" : "Runs", isGroup ? vtxtPtr (tsnp->groupsHeader) : vtxtPtr (tsnp->runsHeader), 1, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Allele Frequency or counts measured by MAGIC if more than 10 covering reads"
    ;

  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      KEYSET rg = 0 ;

      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      if (! keySetMax (tsnp->covers)) /* parse the counts */
	snpBrsParseCounts (tsnp, snp) ;
      
      if (! strcmp (ti->tag, "Runs"))
	rg = tsnp->runs ;
      if (! strcmp (ti->tag, "Groups"))
	rg = tsnp->groups ;
      
      if (rg)
	{
	  aceOut (tsnp->ao, "\t") ;
	  for (int ir = 0 ; ir < keySetMax (rg); ir++)
	    {
	      int r = keySet (rg, ir) ;
	      int c = keySet (tsnp->covers, r) ;
	      int m = keySet (tsnp->mutant, r) ;
	      
	      if (!r)
		continue ;
	      if (isFrequency)
		{
		  if (c >= 10)
		    aceOutf (tsnp->ao, "\t%.2f", 100.0*m/c) ;
		  else
		    aceOutf (tsnp->ao, "\t-10") ;
		}
	      else
		aceOutf (tsnp->ao, "\t%d:%d", m, c) ;
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpBrsFrequencyCounts */

static void snpBrsRunFrequency (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, TRUE, FALSE) ;
}
static void snpBrsRunCounts (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, FALSE, FALSE) ;
}
static void snpBrsGroupFrequency (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, TRUE, TRUE) ;
}
static void snpBrsGroupCounts (TSNP *tsnp, SNP *snp)
{
  return snpBrsFrequencyCounts (tsnp, snp, FALSE, TRUE) ;
}


/*************************************************************************************/

static void snpBrsAllRuns (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Sum", "All mutant\tAll coverage\tAll frequency", 10, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Allele counts and frequency when summing all runs in this table"
    ;

  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      KEYSET rg = 0 ;

      if (snp == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	  continue ;
	}
      if (! keySetMax (tsnp->covers)) /* parse the counts */
	snpBrsParseCounts (tsnp, snp) ;

      if (ti->col == 10)
	{	 
	  int c = tsnp->allC ;
	  int m = tsnp->allM ;

	  aceOutf (tsnp->ao, "\t%d\t%d", m, c) ;
	  if (c >= 10)
	    aceOutf (tsnp->ao, "\t%.2f", 100.0*m/c) ;
	  else
	    aceOutf (tsnp->ao, "\t-10") ;
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }
  ac_free (h) ;
  return;
}


/*************************************************************************************/
/* template for a future chapter */
static void snpOther (TSNP *tsnp, SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Statistics before alignment"
    ;
  if (snp == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (snp == 0)
	aceOutf (tsnp->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  aceOutf (tsnp->ao, "\t%s", ac_name(snp->Snp)) ;
	}
      else  if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    default:
	      break ;
	    }
	}
      else
	snpShowTag (tsnp, snp, ti) ;
    }

   ac_free (h) ;
  return;
}  /* snpOther */

/*************************************************************************************/

static BOOL snpFilter (TSNP *tsnp, SNP *snp)
{
  BOOL ok = TRUE ;
  int minF = tsnp->minSnpFrequency ;
  int minC = tsnp->minSnpCover ;


  if (tsnp->snpType < 3 && minF + minC > 0)
    {
      AC_HANDLE h = ac_new_handle () ;
      AC_TABLE tt = ac_tag_table (snp->Snp, "BRS_counts", h) ;
      int ir ;

      ok = FALSE ;
      if (tt && tt->rows)
	for (ir = 0 ; ir < tt->rows && !ok ; ir++)
	  {
	    int c = ac_table_int (tt, ir, 1, 0) ;
	    int m = ac_table_int (tt, ir, 2, 0) ;

	    if (c >= minC && 100 * m >= minF * c)
	      ok = TRUE ;
	  }

      ac_free (h) ;
    }
  return ok ;
} /* snpFilter */

/*************************************************************************************/

static const char *allMethods = "VISdgrDGR" ;

static MM methods [] = {
  {'V', &snpVCF} ,
  {'I', &snpIdentifiers} ,
  {'d', &snpDanLiFrequency} ,
  {'D', &snpDanLiCounts} ,
  {'r', &snpBrsRunFrequency} , 
  {'R', &snpBrsRunCounts} ,
  {'g', &snpBrsGroupFrequency} , 
  {'G', &snpBrsGroupCounts} ,
  {'S', &snpBrsAllRuns} ,
{ 0, 0 }
} ;

/*************************************************************************************/

static BOOL snpExportSnp (TSNP *tsnp, int iSnp, int type)
{       
  const char *ccp, *lineName ;
  MM *mm ;
  SNP *snp = 0 ;
  SNP snp0 ;
  
  memset (&snp0, 0, sizeof(SNP)) ;

  switch (type)
    {
    case 0: 
      lineName =  "### Caption" ; 
      snp = (void *) 1 ;
      break ;
    case 1: 
      lineName =  "### MAGIC SNP" ;
      snp = (void *) 0 ;
      break ;
    case 2: 
      snp = &snp0 ;
      snp->h = ac_new_handle () ;
      tsnp->covers = arrayHandleCreate (dictMax (tsnp->runDict) + 1, KEY, snp->h) ;
      tsnp->mutant = arrayHandleCreate (dictMax (tsnp->runDict) + 1, KEY, snp->h) ;
      snp->Snp = ac_table_obj (tsnp->snps, iSnp, 0, snp->h) ;
      tsnp->allC = tsnp->allM = 0 ;
      lineName = "" ;

      break ; 
    }
   
  if (type < 2 || snpFilter (tsnp, snp))
    {
      aceOutf (tsnp->ao, "%s", lineName) ;
      
      ccp = tsnp->export - 1 ;
      while (*++ccp)
	{
	  mm = methods - 1 ;
	  while (mm++, mm->cc)
	    if (mm->cc == *ccp)
	      mm->f(tsnp, snp) ;
	}
      aceOutf (tsnp->ao, "\n" ) ;
    }
  ac_free (snp0.h) ;

  return TRUE ;
} /* snpExporTSNP */

/*************************************************************************************/
static void snpExportDoubleDetect (TSNP *tsnp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsnp->outFileName, ".double_counts.txt", tsnp->gzo, h) ;
  aceOutDate (ao, "##", tsnp->dbName) ;
  aceOutDate (ao, "##", tsnp->project) ;

  int N = dictMax (tsnp->runDict) + 1 ;
  for (int r1 = 1 ; r1 < N ; r1++)
    {
      aceOut (ao, "\n") ;
      for (int r2 = r1 ; r2 <  N ; r2++)
	aceOutf (ao, "%s\t%s\t%d\n"
		 , dictName (tsnp->runDict, r1)
		 , dictName (tsnp->runDict, r2)
		 , keySet (tsnp->doubleDetect, N *r1 + r2) 
		 ) ;
    }
  ac_free (h) ;
} /* snpExportDoubleDetect */

/*************************************************************************************/

static int snpGetSnpsRunsGroups (TSNP *tsnp)
{
  int ns = 0 ;
  const char *errors = 0 ;
  char *qq ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = 0 ;
  KEYSET ks ;

  tsnp->runs = keySetHandleCreate (tsnp->h) ;
  tsnp->groups = keySetHandleCreate (tsnp->h) ;
  DICT *dict ;
  vTXT txt ; 

  dict = tsnp->runDict = dictHandleCreate (256, tsnp->h) ;

  if (tsnp->project)
    {
      int r, g, nn = 0 ;

      tsnp->rrs = arrayHandleCreate (200, RR, tsnp->h) ;

      qq = hprintf (h, "select r, t from p in ?project where p == \"%s\" , r in p->run where r ISA runs, t in r->%s", tsnp->project, tsnp->orderBy) ;
      tt = ac_bql_table (tsnp->db
			 , qq
			 , 0
			 , "+2"
			 , &errors
			 , h
			 ) ;

      txt = tsnp->runsHeader = vtxtHandleCreate (tsnp->h) ;
      ks = tsnp->runs ; nn = 0 ;
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tt, ir, 0, "toto") ;
	    dictAdd (dict, ccp, &r) ;
	    vtxtPrintf (txt, "\t%s", ccp) ;
	    keySet (ks, nn++) = r ;
	  }
      else
	vtxtPrintf (txt, "\t%s", "None") ;

      qq = hprintf (h, "select r, t from p in ?project where p == \"%s\", r in p->run where r ISA Groups, t in r->%s", tsnp->project, tsnp->orderBy) ;
      tt = ac_bql_table (tsnp->db
			 , qq
			 , 0
			 , "+2"
			 , &errors
			 , h
			 ) ;
      txt = tsnp->groupsHeader = vtxtHandleCreate (tsnp->h) ;
      ks = tsnp->groups ; nn = 0 ;
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tt, ir, 0, "toto") ;
	    dictAdd (dict, ccp, &g) ;
	    RR *rr = arrayp (tsnp->rrs, g, RR) ;
	    vtxtPrintf (txt, "\t%s", ccp) ;
	    keySet (ks, nn++) = g ;
	    rr->isGroup = TRUE ;
	  }
      else
	vtxtPrintf (txt, "\t%s", "None") ;

      qq = hprintf (h, "select g, r from p in ?project where p == \"%s\", g in p->run where g#union_of, r in g->union_of, p2 in r->project where p2 == p", tsnp->project) ;
      tt = ac_bql_table (tsnp->db
			 , qq
			 , 0
			 , 0
			 , &errors
			 , h
			 ) ;
      
      /* rr->r2g: list of groups g of which r is a member */
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *gNam = ac_table_printable (tt, ir, 0, "toto") ;
	    int g = 0 ;
	    if (dictFind (tsnp->runDict, gNam, &g))
	      {
		int r = 0 ;
		const char *rNam = ac_table_printable (tt, ir, 1, "toto") ;

		if ( dictFind (tsnp->runDict, rNam, &r))
		  {
		    RR *rr = arrayp (tsnp->rrs, r, RR) ;
		    KEYSET ks = rr->r2g ;
		    if (!ks) 
		      ks = rr->r2g = keySetCreate () ;
		    keySet (ks, keySetMax (ks)) = g ;
		  }
	      }
	  }
    }

  switch (tsnp->snpType)
    {
    case 0: 
      qq ="select s from s in ?Variant  " ;
      break ;
    case 3 : 
      qq = "select s from s in ?Variant where s#danli_counts " ;
      break ;
    case 4 : 
      qq = "select s from s in ?Variant, m in s->danli_counts where m, x in m[1], y in m[3] where 100 * x > 98 * y" ;
      break ;
    }

  tsnp->snps = ac_bql_table (tsnp->db
			     , qq
			     , 0
			     , 0
			     , &errors
			     , tsnp->h
			     ) ;
  
  ns = tsnp->snps ? tsnp->snps->rows : 0  ;
  
  fprintf (stderr, " snpGetSnps got %d snps, %d runs , %d groups in project %s\n"
	   , ns
	   , tsnp->runs ? keySetMax (tsnp->runs) : 0
	   , tsnp->groups ? keySetMax (tsnp->groups) : 0 
	   , tsnp->project
	   ) ;
  
  if (tsnp->doub)
    tsnp->doubleDetect = keySetHandleCreate (tsnp->h) ;

  ac_free (h) ;
  return ns ;
} /* snpGetSnpsRunsGroups */

/*************************************************************************************/

static void snpExportHeader (TSNP *tsnp)
{
  aceOutf (tsnp->ao, "## Project %s, quality control report, File %s : %s\n"
	   , tsnp->project
	   , aceOutFileName (tsnp->ao)
	   , timeShowNow () 
	   ) ;

  return ;
} /* snpExportHeader */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// snpsummary: complete SNP summary for the Magic pipeline\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2022, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a project name\n"
	    "//    import all data from TSNP_DB, export a single SNP table\n"
	    "//    one line per snp, with snp-metadata, then groups or runs counts in comulm chapters\n"
	    "//    a smaller table can be obtained using different options\n"
	    "//\n"
	    "// Syntax:\n"
	    "// snpsummary -db <dir> [options]\n"
	    "//   -db <dir>: snp-database directory (usually: tmp/TSNP_DB/$zone)\n"
	    "// All other parameters are optional and may be specified in any order\n"
	    "//   -p\n" 
	    "//   --project <projectName> export counts to runs and groups belonging to the project (usually $MAGIC)\n" 
	    "//      default : only export the identifiers\n"
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file_prefix.action\n" 
            "//   --gzo : gzip all output files\n"
	    "//   -e\n" 
	    "//   --export [TPbafpmg...] : only export some groups of columns, in the requested order\n"
	    "//      default: if -export is not specified, export all columns in default order\n"
	    "//            I: snp identifiers\n"
	    "//            G: groups counts\n"
	    "//            g: groups allele frequencies\n"
	    "//            R: runs counts\n"
	    "//            r: runs allele frequencies\n"
	    "//            D: Dan Li counts\n"
	    "//            d: Dan Li allele frequencies\n"
	    "//   --orderBy <tag> : sort the runs by this tag\n"
	    "//      default: sorting_title\n"
	    "//      Only export in that order the runs belonging to the project, present in this list\n"
	    "//      Insert a blank line if the sorting value looks like A__B and A changes\n"
	    "//   --snpType [0,1,2,3] : which runs\n"
	    "//       0 [default]: any_snp\n"
	    "//       1 : skip the rejected snps\n"
	    "//       2 : only the rejected snps\n"
	    "//       3 : only the DanLi snps\n"
	    "//       9 : no SNP, only the title lines\n"
	    "//    --doubleDetect\n"
	    "//       Count  snp detetced at in 2 libs at 10:4 and at least once at 20:8\n"
	    "// Caveat:\n"
	    "//   Caption lines at the top start with a #\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  snpsummary -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  TSNP tsnp ;
  AC_HANDLE h = 0 ;
  int ns = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&tsnp, 0, sizeof (TSNP)) ;
  tsnp.h = h ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;
  tsnp.dbName = "NULL" ;
  getCmdLineOption (&argc, argv, "-db", &tsnp.dbName) ; 
  getCmdLineOption (&argc, argv, "--db", &tsnp.dbName) ; 

  tsnp.export =  allMethods ;
  tsnp.project = "NULL" ;
  getCmdLineOption (&argc, argv, "-e", &tsnp.export) ;
  getCmdLineOption (&argc, argv, "--export", &tsnp.export) ;
  getCmdLineOption (&argc, argv, "-p", &tsnp.project) ;
  getCmdLineOption (&argc, argv, "--project", &tsnp.project) ;
  getCmdLineInt (&argc, argv, "--snpType", &tsnp.snpType) ;

  tsnp.doub =   getCmdLineBool (&argc, argv, "--doubleDetect") ;
  tsnp.gzo =   getCmdLineBool (&argc, argv, "--gzo") ;

  tsnp.minSnpFrequency = 0 ;
  tsnp.minSnpCover = 0 ;
  getCmdLineInt (&argc, argv, "--minSnpFrequency", &tsnp.minSnpFrequency) ;
  getCmdLineInt (&argc, argv, "--minSnpCover", &tsnp.minSnpCover) ;
  if (! tsnp.dbName)
    {
      fprintf (stderr, "Sorry, missing parameter -db, please try snpsummary -help\n") ;
      exit (1) ;
    }
  if (tsnp.dbName)
    {
      const char *errors ;

      tsnp.db = ac_open_db (tsnp.dbName, &errors);
      if (! tsnp.db)
	{
	  fprintf (stderr, "Failed to open db %s, error %s", tsnp.dbName, errors) ;
	  exit (1) ;
	}
    }
  tsnp.orderBy = "Sorting_title" ;
  getCmdLineOption (&argc, argv, "--orderBy", &tsnp.orderBy) ;
  
  getCmdLineOption (&argc, argv, "-o", &tsnp.outFileName) ; 
  tsnp.ao = aceOutCreate (tsnp.outFileName, ".SNP_summary.txt", tsnp.gzo, h) ;
  aceOutDate (tsnp.ao, "##", tsnp.dbName) ;

  if (! tsnp.ao)
    messcrash ("cannot create output file %s\n", tsnp.outFileName ? tsnp.outFileName : EMPTY ) ;


  if (argc > 1) usage (messprintf ("Unknown parameters %s", argv[1])) ;

  snpExportHeader (&tsnp) ;

  ns = snpGetSnpsRunsGroups (&tsnp) ;
  snpExportSnp (&tsnp, 0, 0) ; /* chapter captions */
  snpExportSnp (&tsnp, 0, 1) ; /* column titles */

  if (ns)
    {
      int iSnp ;

      for (iSnp = 0 ; iSnp < ns ; iSnp++)
	snpExportSnp (&tsnp, iSnp, 2) ;
    }
  
  if (tsnp.doub)
    snpExportDoubleDetect (&tsnp) ;

  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


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

typedef struct snpStruct {
  const char *project ;
  const char *dbName ;
  const char *export ;
  int snpType ; /* 0: any, 1=sublib; 2=runs ; 3=group */
  int minSnpFrequency, minSnpCover ;
  ACEOUT ao ;
  Array runs ;
  DICT *bloomDict ;
  DICT *runsDict ;
  DICT *groupsDict ;
  vTXT runsHeader, groupsHeader ;
  AC_DB db ;
  AC_HANDLE h ;
  AC_TABLE truns ;
  AC_TABLE groups ;
  AC_TABLE snps ;
  const char *Etargets[4] ;
  const char *EtargetsBeau[4] ;
  const char *orderBy ;
  const char *outFileName ;
  BOOL gzi, gzo ;
  Array r2gs ;
} TSNP ;

typedef enum { T_ZERO = 0, T_Raw, T_RawA, T_RawKb,  T_Length, T_Seq, T_Read, T_kb, T_Rejected, T_Unaligned, T_aliFrag, T_cFrag, T_uFrag, T_A, T_T, T_G, T_C, T_MAX } T_TYPE ;
typedef enum { F_ZERO = 0, F_Hide, F_p10, F_perCentRead, T_FMAX } T_FORMAT ;
typedef struct rcStruct {
  AC_HANDLE h ;
  AC_OBJ snp ;
  AC_OBJ run ;
  float var[T_MAX] ;
} RC ;

typedef struct r2gStruct {
  AC_HANDLE h ;
  int run ;
  KEYSET r2g ;
  KEYSET g2r ;
} R2G ;

typedef void (*TSNPFunc)(TSNP *tsnp, RC *rc) ;
typedef struct mmStruct { char cc ; TSNPFunc f ;} MM ;
typedef struct tagtitleStruct { const char *tag, *title ; int col ; T_TYPE setVar ; T_FORMAT format ;} TT ;

#define EMPTY ""

/*************************************************************************************/
/*************************** actual work *********************************************/
/*************************************************************************************/

static void snpChapterCaption (TSNP *snp, TT *tts, const char *caption) 
{
  TT *ti ;
  int n = 0 ;
  char buf[8128] ;
  
  memset (buf, 0, sizeof (buf)) ;

  if (caption)
    aceOutf (snp->ao, "\t\t%s", caption) ;

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
      aceOutf (snp->ao, "%s", buf) ;
    }
  return ;
} /*  snpChapterCaption */

/*************************************************************************************/

static int snpShowMultiTag (TSNP *snp, AC_TABLE tt)
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
      aceOutf (snp->ao, "%s%s"
	       , n == nn - 1 ? "" : ", "
	       , ccp
	       ) ;
      ac_free (txt) ;
    }
  return nn ;
} /* snpShowMultiTag */

/*************************************************************************************/

static int snpShowTable (TSNP *snp, AC_TABLE tt, int nw)
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
	    aceOutf (snp->ao, "%s%s", nw++ ? ", " : "",  ccp) ;			  
	  }
      }
  return nw ;
} /* snpShowTable */

/*************************************************************************************/

static void snpShowMillions (TSNP *snp, float z)
{
  if (z == 0)
    aceOutf (snp->ao, "\t0") ;
  else if (z > 100000 || z < -100000)
    aceOutf (snp->ao, "\t%.3f",  z/1000000) ;
  else
    aceOutf (snp->ao, "\t%.6f",  z/1000000) ;
} /* snpShowMillions */

/*************************************************************************************/

static void snpShowPercent (TSNP *snp, RC *rc, long int z, BOOL p)
{
  if (p)
    {
      float z1 =  rc->var[T_Read] ? 100 * z / rc->var[T_Read] : 0 ;
      
      if (rc->var[T_Read] > 0) 
	{
	  aceOutf (snp->ao, "\t") ;
	  if (0)
	    aceOutf (snp->ao, "%f", z1) ;
	  else
	    aceOutPercent (snp->ao, z1) ;
	}
      else 
	aceOutf (snp->ao, "\t-") ;
    }
  else
    {
      if (z == 0)
	aceOutf (snp->ao, "\t0") ;
      else
	aceOutf (snp->ao, "\t%ld", z) ;
    }
} /*  snpShowPercent */

/*************************************************************************************/
static void snpShowTag (TSNP *snp, RC *rc, TT *ti)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp = EMPTY ;
  float z = 0 ;

  if (strcmp (ti->tag, "Spacer"))
    {
      ccp = EMPTY ; tt = ac_tag_table (rc->snp, ti->tag, h) ;
      if (!tt) tt = ac_tag_table (rc->run, ti->tag, h) ;
      if (tt && tt->rows && tt->cols >= ti->col)
	{
	  ccp = ac_table_printable (tt, 0, ti->col, EMPTY) ;
	  z = ac_table_float (tt, 0, ti->col, 0) ;
	  if (ti->setVar)
	    rc->var[ti->setVar] = z ;
	}
    }

  switch (ti->format)
    {
    default:
      aceOutf (snp->ao, "\t%s", ccp) ;
      break ;
    case F_p10:
      aceOutf (snp->ao, "\t%.1f", z/10) ;
      break ;
    case F_perCentRead:
      aceOutf (snp->ao, "\t%.4f", rc->var[T_Read] ? 100.0 * z/rc->var[T_Read] : 0) ;
      break ;
    case F_Hide:
      break ;
    }
} /* snpShowTag */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* Data characteristics before alignment */
static void snpBeforeAli (TSNP *snp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  int ir ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Average read length (nt)", 1, 0, 0} ,
    { "Compute", "Average fragment multiplicity\tMaximal fragment multiplicity\tMillion distinct sequences\tMillion raw reads\tMegabases sequenced", 2, 0, 0} ,
    { "Compute", "%A\t%T\t%G\t%C\t%N\t%GC", 3, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Insert sizes measured in pairs"
    ;
  if (rc == (void *) 1)
    return  snpChapterCaption (snp, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (snp->ao, "\t%s", ti->title) ;
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 1: /* read length */
	      {
		double mbp = 0, mr = 0, mbp1 = 0, mr1 = 0 ; ;

		aceOutf (snp->ao, "\t") ;
		tt = ac_tag_table (rc->snp, "Letter_profile", h) ;
		if (tt)
		  {
		    int ir, nf = 0, n1 ; float n2 ;
		    char buf[64] ; buf[0] = 0 ;
		    for (ir = 0 ; ir < tt->rows  ; ir++)
		      {
			ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
			if (ir == 0) 
			  strncpy (buf, ccp, 64) ;
			if (strncmp (buf, ccp, 64))
			  {
			    mbp1 += mbp ; mr1 += mr ;
			    if (mr)
			      {
				int uk =  .005 + mbp/mr ;
				float uf =  mbp/mr ;
				int u1 = 100 * uf + 0.5 ;
				if (100 * uk == u1)
				  aceOutf (snp->ao, "%s%s:%d", nf ? ", " : EMPTY, buf, uk) ;
				else
				  aceOutf (snp->ao, "%s%s:%.2f", nf ? ", " : EMPTY, buf, uf) ;
			      }
			    else
			      aceOutf (snp->ao, "%s%s:0",  nf ? ", " : EMPTY, buf) ;
			    nf++ ; mbp = mr = 0 ;
			    strncpy (buf, ccp, 64) ;
			  }
			n1 = ac_table_int (tt, ir, 1, 0) ;
			n2 = ac_table_float (tt, ir, 7, 0) ; 
			mr += n2 ;
			mbp += n1 * n2 ;
			if (0) fprintf(stderr, "n1=%d  n2=%f mr=%f\n",n1,n2, mr); 
		      }
		    mbp1 += mbp ; mr1 += mr ;
		    if (tt && tt->rows > 999)
		      {
			mr = rc->var[T_Read] ;
			mbp = rc->var[T_Length] ;
		      }
		    if (mr)
		      {
			int uk =  .005 + mbp/mr ;
			float uf =  mbp/mr ;
			int u1 = 100 * uf + 0.5 ;
			if (100 * uk == u1)
			  aceOutf (snp->ao, "%s%s:%d", nf ? ", " : EMPTY, buf, uk) ;
			else
			  aceOutf (snp->ao, "%s%s:%.2f", nf ? ", " : EMPTY, buf, uf) ;
		      }
		    else
		      aceOutf (snp->ao, "%s%s:0",  nf ? ", " : EMPTY, buf) ;
		  }
		else
		  aceOutf (snp->ao, "%d", (int)(rc->var[T_Length] + .5)) ;
		break ;
	      }
	    case 2:
	      aceOutf (snp->ao, "\t%.2f", rc->var[T_Seq] ? rc->var[T_Read]/rc->var[T_Seq] : 0) ;
	      aceOutf (snp->ao, "\t%d", ac_tag_int (rc->snp, "Maximal_read_multiplicity", 0)) ;
	      snpShowMillions (snp, rc->var[T_Seq]) ; 
	      snpShowMillions (snp, rc->var[T_Read]) ;
	      aceOutf (snp->ao, "\t%.3f", rc->var[T_kb]/1000) ;   
	      break ;
	    case 3: /* ATGC */
	      {
		float z, za, zt, zg, zc, zn ;
		tt = ac_tag_table (rc->snp, "ATGC_kb", h) ;
		z = za = zt = zg = zc = zn = 0 ;
		for (ir = 0 ; tt && ir < tt->rows && ir < 1 ; ir++)
		  {
		    za = ac_table_float (tt, ir, 5, 0) ; 
		    zt = ac_table_float (tt, ir, 6, 0) ; 
		    zg = ac_table_float (tt, ir, 7, 0) ; 
		    zc = ac_table_float (tt, ir, 8, 0) ; 
		    zn = ac_table_float (tt, ir, 9, 0) ; 
		  } 
		z = za + zt + zg + zc + zn ; if (z == 0) z = 1 ;
		aceOutf (snp->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			 , 100 * za/z, 100 *zt/z, 100 *zg/z, 100 * zc/z, 100 *zn/z
			 , 100 * (zg + zc) /z
			 ) ;
	      }
	    }
	}
      else
	snpShowTag (snp, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpBeforeAli */

/*************************************************************************************/
/* Data characteristics before alignment */
static void snpAli (TSNP *snp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  const char *ccp ;
  float z, zb, zc ;
  int ir ;
  TT *ti, tts[] = { 
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Compute", "Million raw reads", 20, 0, 0} ,  /* copied from snpBeforeAli */
    { "Compute", "Million reads aligned on any target by Magic\tAverage length aligned per read (nt)\tMb aligned on any target\t% Mb aligned on any target before clipping\t% length aligned on average before clipping", 1, 0, 0} , 
    { "Compute", "% length aligned after clipping adaptors and barcodes (nt)", 2, 0, 0 } ,
    { "Compute", "% Reads aligned on any target", 3, 0, 0 } ,
 
    {  0, 0, 0, 0, 0}
  }; 
  const char *caption =
    "Global alignment statistics"
    ;
  if (rc == (void *) 1)
    return  snpChapterCaption (snp, tts, caption) ;

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (snp->ao, "\t%s", ti->title) ;
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
      else if (! strcmp (ti->tag, "Compute"))
	{
	  switch (ti->col)
	    {
	    case 999:
	      z =  rc->var[T_G] +rc->var[T_C] ;
	      if (z > 1000) z = 1000 ;
	      aceOutf (snp->ao, "\t%.1f", z/10) ;
	      break ;
	    case 20: /* Million raw reads */
	      snpShowMillions (snp, rc->var[T_Read]) ;
	      break ;
	    case 1:
	      ccp = EMPTY ; tt = ac_tag_table (rc->snp, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
 		{
 		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
 		  if (! strcasecmp (ccp, "any"))
 		    {
 		      z = ac_table_float (tt, ir, 3, 0) ;
 		      zb = ac_table_float (tt, ir, 5, 0) ;
 		      zc = ac_table_float (tt, ir, 7, 0) ;
 		    }
 		} 

	      aceOutf (snp->ao, "\t%.3f", z/1000000) ;
	      aceOutf (snp->ao, "\t%.2f", zc) ;  /* average length aligned */
	      aceOutf (snp->ao, "\t%.3f", zb/1000) ; /* Mb aligned on any target */

	      aceOutf (snp->ao, "\t%.2f", rc->var[T_kb] ? 100 *zb / rc->var[T_kb] : 0) ; /* % Mb aligned on any target before clipping */
	      {
		float z1 = rc->var[T_kb], z2 = 0 ;
		AC_TABLE tt2 = ac_tag_table (rc->snp, "Unaligned", h) ;
		z2 = tt2 ? ac_table_float (tt2, 0, 4,0) : 0 ; /* kb Unaligned */
		aceOutf (snp->ao, "\t%.2f ", 100*zb/(z1 - z2)) ;  /* % length aligned on average before clipping */
	      }
	   
	 
	      break ;
	    case 2:   /*  Average clipped length */
	      {
		int ir ;
		float z = -1, zClipped = -1 ;
		AC_TABLE tt = ac_tag_table (rc->snp, "nh_Ali", h) ;
		
		for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		  {
		    if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), "any"))
		      {
			z =  ac_table_float (tt, ir, 7, -1) ; 
			zClipped = ac_table_float (tt, ir, 11, -1) ;
			break ;
		      }
		  }
		aceOutf (snp->ao, "\t") ;
		if (z > -1)
		  {
		    aceOutPercent (snp->ao, 100.00 * z/zClipped) ;
		  } 
		  
	      }
	      break ;

	    case 3:
	      ccp = EMPTY ; tt = ac_tag_table (rc->snp, "nh_Ali", h) ;
 	      z = zb = zc = 0 ;
 	      for (ir = 0 ; tt && ir < tt->rows ; ir++)
 		{
 		  ccp = ac_table_printable (tt, ir, 0, EMPTY) ;
 		  if (! strcasecmp (ccp, "any"))
 		    {
 		      z = ac_table_float (tt, ir, 3, 0) ;
 		    }
 		} 

	      aceOutf (snp->ao, "\t%.2f", rc->var[T_Read] ? 100 *z / rc->var[T_Read] : 0) ;
	 
	      break ;
	    }
	}
      else
	snpShowTag (snp, rc, ti) ;
    }

  ac_free (h) ;
  return;
}  /* snpAli */

/*************************************************************************************/

static void snpAvLengthAli (TSNP *snp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tt ;
  int ir, irAny = -1 ;
  float z, avClipped = 1 ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "TITLE", "Title", 10, 0, 0} ,
    { "Jump5", "Bases clipped 5prime in read 1, so that a maximal number of reads alig from base 1", 0, 0, 0} ,
    { "Jump5", "Bases clipped 5prime in read 2, so that a maximal number of reads alig from base 1", 1, 0, 0} ,
    { "Compute", "Average length after clipping adaptors and barcodes (nt)", 1, 0, 0 } ,
    
    { "Compute", "any:Average length aligned (nt)", 2, 0, 0 } ,
    { "Compute", "any1:Average length aligned in read 1 (nt)", 2, 0, 0 } ,
    { "Compute", "any2:Average length aligned in read 2 (paired end)", 2, 0, 0 } ,

    /*  { "Compute", "RT_:Average length aligned on reference transcriptome (nt)", 2, 0, 0 } , :"reference transcriptome used in project like RefSeq or Flybase or Encode (37.104 models) */
    { "Compute", "KT_RefSeq:Average length aligned on RefSeq (nt)", 2, 0, 0 } , /*:"Encode.37.70 Jan2013" */
    { "Compute", "MT_EBI:Average length aligned on EBI (nt)", 2, 0, 0 } , /*:"Encode.37.70 Jan2013" */
    { "Compute", "ET_av:Average length aligned on AceView (nt)", 2, 0, 0 } , /*: AceView" */
    { "Compute", "LT_magic:Average length aligned on magic (nt)", 2, 0, 0 } , /*: Future AceView" */
    { "Compute", "NT_:Average length aligned on Other (nt)", 2, 0, 0 } , /*: other unnoficial  annotation */

    { "Compute", "Z_genome:Average length aligned on Genome (nt)", 2, 0, 0 } ,
    { "Compute", "A_mito:Average length aligned on Mitochondria (nt)", 2, 0, 0 } ,
    { "Compute", "B_rRNA:Average length aligned on rRNA (nt)", 2, 0, 0 } ,

    { "Compute", "C_chloro:Average length aligned on Chloroplast (nt)", 2, 0, 0 } ,
    { "Compute", "QT_smallRNA:Average length aligned on small RNAs (nt)", 2, 0, 0 } ,
    { "Compute", "1_DNASpikeIn:Average length aligned on DNA spikeIn (nt)", 2, 0, 0 } ,
    { "Compute", "0_SpikeIn:Average length aligned on RNA spikeIn (nt)", 2, 0, 0 } ,
    { "Compute", "z_gdecoy:Average length aligned on Imaginary genome specificity control (nt)", 2, 0, 0 } ,
    

    { "Spacer", "", 0, 0, 0} ,

    { "Compute", "any:Average % length aligned (nt)", 2, 0, F_p10 } ,
    { "Compute", "any1:Average % length aligned in read 1 (nt)", 2, 0, F_p10 } ,
    { "Compute", "any2:Average % length aligned in read 2 (paired end)", 2, 0, F_p10 } ,

    { "Compute", "KT_RefSeq:Average % length aligned on RefSeq", 2, 0, F_p10} , /*:"Encode.37.70 Jan2013" */
    { "Compute", "MT_EBI:Average % length aligned on EBI", 2, 0, F_p10} , /*:"Encode.37.70 Jan2013" */
    { "Compute", "ET_av:Average % length aligned on AceView", 2, 0, F_p10} , /*: AceView" */
    { "Compute", "LT_magic:Average % length aligned on magic", 2, 0, F_p10} , /*: Future AceView" */
    { "Compute", "NT_:Average % length aligned on Other", 2, 0, F_p10} , /*: other unnoficial  annotation */

    { "Compute", "Z_genome:Average % length aligned on Genome", 2, 0, F_p10} ,
    { "Compute", "A_mito:Average % length aligned on Mitochondria", 2, 0, F_p10} ,
    { "Compute", "B_rRNA:Average % length aligned on rRNA", 2, 0, F_p10} ,

    { "Compute", "C_chloro:Average % length aligned on Chloroplast", 2, 0, F_p10} ,
    { "Compute", "QT_smallRNA:Average % length aligned on small RNAs", 2, 0, F_p10} ,
    { "Compute", "1_DNASpikeIn:Average % length aligned on DNA spikeIn", 2, 0, F_p10} ,
    { "Compute", "0_SpikeIn:Average % length aligned on RNA spikeIn", 2, 0, F_p10} ,
    { "Compute", "z_gdecoy:Average % length aligned on Imaginary genome specificity control", 2, 0, F_p10} ,
   
    {  0, 0, 0, 0, 0}
  } ; 
  
  const char *caption =
    hprintf (h, "Average length (nt) aligned per target in  %s"
	     "Average length aligned per target in nucleotide (table to the left) and in percentage of the length of the read, after clipping adaptors and barcodes (table to the right)." 
	     "Average lengths are reported only if more than 1000 reads are aligned. "
	     "This is seldom the case for the imaginary decoy genome, used as a mapping specificity control, yet the histogram of length aligned on this target is used to set the thresholds for filtering alignments of low quality." 
	     "Ribosomal RNA (rRNA) is encoded in the genome in many copies with slight sequence variations, hence when we align the archetypic rRNA precursor, usually from RefSeq or GenBank, the quality of the alignments is not expected to be perfect."  
	     "Similarly, the small RNA targets, including tRNA, miRNA and other small RNAs often delineates only the mature product, and the actual genes are larger than the gene models. Hence usually only a fraction of the read lengths align on these targets. This is also true for the imperfect but useful human RNA genes from UCSC."
	     
	     , snp->project
	     ) ;
  if (rc == (void *) 1)
    {
      h = ac_new_handle () ;
      snpChapterCaption (snp, tts, caption) ;
      ac_free (h) ;
      return ;
    }
  
  for (ti = tts ; ti->tag ; ti++)
    {
      int nw = 0 ;
      char buf[256] ;

      if (rc == 0)
	{
	  const char *ccp = strchr (ti->title, ':') ;
	  if (!ccp || ti->col != 2)
	    ccp = ti->title ;
	  else 
	    ccp++ ;
	  aceOutf (snp->ao, "\t%s", ccp) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}
 
      tt = ac_tag_table (rc->snp, "nh_Ali", h) ;

      if (! strcmp (ti->tag, "Compute"))
	{
	  aceOutf (snp->ao, "\t") ;
	  switch (ti->col)
	    {
	    case 1:   /*  Average clipped length */
	      irAny = ir ;
	      for (ir = 0 ; tt &&  ir < tt->rows; ir++)
		{
		  if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), "any"))
		    irAny = ir ;
		}

	      avClipped = z = irAny > 0 ? ac_table_float (tt, irAny, 11, -1) : -1 ; 
	      if (z > -1)
		{
		  nw++ ;
		  if (ti->format == F_p10)
		    aceOutPercent (snp->ao, 100.00 * z/avClipped) ;
		  else
		    aceOutf (snp->ao, "%.2f", z) ;
		}
	      break ;
	      
	    case 2:   /*  Average aligned length */
	      strncpy (buf, ti->title, 255) ;
	      { char *cp = strchr (buf, ':') ;
		if (cp) *cp = 0 ;
	      } 
	      for (ir = 0 ; tt && ir < tt->rows; ir++)
		{
		  if (! strcasecmp (ac_table_printable (tt, ir, 0, ""), buf))
		    {
		      z =  ac_table_float (tt, ir, 3, -1) ; 
		      avClipped = ac_table_float (tt, ir, 11, -1) ;
		      if (z > 1000) /* at least 1000 tags aligned */
			{
			  z = irAny > 0 ? ac_table_float (tt, ir, 7, -1) : -1 ; 
			  if (z > -1)
			    {
			      nw++ ;
			      if (ti->format == F_p10)
				aceOutPercent (snp->ao, 100.00 * z/avClipped) ;
			      else
				aceOutf (snp->ao, "%.2f", z) ;
			    }
			}
		      break ;
		    }
		}
	      break ;
	    }
	  if (! nw) aceOut (snp->ao, EMPTY) ;
	}
      else
	snpShowTag (snp, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
}  /* snpAvLengthAli */
  
/*************************************************************************************/
/*************************************************************************************/

static void snpMismatchTypes (TSNP *snp, RC *rc)
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
  if (rc == (void *) 1)
    return  snpChapterCaption (snp, tts, caption) ;


  memset (zSub, 0, sizeof (zSub)) ;
  /*
Distribution of mismatches in best unique alignments, absolute, observed, counts
*/

  for (ti = tts ; ti->tag ; ti++)
    {
      char buf[256] ;

      if (rc == 0)
	{
	  aceOutf (snp->ao, "\t%s", ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}

      tt = ac_tag_table (rc->snp, "Error_profile", h) ;

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
	      aceOutf (snp->ao, "\t%.0f", zz) ;
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
	      aceOutf (snp->ao, "\t%.0f", zAny) ;
	      aceOutf (snp->ao, "\t%.5f", denominator > 0 ? zAny / (1000 * denominator) : 0) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zTransition, zTransversion) ;

	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion, zDeletion ) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zSlidingInsertion, zSlidingDeletion) ;

	      { 
		int i ; 
		for (i = 0 ; i < 12 ; i++)
		  aceOutf (snp->ao, "\t%.0f", zSub[i]) ;
	      }

	      aceOutf (snp->ao, "\t%.0f", zInsertionA) ;
	      aceOutf (snp->ao, "\t%.0f", zInsertionT) ;
	      aceOutf (snp->ao, "\t%.0f", zInsertionG) ;
	      aceOutf (snp->ao, "\t%.0f", zInsertionC) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionA) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionT) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionG) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionC) ;

	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion1, zDeletion1) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion2, zDeletion2) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion3, zDeletion3) ;

	      aceOutf (snp->ao, "\t%.5f", denominator > 0 ? zTransition / (1000 * denominator) : 0) ;
	      aceOutf (snp->ao, "\t%.5f", denominator > 0 ? zTransversion / (1000 * denominator) : 0) ;
	      aceOutf (snp->ao, "\t%.5f", denominator > 0 ? (zInsertion + zSlidingInsertion)/ (1000 * denominator) : 0) ;
	      aceOutf (snp->ao, "\t%.5f", denominator > 0 ? (zDeletion + zSlidingDeletion)/ (1000 * denominator) : 0) ;

	      if (zAny)
		{
		  aceOutf (snp->ao, "\t%.2f", 100 * zTransition / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zTransversion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * (zInsertion + zDeletion + zSlidingInsertion + zSlidingDeletion) / zAny) ;

		  aceOutf (snp->ao, "\t%.2f", 100 * zInsertion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zDeletion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zSlidingInsertion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zSlidingDeletion / zAny) ; 
		}
	      else
		aceOut (snp->ao, "\t\t\t\t\t\t\t") ;
	      break ;
	    default: /* alread  treated in case 2 */
	      break ;

	    }
	}
      else
	snpShowTag (snp, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* snpMismatchTypes */

/*************************************************************************************/

static void snpSnpTypesDo (TSNP *snp, RC *rc, BOOL isRejected)
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
  if (rc == (void *) 1)
    return  snpChapterCaption (snp, tts, caption) ;


  memset (zSub, 0, sizeof (zSub)) ;
  /*
Distribution of mismatches in best unique alignments, abslute, oberved, counts
*/

  for (ti = tts ; ti->tag ; ti++)
    {
      char buf[256] ;

      if (rc == 0)
	{
	  aceOutf (snp->ao, "\t%s %s", isRejected ? "Rejected " : "",  ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
	}

      tt = ac_tag_table (rc->snp, "SNP_profile", h) ;

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
	      aceOutf (snp->ao, "\t%.0f", zz) ;
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
	      aceOutf (snp->ao, "\t%.0f", zAny) ;
	      aceOutf (snp->ao, "\t%.5f", denominator > 0 ? zAny / (1000 * denominator) : 0) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zTransition, zTransversion) ;

	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion, zDeletion ) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zSlidingInsertion, zSlidingDeletion) ;

	      { 
		int i ; 
		for (i = 0 ; i < 12 ; i++)
		  aceOutf (snp->ao, "\t%.0f", zSub[i]) ;
	      }

	      aceOutf (snp->ao, "\t%.0f", zInsertionA) ;
	      aceOutf (snp->ao, "\t%.0f", zInsertionT) ;
	      aceOutf (snp->ao, "\t%.0f", zInsertionG) ;
	      aceOutf (snp->ao, "\t%.0f", zInsertionC) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionA) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionT) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionG) ;
	      aceOutf (snp->ao, "\t%.0f", zDeletionC) ;

	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion1, zDeletion1) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion2, zDeletion2) ;
	      aceOutf (snp->ao, "\t%.0f\t%.0f", zInsertion3, zDeletion3) ;


	      if (zAny)
		{
		  aceOutf (snp->ao, "\t%.2f", 100 * zTransition / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zTransversion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * (zInsertion + zDeletion + zSlidingInsertion + zSlidingDeletion) / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zInsertion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zDeletion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zSlidingInsertion / zAny) ;
		  aceOutf (snp->ao, "\t%.2f", 100 * zSlidingDeletion / zAny) ; 
		  if ( zSub[1] > 100)
		    aceOutf (snp->ao, "\t%.2f", zSub[0]/zSub[1] -1.0) ;
		  else
		    aceOut (snp->ao, "\t") ;
		}
	      else
		aceOut (snp->ao, "\t\t\t\t\t\t\t\t") ;
	      break ;
	    default: /* alread  treated in case 2 */
	      break ;

	    }
	}
      else
	snpShowTag (snp, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /* snpSnpTypesDo */

/***********/

static void snpSnpTypes (TSNP *snp, RC *rc)
{ 
  snpSnpTypesDo (snp, rc, FALSE) ;
}  /* snpSnpTypes */

/***********/

static void snpSnpRejectedTypes (TSNP *snp, RC *rc)
{ 
  snpSnpTypesDo (snp, rc, TRUE) ;
}  /* snpSnpRejectedTypes */

/*************************************************************************************/

static void snpSnpCoding (TSNP *snp, RC *rc)
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
  if (rc == (void *) 1)
    return  snpChapterCaption (snp, tts, caption) ;
  if (rc)
    {
      ttS = ac_tag_table (rc->snp, "SNP", h) ;
      ttG = ac_tag_table (rc->snp, "Genomic", h) ;
      ttP = ac_tag_table (rc->snp, "Protein_changing", h) ;
    }

  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (snp->ao, "\t%s", ti->title) ;
	  continue ; 
	}
      else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
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
	      nT = ac_tag_int (rc->snp, "Tested_sites", 0) ;
	      nnM = ac_tag_int (rc->snp, "Not_measurable_sites", 0) ; 
	      nR = ac_tag_int (rc->snp, "Rejected_sites", 0) ;
	      if (nT > 0)
		{
		  nC = nT - nnM ;
		  z1 = nT ; z1 = 100.0 * nC / z1 ;
		  z2 = nC ; z2 = nC > 0 ? (100.0 * nR / z2) : 0 ;
		  aceOutf (snp->ao, "\t%d\t%.2f\t%d\t%.2f", nT, z1, nC, z2) ; 
		}
	      else
		aceOut (snp->ao, "\t\t\t\t") ;
	      break ;
	    case 2:
	      tt = ttG ;
	      if (tt)
		{
		  aceOutf (snp->ao, "\t%d\t%d\t%d\t%d\t%d\t%d"
			   , ac_table_int (tt, 0, 0, 0)
			   , ac_table_int (tt, 0, 2, 0)
			   , ac_table_int (tt, 0, 4, 0)
			   , ac_table_int (tt, 0, 6, 0) 
			   , ac_table_int (tt, 0, 8, 0)
			   , ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (snp->ao, "\t\t\t\t\t\t") ;
	     
	      break ;
	    case 3: 
	      tt = ttP ;
	      if (tt)
		{
		  aceOutf (snp->ao, "\t%d\t%d\t%d\t%d"
			   , ac_table_int (tt, 0, 0, 0)
			   , ac_table_int (tt, 0, 2, 0)
			   , ac_table_int (tt, 0, 4, 0) +  ac_table_int (tt, 0, 6, 0) + ac_table_int (tt, 0, 8, 0)
			   , ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (snp->ao, "\t\t\t\t") ;
	      break ;

	    case 4:  
	      nC = nT = 0 ;
	      tt = ac_tag_table (rc->snp, "Genomic", h) ;
	      if (tt)
		{
		  nC = ac_table_int (tt, 0, 6, 0) +  ac_table_int (tt, 0, 8, 0) ;
		  nT = ac_table_int (tt, 0, 10, 0) ;
		}
	      aceOut (snp->ao, "\t") ;  
	      if (nC + nT >= 10)
		{
		  z1 = (nC < 100 * nT ?  nC/(1.0 * nT) : 1000) ;
		  aceOutf (snp->ao, "%.2f", z1) ;
		}
		
	      break ;
#ifdef USELESS
	    case 5:
	      tt = ttG ;
	      nT = tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
	      if (nT > 10)
		{
		  float z = 100.0/nT ;
		  aceOutf (snp->ao, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
			   , z * ac_table_int (tt, 0, 2, 0)
			   , z * ac_table_int (tt, 0, 4, 0)
			   , z * ac_table_int (tt, 0, 6, 0) 
			   , z * ac_table_int (tt, 0, 8, 0)
			   , z * ac_table_int (tt, 0, 10, 0)			   
			   ) ;
		}
	      else
		aceOut (snp->ao, "\t\t\t\t\t") ;
	      break ;
	    case 6:
	      tt = ttG ;
	      nT = tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
	      tt = ttP ;
	      if (nT > 10)
		{
		  float z = 100.0/nT ;
		  aceOutf (snp->ao, "\t%.2f\t%.2f\t%.2f\t%.2f"
			   , z * ac_table_int (tt, 0, 0, 0)
			   , z * ac_table_int (tt, 0, 2, 0)
			   , z * (ac_table_int (tt, 0, 4, 0) +  ac_table_int (tt, 0, 6, 0) + ac_table_int (tt, 0, 8, 0))
			   , z * ac_table_int (tt, 0, 10, 0)
			   ) ;
		}
	      else
		aceOut (snp->ao, "\t\t\t\t") ;
	      break ;
#endif
	    }
	}
      else
	snpShowTag (snp, rc, ti) ;
    }
  
  ac_free (h) ;
  return;
} /*  snpSnpCoding */

/*************************************************************************************/

static void snpIdentifiers (TSNP *tsnp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Line", "Line Number", 0, 0, 0} ,
    { "VCF", "VCF Chromosome\tVCF position\tVCF ID\tVCF Reference\tVCF Variant", 1, 0, 0} ,
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
  if (rc == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "Line"))
	{
	  aceOutf (tsnp->ao, "\t%d", ++line) ;
	}
      else if (! strcmp (ti->tag, "VCF"))
	{
	  AC_TABLE tt = ac_tag_table (rc->snp, "VCF", h) ;
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
	  if (ac_has_tag (rc->snp, "Intergenic"))
	    aceOut (tsnp->ao, "Intergenic\t") ;
	  else if (ac_has_tag (rc->snp, "Intronic"))
	    aceOut (tsnp->ao, "Intronic\t") ;
	  else if (ac_has_tag (rc->snp, "Non_coding_transcript"))
	    aceOut (tsnp->ao, "Non_coding_transcript\t") ;
	  else if (ac_has_tag (rc->snp, "UTR_5prime"))
	    aceOut (tsnp->ao, "UTR_5prime\t") ;
	  else if (ac_has_tag (rc->snp, "UTR_3prime"))
	    aceOut (tsnp->ao, "UTR_3prime\t") ;
	  else if (ac_has_tag (rc->snp, "Synonymous"))
	    aceOutf (tsnp->ao, "Coding Synonymous\t%s", ac_tag_text (rc->snp, "Synonymous", "")) ;
	  else if (ac_has_tag (rc->snp, "AA_substitution"))
	    aceOutf (tsnp->ao, "Coding substitution\t%s", ac_tag_text (rc->snp, "AA_substitution", "")) ;
	  else if (ac_has_tag (rc->snp, "Length_variation"))
	    {
	      AC_TABLE tt = ac_tag_table (rc->snp, "Length_variation", h) ;
	      aceOutf (tsnp->ao, "%s\t%s %s"
		       , ac_table_printable (tt, 0, 0, "")
		       , ac_table_printable (tt, 0, 1, "")
		       , ac_table_printable (tt, 0, 2, "")
		       ) ;
	    }

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
	  tt1 = ac_tag_table (rc->snp, tag1, h) ;
	  tt2 = ac_tag_table (rc->snp, tag2, h) ;
	  aceOut (tsnp->ao, "\t") ;
	  if (tt1 && tt2 && tt1->cols > ic && tt2->cols > ic)
	    aceOutf (tsnp->ao, "%s > %s"
		     , ac_table_printable (tt1, 0, ic, "")
		     , ac_table_printable (tt2, 0, ic, "")
		     ) ;
	}
      else  /* gName rName pName Dan_Li Typ */
	snpShowTag (tsnp, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpIdentifiers */

/*************************************************************************************/

static void snpDanLiCounts (TSNP *tsnp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "AGLR1 m:m+w", 10, 0, 0} ,
    { "Compute", "AGLR2 m:m+w", 11, 0, 0} ,
    { "Compute", "ROCR1 m:m+w", 12, 0, 0} ,
    { "Compute", "ROCR2 m:m+w", 13, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  const char *caption =
    "Dan Li counts"
    ;
  if (rc == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
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
	  tt = ac_tag_table (rc->snp, tag, h) ;

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
		  }
	    }
	}
      else
	snpShowTag (tsnp, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpDanLiCounts */

/*************************************************************************************/

static void snpDanLiFrequency (TSNP *tsnp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    { "Compute", "AGLR1 m/m+w", 1, 0, 0} ,
    { "Compute", "AGLR2 m/m+w", 2, 0, 0} ,
    { "Compute", "ROCR1 m/m+w", 3, 0, 0} ,
    { "Compute", "ROCR2 m/m+w", 4, 0, 0} ,
    { "Inter", "Conflict", 0, 0, 0} ,
    {  0, 0, 0, 0, 0}
  } ; 

  float alpha = 1.5 ;
  float dqq[1024], qq[1024] ;
  int iti, qqi[1024], qqc[1024], qqm[1024] ;

  memset (qq, 0, sizeof(qq)) ;
  memset (dqq, 0, sizeof(dqq)) ;
  memset (qqi, 0, sizeof(qqi)) ;
  memset (qqc, 0, sizeof(qqc)) ;
  memset (qqm, 0, sizeof(qqm)) ;


  const char *caption =
    "Dan Li counts"
    ;
  if (rc == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (iti = 0, ti = tts ; ti->tag ; iti++, ti++)
    {
      if (rc == 0)
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
	  tt = ac_tag_table (rc->snp, tag, h) ;

	  aceOut (tsnp->ao, "\t") ;
	  if (tt)
	    {
	      int ir ;
	      for (ir = 0 ; ir < tt->rows ; ir++)
		if (!strcmp (ac_table_printable (tt, ir, 0, "toto"), run))
		  {
		    int m = ac_table_int (tt, ir, 1, 0) ;
		    int c = ac_table_int (tt, ir, 3, 0) ;
		    if (c > 0)
		      {
			qq[ti->col] = 1.0*m/c ;
			dqq[ti->col] = 1.0/sqrt(c) ;
			qqi[ti->col] = iti ;
			qqc[ti->col] = c ;
			qqm[ti->col] = m ;

			aceOutf (tsnp->ao, "%.2f", 100.0 * m/c) ;
		      }
		  }
	    }
	}
      else if (! strcmp (ti->tag, "Inter"))
	{
	  int kk = 0 ;
	  int ii = ti->col ;
          float ccc = 0, mmm = 0, ppp = 0, dppp ;
	  aceOut (tsnp->ao, "\t") ;
	  
	  for (ii = 1 ; ii < 5 ; ii++)
	    if (qqi[ii])
	      {
		kk++ ;
		ccc += qqc[ii] ;
		mmm += qqm[ii] ;
	      }
	  if (kk >= 4) /* at least 4 measures */
	    {
	      ppp = mmm/ccc ;
	      dppp = 1.0/sqrt (ccc) ;
	      for (ii = 1 ; ii < 5 ; ii++)
		if (qqi[ii])
		  {
		    float z = qq[ii] - ppp ;
		    if (z < 0) z = -z ;
		    if (z > alpha * (dppp + dqq[ii]))
		      {
			char buf[15] ;
			memcpy (buf, tts[qqi[ii]].title, 5) ;
			buf[5] = 0 ;
			aceOutf (tsnp->ao, "YBAD %s:%.1f:%.1f ", buf,100*qq[ii],100*ppp) ;		    
		      }
		  }
	    }
	}
      else
	snpShowTag (tsnp, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpDanLiFrequency */

/*************************************************************************************/

static void snpBrsFrequencyCounts (TSNP *tsnp, RC *rc, BOOL isFrequency, BOOL isGroup)
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
  AC_TABLE tt = 0 ;

  if (rc == (void *) 1)
    return  snpChapterCaption (tsnp, tts, caption) ;
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	{
	  aceOutf (tsnp->ao, "\t%s", ti->title) ;
	}
      else if (! strcmp (ti->tag, "Runs"))
	{
	  const char *tag ;
	  int ir, jr ;
	  tag = "BRS_counts" ;
	  if (! tt)
	    tt = ac_tag_table (rc->snp, tag, h) ;

	  aceOut (tsnp->ao, "\t") ;
	  for (ir = 0 ; ir < tsnp->truns->rows ; ir++)
	    {
	      const char *runNam = ac_table_printable (tsnp->truns, ir, 0, "toto") ; 
	      int m = 0, c = 0 ;
	      
	      aceOut (tsnp->ao, "\t") ;
	      if (tt)
		{
		  for (jr = 0 ; jr < tt->rows ; jr++)
		    if (!strcmp (ac_table_printable (tt, jr, 0, "toto"), runNam))
		      {
			c = ac_table_int (tt, jr, 1, 0) ;
			m = ac_table_int (tt, jr, 2, 0) ;
			break ;
		      }
		}	      
	      if (isFrequency)
		{
		  if (c >= 10)
		    aceOutf (tsnp->ao, "%.2f", 100.0*m/c) ;
		  else
		    aceOutf (tsnp->ao, "-10") ;
		}
	      else
		aceOutf (tsnp->ao, "%d:%d", m,c) ;
	    }
	  aceOut (tsnp->ao, "\t") ;
	}
      else if (! strcmp (ti->tag, "Groups"))
	{
	  const char *tag ;
	  int ir, jr ;
	  tag = "BRS_counts" ;
	  if (! tt)
	    tt = ac_tag_table (rc->snp, tag, h) ;

	  for (ir = 0 ; ir < tsnp->groups->rows ; ir++)
	    {
	      const char *runNam = ac_table_printable (tsnp->groups, ir, 0, "toto") ; 
	      int m = 0, c = 0 ;
	      
	      aceOut (tsnp->ao, "\t") ;
	      if (tt)
		{
		  for (jr = 0 ; jr < tt->rows ; jr++)
		    if (!strcmp (ac_table_printable (tt, jr, 0, "toto"), runNam))
		      {
			c = ac_table_int (tt, jr, 1, 0) ;
			m = ac_table_int (tt, jr, 2, 0) ;
			break ;
		      }   
		}
	      if (isFrequency)
		{
		  if (c >= 10)
		    aceOutf (tsnp->ao, "%.2f", 100.0*m/c) ;
		  else
		    aceOutf (tsnp->ao, "-10") ;
		}
	      else
		aceOutf (tsnp->ao, "%d:%d", m,c) ;
	    }
	  aceOut (tsnp->ao, "\t") ;
	}
      else
	snpShowTag (tsnp, rc, ti) ;
    }
  ac_free (h) ;
  return;
}  /* snpBrsFrequencyCounts */

static void snpBrsRunFrequency (TSNP *tsnp, RC *rc)
{
  return snpBrsFrequencyCounts (tsnp, rc, TRUE, FALSE) ;
}
static void snpBrsRunCounts (TSNP *tsnp, RC *rc)
{
  return snpBrsFrequencyCounts (tsnp, rc, FALSE, FALSE) ;
}
static void snpBrsGroupFrequency (TSNP *tsnp, RC *rc)
{
  return snpBrsFrequencyCounts (tsnp, rc, TRUE, TRUE) ;
}
static void snpBrsGroupCounts (TSNP *tsnp, RC *rc)
{
  return snpBrsFrequencyCounts (tsnp, rc, FALSE, TRUE) ;
}


/*************************************************************************************/
/* template for a future chapter */
static void snpOther (TSNP *snp, RC *rc)
{
  AC_HANDLE h = ac_new_handle () ;
  TT *ti, tts[] = {
    { "Spacer", "", 0, 0, 0} ,
    {  0, 0, 0, 0, 0}
  }; 

  const char *caption =
    "Statistics before alignment"
    ;
  if (rc == (void *) 1)
    return  snpChapterCaption (snp, tts, caption) ;
  
  for (ti = tts ; ti->tag ; ti++)
    {
      if (rc == 0)
	aceOutf (snp->ao, "\t%s", ti->title) ;
       else if (! strcmp (ti->tag, "TITLE"))
	{
	  const char *ccp = ac_tag_printable (rc->run, "Title", 0) ;
	  aceOutf (snp->ao, "\t%s", ccp ? ccp : ac_name(rc->run)) ;
	  continue ;
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
	snpShowTag (snp, rc, ti) ;
    }

   ac_free (h) ;
  return;
}  /* snpOther */

/*************************************************************************************/

static BOOL snpFilter (TSNP *tsnp, RC *rc)
{
  BOOL ok = TRUE ;
  int minF = tsnp->minSnpFrequency ;
  int minC = tsnp->minSnpCover ;


  if (minF + minC)
    {
      AC_HANDLE h = ac_new_handle () ;
      AC_TABLE tt = ac_tag_table (rc->snp, "BRS_counts", h) ;
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

static const char *allMethods = "Idr" ;

static MM methods [] = {
  {'I', &snpIdentifiers} ,
  {'d', &snpDanLiFrequency} ,
  {'D', &snpDanLiCounts} ,
  {'r', &snpBrsRunFrequency} , 
  {'R', &snpBrsRunCounts} ,
  {'g', &snpBrsGroupFrequency} , 
  {'G', &snpBrsGroupCounts} ,
{ 0, 0 }
} ;

/*************************************************************************************/

static BOOL snpExportSnp (TSNP *tsnp, int iSnp, int type)
{       
  const char *ccp, *lineName ;
  MM *mm ;
  RC *rc = 0 ;
  RC rc0 ;
  
  memset (&rc0, 0, sizeof(RC)) ;

  switch (type)
    {
    case 0: 
      lineName =  "### Caption" ; 
      rc = (void *) 1 ;
      break ;
    case 1: 
      lineName =  "### MAGIC SNP" ;
      rc = (void *) 0 ;
      break ;
    case 2: 
      rc = &rc0 ;
      rc->h = ac_new_handle () ;
      rc->snp = ac_table_obj (tsnp->snps, iSnp, 0, rc->h) ;
      lineName = ac_table_printable (tsnp->snps, iSnp, 0, "-") ;

      break ; 
    }
   
  if (type < 2 || snpFilter (tsnp, rc))
    {
      aceOutf (tsnp->ao, "%s", lineName) ;
      
      ccp = tsnp->export - 1 ;
      while (*++ccp)
	{
	  mm = methods - 1 ;
	  while (mm++, mm->cc)
	    if (mm->cc == *ccp)
	      mm->f(tsnp, rc) ;
	}
      aceOutf (tsnp->ao, "\n" ) ;
    }
  ac_free (rc0.h) ;

  return TRUE ;
} /* snpExportSnp */

/*************************************************************************************/
/*************************************************************************************/

static int snpGetSnpsRunsGroups (TSNP *tsnp)
{
  int ns = 0 ;
  const char *errors = 0 ;
  char *qq ;
  AC_HANDLE h = ac_new_handle () ;
  tsnp->runs = 0 ;
  DICT *dict ;
  vTXT txt ; 

  if (tsnp->project)
    {
      qq = hprintf (h, "select r, t from p in ?project where p == \"%s\" , r in p->run where r ISA runs, t in r->%s", tsnp->project, tsnp->orderBy) ;
      tsnp->truns = ac_bql_table (tsnp->db
				  , qq
				  , 0
				  , "+2"
				  , &errors
				  , tsnp->h
				  ) ;
      qq = hprintf (h, "select r, t from p in ?project where p == \"%s\", r in p->run where r ISA Groups, t in r->%s", tsnp->project, tsnp->orderBy) ;
      tsnp->groups = ac_bql_table (tsnp->db
				   , qq
				   , 0
				   , "+2"
				   , &errors
				   , tsnp->h
				   ) ;
      txt = tsnp->runsHeader = vtxtHandleCreate (tsnp->h) ;
      dict = tsnp->runsDict = dictHandleCreate (256, tsnp->h) ;
      if (tsnp->truns)
	for (int ir = 0 ; ir < tsnp->truns->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tsnp->truns, ir, 0, "toto") ;
	    dictAdd (dict, ccp, 0) ;
	    vtxtPrintf (txt, "\t%s", ccp) ;
	  }
      vtxtPrint (txt, "\t") ;
      txt = tsnp->groupsHeader = vtxtHandleCreate (tsnp->h) ;
      dict = tsnp->groupsDict = dictHandleCreate (256, tsnp->h) ;
      if (tsnp->groups)
	for (int ir = 0 ; ir < tsnp->groups->rows ; ir++)
	  {
	    const char *ccp = ac_table_printable (tsnp->groups, ir, 0, "toto") ;
	    dictAdd (dict, ccp, 0) ;
	    vtxtPrintf (txt, "\t%s", ccp) ;
	  }
      vtxtPrint (txt, "\t") ;

      qq = hprintf (h, "select g, r from p in ?project where p == \"%s\", g in p->run where g#union_of, r in g->union_of, p2 in r->project where p2 == p", tsnp->project) ;
      AC_TABLE tt = ac_bql_table (tsnp->db
				  , qq
				  , 0
				  , 0
				  , &errors
				  , h
				  ) ;
      if (tt)
	for (int ir = 0 ; ir < tt->rows ; ir++)
	  {
	    const char *gNam = ac_table_printable (tsnp->groups, ir, 0, "toto") ;
	    const char *rNam = ac_table_printable (tsnp->groups, ir, 1, "toto") ;
	    int g = 0 ;
	    dictFind (tsnp->groupsDict, gNam, &g) ;
	    if (g >= 1)
	      {
		R2G* g2r ;
		int r = 0 ;

		if (!tsnp->r2gs)
		  tsnp->r2gs = arrayHandleCreate (200, R2G, tsnp->h) ;
		g2r = arrayp (tsnp->r2gs, g, R2G) ; 

		if (! dictFind (tsnp->runsDict, rNam, &r))
		  dictFind (tsnp->groupsDict, rNam, &r) ;
		if (r >= 1)
		  {
		    R2G *r2g ; 
		    r2g = arrayp (tsnp->r2gs, r, R2G) ; 
		  }
	      }
	  }
    }

  tsnp->snps = ac_bql_table (tsnp->db
			     , tsnp->snpType == 3 ? "select s from s in ?Variant where s#danli_counts " : "select s from s in ?Variant  "
			     , 0
			     , 0
			     , &errors
			     , tsnp->h
			     ) ;
  
  ns = tsnp->snps ? tsnp->snps->rows : 0  ;
  
  fprintf (stderr, " snpGetSnps got %d snps, %d runs , %d groups in project %s\n"
	   , ns
	   , tsnp->truns ? tsnp->truns->rows : 0 
	   , tsnp->groups ? tsnp->groups->rows : 0 
	   , tsnp->project
	   ) ;
  
  ac_free (h) ;
  return ns ;
} /* snpGetSnpsRunsGroups */

/*************************************************************************************/

static void snpExportHeader (TSNP *snp)
{
  aceOutf (snp->ao, "## Project %s, quality control report, File %s : %s\n"
	   , snp->project
	   , aceOutFileName (snp->ao)
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

  tsnp.minSnpFrequency = 20 ;
  tsnp.minSnpCover = 20 ;
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
  
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


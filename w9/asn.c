/*  Last edited: Dec  4 14:20 1998 (fw) */
/*****************************************************************************
 *
 *   asn.c
 *   	by LaDeana Hillier, WashU, and help from Jim Ostell, NCBI.
 *   original version: 11/15/93
 *   modified version: 9/95 (update for ACEDB 4) Dave Ficenec, WashU.
 *   modified version: 3/96 (update for human) Dave Ficenec, WashU.
 *   modified version: 10/96 (update for human) Suemee Shin, WashU.
 *
 *   This is a modification of the embl.c file (comment below) to produce
 *   an ASN.1 format direct submission for NCBI. It assumes the NCBI
 *   toolkit is installed and that either the environment variable "NCBI" is
 *   set to a directory in which the file ".ncbirc" exists or that ".ncbirc"
 *   is in the current directory.
 *
 *   ".ncbirc" contains a path to /ncbi/data which contains the data files
 *   seqcode.val and gc.val, sequence alphabets and genetic codes, respectively
 *   used by the toolkit. If you have Entrez installed, then this is already
 *   properly configured.
 *
 *   These instructions apply only to UNIX machines.
 *
 ****************************************************************************
 *
 *  File: embl.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *
 *  This file is part of the ACEDB genome database package, written by
 *  Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *  Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 ***************************************************************************/

/* $Id: asn.c,v 1.3 2009/06/22 19:36:29 mieg Exp $ */

/*
WRG 961004 Fixed partial interval problems
WRG 960929 Sequence similarity comments will appear at beginning of /note=
WRG 960929 Added calls to new NCBI subutil functions for WASHU and htgs_3 tags
WRG 960917 Added conditional code related to HUMAN.
*/

#include "acedb.h"

#include "display.h"
#include "bs.h"
#include "classes.h"
#include "tags.h"
#include "systags.h"
#include "dna.h"
#include "lex.h"
#include "gmap.h"
#include "subutil.h"    /* this file includes NCBI stuff */

enum {NOPROJECT, ELEGANS, YEAST, HUMAN, BRIGGSAE, FUGU};

static char   Institution[] = "WUGSC"; /* submittor tag */
static int    gProject;          /* Project id                      */
static int    gAsnDebug;         /* Debgging output off/on.         */
static int    gDay,gMon,gYear;   /* Submit date.                    */
static char   gChromosome[10];   /* Chromosome name                 */
static char   gSequence[128];    /* Current sequence name           */
static char   gRawClone[128];    /* Unadulterated clone name        */
static char   gClone[128];       /* Current clone name              */
static char   gClonePat[128];    /* Current clone+gene name pattern */
static char   gRawClonePat[128];   /* Current clone+gene name pattern */
static char   gCloneType[128];   /* Current clone type              */
static char   gLibrary[128];     /* Library                         */
static char   keyword[10][256];  /* keywords                        */
static char   gMapPosition[128]; /* Map position, xq28              */
static char   gLocus[256];       /* Current locus name              */
static int    gLength;           /* Length of current DNA sequence  */
static char   gBuf[32000];       /* All-purpose char buffer         */
static int    gOffset;           /* Subseq in superseq coords       */
static BOOL   gOffForward;       /* Subseq codirect. with superseq  */
static int    gMiss1,gMiss2;     /* Amount clipped from subseq      */

static NCBISubPtr  nsp;          /* NCBI submission pointer         */
static SeqEntryPtr nucprot_entry;/* NCBI nucleotide & protein entry */
static SeqEntryPtr nuc_entry;    /* NCBI nucleotide entry           */
static SeqEntryPtr prot_entry;   /* NCBI protein entry              */
static SeqEntryPtr the_entry;    /* NCBI entry (one of above)       */

static KEY  _Database;
static KEY  _scRNA;
static KEY  _misc_RNA;
static KEY  _pseudogene;
static KEY  _Feature;
static KEY  _Clone_Library;
static KEY  _EMBL_feature;
static KEY  _EMBL_threshold;
static KEY  _EMBL_qualifier;
static KEY  _EMBL_dump_YES;
static KEY  _EMBL_dump_NO;

char *catnames(char *, char *);
unsigned long get_cDNA_NID(char *);


/*********************************************************************
 * checkFileExists:
 *
 * Check if file in current directory; prompt user to continue
 * if it does; exit at user request.
 ********************************************************************/

static int checkFileExists(char *filename)
{
	struct stat	sbuf;
	char	a[256];

	strcpy(a, "");
	if (stat(filename, &sbuf) == 0) {
		printf("  ASN file %s already exists.  Abort this ASN dump? [n] ",
				filename);
		fgets(a, sizeof(a), stdin);
    }
	return (a[0] == 'Y' || a[0] == 'y');
}

/*********************************************************************
 * sortExons:
 *
 * Swap exon limits if start > end; then sort exons by start.
 ********************************************************************/

static void sortExons(Array a)
{
  int i,j,lo,hi;

  for (i=0; i<arrayMax(a); i+=2) {
    if ( arr(a,i,BSunit).i > arr(a,i+1,BSunit).i ) {
      j = arr(a,i+1,BSunit).i;
      arr(a,i+1,BSunit).i = arr(a,i,BSunit).i;
      arr(a,i,BSunit).i = j;
    }
  }

  for (j=2; j+1<arrayMax(a); j+=2) {
    lo = arr(a,j  ,BSunit).i;
    hi = arr(a,j+1,BSunit).i;
    i=j-2;
    while (i>0 && arr(a,i,BSunit).i > lo) {
      arr(a,i+2,BSunit).i = arr(a,i+0,BSunit).i;
      arr(a,i+3,BSunit).i = arr(a,i+1,BSunit).i;
      i-=2;
    }
    arr(a,i+2,BSunit).i = lo;
    arr(a,i+3,BSunit).i = hi;
  }
}


/*********************************************************************
 * offAdjust:
 *
 * Converts limits to dump sequence coordinates.
 ********************************************************************/

static void offAdjust (int *x1, int *x2)
{
  if (gOffForward) {
    *x1 -= gOffset - 1 ;
    *x2 -= gOffset - 1 ;
  }
  else {
    *x1 = gOffset + 1 - *x1 ;
    *x2 = gOffset + 1 - *x2 ;
  }
}

/*********************************************************************
 * clipAdd:
 *
 * Accumulate counts of nucleotides outside superseq.
 * Modified from EMBL version: does not build exon array.
 ********************************************************************/

static void clipAdd (int x1, int x2)
{
  if (x2 < 1) {
    gMiss1 += x2 - x1 + 1 ;
    return ;
  }
  if (x1 > gLength) {
    gMiss2 += x2 - x1 + 1 ;
    return ;
  }
  if (x1 < 1) {
    gMiss1 += 1 - x1 ;
    x1 = 1 ;
  }
  if (x2 > gLength) {
    gMiss2 += x2 - gLength ;
    x2 = gLength ;
  }
}

/*********************************************************************
 * FixCloneName:
 *
 * Makes character string uppercase and returns a pointer past "H_".
 ********************************************************************/

static char * FixCloneName(char *name)
{
  int i;
  int l = strlen(name);
  char *newname;

  for (i=0; i<l; i++) {
    name[i] = toupper(name[i]);
  }

  if (strncasecmp(name,"H_",2) != 0) {
    newname = name;
  }
  else {
    newname = name+2;
  }

  return newname;
}

/*********************************************************************
 * FixCitationName:
 *
 * Best guess as to first initials and last name.
 * Returns 1 if arguments look reversed (need fixing), 0 otherwise.
 ********************************************************************/

static int FixCitationName(char *initials, char *lastname)
{
  if (strchr(initials,'.') != NULL) {
    return 0;
  }
  else if (strchr(lastname,'.') != NULL) {
    return 1;
  }
  else if ( strlen(initials) > strlen(lastname) ) {
    return 1;
  }
  else {
    return 0;
  }
}

/*********************************************************************
 * AddFinishers:
 *
 * Adds finisher reference from key to NCBI submission.
 ********************************************************************/

static void AddFinishers(OBJ this, NCBISubPtr sp, SeqEntryPtr ep, char *title)
{
  OBJ	  obj;
  PubPtr  pub;
  static  Array a;
  int	  i;
  char    *finisher;
  char    lastname[80];
  char    firstini[80];

  a = arrayReCreate(a, 32, BSunit);

  if (bsFindTag(this, _From_Author) && bsFlatten(this,1,a)) {
    if (gAsnDebug) printf("  Adding finisher reference...\n");
    pub = CitArtBuild(sp, title, NULL, NULL, NULL, NULL,
		      gMon, gDay, gYear, PUB_STATUS_UNPUBLISHED);
    for (i=0; i<arrayMax(a); i++) {
      finisher = name(arr(a,i,BSunit).k);
      sscanf(finisher,"%s %s",lastname,firstini);
      if (gAsnDebug) printf("  Add finisher: %s\n",finisher);
      if ( FixCitationName(firstini,lastname) )
	AddAuthorToPub(sp, pub, firstini, NULL, NULL, lastname, NULL);
      else
	AddAuthorToPub(sp, pub, lastname, NULL, NULL, firstini, NULL);
    }
    AddPubToEntry(sp, ep, pub);
  }
  else {
    printf("  WARNING: No finisher name...skipping reference.\n");
  }

}


/*********************************************************************
 * AddElegansReference:
 *
 * Adds C. elegans reference from OBJ to NCBI submission.
 * PLEASE remove this routine when C. elegans papers
 * become part of the ACeDB database; the general routine
 * AddReferences which is already called will add any reference
 * tags from main sequence entry.
 ********************************************************************/

static void AddElegansReference(NCBISubPtr sp, SeqEntryPtr ep)
{
  PubPtr  pub;

  strcpy(gBuf,"2.2 Mb of contiguous nucleotide sequence from chromosome III of C. elegans");
  pub = CitArtBuild(sp, gBuf, "Nature", "368", NULL, "32-38", 0, 0, 1994, PUB_STATUS_PUBLISHED);

  AddAuthorToPub(sp, pub, "Wilson",NULL,NULL,"R.",NULL);
  AddAuthorToPub(sp, pub, "Ainscough",NULL,NULL,"R.",NULL);
  AddAuthorToPub(sp, pub, "Anderson",NULL,NULL,"K.",NULL);
  AddAuthorToPub(sp, pub, "Baynes",NULL,NULL,"C.",NULL);
  AddAuthorToPub(sp, pub, "Berks",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Bonfield",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Burton",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Connell",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Copsey",NULL,NULL,"T.",NULL);
  AddAuthorToPub(sp, pub, "Cooper",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Coulson",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Craxton",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Dear",NULL,NULL,"S.",NULL);
  AddAuthorToPub(sp, pub, "Du",NULL,NULL,"Z.",NULL);
  AddAuthorToPub(sp, pub, "Durbin",NULL,NULL,"R.",NULL);
  AddAuthorToPub(sp, pub, "Favello",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Fulton",NULL,NULL,"L.",NULL);
  AddAuthorToPub(sp, pub, "Gardner",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Green",NULL,NULL,"P.",NULL);
  AddAuthorToPub(sp, pub, "Hawkins",NULL,NULL,"T.",NULL);
  AddAuthorToPub(sp, pub, "Hillier",NULL,NULL,"L.",NULL);
  AddAuthorToPub(sp, pub, "Jier",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Johnston",NULL,NULL,"L.",NULL);
  AddAuthorToPub(sp, pub, "Jones",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Kershaw",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Kirsten",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Laister",NULL,NULL,"N.",NULL);
  AddAuthorToPub(sp, pub, "Latreille",NULL,NULL,"P.",NULL);
  AddAuthorToPub(sp, pub, "Lightning",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Lloyd",NULL,NULL,"C.",NULL);
  AddAuthorToPub(sp, pub, "McMurray",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Mortimore",NULL,NULL,"B.",NULL);
  AddAuthorToPub(sp, pub, "O'Callaghan",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Parsons",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Percy",NULL,NULL,"C.",NULL);
  AddAuthorToPub(sp, pub, "Rifken",NULL,NULL,"L.",NULL);
  AddAuthorToPub(sp, pub, "Roopra",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Saunders",NULL,NULL,"D.",NULL);
  AddAuthorToPub(sp, pub, "Shownkeen",NULL,NULL,"R.",NULL);
  AddAuthorToPub(sp, pub, "Smaldon",NULL,NULL,"N.",NULL);
  AddAuthorToPub(sp, pub, "Smith",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Sonnhammer",NULL,NULL,"E.",NULL);
  AddAuthorToPub(sp, pub, "Staden",NULL,NULL,"R.",NULL);
  AddAuthorToPub(sp, pub, "Sulston",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Thierry-Mieg",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Thomas",NULL,NULL,"K.",NULL);
  AddAuthorToPub(sp, pub, "Vaudin",NULL,NULL,"M.",NULL);
  AddAuthorToPub(sp, pub, "Vaughan",NULL,NULL,"K.",NULL);
  AddAuthorToPub(sp, pub, "Waterston",NULL,NULL,"R.",NULL);
  AddAuthorToPub(sp, pub, "Watson",NULL,NULL,"A.",NULL);
  AddAuthorToPub(sp, pub, "Weinstock",NULL,NULL,"L.",NULL);
  AddAuthorToPub(sp, pub, "Wilkinson-Sproat",NULL,NULL,"J.",NULL);
  AddAuthorToPub(sp, pub, "Wohldman",NULL,NULL,"P.",NULL);

  AddPubToEntry(sp, ep, pub);
}

/*********************************************************************
 * AddBriggsaeReference:
 *
 * Adds C. briggsae reference from OBJ to NCBI submission.
 * PLEASE remove this routine when C. broggsae papers
 * become part of the ACeDB database; the general routine
 * AddReferences which is already called will add any reference
 * tags from main sequence entry.
 ********************************************************************/

static void AddBriggsaeReference(NCBISubPtr sp, SeqEntryPtr ep)
{
  PubPtr  pub;

  strcpy(gBuf,"The C. briggsae Genome Sequencing Project");
  pub = CitArtBuild(sp, gBuf, NULL, NULL, NULL, NULL, 0, 0, 1996,
		    PUB_STATUS_UNPUBLISHED);
  AddAuthorToPub(sp, pub, "Washington University Genome Sequencing Center",
		 NULL,NULL,NULL,NULL);

  AddPubToEntry(sp, ep, pub);
}

/*********************************************************************
 * AddReferences:
 *
 * Adds reference(s) from OBJ to NCBI submission.
 ********************************************************************/

static void AddReferences(OBJ this, NCBISubPtr sp, SeqEntryPtr ep)
{
  OBJ	  obj;
  BOOL    status;
  KEY     rkey,gpkey;
  KEY     title,journal;
  PubPtr  pub;
  static  Array a;
  int	  i,yr;
  char    *volume,*yt;
  char    *author;
  char    *page1,*page2;
  char    lastname[30];
  char    firstini[30];
  char    journbuf[128];
  char    pagesbuf[128];
  char    volumbuf[128];
  char    titlebuf[512];

  a = arrayReCreate(a, 32, BSunit);

  status = bsGetKey(this, _Reference, &rkey);

  while (status) {

    if (!(obj = bsCreate(rkey))) {
      printf("  Cannot create obj from key %s\n",name(rkey));
      return;
    }

    if (gAsnDebug) printf("  Adding reference %s\n",name(rkey));

    if (bsGetKey(obj, _Title, &title)) {
      strcpy(titlebuf,name(title));
    }
    else {
      strcpy(titlebuf,"");
      printf("  WARNING: Untitled from %s --- assuming unpublished\n",name(rkey));
    }
    if (gAsnDebug) printf("  Title is %s\n",titlebuf);

    if (bsGetKey(obj, _Journal, &journal)) {
      strcpy(journbuf,name(journal));
    }
    else {
      strcpy(journbuf,"");
      printf("  WARNING: No journal for %s --- assuming unpublished.\n",name(rkey));
    }
    if (gAsnDebug) printf("  Journal is %s\n",journbuf);

    if (bsGetData(obj, _Volume, _Text, &volume)) {
      strcpy(volumbuf,volume);
    }
    else {
      strcpy(volumbuf,"");
      printf("  WARNING: No volume for %s --- assuming unpublished.\n",name(rkey));
    }
    if (gAsnDebug) printf("  Volume is %s\n",volumbuf);

    if (!bsGetData(obj, _Year, _Int, &yr)) {
      if (!bsGetData(obj, _Year, _Text, &yt)) {
	printf("  WARNING: Could not get a year from %s --- using today\n",name(rkey));
	yr = gYear;
      }
      sscanf(yt,"%d",&yr);
    }
    if (gAsnDebug) printf("  Year is %d\n",yr);

    if (!bsGetData(obj, _Page, _Text, &page1)) {
      page1 = "";
      printf("  WARNING: No lower page from %s\n",name(rkey));
    }
    if (gAsnDebug) printf("  Page1 is  %s\n",page1);

    if (!bsGetData(obj, _bsRight, _Text, &page2)) {
      page2 = "";
      printf("  WARNING: No lower page from %s\n",name(rkey));
    }
    if (gAsnDebug) printf("  Page2 is  %s\n",page2);

    if (page1 == "" || page2 == "") {
      strcpy(gBuf,"");
    }
    else {
      sprintf(gBuf,"",page1,page2);
    }

    if (strcmp(volumbuf,"") == 0 || strcmp(journbuf,"") == 0) {
      pub = CitArtBuild(sp, titlebuf, NULL, NULL, NULL,
			NULL, 0, 0, yr, PUB_STATUS_UNPUBLISHED);
    }
    else {
      pub = CitArtBuild(sp, titlebuf, journbuf, volumbuf, NULL,
			gBuf, 0, 0, yr, PUB_STATUS_PUBLISHED);
    }

    if (gAsnDebug) printf("  CitArtBuild successful\n");

    if (bsFindTag(obj, _Author) && bsFlatten(obj,1,a)) {
      for (i=0; i<arrayMax(a); i++) {
	author = name(arr(a,i,BSunit).k);
	sscanf(author,"%s %s",lastname,firstini);
	if ( FixCitationName(firstini,lastname) )
	  AddAuthorToPub(sp, pub, firstini, NULL, NULL, lastname, NULL);
	else
	  AddAuthorToPub(sp, pub, lastname, NULL, NULL, firstini, NULL);
	if (gAsnDebug) printf("  Adding author: %s\n",author);
      }
    }

    AddPubToEntry(sp, ep, pub);

    bsDestroy(obj);

    status = bsGetKey(this, _bsDown, &rkey);
  }
}


/*********************************************************************
 * subseqFeatures:
 *
 * Adds subsequence features to NCBI submission.
 ********************************************************************/

static int subseqFeatures (KEY key, int subBeg, int subEnd)
{
  OBJ	obj;                /* general purpose object               */
  BOOL	isPlus;             /* subsequence is on plus strand        */
  Nlm_Boolean	start_not_found, end_not_found;
  Nlm_Boolean	snf, enf;
  int	i;                  /* general purpose int                  */
  int   cdsBeg,cdsEnd;      /* left and right ends of CDS subseq    */
  int   has_cDNA = 0;         /* cdna hit for coding region           */
  int   has_Locus = 0;        /* identified protein for coding region */
  int	has_DBRemark = 0;
  int	has_BriefID = 0;
  KEY   gpkey;              /* general purpose key                  */
  char  *text;              /* pointer to text items                */
  char  *lcp;
  char  trna_codon[32];     /* trna codon (gca, etc.)               */
  char  rnaname[32];        /* rna/trna name                        */
  char  trna_AA[32];        /* trna amino acid (A, etc.)            */
  char  gene[256];          /* name of gene (locus or subsequence)  */
  CharPtr name1 = NULL;     /* for protein names                    */
  CharPtr descname = NULL;  /* for protein names                    */
  CharPtr	gbp;            /* pointer to tail of gBuf              */
  CharPtr	cDNA_description = NULL;
  SeqFeatPtr sfp;           /* general feature NCBI pointer         */
  static Array a;	    /* work array for bsFlatten             */

  strcpy(gene, "");

  a = arrayReCreate(a, 32, BSunit);

  offAdjust(&subBeg,&subEnd);
  isPlus = (subBeg < subEnd);

  if (isPlus && (subEnd < 1 || subBeg > gLength) ||
      !isPlus && (subBeg < 1 || subEnd > gLength)) {
    return(1);
  }

  if (!(obj = bsCreate(key))) {
    printf("  Cannot create obj from key %s\n",name(key));
    return(1);
  }

  if (!bsFindTag(obj, _CDS) && !bsFindTag(obj, _Transposon) && !bsFindTag(obj, _Transcript)) {
    bsDestroy(obj);
    return(0);
  }

  fflush(stdout);
  fflush(stderr);
  if (gAsnDebug) printf("  Entering subsequence %s from %d to %d.\n",name(key),subBeg,subEnd);

  /* Add transposon feature: assumes text field is transposon type. */

  if (bsGetData(obj, _Transposon, _Text, &text)) {
    if (gAsnDebug) printf("    Is a transposon.\n");
    sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, NULL);
    AddIntervalToFeature(nsp,sfp,nuc_entry,NULL,subBeg,subEnd,isPlus,FALSE,FALSE);
    MakeImpFeature (nsp,sfp,"source");
    switch (gProject) {
    case YEAST:
      AddQualToImpFeature(nsp,sfp,"organism","Saccharomyes cerevisiae");
      break;
    case ELEGANS:
      AddQualToImpFeature(nsp,sfp,"organism","Caenorhabditis elegans");
      break;
    case BRIGGSAE:
      AddQualToImpFeature(nsp,sfp,"organism","Caenorhabditis briggsae");
      break;
    case HUMAN:
      AddQualToImpFeature(nsp,sfp,"organism","Homo sapiens");
      break;
    default:
      break;
    }
    AddQualToImpFeature(nsp,sfp,"transposon",text);
  }


  /* Add RNA feature: assume tRNA text is "codon name letter" */

  if (bsFindTag(obj, _Transcript) && bsFlatten(obj, 2, a)) {
    if (gAsnDebug) printf("    Is a transcript.\n");
    gpkey = arr(a,0,BSunit).k;
    if (gpkey == _mRNA) {
      if (gAsnDebug) printf("    Is a mRNA.\n");
      strcpy(rnaname,name(arr(a,1,BSunit).k));
    }
    else if (gpkey == _tRNA) {
      if (gAsnDebug) printf("    Is a tRNA.\n"); fflush(stdout);
      if (arr(a,1,BSunit).s == NULL) {
	printf("  FATAL: ===================================\n");
	printf("         Could not find the text field for %s.\n",name(key));
	printf("         Please enter codon, rna name, and rna aa \n");
	printf("         in tRNA text field and try again.\n");
	printf("         Example text: \"TAC Tyr Y\"\n");
	printf("  FATAL: ===================================\n");
	return(1);
      }
      else {
	strcpy(gBuf,arr(a,1,BSunit).s);
      }
      if (gAsnDebug) printf("   tRNA info: %s\n",gBuf); fflush(stdout);
      if (sscanf(gBuf,"%s %s %s",trna_codon,rnaname,trna_AA) != 3) {
	printf("  FATAL: Error in tRNA text \"%s\".\n",gBuf);
	return(1);
      }
      if ( strcmp(rnaname,"Sup") == 0 ) {
	printf("  Comment: adding tRNA-Sup (putative) entry.\n");
	strcpy(rnaname,"tRNA-Sup (putative)");
      }
    }
    else {
      if (gAsnDebug) printf("    Is other RNA.\n");
      if (arr(a,1,BSunit).s != NULL) {
	strcpy(rnaname,arr(a,1,BSunit).s);
      }
      else if (bsGetKey(obj, _Locus, &gpkey)) {
	if (gAsnDebug) printf("    Locus found; using %s as name of rna.\n",name(gpkey));
	strcpy(rnaname, name(gpkey));
      }
      else {
	printf("  Error in RNA: no text or locus.\n");
      }
    }

    sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, rnaname);
    AddIntervalToFeature(nsp, sfp, nuc_entry, NULL, subBeg, subEnd, isPlus, FALSE, FALSE);

    if (gpkey == _mRNA) {
      MakeRNAFeature(nsp,sfp,RNA_TYPE_mRNA,FALSE,NULL,NULL,NULL);
    }
    else if (gpkey == _tRNA) {
      MakeRNAFeature(nsp,sfp,RNA_TYPE_tRNA,FALSE,NULL,trna_AA,trna_codon);
    }
    else if (gpkey == _snRNA) {
      MakeRNAFeature(nsp,sfp,RNA_TYPE_snRNA,FALSE,NULL,NULL,NULL);
    }
    else if (gpkey == _scRNA) {
      MakeRNAFeature(nsp,sfp,RNA_TYPE_scRNA,FALSE,NULL,NULL,NULL);
    }
    else {
      MakeRNAFeature(nsp,sfp,RNA_TYPE_other,FALSE,NULL,NULL,NULL);
    }
  }


  /* Return if subsequence is not CDS. */

  if (!bsFindTag(obj, _CDS) ) {
    bsDestroy(obj);
    if (gAsnDebug)
      printf("  Exiting subsequence:  no CDS found.\n");
    fflush(stdout);
    fflush(stderr);
    return 0;
  }


  /*Add CDS;
    check for matching cDNA first and add to feature comment (start buf)
    */

  strcpy(gbp = gBuf, "");
  has_cDNA = (bsFindTag(obj, _Matching_cDNA) && bsFlatten(obj,1,a));
  if (has_cDNA) {
    char	*cp;
    unsigned long	nid;

    if (gAsnDebug)
      printf("    Found matching cDNA.\n");
    switch (gProject) {
    case ELEGANS:
      for (i = 0; i < arrayMax(a); ++i) {
	if (i > 0)
	  gbp = Nlm_StrMove(gbp, "; ");
	gbp += sprintf(gbp, "coded for by C. elegans cDNA %s",
		       name(arr(a,i,BSunit).k));
      }
      break;
    case HUMAN:
      gbp += sprintf(gbp, "coded for by human cDNA%s %s",
		     (arrayMax(a) > 1 ? "s" : ""),
		     cp = name(arr(a,0,BSunit).k));
      if ((nid = get_cDNA_NID(cp)) != 0) {
	gbp += sprintf(gbp, " (NID:g%lu)", nid);
      }
      if (arrayMax(a) > 1) {
	for (i = 1; i + 1 < arrayMax(a) ; ++i) {
	  gbp += sprintf(gbp, ", %s", cp = name(arr(a,i,BSunit).k));
	  if ((nid = get_cDNA_NID(cp)) != 0) {
	    gbp += sprintf(gbp, " (NID:g%lu)", nid);
	  }
	}
	gbp += sprintf(gbp, " and %s", cp = name(arr(a,i,BSunit).k));
	if ((nid = get_cDNA_NID(cp)) != 0) {
	  gbp += sprintf(gbp, " (NID:g%lu)", nid);
	}
      }
      cDNA_description = StrSave(gBuf);
      strcpy(gbp = gBuf, "");
      break;
    default:
      break;
    }
  }


  /* Add TSL_site to feature comment. */

  if (bsGetData(obj, _TSL_site, _Int, &i)) {
    if (gAsnDebug)
      printf("    Adding TSL site in comment field.\n");
    if (gBuf[0] != NULLB)
      gbp = StrMove(gbp, "; possible");
    else
      gbp = StrMove(gbp, "Possible");
    gbp += sprintf(gbp,
		   " TSL site at %d",isPlus?i+subBeg-1:subBeg+1-i);
  }


  /* Check for locus and set gene name (use subsequence if no locus). */

  has_Locus = bsGetKey(obj, _Locus, &gpkey);
  lcp = MemNew(strlen(name(gpkey)) + 2);
  sprintf(lcp, "%s", name(gpkey));
  if (has_Locus && (*lcp != '*')) {
    if (gAsnDebug)
      printf("    Locus found; using \"%s\" as name of gene.\n", name(gpkey));
    strcpy(gene, name(gpkey));       /*ss --where locus is set to gene*/
  }
  else {
    if (gAsnDebug)
      printf("    No locus found; using \"%s\" as name of gene.\n", name(key));
    strcpy(gene, "");
    switch (gProject) {
    case HUMAN:
      if (strncasecmp(name(key), gClonePat, strlen(gClonePat)) == 0
	  || strncasecmp(name(key), gRawClonePat, strlen(gRawClonePat)) == 0) {
	strcpy(gene, Institution);
	strcat(gene, ":");
      }
      break;
    default:
      break;
    }
    strcat(gene, name(key));
  }
  lcp = MemFree(lcp);

  /* Build CDS feature. */

  /* Get source exons and gene limits (cdsBeg,cdsEnd) */

  if (gAsnDebug)
    printf("    Start analysis of source exons...\n");

  cdsBeg = subEnd;
  cdsEnd = subBeg;
  gMiss1 = gMiss2 = 0;
  if (bsFindTag(obj, _Source_Exons) && bsFlatten(obj, 2, a)) {
    sortExons(a);
    for (i = 0; i < arrayMax(a) ; i += 2) {
      if (gAsnDebug)
	printf("    >>> %5d %5d\n",
	       arr(a,i,BSunit).i,arr(a,i+1,BSunit).i);
      if (isPlus) {
	arr(a,i  ,BSunit).i += subBeg - 1;
	arr(a,i+1,BSunit).i += subBeg - 1;
	clipAdd(arr(a,i,BSunit).i, arr(a,i+1,BSunit).i);
      }
      else {
	arr(a,i  ,BSunit).i = subBeg + 1 - arr(a,i  ,BSunit).i;
	arr(a,i+1,BSunit).i = subBeg + 1 - arr(a,i+1,BSunit).i;
	clipAdd(arr(a,i+1,BSunit).i, arr(a,i,BSunit).i);
      }
    }
    cdsBeg = arr(a,0            ,BSunit).i;
    cdsEnd = arr(a,arrayMax(a)-1,BSunit).i;
  }
  else {
    if (isPlus)
      clipAdd(subBeg, subEnd);
    else
      clipAdd(subEnd, subBeg);
  }

  if (isPlus) {
    if (!gMiss1 && bsFindTag (obj, _Start_not_found))
      if (bsGetData (obj, _bsRight, _Int, &i))
	gMiss1 = 2 + i;
      else
	gMiss1 = 3;		/* no frame change */
    if (!gMiss2 && bsFindTag (obj, _End_not_found))
      gMiss2 = 3;
  }
  else {
    if (!gMiss2 && bsFindTag (obj, _Start_not_found))
      if (bsGetData (obj, _bsRight, _Int, &i))
	gMiss2 = 2 + i;
      else
	gMiss2 = 3;		/* no frame change */
    if (!gMiss1 && bsFindTag (obj, _End_not_found))
      gMiss1 = 3;
  }

  start_not_found = FALSE;
  if (bsFindTag(obj, _Start_not_found)) {
    start_not_found = TRUE;
    if (gAsnDebug)
      printf("  WARNING: Start was NOT found\n");
  }

  end_not_found = FALSE;
  if (bsFindTag(obj, _End_not_found)) {
    end_not_found = TRUE;
    if (gAsnDebug)
      printf("  WARNING: End was NOT found\n");
  }

  if (has_Locus || has_cDNA)
    sfp = FeatureBuild(nsp, the_entry, (start_not_found || end_not_found), EVIDENCE_NOT_SET, FALSE, gBuf);
  else
    sfp = FeatureBuild(nsp, the_entry, (start_not_found || end_not_found), EVIDENCE_NOT_EXPERIMENTAL, FALSE, gBuf);

  for (i = 0; i + 1 < arrayMax(a) ; i += 2) {
    snf = enf = FALSE;
    if (i == 0)
      snf = start_not_found;
    if (i + 2 == arrayMax(a))
      enf = end_not_found;
    if (gAsnDebug && snf)
      printf("    AddIntervalToFeature, start_not_found.\n");
    if (gAsnDebug && enf)
      printf("    AddIntervalToFeature, end_not_found.\n");
    AddIntervalToFeature(nsp, sfp, nuc_entry, NULL,
			 arr(a,   i, BSunit).i,
			 arr(a, i+1, BSunit).i, isPlus, snf, enf);
  }

  if (gAsnDebug)
    printf("    Source exons added.\n");

  MakeCdRegionFeature(nsp, sfp, 1, 0, NULL, NULL);
  prot_entry = TranslateCdRegion(nsp, sfp, nucprot_entry, name(key), NULL, NULL, 0);
  AddCompleteness(nsp, prot_entry, sfp);

  if (gAsnDebug)
    printf("    Coding translation successful.\n");

  if (prot_entry != NULL) {
    if (gAsnDebug)
      printf("    Adding protein entry.\n");

    /* Only give protein name1 if we have synonyms for the protein and  */
    /* and we actually know the product: only know what the protein     */
    /* encodes if know the locus (at least in elegans).                 */
    /*
      For HUMAN, if there is no Locus or it starts with a '*',
      then it is assumed the protein has been positively identified,
      but we just don't have a gene name.  In this case, the '*' is
      skipped over and not dumped, and what will be seen in the GenBank
      flat file is a /product="{briefID}" and /gene="{pseudoname}". */

    has_DBRemark = bsGetKey(obj, _DB_remark, &gpkey);
    if (has_DBRemark) {
      switch (gProject) {
      case HUMAN:
	if (bsGetKey(obj, _Brief_identification, &gpkey)) {
	  while (bsGetKey(obj, _bsDown, &gpkey))      /* goto last Brief_identification */
	    ;
	  name1 = StringSaveNoNull(name(gpkey));      /*ss --such as /product="semaphorin IV" from BI*/
	}
      default:
	break;
      }
      bsGetKey(obj, _DB_remark, &gpkey);
      while (bsGetKey(obj, _bsDown, &gpkey))	      /* goto last db remark */
	;
      if (has_Locus) {
	if (gAsnDebug)
	  printf("    Locus taken from DB remark: %s.\n", name(gpkey));
	if (name1 == NULL)
	  name1 = StringSaveNoNull(name(gpkey));
	switch (gProject) {
	  /*char	*cp;*/
	case HUMAN:
	  if (strncasecmp(name(key), gRawClonePat, strlen(gRawClonePat)) == 0) {
	    /*cp = MemNew(strlen(name(key)) + strlen(Institution) + 2);
	    sprintf(cp, "%s:%s", Institution, name(key));
	    descname = catnames(cp, name(gpkey));  
	    cp = MemFree(cp);*/
	    descname = StringSaveNoNull(name(gpkey));         /*ss --such as /note="" from DB remark*/
	  }
	  else
	    descname = catnames(name(key), name(gpkey));
	  break;
	default:
	  descname = StringSaveNoNull(name(key));
	  break;
	}
      }
      else {
	if (gAsnDebug)
	  printf("    No locus found from DB remark.\n");
	if (gProject != HUMAN)
	  name1 = MemFree(name1);
	if (gAsnDebug)
	  printf("   name1 = \"%s\"\n", (name1 != NULL ? name1 : ""));
	if (name1 == NULL || *name1 != '*') {
	  descname = catnames(name1, name(gpkey));
	  name1 = MemFree(name1);
	}
	else {
	  descname = StringSaveNoNull(name(gpkey));
	}
      }
    } /* has_DBRemark */
    else {
      if (gAsnDebug)
	printf("    No DBremark...checking for briefID\n");
      if (bsGetKey(obj, _Brief_identification, &gpkey)) {
	while (bsGetKey(obj, _bsDown, &gpkey))    /* goto last Brief_identification */ 
	  ;
	if (gAsnDebug && name1 != NULL)
	  printf("    Chucking name1=\"%s\"\n", name1);
	name1 = MemFree(name1);
	if (has_Locus) {
	  if (gAsnDebug)
	    printf("    Locus found from brief id: %s.\n",name(gpkey));
	  name1 = StringSaveNoNull(name(gpkey));
	  descname = StringSaveNoNull(name(key));
	}
	else {
	  char	*cp;
	  if (gAsnDebug)
	    printf("    No locus found from brief id.\n");
	  cp = MemNew(strlen(name(gpkey)) + 100);
	  switch (gProject) {
	  case HUMAN:
	    sprintf(cp, "%s", name(gpkey));
	    break;
	  default:
	    sprintf(cp, "Similar to %s", name(gpkey));
	    break;
	  }
	  /* Prepend cp to the growing gBuf */
	  if (gAsnDebug)
	    printf("    concatenating \"%s\" & \"%s\"\n", cp, gBuf);
	  /* if (*cp == '*') {
	    name1 = cp;
	    descname = gBuf;
	  }
	  else {*/
	  descname = catnames(cp, gBuf);
	  cp = MemFree(cp);
	}
      }
    }
 

    descname = catnames(descname, cDNA_description);   /*ss --where cDNA gets added*/
    
    sfp = FeatureBuild(nsp, prot_entry, (start_not_found || end_not_found), FALSE, FALSE, NULL);
    AddIntervalToFeature(nsp, sfp, prot_entry, NULL, 1, -1, TRUE, start_not_found, end_not_found);
    
    /* Argument name1 is the thing that goes in for product and  */
    /* should only be used when we KNOW what the product is. */
    
    if (gAsnDebug) {
      if (name1 != NULL)
	printf("    Protein name1 is %s\n", name1);
      else
	printf("    Protein name1 is NULL.\n");
      if (descname != NULL)
	printf("    Protein descriptive name is %s\n", descname);
      else
	printf("    Protein descriptive name is NULL.\n");
    }
  
    MakeProteinFeature(nsp, sfp,
		       (name1 == NULL || *name1 != '*' ? name1 : name1+1),
		       NULL, NULL, descname, NULL, NULL,
		       NULL, NULL, NULL, NULL);
    
    name1 = MemFree(name1);
    descname = MemFree(descname);
    
    if (gAsnDebug)
      printf("    Adding gene %s from %d to %d.\n", gene, cdsBeg, cdsEnd);
    
    sfp = FeatureBuild(nsp, nuc_entry, (start_not_found || end_not_found), FALSE, FALSE, NULL);
    
    /* MakeGeneFeature args: nsp, sfp, gene_name, allele,   */
    /* descriptive name, map_location, is_pseudogene,       */
    /* genetic_database, gene_id, synonym1, syn2, syn3)     */
    
    MakeGeneFeature(nsp,sfp,gene,NULL,NULL,NULL,FALSE,NULL,NULL,NULL,NULL,NULL);
    
    AddIntervalToFeature(nsp,sfp,nuc_entry,NULL,cdsBeg,cdsEnd,isPlus,start_not_found,end_not_found);
  }

  if (gAsnDebug)
    printf("  All done with subsequence.\n");
  bsDestroy(obj);
  arrayDestroy(a);
  return 0;
}

/*********************************************************************
 * emblFeature:
 *
 * Adds EMBL_feature to NCBI submission.
 ********************************************************************/

static int emblFeature(KEY tag, int beg, int end, BSunit u)
{
  SeqFeatPtr sfp;
  BOOL isPlus = (beg <= end);

  if (u.s)
    strcpy(gBuf,u.s);
  else
    strcpy(gBuf," ");

  if (gAsnDebug)
    printf("  Adding feature %s from %d to %d.\n",name(tag),beg,end);

  sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, gBuf);
  AddIntervalToFeature(nsp, sfp, nuc_entry, NULL, beg, end, isPlus, FALSE, FALSE);
  MakeImpFeature(nsp, sfp, name(tag));

  return(0);
}

/*********************************************************************
 * parseQualifierString:
 *
 * Parses text string from ACeDB Feature type GB_*. Text string
 * contains multiple GenBank qualifers. See dumpGBFeatures for
 * format of the qualifier string.
 *
 * D. Ficenec WUGSC.
 *
 ********************************************************************/

static int parseQualString(char *qstr, char *qbuf, char *vbuf)
{
  static char *qhere = NULL;
  static char *nextq = NULL;

  char *currq;
  int   nscan;

  if (qstr != NULL) {
    if (strlen(qstr) == 0) {
      printf("  WARNING: Empty qualifer string.\n");
      qhere = NULL;
      return 0;
    }
    if (strncmp(qstr,"/",1) != 0) {
      printf("  WARNING: Qualifer string must start with /\n");
      qhere = NULL;
      return 0;
    }
    qhere = (char *) malloc(strlen(qstr)+1);
    strcpy(qhere,qstr);
    nextq = qhere + 1;
  }

  if (nextq == NULL) {
    qhere = NULL;
    return 0;
  }

  currq = nextq;
  if ((nextq = (char *) strstr(currq+2," /")) != NULL) {
    nextq[0] = '\0';
    nextq += 2;
  }

  nscan = sscanf(currq,"%[^= ]%*[= ]%[^\f]",qbuf,vbuf);
  if (nscan == 1) {
    strcpy(vbuf,"");
  }

  return 1;
}

/*********************************************************************
 * dumpGBFeatures:
 *
 * Adds other features (attached to main sequence only) to NCBI
 * submission. Only features which start with GB_ and have the
 * EMBL_dump_YES flag set are dumped. The text contains the feature
 * note and/or qualifier(s); each one starting with a blank followed
 * by a / (/ in qualifier must be preceded by a non-space character).
 * If the qualifier has a value, then use "/qualifer=value" with
 * optional spaces between the qualifier and the the equals sign.
 *
 * D. Ficenec WUGSC.
 *
 ********************************************************************/

static int dumpGBFeatures(KEY tag)
{
  OBJ mobj,sobj;
  KEY mkey;
  BOOL isPlus;
  BSMARK mmark = NULL;
  BSMARK smark = NULL;
  SeqFeatPtr sfp;
  int dumpYes;
  int hasEMBLfeature;
  int mdone,sdone;
  int mbeg,mend;
  float msig;
  /*  char qbuf[128];
  char vbuf[2048];*/
  char ftype[128];
  char ntxt[256];
  char	slashname[256];
  char *mtxt, *cp;
  char *etype;

  sobj = bsCreate(tag);

  if (gAsnDebug) printf("  Adding GB_ features...\n");

  if (bsFindTag(sobj, _Feature)) {
    bsGetKeyTags(sobj,_bsRight,&mkey);
    mdone = 0;
    while (!mdone) {
      mobj = bsCreate(mkey);
      hasEMBLfeature = bsGetData(mobj, _EMBL_feature, _Text, &etype);
      if (hasEMBLfeature) strcpy(ftype,etype);
      bsDestroy(mobj);
      if (strncmp(name(mkey),"GB_",3) == 0 && hasEMBLfeature) {
	mmark = bsMark(sobj,mmark);
	bsGetData(sobj,_bsRight,_Int,&mbeg);
	sdone = 0;
	while (!sdone) {
	  smark = bsMark(sobj,smark);
	  bsGetData(sobj,_bsRight,_Int,&mend);
	  bsGetData(sobj,_bsRight,_Float,&msig);
	  bsGetData(sobj,_bsRight,_Text,&mtxt);
	  strcpy(ntxt, mtxt);
	  slashname[0] = NULLB;
	  if (mtxt[0] == '/' && isalpha(mtxt[1]) && (cp = strchr(mtxt,'=')) != NULL) {
		memcpy(slashname, mtxt+1, cp-mtxt-1);
		slashname[cp-mtxt-1] = NULLB;
		strcpy(ntxt, cp+1);
	  }
	  dumpYes = (bsPushObj(sobj) && bsFindTag(sobj,_EMBL_dump_YES)) ? 1 : 0;
	  if (dumpYes) {
	    isPlus = (mbeg <= mend);
	    if (strncasecmp(ntxt,"Grail",5) ==0) {
	      sfp = FeatureBuild(nsp, nuc_entry, FALSE, EVIDENCE_NOT_EXPERIMENTAL, FALSE, NULL);
	    }
	    else {
	      sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, NULL);
	    }
	    MakeImpFeature(nsp, sfp, ftype);
	    AddIntervalToFeature(nsp, sfp, nuc_entry, NULL, mbeg, mend, isPlus, FALSE, FALSE);
      	    /*if (parseQualString(mtxt,qbuf,vbuf)) {
	      AddQualToImpFeature(nsp, sfp, qbuf, vbuf);
	      while (parseQualString(NULL,qbuf,vbuf)) {
		AddQualToImpFeature(nsp, sfp, qbuf, vbuf);
	      }
	    }*/
	    if (strcmp(ftype, "exon") == 0) {
		  if (slashname[0] == NULLB)
			strcpy(slashname, "note");
	      AddQualToImpFeature(nsp,sfp,slashname,ntxt);
		    if (gAsnDebug) {
		      printf("  >%s as %s [%6d,%6d] /%s=\"%s\"\n",name(mkey),ftype,mbeg,mend,slashname,ntxt);
		    }
	    }
	    else {
			if (slashname[0] == NULLB)
				strcpy(slashname, "rpt_family");
	      AddQualToImpFeature(nsp,sfp,slashname,ntxt);
		    if (gAsnDebug) {
		      printf("  >%s as %s [%6d,%6d] /%s=\"%s\"\n",name(mkey),ftype,mbeg,mend,slashname,ntxt);
		    }
	    }
	  }
	  bsGoto(sobj,smark);
	  sdone = !bsGetData(sobj,_bsDown,_Int,&mbeg);
	}
	bsGoto(sobj,mmark);
      }
      mdone = !bsGetKeyTags(sobj,_bsDown,&mkey);
    }
  }

  bsDestroy(sobj);
  bsMarkFree(mmark);
  bsMarkFree(smark);

}

/*********************************************************************
 * dumpMiscFeatures:
 *
 * Adds other features (attached to main sequence only) to NCBI
 * submission.
 *
 * S. Shin WUGSC.
 *
 ********************************************************************/

static int dumpMiscFeatures(KEY tag)
{
  OBJ mobj,sobj;
  KEY mkey;
  BOOL isPlus;
  BSMARK mmark = NULL;
  BSMARK smark = NULL;
  SeqFeatPtr sfp;
  int dumpYes;
  int hasEMBLfeature;
  int mdone,sdone;
  int mbeg,mend;
  Boolean	cpg_island;
  float msig;
  /*  char qbuf[128];
  char vbuf[2048];*/
  char qtxt[2176];
  char ftype[128];
  char *mtxt;
  char *etype;

  sobj = bsCreate(tag);

  if (gAsnDebug) printf("  Adding Misc. features...\n");

  if (bsFindTag(sobj, _Feature)) {
    bsGetKeyTags(sobj,_bsRight,&mkey);
    mdone = 0;
    while (!mdone) {
      mobj = bsCreate(mkey);
      hasEMBLfeature = bsGetData(mobj, _EMBL_feature, _Text, &etype);
      if (hasEMBLfeature)
		strcpy(ftype,etype);
      bsDestroy(mobj);
      if (strncasecmp(name(mkey),"Predicted_CpG_island",20) == 0 && hasEMBLfeature) {
	mmark = bsMark(sobj,mmark);
	bsGetData(sobj,_bsRight,_Int,&mbeg);
	sdone = 0;
	while (!sdone) {
	  smark = bsMark(sobj,smark);
	  bsGetData(sobj,_bsRight,_Int,&mend);
	  bsGetData(sobj,_bsRight,_Float,&msig);
	  bsGetData(sobj,_bsRight,_Text,&mtxt);
	  strcpy(qtxt, "CpG_island ");
	  /* enclose mtxt between matching parentheses */
	  if (mtxt[0] != '(') /* open parenthesis (if not already present) */
		strcat(qtxt, "(");
	  strcat(qtxt,mtxt);
	  if (mtxt[strlen(mtxt)-1] != ')') /* close parenthesis (if not present) */
		strcat(qtxt, ")");

	  /*dumpYes = (bsPushObj(sobj) && bsFindTag(sobj,_EMBL_dump_YES)) ? 1 : 0;*/
	  dumpYes = 1;
	  if (dumpYes) {
	    isPlus = (mbeg <= mend);
	    sfp = FeatureBuild(nsp, nuc_entry, FALSE, EVIDENCE_NOT_SET, FALSE, NULL);
	    MakeImpFeature(nsp, sfp, ftype);
	    AddIntervalToFeature(nsp, sfp, nuc_entry, NULL, mbeg, mend, isPlus, FALSE, FALSE);
	    /*	    if (parseQualString(mtxt,qbuf,vbuf)) {
	      AddQualToImpFeature(nsp, sfp, qbuf, vbuf);
	      while (parseQualString(NULL,qbuf,vbuf)) {
		AddQualToImpFeature(nsp, sfp, qbuf, vbuf);
	      }
	    }*/
	    AddQualToImpFeature(nsp,sfp,"note",qtxt);
     	    if (gAsnDebug) {
	      printf("  >%s as %s [%6d,%6d] %s\n",name(mkey),ftype,mbeg,mend,mtxt);
	    }
	  }
	  bsGoto(sobj,smark);
	  sdone = !bsGetData(sobj,_bsDown,_Int,&mbeg);
	}
	bsGoto(sobj,mmark);
      }
      mdone = !bsGetKeyTags(sobj,_bsDown,&mkey);
    }
  }

  bsDestroy(sobj);
  bsMarkFree(mmark);
  bsMarkFree(smark);

}

/*********************************************************************
 * asnDoDump:
 *
 * Called from asnDumpKey.
 *
 * Process:
 *         [1] Check sequence (DNA, Genomic_Canonical);
 *         [2] Get sequence info (locus, finisher, etc);
 *         [3] Add sequence info and publications;
 *         [4] Add EMBL features;
 *         [5] Add method features;
 *         [6] Add subsequence features;
 ********************************************************************/

static int asnDoDump (KEY seq)
{
	KEY	key;                     /* General purpose key               */
	KEY	clone;                   /* Key for clone                     */
	OBJ	Seq;                     /* Sequence object (key seq)         */
	OBJ	obj;                     /* General purpose object            */
	SeqFeatPtr	sfp;              /* General purpose seq feature ptr   */
	PubPtr	pub;                  /* General purpose publication ptr   */
	static Array	a;	       /* Working array for bsFlatten       */
	Array	dna;                   /* Holds nucleotide sequence         */
	int	i,ii;                       /* General purpose int               */
	int	isCDS=0;                 /* Flag for CDS in subsequences      */
	int	status;                  /* Subroutine return status          */
	char	*message = 0;           /* Error message                     */
	char	*dnastr;                /* DNA string.                       */
	char	namebuf1[64];           /* Last name or first initial        */
	char	namebuf2[64];           /* Last name or first initial        */
	char	*asnbr = NULL;          /* Pointer to accession number       */
	char	*cp;
	CharPtr	gbp;

  /* Initialize globals and set working array. */

  if (! _scRNA) {
    lexaddkey("scRNA", &_scRNA, 0);
    lexaddkey("misc_RNA", &_misc_RNA, 0);
    lexaddkey("pseudogene", &_pseudogene, 0);
    lexaddkey("Database", &_Database, 0);
    lexaddkey("Clone_Library", &_Clone_Library, 0);
    lexaddkey("Feature", &_Feature, 0);
    lexaddkey("EMBL_feature", &_EMBL_feature, 0);
    lexaddkey("EMBL_threshold", &_EMBL_threshold, 0);
    lexaddkey("EMBL_qualifier", &_EMBL_qualifier, 0);
    lexaddkey("EMBL_dump_YES", &_EMBL_dump_YES, 0);
    lexaddkey("EMBL_dump_NO", &_EMBL_dump_NO, 0);
  }

  nucprot_entry = NULL;
  nuc_entry = NULL;

  a = arrayReCreate(a, 32, BSunit);

  /* Check sequence. */

  if (!(class(seq) == _VSequence) || !(Seq = bsCreate(seq))) {
    message = "  Sequence object missing.";
    goto end_asnDoDump;
  }

  if (!bsFindTag(Seq, _Genomic_Canonical)) {
    gProject = NOPROJECT;
    message = "  Sequence must be Genomic_Canonical to do ASN.1 dump.";
    goto end_asnDoDump;
  }

  /* Get DNA array. */

  if (!bsGetKey(Seq, _DNA, &key) || !(dna = dnaGet(key))) {
    message = "  No DNA attached to Sequence object.";
    goto end_asnDoDump;
  }

  /* Get info: name, clone, length, accession number, locus, chromosome. */

  strcpy(gSequence,FixCloneName(name(seq)));

  bsGetKey(Seq, _Clone, &clone);

  strcpy(gRawClone, name(clone));
  strcpy(gClone, FixCloneName(gRawClone));

  strcpy(gClonePat, gClone);
  strcat(gClonePat, ".");
  strcpy(gRawClonePat, gRawClone);
  strcat(gRawClonePat, ".");

  obj = bsCreate(clone);

  if (bsFindTag(obj, _Type) && bsFlatten(obj, 1, a)) {
    strcpy(gCloneType,name(arr(a,0,BSunit).k));
  }
  else {
    printf("  WARNING: Could not find Clone name.\n");
    strcpy(gCloneType,"clone");
  }


  if (bsFindTag(obj, _Clone_Library) && bsFlatten(obj, 1, a)) {
    strcpy(gLibrary,name(arr(a,0,BSunit).k));
  }
  else {
    strcpy(gLibrary,"library");
    printf("  WARNING: Could not find Clone_Library name.\n");
  }

  if (bsFindTag(obj, _Map) && bsFlatten(obj, 1, a)) {
    strcpy(gMapPosition,name(arr(a,0,BSunit).k));
  }
  else {
    printf("  WARNING: Could not find library map position under clone.\n");
    strcpy(gMapPosition,"map_position");
  }

  if (gAsnDebug) {
    printf("  Sequence is %s\n",gSequence);
    printf("  Clone is %s of type %s\n",gClone,gCloneType);
    printf("  Library is %s\n",gLibrary);
    printf("  Map position is %s\n",gMapPosition);
  }

  gLength = arrayMax(dna);

  if (gAsnDebug) printf("  DNA length is %d.\n", gLength);

  if (bsFindTag(Seq, _Database) && bsFlatten(Seq,3,a) && getenv("ASNACCESSION") == NULL) {
    asnbr = arr(a,2,BSunit).s;
  }

    if (bsFindTag(Seq, _Keyword) && bsFlatten (Seq, 1, a)) {
      ii=0;
      for (i = 0; i < arrayMax(a) ; ++i) {
	sprintf(keyword[ii++],"%s",name(arr(a,i,BSunit).k));
      }
    }


   if (gAsnDebug) printf("keyword[0] is %s\n",keyword[0]);

  switch (gProject) {
  case YEAST:
    sprintf(gLocus, "YSCL%s", gSequence);
    strcpy(gChromosome,"XII");
    if (gAsnDebug) printf("  Yeast chromosome XII.\n");
    break;

  case BRIGGSAE:
    sprintf(gLocus, "CBR%s", gSequence);
    if (bsGetKey(Seq, _Map, &key)) {
      sscanf(name(key),"Sequence-%s",gChromosome);
      if (gAsnDebug) printf("  Briggsae chromosome %s.\n",gChromosome);
    }
    else {
      printf("  WARNING: Could not find chromosome for C. briggsae.\n");
      strcpy(gChromosome,"unknown");
    }
    break;

  case ELEGANS:
    sprintf(gLocus, "CEL%s", gSequence);
    if (bsGetKey(Seq, _Map, &key)) {
      sscanf(name(key),"Sequence-%s",gChromosome);
      if (gAsnDebug) printf("  Elegans chromosome %s.\n",gChromosome);
    }
    else {
      message = " Could not find chromosome (no map).";
      goto end_asnDoDump;
    }
    break;

  case HUMAN:
	if (asnbr != NULL && asnbr[0] != NULLB)
    	sprintf(gLocus, "HS%s", asnbr);
	/*sprintf(gLocus, "HUM%s", gSequence);*/
    if (bsGetKey(Seq, _Map, &key)) {
      sscanf(name(key),"Sequence-%s",gChromosome);
      if (gAsnDebug) printf("  Human chromosome %s.\n",gChromosome);
    }
    else {
      message = " Could not find chromosome (no map).";
      goto end_asnDoDump;
    }
    break;

  default:
    message = "  Project not implemented yet.";
    goto end_asnDoDump;
  }


	/*
	 * Determine is any CDS subsequences. Create appropriate
	 * sequence entry type. Add DNA.
	 */

	if (bsFindTag(Seq, _Subsequence) && bsFlatten(Seq, 3, a)) {
		for (i=0; i<arrayMax(a); i+=3) {
			if (obj = bsCreate(arr(a,i,BSunit).k)) {
				if (bsFindTag(obj, _CDS)) isCDS++;
				bsDestroy(obj);
			}
		}
	}
	if (gAsnDebug) printf("  Found %d CDS subsequences.\n",isCDS);

	if (isCDS) {
		nucprot_entry = AddNucProtToSubmission(nsp);
		nuc_entry = AddSeqToNucProtEntry
				(nsp, nucprot_entry, gRawClone,
				(gLocus[0] != '\0' ? gLocus : NULL), asnbr, 0,
				MOLECULE_CLASS_DNA, MOLECULE_TYPE_GENOMIC, (Int4)gLength,
				TOPOLOGY_LINEAR, STRANDEDNESS_DOUBLE);
		the_entry = nucprot_entry;
	}
	else {
		printf("  WARNING: No CDS subsequences found.\n");
		nuc_entry = AddSeqOnlyToSubmission(nsp, gRawClone,
				(gLocus[0] != '\0' ? gLocus : NULL), asnbr, 0,
				MOLECULE_CLASS_DNA, MOLECULE_TYPE_GENOMIC, (Int4)gLength,
				TOPOLOGY_LINEAR, STRANDEDNESS_DOUBLE);
		the_entry = nuc_entry;
	}

	switch (gProject) {
	case HUMAN:
		if (!AddTechToEntry(nsp, nuc_entry, MI_TECH_htgs_3))
		/*if (!AddTechToEntry(nsp, the_entry, MI_TECH_htgs_3))*/
			printf("  AddTechToEntry() FAILED\n");
		break;
	default:
		break;
	}

	dnastr = (char *) malloc(gLength+1);
	dnastr[gLength] = '\0';
	for (i = 0; i < gLength; ++i) {
		dnastr[i] = TO_UPPER(dnaDecodeChar[arr(dna,i,char)]);
	}
	AddBasesToBioseq(nsp, nuc_entry, dnastr);
	free(dnastr);

  /*
   * Master switch statement for adding entry annotation:
   * organism type, taxonomy, refrences, comments, etc.
   */

	switch (gProject) {
	case YEAST:

    AddOrganismToEntry(nsp, the_entry, "Saccharomyces cerevisiae",
		       "S. cerevisiae", NULL, "S288C (AB972)", NULL,
		       NULL, NULL );

    strcpy(gBuf,"Eukaryota; Plantae; Thallobionta; Eumycota; ");
    strcat(gBuf,"Hemiascomycetes; Endomycetales; Saccharomycetaceae");
    /* NULL NULL NULL is where I can enter keywords */
    AddGenBankBlockToEntry(nsp, the_entry, gBuf, "PLN",NULL,NULL,NULL);

    sprintf(gBuf, "Saccharomyces cerevisiae chromosome XII cosmid %s", gSequence);
    AddTitleToEntry(nsp, nuc_entry, gBuf);
    AddReferences(Seq, nsp, the_entry);
    sprintf(gBuf,"The sequence of S. cerevisiae cosmid %s", gSequence);
    AddFinishers(Seq, nsp, the_entry, gBuf);

    strcpy(gBuf,"            Submitted by: ~");
    strcat(gBuf,"            Genome Sequencing Center ~");
    strcat(gBuf,"            Department of Genetics, Washington University, ~");
    strcat(gBuf,"            St. Louis, MO 63110, USA~");
    strcat(gBuf,"            e-mail: mj@sequencer.wustl.edu~");

    if (bsFindTag (Seq, _DB_remark) && bsFlatten (Seq, 1, a)) {
      strcat(gBuf, "~\n            NEIGHBORING COSMID INFORMATION:\n");
      for (i = 0 ; i < arrayMax(a) ; ++i) {
	sprintf(gBuf+strlen(gBuf), "~~%s",name(arr(a,i,BSunit).k));
      }
    }
    if (bsFindTag (Seq, _TSL_site)) {
      sprintf(gBuf+strlen(gBuf), "~\"TSL\" = trans-spliced leader.");
    }
    if (strlen(gBuf) > (size_t) 0) {
      AddCommentToEntry(nsp, nuc_entry, gBuf);
    }

    break;

  case BRIGGSAE:

    AddOrganismToEntry(nsp, the_entry, "Caenorhabditis briggsae",
		       "C. briggsae hermaphrodite mixed whole animal DNA.",
		       NULL, "GujArat G16", NULL, NULL, NULL );

    sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, NULL);
    AddIntervalToFeature(nsp,sfp,nuc_entry,NULL,1,gLength,TRUE,FALSE,FALSE);
    MakeImpFeature (nsp,sfp,"source");
    AddQualToImpFeature(nsp,sfp,"organism","Caenorhabditis briggsae");
    if (strcmp(gChromosome,"unknown") != 0) {
      AddQualToImpFeature(nsp,sfp,"chromosome",gChromosome);
    }
    AddQualToImpFeature(nsp,sfp,"clone",gSequence);
    AddQualToImpFeature(nsp,sfp,"strain","GujArat G16");

    strcpy(gBuf,"Eukaryota; Animalia; Eumetazoa; Nematoda; Secernentea; ");
    strcat(gBuf,"Rhabditida; Rhabditina; Rhabditoidea; Rhabditidae.");
    AddGenBankBlockToEntry(nsp, the_entry, gBuf, "INV", NULL, NULL, NULL);
    sprintf(gBuf, "Caenorhabditis briggsae cosmid %s", gSequence) ;
    AddTitleToEntry(nsp, nuc_entry, gBuf);
    AddBriggsaeReference(nsp, the_entry);  /* Remove when briggsae papers in ACeDB. */
    AddReferences(Seq, nsp, the_entry);
    sprintf(gBuf, "The sequence of C. briggsae cosmid %s", gSequence) ;
    AddFinishers(Seq, nsp, the_entry, gBuf);

    strcpy(gBuf,"Submitted by: ~");
    strcat(gBuf,"         Genome Sequencing Center ~");
    strcat(gBuf,"         Department of Genetics, Washington University, ~");
    strcat(gBuf,"         St. Louis, MO 63110, USA ~");
    strcat(gBuf,"         e-mail: mmarra@watson.wustl.edu ~");

    if (bsFindTag(Seq, _DB_remark) && bsFlatten (Seq, 1, a)) {
      strcat(gBuf, "~~            NEIGHBORING COSMID INFORMATION:~");
      for (i = 0; i < arrayMax(a) ; ++i) {
	sprintf(gBuf+strlen(gBuf), "~~%s",name(arr(a,i,BSunit).k));
      }
    }

    if (isCDS) {
      if (gAsnDebug) printf("  Adding comment about GeneFinder predictions.\n");
      sprintf(gBuf+strlen(gBuf), " ~ ~NOTES:~~Coding sequences below are predicted from");
      sprintf(gBuf+strlen(gBuf), " computer analysis, using the program Genefinder");
      sprintf(gBuf+strlen(gBuf), "(P. Green and L. Hillier, ms in preparation).");
    }

    if (bsFindTag (Seq, _TSL_site)) {
      sprintf(gBuf+strlen(gBuf), "~\"TSL\" = trans-spliced leader.");
    }

    if (strlen(gBuf) > (size_t) 0) AddCommentToEntry(nsp, nuc_entry, gBuf);

    break;

  case ELEGANS:

    AddOrganismToEntry(nsp, the_entry, "Caenorhabditis elegans",
		       "C. elegans hermaphrodite mixed whole animal DNA.",
		       NULL, "Bristol N2", NULL, NULL, NULL );

    sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, NULL);
    AddIntervalToFeature(nsp,sfp,nuc_entry,NULL,1,gLength,TRUE,FALSE,FALSE);
    MakeImpFeature (nsp,sfp,"source");
    AddQualToImpFeature(nsp,sfp,"organism","Caenorhabditis elegans");
    AddQualToImpFeature(nsp,sfp,"chromosome",gChromosome);
    AddQualToImpFeature(nsp,sfp,"clone",gSequence);
    AddQualToImpFeature(nsp,sfp,"strain","Bristol N2");

    strcpy(gBuf,"Eukaryota; Animalia; Eumetazoa; Nematoda; Secernentea; ");
    strcat(gBuf,"Rhabditida; Rhabditina; Rhabditoidea; Rhabditidae.");
    AddGenBankBlockToEntry(nsp, the_entry, gBuf, "INV", NULL, NULL, NULL);
    sprintf(gBuf, "Caenorhabditis elegans cosmid %s", gSequence) ;
    AddTitleToEntry(nsp, nuc_entry, gBuf);
    AddElegansReference(nsp, the_entry);  /* Remove when elegans papers in ACeDB. */
    AddReferences(Seq, nsp, the_entry);
    sprintf(gBuf, "The sequence of C. elegans cosmid %s", gSequence) ;
    AddFinishers(Seq, nsp, the_entry, gBuf);

    strcpy(gBuf,"Submitted by: ~");
    strcat(gBuf,"         Genome Sequencing Center ~");
    strcat(gBuf,"         Department of Genetics, Washington University, ~");
    strcat(gBuf,"         St. Louis, MO 63110, USA, and ~");
    strcat(gBuf,"         Sanger Centre, Hinxton Hall~");
    strcat(gBuf,"         Cambridge CB10 IRQ, England ~");
    strcat(gBuf,"         e-mail: rw@nematode.wustl.edu and jes@sanger.ac.uk~");

    strcat(gBuf,"~~NOTICE:  This sequence may not be the entire insert of this clone.~");
    strcat(gBuf,"It may be shorter because we only sequence overlapping sections~");
    strcat(gBuf,"once, or longer because we provide a small overlap between ~");
    strcat(gBuf,"neighboring submissions.~");

    strcat(gBuf,"~This sequence was finished as follows unless otherwise noted:~");
    strcat(gBuf,"all regions were double stranded or sequenced with an alternate ~");
    strcat(gBuf,"chemistry; an attempt was made to resolve all sequencing problems,~");
    strcat(gBuf,"such as compressions and repeats; all regions were covered by ~");
    strcat(gBuf,"sequence from more than one subclone ~");

    if (bsFindTag(Seq, _DB_remark) && bsFlatten (Seq, 1, a)) {
      strcat(gBuf, "~~            NEIGHBORING COSMID INFORMATION:~");
      for (i = 0; i < arrayMax(a) ; ++i) {
	sprintf(gBuf+strlen(gBuf), "~~%s",name(arr(a,i,BSunit).k));
      }
    }



    if (isCDS) {
      if (gAsnDebug) printf("  Adding comment about GeneFinder predictions.\n");
      sprintf(gBuf+strlen(gBuf), " ~ ~NOTES:~~Coding sequences below are predicted from");
      sprintf(gBuf+strlen(gBuf), " computer analysis, using the program Genefinder");
      sprintf(gBuf+strlen(gBuf), "(P. Green and L. Hillier, ms in preparation).");
    }

    if (bsFindTag (Seq, _TSL_site)) {
      sprintf(gBuf+strlen(gBuf), "~\"TSL\" = trans-spliced leader.");
    }

    if (strlen(gBuf) > (size_t) 0) AddCommentToEntry(nsp, nuc_entry, gBuf);

    break;

  case HUMAN:

    AddOrganismToEntry(nsp, the_entry, "Homo sapiens", "H. sapiens",
		       NULL, NULL, NULL, NULL, NULL );

    sfp = FeatureBuild(nsp, nuc_entry, FALSE, FALSE, FALSE, NULL);
    AddIntervalToFeature(nsp,sfp,nuc_entry,NULL,1,gLength,TRUE,FALSE,FALSE);
    MakeImpFeature (nsp,sfp,"source");
    AddQualToImpFeature(nsp,sfp,"organism","Homo sapiens");
    AddQualToImpFeature(nsp,sfp,"chromosome",gChromosome);
    AddQualToImpFeature(nsp,sfp,"clone",gRawClone);
    AddQualToImpFeature(nsp,sfp,"clone_lib",gLibrary);
    AddQualToImpFeature(nsp,sfp,"map",gMapPosition);

    gbp = StrMove(gBuf,"Eukaryota; Animalia; Chordata; Vertebrata; Mammalia; Theria; ");
    gbp = StrMove(gbp,"Eutheria; Primates; Haplorhini; Catarrhini; Hominidae.");

    /* keywords go in last three fields */

    /*    AddGenBankBlockToEntry(nsp, the_entry, gBuf, "PRI", NULL, NULL, NULL);*/

    if (ii==0)
        AddGenBankBlockToEntry(nsp, the_entry, gBuf, "PLN",NULL,NULL,NULL);
    else if (ii==1)
      AddGenBankBlockToEntry(nsp, the_entry, gBuf, "PLN", keyword[0],NULL,NULL);
    else if (ii==2)
      AddGenBankBlockToEntry(nsp, the_entry, gBuf, "PLN", keyword[0], keyword[1], NULL);
    else if (ii==3)
      AddGenBankBlockToEntry(nsp, the_entry, gBuf, "PLN", keyword[0], keyword[1], keyword[2]);

	gbp = gBuf;
    gbp += sprintf(gBuf, "Human %s clone %s", gCloneType, gClone);
	if (gMapPosition[0] != '\0') {
    	gbp = StrMove(gbp, " from ");
		gbp = StrMove(gbp, gMapPosition);
	}
    AddTitleToEntry(nsp, nuc_entry, gBuf);

    AddReferences(Seq, nsp, the_entry);

    sprintf(gBuf, "The sequence of H. sapiens %s clone %s", gCloneType, gClone);
    AddFinishers(Seq, nsp, the_entry, gBuf);

    strcpy(gbp = gBuf,"Submitted by:~");
    gbp = StrMove(gbp,"   Genome Sequencing Center~");
    gbp = StrMove(gbp,"   Department of Genetics, Washington University~");
    gbp = StrMove(gbp,"   St. Louis, MO 63108, USA~");
    gbp = StrMove(gbp,"   e-mail: sapiens@watson.wustl.edu~");
    gbp = StrMove(gbp,"~NOTICE:  This sequence may not represent the entire insert of this clone.  ");
    gbp = StrMove(gbp,"It may be shorter because we only sequence overlapping sections ");
    gbp = StrMove(gbp,"once, or longer because we provide a small overlap between ");
    gbp = StrMove(gbp,"neighboring submissions.");
    gbp = StrMove(gbp,"~~This sequence was finished as follows unless otherwise noted:~");
    gbp = StrMove(gbp,"all regions were double stranded or sequenced with an alternate ");
    gbp = StrMove(gbp,"chemistry; an attempt was made to resolve all sequencing problems, " );
    gbp = StrMove(gbp,"such as compressions and repeats; all regions were covered by ");
    gbp = StrMove(gbp,"sequence from more than one subclone; and the assembly was ");
    gbp = StrMove(gbp,"confirmed by restriction digest.");

    if (strcmp(gLibrary,"LL0XNCC01-U")==0) {

      gbp = StrMove(gbp,"~~SOURCE INFORMATION: This clone is from a chromosome X specific cosmid library LL0XNCC01 \"U\".  ");
      gbp = StrMove(gbp,"The source of the chromosomes was a human/hamster hybrid, GM07297-F, from ");
      gbp = StrMove(gbp,"Robert Nussbaum at University of Pennsylvania School of Medicine.");
	  gbp = StrMove(gbp, "~VECTOR:  Lawrist16");

	gbp = StrMove(gbp,"~~Clone reference:~");
	gbp = StrMove(gbp,"M. Grieff, R. Mazzarella, M.P. Whyte, R.V. Thakker, R. Wilson, S. Chissoe, and D. Schlessinger.  ");
	gbp = StrMove(gbp,"\"X-linked hypophosphatemia candidate gene region: sequence data on 80 kb ");
	gbp = StrMove(gbp,"containing spermine synthase and the 5' region of PEX\", in preparation.");

    }

	if (strcmp(gLibrary,"LL22NC03")==0) {
		gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
		gbp = StrMove(gbp,
		"This clone is from the human chromosome 22-specific cosmid library ");
		gbp = StrMove(gbp,
		"LL22NCO3, constructed at the Biomedical Sciences Division, Lawrence ");
		gbp = StrMove(gbp,
		"Livermore National Laboratory, Livermore, CA 94550 under the ");
		gbp = StrMove(gbp,
		"auspices of the National Laboratory Gene Library Project sponsored ");
		gbp = StrMove(gbp,
		"by the US Department of Energy.  The source of the flow sorted ");
		gbp = StrMove(gbp,
		"chromosomes was a human/hamster hybrid containing chromosomes Y, 22 and 9.  ");
		gbp = StrMove(gbp,
		"This clone is part of a cosmid contig isolated using YACs from the ");
		gbp = StrMove(gbp,
		"Sanger Centre chromosome 22 YAC contig described by ");
		gbp = StrMove(gbp,
		"Collins et al., Nature 377 Suppl., 367-379.");
		gbp = StrMove(gbp,"~VECTOR:  Lawrist16");
	}

    if (strcmp(gLibrary,"LLNL3")==0) {
	gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
        gbp = StrMove(gbp,"This clone is from a chromosome 3 specific library described by ");
        gbp = StrMove(gbp,"Wei et al., Cancer Research 56:1487-1492 (1996).");
		gbp = StrMove(gbp, "~VECTOR:  pWE15");
    }

    if (strcmp(gLibrary,"Stratagene_placenta_951202")==0) {
	gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
        gbp = StrMove(gbp,"This human genomic clone is from a male, caucasian placental cosmid library purchased from Stratagene.");
		gbp = StrMove(gbp, "~VECTOR:  pWE15");

    }

    if (strcmp(gLibrary,"GSBAC1")==0) {
	gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
	gbp = StrMove(gbp,"This clone is from Genome Systems first BAC library.");
	gbp = StrMove(gbp,"~Cell line:  lymphoblastoid");
	gbp = StrMove(gbp,"~Haplotypes:  two");
        gbp = StrMove(gbp,"~VECTOR:  pBELO");
        gbp = StrMove(gbp,"~Selection:  chloramphenicol");
    }

    if (strcmp(gLibrary,"RPCI-1")==0 || strcmp(gLibrary,"RPCI-3")==0 || strcmp(gLibrary,"RPCI-4")==0 || strcmp(gLibrary,"RPCI-5")==0 || strcmp(gLibrary,"RPCI-6")==0) {
		gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
		gbp = StrMove(gbp,"The sequence for this clone was derived from a human PAC library ");
		gbp = StrMove(gbp,"prepared by Pieter de Jong and coworkers at Roswell Park ");
		gbp = StrMove(gbp,"Cancer Institute using the method described by Ioannou et al., ");
		gbp = StrMove(gbp,"Nature Genetics 6:84-89 (1994).");

		if (strcmp(gLibrary,"RPCI-6")==0) {
			gbp = StrMove(gbp,"  The library is from one female donor.");
			gbp = StrMove(gbp, "~VECTOR:  pPAC4");
		}
		else {
			gbp = StrMove(gbp,"  The library is from one male donor.");
			gbp = StrMove(gbp, "~VECTOR:  pCYPAC2");
		}
    }


    if (strcmp(gLibrary,"CITB-978SK-B")==0) {
	gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
	gbp = StrMove(gbp,"This clone is from the first release of the human BAC library.  The library contains cloned DNA ");
	gbp = StrMove(gbp,"from a human male fibroblast cell line 978SK.  For references see:  ");

        gbp = StrMove(gbp,"Shizuya et al., Proc. Natl. Acad. Sci. 89:8794-8797 (1992); ");

        gbp = StrMove(gbp,"Kim et al., Genomics 34:213-218 (1996).");

	gbp = StrMove(gbp,"~VECTOR:  pBELO~Selection:  chloramphenicol");

    }

    if (strcmp(gLibrary,"CITB-HS-A")==0) {
		gbp = StrMove(gbp,"~~SOURCE INFORMATION:~");
	gbp = StrMove(gbp,"This clone is from a release of the human BAC library. The library contains cloned DNA ");
	gbp = StrMove(gbp,"from a human sperm.  For references see:  ");

        gbp = StrMove(gbp,"Shizuya et al., Proc. Natl. Acad. Sci. 89:8794-8797 (1992); ");

        gbp = StrMove(gbp,"Kim et al., Genomics 34:213-218 (1996).");

	gbp = StrMove(gbp,"~VECTOR:  pBELO~Selection:  chloramphenicol");
    }


    if (bsFindTag(Seq, _DB_remark) && bsFlatten (Seq, 1, a)) {
      gbp = StrMove(gbp, "~~NEIGHBORING SEQUENCE INFORMATION:");
      for (i = 0; i < arrayMax(a) ; ++i) {
			gbp += sprintf(gbp, "~%s", name(arr(a,i,BSunit).k));
			if (i + 1 < arrayMax(a) )
				gbp = StrMove(gbp, "~");
      }
    }
    if (gBuf[0] != NULLB)
		AddCommentToEntry(nsp, nuc_entry, gBuf);

    break;

  default:

    message = "  Project not implemented yet.";
    goto end_asnDoDump;

  }


  /* Now dump EMBL_feature, Subsequence, and Feature (human only). */
  /* Recurse through parent. */

  key = seq;
  obj = bsCreate(key);
  gOffset = 1;
  gOffForward = TRUE;

  while (TRUE) {
    KEY new;
    int i,x1=0,x2=0;

    if (gAsnDebug) {
      printf("  Analyzing subseq/features from %s.\n",name(key));
    }

    if (bsFindTag(obj, _EMBL_feature) && bsFlatten(obj, 4, a)) {
      for (i = 0; i < arrayMax(a); i += 4) {
	status = emblFeature(arr(a,i+0,BSunit).k,
			     arr(a,i+1,BSunit).i,
			     arr(a,i+2,BSunit).i,
			     arr(a,i+3,BSunit));
      }
      if (status) {
	bsDestroy(obj);
	goto end_asnDoDump;
      }
    }

    if (bsFindTag(obj, _Subsequence) && bsFlatten(obj, 3, a)) {
      for (i = 0; i < arrayMax(a); i += 3) {
	status = subseqFeatures(arr(a,i+0,BSunit).k,
				arr(a,i+1,BSunit).i,
				arr(a,i+2,BSunit).i);
      }
      if (status) {
	bsDestroy(obj);
	goto end_asnDoDump;
      }
    }


    if (gProject == HUMAN) {
      status = dumpGBFeatures(key);
      status = dumpMiscFeatures(key);
      if (status) {
	bsDestroy(obj);
	goto end_asnDoDump;
      }
    }

    if (gProject == HUMAN) {
      status = dumpMiscFeatures(key);
      if (status) {
	bsDestroy(obj);
	goto end_asnDoDump;
      }
    }

    if (!bsGetKey(obj, _Source, &new)) {
      bsDestroy (obj);
      break;
    }

    if (gAsnDebug) {
      printf("  Using %s found _Source %s.\n",name(key),name(new));
    }

    bsDestroy(obj);

    if (!(obj = bsCreate (new))) {
      break;
    }

    if (!bsFindKey (obj, _Subsequence, key) ||
	!bsGetData (obj, _bsRight, _Int, &x1) ||
	!bsGetData (obj, _bsRight, _Int, &x2)) {
      break;
    }

    if (gAsnDebug) {
      printf("  New key %s has subseq with left/right limits.\n",name(new));
    }

    if (x2 > x1)
	gOffset += x1 - 1;
    else {
      gOffForward = !gOffForward;
      gOffset = x1 + 1 - gOffset;
    }

    key = new;
  }

 end_asnDoDump:

  if (message) printf("  Fatal error dumping %s: %s\n", name(seq), message);

  bsDestroy(Seq);
  arrayDestroy(dna);

  return (message) ? 1 : 0;
}


/*********************************************************************
 * asnDumpKey:
 *
 * Called from displayBlock (pointer to function passed from asnDump).
 *
 * Process:
 *         [0] Check for valid ACEDB_PROJECT;
 *         [1] Check for valid ASN output file;
 *         [2] Get submission date and call NCBISubCreate;
 *         [3] Add author for NCBI reference;
 *         [4] Call asnDoDump to do the rest;
 *         [5] Validate NCBI submission and report errors.
 ********************************************************************/

void asnDumpKey (KEY key)
{
	struct tm	*date = NULL;	/* MDH 11-14-94 for submit date.             */
	time_t	clock;			/* MDH 11-14-94 for submit date.             */
	BOOL	PubHold = FALSE;	/* Do not hold submission until publication. */
	PubPtr	pub;			/* Pointer to NCBI publication feature.      */
	OBJ		Seq;			/* Object for key seq.                       */
	char	asnfile[1024];	/* Absolute full file name for ASN.1 output. */
	char	dname[DIR_BUFFER_SIZE] = {0};	/* Path name of ASN.1 output.*/
	char	fname[FIL_BUFFER_SIZE] = {0};	/* ASN.1 output file (no extension). */
	char	*ptype;			/* Pointer to project name (from getenv).    */
	int		i;				/* Loop variable.                            */
	int		asnauto = 0;	/* Is set to 1 if ASNAUTO defined in env.    */
	int		asncase = 0;	/* Is set to 1 if ASNCASE defined in env.    */
	int		asnprompt = 0;	/* Is set to 1 if ASNPROMPT defined in env.  */
	Int2	nasnerrors;		/* Number of errors in the ASN.1 submission  */
	FILE	*fil = NULL;	/* Temp file pointer to test file open.      */
	char	*sub_email;

	printf("Begin ASN.1 dump for sequence %s.\n",name(key));

	/* Determine if debug output on.  */
	gAsnDebug = (getenv("ASNDEBUG") != NULL);

	/* Determine if output file determined from sequence name.  */
	asnauto = (getenv("ASNAUTO") != NULL);

	/* Determine if output file should be lower case.  */
	asncase = (asnauto && getenv("ASNCASE") != NULL);

	/* Determine if prompt before overwriting ASN file. */
	asnprompt = (asnauto && getenv("ASNPROMPT") != NULL);

	/* Determine project type from ACEDB_PROJECT; return if invalid project */
	ptype = getenv("ACEDB_PROJECT");

  if (ptype == NULL)
    gProject = NOPROJECT;
  else if (!strcmp(ptype,"elegans"))
    gProject = ELEGANS;
  else if (!strcmp(ptype,"yeast"))
    gProject = YEAST;
  else if (!strcmp(ptype,"human"))
    gProject = HUMAN;
  else if (!strcmp(ptype,"briggsae"))
    gProject = BRIGGSAE;
  else if (!strcmp(ptype,"fugu"))
    gProject = FUGU;
  else
    gProject = NOPROJECT;

  if (gProject == NOPROJECT) {
    printf("  Fatal error: Must set ACEDB_PROJECT to valid project name.\n");
    printf("Abort ASN.1 dump for sequence %s.\n",name(key));
    return;
  }

  /* Get dname/fname.asn and create file name string for NCBI toolkit. */

  if (asnauto) {
    sprintf(asnfile, "%s.asn", name(key));
    if (asncase > 0) {
      for (i=0; i<strlen(asnfile); i++) {
		if (isupper(asnfile[i]))
			asnfile[i] = tolower(asnfile[i]);
      }
    }
    if (asnprompt) {
        if (checkFileExists(asnfile)) {
	  printf("Abort ASN.1 dump for sequence %s.\n",name(key));
	  return;
	}
    }
  }
  else {
    if ( !(fil = filqueryopen(dname,fname,"asn","w","NCBI Submission file")) ) {
      printf("  Fatal error: File open failed for %s/%s.asn.\n",dname,fname);
      printf("Abort ASN.1 dump for sequence %s.\n",name(key));
      return;
    }
    fclose(fil);
    sprintf(asnfile, "%s/%s.asn", dname, fname);
  }

  /* Create submit date. MDH 11-14-94 */

  clock = time(NULL);
  date  = localtime(&clock);
  gDay   = date->tm_mday;         /* days of month (1-31) */
  gMon   = 1 + date->tm_mon;      /* months since jan (0-11) */
  gYear  = 1900 + date->tm_year;  /* years since 1900 */

 	/* Create submission. */

	switch (gProject) {
	case HUMAN:
		sub_email = "humansub@watson.wustl.edu";
		break;
	default:
		sub_email = "lhillier@watson.wustl.edu";
		break;
	}

	nsp = NCBISubCreate("Hillier","LaDeana",NULL,"L.",NULL,
		      "Washington University","Department of Genetics",
		      "4444 Forest Park Avenue","St. Louis","Missouri",
		      "USA","63108","314-286-1800","314-286-1810",
		      sub_email,PubHold,gMon,gDay,gYear);

	if (nsp == NULL) {
		printf("  Fatal Error: NCBISubCreate returned NULL pointer.\n");
		printf("Abort ASN.1 dump for sequence %s.\n",name(key));
		return;
	}

	switch (gProject) {
	case HUMAN:
		DefineSubmittorKey(nsp, Institution);
		AddTechToEntry(nsp, the_entry, MI_TECH_htgs_3);
		break;
	default:
		break;
	}

	pub = CitSubBuild(nsp, gMon, gDay, gYear, MEDIUM_EMAIL);
	CitSubForSubmission(nsp, pub);
	AddAuthorToPub(nsp, pub, "Waterston", "Robert", NULL, "R.", NULL);

	if (asnDoDump(key)) {
		printf("Abort ASN.1 dump for sequence %s.\n",name(key));
		return;
	}

	/* Validate submission and report errors; write ASN.1 file. */
	if ((nasnerrors = NCBISubValidate(nsp, stderr)) == 0)
		printf("  NCBI validate:  SUCCESS!\n");
	else
		printf("  NCBI validate:  %d ERRORS in ASN.1 submission.\n",
				(int)nasnerrors);

	NCBISubWrite(nsp, asnfile);
	NCBISubFree(nsp);

	printf("End ASN.1 dump for sequence %s.\n",name(key));
}

/*********************************************************************
 * asnDump:
 *
 * Entry point for ASN.1 dump; calling routine in w7/fmapcontrol.c
 ********************************************************************/

void asnDump (void)
{
  displayBlock(asnDumpKey,"I will NCBI dump the corresponding sequence ") ;
}

/*********************************************************************
 * catnames:
 *
 * Return a malloced string that is the concatenation of name1 & name2
 ********************************************************************/
char *
catnames (char *name1, char *name2)
{
	char	*cp, *cp0;
	size_t	len1, len2;
	size_t	len;

	name1 = (name1 == NULL ? "" : name1);
	name2 = (name2 == NULL ? "" : name2);

	len1 = strlen(name1);
	len2 = strlen(name2);

	len = len1 + len2 + 3; /* make room for "; " + NUL byte */
	cp0 = cp = Nlm_MemNew(len);

	cp = Nlm_StrMove(cp, name1);
	if (*name1 != NULLB && *name2 != NULLB)
		cp = Nlm_StrMove(cp, "; ");

	MemCpy(cp, name2, len2 + 1);

	return cp0;
}

/****************************************************************************
 * get_acc_nid:
 *
 * Return the NID associated with the given accession
 ****************************************************************************/
static unsigned long
get_acc_nid(char *fname, char *acc)
{
	FILE	*pp;
	char	*cp;
	char	buf[4096];
	unsigned long	nid = 0;

	sprintf(buf, "getfa %s %s", fname, acc);

	if ((pp = popen(buf, "r")) == NULL)
		return nid;

	if (fgets(buf, sizeof buf, pp) == NULL || buf[0] != '>')
		goto Return;

	if ((cp = StringStr(buf, "gi|")) == NULL)
		goto Return;

	sscanf(cp + 3, "%lu", &nid);

Return:
	pclose(pp);
	return nid;
}

/****************************************************************************
 * get_cDNA_NID:
 *
 * Return the NID associated with the given accession
 * -- looking in PRImate and EST divisions only.
 ****************************************************************************/
unsigned long
get_cDNA_NID(char *acc)
{
	unsigned long	nid;

	if ((nid = get_acc_nid("pri.nt", acc)) != 0)
		return nid;

	if ((nid = get_acc_nid("est.nt", acc)) != 0)
		return nid;

	return nid;
}


/****************************************************************************

    Commented out for now since info is in DB_remark; at
    least for C. elegans...This goes in the build of the
    main comment; do not delete. Save for future reference.
    Note that this code has never been tested.

    if (clone && (obj = bsCreate(clone))) {
      KEY sseq;
      OBJ Sseq;
      if (bsFindKey(Seq, _Clone_left_end, clone) && bsGetData(Seq, _bsRight, _Int, &i))
        sprintf(gBuf+strlen(gBuf),"  The true left end of clone %s is at %d in this sequence.~",
		name(clone), i);
      else if (bsGetKey(obj, _Clone_left_end, &sseq) && (Sseq = bsCreate (sseq))) {
	if (bsFindKey(Sseq, _Clone_left_end, clone) && bsGetData(Sseq, _bsRight, _Int, &i))
	  sprintf(gBuf+strlen(gBuf),"  The true left end of clone %s is at %d in sequence CEL%s~",
		  name(clone), i, name(sseq));
	bsDestroy(Sseq);
      }
      if (bsFindKey(Seq, _Clone_right_end, clone) && bsGetData(Seq, _bsRight, _Int, &i))
        sprintf(gBuf+strlen(gBuf),"  The true right end of clone %s is at %d in this sequence.~",
	      name(clone), i);
      else if (bsGetKey(obj, _Clone_right_end, &sseq) && (Sseq == bsCreate (sseq))) {
	if (bsFindKey(Sseq, _Clone_right_end, clone) && bsGetData(Sseq, _bsRight, _Int, &i))
	sprintf(gBuf+strlen(gBuf),"  The true right end of clone %s is at %d in sequence CEL%s~",
	       name(clone), i, name(sseq));
	bsDestroy (Sseq);
      }
      if (bsFindTag(Seq, _Clone_left_end) && bsFlatten(Seq, 2, a))
        for (i = 0 ; i < arrayMax(a) ; i += 2)
	  if (arr(a, i, BSunit).k != clone && arr(a, i+1, BSunit).i)
	    sprintf(gBuf+strlen(gBuf),"  The true left end of clone %s is at %d in this sequence.~",
		    name(arr(a, i, BSunit).k), arr(a, i+1, BSunit).i);
      if (bsFindTag(Seq, _Clone_right_end) && bsFlatten(Seq, 2, a))
        for (i = 0 ; i < arrayMax(a) ; i += 2)
          if (arr(a, i, BSunit).k != clone && arr(a, i+1, BSunit).i)
            sprintf(gBuf+strlen(gBuf),"  The true right end of clone %s is at %d in this sequence.~",
		    name(arr(a, i, BSunit).k), arr(a, i+1, BSunit).i);
      bsDestroy(obj);
    }
    if (bsGetKey(Seq, _Overlap_left, &key)) {
      strcat(gBuf, "The start of this sequence ");
      if (obj = bsCreate(key)) {
	int len, off;
	KEY dna;
	if (bsGetKey (obj, _DNA, &dna) &&
	    bsGetData (obj, _bsRight, _Int, &len) &&
	    bsFindKey (obj, _Overlap_right, seq) &&
	    bsGetData (obj, _bsRight, _Int, &off))
	sprintf(gBuf+strlen(gBuf), "(1..%d) ", len - off + 1);
	bsDestroy (obj);
      }
      sprintf(gBuf+strlen(gBuf),"overlaps with the end of sequence CEL%s.~",name (key));
    }
    if (bsGetKey(Seq, _Overlap_right, &key)) {
      strcat (gBuf, "The end of this sequence ");
      if (bsGetData(Seq, _bsRight, _Int, &i))
      sprintf(gBuf+strlen(gBuf), "(%d..%d) ", i, gLength);
      sprintf(gBuf+strlen(gBuf), "overlaps with the start of sequence CEL%s.~",name (key));
    }

*****************************************************************************/

 

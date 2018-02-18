/*  File: UCSCtrackCompare.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2004
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file was developped at the occasion of the EGASP/Gencode project
 * Please quote:
 *   D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov) 
 *     Genomewide AceView annotation...
 *     pubmedcentral 2006
 *
 *  This code should be linked with the acedb libraries
 *  available from http://www.aceview.org
 *        ftp://ftp.ncbi.nlm.nih.gov/repository/acedb/AceView
 *  The present code and the acedb schema for this project
 *  are available from mieg@ncbi.nlm.nih.gov
 */

#include "../wac/ac.h"
#include "keyset.h"
#include <errno.h>
#include <math.h>
#include "bitset.h"
#include "freeout.h"
#include "vtxt.h"
#include "dict.h"
#include "aceio.h"

#define MAXEXON 128

static void usage (void) ;

typedef enum {ZERO=0, HAV, HAV_no_CDS, HPV, ACV, DOG, ENS, EXO, EXH, FGH, GID, GNM, JGS, PAA, PAM, PAU, PAN,
	     SG2, SPD, TWS, PAC, PAO, AGI, AGA, AGD, AGE, GNZ, SAG,
	     KGP, CCP, RGP, MGP, VGP, EGP, AVP, ECP, TWP, NSP, NEP, SGP, 
	     GIP, GSP, EXP, AUP, EWP, KKG, XXX} LABID ;

typedef struct mmStruct {
  int gene ;
  KEY mrna ;
  KEY chrom ;
  int labId ;
  int z1, m1, m2, c1, c2, quality, score, twin ;
  int nExons, nIntrons ;
  BOOL hasCds, isDown, atgComplete, stopComplete ;
  int regionName ;
  int a1 [MAXEXON ] ;  
  int a2 [MAXEXON ] ;  
  int a1_original [MAXEXON ] ;  
  int a2_original [MAXEXON ] ;  
} MM ;


typedef struct labostruct { 
  int labId ;
  BOOL hasStop, hasBinSection ;
  int filName ;
  int nbp, nvp, nfpi, nfpnoi, nfn, badStop, badAtg, badAtgAndStop, badElseWhere, badDenom ;
  int nex, extp, exfpKnown, exfpNew, exfpNew1, exfpCds, exfn ;
  int nTranscript,  nIntronlessTranscript, nExon, nTerminalExon ;
  int nExactTranscript, nSimilarTranscript, nNewVariant, nInNewGene, nMissedVariant, nInMissedGene ;
  int nTrueIntron, nFalseIntron, nTrueIntronInNewVariant, nFalseIntronInNewVariant, nFalseIntronInNewGene ;
  int nMissedIntronInMissedGene, nMissedIntronInTouchedGene ;
  int nSimilarIntronlessVariant, nNewIntronlessVariant, nIntronlessInNewGene ;

  BitSet bbDi, bbRi, bbDnoi, bbRnoi ;
  Array exons ;
} LAB ;

typedef struct utcStruct {
  AC_HANDLE h ;
  DICT *dict ;
  DICT *chromDict ; 
  DICT *geneDict ;
  char *template ;
  char *keysetName ;
  BOOL cds, bp, introns, exons, regions, global, multi, filter, nosort, pairs, aceExport ;
  int makeConcensus ;
  const char *chromConfigFile ;
  const char *trackConfigFile ;
  Array aa, chromArray, labs ;
  LAB *refLab ;
} UTC ;

static int REF = 0 ;
static int ACVV = 0 ;

LAB labosNoCDS[1] ;

#ifdef JUNK

typedef struct laboOldstruct { 
  LABID labId ;
  BOOL hasStop, hasBinSection ;
  char colour[64] ;
  char labName[32] ; 
  char filName[64] ;
  int nbp, nvp, nfpi, nfpnoi, nfn, badStop, badAtg, badAtgAndStop, badElseWhere, badDenom ;
  int nex, extp, exfpKnown, exfpNew, exfpNew1, exfpCds, exfn ;
  int nTranscript,  nIntronlessTranscript, nExon, nTerminalExon ;
  int nExactTranscript, nSimilarTranscript, nNewVariant, nInNewGene, nMissedVariant, nInMissedGene ;
  int nTrueIntron, nFalseIntron, nTrueIntronInNewVariant, nFalseIntronInNewVariant, nFalseIntronInNewGene ;
  int nMissedIntronInMissedGene, nMissedIntronInTouchedGene ;
  int nSimilarIntronlessVariant, nNewIntronlessVariant, nIntronlessInNewGene ;

  BitSet bbDi, bbRi, bbDnoi, bbRnoi ;
  Array exons ;
} LABold ;

LABold labos [] = {
  /*
  { KKG, TRUE,  "BLUE",          "King-Kong " , "DATA_egasp/kingkong" } , 
    * kingkong == exhon-hunter + sgp2 + genemark + genzilla 
    */
  { HAV, TRUE,  "BLUE",          "Utc " , "DATA_egasp/encodeGencodeGeneKnownOrPutativeOct05" } , 
  { ACV, FALSE, "MAGENTA",       "AceView ", "DATA_egasp/encodeEgaspFullAceview" } , 
  /*
  { ZERO, TRUE, "", "", ""} ,
  */
  { DOG, FALSE, "lightblue",     "UP Dogfish",          "DATA_egasp/encodeEgaspFullDogfish" } ,
  { ENS, TRUE,  "paleyellow",    "Ensembl ",             "DATA_egasp/encodeEgaspFullEnsembl"  } ,
  { EXO, FALSE, "palegreen" ,    "Exogean " ,            "DATA_egasp/encodeEgaspFullExogean" } ,

  { EXH, TRUE,  "paleviolet" ,   "UP ExonHunter " ,     "DATA_egasp/encodeEgaspFullExonhunter" } ,
  { FGH, FALSE, "lightblue",     "Fgenesh ",             "DATA_egasp/encodeEgaspFullFgenesh"  } ,

  { GID, TRUE,  "yellow",        "UP GeneID" ,          "DATA_egasp/encodeEgaspFullGeneId" } ,
  { GNM, TRUE,  "orange",        "UP GeneMark",         "DATA_egasp/encodeEgaspFullGenemark" } , 
  { JGS, FALSE, "violet",        "UP Jigsaw",           "DATA_egasp/encodeEgaspFullJigsaw" } ,

  { PAA, FALSE, "lightcyan",     "Pairagon",            "DATA_egasp/encodeEgaspFullPairagonAny"  } ,

  { SG2, TRUE,  "yellow",        "UP SGP2 ",             "DATA_egasp/encodeEgaspFullSgp2" } ,
  { TWS, TRUE,  "lightgreen",    "P Twinscan" ,         "DATA_egasp/encodeEgaspFullTwinscan" } ,

  { AGA, FALSE, "paleorange" ,   "UP Augustus" ,        "DATA_egasp/encodeEgaspPartAugustusAny" } ,

  { GNZ, FALSE, "brown",         "UP GeneZilla",        "DATA_egasp/encodeEgaspPartGenezilla"} ,
  { SAG, FALSE, "brown",         "UP Saga ",             "DATA_egasp/encodeEgaspPartSaga" } ,
  /*
   * PUBLIC
  */
  { KGP, TRUE,  "lightblue",     "*KnownGene" ,         "DATA_public/knownGene" } ,
  { CCP, TRUE,  "darkblue",      "*P CCDS ",             "DATA_public/ccdsGene" } ,
  { RGP, TRUE,  "cyan",          "*RefSeq " ,            "DATA_public/refGene" } ,

  { MGP, TRUE,  "lightcyan",     "*MGC    " ,               "DATA_public/mgcGenes" } ,
  { EGP, TRUE,  "green",         "*Ensembl" ,           "DATA_public/ensGene" } ,
  { AVP, TRUE,  "cerise",        "*AceView",            "DATA_public/acembly"  } ,
  { ECP, FALSE, "purple",        "*ECgene " ,            "DATA_public/ECgene" } ,

  { NEP, TRUE,  "darkgray",      "*U NscanEst",         "DATA_public/nscanEstGene" } ,
  { GSP, TRUE,  "lightgreen",    "*UP GenScan" ,        "DATA_public/genscan" } ,
  { EWP, TRUE,  "midblue",       "*ExonWalk" ,          "DATA_public/exonWalk" } ,

  { ZERO, TRUE, "", "", ""} ,

  { ZERO, TRUE, "", "", ""} 
} ;

#endif
	
/* summary of table 4 5 6 of guigo et al for comparison to UCSCtrackCompare */

typedef struct guigoStruct { int labId ;  char labNam[128]; 
  float bpSn, bpSp, exonSn, exonSp ;
  float cdsBpSn, cdsBpSp, cdsExonSn, cdsExonSp, cdsTrSn, cdsTrSp ;
} GUIGO ;
GUIGO guigos[] = {
{ 0, "UP Augustus", 40.46, 84.50, 33.60, 54.80, 94.42, 82.43, 74.67, 76.76, 22.65, 35.59 }, /* Augustus-any */
{ 0, "Pairagon", 56.31, 89.36, 41.23, 74.93, 87.77, 92.78, 76.85, 88.91, 39.29, 60.34 }, /* Pairagon/NSCAN */
{ 0, "UP Jigsaw" , 40.02, 93.35, 36.10, 64.56, 94.56, 92.19, 80.61, 89.33, 34.05, 65.95 }, /* Jigsaw */
{ 0, "Fgenesh", 48.87, 81.16, 35.84, 58.41, 91.09, 76.89, 75.18, 69.31, 36.21, 41.61 }, /* fgenesh++ */
{ 0, "XXX", 33.73, 77.25, 23.76, 43.97, 78.65, 75.29, 52.39, 62.93, 11.09, 17.22 }, /* Augustus-abinit */
{ 0, "UP GeneMark", 33.09, 65.28, 22.74, 40.61, 78.43, 37.97, 50.58, 29.01, 6.93, 3.24 }, /* GenemarkHMM */
{ 0, "UP GeneZilla", 37.42, 52.08, 29.60, 36.58, 87.56, 50.93, 62.08, 50.25, 9.09, 8.84 }, /* Genezilla */
{ 0, "UP ExonHunter", 38.93, 61.44, 30.93, 29.46, 90.46, 59.67, 64.44, 41.77, 10.48, 6.33 }, /* Exonhunter */
{ 0, "Exogean", 60.58, 94.73, 48.87, 76.29, 84.18, 94.33, 79.34, 83.45, 42.53, 52.44 }, /* Exogean */
{ 0, "XXX", 39.65, 85.47, 33.38, 55.54, 92.62, 83.45, 74.10, 77.40, 22.50, 37.01 }, /* Augustus-EST */
{ 0, "AceView", 88.08, 79.47, 64.16, 61.18, 90.94, 79.14, 85.75, 56.98, 44.68, 19.31 }, /* Aceview */
{ 0, "PairagonNovel", 56.22, 89.35, 41.11, 74.98, 87.56, 92.77, 76.63, 88.95, 39.29, 60.64 }, /* Pairagon */
{ 0, "Ensembl", 61.61, 95.26, 41.61, 73.41, 90.18, 92.02, 77.53, 82.65, 39.75, 54.64 }, /* Ensembl */

{ 0, "XXX", 38.12, 82.26, 28.97, 49.71, 88.86, 80.15, 63.06, 69.14, 12.33, 18.64 }, /* Augustus-dual */
{ 0, "XXX", 39.55, 78.69, 32.41, 65.25, 85.38, 89.02, 67.66, 82.05, 16.95, 36.71 }, /* NSCAN */
{ 0, "UP Saga", 22.00, 81.54, 19.36, 39.55, 52.54, 81.39, 38.82, 50.73, 2.16, 3.44 }, /* Saga */

{ 0, "UP Dogfish", 27.49, 87.44, 26.36, 63.97, 64.81, 88.24, 53.11, 77.34, 5.08, 14.61 }, /* Dogfish */
{ 0, "P Twinscan", 36.16, 76.04, 30.48, 45.30, 84.25, 74.13, 65.56, 61.65, 15.87, 15.11 }, /* MARS */

{ 0, "UP GeneID", 31.80, 79.75, 24.41, 43.65, 75.03, 78.83, 51.41, 63.92, 5.39, 10.69 }, /* GeneID_u12 */
{ 0, "UP SGP2", 35.40, 68.04, 28.01, 29.77, 82.84, 66.80, 59.73, 49.20, 9.71, 8.44 }, /* SGP-u12 */
{ 0, "XXX", 18.94, 82.64, 16.96, 30.23, 35.99, 94.25, 29.81, 35.09, 0.00, 0.00 }, /* Spida */
{ 0, "XXX", 3.47, 97.84, 0.32, 6.32, 8.065, 95.77, 1.66, 27.22, 0.00, 0.00 }, /* Dogfish-exon */
{ 0, "XXX", 40.46, 84.50, 32.68, 51.58, 94.42, 82.43, 39.80, 40.89, 0.00, 0.00 }, /* Augustus-exon */
{ 0, "*ECgene", 93.00, 38.68, 58.17, 34.83, 96.36, 47.30, 86.22, 35.08, 56.86, 8.84 }, /* ECgene */
{ 0, "XXX", 91.94, 53.98, 65.51, 44.28, 96.43, 58.47, 84.66, 38.32, 33.90, 7.96 }, /* acembly AVP non comparable */
{ 0, "*P CCDS", 23.86, 99.68, 22.45, 76.80, 56.87, 99.52, 51.95, 97.75, 28.97, 85.58 }, /* ccdsGene */
{ 0, "*Ensembl", 62.43, 95.27, 41.71, 72.65, 91.39, 91.92, 77.71, 82.39, 40.52, 54.09 }, /* ensGene */
{ 0, "XXX", 32.56, 77.42, 25.75, 53.16, 76.77, 76.48, 53.84, 61.08, 4.78, 8.78 }, /* geneid */
{ 0, "*UP GenScan", 36.76, 63.20, 29.40, 42.31, 84.17, 60.60, 58.65, 46.37, 7.40, 10.13 }, /* genscan */
{ 0, "XXX", 65.74, 91.82, 43.84, 74.57, 89.10, 93.61, 78.11, 82.28, 43.45, 46.93 }, /* knownGene KGP non comparable */
{ 0, "XXX", 29.17, 96.73, 21.21, 74.10, 44.06, 97.56, 42.95, 93.61, 23.73, 78.24 }, /* mgcGenes MGP non comparable  */
{ 0, "XXX", 57.51, 97.07, 38.35, 83.51, 85.34, 98.50, 73.23, 94.67, 41.91, 75.21 }, /* refGene RGP non comparable */
{ 0, "XXX", 34.99, 82.86, 28.83, 56.44, 82.81, 82.20, 60.56, 65.16, 8.17, 12.59 }, /* sgpGene */
{ 0, "XXX", 33.40, 86.28, 27.79, 63.26, 78.16, 84.59, 58.43, 73.11, 10.63, 20.25 }, /* twinscan */
{ 0, "", 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  }
} ;


typedef struct badStruct {char *chromNam ; int chrom,  a1,  a2 ; int regionName ;} REGION ;
#ifdef JUNK
REGION baddies[] = {
  {"chr22", 0,      30128508  ,   31828507, "ENm004"}, 
  {"chrX",  0,     152635145  ,  153973591, "ENm006"}, 
  {"chr13", 0,      29418016  ,   29918015, "ENr111"}, 
  {"chr10", 0,      55153819  ,   55653818, "ENr114"}, 
  {"chr13", 0,     112338065  ,  112838064, "ENr132"}, 
  {"chr6", 0,      132218540  ,  132718539, "ENr222"}, 
  {"chr6", 0,       73789953  ,   74289952, "ENr223"}, 
  {"chr1", 0,      147971134  ,  148471133, "ENr231"}, 
  {"chr9", 0,      128764856  ,  129264855, "ENr232"}, 
  {"chr6", 0,      108371397  ,  108871396, "ENr323"}, 
  {"chrX",  0,     122507850  ,  123007849, "ENr324"}, 
  {"chr20", 0,      33304929  ,   33804928, "ENr333"}, 
  {"chr6", 0,       41405895  ,   41905894, "ENr334"}, 
  { 0, 0, 0, 0, 0}
}  ;

REGION gooddies[] = {
 /* 
  {"chr18", 0,   59412301, 59912300 , "ENr122"}, 
  { 0, 0, 0, 0, 0} ,

  {"chr16", 0,          1,    500000, "ENm008"}, 
  { 0, 0, 0, 0, 0} ,

  {"chr2", 0,   234273825, 234773888, "ENr131"}, 
  { 0, 0, 0, 0, 0} ,
  {"chr2", 0,   234733888, 234773888, "ENr131"}, 
  { 0, 0, 0, 0, 0} ,
  */
  {"chr7", 0,   115404472, 117281897, "ENm001"}, 
  {"chr5", 0,   131284314, 132284313, "ENm002"}, 
  {"chr11", 0,  115962316, 116462315, "ENm003"}, 
  {"chr21", 0,   32668237,  34364221, "ENm005"}, 
  {"chr19", 0,   59023585,  60024460, "ENm007"}, 
  {"chr16", 0,          1,    500000, "ENm008"}, 
  {"chr11", 0,    4730996,   5732587, "ENm009"}, 
  {"chr7", 0,    26730761,  27230760, "ENm010"}, 
  {"chr11", 0,    1699992,   2306039, "ENm011"}, 
  {"chr7", 0,   113527084, 114527083, "ENm012"}, 
  {"chr7", 0,    89428340,  90542763, "ENm013"}, 
  {"chr7", 0,   125672607, 126835803, "ENm014"}, 
  {"chr2", 0,    51570356,  52070355, "ENr112"}, 
  {"chr4", 0,   118604259, 119104258, "ENr113"}, 
  {"chr2", 0,   118010804, 118510803, "ENr121"}, 
  {"chr18", 0,   59412301,  59912300, "ENr122"}, 
  {"chr12", 0,   38626477,  39126476, "ENr123"}, 
  {"chr2", 0,   234273825, 234773888, "ENr131"}, 
  {"chr21", 0,   39244467,  39744466, "ENr133"}, 
  {"chr16", 0,   25780428,  26280428, "ENr211"}, 
  {"chr5", 0,   141880151, 142380150, "ENr212"}, 
  {"chr18", 0,   23719232,  24219231, "ENr213"}, 
  {"chr5", 0,    55871007,  56371006, "ENr221"}, 
  {"chr15", 0,   41520089,  42020088, "ENr233"}, 
  {"chr14", 0,   52947076,  53447075, "ENr311"}, 
  {"chr11", 0,  130604798, 131104797, "ENr312"}, 
  {"chr16", 0,   60833950,  61333949, "ENr313"}, 
  {"chr8", 0,   118882221, 119382220, "ENr321"}, 
  {"chr14", 0,   98458224,  98958223, "ENr322"}, 
  {"chr2", 0,   220102851, 220602850, "ENr331"}, 
  {"chr11", 0,   63940889,  64440888, "ENr332"}, 
  { 0, 0, 0, 0, 0}
}  ;

REGION gooddies2[] = {
  {"chr7", 0,   115404472, 117281897, "ENm001"}, 
  {"chr5", 0,   131284314, 132284313, "ENm002"}, 
  {"chr11", 0,  115962316, 116462315, "ENm003"}, 
  {"chr21", 0,   32668237,  34364221, "ENm005"}, 
  {"chr19", 0,   59023585,  60024460, "ENm007"}, 
  {"chr16", 0,          1,    500000, "ENm008"}, 
  {"chr11", 0,    4730996,   5732587, "ENm009"}, 
  {"chr7", 0,    26730761,  27230760, "ENm010"}, 
  {"chr11", 0,    1699992,   2306039, "ENm011"}, 
  {"chr7", 0,   113527084, 114527083, "ENm012"}, 
  {"chr7", 0,    89428340,  90542763, "ENm013"}, 
  {"chr7", 0,   125672607, 126835803, "ENm014"}, 
  {"chr2", 0,    51570356,  52070355, "ENr112"}, 
  {"chr4", 0,   118604259, 119104258, "ENr113"}, 
  {"chr2", 0,   118010804, 118510803, "ENr121"}, 
  {"chr18", 0,   59412301,  59912300, "ENr122"}, 
  {"chr12", 0,   38626477,  39126476, "ENr123"}, 
  {"chr2", 0,   234273825, 234773888, "ENr131"}, 
  {"chr21", 0,   39244467,  39744466, "ENr133"}, 
  {"chr16", 0,   25780428,  26280428, "ENr211"}, 
  {"chr5", 0,   141880151, 142380150, "ENr212"}, 
  {"chr18", 0,   23719232,  24219231, "ENr213"}, 
  {"chr5", 0,    55871007,  56371006, "ENr221"}, 
  {"chr15", 0,   41520089,  42020088, "ENr233"}, 
  {"chr14", 0,   52947076,  53447075, "ENr311"}, 
  {"chr11", 0,  130604798, 131104797, "ENr312"}, 
  {"chr16", 0,   60833950,  61333949, "ENr313"}, 
  {"chr8", 0,   118882221, 119382220, "ENr321"}, 
  {"chr14", 0,   98458224,  98958223, "ENr322"}, 
  {"chr2", 0,   220102851, 220602850, "ENr331"}, 
  {"chr11", 0,   63940889,  64440888, "ENr332"}, 
  { 0, 0, 0, 0, 0}
}  ;
#endif
/*************************************************************************************/

static BOOL utcConfigRegions (UTC *utc)
{
  AC_HANDLE h = ac_new_handle () ;
  REGION *rg ;
  int ir = 0 ;
  Array regions = 0 ;
  BOOL ok = TRUE ;
  ACEIN ai = 0 ;
  const char *cp ;

  if (utc->chromConfigFile)
    ai = aceInCreateFromFile (utc->chromConfigFile, "r", 0, h) ;
  if (ai)
    {
      regions = arrayHandleCreate (128, REGION, utc->h) ;
      while (ok && aceInCard (ai))
	{
	  /* parse the file */
	  cp = aceInWord (ai) ;

	  /* jump empty lines and comment lines */
	  if (!cp || *cp == '#')
	    continue ;
	  /* break on ZERO ascii letter */
	  if (*cp == '0')
	    break ;
	  rg = arrayp (regions, ir++, REGION) ;
	  dictAdd (utc->chromDict, cp, &(rg->chrom)) ;
	  
	  if (aceInInt (ai, &(rg->a1)) && aceInInt (ai, &(rg->a2)))
	    {
	      if (rg->a1 < 3) rg->a1 = 3 ; 
	    }
	  else
	    { ok = FALSE ; break ; }
	  if ((cp = aceInWord (ai)))
	    {
	      dictAdd (utc->chromDict, cp, &(rg->regionName)) ;
	    }
	  else
	    { ok = FALSE ; break ; }
	}
    }

  if (ok && regions && arrayMax(regions))
    utc->chromArray = regions ;
  else
    arrayDestroy (regions) ;
  ac_free (h) ;
  return  utc->chromArray ? arrayMax(utc->chromArray) : 0 ;
} /* utcConfig */

/*************************************************************************************/

static BOOL utcFilter (UTC *utc, int chrom, int m1, int m2, int *z1p, int *z2p, int *snamp)
{
  REGION *rg ;
  int ir ;

  if (utc->chromArray && arrayMax(utc->chromArray))
    {
      for (ir = 0 ; ir < arrayMax (utc->chromArray) ; ir++)
	{
	  rg = arrayp (utc->chromArray, ir, REGION) ;
	  if (rg->chrom == chrom && rg->a1 - 1 <= m2 && rg->a2 >= m1)
	    { 
	      *z1p = rg->a1 - 1 ; *z2p = rg->a2 ;	 
	      *snamp = rg->regionName ; 
	      return TRUE ;
	    }
	}
      return FALSE ;
    }
  else
    {
      const char *cNam = dictName (utc->chromDict, chrom) ;
      dictAdd (utc->dict, cNam, snamp) ;
      *z1p = 0 ; *z2p = m1 < m2 ? m2 : m1 ;
    }
  return TRUE ;
} /* utcFilter */

/*************************************************************************************/

static BOOL utcConfigTracks (UTC *utc)
{
  AC_HANDLE h = ac_new_handle () ;
  LAB *lab ;
  int ir = 0, iLab ;
  Array labs = 0 ;
  BOOL ok = TRUE ;
  ACEIN ai = 0 ;
  const char *cp ;

  if (utc->trackConfigFile)
    ai = aceInCreateFromFile (utc->trackConfigFile, "r", 0, h) ;
  if (ai)
    {
      labs = arrayHandleCreate (128, LAB, utc->h) ;
      while (ok && aceInCard (ai))
	{
	  /* parse the file */
	  cp = aceInWord (ai) ;

	  /* jump empty lines and comment lines */
	  if (!cp || *cp == '#')
	    continue ;
	  /* break on ZERO ascii letter */
	  if (*cp == '0')
	    break ;
	  if (! strcmp (cp, "REF"))
	    {
	      if ((cp = aceInWord (ai)))
		dictAdd (utc->chromDict, cp, &REF) ;	
	      else
		{ ok = FALSE ; break ; }
	      continue ;
	    }
	  lab = arrayp (labs, ir++, LAB) ;
	  dictAdd (utc->chromDict, cp, &(lab->labId)) ;

	  if (! aceInInt (ai, (int*)(&lab->hasStop)) ||
	      (lab->hasStop != 0 && lab->hasStop != 1)
	      )
	    { 
	      fprintf (stderr, "ERROR: missing 1/0 hasstop flag in track %s in file %s\n"
		       , dictName (utc->chromDict, lab->labId)
		       , utc->trackConfigFile
		       );
	      ok = FALSE ; break ; 
	    }
	  if (! aceInInt (ai, (int*)(&lab->hasBinSection)) ||
	      (lab->hasBinSection != 0 && lab->hasBinSection != 1)
	      )
	    { 
	      fprintf (stderr, "ERROR: missing 1/0 hasBinSection flag in track %s in file %s\n"
		       , dictName (utc->chromDict, lab->labId)
		       , utc->trackConfigFile
		       );
	      ok = FALSE ; break ; 
	    }
	  if ((cp = aceInWord (ai)))
	    dictAdd (utc->chromDict, cp, &(lab->filName)) ;
	  else
	    { ok = FALSE ; break ; }
	}
    }

  if (ok && REF && labs && arrayMax(labs))
    for (iLab = 0, lab = arrp (labs, 0, LAB) ; iLab < arrayMax (labs) ; lab++, iLab++)
      if (lab->labId == REF)
	utc->refLab = lab ;
  if ((!REF || !utc->refLab) && !utc->aceExport)
    {
      fprintf (stderr, 
	       "ERROR: The reference track REF = %s defined in the file %s\n"
	       "is not present in the track list defined furtherdown in the same file\n"
	       , REF ? dictName (utc->chromDict, REF) : "\"this MANDATORY value is MISSING\""
 	       , utc->trackConfigFile ? utc->trackConfigFile : "NULL"
	       ) ;
      return 0 ;
    }
  
  /* try to locate the guigos track names */
  {
    GUIGO *guigo = 0 ;
    
    for (guigo = guigos ; *guigo->labNam ; guigo++)
      dictFind (utc->chromDict, guigo->labNam, &guigo->labId) ;
  }

  /* try to locate the AceView program */
  dictFind (utc->chromDict, "AceView", &ACVV) ;
  if (ok && (utc->aceExport || utc->refLab) && labs && arrayMax(labs))
    utc->labs = labs ;
  else
    arrayDestroy (labs) ;
  ac_free (h) ;
  return  utc->labs ? TRUE : FALSE ;
} /* utcConfigTracks */

/*************************************************************************************/
/*************************************************************************************/

static BOOL utcGetGtfData (UTC *utc, LAB *lab, BOOL utcNoCDS)
{
  BOOL ok = FALSE, isDown = TRUE ;
  char *cp ;
  int regionName ;
  int i, j, level, line = 0 ;
  int nn, m1, m2, c1, c2, x, z1, z2 ;
  int mrna, chrom, nmrna, nNewMrna = 0, nFiltered = 0, nClipped = 0, nExons ;
  MM *mm = 0 ;
  FILE *f = 0 ;

  /*
    if (utc->cds)
    f = filopen ("allMrnaCDS","txt","r") ;
    else
    f = filopen ("allMrnaExons","txt","r") ;
  */
  if (! filName (dictName (utc->chromDict, lab->filName),"txt","r") ||
      ! (f = filopen (dictName (utc->chromDict, lab->filName),"txt","r")))
    {
      fprintf (stderr, "ERROR: cannot open the track file %s : %s\n"
	       , dictName (utc->chromDict, lab->labId)
	       , dictName (utc->chromDict, lab->filName)
	       ) ;
      return FALSE ;
    }

  if (0) fprintf (stderr, "Parsing %s.txt\n", dictName (utc->chromDict, lab->filName)) ;
  if (! utc->dict)
    {
      nmrna = 0 ;
      utc->dict = dictHandleCreate (10000, utc->h) ;
      dictAdd (utc->dict, "__toto__", 0) ; /* avoid 0 */
      utc->aa = arrayHandleCreate (10000, MM, utc->h) ;
    }
  nmrna = arrayMax (utc->aa) ;
  level = freesetfile (f, 0) ;

  ok = TRUE ; nn = 0 ;

  while (freecard (level))
    {
      cp = freepos () ;
      while (*cp++) if (*cp == ',') *cp = ' ' ;
      ok = FALSE ;
      line++ ;
      if (lab->hasBinSection) 
	{
	  cp = freeword () ;
	  if (!cp) continue ;
	}   
      /* parse the mrna/chrom */
      cp = freeword () ;
      if (!cp) continue ;
      dictAdd (utc->dict, cp, &mrna) ;
      cp = freeword () ;
      if (!cp) continue ;
      dictAdd (utc->chromDict, cp, &chrom) ;
      cp = freeword () ;
      if (!cp) continue ;
      switch ((int) *cp)
	{
	case '+': isDown = TRUE ; break ;
	case '-': isDown = FALSE ; break ;
	default: continue ;
	}
      freeint (&m1) ; /* mrna start stop */
      freeint (&m2) ;
      c1 = c2 = 0 ;
      freeint (&c1) ; /* cds  start stop */
      freeint (&c2) ;
      if (utc->cds && !lab->hasStop && c1 != 0 && c2 != 0)
	{
	  if (isDown) c2 += 3 ; else c1 -= 3 ;
	} 
      freeint (&nExons) ;
      if (nExons < 1) continue ;
      if (utcNoCDS && (c1 || c2))
	continue ;
      if (! utcNoCDS && utc->cds)
	{ m1 = c1 ; m2 = c2 ; }
      if (!utcFilter (utc, chrom, m1, m2, &z1, &z2, &regionName))
	{ nFiltered++ ; continue ; }
      if (z1 > m1 || z2 < m2) nClipped++ ;
      nNewMrna++ ;
      /* store in new  mm */
	{
	  mm = arrayp (utc->aa, nmrna++, MM) ;
	  mm->chrom = chrom ;
	  mm->regionName = regionName ;
	  mm->isDown = isDown ;
	  mm->mrna = mrna ;
	  mm->m1 = m1 > z1 ? m1 : z1 ;
	  mm->m2 = m2 < z2 ? m2 : z2 ;
	  if (c1 || c2) mm->hasCds = TRUE ;
	  mm->c1 = !c1 || c1 > z1 ? c1 : z1 ;
	  mm->c2 = !c2 || c2 < z2 ? c2 : z2 ;
	  mm->z1 = z1 ;
	  mm->nExons = 0 ;
	  mm->labId = lab->labId ;
	}
      if (! utcNoCDS && utc->cds) /* we need that to reclip the exons */
	{
	  z1 = z1 < mm->c1 ? mm->c1 - 2 : z1 ;
	  z2 = z2 > mm->c2 ? mm->c2 + 2 : z2 ; 
	}
      else
	{
	  z1 = z1 < mm->m1 ? mm->m1 - 2 : z1 ;
	  z2 = z2 > mm->m2 ? mm->m2 + 2 : z2 ; 
	}
      mm->nExons = nExons ;
      for (i = 0 ; i < mm->nExons ; i++)
	{ 
	  freeint (&x) ; 
	  if (x < z1 + 1) x = z1 + 2 ; 
	  if (x > z2) x = z2 ; 
	  mm->a1[i] = x ;
	}
      if (utc->cds && !isDown && !lab->hasStop && mm->a1[0] == c1 + 3)
	mm->a1[0] -= 3;
      for (i = 0 ; i < mm->nExons ; i++)
	{ 
	  freeint (&x) ; 
	  if (x < z1) x = z1 ; 
	  if (x > z2 - 1) x = z2 - 2 ; 
	  mm->a2[i] = x ;
	}
      freeword() ;
      cp = freeword() ;
      mm->gene = 0 ;
      if (cp && *cp) dictAdd (utc->geneDict, cp, &mm->gene) ;
      cp = freeword () ;
      if (cp && !strcasecmp (cp, "cmpl"))
	{
	  if (isDown) 
	    mm->atgComplete = TRUE ;
	  else
	    mm->stopComplete = TRUE ;
	}
      cp = freeword () ;
      if (cp && !strcasecmp (cp, "cmpl"))
	{
	  if (! isDown) 
	    mm->atgComplete = TRUE ;
	  else
	    mm->stopComplete = TRUE ;
	}
      if (!utc->refLab || ! strcmp (dictName(utc->chromDict, utc->refLab->labId), "HAV"))
	mm->atgComplete = mm->stopComplete = TRUE ;
      if (utc->cds && isDown && !lab->hasStop && mm->a2[mm->nExons-1] == c2 - 3)
	mm->a2[mm->nExons-1] += 3;
      for (i = x = 0 ; i < mm->nExons ; i++)
	{
	  if (mm->a1[i] <= mm->a2[i])
	    {
	      if (x < i) 
		{ mm->a1[x] = mm->a1[i] ; mm->a2[x] = mm->a2[i] ;}
	      x++ ;
	    }
	}
      mm->nExons = x ; /* eliminate out of range exons */
      /* remember the original exons */
       for (i = 0 ; i < mm->nExons ; i++)
	{
	  mm->a1_original[i] = mm->a1[i] ;
	  mm->a2_original[i] = mm->a2[i] ;
	}
     
      if (!isDown)
	{
	  m1 = mm->m1 ; mm->m1 = mm->m2 ; mm->m2 = m1 ;
	  m1 = mm->c1 ; mm->c1 = mm->c2 ; mm->c2 = m1 ;

	  /* swap the exon order */
	  for (i = 0 ; i < mm->nExons/2 ; i++)
	    {
	      j = mm->nExons - i - 1 ;
	      m1 = mm->a1[i] ; mm->a1[i] = mm->a1[j] ; mm->a1[j] = m1 ;
	      m2 = mm->a2[i] ; mm->a2[i] = mm->a2[j] ; mm->a2[j] = m2 ;
	    } 
	  /* swap the exon coordinates */
	  for (i = 0 ; i < mm->nExons ; i++)
	    {
	      m1 = mm->a1[i] ; mm->a1[i] = mm->a2[i] ; mm->a2[i] = m1 ;
	    }
	}
      if (utc->introns)
	{
	  for (i = 0 ; i < mm->nExons - 1 ; i++)
	    { mm->a1[i] = mm->a2[i] ; mm->a2[i] = mm->a1[i + 1] ; }
	  mm->nExons-- ;
	}
      if (0)
	for (i = 0 ; i < mm->nExons ; i++)
	  fprintf (stderr, "exons %3d\t%d\t%d\n", i, mm->a1[i], mm->a2[i]) ;
    }
  freeclose (level) ; /* will close f */

  if (0) fprintf (stderr, "Found %d mrna; filtered %d mrna,; clipped %d mrna aaMax = %u nn = %d ok = %d\n"
		  , nNewMrna, nFiltered, nClipped, arrayMax (utc->aa), nn, ok) ;
  return nmrna ? TRUE : FALSE ;
} /* utcGetGtfData */

/*************************************************************************************/

static BOOL utcGetData (UTC *utc)
{
  LAB *lab, *lab0 ;
  int ii ;

  for (ii = 0, lab = arrp (utc->labs, 0, LAB) ; ii < arrayMax (utc->labs) ; lab++, ii++)
    utcGetGtfData (utc, lab, FALSE) ;
  if (utc->cds)
    {
      lab = labosNoCDS ;
      lab0 = utc->refLab ;
      *lab = *lab0 ;
      dictAdd (utc->chromDict, "HAV_no_CDS", &lab->labId) ;
      utcGetGtfData (utc, lab, TRUE) ;
    }
  return TRUE ;
} /* utcGetData */

/*************************************************************************************/
/*************************************************************************************/

static int utcBpCumulate (LAB *lab, MM *mm)
{
  int i, ii, nn ;
  BitSet bb ;

  if (mm->nExons > 1)
    bb = (mm->m1 < mm->m2 ? lab->bbDi : lab->bbRi) ;
  else
    bb = (mm->m1 < mm->m2 ? lab->bbDnoi : lab->bbRnoi) ; 

  for (ii = 0 ; ii < mm->nExons ; ii++)
    if (mm->m1 < mm->m2)
      {
	if (mm->a1[ii] < mm->z1)
	  messcrash ("mm->a1[%d] = %d < mm->z1 = %d", ii,  mm->a1[ii],  mm->z1) ;
	for (i = mm->a1[ii] ; i <= mm->a2[ii] ; i++)
	  { nn++; bitSet (bb, i - mm->z1) ; }
      }
    else
      {
	if (mm->a2[ii] < mm->z1)
	  messcrash ("mm->a2[%d] = %d < mm->z1 = %d", ii,  mm->a2[ii],  mm->z1) ;
	for (i = mm->a1[ii] ; i >= mm->a2[ii] ; i--)
	  { nn++; bitSet (bb, i - mm->z1) ; }
      }
      

  return nn ;
} /* utcBpCumulate */

/*************************************************************************************/

static int utcBpCount (UTC *utc, BOOL isDown)
{
  int ii, nn = 0, xx ;
  BitSet b, b0, b1i, b1noi ; 
  LAB *lab ;

  /* havana */
  lab = utc->refLab ;
  /* in havana, merge all bases */
  bitSetOR (lab->bbDi, lab->bbDnoi) ;
  bitSetOR (lab->bbRi, lab->bbRnoi) ;
  bitSetReCreate (lab->bbDnoi, 100) ;   bitSetReCreate (lab->bbRnoi, 100) ;

  b0 = isDown ? lab->bbDi : lab->bbRi ;
  for (ii = 0, lab = arrp (utc->labs, 0, LAB) ; ii < arrayMax (utc->labs) ; lab++, ii++)
    {
      b1i = isDown ? lab->bbDi : lab->bbRi ;
      b1noi = isDown ? lab->bbDnoi : lab->bbRnoi ;
      b = bitSetCreate (1000000, 0) ;
      bitSetOR (b, b1i) ;
      bitSetOR (b, b1noi) ;
      lab->nbp += bitSetCount (b) ;
      xx = bitSetCount (b) ;
      bitSetDestroy (b) ;  

      b = bitSetCreate (1000000, 0) ;
      bitSetOR (b, b1i) ;
      bitSetOR (b, b1noi) ;
      bitSetAND (b, b0) ;
      lab->nvp += bitSetCount (b) ;
      xx -= bitSetCount (b) ;
      bitSetDestroy (b) ;

      b = bitSetCreate (1000000, 0) ;
      bitSetOR (b, b0) ;
      bitSetMINUS (b, b1i) ;
      bitSetMINUS (b, b1noi) ;
      lab->nfn += bitSetCount (b) ;
      bitSetDestroy (b) ;
      
      b = bitSetCreate (1000000, 0) ;
      bitSetOR (b, b1i) ;
      bitSetMINUS (b, b0) ;
      lab->nfpi += bitSetCount (b) ;
      xx -= bitSetCount (b) ;
      bitSetDestroy (b) ;

      b = bitSetCreate (1000000, 0) ;
      bitSetOR (b, b1noi) ;
      bitSetMINUS (b, b0) ;
      bitSetMINUS (b, b1i) ;
      lab->nfpnoi += bitSetCount (b) ;
      xx -= bitSetCount (b) ;
      bitSetDestroy (b) ;
      
      if (xx) messcrash ("bad non zero xx = %d", xx) ;
    }
     
  return nn ;
} /* utcBpCount */

/*************************************************************************************/
/* count false positive/negative in bp */
static BOOL utcBpChrom (UTC *utc, int regionName)
{
  int ii, iLab, nn ;
  MM* mm ;
  LAB *lab ;

  /* loop on all havana */
  nn = arrayMax (utc->aa) ;
  for (ii = 0, mm = arrp (utc->aa, 0, MM) ; ii < nn ; ii++, mm++)
    if (mm->regionName == regionName)
      for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
	if (lab->labId == mm->labId)
	  { utcBpCumulate (lab, mm) ; break ; }
  
  return TRUE ;
} /* utcBpChrom */

/*************************************************************************************/

static int truePositivesLabsBpOrder (const void *a, const void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = la->nvp - lb->nvp ;
  if (n) return -n ;              /* best true-positives = best sensibility */
  n = la->nfpi - lb->nfpi ;
  if (n) return n ;               /* low false negatives */
  n = la->labId - lb->labId ;
  return n ;
} /* truePositivesLabsBpOrder */

/*************************************************************************************/
/* count false positive/negative in bp */
static BOOL utcBpRegion (UTC *utc, int regionName)
{
  LAB *lab ;
  int iLab ;
  AC_HANDLE h = ac_new_handle () ;

  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      lab->bbDi= bitSetCreate (1000000, h) ;
      lab->bbRi = bitSetCreate (1000000, h) ;
      lab->bbDnoi = bitSetCreate (1000000, h) ;
      lab->bbRnoi = bitSetCreate (1000000, h) ;
    }
  lab->bbDi= bitSetCreate (1000000, h) ;
  lab->bbRi = bitSetCreate (1000000, h) ;
  lab->bbDnoi = bitSetCreate (1000000, h) ;
  lab->bbRnoi = bitSetCreate (1000000, h) ;
  
  utcBpChrom (utc, regionName) ;
  utcBpCount (utc, TRUE) ;
  utcBpCount (utc, FALSE) ;
  ac_free (h) ;

  return TRUE ;
} /* utcBpRegion */

/*************************************************************************************/
/* count false positive/negative in bp */
static BOOL utcBp (UTC *utc)
{
  int nn, ir, iLab ;
  LAB *lab0, *lab ;
  Array  sortedLabs = 0 ;
  REGION *bb ;

   for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    lab->nbp = lab->nvp = lab->nfpi = lab->nfpnoi = lab->nfn = 0 ;
  lab->nbp = lab->nvp = lab->nfpi = lab->nfpnoi = lab->nfn = 0 ;

  /* loop on all havana */
  for (ir = 0 ; ir < arrayMax (utc->chromArray) ; ir++)
    {
      bb = arrayp (utc->chromArray, ir, REGION) ;
      utcBpRegion (utc, bb->regionName) ; 
    }

  nn = 0 ; sortedLabs = arrayHandleCreate (100, LAB, utc->h) ; 
  lab0 = utc->refLab ;
  lab0->nvp = lab0->nbp ; 
  lab0->nfpnoi = lab0->nfpi = 0 ;
  for (iLab = 0, lab0 = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab0++, iLab++)
    {
      lab = arrayp (sortedLabs, nn++, LAB) ;
      *lab = *lab0 ;
    }
  if (! utc->nosort)
    arraySort (sortedLabs, truePositivesLabsBpOrder) ;
  freeOut ("Program"
	   "\tTrue positive"
	   "\tFalse positive in transcript with intron"
	   "\tFalse positive in intronless transcript"
	   "\tFalse negative"
	   "\tProgram"
	   "\tSensitivity"
	   ) ;
  if (! utc->multi)
    freeOut ("\tSensitivity eval" 
	     "\tProgram"
	     "\tRaw specificity"
	     "\tSpecificity eval"
	     ) ;
  else
    freeOut ("\tRaw specificity"
	     ) ;
  
  freeOut ("\tSpecificity, ignoring single exon transcripts"
	   "\tNumber bp"
	   "\tVerif\n"
	   ) ;

  for (lab = lab0 = arrp (sortedLabs, 0, LAB), nn = 0 ; nn < arrayMax (sortedLabs) ; nn++, lab++)
    {
      GUIGO *guigo = 0 ;

      if (! utc->multi)
	for (guigo = guigos ; guigo->labId ; guigo++)
	  if (guigo->labId == lab->labId)
	    break ;
      if (guigo && guigo->labId != lab->labId)
	guigo = 0 ;

      if (!lab->nbp)
	{
	  freeOutf ("%s nbp == 0\n", dictName (utc->chromDict, lab->labId)) ;
	  continue ;
	}

      freeOutf ("%s\t%d\t%d\t%d\t%d\t%s\t%.2f"
		, dictName (utc->chromDict, lab->labId)
		/* #bp = vrai-positif + faux-positifs ; faux-nagatifs */
		, lab->nvp, lab->nfpi, lab->nfpnoi, lab->nfn 
		, dictName (utc->chromDict, lab->labId)
		, 100.0 * lab->nvp/lab0->nbp
		) ;

      if (! utc->multi)
	freeOutf ("\t%.2f\t%s"
		  , guigo ? (utc->cds ? guigo->cdsBpSn : guigo->bpSn) : 0
		  , dictName (utc->chromDict, lab->labId)
		  ) ;

      freeOutf ("\t%.2f"
		, 100.0 * lab->nvp/(lab->nvp + lab->nfpnoi + lab->nfpi)
	      ) ;
      if (! utc->multi)
	freeOutf ("\t%.2f", guigo ? (utc->cds ? guigo->cdsBpSp : guigo->bpSp) : 0) ;
      freeOutf ("\t%.2f"
		, 100.0 * lab->nvp/(lab->nvp + lab->nfpi)
	      ) ;
 
      freeOutf ("\t%d\t%d\n"
		, lab->nbp
		, lab->nbp - (lab->nvp + lab->nfpnoi + lab->nfpi)
	      ) ;
    }
  return TRUE ;
} /* utcBp */

/*************************************************************************************/
/*************************************************************************************/
/* count false positive/negative in bp */
typedef struct exonStruct { BOOL isDown, isFp ; int a1, a2 ; int regionName ; int isHavana, isAceView, count, nExons, labId ; } EX ;
static int exonOrder (const void *va, const void *vb)
{
  const EX *ea = (const EX*)va, *eb = (const EX*)vb ;
  int nn ;

  nn = ea->regionName - eb->regionName ;
  if (nn) return nn ;

  nn = ea->isDown - eb->isDown ;
  if (nn) return nn ;
 
  nn = ea->a1 - eb->a1 ;
  if (nn) return nn ;

  nn = ea->a2 - eb->a2 ;
  if (nn) return nn ;
  return 0 ;
} /* exonOrder */

/*************************************************************************************/

static int exonOrderLab (const void *va, const void *vb)
{
  const EX *ea = (const EX*)va, *eb = (const EX*)vb ;
  int nn ;

  nn = ea->regionName - eb->regionName ;
  if (nn) return nn ;

  nn = ea->isDown - eb->isDown ;
  if (nn) return nn ;
 
  nn = ea->a1 - eb->a1 ;
  if (nn) return nn ;

  nn = ea->a2 - eb->a2 ;
  if (nn) return nn ;

  nn = ea->labId - eb->labId ;
  if (nn) return nn ;

  return 0 ;
} /* exonOrderLab */

/*************************************************************************************/

static int utcExonCountLab (LAB *lab0, LAB *lab)
{
  Array exons0 = lab0->exons ;
  Array exons = lab->exons ;
  int ii0, nn0 = arrayMax (exons0), oldA1 = -9999, oldA2 = 0;
  int ii, nn = arrayMax (exons) ;
  int k ;
  EX *ee0, *ee ;

  lab0->nex = lab0->extp = nn0 ;
  lab->nex = nn ;
  if (!nn && !nn0)
    return 0 ;

  ii = ii0 = 0 ; ee = ee0 = 0 ;
  if (nn0)
    ee0 = arrp (exons0, 0, EX) ;
  if (nn)
    ee = arrp (exons, 0, EX) ;
  while (ii0 < nn0 && ii < nn)
    {
      k = exonOrder (ee0, ee) ;
      if (k < 0)
	{ 
	  oldA1 = -9999 ;
	  ee0++ ; ii0++ ; lab->exfn++ ;
	}
      else if (k > 0)
	{ 
	  if (ee->a1 == oldA1 && ee->a2 == oldA2)
	    lab->exfpKnown++ ; 
	  else
	    {
	      if (ee->nExons > 1)
		lab->exfpNew++ ; 
	      else
		lab->exfpNew1++ ; 
	      oldA1 = -9999 ;
	    }
	  ee->isFp = TRUE ;
	  ee++ ; ii++ ; 
	}
      else
	{ 
	  lab->extp++ ; 
	  oldA1 = ee->a1 ; oldA2 = ee->a2 ;
	  ee0++ ; ii0++ ; ee++ ; ii++ ; 
	}
    }  
  lab->exfn += nn0 - ii0 ;

  while (ii < nn)
    {
      if (ee->nExons > 1)
	lab->exfpNew++ ; 
      else
	lab->exfpNew1++ ; 
      ee->isFp = TRUE ;
      ee++ ; ii++ ; 
    }

  return nn ;
} /* utcExonCountLab */

/*************************************************************************************/
/* split the false positive exons in cds case according to
 * if they appear in havana devoid of CDS annot
 */

static int utcExonSplitFp (LAB *lab0, LAB *lab)
{
  Array exons0 = lab0->exons ;
  Array exons = lab->exons ;
  int ii0, nn0 = arrayMax (exons0) ;
  int ii, nn = arrayMax (exons) ;
  int k ;
  EX *ee0, *ee ;

  lab->exfpCds = 0 ;
  if (!nn && !nn0)
    return 0 ;
  if (!lab->exfpNew1 && !lab->exfpNew1)
     return 0 ;
  ii = ii0 = 0 ; ee = ee0 = 0 ;
  if (nn0)
    ee0 = arrp (exons0, 0, EX) ;
  if (nn)
    ee = arrp (exons, 0, EX) ;
  while (ii0 < nn0 && ii < nn)
    {
      k = exonOrder (ee0, ee) ;
      if (k < 0)
	{ 
	  ee0++ ; ii0++ ; 
	}
      else if (k > 0)
	{
	  ee++ ; ii++ ; 
	}
      else
	{ 
	  if (ee->isFp)
	    {
	      if (ee->nExons > 1)
		lab->exfpNew-- ; 
	      else
		lab->exfpNew1-- ; 
	      lab->exfpCds++ ;
	    }
	  ee0++ ; ii0++ ; ee++ ; ii++ ; 
	}
    }  
  lab->exfn += nn0 - ii0 ;

  return lab->exfpCds ;
} /* utcExonSplitFp */

/*************************************************************************************/

static int utcExonCumulate (UTC *utc, LAB *lab, MM *mm, Array allExons)
{
  int ii ;
  Array exons = lab->exons ;
  int nn = arrayMax (exons) ; 
  int nnAll = allExons ? arrayMax (allExons) : 0 ;
  EX *ee, *ff ;

  for (ii = 0 ; ii < mm->nExons ; ii++)
    {
      ee = arrayp (exons, nn++, EX) ;
      ff = allExons ? arrayp (allExons, nnAll++, EX) : 0 ;
      if (mm->m1 < mm->m2)
	{
	  ee->isDown = TRUE ; 
	  ee->regionName = mm->regionName ;
	  ee->a1 =  mm->a1[ii] ; ee->a2 =  mm->a2[ii] ;
	}
      else
	{
	  ee->isDown = FALSE ;
	  ee->regionName = mm->regionName ;
	  ee->a1 =  mm->a2[ii] ; ee->a2 =  mm->a1[ii] ;
	}
      ee->nExons = utc->introns ? 2 : (mm->nExons > 1 ? 2 : 1) ;
      if (ff)
	{
	  *ff = *ee ;
	  if (lab->labId == REF) ff->isHavana = TRUE ;
	  if (lab->labId == ACVV) ff->isAceView = TRUE ;
	  ff->labId = lab->labId ;
	}
    }
  return nn ;
} /* utcExonCumulate */

/*************************************************************************************/

static int utcExonCount (UTC *utc)
{
  LAB *lab, *lab0 ;
  int ii, iLab, nn ;
  MM *mm1 ;

  lab0 = utc->refLab ;
  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    if (lab != lab0)
      utcExonCountLab (lab0, lab) ;

  lab0 = labosNoCDS ;
  if (utc->exons && utc->cds)
    {
      nn = arrayMax (utc->aa) ;
      lab0->exons = arrayCreate (10000, EX) ;
      for (ii = 0, mm1 = arrp (utc->aa, 0, MM) ; ii < nn ; ii++, mm1++)
	if (mm1->labId == lab0->labId)
	  utcExonCumulate (utc, lab0, mm1, 0) ; 
      arraySort (lab0->exons, exonOrder) ;
      if (!utc->multi)
	arrayCompress (lab0->exons) ;
      
      for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
	utcExonSplitFp (lab0, lab) ;
    }

  return 1 ;
} /* utcExonCount */

/*************************************************************************************/

static int truePositivesLabsOrder (const void *a, const  void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = la->extp - lb->extp ;
  if (n) return -n ;              /* best true-positives = best sensibility */
  n = la->exfn - lb->exfn ;
  if (n) return n ;               /* low false negatives */
  n = la->exfpNew - lb->exfpNew + la->exfpNew1 - lb->exfpNew1 + la->exfpKnown - lb->exfpKnown ;
  if (n) return n ;               /* low false positives */
  n = la->labId - lb->labId ;
  return n ;
}

/*************************************************************************************/

static BOOL utcExons (UTC *utc)
{
  AC_HANDLE h = 0 ;
  int ii, iLab, nn ;
  LAB *lab0, *lab ;
  MM* mm1 ;
  Array sortedLabs = 0, allExons = 0 ;

  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      lab->exons = arrayCreate (10000, EX) ;
      lab->nex = lab->extp = lab->exfpNew = lab->exfpNew1 = lab->exfpCds = lab->exfpKnown = lab->exfn = 0 ;
    }
  lab->exons = arrayCreate (10000, EX) ;
  lab->nex = lab->extp = lab->exfpNew = lab->exfpNew1 = lab->exfpCds = lab->exfpKnown = lab->exfn = 0 ;

  /* loop on all havana */
  allExons = arrayHandleCreate (10000, EX, h) ;
  nn = arrayMax (utc->aa) ;
  for (ii = 0, mm1 = arrp (utc->aa, 0, MM) ; ii < nn ; ii++, mm1++)
    for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
      if (lab->labId == mm1->labId)
	{  utcExonCumulate (utc, lab, mm1, allExons) ; break ; }

  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      arraySort (lab->exons, exonOrder) ;
      if (!utc->multi) arrayCompress (lab->exons) ;
    }
  utcExonCount (utc) ;
  if (utc->introns)
    {
      freeOutf ("Program"
	       "\tTrue positive, identical to %s"
		, dictName (utc->chromDict, utc->refLab->labId)
	       ) ;
      if (utc->multi)
	freeOutf ("\tFalse positive, present in %s but over-used"
		  , dictName (utc->chromDict, utc->refLab->labId)
		  ) ;
      freeOutf ("\tFalse positive, new"
	       "\tFalse negative"
	       "\tProgram"
	       "\tSensitivity"
	       "\tProgram"
	       "\tSpecificity"
	       "\t%% missed introns"
	       "\t%% wrong introns"
	       "\tAll introns"
		, dictName (utc->chromDict, utc->refLab->labId)
	       ) ;
    }
  else
    {
      freeOutf ("Program"
	       "\tTrue positive, identical to %s"
		, dictName (utc->chromDict, utc->refLab->labId)
	       ) ;
      if (utc->multi)
	freeOutf ("\tFalse positive, present in %s but over-used" 
		  , dictName (utc->chromDict, utc->refLab->labId)
		 ) ;
      if (utc->cds)
	freeOutf ( "\tFalse positive, non-coding in %s"
		  , dictName (utc->chromDict, utc->refLab->labId)
		  ) ;
      freeOutf ("\tFalse positive, new in multi-exon transcript"
		"\tFalse positive, new in single exon transcript"
		"\tFalse negative"
		) ;
      if (! utc->multi)
	freeOut ("\tProgram"
		 "\tSensitivity"
		 "\tSensitivity eval"
		 "\tProgram"
		 "\tRaw specificity"
		 "\tSpecificity eval"
		 ) ;
      else
	freeOut ("\tProgram"
		 "\tSensitivity"
		 "\tProgram"
		 "\tRaw specificity"
		 ) ;
      
      freeOut ("\tSpecificity, ignoring single exon transcripts" 
	       "\t% missed exons"
	       "\t% wrong exons"
	       "\tAll exons"
	       ) ;
    }
  if (!utc->multi)
    freeOut ("\tVerif2") ;
  freeOut ("\n") ;

  nn = 0 ; sortedLabs = arrayHandleCreate (100, LAB, utc->h) ; 
  for (iLab = 0, lab0 = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab0++, iLab++)
    {
      lab = arrayp (sortedLabs, nn++, LAB) ;
      *lab = *lab0 ;
    }
  if (! utc->nosort)
    arraySort (sortedLabs, truePositivesLabsOrder) ;
  for (lab = lab0 = arrp (sortedLabs, 0, LAB), nn = 0 ; nn < arrayMax (sortedLabs) ; nn++, lab++)
    {
      GUIGO *guigo = 0 ;

      if (! utc->multi && utc->exons)
	for (guigo = guigos ; guigo->labId ; guigo++)
	  if (guigo->labId == lab->labId)
	    break ;
      if (guigo && guigo->labId != lab->labId)
	guigo = 0 ;

      if (!lab->nex)
	{
	  freeOutf ("%s n%s == 0\n", dictName (utc->chromDict, lab->labId), utc->introns ? "introns" : "exons" ) ;
	  continue ;
	}

      freeOutf ("%s\t%d", dictName (utc->chromDict, lab->labId), lab->extp) ;
      if (utc->multi)
	freeOutf ("\t%d",lab->exfpKnown) ;
      if (utc->cds && utc->exons)
	freeOutf ("\t%d", lab->exfpCds) ;
      freeOutf ("\t%d", lab->exfpNew) ;
      if (utc->exons)
	freeOutf ("\t%d", lab->exfpNew1) ;
      freeOutf ("\t%d", lab->exfn) ;
      /* SN / SP */
      if (utc->introns)
	freeOutf ("\t%s\t%.2f\t%s\t%.2f"
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab0->nex 
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab->nex
		  ) ;
      else if (! utc->multi && utc->exons && utc->cds)
	freeOutf ("\t%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f\t%.2f"
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab0->nex 
		  , guigo ? guigo->cdsExonSn : 0
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab->nex
		  , guigo ? guigo->cdsExonSp : 0
		  , 100.0 * lab->extp/(lab->nex - lab->exfpNew1)
		  ) ;
      else if (! utc->multi && utc->exons && ! utc->cds)
	freeOutf ("\t%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f\t%.2f"
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab0->nex 
		  , guigo ? guigo->exonSn : 0
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab->nex
		  , guigo ? guigo->exonSp : 0
		  , 100.0 * lab->extp/(lab->nex - lab->exfpNew1)
		  ) ;
      else
	freeOutf ("\t%s\t%.2f\t%s\t%.2f\t%.2f"
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab0->nex 
		  , dictName (utc->chromDict, lab->labId)
		  , 100.0 * lab->extp/lab->nex
		  , 100.0 * lab->extp/(lab->nex - lab->exfpNew1)
		  ) ;
      /*  missed introns wrong introns */
      freeOutf ("\t%.2f\t%.2f"
		, 100.0 * lab->exfn/lab0->nex
		, 100.0 * (lab->exfpKnown + lab->exfpNew + lab->exfpNew1 + lab->exfpCds) /lab->nex
		) ;
      /* all introns */  
      if (!utc->multi)
	freeOutf ("\t%d\t%d\n"
		  , lab->nex - (utc->multi ? 0 : lab->exfpKnown)
		  , lab->exfpKnown
		  ) ;
      else
	freeOutf ("\t%d\n", lab->nex) ;
    }

  /* histogramme des supports d'exons */
  if (arrayMax (allExons))
    {
      EX *ee, *ff ;
      int ii, jj, countMax = 0, nJustHavana ;
      KEYSET acvHisto, otherHisto ;

      arraySort (allExons, exonOrderLab) ;
      for (ii = jj = 0, ee = arrp (allExons, 0, EX), ff = ee ;
	   ii < arrayMax (allExons) ; ii++, ee++)
	{
	  if (! exonOrder (ee, ff)) /* repetition of same coordinates */
	    {
	      if (ee == ff || ee->labId != ff->labId)
		ff->count++ ;
	      ff->labId = ee->labId ;
	      ff->isHavana |= ee->isHavana ;
	      ff->isAceView |= ee->isAceView ;
	      if (ff->count > countMax)
		countMax = ff->count ;
	      if (0 && (! ff->isAceView) && ff->count >=14)
		printf ("Missed n=%d in aceview %s %d %d\n"
			, ff->count
			, dictName(utc->dict, ee->regionName)
			, ee->a1, ee->a2
			) ;
	    }
	  else
	    {
	      ff++ ; jj++ ;
	      if (jj < ii) 
		*ff = *ee ; 
	      ff->count = 1 ;
	    }
	}
      arrayMax (allExons) = jj + 1 ;
      acvHisto = keySetHandleCreate (h) ; keySet (acvHisto, countMax +1) = 0 ;
      otherHisto = keySetHandleCreate (h) ; keySet (otherHisto, countMax +1) = 0 ;
      countMax = 0 ;
      for (ii = nJustHavana = 0, ee = arrp (allExons, 0, EX) ; ii < arrayMax (allExons) ; ii++, ee++)
	{
	  if (ee->isHavana)
	    { if (ee->count == 1) nJustHavana++ ; }
	  else
	    {
	      if (ee->isAceView)
		keySet (acvHisto, ee->count) += 1 ; 
	      else
		keySet (otherHisto, ee->count) += 1 ; 
	      if (ee->count > countMax)
		countMax = ee->count ;
	    }	  
	}
      freeOut ("\n\nSupport histogram: how many features are supported by n methods\n") ;
      freeOutf( "%d are specific of %s\n", nJustHavana, dictName (utc->chromDict, utc->refLab->labId)) ;
      freeOut ("Number of supporting methods"
	       "\tPresent in AceView"
	       "\tAbsent in AceView"
	       "\n"
	       ) ;
      for (ii = 0 ; ii <= countMax ; ii++)
	freeOutf("%d\t%d\t%d\n"
		 , ii, keySet (acvHisto, ii), keySet (otherHisto,ii)) ;
    }
      
  ac_free (h) ;
  return TRUE ;
} /* utcExons */

/*************************************************************************************/
/*************************************************************************************/
/* export number of transcript per region per method */
static int  utcRegionsDo (UTC *utc, LABID labId, REGION *rg)
{
  int k = 0, ii, nn ;
  MM *mm ;

  nn = arrayMax (utc->aa) ;
  for (ii = 0, mm = arrp (utc->aa, 0, MM) ; ii < nn ; ii++, mm++)
    if (mm->regionName == rg->regionName && mm->labId == labId )
      { 
	k++ ; 
	if (0) printf ("ZZZZ\t%s\n", dictName(utc->dict, mm->mrna)) ;
      }

  return k ;
} /* utcRegionsDo */
      
/*************************************************************************************/
/* export number of transcript per region per method */
static BOOL  utcRegions (UTC *utc)
{
  REGION *rg ;
  LAB *lab ;
  int ir, iLab, nn ;

  freeOutf ("Program") ;
  for (ir = 0 ; ir < arrayMax (utc->chromArray) ; ir++)
    {
      rg = arrayp (utc->chromArray, ir, REGION) ;
      freeOutf ("\t%s", dictName (utc->chromDict,rg->regionName)) ;
    }
  /* loop on all methods */
  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      lab->nbp = lab->nvp = 0 ;
      freeOutf ("\n%s", dictName (utc->chromDict, lab->labId)) ;
      /* loop on all chroms */
      for (ir = 0 ; ir < arrayMax (utc->chromArray) ; ir++)
	{
	  rg = arrayp (utc->chromArray, ir, REGION) ;
	  nn =  utcRegionsDo (utc, lab->labId, rg) ;
	  freeOutf ("\t%d", nn) ;
	}
    }
  freeOut ("\n\n") ;

  return TRUE ;
} /* utcRegions */
      
/*************************************************************************************/
/*************************************************************************************/
/* export global statistics on all methods */
static BOOL  utcGlobal (UTC *utc)
{
  MM *mm ;
  EX *exp ;
  LAB *lab ;
  LABID labId ;
  int i, ii, iLab, nbp, ncds, nTr, nTr1, nIntrons, nExons, nExonsTerminal ;
  int nGenes, nGenesWithIntrons, nCodingGenesWithIntrons  ;
  Array allExons = 0, allIntrons = 0 ;
  KEYSET genes = 0, genesWithIntrons = 0, codingGenesWithIntrons = 0 ;
  
  freeOut ("Program"
	   "\tNucleotides"
	   "\tExons"
	   "\t%% Terminal exons"
	   "\tIntrons"
	   "\tmRNA with intron"
	   "\tSingle exon mRNA"
	   "\t% all mRNA with CDS"
	   "\tGenes"
	   "\tGenes with introns"
	   "\tCoding Genes with introns"
	   ) ;
  
  /* loop on all methods */
  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      labId = lab->labId ;
      
      nbp = ncds = nTr = nTr1 = nIntrons = nExons = nExonsTerminal = 0 ;
      allExons = arrayReCreate (allExons, 10000, EX) ;
      allIntrons = arrayReCreate (allIntrons, 10000, EX) ;
      freeOutf ("\n%s", dictName (utc->chromDict, lab->labId)) ;
      genes = keySetReCreate (genes) ;   
      genesWithIntrons = keySetReCreate (genesWithIntrons) ;
      codingGenesWithIntrons = keySetReCreate (codingGenesWithIntrons) ;
      nGenes = nGenesWithIntrons = nCodingGenesWithIntrons = 0 ;
      for (ii = 0, mm = arrp (utc->aa, ii, MM) ; ii < arrayMax (utc->aa); mm++, ii++)
	{
	  if (mm->labId != labId)
	    continue ;
	  if (mm->nExons < 1)
	    continue ;
	  else if (mm->nExons == 1)
	    nTr1++ ;
	  else 
	    {
	      keySet (genesWithIntrons, nGenesWithIntrons++) = mm->gene ;
	      if (mm->hasCds)
		keySet (codingGenesWithIntrons, nCodingGenesWithIntrons++) = mm->gene ;
	    }
	  nTr++ ;
	  keySet (genes, nGenes++) = mm->gene ;

	  nExonsTerminal += (mm->nExons == 1 ? 1 : 2) ;
	  if (mm->hasCds) ncds++ ;
	  for (i = 0 ; i < mm->nExons ; i++)
	    {
	      exp = arrayp (allExons, nExons++, EX) ;
	      exp->regionName = mm->regionName ;
	      exp->a1 = mm->a1[i] ;
	      exp->a2 = mm->a2[i] ;
	    }
	  for (i = 0 ; i < mm->nExons - 1 ; i++)
	    {
	      exp = arrayp (allIntrons, nIntrons++, EX) ;
	      exp->regionName = mm->regionName ;
	      exp->a1 = mm->a2[i] ;
	      exp->a2 = mm->a1[i+1] ;
	    }
	}
      arraySort (allExons, exonOrder) ; 
      arraySort (allIntrons, exonOrder) ;
      keySetSort (genes) ;
      keySetCompress (genes) ;
      keySetSort (genesWithIntrons) ;
      keySetCompress (genesWithIntrons) ;
      keySetSort (codingGenesWithIntrons) ;
      keySetCompress (codingGenesWithIntrons) ;
      if (!utc->multi)
	{
	  arrayCompress (allExons) ;
	  arrayCompress (allIntrons) ;
	  nExonsTerminal = 0 ;
	}
      nExons = arrayMax (allExons) ;
      nIntrons = arrayMax (allIntrons) ;
       nbp = 0 ;
       if (utc->multi) 
	for (i = 0, exp = arrp (allExons, 0, EX) ; i < arrayMax (allExons) ; i++, exp++)
	  nbp += (exp->a2 > exp->a1 ? exp->a2 - exp->a1 + 1 : exp->a1 - exp->a2 + 1) ;
      freeOutf ("\t%d\t%d\t%.0f\t%d\t%d\t%d\t%.0f\t%d\t%d\t%d"
		, nbp
		, nExons
		, 100.0 * nExonsTerminal/(nExons ? nExons : 1)
		, nIntrons
		, nTr - nTr1
		, nTr1
		, 100.0 * ncds/(nTr ? nTr : 1)
		, keySetMax (genes) 
		, keySetMax (genesWithIntrons) 
		, keySetMax (codingGenesWithIntrons) 
		) ;
    }
  freeOut ("\n\n") ;
  arrayDestroy (allExons) ;
  arrayDestroy (allIntrons) ;
  keySetDestroy (genes) ;  
  keySetDestroy (genesWithIntrons) ;
  keySetDestroy (codingGenesWithIntrons) ;
  return TRUE ;
} /* utcGlobal */
      
/*************************************************************************************/
/*************************************************************************************/
/* export in .ace format the genes of all methods mapping in the acceptable regions */
static BOOL  utcAceExport (UTC *utc)
{
  char buf[1000] ;
  char *cdsBuff = utc->cds ? "_cds" : "" ;
  MM *mm ;
  LAB *lab ;
  int i, ii, iLab, a1, a2, x1, x2, c1, c2, u1, u2, y1, y2 ;
  /*
    LABID ll2[] = { HAV, ACV, ENS, EXO, EXH, FGH, KGP, ECP, EWP, ZERO} ;
    LABID ll[] = { HAV, ACV, AVP, ZERO} ;
  */
  if (utc->aa)
    for (ii = 0, mm = arrp (utc->aa, ii, MM) ; ii < arrayMax (utc->aa); mm++, ii++)
      {
	/*
	  for (i = 0 ; ll[i] ; i++)
	  if (mm->labId == ll[i]) 
	  break ;
	  if (!ll[i])
	  continue ;
	*/
	for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
	  if (lab->labId == mm->labId)
	    {
	      sprintf (buf, "\"%s%s_%s"
		       , dictName (utc->chromDict, lab->labId)
		       , cdsBuff,  1 + freeprotect(dictName (utc->dict, mm->mrna))) ;
	      break ;
	    }
	freeOutf ("Sequence %s\nmRNAs %s %d %d\n\n"
		  , dictName (utc->dict, mm->regionName)
		  , buf
		  , mm->m1 - mm->z1 + 1 
		  , mm->m2 - mm->z1 + 1 
		  ) ;
	freeOutf ("mRNA %s\n", buf) ;
	if (mm->gene)
	  freeOutf ("Model_of_gene  \"%s\"\n", dictName (utc->geneDict, mm->gene)) ;
	freeOutf ("Covers %d \"bp from\" %d to %d\n"
		  , mm->m1 < mm->m2 ? mm->m2 - mm->m1 + 1 :  mm->m1 - mm->m2 + 1
		  , mm->m1 - mm->z1 + 1 
		  , mm->m2 - mm->z1 + 1 
		  ) ;
	c1 = mm->m1 < mm->m2 ? mm->c1 - mm->m1 + 1 : mm->m1 - mm->c1 + 1 ;
	c2 = mm->m1 < mm->m2 ? mm->c2 - mm->m1 + 1 : mm->m1 - mm->c2 + 1 ;
	c1 = 1 ; c2 =  mm->m1 < mm->m2 ? mm->m2 - mm->m1 + 1 :  mm->m1 - mm->m2 + 1 ;
	if (c1 < c2)
	  freeOutf ("Product %d %d\n", c1, c2) ;
	freeOutf ("IntMap %s %d %d\n"
		  , dictName (utc->chromDict, mm->chrom)
		  , mm->m1+ 1, mm->m2
		  ) ;
	freeOutf ("Nb_exons %d\n",  mm->nExons) ;
	for (x1 = 1, i = 0 ; i < mm->nExons ; i++)	
	  { 
	    a1 = mm->m1 < mm->m2  ? mm->a1[i] - mm->m1 + 1 : mm->m1 - mm->a1[i] + 1 ;
	    a2 = mm->m1 < mm->m2  ? mm->a2[i] - mm->m1 + 1 : mm->m1 - mm->a2[i] + 1 ;
	    x2 = x1 + a2 - a1 ;
	    freeOutf ("Splicing %d %d %d %d Exon %d Length %d bp\n"
		      , a1, a2, x1, x2, i, x2 - x1 + 1
		      ) ;
	    if (c1 < c2)
	      {
		u1 = a1 > c1 ? a1 : c1 ;
		u2 = a2 < c2 ? a2 : c2 ; 
		y1 = x1 + u1 - a1 ;
		y2 = x2 + u2 - a2 ;
		if (u1 <= u2)
		  freeOutf ("Coding %d %d %d %d A\n", u1, u2, y1, y2) ;
	      }
	    x1 = x2 + 1 ;
	    if (i == mm->nExons - 1)
	      continue ;
	    a1 = a2+1 ;
	    a2 = mm->m1 < mm->m2  ? mm->a1[i+1] - mm->m1  : mm->m1 - mm->a1[i+1]  ;
	    freeOutf ("Splicing %d %d %d %d Intron gt_ag Length %d bp\n"
		      , a1, a2, x2, x2+1, a2 - a1 + 1
		      ) ;
	  }
	freeOutf ("\n") ;
      }
  freeOutf ("\n\n") ;
  return TRUE ;
} /* utcAceExport */
      
/*************************************************************************************/
/*************************************************************************************/

static BOOL utcMakeConcensusExport (UTC *utc, MM *mm)
{
  static int nModel = 0 ;
  char *prefix = "KK" ;
  int i ;
  
  /* export this model in gff format */
  freeOutf ("%s_%d", prefix, ++nModel) ;
  freeOutf ("\t%s", dictName (utc->chromDict, mm->chrom)) ;
  if (mm->isDown)
    freeOutf ("\t+\t%d\t%d", mm->m1, mm->m2) ;
  else
    freeOutf ("\t-\t%d\t%d", mm->m2, mm->m1) ;
  if (! utc->cds)
    freeOutf ("\t0\t0") ;
  else
    freeOutf ("\t%d\t%d", mm->c1, mm->c2) ;
  freeOutf ("\t%d\t", mm->nExons) ;
  for (i = 0 ; i <  mm->nExons ; i++)
    freeOutf ("%d,", mm->a1_original[i]) ;
  freeOutf ("\t") ;
  for (i = 0 ; i <  mm->nExons ; i++)
    freeOutf ("%d,", mm->a2_original[i]) ;
  freeOutf ("\t-\t-") ;
  if ((mm->isDown && mm->atgComplete) || (!mm->isDown && mm->stopComplete) )
    freeOutf ("\tcmpl") ;
  else
    freeOutf ("\t-") ;
  if ((!mm->isDown && mm->atgComplete) || (mm->isDown && mm->stopComplete) )
    freeOutf ("\tcmpl") ;
  else
    freeOutf ("\t-") ;
  freeOutf ("\n") ;
  return TRUE ;
} /* utcMakeConcensusExport */

/*************************************************************************************/
/*************************************************************************************/

static int utcMakeConcensusOrder2 (const void *a, const void *b)
{
  const MM *ma = (const MM *) a , *mb = (const MM *) b ;
  int i, nn ;

  nn = ma->chrom - mb->chrom ;
  if (nn) return nn ;
  nn = ma->nExons - mb->nExons ; 
  if (nn) return nn ;
  for (i = 0 ; i < ma->nExons - 1 ; i++)
    {
      nn = ma->a1_original[i+1] - mb->a1_original[i+1] ;
      if (nn) return nn ;
      nn = ma->a2_original[i] - mb->a2_original[i] ;
      if (nn) return nn ;
    }
  return 0 ;
} /* utcMakeConcensusOrder */

/*************************************************************************************/

static int utcMakeConcensusOrder (const void *a, const void *b)
{
  const MM *ma = (const MM *) a , *mb = (const MM *) b ;
  int nn ;

  nn = utcMakeConcensusOrder2 (a, b) ;
  if (nn) return nn ;

  nn = ma->labId - mb->labId ; 
  if (nn) return nn ;
  return 0 ;  
} /* utcMakeConcensusOrder */

/*************************************************************************************/
/* export in .gff format the models present in at least MakeConcensus tracks */
static BOOL  utcMakeConcensus (UTC *utc)
{
  MM *mm, *mm2 ;
  int oldId, ii, jj, nn ;

  if (!utc->aa)
    return FALSE ;
  arraySort (utc->aa, utcMakeConcensusOrder) ;
  
  for (ii = 0, mm = arrp (utc->aa, ii, MM) ; ii < arrayMax (utc->aa); mm++, ii++)
    {
      if (utc->cds && !mm->hasCds)
	continue ;
      nn = 1 ; 
      oldId = mm->labId ;
      for (jj = ii + 1, mm2 = mm + 1 ;  jj < arrayMax (utc->aa); mm2++, jj++)
	{
	  if (utcMakeConcensusOrder2 (mm, mm2))
	    break ;
	  if (mm2->labId == oldId)
	    continue ;
	  nn++ ;  /* i found an identical model from a different track */
	  oldId = mm2->labId ;
	}
      if (nn >= utc->makeConcensus)
	utcMakeConcensusExport (utc, mm) ;
      ii += mm2 - mm - 1 ;
      mm = mm2 - 1 ;
    }
  
  return TRUE ;
} /* utcMakeConcensus */
      
/*************************************************************************************/
/*************************************************************************************/

static BOOL utcScore (MM *mm1, MM *mm2, int *scorep, BOOL *exactp
			  , int *nTrueIntronp, int *nFalseIntronp, int *nMissedIntronp)
{
  int x, y, z, t, i, j, score = 0 ;
  BOOL ok = FALSE, exact = TRUE ;
  
  *nTrueIntronp = 0 ;
  *nFalseIntronp = 0 ;
  *nMissedIntronp = 0 ;
  if (mm1->nExons > 1 &&  mm2->nExons > 1)
    {
      /* study exon starts */
      i = j = 1 ;
      while (i < mm1->nExons && j < mm2->nExons)
	{
	  x = mm1->a1[i] ; y = mm2->a1[j] ;
	  if (x == y)
	    { i++ ; j++ ; score ++ ; ok = TRUE ; }
	  else if (x < y)
	    { i++ ; score-- ; exact = FALSE ; }
	  else if (x > y)
	    { j++ ; score-- ; exact = FALSE ; }
	}
      while (i < mm1->nExons) { i++ ; score-- ; }
      while (j < mm2->nExons) { j++ ; score-- ; }
 
      /* study exon stop */
      i = j = 0 ;
      while (i < mm1->nExons - 1 && j < mm2->nExons - 1)
	{
	  x = mm1->a2[i] ; y = mm2->a2[j] ;
	  if (x == y)
	    { i++ ; j++ ; score ++ ;  ok = TRUE ; }
	  else if (x < y)
	    { i++ ; score-- ; exact = FALSE ; }
	  else if (x > y)
	    { j++ ; score-- ; exact = FALSE ;  }
	}
      while (i < mm1->nExons - 1) { i++ ; score-- ; }
      while (j < mm2->nExons - 1) { j++ ; score-- ; }
    }
  else if (mm1->nExons > 1 &&  mm2->nExons == 1)
    {
      for (i = 0 ; i < mm1->nExons ; i++)
	if (  /* they intersect */
	    (mm1->regionName == mm2->regionName  && mm1->m1 < mm1->m2 && mm2->a1[0] <  mm2->a2[0] && mm1->a1[i] < mm2->a2[0] && mm1->a2[i] > mm2->a1[0]) ||
	    (mm1->regionName == mm2->regionName  && mm1->m1 > mm1->m2 && mm2->a1[0] >  mm2->a2[0] && mm1->a1[i] > mm2->a2[0] && mm1->a2[i] < mm2->a1[0])
	    )
	  { ok = TRUE ; score = -1000 - mm1->nExons ; } /* prefer a hit to a havana with few exons */
    }
  else if (mm1->nExons == 1 &&  mm2->nExons == 1)
    {
      if (
	  (mm1->regionName == mm2->regionName && mm1->m1 < mm1->m2 && mm1->m1 < mm2->m2 && mm1->m2 > mm2->m1) ||
	  (mm1->regionName == mm2->regionName && mm1->m1 > mm1->m2 && mm1->m1 > mm2->m2 && mm1->m2 < mm2->m1)
	  ) /* they intersect */
	{ ok = TRUE ; score = 1000 ; }
    }
  /* count the true and false introns */
  if (! ok)
    (*nFalseIntronp) += mm2->nExons - 1 ;
  else
    {
      /* study introns */
      i = j = 0 ;
      while (i < mm1->nExons - 1 && j < mm2->nExons - 1)
	{
	  x = mm1->a2[i] ; y = mm2->a2[j] ;
	  z = mm1->a1[i+1] ; t = mm2->a1[j+1] ;
	  
	  if (x == y && z == t)
	    { i++ ; j++ ;  (*nTrueIntronp)++ ; }
	  else if (x <= y)
	    {  (*nMissedIntronp)++ ; i++ ; }
	  else if (x > y)
	    {  (*nFalseIntronp)++ ; j++ ; }
	}
      while (i < mm1->nExons - 1) { i++ ;  (*nMissedIntronp)++ ; }
      while (j < mm2->nExons - 1) { j++ ;  (*nFalseIntronp)++ ; }
    }
  if (ok)
    { *scorep = score ; *exactp = exact ; }
  return ok ;
} /* utcScore */

/*************************************************************************************/
/* inside one region, havana = jj0->jj1, ilab = jj2->jj3
 * compute all scalar products and create the best pairs
 */
typedef struct paireStruct { int score, ii, jj, nTrueIntron, nFalseIntron, nMissedIntron, niiExon, njjExon , deltaC2 ; BOOL bad, exact ; } PAIRE ;

static int utcScoreOrder (const void *a, const void *b)
{
  const PAIRE *pa = (const PAIRE *)a, *pb = (const PAIRE *)b ;
  int n ;
  
  n = pa->score - pb->score ;
  if (n) return - n ;
  n = pa->deltaC2 - pb->deltaC2 ;
  return n ;
} /*  utcScoreOrder */

static int utcA3 (UTC *utc, LAB *lab, Array newIntrons, int jj0, int jj1, int jj2, int jj3)
{
  int score = 0 ;
  BOOL exact = FALSE ;
  int ii, jj, nn, nii, njj ;
  int nTrueIntron, nFalseIntron, nMissedIntron ;
  MM *mm1, *mm2 ;
  Array scores = arrayCreate (1000, PAIRE) ;
  BitSet iUsed, jUsed ;
  PAIRE *pp ;

  /* for each lab transcript find best havana twin */
  for (nn = 0, jj = jj2, mm2 = arrp (utc->aa, jj, MM) ; jj < jj3 ; jj++, mm2++) 
    for (ii = jj0, mm1 = arrp (utc->aa, ii, MM) ; ii < jj1 ; ii++, mm1++)
      { 
	if (utcScore (mm1, mm2, &score, &exact, &nTrueIntron, &nFalseIntron, &nMissedIntron))
	  {
	    pp = arrayp (scores, nn++, PAIRE) ;
	    pp->score = score ; pp->exact = exact ;
	    pp->ii = ii ; pp->jj = jj ;
	    pp->nTrueIntron = nTrueIntron ;
	    pp->nFalseIntron = nFalseIntron ;
	    pp->nMissedIntron = nMissedIntron ;
	    pp->niiExon = mm1->nExons ;
	    pp->njjExon = mm2->nExons ;
            pp->deltaC2 = mm1->c2 >= mm2->c2 ? mm1->c2 - mm2->c2 :  mm2->c2 - mm1->c2 ;
	  }
      }
  nii = jj1 - jj0 ;
  njj = jj3 - jj2 ;
  iUsed = bitSetCreate (jj1, 0) ;
  jUsed = bitSetCreate (jj3, 0) ;
  if (nn)
    {
      arraySort (scores, utcScoreOrder) ;
      /* first look for best reciprocal pairs */
      for (ii = 0, pp = arrayp (scores, 0, PAIRE) ; ii < arrayMax (scores) ; ii++, pp++)
	{
	  if (bit (iUsed, pp->ii) || bit (jUsed, pp->jj))
	    continue ;
	  bitSet (iUsed, pp->ii) ;
	  bitSet (jUsed, pp->jj) ;
	  nii-- ; njj-- ;
	  if (pp->niiExon == 1 && pp->njjExon == 1)
	    {
	      lab->nSimilarIntronlessVariant++ ;
	    }
	  else if (pp->niiExon == 1 && pp->njjExon > 1)
	    {
	      lab->nNewVariant++ ;
	      lab->nFalseIntronInNewVariant += pp->njjExon - 1 ;
	    }
	  else if (pp->niiExon > 1 && pp->njjExon == 1)
	    {
	      lab->nNewIntronlessVariant++ ;
	      lab->nMissedVariant++ ; 
	      lab->nMissedIntronInTouchedGene += pp->nMissedIntron ;
	    }
	  else if (pp->niiExon > 1 && pp->njjExon > 1)
	    {
	      if (pp->exact)
		{
		  int myc1, myc2 ;
		  int mybad = 0 ;

		  lab->nExactTranscript++ ;
		  mm1 = arrp (utc->aa, pp->ii, MM) ;
		  mm2 = arrp (utc->aa, pp->jj, MM) ;
		  myc1 = mm2->c1 ; myc2 = mm2->c2 ;
		  if (myc1 && mm1->c1 && mm1->c1 < mm1->c2)
		    {
		      myc1 = mm2->c1 ;
		      myc2 = mm2->c2 + (! utc->cds && ! lab->hasStop ? 3 :0) ;
		      if (mm1->atgComplete && mm1->c1 != myc1)
			mybad |= 1 ; 
		      if (mm1->stopComplete && mm1->c2 != myc2)
			mybad |= 2 ; 
		      if (mm1->c1 < myc2 && mm1->c2 > myc1)
			mybad |= 4 ;
		    }
		  else if (myc1 && mm1->c1 && mm1->c1 > mm1->c2)
		    { 
		      myc1 = mm2->c1 ;
		      myc2 = mm2->c2 - (! utc->cds && ! lab->hasStop ? 3 :0) ; 
		      if (mm1->atgComplete && mm1->c1 != myc1)
			mybad |= 1 ; 
		      if (mm1->stopComplete && mm1->c2 != myc2)
			mybad |= 2 ; 
		      if (mm1->c1 > myc2 && mm1->c2 < myc1)
			mybad |= 4 ;
		    } 
		  if (myc1 && mm1->c1 && (mm1->atgComplete  || mm1->stopComplete))
		    lab->badDenom++ ;
		  if (mybad == 5) lab->badAtg++ ;
		  if (mybad == 6) lab->badStop++ ;
		  if (mybad == 3) lab->badElseWhere++ ;
		  if (mybad == 7) lab->badAtgAndStop++ ;

		  if (0 && (mybad  & 2) && (lab->labId == ACVV || lab->labId == REF))
		    printf ("%s\t%s%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n"
			    , ! mybad ? "__OK__" : "__BAD__"
			    , mm1->isDown ? "++++" : "----"
			    , mybad == 3 ? "***" : "..."
			    , mybad
			    , dictName (utc->dict, mm1->regionName)
			    , dictName (utc->chromDict, mm1->chrom)
			    , mm1->c1, myc1
			    , mm1->c2, myc2
			    , dictName (utc->dict, mm1->mrna)
			    , dictName (utc->dict, mm2->mrna)
			    ) ;
		}
	      else
		lab->nSimilarTranscript++ ;
	      
	      lab->nTrueIntron += pp->nTrueIntron ;
	      lab->nFalseIntron += pp->nFalseIntron ;
	      lab->nMissedIntronInTouchedGene += pp->nMissedIntron ;
	    }
	}
      /* now look for matches from havana to other from method */
      for (ii = 0, pp = arrayp (scores, 0, PAIRE) ; ii < arrayMax (scores) ; ii++, pp++)
	{
	  if (bit (jUsed, pp->jj))
	    continue ;
	  bitSet (jUsed, pp->jj) ; njj-- ; 
	  if (pp->njjExon > 1)
	    {
	      lab->nNewVariant++ ;
	      lab->nFalseIntronInNewVariant += pp->njjExon - 1 ;
	    }
	  else
	    {
	      lab->nNewIntronlessVariant++ ;
	    }
	}
      /* now look for matches to havana from other from method */
      for (ii = 0, pp = arrayp (scores, 0, PAIRE) ; ii < arrayMax (scores) ; ii++, pp++)
	{
	  if (bit (iUsed, pp->ii))
	    continue ;
	  bitSet (iUsed, pp->ii) ; nii-- ;
	  mm1 = arrp (utc->aa, pp->ii, MM) ;
	  if (mm1->nExons > 1)
	    {
	      lab->nMissedVariant++ ; 
	      lab->nMissedIntronInTouchedGene += mm1->nExons - 1 ;
	    }
	}
    }
  for (jj = jj2 ; jj < jj3 ; jj++)
    {
      if (bit (jUsed, jj))
	continue ;
      bitSet (jUsed, jj) ; njj-- ;
      mm2 = arrp (utc->aa, jj, MM) ;
      if (mm2->nExons > 1)
	{
	  int i ;
	  EX *ee ;

	  lab->nInNewGene++ ;
	  lab->nFalseIntronInNewGene += mm2->nExons - 1 ;
	  switch (lab->labId)
	    {
	    default: 
	      for (i = 0 ; i < mm2->nExons - 1 ; i++)
		{
		  ee = arrayp (newIntrons, arrayMax(newIntrons), EX) ;      
		  if (mm2->m1 < mm2->m2)
		    {
		      ee->isDown = TRUE ; 
		      ee->regionName = mm2->regionName ;
		      ee->a1 =  mm2->a2[i+1] ; ee->a2 =  mm2->a1[i] ;
		    }
		  else
		    {
		      ee->isDown = FALSE ;
		      ee->regionName = mm2->regionName ;
		      ee->a1 =  mm2->a2[i+1] ; ee->a2 =  mm2->a1[i] ;
		    }
		}
	      break ;
	    }
	}
      else
	{
	  lab->nIntronlessInNewGene++ ;
	  if (0) printf (">nIntronlessInNewGene: %d, NOP<->%d\n", lab->nIntronlessInNewGene, jj) ;
	}
    }
  for (ii = jj0 ; ii < jj1 ; ii++)
    {
      if (bit (iUsed, ii))
	continue ;
      bitSet (iUsed, ii) ; nii-- ;
      mm1 = arrp (utc->aa, ii, MM) ;
      if (mm1->nExons > 1)
	{
	  lab->nInMissedGene++ ;
	  lab->nMissedIntronInMissedGene += mm1->nExons - 1 ;
	}
    }

  bitSetDestroy (iUsed) ;
  bitSetDestroy (jUsed) ;

  arrayDestroy (scores) ;

  return TRUE ;
} /* utcA3 */
      
/*************************************************************************************/
/* inside one region, collect the havana, then collect each ilab in turn */
static BOOL utcA2 (UTC *utc, Array newIntrons, int jj0, int jj1)
{
  int ii, j0, j1, j2, j3, iLab ;
  MM *mm ;
  LAB *lab ;

  /* locate the reference set */
  j0 = jj1 + 1 ;
  for (ii = jj0, mm = arrp (utc->aa, ii, MM) ; ii < jj1 ; ii++, mm++)
    if (mm->labId == REF)
      { j0 = ii ; break ; }
   for (; ii <= jj1 ; ii++, mm++)
    if (mm->labId != REF)
      break ;
  j1 = ii ;

  /* locate each ilab */
  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      j2 = j3 = 0 ;
      for (ii = jj0, mm = arrp (utc->aa, ii, MM) ; ii < jj1 ; ii++, mm++)
	if (mm->labId == lab->labId)
	  {
	    j2 = ii ; 
	    for (; ii <= jj1 ; ii++, mm++)
	      if (mm->labId != lab->labId)
		break ;
	    j3 = ii ;
	    break ;
	  }
      utcA3 (utc, lab, newIntrons, j0, j1, j2, j3) ;
    }
  
  return TRUE ;
} /* utcA2 */
      
/*************************************************************************************/

static int utcGlobalOrder (const void *a, const void *b)
{
  const MM *ma = (const MM *) a, *mb = (const MM *) b ;
  int n ;

  n = ma->regionName - mb->regionName ;
  if (n) return n ;
  
  n = ma->labId - mb->labId ;
  if (n) return n ;

  n = ma->m1 - mb->m1 ;
  if (n) return n ;
  
  n = ma->m2 - mb->m2 ;
  if (n) return n ;

  return 0 ;
} /* utcGlobalOrder */

/*************************************************************************************/

static int nTranscriptIntronsLabsOrder (const void *a, const void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = la->nTranscript - la->nIntronlessTranscript
    - lb->nTranscript + lb->nIntronlessTranscript ;
  if (n) return -n ;              /* most transcripts with introns */
  n = la->nTranscript - lb->nTranscript ;
  if (n) return -n ;              /* most transcripts */
  n = la->labId - lb->labId ;
  return n ;
} /* nTranscriptIntronsLabsOrder */

/*************************************************************************************/

static int nExactLabsOrder (const void *a, const void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = la->nExactTranscript - lb->nExactTranscript ;
  if (n) return -n ;              /* best nExactTranscript */
  n = la->nSimilarTranscript - lb->nSimilarTranscript ;
  if (n) return -n ;              /* best nSimilarTranscript */
  n = la->labId - lb->labId ;
  return n ;
} /* nExactLabsOrder */

/*************************************************************************************/

static int nExactIntronsLabsOrder (const void *a, const void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = la->nTrueIntron - lb->nTrueIntron ;
  if (n) return -n ;              /* best nExactTranscript */
  n = la->nFalseIntron - lb->nFalseIntron ;
  if (n) return -n ;              /* best nSimilarTranscript */
  n = la->labId - lb->labId ;
  return n ;
} /* nExactLabsOrder */

/*************************************************************************************/

static int labCdsOrder (const void *a, const void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = 
    (la->badDenom - la->badAtgAndStop - la->badStop - la->badElseWhere) -
    (lb->badDenom - lb->badAtgAndStop - lb->badStop - lb->badElseWhere) ;
  if (n) return -n ;              /* best nExactTranscript */
  n = la->badDenom - lb->badDenom ;
  return -n ;
} /* labCdsOrder */

/*************************************************************************************/

static int labNoIntronOrder (const void *a, const void *b)
{
  const LAB *la = (const LAB *)a, *lb = (const LAB *) b ;
  int n ;

  if (la->labId == REF) return -1 ; /* havana first */
  if (lb->labId == REF) return 1 ;
  n = la->nSimilarIntronlessVariant - lb->nSimilarIntronlessVariant ;
  if (n) return -n ;              /* best nExactTranscript */
  n = la->nNewIntronlessVariant - lb->nNewIntronlessVariant ;
  return -n ;
} /* labNoIntronOrder */

/*************************************************************************************/
/* analyse recoprocal matching pairs */
static BOOL  utcA1 (UTC *utc)
{
  int nn ;
  int iLab, ii, jj0 ;
  int oldRegionName = 0 ;
  MM *mm ;
  LAB *lab, *lab0 ;
  Array sortedLabs = 0 ;
  Array newIntrons = arrayHandleCreate (30000, EX, utc->h) ;

  arraySort (utc->aa, utcGlobalOrder) ;
  nn = arrayMax (utc->aa) ;
  array (utc->aa, nn, MM).chrom = -2 ; /* make room */
  arrayMax (utc->aa) = nn ;
  freeOutf ("Found %d mrna\n", nn) ;
  if (!nn) return FALSE ;
 
  for (iLab = 0, lab = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab++, iLab++)
    {
      lab->nTranscript = 0 ;
      lab->nIntronlessTranscript = 0 ;
      lab->nExon  = 0 ;
      lab->nTerminalExon = 0 ;
      lab->nExactTranscript = 0 ;
      lab->nSimilarTranscript = 0 ;
      lab->badAtg = lab->badStop = lab->badAtgAndStop = lab->badDenom = lab->badElseWhere = 0 ;
      lab->nNewVariant = 0 ;
      lab->nInNewGene = 0 ;
      lab->nInMissedGene = 0 ;
      lab->nMissedVariant = 0 ;
      lab->nTrueIntron = 0 ;
      lab->nTrueIntronInNewVariant = 0 ;
      lab->nFalseIntronInNewVariant = 0 ;
      lab->nFalseIntronInNewGene = 0 ;
      lab->nMissedIntronInMissedGene = 0 ;
      lab->nMissedIntronInTouchedGene = 0 ;
      lab->nSimilarIntronlessVariant = 0 ;
      lab->nNewIntronlessVariant = 0 ;
      lab->nIntronlessInNewGene = 0 ;
    }
  /* loop on all chroms */
  for (ii = jj0 = 0, mm = arrp (utc->aa, 0, MM) ; ii <= nn ; ii++, mm++)
    {
      if (mm->regionName != oldRegionName)
	{
	  oldRegionName = mm->regionName ;
	  if (ii) utcA2 (utc, newIntrons, jj0, ii) ;
	  jj0 = ii ;
	}
    }
  nn = arrayMax (utc->aa) ;
  for (iLab = 0, lab0 = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab0++, iLab++)
    {
      lab0->nTranscript = 0 ;  
      for (ii = 0, mm = arrp (utc->aa, 0, MM) ; ii < nn ; ii++, mm++)
	if (mm->labId == lab0->labId)
	  {
	    lab0->nTranscript++ ;
	    lab0->nIntronlessTranscript += mm->nExons > 1 ? 0 : 1 ;
	    lab0->nExon += mm->nExons ;
	    lab0->nTerminalExon += mm->nExons > 1 ? 2 : 1 ;
	  }
    }

  lab0 = utc->refLab ;
  lab0->nExactTranscript = lab0->nTranscript - lab0->nIntronlessTranscript ;
  lab0->nTrueIntron = lab0->nExon - lab0->nTranscript ;
  lab0->nSimilarIntronlessVariant = lab0->nIntronlessTranscript ;
  lab0->badAtg = lab0->badStop = lab0->badAtgAndStop = lab->badElseWhere = 0 ;
  nn = 0 ; sortedLabs = arrayHandleCreate (100, LAB, utc->h) ; 
  for (iLab = 0, lab0 = arrp (utc->labs, 0, LAB) ; iLab < arrayMax (utc->labs) ; lab0++, iLab++)
    {
      lab = arrayp (sortedLabs, nn++, LAB) ;
      *lab = *lab0 ;
    }
  if (! utc->nosort)
    arraySort (sortedLabs, nTranscriptIntronsLabsOrder) ;

  lab0 = utc->refLab ;
  freeOutf ("\n\nCount the transcripts in each method") ;
  freeOutf ("\nProgram"
	    "\tTranscript with intron"
	    "\tIntronless transcript"
	    "\tAll transcript"
	    "\tExon (multi)"
	    "\tTerminal exon"
	    "\tIntron (multi)"
	    "\t%% Terminal exon"
	    ) ;
   for (lab = arrp (sortedLabs, 0, LAB) ; lab->labId ; lab++)
     {
       freeOutf ("\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f" 
		 , dictName (utc->chromDict, lab->labId)
		 , lab->nTranscript - lab->nIntronlessTranscript
		 , lab->nIntronlessTranscript
		 , lab->nTranscript
		 , lab->nExon
		 , lab->nTerminalExon 
		 , lab->nExon - lab->nTranscript
		 , 100.0 * lab->nTerminalExon/ (lab->nExon ? lab->nExon : 1)
		 ) ;
     }

   if (! utc->nosort)
     arraySort (sortedLabs, nExactLabsOrder) ;
   freeOut ("\n\nAnalyse the transcript with introns") ;
   freeOutf ("\nProgram"
	    "\tTranscript identical to %s"
	    "\tTranscript best matching %s"
	    "\tNew transcript in %s gene"
	    "\tNew transcript in new gene"
	    "\tProgram"
	    "\tMissed %s transcript in touched gene"
	    "\tMissed %s transcript in missed gene"
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	    ) ;
   if (utc->cds)
     freeOut ("\tSensitivity"
	      "\tSensitivity eval"
	      "\tSpecificity"
	      "\tSpecificity eval" 
	      "\tVerif"
	      ) ;
   else
     freeOut ("\tSensitivity"
	      "\tSpecificity" 
	      "\tVerif"
	      ) ;

   for (lab = arrp (sortedLabs, 0, LAB) ; lab->labId ; lab++)
     {
       GUIGO *guigo = 0 ;
       
       if (utc->cds)
	 for (guigo = guigos ; guigo->labId ; guigo++)
	   if (guigo->labId == lab->labId)
	     break ;
       if (guigo && guigo->labId != lab->labId)
	 guigo = 0 ;

       freeOutf ("\n%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d" 
		 , dictName (utc->chromDict, lab->labId)
		 , lab->nExactTranscript
		 , lab->nSimilarTranscript
		 , lab->nNewVariant
		 , lab->nInNewGene
		 , dictName (utc->chromDict, lab->labId)
		 , lab->nMissedVariant
		 , lab->nInMissedGene 
		 ) ; 
       if (utc->cds)
	 freeOutf ("\t%.2f\t%.2f\t%.2f\t%.2f"
		   , 100.0 * lab->nExactTranscript/(lab0->nTranscript - lab0->nIntronlessTranscript)
		   , guigo ? guigo->cdsTrSn : 0
		   , 100.0 * lab->nExactTranscript/(lab->nTranscript - lab->nIntronlessTranscript)
		   , guigo ? guigo->cdsTrSp : 0
		   ) ;
       else
	 freeOutf ("\t%.2f\t%.2f"
		   , 100.0 * lab->nExactTranscript/(lab0->nTranscript - lab0->nIntronlessTranscript)
		   , 100.0 * lab->nExactTranscript/(lab->nTranscript - lab->nIntronlessTranscript)
		   ) ;

       freeOutf ("\t%d" 
		 , lab->nExactTranscript + lab->nSimilarTranscript + lab->nMissedVariant + lab->nInMissedGene - lab0->nExactTranscript
		 ) ;
     }

   { 
     int n1, n2 ;
     arraySort (newIntrons,  exonOrder) ;
     n1 = arrayMax (newIntrons) ;
     arrayCompress (newIntrons) ;
     n2 = arrayMax (newIntrons) ;
     freeOutf ("\n\n\nThere are %d unique introns ammong the %d introns in new gene\n\n"
	       , n2, n1) ;
     
   }

   freeOut ("\n\n\nAnalysis of the Start and Stop among CDS with the same intron structure as the reference\n") ;
   freeOut ("\nProgram"
	    "\tSame Start and End"
	    "\tDifferent start"
	    "\tDifferent end"
	    "\tDifferent start and end"
	    "\tElsewhere"
	    "\tReference CDS model with same intron structure"
	    "\tRate"
	    ) ;
   arraySort (sortedLabs, labCdsOrder) ;
   for (lab = arrp (sortedLabs, 0, LAB) ; lab->labId ; lab++)
     {
       freeOutf ("\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f"
		 , dictName (utc->chromDict, lab->labId)
		 , lab->badDenom - lab->badAtgAndStop - lab->badAtg - lab->badStop - lab->badElseWhere
		 , lab->badAtg
		 , lab->badStop
		 , lab->badAtgAndStop
		 , lab->badElseWhere
		 , lab->badDenom
		 , 100.0 * (lab->badAtgAndStop + lab->badElseWhere + lab->badStop) / (lab->badDenom ? lab->badDenom : 1)
		 ) ; 
     }


   freeOutf ("\n\nAnalysis of the transcripts without intron") ;
   freeOutf ("\nProgram"
	     "\tIntronless transcript also found in %s"
	     "\tNew intronless transcript in %s gene"
	     "\tNew intronless transcript in new gene"
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     ) ;
   arraySort (sortedLabs, labNoIntronOrder) ;
   for (lab = arrp (sortedLabs, 0, LAB) ; lab->labId ; lab++)
     {
       freeOutf ("\n%s\t%d\t%d\t%d"
		 ,  dictName (utc->chromDict, lab->labId)
		 , lab->nSimilarIntronlessVariant
		 , lab->nNewIntronlessVariant
		 , lab->nIntronlessInNewGene
		 ) ;
     }

   if (! utc->nosort)
     arraySort (sortedLabs, nExactIntronsLabsOrder) ;
   freeOutf ("\n\nAnalysis of the introns") ;
   freeOutf ("\nProgram"
	     "\tCorrect intron in transcript exact or similar to %s"
	     "\tDifferent intron in transcript exact or similar to %s"
	     "\tIntron in additional variant"
	     "\tIntron in new gene"
	     "\tMissed intron in touched %s gene"
	     "\tMissed intron in missed %s gene" 
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     , dictName (utc->chromDict, utc->refLab->labId)
	     ) ;
   
   for (lab = arrp (sortedLabs, 0, LAB) ; lab->labId ; lab++)
     {
       freeOutf ("\n%s\t%d\t%d\t%d\t%d\t%d\t%d"
		 , dictName (utc->chromDict, lab->labId)
		 , lab->nTrueIntron
		 , lab->nFalseIntron
		 /* 	 , lab->nTrueIntronInNewVariant */
		 , lab->nFalseIntronInNewVariant
		 , lab->nFalseIntronInNewGene 
		 , lab->nMissedIntronInTouchedGene 
		 , lab->nMissedIntronInMissedGene 
		 ) ;
     }
  freeOutf ("\n") ;
      
  return TRUE ;
} /* utcA1 */
      
/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: UCSCtrackCompare [-cds] [-transcripts  | -exons | -introns | -projectedExons | -projectedIntrons | -bp | -global | -projectedGlobal | -regions | -aceExport] [-nosort] \n") ;
  fprintf (stderr, "// Example:  utc -filter -transcripts \n") ;  
  fprintf (stderr, "// -c chromFile : name of the chromosome regions  configuration file\n") ;
  fprintf (stderr, "//       Format :\n") ;
  fprintf (stderr, "//          # lines are comments\n") ;
  fprintf (stderr, "//          Otherwise there are 4 fields per line\n") ;
  fprintf (stderr, "//            1: chrom as spelled in the gff files\n") ;
  fprintf (stderr, "//            2: start of region on that chromosome, in bp\n") ;
  fprintf (stderr, "//            3: end of region on that chromosome, in bp\n") ;
  fprintf (stderr, "//               all data are clipped at the regions boundaries\n") ;
  fprintf (stderr, "//            4: region name used in some of the reports\n") ;
  fprintf (stderr, "//          0 (Zero in ascii) as first char on line stops the parsing\n") ;
  fprintf (stderr, "//       Example :\n") ;
  fprintf (stderr, "//          ENm004  chr22  30128508  31828507\n") ;
  fprintf (stderr, "// -t trackConfigFile : name of the track configuration file\n") ;
  fprintf (stderr, "//       Format :\n") ;
  fprintf (stderr, "//          # lines are comments\n") ;
  fprintf (stderr, "//          Otherwise there are 3 fields per line\n") ;
  fprintf (stderr, "//            1: track name used in some of the reports\n") ;
  fprintf (stderr, "//            2: 0 or 1, indicating if the stop codon is included or not in the CDS\n") ;
  fprintf (stderr, "//            3: name of the track file in gff format downloaded from UCSC\n") ;
  fprintf (stderr, "//          0 (Zero in ascii) as first char on line stops the parsing\n") ;
  fprintf (stderr, "//       Example :\n") ;
  fprintf (stderr, "//          Exonify\t1\texoniphy.txt\n") ;
  fprintf (stderr, "//          0\tDo not analysize the other tracks\n") ;
  fprintf (stderr, "//          Vega\t1\tvegaGeneJune08.txt\n") ;
  fprintf (stderr, "// -cds : restrict to coding part\n") ;
  fprintf (stderr, "// -aceExport : export the models in AceView .ace format\n") ;
  fprintf (stderr, "// -makeConcensus n : exports a gff file with the models shared by n of the tracks\n") ;
  fprintf (stderr, "// -nosort : sort the lab as in the trackConfig file, rather than by score\n") ;
  fprintf (stderr, "// -global :  global statistics for each method\n") ;
  fprintf (stderr, "// -projectedGlobal :  global statistics after projection on the genome\n") ;
  fprintf (stderr, "// -regions : count the transcript per method per region\n") ;
  fprintf (stderr, "// -bp : count the projected base pairs\n") ;
  fprintf (stderr, "// -projectedExons : count the exons after projection on the genome\n") ;
  fprintf (stderr, "// -exons : count the exons\n") ;
  fprintf (stderr, "// -projecteIintrons : count the introns after projection on the genome\n") ;
  fprintf (stderr, "// -introns : count the introns\n") ;
  fprintf (stderr, "// -transcripts : analyse the complete transcripts\n") ;
  fprintf (stderr, "//   no options --> this on line help\n") ;
  fprintf (stderr, "//   full documentation and examples are available from http://aceview.org/Software\n") ;
  exit (1) ;
} /* usage */

/*************************************************************************************/

static void utcTitle (UTC *utc, char *title)
{
  freeOutf ("######################################################\n") ;
  freeOutf (title) ; 
  freeOutf ("%s, %d regions%s\n"
	    , utc->multi	? "multi" : "projected"
	    , arrayMax (utc->chromArray)
	    , utc->cds ? ", restricted to coding part" : ""
	    ) ;
} /* utcTitle */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  int outlevel = 0 ;
  UTC utc ;
  AC_HANDLE h = 0 ;
  const char *cp ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&utc, 0, sizeof (UTC)) ;
  utc.h = h ;

  /* filters on the data */
  utc.cds = getCmdLineOption (&argc, argv, "-cds", 0) ;
  
  /* control the output format */
  utc.nosort = getCmdLineOption (&argc, argv, "-nosort", 0) ;
  
  /* locate the configuration files */
  getCmdLineOption (&argc, argv, "-c", &(utc.chromConfigFile)) ;
  getCmdLineOption (&argc, argv, "-t", &(utc.trackConfigFile)) ;

  /* what do we do */
  if (getCmdLineOption (&argc, argv, "-aceExport", 0))
    {
      utc.aceExport = TRUE ;
    }
  else if (getCmdLineOption (&argc, argv, "-makeConcensus", &cp))
    {
      int n = 0 ;
      char cc ;

      if (cp && sscanf (cp, "%d%c", &n, &cc) == 1 && n>0)
	utc.makeConcensus = n ;
      else
	{
	  fprintf (stderr
		   , "Could not read the integer \'%s\' needed after -makeConcensus"
		   , cp ? cp : "NULL"
		   ) ;
	  exit (1) ;
	}
    }
  else if (getCmdLineOption (&argc, argv, "-projectedExons", 0))
    {
      utc.exons = TRUE ;
      utc.multi = FALSE ;
    }
  else if (getCmdLineOption (&argc, argv, "-projectedIntrons", 0))
    {
      utc.introns = TRUE ;
      utc.multi =  FALSE ;
    }
  else if (getCmdLineOption (&argc, argv, "-exons", 0))
    {
      utc.exons = TRUE ;
      utc.multi = TRUE ;
    }
  else if (getCmdLineOption (&argc, argv, "-introns", 0))
    {
      utc.introns = TRUE ;
      utc.multi = TRUE ;
    }
  else if (getCmdLineOption (&argc, argv, "-global", 0))
    {
      utc.global = TRUE ;
      utc.multi =  TRUE ;
    }
  else if (getCmdLineOption (&argc, argv, "-projectedGlobal", 0))
    {
      utc.global = TRUE ;
      utc.multi =  FALSE ;
    }
  else if (getCmdLineOption (&argc, argv, "-regions", 0))
    utc.regions = TRUE ;
  else if (getCmdLineOption (&argc, argv, "-transcripts", 0))
    {
      utc.pairs = TRUE ; 
      utc.multi =  TRUE ;
    }
  else if (getCmdLineOption (&argc, argv, "-bp", 0))
    utc.bp = TRUE ; 
  else usage () ;

  /* parse the config files : lists of tracks and genome regions */
  utc.chromDict = dictHandleCreate (300, utc.h) ;
  utc.geneDict = dictHandleCreate (300, utc.h) ;
  dictAdd (utc.chromDict, "__toto__", 0) ; /* avoid 0 */
  utc.labs = arrayHandleCreate (64, LAB, utc.h) ;
  utc.chromArray = arrayHandleCreate (64, REGION, utc.h) ;

  if (! utcConfigRegions (&utc) && !utc.aceExport)
    {
      fprintf (stderr
	       , "Could not read the chromConfigFile %s\n"
	       , utc.chromConfigFile) ; 
      exit (1) ;
    }
  if (! utcConfigTracks (&utc))
    {
      fprintf (stderr
	       , "Could not read the trackConfigFile %s\n"
	       , utc.trackConfigFile) ; 
      exit (1) ;
    }

  outlevel = freeOutSetFile (stdout) ;	
  freeOutf ("// start: %s\n", timeShowNow()) ;

  utcGetData (&utc) ;
  freeOutf ("// data parsed, start statistics: %s\n", timeShowNow()) ;

  if (utc.regions)  
    { 
      freeOut ("######################################################\n") ;
      if (utc.cds)
	freeOut ("## Coding transcripts only") ;
      else
	freeOut ("## Coding or non coding transcripts") ;

      freeOut (", 31 test regions followed by 13 training regions\n") ;
      utcRegions (&utc) ;
    }
  else if (utc.global)  
    { 
      utcTitle (&utc, "## Global statistics on all methods") ;
      utcGlobal (&utc) ;
    }
  else if (utc.pairs)  
    { 
      utcTitle (&utc, "## Intron chaining (all introns and exons count multiply)") ;
      utcA1 (&utc) ;
    }
  else if (utc.aceExport)  
    { 
      utcAceExport (&utc) ;
    }
 else if (utc.makeConcensus)
    { 
      utcMakeConcensus (&utc) ;
    }
  else if (utc.bp)  
    { 
      utcTitle (&utc, "## BasePairs ") ;
      utcBp (&utc) ;
    }
  else if (utc.exons) 
    {
      utcTitle (&utc, "## Exons ")  ;
      utcExons (&utc) ;
    }
  else if (utc.introns) 
    {
      utcTitle (&utc, "## Introns ")  ;
      utcExons (&utc) ;
    }
    
  freeOutf ("// done: %s\n", timeShowNow()) ;
  freeOutf ("######################################################\n") ;
  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;
  
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


#include "ac.h"
#include "bitset.h"
#include "dna.h"
#include "mytime.h"
#include "freeout.h"
#include "dict.h"
#include "../wacext/golay.h"

typedef struct probeAlignStruct {
  AC_HANDLE h ;
  const char *probeFileName ;
  const char *targetFileName ;
  const char *ignoredProbeFileName ;
  const char *selectedProbeFileName ;
  const char *previousScoreFileName ;
  const char *outFileName ;
  const char *hitFileName ;
  const char* wordFrequencyConstructTable ;
  const char* wordFrequencyTable ;
  FILE *probeFile, *targetFile, *ignoredProbeFile, *selectedProbeFile, *hitFile;
  Stack probeStack, qualityStack ;
  unsigned long int mask, mask1, mask2 ;
  int dx1, dx2, gap ;
  int target, targetGene ;
  int wffSize ;
  unsigned char *wff ; /* word frequency table of words of length wffSize */
  int jump5 ; /* jump that many base at the 5' or 3' end of the probes */
  int clipAt ; /* clip all probes at that maximal length */
  int exactTargetStart, exactTargetStop; /* no error in this range on the target, used to search exon juntions */
  int overhangLength ; /* used in -splice mode to see that the overhanging prefix exists in the genome */
  int intronMaxLength ; /*  used in -splice mode when veryfying the overhand matches the genome */
  BOOL showTargetPrefix ; /* export 30bp upstream of the alignment */
  BOOL probeSetName ; /* if provided, export the column */
  BOOL clipPolyA ; /* clip terminal polyA */
  BOOL clipPolyT ; /* clip initial polyT */
  BOOL solid ;  /* we are aligning in color space, apply the error correcting procedure */
  BOOL fastQ ; /* probe file is in fastQ format */
  int probeQquality ; /* ignore, if fastq format, right clip until a letter of given quality is reached */
  int fastQbase ; /* base of the quality factors, may be 33 (NCBI) or 64 (ILM/solexa), if set, report the quality of the SNP */
  int nbpQclipped ;   /* how many bp where clipped on the right */
  int nseqQclipped ;  /* in how many sequence tags */
  BOOL stranded ;
  int strandBonus ; /* malus if the the tag is antisense relative to the target */
  BOOL strand, antistrand, antiprobe, aceOutGenomic, aceOutMrna, aceOutNM, silent ;
  BOOL bestHit, slide, getTm, golay, multi, preScanGenome, preScanProbes, strandedTarget ;
  float minTM ;
  unsigned char **sl, **slR, **exitVector ;
  int aliBonusForThisAli, rightVectorClipForThisAli ; /* belong to mm, but it saves memory to have them as globals */
  int errMin, errMax, errRateMax, errCost, bonus, slBonus, 
    maxHit, maxHitR, minEntropy, nExportedHits,
    probeMinLength, probeLengthMin, probeLengthMax, errPosMin, errPosMax, minAli,
    seedOffset, seedShift, seedOffsetUsed, seedError, seedLength, seedLengthRequired ;
  DICT *previousScoreDict, *ignoredProbeDict, *selectedProbeDict, *exportDict, *targetDict ;
  Array probes, letters, dna0, previousScores, exportHits, exportDonors, exportAcceptors, exportIntrons ;
  Associator ass, ass1, ass2, assMulti[32] ;
  Array scores[1024] ;
  int hashPhase, nOligo ;
  Array oligoFrequency ;
  BOOL isShortTarget ;
  BOOL gzi, gzo ;
  ACEIN ai ; ACEOUT ao ;
} PROBEALIGN ;

typedef struct mmStruct {
  int probeLength ;
  int probeLeftClip ;
  int probeRightClip ;
  int probeName ;     /* offsett in PROBEALIGN->probeStack */
  int probeSetName ;  /* offsett in PROBEALIGN->probeStack */
  int probeSequence ; /* offsett in PROBEALIGN->probeStack */
  int probeQuality ;  /* offsett in PROBEALIGN->qualityStack */
  int nHit ;          /* number of hits reported */
  int nPerfectHit ;          /* number of perfect hits reported */
  int bestHit ;       /* minimal N +err+unaligned so far */
  int bestScore ;     /* best score so far */
  int nScore ;        /* offset in scores array for that loop */
  int lastTarget, lastA1, lastA2, lastX1, lastX2 ;
  int entropy ;       /* complexity in bp equivalent */
  int mult ;          /* multiplicity in fastc format */
  int same ;          /* number of letters in common with (mm+1) */ 
  float tm ;          /* TM Maniatis formula */
  int ipx ;
} MM ;

BOOL probeAlignVerifyProbeHit (PROBEALIGN *pp, MM *mm, Array dna, Array dnaR, int pos, int nLoop
			       , const unsigned char *probeDna, int probeLength
			       , BOOL isUp
			       , int *nerrp, int *nNp
			       , int *a1p, int *a2p
			       , int *x1p, int *x2p
			       , int *s1p, int *s2p
			       , int *errLeftPos, int *errLeftPosA
			       , char *errLeftVal, unsigned char *prefix, unsigned char *targetPrefix
			       , int *errRightPos, int *errRightPosA
			       , char *errRightVal, unsigned char *suffix
			       ) ;

void probeAlignExportHit (PROBEALIGN *pp, MM *mm, int nLoop
			  , int a1, int a2, int x1, int x2, int s1, int s2
			  , int nErr, int nN
			  , int errLeftPos, int errLeftPosA, char *errLeftVal, unsigned char *prefix, unsigned char *targetPrefix
			  , int errRightPos, int errRightPosA, char *errRightVal, unsigned char *suffix
			  ) ;

int golayHashProbes (PROBEALIGN *pp) ; 
int golayGetProbeHits (PROBEALIGN *pp, const char *cosmidName, Array dna, BOOL isUp) ;

#ifndef ATITLE_DEF
#define ATITLE_DEF

#include    "../wac/ac.h"
#include "biolog.h"

char* gtCleanUp (const char *ptr) ;
char* gtCleanQuotes (const char *ptr) ; /* keep all upper/lower */
char* gtSquareCleanUp (const char *ptr) ; /* avoids doubling brackets [[...]] */
char* gtLowerCleanUp (const char *ptr) ;
char* gtSetUpper (const char *buf) ;
char* gtSetAllUpper (const char *buf) ;

char *gtGfName (AC_OBJ oGF) ;
char *gtYkName (const char *old) ;
char *gtGeneName (vTXT blkp, GMP *gmp) ;
char *gtGeneTitle (vTXT blkp, GMP *gmp, BOOL showSpecies) ;
char *gtGeneDescriptor (vTXT blkp, GMP *gmp) ;
char *gtGeneFunctionalName (vTXT blkp, GMP *gmp) ;
char *gtGeneFunctionalDescriptor (vTXT blkp, GMP *gmp, AC_OBJ oProduct, BOOL doQuote) ;
AC_OBJ gtPredictedMrna2Product (AC_OBJ oGF, AC_HANDLE h) ;
AC_OBJ gtGene2Product (AC_OBJ oGene, AC_HANDLE h) ;

AC_OBJ gtMrna2Product (AC_OBJ oMrna, AC_HANDLE h) ;
char *gtMrnaSuffix (const char *geneName, const char *mrnaName, AC_HANDLE h)  ;
char *gtMrnaName (vTXT blkp, GMP *gmp) ;
char *gtMrnaTitle (vTXT blkp, GMP *gmp) ;
char *gtProductName (vTXT blkp, GMP *gmp, AC_OBJ oProduct) ;
char *gtProductTitle (vTXT blkp, GMP *gmp) ;
int   gtPredictedMrnaSupport (AC_OBJ oGF, int *nclp) ;
char *gtPredictedMrnaName (vTXT blkp, GMP *gmp, AC_OBJ oMrna) ;
char *gtPredictedMrnaTitle (vTXT blkp, GMP *gmp) ;
char *gtPredictedTrnaTitle (vTXT blkp, GMP *gmp, BOOL showSpecies, char **typep, AC_HANDLE h) ;
char *gtPredictedProductName (vTXT blkp, GMP *gmp, AC_OBJ oProduct, BOOL insist) ; 
char *gtPredictedProductTitle (vTXT blkp, GMP *gmp) ;
char *gtProductMolecularName (vTXT blkp, GMP *gmp, AC_OBJ oProduct) ;

DICT *gtGeneAliases (GMP *gmp, BOOL isLink) ;
char *gtCloneGroup (AC_OBJ oGene, BOOL isLink) ;
char *gtOst (AC_OBJ oGene) ;
char *gtHyman (AC_OBJ oGene) ;
char *gtPiano (AC_OBJ oGene) ;
BOOL  gtIsExperimentalGene (AC_OBJ oGene, AC_OBJ oGF) ;
BOOL  gtIsEssentialGene (AC_OBJ oGene) ;

char *gtReportTitles (AC_DB db, AC_OBJ oProduct) ; /* returns ptr to .ace buffer */

BOOL gmpSuper (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *title) ;
BOOL gmpBlockInit (vTXT blkp, GMP *gmp, BOOL showIt, int page) ;
BOOL gmpBlockStart (vTXT blkp, GMP *gmp, const char *header) ;
BOOL gmpBlockClose (vTXT blkp, GMP *gmp, const char *header) ;
BOOL gmpImportRemainder (vTXT blkp, GMP *gmp, AC_OBJ obj) ;
BOOL gmpChapter (vTXT blkp, GMP *gmp, const char *header, const char *title) ;
BOOL gmpChapter2 (vTXT blkp, GMP *gmp, const char *header, const char *title, const char *bubble) ;
BOOL gmpChapterClose (vTXT blkp, GMP *gmp, const char *header, BOOL register) ;
BOOL gmpSection (vTXT blkp, GMP *gmp, char *header, const char *title) ;
BOOL gmpSection2 (vTXT blkp, GMP *gmp, char *header, const char *title, const char *subTitle) ;
BOOL gmpSubSection (vTXT blkp, GMP *gmp, char *header, const char *title) ;
BOOL gmpCaption (vTXT blkp, GMP *gmp, char *header, const char *title) ;
BOOL gmpHelp (vTXT blkp, GMP *gmp, char *file, const char *title) ;

void gmpJumpPointInit (void) ;
void gmpJumpPointShow (vTXT blkp, GMP *gmp, BOOL showNow, BOOL verbose) ;
void gmpJumpPointDestroy (void) ;

int gmpURL (vTXT blkp, GMP *gmp, char *link, const char *text) ;
BOOL gmpAction (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *action, const char *text) ;
int gmpObjLink (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *text) ;
int gmpObjLinkAnchor (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *anchor, const char *text) ;
int gmpFakeObjLink (vTXT blkp, GMP *gmp, const char *cl, AC_OBJ obj, const char *text) ;

#endif

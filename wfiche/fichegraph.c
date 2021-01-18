/*  File: fiche.c
 *  Author: Vahan Simonyan
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 *  Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk and
 *  Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: View multiple trees in a user configurable way
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 17 03:17 1998 (rd)
 * Created:
 *-------------------------------------------------------------------
 */

/* %W% %G% */

#include "acedb.h"
#include "display.h"
#include "pick.h"
#include "query.h"
#include "dna.h"
#include "peptide.h"
#include "parse.h"

#include    <vstd.h>
#include    <gmleditor.h>
#include    "../wfiche/biolog.h"
#include "../wfiche/gtitle.h"
#include    <vtxt.h>
#include    <vtxt_.h>

/* #define FICHESAVEFEATUREON */

#define clnVar(v_var)   vstrCleanEnds((v_var),(v_var),0," \t\n\r",1)

typedef struct  ficheStruct
{
  KEY key ;
  char    keyName[256] ;
  Graph graph ;
  void *  ficheEditor ;
  int     winW, winH ;
  
  char * ficheBufr ;
  int curMode ;
} ficheWINDOW ;

static vMEX ficheWindows ;

static int     shiftX = 1, shiftY = 5 ;
static float   fichH = 0.9 ;
static float   fichW = 0.95 ;

static gmlEDITCONFIG ficheEditorCfg = 
{
  gmlEDIT_VSCROLL,        /* style */
  PALECYAN,                   /* bgColor */
  PALEYELLOW, DARKBLUE,    /* selectBgColor, selectColor */
  RED, PALEYELLOW,         /* cursorBgColor, cursorColor */
  PALECYAN, DARKBLUE,          /* textBgColor, textColor */
  PALEYELLOW, BLACK, PALEBLUE, /* scrlArrColor, scrlLineColor, scrlElevatorColor */
} ;

void    ficheActionButton (void * arg) ;

enum    {fModeFiche = 0, fSaveFiche, fSaveGenbank, fSaveRefseq, fRemove, fLoadFiche, fPrintFiche, fModeAsnGenbank, fModeAsnRefseq, fModeFlatGB, fModeFlatRS } ;
static COLOUROPT colButtons[] = 
{
  { (ColouredButtonFunc) graphDestroy, assVoid (0) , BLACK, PALEBLUE, "Quit"}, 
  {ficheActionButton, assVoid (fPrintFiche) , BLACK, PALEBLUE, "Print"}, 
  
#ifdef FICHESAVEFEATUREON
  {ficheActionButton, assVoid (fSaveFiche) , BLACK, PALEBLUE, "Save"}, 
#endif
  {ficheActionButton, assVoid (-1) , WHITE, WHITE, "            "}, 
  {ficheActionButton, assVoid (fSaveGenbank) , BLUE, PALEYELLOW, "GenBank"}, 
  {ficheActionButton, assVoid (fSaveRefseq) , BLUE, PALEYELLOW, "RefSeq"}, 
  
  {ficheActionButton, assVoid (-1) , WHITE, WHITE, "                        "}, 
  {ficheActionButton, assVoid (-1) , WHITE, WHITE, "                        "}, 
  {ficheActionButton, assVoid (fModeFiche) , WHITE, BLUE, "Fiche Mode"}, 
  {ficheActionButton, assVoid (-1) , WHITE, WHITE, " "}, 
  {ficheActionButton, assVoid (fModeAsnGenbank) , WHITE, BLUE, "ASN GB Mode"}, 
  {ficheActionButton, assVoid (fModeFlatGB) , WHITE, BLUE, "FF GB"}, 
  {ficheActionButton, assVoid (-1) , WHITE, WHITE, " "}, 
  {ficheActionButton, assVoid (fModeAsnRefseq) , WHITE, BLUE, "ASN RS Mode"}, 
  {ficheActionButton, assVoid (fModeFlatRS) , WHITE, BLUE, "FF RS"}, 
  {ficheActionButton, assVoid (-1) , WHITE, WHITE, " "}, 
  
  {0, assVoid (0) }
} ;

static int ficheGetActiveWindow (void) 
{
  int iw ;
  ficheWINDOW  *fw =  (ficheWINDOW * ) (ficheWindows.buf) ;
  Graph graph  =  graphActive () ;
  
  /* look for our window in out window list */    
  for (iw = 0 ; iw < ficheWindows.cnt ; iw++) 
    {
      if (graph  ==  fw[iw].graph) return iw ;
    }
  return -1 ;
}

static int ficheGetKeyWindow (char * keyName) 
{
  int iw ;
  ficheWINDOW  * fw =  (ficheWINDOW * ) (ficheWindows.buf) ;
  
  for (iw = 0 ; iw<ficheWindows.cnt ; iw++) 
    {
      if (!strcmp (keyName, fw[iw].keyName)) return iw ;
    }
  return -1 ;
}

/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/                                          _/
  _/  DRAWING FUNCTIONS                       _/
  _/                                          _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static void cFicheDraw (void) 
{
  ficheWINDOW  * fw =  (ficheWINDOW * ) (ficheWindows.buf) ;
  int iw = ficheGetActiveWindow () ;
  float y = 1;
  if ( iw  ==  -1) return ;
  
  graphColouredButtons (colButtons, 2, &y, 100) ;
  gmlEditorDraw (fw[iw].ficheEditor) ;
  
  graphTextBounds (fw[iw].winW, fw[iw].winH) ;
  graphRedraw () ;
  pickRememberDisplaySize (CFICHE) ;
}


/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  Contruction and destruction             _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static void cFicheDestroy (void) 
{
  ficheWINDOW  * fw =  (ficheWINDOW * ) (ficheWindows.buf) ;
  int iw = ficheGetActiveWindow () ;
  if (iw == -1) return ;
  
#ifdef FICHESAVEFEATUREON
  if (gmlEditorDirty (fw[iw].ficheEditor, -1) && messQuery ("Do you want to save your changes ?")) 
    ficheActionButton (assVoid (fSaveFiche)) ;
#endif
  
  /* mieg jan2, 2003
     the if (0) introduces i think a memory leak
     but without it i crash if i open a fiche go to gb mode and exit immediatly */

  if (0) messfree (fw[iw].ficheBufr) ; fw[iw].ficheBufr = 0 ;
  fw[iw].graph = 0 ;
  gmlEditorDestroy (fw[iw].ficheEditor) ; fw[iw].ficheEditor = 0 ;

  
  
  /*
    memExclude (&ficheWindows, iw*sizeof (ficheWINDOW) , sizeof (ficheWINDOW)) ;
  */
}

BOOL cFicheDisplay (KEY key, KEY from, BOOL isOldGraph) 
{
  ficheWINDOW sfw, * fw ;
  int iw = ficheGetKeyWindow (name (key)) ; /* get the window for this key */
  int xPos, yPos, xLen, yLen ;
  
  if (iw == -1) 
    { /* the first time we open this fiche */
      fw = &sfw ; memset (fw, 0, sizeof (ficheWINDOW)) ;
      fw->key = key ; /* create a ficheWindow object in the ficheWindow list */
      strcpy (fw->keyName, name (key)) ;
    }
  else 
    {
      fw = & ((ficheWINDOW * ) (ficheWindows.buf)) [iw] ;
      /* do not open if it is already open */
      if (fw->graph) 
	{ graphActivate (fw->graph) ; graphPop () ; return TRUE ; }
    }
  
  /* create a display and register treatment functions  */
  if (! displayCreate (CFICHE)) return FALSE ;
  graphRegister (RESIZE, cFicheDraw) ;
  graphRegister (DESTROY, cFicheDestroy) ;    
  
  graphRetitle (messprintf ("%s: %s", className (key) , name (key)) ) ;
  graphFitBounds (&fw->winW, &fw->winH) ; /* determine the size of Fiche window */
  fw->graph = graphActive () ;  
  graphPop () ;
  
  /* create and configure the fiche editor */    
  yPos = shiftY ; yLen = fw->winH * fichH ; xPos = shiftX ; xLen = fw->winW * fichW ;
  fw->ficheEditor = gmlEditorInit ("Fiche Content", xPos, yPos, xLen, yLen, -1, gmlEDIT_VSCROLL, &ficheEditorCfg) ;
  
  /* add this to the list if it didn't exist before */
  if (iw == -1) 
    {
      vmexAppend (&ficheWindows, fw, sizeof (ficheWINDOW)) ;
      iw = ficheWindows.cnt-1 ;
    }
  
  /* set and redraw the Fiche */
  ficheActionButton (assVoid (fLoadFiche)) ;
  gmlEditorDirty (fw->ficheEditor, 0) ;
  
  fw->curMode = fModeFiche ;
  
  return TRUE ;
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/                                          _/
_/  FICHE-ACTION BUTTONS                    _/
_/                                          _/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

void ficheSaveFields (char * keyName, char * ficheBufr) 
{
  char * titleFromFiche, *commentSection ;
  vSTR bfr ;
  int mrnaSig = 0 ;
  Array dna = 0 ;
  KEYSET ks = 0 ;
  KEY mrna = 0, tg = 0 , gene = 0, gbId = 0 ;
  BOOL res = TRUE ;
   
  vSet0 (bfr) ;
  
  /* prepare to parse changes after edition */
  /* get the mRNA signature */
  ks = query (0, messprintf( "query find mRNA IS \"%s\" AND DNA", keyName)) ;
  mrna = keySetMax (ks) ? keySet (ks, 0) : 0 ;
  if (mrna &&
      (dna = dnaGet (mrna)))
    {
      mrnaSig = hashArray(dna) ;
      arrayDestroy (dna) ;
      tg = keyGetKey (mrna, str2tag("Transcribed_gene")) ;
      if (tg) gene = keyGetKey (tg, str2tag("Gene")) ;
      if (gene) gbId = keyGetKey (gene, str2tag("GenBank_id")) ;
    }
  if (gbId)
    res = messQuery ("This MRNA was already  submitted\nDo you want to resubmit ?") ;
  
  if (res) 
    {
      ficheSectionizeFiche (ficheBufr, &titleFromFiche,  &commentSection) ;
      
      vstrPrintf (&bfr, "Genbank : %s\n", keyName) ;
      vstrPrintf (&bfr, "mRNA_signature %d\n", mrnaSig) ;
      if (gene) vstrPrintf (&bfr, "Gene_id %s\n", name (gene)) ;
      if (titleFromFiche && *titleFromFiche) 
	vstrPrintf (&bfr, "Info Title \"%s\"\n", titleFromFiche) ;
     
      if (commentSection) 
	{
	  vstrPrintf (&bfr, "Fiche \"mRNAfiche:%d\"\n\n"
		     "LongText : \"mRNAfiche:%d\"\n"
		     , mrnaSig, mrnaSig) ;
	  vstrPrintfWrapped (&bfr, commentSection, " \t\n\r", 60, 9, 0) ;
	  vstrPrintf (&bfr, "\n***LongTextEnd***\n") ;
	}

      if (commentSection) * (commentSection-2)  = '\n' ;
    }
    
  parseBuffer (vstrPtr (&bfr), 0) ;
  vstrEmpty (&bfr) ;
}
	
	
char * ficheMRNAConvertFiche (char *keyName, char *ficheBufr, int toMode) 
{
  char * asnBfr = 0, * ffBfr = 0, style, * ffBfr1 = 0 ;
  AC_DB db = 0 ;
  AC_OBJ lObj = 0 ;
  vTXT asnStr = vtxtCreate () ;
  
  if (!ficheBufr) return 0 ;
  /* translate from Fiche mode to ... */
  if (toMode == fModeFiche) return strnew (ficheBufr, 0) ;
  
  style  =  (toMode == fModeAsnRefseq || toMode == fModeFlatRS) ? 'r' : 's' ;
  
  if ((db = ac_open_db ("local", 0)) && 
      (lObj = ac_get_obj (db, "mrna", keyName, 0)))
    { /* work with the copy since it gets sectorized */
      char *ptr, *cpyBufr ;
      if((cpyBufr=strnew (ficheBufr, 0)))
	{
	  GMP *gmp = gmpCreate (db, 0, 0, lObj, 0, 0, style, 'm') ;
	  if ((ptr = strstr(cpyBufr,"\nSEQUENCES")) != 0) 
	    *ptr=0 ;/* cut out the sequence information */
	  asnBfr = fAsnGenerateMRNA (asnStr, gmp, cpyBufr, 1) ;
	  messfree (cpyBufr) ;
	  gmpDestroy (gmp) ;
	}
    }

  /* close the object */
  ac_free (lObj) ; ac_db_close (db) ;
  
  if (toMode == fModeAsnGenbank || toMode == fModeAsnRefseq ) 
    {
      ffBfr = vtxtFileExecution (keyName, "asnbeautifier -s", asnBfr, 0, 0) ;
    }
  else if (toMode == fModeFlatGB || toMode == fModeFlatRS ) 
    {
      /*
      * you might think it was silly to run the beautifier just to hand
      * the input to another program, but asn2gb is not ours and it can't
      * handle input lines with more than 1024 characters.  The beautifier
      * reformats the valid ASN.1 data into a form that asn2gb can read.
      */
      ffBfr1 = vtxtFileExecution (keyName, "asnbeautifier -s", asnBfr, 0, 0) ;
      if (! ffBfr1)
	messout ("asnbeautifier failed\n");
      else
        {
	  ffBfr = vtxtFileExecution (keyName, "asn2gb", ffBfr1, 0, 0) ;
	  messfree (ffBfr1) ;
 	}
    }
  vtxtDestroy (asnStr) ;
  return ffBfr ;
}

static void ficheSynchronizeContent (ficheWINDOW *fwp)
{
  if(!fwp->curMode || fwp->curMode==fModeFiche)
    {
      messfree (fwp->ficheBufr) ;
      fwp->ficheBufr = gmlEditorGetBuffer (fwp->ficheEditor,0) ;
    }
} /* ficheSynchronizeContent	  */
  
#define FICHE_ROOT_DIR "/net/vesta/a/mieg"
void ficheActionButton (void * arg) 
{
  char            fileName[MAXPATHLEN] ;
  ficheWINDOW  *  fw =  (ficheWINDOW * ) (ficheWindows.buf) ;
  int             iw = ficheGetActiveWindow () ;
  char         *  tmpBufr, * ptr ;
  int             toDo = assInt (arg) ;

  if (iw  ==  -1) return ;
  
  switch (toDo) 
    {
    case    fLoadFiche:
#ifdef FICHESAVEFEATUREON
      sprintf (fileName, "%s/YkFiches/%s.fiche", FICHE_ROOT_DIR, fw[iw].keyName) ;
      fw[iw].ficheBufr = vtxtFileGetContent (fileName, 0) ;
#endif
      if (!fw[iw].ficheBufr)
	{ 
	  AC_DB db = 0 ;
	  AC_OBJ lObj = 0 ;
	  GMP *gmp = 0 ;
	  char *ptr ;
	  vTXT buf2 = vtxtCreate () ;
	 
	  if ((db = ac_open_db ("local", 0)) && 
	      (lObj =  ac_get_obj (db, "MRNA",  fw[iw].keyName, 0))) 
	    {
	      gmp = gmpCreate (db, 0, 0, lObj, 0, 0, 's', 'm') ;
	      gtMrnaTitle (buf2, gmp) ;
	      vtxtPrintf (buf2, "\n\n") ;
	      ficheNewMrnaSubmissionComment (buf2, gmp) ;
	      ptr =  vtxtPtr (buf2) ;
	      fw[iw].ficheBufr = strnew (ptr, 0) ;
	      gmpDestroy (gmp) ;
	    }
	  ac_free (lObj) ; ac_db_close (db) ;
	    /* swormView (0, "local", "MRNA", fw[iw].keyName, 's', "FICHE", 0) ; */
	  vtxtDestroy (buf2) ;
	}
				/* set the content of fiche window */
      gmlEditorSetBuffer (fw[iw].ficheEditor, fw[iw].ficheBufr, 0) ;
      gmlEditorDirty (fw[iw].ficheEditor, 0) ;
      
      cFicheDraw () ;
      
      break ;
      
    case    fRemove:
      {
#ifdef FICHESAVEFEATUREON
	sprintf (fileName, "%s/YkFiches/%s.fiche", FICHE_ROOT_DIR, fw[iw].keyName) ;
 	remove (fileName) ; /* remove the ficheFile*/
#endif
      }
      break ;
      
    case    fPrintFiche:
      
      if ((tmpBufr = gmlEditorGetBuffer (fw[iw].ficheEditor, 0)) ) 
	{
	  vSTR blk ;
	  memset (&blk, 0, sizeof (vSTR)) ;
	  vstrPrintfWrapped (&blk, tmpBufr, " \t", 60, 0, 0) ;
	  if ((ptr = vtxtFileExecution (fw[iw].keyName, "lpr", vstrPtr (&blk) , 0, 0)) ) 
	    gmlEditorSetBuffer (fw[iw].ficheEditor, ptr, 0) ; 					
	  messfree (ptr) ; messfree (vstrPtr (&blk)) ; messfree (tmpBufr) ;
	}
      fw[iw].curMode = toDo ;
      
      break ;
      
    case    fSaveFiche:
      ficheSynchronizeContent (&(fw[iw])) ;
      sprintf (fileName, "%s/YkFiches/%s.fiche", FICHE_ROOT_DIR, fw[iw].keyName) ;
      vtxtFileSetContent (fileName, fw[iw].ficheBufr) ;
      gmlEditorDirty (fw[iw].ficheEditor, 0) ;
      
      break ;
      
    case    fSaveGenbank:
      ficheSynchronizeContent (&(fw[iw])) ;
      sprintf (fileName, "%s/YkGenBank/%s.asn", FICHE_ROOT_DIR, fw[iw].keyName) ;

      ficheSaveFields (fw[iw].keyName, fw[iw].ficheBufr) ;
      if ((tmpBufr = ficheMRNAConvertFiche (fw[iw].keyName, fw[iw].ficheBufr, fModeAsnGenbank)) ) 
	{
	  vtxtFileSetContent (fileName, tmpBufr) ;

	  messfree (tmpBufr) ;
	  ficheActionButton (assVoid (fModeAsnGenbank)) ;
	}
      break ;
    
    case    fSaveRefseq:
      ficheSynchronizeContent (&(fw[iw])) ;
      sprintf (fileName, "%s/YkRefSeq/%s.asn", FICHE_ROOT_DIR, fw[iw].keyName) ;
      
      /*				ficheSaveFields (fw[iw].keyName, fw[iw].ficheBufr) ; */
      if ((tmpBufr = ficheMRNAConvertFiche (fw[iw].keyName, fw[iw].ficheBufr, fModeAsnRefseq)) ) 
	{
	  vtxtFileSetContent (fileName, tmpBufr) ;
	  messfree (tmpBufr) ;
	  ficheActionButton (assVoid (fModeAsnRefseq)) ;
      }
      break ;
      
    case    fModeFiche:
    case    fModeAsnGenbank:
    case    fModeAsnRefseq:
    case    fModeFlatGB:
    case    fModeFlatRS: 
      {
	if (fw[iw].curMode == toDo) break ;
	gmlEditorSetStyle (fw[iw].ficheEditor, gmlEDIT_READONLY, (toDo == fModeFiche) ? 0 : 1) ;
	
	ficheSynchronizeContent (&(fw[iw])) ;
	
	if ((tmpBufr = ficheMRNAConvertFiche (fw[iw].keyName, fw[iw].ficheBufr, toDo)) ) 
	  {
	    gmlEditorSetBuffer (fw[iw].ficheEditor, tmpBufr, 0) ;
	    messfree (tmpBufr) ;
	    gmlEditorUpdate (fw[iw].ficheEditor) ;
	  }
	fw[iw].curMode = toDo ;
      }
      break ;
      
    default :break ;
    }
}



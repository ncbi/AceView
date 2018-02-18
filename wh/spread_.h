/*  File: spread_.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
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
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: private header for spreadsheet operations
 *              it completes the opaque SPREAD type
 *              and declares all sprdop + sprddef function used within
 *              the spreadPackage (basic ops + display ops)
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 25 14:09 1999 (fw)
 * Created: Tue Dec 22 16:25:35 1998 (fw)
 *-------------------------------------------------------------------
 */


#ifndef ACEDB_SPREAD__H
#define ACEDB_SPREAD__H

#include "bump.h"
#include "query.h"
#include "table.h"
#include "bs.h"
#include "spread.h"		/* public header */

/*************************************************************/
/* spread - structure containing the spreadsheet definitions */
/*************************************************************/

typedef enum { SHOW_ALL=0, SHOW_MIN, SHOW_MAX, SHOW_AVG , SHOW_VAR, SHOW_SUM, SHOW_MULTI, SHOW_COPY, SHOW_DNA, SHOW_PEP, SHOW_COMPUTE } SpShowType ;

typedef struct colonne_struct
{ 
  BOOL hidden, map ; 
  int extend ;
  SpShowType showType ;
  int colonne ;
  int 
    mandatory , /* 2:required,  1:optional, 0:must be absent */
    realType , /* one of k:KEY, b:BOOL, i:Int, f:Float, t:Text D:DNA, P:Peptide*/
    type ;  /* idem or  count */
  KEY
    classe , /* class if type == 'k' */
    tag ;   /* tag in object I build from */
  BOOL nonLocal ; /* true if this column !has a correct objMark, i.e. average */
  Stack tagStack ;
  Array tagMenu ;
  const char *hiddenp , *optionalp , *typep , *showtypep, *classp , *tagp, *extendp ;

  int dna1, dna2, pep1, pep2 ;
  int position;
  int width ;
  char widthBuffer[12], dna1Buffer[256], dna2Buffer[256] ;

  char subtitleBuffer[60] ;
  char legendBuffer[1024] ;
  
  int from ; /* If 0, this column is the fundamental column of the table */
  char fromBuffer[8] ;
  
  Stack text ;      /* to stack the Text data */
  Array dnaD, dnaR, pep ;

  Stack dnaStack, pepStack ;      /* to stack the Dna or peptide data */
  char conditionBuffer[360] ; /* additional restriction on the new object */
  BOOL condHasParam ;        /* so must reevaluate each time */
  char tagTextBuffer[360] ; /* To edit by hand the tag filed */
  int mark ;
  float centre, mag, min, max ;
  BUMP bump ; 
  Array segs , segps ;  /* segps are used to reoder the segs without loosing the friend info */
  BOOL flip ;
  int mapColBox, subTitleBox, legendBox ;
  QueryRegExp *dnaCondition ;
  COND tagCondition, countCondition, conditionCondition ;
} COL ;


/***************** SPREAD structure *************************/
extern magic_t SPREAD_MAGIC;	/* to verify a SPREAD pointer */

struct SpreadStruct {
  magic_t *magic;        /* == &SPREAD_MAGIC*/

  /*#ifndef NON_GRAPHIC*/
  /*  Graph graph , mapGraph ;*/
  int graph, mapGraph;		/* needs attention - HACK HACK HACK */
  /*#endif*/

  char *fileName ;
  char *dirName ;
  char style ; /* a: ace, p: perl, j: java */
  BOOL showTitle ; /* dump title when printing */
  int 
    tableauBox , activeBox , editModeBox,
    definitionBox , defBoxPerCol , activeDefBox ,
    hideBox, optionalBox, typeBox, paramBox, titleBox, subTitleBox, legendBox,
    widthBox, extendBox, fromBox, dnaBox, pepBox, dna1Box, dna2Box, tagBox, tagTextBox, conditionBox,
    sortColonneBox ;
  int activeColonne , numberDisplayedColonnes, sortColonne, nCol ;
  Array pos2col ;
  char colBuffer[12] ;
  char titleBuffer[301] ;
  char paramBuffer[180] ;
  char sortBuffer[8] ;
  Array activeLine ;
  int curr ;
  FILE *fil ;
  BOOL modif , modified , quitWithoutConfirmation , editMode ;
  Array colonnes ; /* array of definitions for each column */
  Array tableau ; /* array of Arrays of BS-cells */
  Array flags ; /* array of flags for each line */
  Stack comments ;
  BOOL interrupt ;
  /* For spreadMap */
  BOOL zoomAll ;
  unsigned int flag ;
  COL * activeMap ;
  int activeMapBox ;
  Array mapBoxes ;
  int  cursorBox , chromoBox , allButton ;
  /* precompute system, stored as a TABLE* under tKey */
  BOOL precompute ;
  BOOL precomputing ;
  TABLE *table ;
  KEY tKey , tmKey ; 
  Array oldClass ; /* to export class names in -a/-J style */
  BOOL isActiveKeySet ;
  int exportedColumn ;  /* if > 0, export that visible column into exportedKeySet */
  KEYSET exportedKeySet ;
  char *href ;  /* used in -X mode for call backs to the aceview server */
} ;

/************************************************************/

typedef struct sprd1 { BSunit u ;
		       KEY parent, grandParent ; BOOL empty ;
		     } SPCELL ;

typedef struct spread_cell_struct *SC ;
struct spread_cell_struct
{ BOOL empty ;
  KEY key, parent, grandParent ;
  int iCol ;  COL *col ; 
  SC up, /* previous column */
   scFrom, /* from/rightof: has the correct mark */
   scGrandParent ; /* has the correct obj */
  BSunit u ;
  OBJ obj ;
  BSMARK mark ;
} ;

/************************************************************/
/**** functions private within the spreadsheet package ******/
/************************************************************/

/* sprdop.c */
void spreadDestroyCol (COL* c) ;
int spreadOrder(const void *x, const void *y);
void spreadReorder (SPREAD spread);
BOOL spreadCheckConditions(SPREAD spread);
BOOL spreadDoRecompute (SPREAD spread) ;

/* sprddef.c */
void spreadDoSaveDefinitions (SPREAD spread, FILE *f) ;
void spreadDoExportDefinitions (SPREAD spread) ;

#endif /* !ACEDB_SPREAD__H */

/*  Last edited: Dec  4 14:19 1998 (fw) */

/* $Id: dotter.c,v 1.8 2020/05/30 16:50:34 mieg Exp $ */

/*
   DOTTER - Sequence-sequence dotplot in a pixel background image

-------------------------------------------------------------
|  File: dotter.c                                           |
|  Author: Erik Sonnhammer (esr@sanger.ac.uk)               |
|  Copyright (C) E Sonnhammer and R Durbin, 1994            |
-------------------------------------------------------------

   Memory requirements:
   The pixelmap + (alphabetsize+4) x qlen + 2 x slen


   Date   Modification
--------  ---------------------------------------------------
93-09-04  Created
...
94-10-09  Added DNA support, changed crosshair to BoxShift calls
94-10-10  Got rid of *rawdata since we can't store the full n*m matrix
	  Implemented zoomfactor everywhere
	  ScoreVector speedup trick
94-10-12  Reverse-complement for DNA (watson = top/forward strand)
94-10-13  Saving, loading and batch running
94-11-01  Re-calling dotter with middle mouse rectangle.
          Tidy scales even with offset.
	  Zoom factor according to qlen*slen
94-11-11  Introduced header in dotter files [variable(bytes)]:
             version (1), zoom(4), qlen4(4), slen4(4)
          Fixed signed char bug (-1 !> 127) for Sun, Alpha.
          Rewrote inner loop - 3 times faster.
94-11-15  Changed *data to unsigned char.  char* doesn't work for Sun/Sol/OSF.
          Reversed blastx mode
94-12-05  Display of Exons and Introns.
94-12-13  Speedup by better use of scoreVec rows.
95-07-12  [2.0] user-features with start == end drawn as lines.
95-07-13  Finally tracked down malloc problem of graphPixelsRaw().
          This was incorrect in graphRampTool and graphPixles too.
95-08-24  [2.1] Added from command line choosable score matrix file (Both Protein & DNA)
95-08-24  [2.1] Crosshair coordinates in dot-matrix.
95-09-06  [2.2] Karlin-Altschul statistics for estimating best windowsize.
95-09-07  [2.2] Calculation of score -> pixel factor, given MATRIX, and sequences.
95-09-08  [2.2] Dotplot files now with: win, pixelFac, Matrix, Matrixname.
95-11-01  [2.2] Added X options capability.
          Emergency workaround if Karlin/Altschul statistics fail. Limit range to 10-50.
	  Improved zooming in with findCommand().
	  Added Xoptions at zooming.
	  Added MSPlist at zooming.
96-04-16  [2.3] Default usage of 'installed' (private) colormap, by adding "-install"
          to the argument list.  Disable with -i. New graphlib routines by Darren Platt.
96-04-23  [2.3] Added dotterBinary filename to zoom-parameters, since "which" always seems
          to fail from dotter started in a popen, which blocked double zooming.
96-04-24  [2.3] rd changed graphColorSquares to use explicit tints by force majeure.
96-07-18  [2.3] fixed bug (initAlignment called before initCrosshair) that caused crashes.
96-07-24  [2.3] fixed bug in LeftDown that caused box crash when Crosshair wasn't on.
96-07-28  [2.3] Changed to checkmark menus.
                Fixed bug that HSPpixel map didn't get reset before calculation.
96-09-20 (rbrusk): in dotterRedraw(), (re-)initialize  vLineBox and hLineBox for all
			new dotgraph's
97-02-28  [2.4] Fixed bugs in feature drawing routine (incorrect resfac in DNA-protein),
          incorrect parsing of multiple sequence files.
97-03-19  [2.4] Changed findCommand to search path self instead of relying on 'which dotter'.
97-11-19  [2.5] Added featurefile on command line.
          Added series and width drawing of feature segments.
	  Added crosshair-full-window option.
	  Fixed annotation drawing of segment boxes (used to only work on lines).
97-11-19  [2.6] For selfcomparisons, duplication of the mirror image was made default.
                (full selfcomp matrix is now never calculated.  -D has reversed effect).
97-11-19  [2.6] Added selectFeatures Tool to hide/unhide feature series.
98-01-15  [2.7] Changed Feature annotation to msp->desc (no length limit).
                Fixed some bugs in Feature drawing.
		Window sizes automatically to accommodated feature files.


          Pending:

	  change filqueryopen to filqueryeditopen (.dotter) when getting new code.

	  Fix revcomp when zooming in.

	  Add to alignment tool the extent and score of the current window.

	  (Raise initial cutoffs when compressed. ?)

	  HSP drawing bugged for reverse matches and for reversed scale
	  (Huge job to fix ... do only if really necessary)

      	  Score matrix for DNA w/ transversions/transitions...? literature?

-------------------------------------------------------------------------------- */

/* CPU usage profiling on SGI:

   cc dotter.c, no optimization
   pixie -o dotter.pixie dotter
   dotter.pixie -wD -b t seq seq
   prof -pixie -h -only calcWindow dotter dotter.Addrs dotter.Counts
*/

#include "regular.h"

#include "graph.h"
#include "key.h"
#include "../wh/menu.h"
#include "dotter_.h"

#include <ctype.h>
#include <assert.h>

typedef struct featureSeries_ {
    char *name;
    int nr;
    int on;
    int drawOrder;
} FEATURESERIES;

/* tint stuff used to be in graph.h, now local - rd 960524
   NB colours as Erik liked them before Jean tinkered!
   could rename #define's more sensibly now
*/

#define TINT_WHITE      0x00
#define TINT_HIGHLIGHT1 0x01	/* highest priority, dna highlighting */
#define TINT_HIGHLIGHT2 0x02	/* highlight friends */
#define TINT_RED        0x04 
#define TINT_LIGHTGRAY  0x08 
#define TINT_MAGENTA    0x10 
#define TINT_CYAN       0x20 
#define TINT_LIGHTGREEN 0x40 
#define TINT_YELLOW     0x80 

static int tints[8] = { LIGHTRED, MIDBLUE, RED, LIGHTGRAY, 
			MAGENTA, CYAN, LIGHTGREEN, YELLOW } ;

#define LeftBorder 70		/* x-pos of y-axis' 0-position */
#define TopBorder 65		/* y-pos of x-axis' 0-position */
#define MAXALIGNLEN 501
#define MAXLINE 1024

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

char *dotterVersion = "2.6 LICR-CI",
     *Xoptions=0;


extern void colorMap (void) ;
static void setWindow (void) ;
static void initAlignment(void);
static void keyboard (int key);
static void savePlot(void);
static void togglePrintColors(void);
static void setAlnLen(void);
static void drawBlastHSPgray(void);
static void drawBlastHSPline(void);
static void drawBlastHSPlinef(void);
static void dotterRedraw(void);
static void togglePixelmap(void);
static void toggleCrosshair(void);
static void toggleCrosshairPos(void);
static void toggleCrosshairFullscreen(void);
static void toggleGrid(void);
static void clearHSPs(void);
static void initCrosshair(void);
static void loadFeaturesPrompt(void);
static void callDotterParams(void);
static void readmtx(int mtx[24][24], char *mtxfile);
static void mtxcpy(int mtx[24][24], int BLOSUM62[24][24]);
static void DNAmatrix(int mtx[24][24]);
static void dotterPrint(void);
static void Help(void);
static void selectFeatures(void);

#define toggleCrosshairStr    "Crosshair"
#define toggleCrosshairPosStr "Crosshair coordinates"
#define toggleCrosshairFullscreenStr "Crosshair over whole window"
#define toggleGridStr         "Grid"
#define drawBlastHSPgrayStr   "Draw Blast HSPs (gray pixels)"
#define drawBlastHSPlineStr   "Draw Blast HSPs (red lines)"
#define drawBlastHSPlinefStr  "Draw Blast HSPs (colour = f(score))"
#define togglePixelmapStr     "Pixelmap"
#define selectFeaturesStr     "Feature series selection tool"
static MENU dotterMenu ;
static MENUOPT mainMenu[] = {  
  { graphDestroy,     "Quit"} ,
  { Help,             "Help"} ,
  { graphRampTool,    "Greyramp Tool"} ,
  { initAlignment,    "Alignment Tool"} ,
  { dotterPrint,      "Print"} ,
  { toggleCrosshair,   toggleCrosshairStr} ,
  { toggleCrosshairPos,toggleCrosshairPosStr} ,
  { toggleCrosshairFullscreen, toggleCrosshairFullscreenStr} ,
  { toggleGrid,        toggleGridStr} ,
  { menuSpacer,       ""} ,
  { callDotterParams, "Zoom in with parameter control"} ,
  { savePlot,         "Save current plot"} ,
  { loadFeaturesPrompt,"Load features from file"} ,
  { selectFeatures    ,selectFeaturesStr} ,
  { setWindow,        "Change size of sliding window"} ,
  { menuSpacer,       ""} ,
  { drawBlastHSPgray, drawBlastHSPgrayStr} ,
  { drawBlastHSPline, drawBlastHSPlineStr} ,
  { drawBlastHSPlinef,drawBlastHSPlinefStr} ,
  { clearHSPs,        "Remove HSPs"} ,
  { togglePixelmap,   togglePixelmapStr} ,
  { 0, 0 }
} ;

static MENUOPT alnmenu[] = {
  { graphDestroy,  "Quit" },
  {  graphPrint,    "Print" },
  {  togglePrintColors,   "Toggle colours for printing" },
  {   setAlnLen,     "Set Alignment length" },
  {  0, 0 }
};

enum { BLASTNOTHING, BLASTRED, BLASTFUNC, BLASTGREY };

/* Global variables */

static int    MATRIX[24][24],
              i, 
              qlen,	/* query residues */
              qlen4,	/* query Pixels (pixelmap length) */
              slen,	/* subject residues */
              slen4,	/* subject Pixels (pixelmap height) */
              qoffset,  /* Difference between displayed q coord and position in qseq */
              soffset,	/* Difference between displayed s coord and position in sseq  */
              qseqbox, xqseqbox[3], sseqbox, qposbox, sposbox,
              qseqboxCrick, sseqboxCrick, qposboxCrick, sposboxCrick, 
              oldcolor, 
              backgColor=LIGHTGRAY, alnBackgColor=TINT_LIGHTGRAY, 
              vLineBox, hLineBox, CrosshairPosBox,
              win,		/* The length of the sliding window */
              abetsize,		/* The alphabet size (20 or 4) */
              blastp, blastn, blastx,
              zoom,		/* Zoomfactor = 1 ... */
              oldx, oldy,
              selfcomp,
              Display_selfcomp,
              watsonOnly,
              crickOnly,
              rawMapRev[256], 
              resfac,	/* Residue factor. 3 for DNA-Protein, 1 otherwise */
              reversedScale,
              plusmin,
              RightBorder,
              MSPoffset,	/* Difference between real MSP coord and coord stored in MSP */
              CrosshairON = 1,
              CrosshairPosON = 1,
              CrosshairFullscreenON = 0,
              PixelmapON,
              pixelmap_done,
              pixelFac,
              datalen, 
              ALIGNLEN = 125,      /* use an odd number please */
              gridOn = 0,
              BlastHSPsOn = 0,
              BlastMode = BLASTNOTHING,
              printMode = 0,
              HSPgaps = 0,
              fsBoxStart,
              fsAnnRightOn = 0,
              fsAnnBottomOn = 1,
              alignmentInitialized = 0;

Graph  dotGraph=0; /* CI - export for redrawing.  */
static Graph  alnGraph=0, fsGraph=0;

extern Graph  rampGraph;
static UCHAR *data, *HSPpixels=0, rawMap[256];
static char  *qseq, *sseq, qname[NAMESIZE+1], sname[NAMESIZE+1],
              qseqDisp[MAXALIGNLEN+1], xqseqDisp[3][MAXALIGNLEN+1], sseqDisp[MAXALIGNLEN+1],
              qseqDispCrick[MAXALIGNLEN+1], sseqDispCrick[MAXALIGNLEN+1],
              qpos[10], spos[10],
              qcolors[MAXALIGNLEN+1], xqcolors[3][MAXALIGNLEN+1], scolors[MAXALIGNLEN+1],
              qcolorsCrick[MAXALIGNLEN+1], scolorsCrick[MAXALIGNLEN+1],
             *qrevcompl=0, *pepqseq[3],
              CrosshairPosText[100],
              MATRIX_NAME[81] = "",
             *banner;

char         *dotterBinary=0;

static double rectx, oldrectx, recty, oldrecty, exp_res_score, Lambda;
static float  crossx, crossy;
static FILE  *saveFil;
static MSP   *MSPlist=0;	/* List of MSPs - the first object contains data */
static AC_HANDLE  handle;

static char *colorNames[NUM_TRUECOLORS] = {
    "WHITE", 
    "BLACK", 
    "LIGHTGRAY", 
    "DARKGRAY",
    "RED", 
    "GREEN", 
    "BLUE",
    "YELLOW", 
    "CYAN", 
    "MAGENTA",
    "LIGHTRED", 
    "LIGHTGREEN", 
    "LIGHTBLUE",
    "DARKRED", 
    "DARKGREEN", 
    "DARKBLUE",
    "PALERED", 
    "PALEGREEN", 
    "PALEBLUE",
    "PALEYELLOW", 
    "PALECYAN", 
    "PALEMAGENTA",
    "BROWN", 
    "ORANGE", 
    "PALEORANGE",
    "PURPLE", 
    "VIOLET", 
    "PALEVIOLET",
    "GRAY", 
    "PALEGRAY",
    "CERISE", 
    "MIDBLUE"
};

int atob_0[]	/* NEW (starting at 0) ASCII-to-binary translation table */
	= {
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,23,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
	14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,
	NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
	14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,

	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR 
};

int atob[]	/* OLD (starting at 1) ASCII-to-binary translation table  (Inherited from blast) */
	= {
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,24,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
	15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,
	NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
	15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,

	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
};
char aa_btoa[]	/* binary-to-ASCII translation table */
	= "-ARNDCQEGHILKMFPSTWYVBZX*" ;

/*  BLOSUM62 930809

  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  \* */ 
int BLOSUM62[24][24] = {
  { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
 {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
 {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
 {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
 { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
 {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
 {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
 { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
 {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
 {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
 {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
 {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
 {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
 {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
 {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
 { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
 { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
 {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
 {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
 { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
 {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
 {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
 { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
  {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1}
 };



double aafq[20]	/* Amino acid residue frequencies used by S. Altschul */
	= {.081, .057, .045, .054, .015, .039, .061, .068, .022, .057,
	   .093, .056, .025, .040, .049, .068, .058, .013, .032, .067 } ;

#define NN 5

int ntob[] = {
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

int ntob_compl[] = {
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 3,NN, 2,NN,NN,NN, 1,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 0,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 3,NN, 2,NN,NN,NN, 1,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 0,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

/* binary-to-ASCII translation table */
char bton[] = "ACGTN*";

Array fsArr = 0;


/************************/

/* RMEXP will subtract the expected score to remove background noise
   due to biased composition. Gos's idea - not explored yet.

static void rmExp(void){}
*/


static void menuCheck(MENU menu, int mode, int thismode, char *str)
{
    if (mode == thismode)
	menuSetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
    else
	menuUnsetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
}
static void setMenuCheckmarks(void)
{
    menuCheck(dotterMenu, 1, CrosshairON, toggleCrosshairStr);
    menuCheck(dotterMenu, 1, CrosshairPosON, toggleCrosshairPosStr);
    menuCheck(dotterMenu, 1, CrosshairFullscreenON, toggleCrosshairFullscreenStr);
    menuCheck(dotterMenu, 1, gridOn, toggleGridStr);

    menuCheck(dotterMenu, BlastMode, BLASTGREY, drawBlastHSPgrayStr);
    menuCheck(dotterMenu, BlastMode, BLASTRED, drawBlastHSPlineStr);
    menuCheck(dotterMenu, BlastMode, BLASTFUNC, drawBlastHSPlinefStr);

    menuCheck(dotterMenu, 1, PixelmapON, togglePixelmapStr);

    graphNewMenu(dotterMenu);
}


void fatal(char *format, ...)
{
    va_list  ap;

    printf("\nFATAL ERROR: "); 

    va_start(ap, format);
    vprintf(format, ap);
    va_end(ap);

    printf("\n"); 
    exit(1);
}


static void Help(void)
{
    graphMessage (messprintf("\
Dotter, version %s\n\
Copyright (C) Erik Sonnhammer\n\
and Richard Durbin, 1995\n\
\n\n\
Sliding window length = %d\n\
Pixel values = %d x score/residue\n\
Matrix = %s\n\
Zoom (compression) factor = %d\n\
\n\
LEFT MOUSE BUTTON:\n\
  Position crosshair.\n\
\n\
MIDDLE MOUSE BUTTON:\n\
  Zoom in region.\n\
\n\
RIGHT MOUSE BUTTON:\n\
  Menu.\n\
\n\n\
KEYSTROKES:\n\
  Arrow keys: move up/down/left/right\n\
  < > : move along diagonals\n\
  { } : move along reverse diagonals\n\
\n\n\
RESIDUE COLOURS in alignment tool:\n\
  Cyan      = Identical Residue.\n\
  DarkBlue  = Positive Score.\n\
  No colour = Negative score.\n", dotterVersion, win, pixelFac, MATRIX_NAME, zoom));
}


/* Convert pixel values to colormap values */
static void pixels_real2col(UCHAR *pixels)
{
    register int i;
    for (i = 0; i < datalen; i++) pixels[i] = rawMap[pixels[i]];
}

/* Convert colormap values to real values */
static void pixels_col2real(UCHAR *pixels)
{
    /* Note that the real values may not be exactly the same as
       the original, since many real values are converted into
       the same color value */
    register int i;
    for (i = 0; i < datalen; i++) pixels[i] = (UCHAR)rawMapRev[(int)pixels[i]];
}


/* REVERSEBYTES changes order of n bytes at location ptr
*/
void reversebytes(void *ptr, int n)
{ 
  static char copy[256], *cp;
  int  i;

  cp = ptr;
  memcpy(copy, ptr, n);  /* Note: strcpy doesn't work - stops at \0 */
  
  for(i=0; i<n; i++) *cp++ = copy[n-i-1];
}


/* SAVEPLOT writes a 1 byte fileformat first, various params and then the dot-matrix
*/
static void savePlot(void)
{
    static char 
	dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE] ;
    int 
	i, j,
	batch=1,
	MNlen, MNlenSave,
	mtx;
    UCHAR 
	format = 2;


    if (!saveFil) {
	if (!(saveFil = filqueryopen(dirName, fileName, "dotter","w", 
				     "Save as file: (adds .dotter to filename)")))
	    return ;

	pixels_col2real(data);
	batch = 0;
    }

    MNlen = MNlenSave = strlen(MATRIX_NAME);

#ifdef ALPHA
    reversebytes(&zoom, 4);
    reversebytes(&qlen4, 4);
    reversebytes(&slen4, 4);
    reversebytes(&pixelFac, 4);
    reversebytes(&win, 4);
    reversebytes(&MNlenSave, 4);
#endif

    fwrite(&format,    1, 1, saveFil);
    fwrite(&zoom,      1, 4, saveFil);
    fwrite(&qlen4,     1, 4, saveFil);
    fwrite(&slen4,     1, 4, saveFil);
    fwrite(&pixelFac,  1, 4, saveFil); /* New feature of format 2  */
    fwrite(&win,       1, 4, saveFil); /* New feature of format 2  */
    fwrite(&MNlenSave, 1, 4, saveFil); /* New feature of format 2  */
    fwrite(&MATRIX_NAME, 1, MNlen, saveFil); /* New feature of format 2  */
    for (i = 0; i < 24; i++)
	for (j = 0; j < 24; j++) {
	    mtx = MATRIX[i][j];
#ifdef ALPHA
	    reversebytes(&mtx, 4);
#endif
	    fwrite(&mtx, 1, 4, saveFil); /* New feature of format 2  */
	}

#ifdef ALPHA
    reversebytes(&zoom, 4);
    reversebytes(&qlen4, 4);
    reversebytes(&slen4, 4);
    reversebytes(&pixelFac, 4);
    reversebytes(&win, 4);
#endif
    
    fwrite(data, 1, qlen4*slen4, saveFil);

    if (!batch) pixels_real2col(data);

    fclose(saveFil);
    saveFil = 0;
}


static void loadPlot(char *loadfile)
{
    FILE 
	*fil;
    int 
	i, j, n, 
	dotstart,
	MNlen,
	mtx;
    UCHAR 
	format;

    if (!(fil = fopen (loadfile, "r")))
	fatal("Failed to open %s", loadfile);

    if ((fread(&format, 1, 1, fil)) != 1) fatal("reading file %s", loadfile);
    if ((fread(&zoom,   1, 4, fil)) != 4) fatal("reading file %s", loadfile);
    if ((fread(&qlen4,  1, 4, fil)) != 4) fatal("reading file %s", loadfile);
    if ((fread(&slen4,  1, 4, fil)) != 4) fatal("reading file %s", loadfile);
#ifdef ALPHA
    reversebytes(&zoom, 4);
    reversebytes(&qlen4, 4);
    reversebytes(&slen4, 4);
#endif

    if (format == 1) {
	dotstart = 13;
	/* Don't actually know these variables for sure - guess the most common */
	pixelFac = 50;
	win = 25;
    }
    else if (format == 2) 
    {
	if ((fread(&pixelFac, 1, 4, fil)) != 4) fatal("reading file %s", loadfile);
	if ((fread(&win,      1, 4, fil)) != 4) fatal("reading file %s", loadfile);
	if ((fread(&MNlen, 1, 4, fil)) != 4) fatal("reading file %s", loadfile);
#ifdef ALPHA
	reversebytes(&pixelFac, 4);
	reversebytes(&win, 4);
	reversebytes(&MNlen, 4);
#endif
	if ((fread(&MATRIX_NAME, 1, MNlen, fil)) != MNlen) fatal("reading file %s", loadfile);
	MATRIX_NAME[MNlen] = 0;
	for (i = 0; i < 24; i++)
	    for (j = 0; j < 24; j++) {
		if ((fread(&mtx, 1, 4, fil)) != 4) fatal("reading file %s", loadfile);
#ifdef ALPHA
		reversebytes(&mtx, 4);
#endif
		MATRIX[i][j] = mtx;
	    }
	dotstart = MNlen + 2329;
    }
    else 
	fatal("Unknown dotter file format version: %d", format);

    fseek(fil, 0, SEEK_END);
    n = ftell(fil);

    if (n-dotstart != qlen4*slen4)
	fatal("Wrong number of pixels in %s: %d. Expected %d * %-d = %d\n", 
	      loadfile, n, qlen4, slen4, qlen4*slen4);

    datalen = slen4*qlen4;
    data = (UCHAR *)malloc(datalen);
    fseek(fil, dotstart, SEEK_SET);

    if ((n = fread(data, 1, qlen4*slen4, fil)) != qlen4*slen4)
	fatal("Read wrong number of pixels from %s: %d. Expected %d * %-d = %d\n", 
	      loadfile, n, qlen4, slen4, qlen4*slen4);
    
    fclose(fil);

    printf("I read your dotter dotplot file %s.\n", loadfile);
    if (format == 1) 
	printf("It was in the old format 1, so the windowsize and pixel factor had to be guessed.\n");
    fflush(stdout);
}

static int fsorder(FEATURESERIES *x, FEATURESERIES *y)
{
    /*printf("%s - %s : %d\n", x->name, y->name,  strcmp(x->name, y->name));*/
    return strcmp(x->name, y->name);
}

static int fsnrorder(FEATURESERIES *x, FEATURESERIES *y)
{
    if (x->nr < y->nr)
	return -1;
    else if (x->nr > y->nr)
	return 1;
    else return 0;
}

/* Format of feature files:
 *
 * sequence(1=horiz, 2=vert) start end colour annotation....
 *
 * Note: msp->send is used both as indicator of score and that it belongs to a series
 */
static void loadFeatures(FILE* fil)
{
    char line[MAXLINE+1], *cp, color[MAXLINE+1], series[MAXLINE+1];
    MSP *msp = 0 ;
    int format=1;
    static int curnr=0;
    FEATURESERIES fs;

    fs.on = 1;
    if (fsArr)arraySort(fsArr, (void*)fsorder);

    if (MSPlist) {
	msp = MSPlist;  
	while(msp->next) msp = msp->next;
    }

    while (!feof(fil))
    {
	if (!fgets(line, MAXLINE, fil)) break;
	if ((cp = strchr(line, '\n'))) *cp = 0;
	
	if (!strncasecmp(line, "# dotter feature format 2", 25)) {
	    format = 2;
	    continue;
	}

	if (!strlen(line) || *line == '#')
	    continue;
	
	if (!MSPlist) {
	    msp = MSPlist = (MSP *)handleAlloc(0, handle, sizeof(MSP));
	}
	else {
	    msp->next = (MSP *)handleAlloc(0, handle, sizeof(MSP));
	    msp = msp->next;
	}
	
	if (format == 1) {
	    /* frame=sequence nr */
	    if (sscanf(line, "%s%d%d%s", 
		       msp->frame, &msp->qstart, &msp->qend, color) != 4) {
		messout("Error parsing Feature file on line: \"%s\"\n", line);
		return;
	    }
	    fs.name = messalloc(2);
	    strcpy(fs.name, "0");
	    msp->send = 100;	/* Fake series 0 */
	}
	else if (format == 2) {
	    /* frame=sequence nr; send=score (1-100); id=color */
	    if (sscanf(line, "%s%s%d%d%d%s", 
		       msp->frame, series, &msp->send, &msp->qstart, &msp->qend, color) != 6) {
		messout("Error parsing Feature file on line: \"%s\"\n", line);
		return;
	    }
	    fs.name = messalloc(strlen(series)+1);
	    strcpy(fs.name, series);
	}

	
	if (arrayFind(fsArr, &fs, &i, (void*)fsorder)) {
	    msp->sstart = arrp(fsArr, i, FEATURESERIES)->nr;
	    /* printf("Old: %s %d\n", fs.name, msp->sstart);*/
	    messfree(fs.name);
	}
	else {
	    msp->sstart = fs.nr = curnr++;
	    arrayInsert(fsArr, &fs, (void*)fsorder);
	    /* printf("New: %s %d\n", fs.name, fs.nr);*/
	}

	/* Stick color in ->id */
	for (i = 0; i < NUM_TRUECOLORS; i++)
	    if (!strcasecmp(colorNames[i], color)) {
		msp->id = i;
		break;
	    }
	if (i == 16) {
	    messout("Unrecognised colour: %s", color);
	    msp->id = WHITE;
	}
	
	/* Put annotation in desc */
	cp = strstr(line, color);
	for (; *cp && !isspace(*cp); cp++);
	if (isspace(*cp)) cp++;
	msp->desc = messalloc(strlen(cp)+1);
	strcpy(msp->desc, cp);
	
	/* score -3 -> feature */
	msp->score = -3;
    }

    fclose(fil);

    /* Sort feature segment array by number */
    arraySort(fsArr, (void*)fsnrorder);

    menuUnsetFlags(menuItem(dotterMenu, selectFeaturesStr), MENUFLAG_HIDE);

    CrosshairFullscreenON = 1;

    dotterRedraw();
    if (fsGraph) selectFeatures();
}


static void loadFeaturesPrompt(void)
{
    FILE *fil;
    static char dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE];

    if (!(fil = filqueryopen(dirName, fileName, "","r", "Load file: "))) return;

    loadFeatures(fil);
}


/* INITMATRIX calculates the min, max and mean score of the given matrix.
   Archaic and not needed any more.
* /
static void initMatrix (void)
{
    int i, j, x, sum ;

    minScore = maxScore = sum = 0 ;
    for (i = 0 ; i < abetsize ; i++)
	for (j = 0 ; j < abetsize ; j++)
	    { 
		x = MATRIX[i][j];
		if (x < minScore) minScore = x;
		if (x > maxScore) maxScore = x;
		sum += x;
	    }
	
    meanScore = (float)sum/(abetsize*abetsize);

    if (blastn)
	printf("DNA");
    else
	printf("Protein");

    printf(" score matrix: Max score= %d, Min score= %d, Mean score= %.2f\n", 
	   maxScore, minScore, meanScore);
}
*/


/* CALCWINDOW puts the max diagonal for each pixel into *data
 ************************************************************/
static void calcWindow(void)
{
    AC_HANDLE calcHandle;

    register int 
	q, s, qmax,     /* Loop variables */
	*newsum,	/* New sum pointer */
	*oldsum,	/* Old sum pointer */
	*delrow,	/* Pointer to scoreVec, row to remove */
	*addrow,	/* Pointer to scoreVec, row to add */
	dotpos, dotposq, dotposs;
    
    int 
	**scoreVec = 0 ,     /* Array of precalculated scores for qseq residues */
        *sIndex,	/* Array of binary coded sseq residues */
	ql,	        /* query position in local submatrix (of one pixel) */
	sl,	        /* subject position in local submatrix (of one pixel) */
	i, min, sec, frame, pepqseqlen, 
	*zero = 0, *sum1 = 0, *sum2 = 0,
	win2 = win/2, 
	val;

    float 
	dots,	        /* total number of dots (millions) */
	speed;          /* speed in Mdots/seconds */

    UCHAR 
	*dot, dotValue;

    calcHandle = handleCreate();
    
    speed = 17.2;  /* SGI MIPS R10000 (clobber) */
    /* speed = 5.7;  DEC Alpha AXP 3000/700 */
    /* speed = 3.7;  SGI R4400: */
    dots = qlen/1e6*slen;
    if (selfcomp) dots /= 2;
    if (blastn && !(watsonOnly || crickOnly)) dots *= 2;
    if (blastx) dots *= 3;

    min = (int)(dots/speed/60);
    sec = (int)(dots/speed) - min*60;
    printf("%d vs. %d residues => %.2f million dots. ",
	   qlen, slen, dots);

    if (min+sec >= 2) {
	printf("(Takes ");

	if (min)
	    printf("%d:%.2d minutes", min, sec);
	else 
	    printf("%d seconds", sec);
    
	printf(" on an SGI MIPS R10000)");
    }

    printf("\n");
    fflush(stdout);

    /* Initialize lookup tables for faster execution */
    sIndex = (int *) handleAlloc(0, calcHandle, slen * sizeof(int)) ;

    if (blastp)
    {
	zero = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum1 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum2 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));

	scoreVec = (int**) handleAlloc(0, calcHandle, 25*sizeof(int*));
	for (i = 0; i < 25; i++)
	    scoreVec[i] = (int*) handleAlloc(0, calcHandle, qlen*sizeof(int));
	
	for (i = 0; i < 24; i++)
	    for (q = 0; q < qlen; q++)
	      scoreVec[i][q] = MATRIX[i][(int)atob[(int)qseq[q]]-1];
	
	/* Non-protein symbols in scorevector */
	for (q = 0; q < qlen; q++) scoreVec[24][q] = MATRIX[23][23];
	
	for (s = 0; s < slen; s++) sIndex[s] = atob[(int)sseq[s]]-1;
    }
    else if (blastn) 
    {
	zero = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum1 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum2 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));

	scoreVec = (int**) handleAlloc(0, calcHandle, 6*sizeof(int*));
	
	for (i = 0; i < 6; i++)
	    scoreVec[i] = (int*) handleAlloc(0, calcHandle, qlen*sizeof(int));
	/* Fill the scoreVec below depending on which strand */

	    
	for (s = 0 ; s < slen ; s++) sIndex[s] = ntob[(int)sseq[s]];
    }
    
    /* Reset background */
    for (i = 0; i < slen4*qlen4 ; i++) data[i] = 0;


    if (blastx)
    {
	zero = (int *)handleAlloc(0, calcHandle, qlen/3*sizeof(int));
	sum1 = (int *)handleAlloc(0, calcHandle, qlen/3*sizeof(int));
	sum2 = (int *)handleAlloc(0, calcHandle, qlen/3*sizeof(int));

	scoreVec = (int**) handleAlloc(0, calcHandle, 25*sizeof(int *));
	for (i = 0; i < 25; i++)
	    scoreVec[i] = (int*) handleAlloc(0, calcHandle, qlen/3*sizeof(int));
	
	/* Non-protein symbols in scorevector */
	for (q = 0; q < qlen/3; q++) scoreVec[24][q] = MATRIX[23][23];
	
	for (s = 0; s < slen; s++) sIndex[s] = atob[(int)sseq[s]]-1;

	for (frame = 0; frame < 3; frame++)
	{
	    pepqseqlen = strlen(pepqseq[frame]);
	    for (i = 0; i < 24; i++)
		for (q = 0; q < pepqseqlen; q++)
		  scoreVec[i][q] = MATRIX[i][atob[(int)pepqseq[(int)frame][q]]-1];

	    for (i = 0; i<pepqseqlen; i++) {
		sum1[i] = 0;
		sum2[i] = 0;
	    }

	    for (s = 0; s < slen; ++s)
	    {   
		if (s & 1) {
		    newsum = sum1 ;
		    oldsum = sum2 ;
		}
		else {
		    newsum = sum2 ;
		    oldsum = sum1 ;
		}

		if (s >= win) delrow = scoreVec[sIndex[s-win]];
		else delrow = zero;

		addrow = scoreVec[sIndex[s]];
		*newsum = *addrow++;
		
		qmax = min(win, pepqseqlen);
		for (q = 1; q < qmax ; ++q)
		    *++newsum = *oldsum++ + *addrow++;
		
		qmax = pepqseqlen;
		for ( ; q < qmax ; ++q) {
		    *++newsum = *oldsum++ + *addrow++ - *delrow++ ;
		    if (*newsum > 0 && s >= win) 
		    {
			dotposq = (q-win2)/zoom;
			dotposs = (s-win2)/zoom;
			    
			/* Only fill half the submatrix */
			ql = q-win2 - dotposq*zoom;
			sl = s-win2 - dotposs*zoom;
			if (sl >= ql)
			{
			    dotpos = qlen4*dotposs + dotposq;

			    if (dotpos < 0 || dotpos >= datalen) {
				messerror ( "Pixel out of bounds (%d) in blastx: %d\n",
					datalen-1, dotpos);
			    }
			    else {
				val = *newsum * pixelFac / win;
				dotValue = (val > 255 ? 255 : (UCHAR)val);
				dot = &data[dotpos];
				if (dotValue > *dot) *dot = dotValue;
			    }
			}
		    }
		}
	    }
	}
    }

    if (blastp || (blastn && !crickOnly)) 
    {
	if (blastn)
	    for (i = 0; i < 6; i++)
		for (q = 0; q < qlen; q++)
		  scoreVec[i][q] = MATRIX[i][ntob[(int)qseq[q]]];

	for (s = 0; s < slen; ++s)
	{ 
	    if (s & 1) {
		newsum = sum1 ;
		oldsum = sum2 ;
	    }
	    else {
		newsum = sum2 ;
		oldsum = sum1 ;
	    }

	    if (s >= win) delrow = scoreVec[sIndex[s-win]];
	    else delrow = zero;

	    addrow = scoreVec[sIndex[s]];
	    *newsum = *addrow++;

	    qmax = min(win, qlen);
	    for (q = 1; q < qmax ; ++q)
		*++newsum = *oldsum++ + *addrow++;

	    qmax = (selfcomp ? s+1 : qlen);
	    for ( ; q < qmax ; ++q) 
	    {
		*++newsum = *oldsum++ + *addrow++ - *delrow++;
		if (*newsum > 0 && s >= win) 
		{
		    dotposq = (q-win2)/zoom;
		    dotposs = (s-win2)/zoom;
		    
		    /* Only fill half the submatrix */
		    ql = q-win2 - dotposq*zoom;
		    sl = s-win2 - dotposs*zoom;
		    if (sl >= ql)
		    {
			dotpos = qlen4*dotposs + dotposq;

			if (dotpos < 0 || dotpos > datalen-1) {
			    messerror ( "Pixel out of bounds (%d) in blastp/blastn-forw: %d\n", 
				    datalen-1, dotpos);
			}
			else {
			    /* Keep the max dot value of all diagonals in this pixel */
			    val = *newsum * pixelFac / win;
			    dotValue = (val > 255 ? 255 : (UCHAR)val);
			    dot = &data[dotpos];
			    if (dotValue > *dot) *dot = dotValue;
			}
		    }
		}
	    }
	}
    }


    if (blastn && !watsonOnly) 
    {
	if (blastn)
	    for (i = 0; i < 6; i++)
		for (q = 0; q < qlen; q++)
		  scoreVec[i][q] = MATRIX[i][ntob_compl[(int)qseq[q]]];

	for (i = 0; i<qlen; i++) {
	    sum1[i] = 0;
	    sum2[i] = 0;
	}

	for (s = slen-1; s >= 0; --s)
	{ 
	    if (s & 1) {
		newsum = sum1 ;
		oldsum = sum2 ;
	    }
	    else {
		newsum = sum2 ;
		oldsum = sum1 ;
	    }

	    if (s < slen-win) delrow = scoreVec[sIndex[s+win]];
	    else delrow = zero;

	    addrow = scoreVec[sIndex[s]];
	    *newsum = *addrow++;

	    qmax = min(win, qlen);
	    for (q = 1; q < qmax ; ++q)
		*++newsum = *oldsum++ + *addrow++;

	    qmax = (selfcomp ? s+1 : qlen);
	    for ( ; q < qmax ; ++q) {
		*++newsum = *oldsum++ + *addrow++ - *delrow++ ;
		if (*newsum > 0 && s <= slen-win) 
		{
		    dotposq = (q-win2)/zoom;

		    dotposs = (s+win2)/zoom;

		    /* Only fill half the submatrix */
		    ql = q-win2 - dotposq*zoom;

		    /* Set the origin (0,0) to the bottom left corner of submatrix
		       Ugly but correct. Zoom = pixels/submatrix */
		    sl = zoom-1 - (s+win2 - dotposs*zoom);
		    
		    if (sl >= ql)
		    {
			dotpos = qlen4*dotposs + dotposq;
			
			if (dotpos < 0 || dotpos >= datalen) {
			    messerror ( "Pixel out of bounds (%d) in blastn-rev: %d\n",
				    datalen-1, dotpos);
			}
			else {
			    val = *newsum * pixelFac / win;
			    dotValue = (val > 255 ? 255 : (UCHAR)val);
			    dot = &data[dotpos];
			    if (dotValue > *dot) *dot = dotValue;
			}
		    }
		}
	    }
	}
    }

    if (selfcomp && Display_selfcomp) {
	/* Copy mirror image */

	int dotposCopy;
	
	for (s = 0; s < slen4; ++s) { 
	    for (q = 0; q < s ; ++q) {
		
		dotpos = qlen4*s + q;
		dotposCopy = qlen4*q + s;

		if (dotpos < 0 || dotpos >= datalen)
		    messerror ( "Source pixel out of bounds (%d) in mirrorCopy: %d\n",
				datalen-1, dotpos);
		if (dotposCopy < 0 || dotposCopy >= datalen)
		    messerror ( "Destination pixel out of bounds (%d) in mirrorCopy: %d\n",
				datalen-1, dotposCopy);
		data[dotposCopy] = data[dotpos];
	    }
	}
    }

    handleDestroy(calcHandle);
    pixelmap_done = 1;
}


void findScaleUnit (float cutoff, float *u, float *sub)
{
  float unit = *u ;
  float subunit = *u ;

  if (cutoff < 0)
    cutoff = -cutoff ;

  while (unit < cutoff)
    { unit *= 2 ;
      subunit *= 5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2.5000001 ;	/* safe rounding */
      if (unit >= cutoff)
	break ;
      unit *= 2 ;
      subunit *= 2 ;
    }
  subunit /= 10 ;
  if (subunit > *sub)
    *sub = subunit ;
  *u = unit ;
}


static void drawScale(void)
{
    int i, x, y,
    TickStart;			/* First minor tick   */

    float unit=1,		/* Major tick unit */
          subunit=1;		/* First Major tick   */

    graphRectangle(LeftBorder-1, TopBorder-1, LeftBorder-1+qlen4, TopBorder-1+slen4);

    /* QSEQ */
    /* min 100 pixels between major ticks */
    findScaleUnit(100*zoom*resfac, &unit, &subunit);
    y = TopBorder -1;

    graphTextFormat(FIXED_WIDTH);

    if (!qoffset) 
	TickStart = 0;
    else 
	TickStart = subunit*ceil(qoffset/subunit);
    for (i = TickStart; i < qlen+qoffset; i += subunit)
    {
	if (reversedScale) 
	    x = RightBorder - (i-qoffset)/zoom/resfac;
	else 
	    x = LeftBorder-1 + ceil((float)(i-qoffset)/zoom/resfac);

	graphLine(x, y, x, y-5); /* Minor tick */

	if (!(i % (int)unit)) { /* Major tick */
	    graphLine(x, y, x, y-10);
	    graphText(messprintf("%d", i), x-4, y-25);
	}

	if (gridOn && x != LeftBorder-1) {
	    graphColor(LIGHTRED);
	    graphLine(x, y+1, x, y+slen4-1);
	    graphColor(BLACK);
	}
    }
    
    /* SSEQ */
    /* min 50 pixels between major ticks */
    unit = subunit = 1;
    findScaleUnit(50*zoom, &unit, &subunit);
    x = LeftBorder - 1;

    if (!soffset) 
	TickStart = 0;
    else 
	TickStart = subunit*ceil(soffset/subunit);

    for (i = TickStart; i < slen+soffset; i += subunit) 
    {
	y = TopBorder-1 + (i-soffset)/zoom;

	graphLine(x, y, x-5, y); /* Minor tick */

	if (!(i % (int)unit)) { /* Major tick */
	    graphLine(x, y, x-10, y);
	    graphText(messprintf("%7d", i), 0, y-6);
	}

	if (gridOn && y != TopBorder-1) {
	    graphColor(LIGHTRED);
	    graphLine(x+1, y, x+qlen4-1, y);
	    graphColor(BLACK);
	}
    }

    graphTextFormat(PLAIN_FORMAT);
}


/* Return the drawing order of this series */
int fs2order(MSP *msp, int *next)
{
    int i;

    if (!msp->send) return 0;

    if ((i = arrp(fsArr, msp->sstart, FEATURESERIES)->drawOrder)) {
	return i-1;
    }
    else {
	(*next)++;
	arrp(fsArr, msp->sstart, FEATURESERIES)->drawOrder = *next;
	return *next-1;
    }
}


static void drawGenes(MSP *msp)
{
    float 
	sx, ex, sy, ey, midx, midy, x, y, 
	height,		/* box height/width */
	boxHeight=10,
	textHeight,
	oldLinew;
    char strand;
    int i, 
	nextSeries=0;		/* Next series to be drawn */
    float fh; int ph;


    if (fsArr) {
	for (i = 0; i < arrayMax(fsArr); i++)
	    arrp(fsArr, i, FEATURESERIES)->drawOrder = 0;
    }
    graphScreenSize(0, 0, 0, &fh, 0, &ph);
    textHeight = ph/fh;
    boxHeight = ph/fh;

    oldLinew = graphLinewidth(1);


    if (selfcomp || !reversedScale) 
	strand = '+';
    else 
	strand = '-';

    for (; msp; msp = msp->next)
    {    
	height = boxHeight;

	if (msp->score < 0 /* && msp->frame[1] == strand  &&
	    msp->qstart+MSPoffset > qoffset && msp->qend+MSPoffset < qoffset+qlen */ )
	{

	    sx = ceil((float)(msp->qstart+MSPoffset - qoffset)/zoom/resfac);
	    ex = ceil((float)(msp->qend+MSPoffset - qoffset)/zoom/resfac);

	    if (reversedScale) {
		sx = RightBorder - sx;
		ex = RightBorder - ex;
	    }
	    else {
		sx += LeftBorder-1;
		ex += LeftBorder-1;
	    }

	    y = TopBorder + slen4 + 10;
	    if (msp->frame[1] != strand) y += 20;

	    if (msp->score == -1) /* EXON */
	    {
		oldcolor = graphColor(BLUE); 
		graphRectangle(sx, y, ex, y + height);
	    }
	    else if (msp->score == -2) /* INTRON */
	    {
		oldcolor = graphColor(BLUE); 
		midx = 0.5 * (sx + ex) ;
		graphLine (sx, y + height/2, midx, y) ;
		graphLine (ex, y + height/2, midx, y) ;
	    }
	    else if (msp->score == -3) /* FEATURE SEGMENT - COLOURED BOX OR LINE */
	    {

		/* Adjust height to score */
		if (msp->send) {
		    if (!arrp(fsArr, msp->sstart, FEATURESERIES)->on) continue;
		    height = boxHeight*(float)msp->send/100.0;
		}
		
		if (selfcomp || *msp->frame == '2') { /* VERTICAL */
		    x = LeftBorder + qlen4 + 10 + 
			fs2order(msp, &nextSeries)*boxHeight + 
			fs2order(msp, &nextSeries)*textHeight;

		    /* if (msp->qstart > msp->qend) x += 20; */

		    sy = ceil((float)(msp->qstart - soffset)/zoom);
		    ey = ceil((float)(msp->qend - soffset)/zoom);

		    sy += TopBorder-1;
		    ey += TopBorder-1;

		    graphColor(msp->id);
		    if (sy == ey) {
			graphLine(LeftBorder, sy, LeftBorder-2+qlen4, sy);
		    }
		    else {
			graphFillRectangle(x, sy, x + height, ey);
			graphColor(BLACK);
			graphRectangle(x, sy, x + height, ey);
		    }
		    if (fsAnnRightOn && msp->desc) {
			graphColor(BLACK);
			graphText(msp->desc, x+boxHeight+2, sy);
		    }
		}
		
		if (selfcomp || *msp->frame != '2') { /* HORIZONTAL */
		    y = TopBorder + slen4 + 10 +
			fs2order(msp, &nextSeries)*boxHeight + 
			fs2order(msp, &nextSeries)*textHeight;
		    /* if (msp->qstart > msp->qend) y += 20; */
		    
		    graphColor(msp->id);
		    if (sx == ex) {
			graphLine(sx, TopBorder, sx, TopBorder-2+slen4);
		    }
		    else {
			graphFillRectangle(sx, y, ex, y + height);
			graphColor(BLACK);
			graphRectangle(sx, y, ex, y + height);
		    }
		    if (fsAnnBottomOn && msp->desc) {
			graphColor(BLACK);
			graphText(msp->desc, sx, y+boxHeight);
		    }
		}
	    }


	    if (selfcomp) /* Draw the boxes on vertical axes too */
	    {
		sy = ceil((float)(msp->qstart+MSPoffset - qoffset)/zoom);
		ey = ceil((float)(msp->qend+MSPoffset - qoffset)/zoom);
		
		sy += TopBorder-1;
		ey += TopBorder-1;
		
		x = LeftBorder + qlen4 + 10;
		if (msp->frame[1] != strand) x += 20;
		
		if (msp->score == -1) /* EXON */
	        {
		    oldcolor = graphColor(BLUE); 
		    graphRectangle(x, sy, x + height, ey);
		}
		else if (msp->score == -2) /* INTRON */
		{
		    oldcolor = graphColor(BLUE); 
		    midy = 0.5 * (sy + ey) ;
		    graphLine (x + height/2, sy, x, midy) ;
		    graphLine (x + height/2, ey, x, midy) ;
		}
	    }

	    graphColor(oldcolor); 
	    graphLinewidth(oldLinew);
	}
    }
}


static void graphPixelLine(int strength, int sx, int sy, int ex, int ey)
{
    int i, inc, len, dotpos, x, y;
    UCHAR dotValue;

    if (sx < ex) 
	inc = 1;
    else
	inc = -1;

    len = abs(ex - sx +1);
    
    if (strength < 256)
	dotValue = (UCHAR)strength;
    else 
	dotValue = 255;

    for (i = 0; i < len; i++)
    {
	x = sx + i*inc;
	y = sy + i;

	if (x >= 0 && x < qlen4 && sy >= 0 && sy < slen4) 
	{
	    dotpos = qlen4*(y) + x;

	    if (dotpos < 0 || dotpos > datalen-1) {
		messout("Pixel out of bounds (0-%d) in graphPixelLine: %d."
			"Crash imminent.", datalen-1, dotpos);
	    }
	    else
	    {
		if (dotValue > HSPpixels[dotpos]) HSPpixels[dotpos] = dotValue;
	    }
	}
    }
}


static int score2color(int score)
{
    if (score < 75) return DARKRED;
    if (score < 100) return MAGENTA;
    return RED;
}


/*
static int score2width(int score)
{
    if (score < 75) return 1;
    if (score < 100) return 2;
    return 3;
}
*/

int illegalSubmatrix(int sx, int sy, int shift)
{
    int dotposq, dotposs, ql, sl;

    dotposq = (sx+shift)/zoom;
    dotposs = (sy+shift)/zoom;
    
    ql = (sx+shift) - dotposq*zoom;
    sl = (sy+shift) - dotposs*zoom;

    if (sl >= ql) 
	return 0;		/* legal */
    else 
	return 1;		/* illegal */
}


int illegalSubmatrixRev(int sx, int sy, int shift)
{
    int dotposq, dotposs, ql, sl;

    dotposq = (sx-shift)/zoom;
    dotposs = (sy+shift)/zoom;
    
    ql = (sx-shift) - dotposq*zoom;

    /* Set the origin (0,0) to the bottom left corner of submatrix
       Zoom = pixels/submatrix */
    sl = zoom-1 - (sy - dotposs*zoom);

    if (sl >= ql) 
	return 0;		/* legal */
    else 
	return 1;		/* illegal */
}


/* Notes:
 * msp->[qs][start|end] coordinates go [1..n+1] while real dots go [0..n]
 * In other words, msp-coords [6,6] will make a dot at [5,5].
 */
static void drawBlastHSPs(void)
{
    int 
	sx, ex, sy, ey, /* Screen coordinates [0..n] */
	strength, matchlen, shift;
    char 
	*MSPsname;
    MSP 
	*msp;

    /* char   strand = reversedScale ? '-' : '+'; */

    if (BlastMode == BLASTNOTHING) BlastMode = BLASTRED;

    if (BlastMode == BLASTGREY) /* Reset background */
	for (i=0; i < datalen; i++) HSPpixels[i] = 0;

    for (msp = MSPlist; msp;  msp = msp->next)
      {    
	if ((MSPsname = strchr(msp->sname, ':')))
	  MSPsname++;
	else
	  MSPsname = msp->sname;
	
	if (!strcmp(MSPsname, sname))
	  {
	    matchlen = (abs(msp->send - msp->sstart)+1)/zoom;
	    

/* printf("\n%s: %d,%d - %d,%d\n", 
       msp->sname, msp->qstart, msp->sstart, msp->qend, msp->send);
*/

	    if (msp->qstart < msp->qend) {

		sx = ceil((float)(msp->qstart+MSPoffset - qoffset)/resfac) -1;
		sy = msp->sstart - soffset -1;

		/* Check if we're in an illegal part of submatrix
		 * - We only want to compress in legal submatrix parts */
		shift = 0;
		while(illegalSubmatrix(sx, sy, shift)) shift++;

		sx = (sx + shift)/zoom -shift;
		sy = (sy + shift)/zoom -shift;
		
		ex = sx + (matchlen-1);
	    }
	    else {
		sx = ceil((float)(msp->qstart+MSPoffset - qoffset)/resfac) -1;
		sy = msp->sstart - soffset -1;

		/* Check if we're in an illegal part of submatrix
		 * - We only want to compress in legal submatrix parts */
		shift = 0;
		while(illegalSubmatrixRev(sx, sy, shift)) shift++;

		sx = (sx - shift)/zoom + shift;
		sy = (sy + shift)/zoom - shift;
		
		ex = sx - (matchlen-1);

		/* sx = ceil((float)(msp->qstart+MSPoffset - qoffset)/resfac);
		sx = sx/zoom - 1;
		ex = sx - (matchlen-1);
		sy = (msp->sstart - soffset)/zoom - 1;*/

	    }

	    if (reversedScale) {
		sx = RightBorder-LeftBorder - sx;
		ex = RightBorder-LeftBorder - ex;
	    }

	    ey = sy + matchlen-1;
	    
	    if (BlastMode == BLASTRED) {
		graphColor(RED);
		graphLinewidth(1 /* score2width(msp->score) */);

		graphLine(LeftBorder + sx, TopBorder + sy,
			  LeftBorder + ex, TopBorder + ey);
		graphColor(BLACK);
	    }
	    else if (BlastMode == BLASTFUNC) {
		graphColor(score2color(msp->score));
		graphLinewidth(1);

		graphLine(LeftBorder + sx, TopBorder + sy,
			  LeftBorder + ex, TopBorder + ey);
		graphColor(BLACK);
	    }
	    else if (BlastMode == BLASTGREY) {
		/* strength = win*msp->score/abs(msp->send - msp->sstart); */
		strength = msp->score;
		graphPixelLine(strength, sx, sy, ex, ey);
	    }
	    else {
		messout("Unknown BlastMode = bug -> contact Erik");
	    }
	}
    }
    if (BlastMode == BLASTGREY) pixels_real2col(HSPpixels);
}

static void drawBlastHSPgray()
{
    BlastMode = BLASTGREY;

    BlastHSPsOn = 1;
    dotterRedraw();
}
static void drawBlastHSPline()
{
    BlastMode = BLASTRED;

    BlastHSPsOn = 1;
    dotterRedraw();
}
static void drawBlastHSPlinef(void)
{
    BlastMode = BLASTFUNC;

    BlastHSPsOn = 1;
    dotterRedraw();
}
static void clearHSPs(void)
{
    BlastMode = BLASTNOTHING;

    BlastHSPsOn = 0;
    dotterRedraw();
}


static void setQposSpos(int x, int y) /* x and y real sequence coordinates */
{
    if (reversedScale)
	sprintf(qpos, "%d", qlen - x*resfac - 1 + qoffset);
    else 
	sprintf(qpos, "%d", x*resfac + 1 + qoffset);
    
    sprintf(spos, "%d", y+1 + soffset);
}


/* DRAWALIGNMENT draws the alignment of the diagonal at the crosshair
 * Note: x and y are in sequence positions (Real sequence, starting at 0)
 */
static void drawAlignment(int x, int y)
{
    int    i, frame;
    char  *q, *s, *qrc;

    if (x < 0) x = 0;
    else if (x > qlen-1) x = qlen-1;
    if (y < 0) y = 0;
    else if (y > slen-1) y = slen-1;

    oldx = x;
    oldy = y;

    if (!graphActivate(alnGraph)) return;

    q = qseq + x - ALIGNLEN/2;
    s = sseq + y - ALIGNLEN/2;

    setQposSpos(x, y);

    if (blastp || (blastn && !crickOnly)) 
    {
	for (i = 0; i < ALIGNLEN; i++) 
	{
	    qcolors[i] = scolors[i] = alnBackgColor;

	    qseqDisp[i] = sseqDisp[i] = 0;
	    
	    /* Handle sticky ends */
	    if (q+i < qseq || q+i > qseq+qlen-1) {
		qseqDisp[i] = ' ';
		if (s+i >= sseq && s+i <= sseq+slen-1) { 
		    sseqDisp[i] = s[i];
		}
	    }
	    if (s+i < sseq || s+i > sseq+slen-1) {
		sseqDisp[i] = ' ';
		if (q+i >= qseq && q+i <= qseq+qlen-1) { 
		    qseqDisp[i] = q[i];
		}
	    }
	    
	    if (qseqDisp[i] || sseqDisp[i]) continue;
	    
	    qseqDisp[i] = q[i];
	    sseqDisp[i] = s[i];
	    
	    /* Highlight matching residues */
	    if (blastp) {
		if (q[i] == s[i])
		    qcolors[i] = scolors[i] = TINT_CYAN;
		else if ( MATRIX[(int)atob[(int)q[i]]-1 ][(int) atob[(int)s[i]]-1 ] > 0)
		    qcolors[i] = scolors[i] = TINT_HIGHLIGHT2;
	    }
	    else if (blastn && q[i] == s[i])
		qcolors[i] = scolors[i] = TINT_CYAN;

	}

	graphBoxDraw(qseqbox, BLACK, backgColor);
	graphBoxDraw(sseqbox, BLACK, backgColor);
	
	graphBoxDraw(qposbox, BLACK, backgColor);
	graphBoxDraw(sposbox, BLACK, backgColor);
    }

    /* Crick strand - qseq reversed and complemented */

    if (blastn && !watsonOnly)
    {

	qrc = qrevcompl + qlen-1 - x - ALIGNLEN/2;

	for (i = 0; i < ALIGNLEN; i++) 
	{
	    qcolorsCrick[i] = scolorsCrick[i] = alnBackgColor; /* i.e. Background */

	    qseqDispCrick[i] = sseqDispCrick[i] = 0;

	    /* Handle sticky ends */
	    if (qrc+i < qrevcompl || qrc+i > qrevcompl+qlen-1) {
		qseqDispCrick[i] = ' ';
		if (s+i >= sseq && s+i <= sseq+slen-1) { 
		    sseqDispCrick[i] = s[i];
		}
	    }
	    if (s+i < sseq || s+i > sseq+slen-1) {
		sseqDispCrick[i] = ' ';
		if (qrc+i >= qrevcompl && qrc+i <= qrevcompl+qlen-1) { 
		    qseqDispCrick[i] = qrc[i];
		}
	    }
	    
	    if (qseqDispCrick[i] || sseqDispCrick[i]) continue;
	    
	    qseqDispCrick[i] = qrc[i];
	    sseqDispCrick[i] = s[i];
	    
	    /* Highlight matching residues */
	    if (qrc[i] == s[i])
		qcolorsCrick[i] = scolorsCrick[i] = TINT_CYAN;
	}

	graphBoxDraw(qseqboxCrick, BLACK, backgColor);
	graphBoxDraw(sseqboxCrick, BLACK, backgColor);
	
	graphBoxDraw(qposboxCrick, BLACK, backgColor);
	graphBoxDraw(sposboxCrick, BLACK, backgColor);

    }
    
    if (blastx)
    {
	for (i = 0; i < ALIGNLEN; i++) scolors[i] = alnBackgColor;

	for (frame = 0; frame < 3; frame++)
	{
	    for (i = 0; i < ALIGNLEN; i++) 
	    {
		q = pepqseq[frame] + x - ALIGNLEN/2;

		xqcolors[frame][i] = alnBackgColor;

		xqseqDisp[frame][i] = sseqDisp[i] = 0;
	    
		/* Handle sticky ends */
		if (q+i < pepqseq[frame] || q+i > pepqseq[frame]+qlen/3-1) {
		    xqseqDisp[frame][i] = ' ';
		    if (s+i >= sseq && s+i <= sseq+slen-1) { 
			sseqDisp[i] = s[i];
		    }
		}
		if (s+i < sseq || s+i > sseq+slen-1) {
		    sseqDisp[i] = ' ';
		    if (q+i >= pepqseq[frame] && q+i <= pepqseq[frame]+qlen/3-1) { 
			xqseqDisp[frame][i] = q[i];
		    }
		}
	    
		if (xqseqDisp[frame][i] || sseqDisp[i]) continue;
	    
		xqseqDisp[frame][i] = q[i];
		sseqDisp[i] = s[i];
	    
		/* Highlight matching residues */
		if (q[i] == s[i])
		    xqcolors[frame][i] = scolors[i] = TINT_CYAN;
		else if ( MATRIX[(int)atob[(int)q[i]]-1 ][(int) atob[(int)s[i]]-1 ] > 0 ) {
		    xqcolors[frame][i] = TINT_HIGHLIGHT2;
		    if (scolors[i] != TINT_CYAN) scolors[i] = TINT_HIGHLIGHT2;
		}
	    }
	    graphBoxDraw(xqseqbox[frame], BLACK, backgColor);
	}
	graphBoxDraw(qposbox, BLACK, backgColor);

	graphBoxDraw(sseqbox, BLACK, backgColor);
	graphBoxDraw(sposbox, BLACK, backgColor);
    }

    graphActivate(dotGraph);
}


static void togglePrintColors(void)
{
    backgColor = (backgColor == LIGHTGRAY ?  WHITE : LIGHTGRAY);
    alnBackgColor = (alnBackgColor == TINT_LIGHTGRAY ?  TINT_WHITE : TINT_LIGHTGRAY);
    initAlignment();
}


static void dotterPrint(void)
{
    printMode = 1;
    dotterRedraw();

    graphPrint();

    printMode = 0;
    dotterRedraw();
}


static void setAlnLen(void)
{
    if (!graphPrompt ("Give Alignment length", 
		      messprintf("%d", ALIGNLEN), 
		      "iz"))
	return;

    freeint (&ALIGNLEN);

    if (!(ALIGNLEN % 2)) ALIGNLEN++;

    if (ALIGNLEN > MAXALIGNLEN) ALIGNLEN = MAXALIGNLEN;

    initAlignment();
}


static void fsSelFinish(void)
{
    dotterRedraw();
    graphActivate(fsGraph);
    graphPop();
}
static void fsSelAll(void)
{
    graphActivate(fsGraph);
    for (i = 0; i < arrayMax(fsArr); i++) {
	arrp(fsArr, i, FEATURESERIES)->on = 1;
	graphBoxDraw(fsBoxStart+i, WHITE, BLACK);
    }
    fsSelFinish();
}
static void fsSelNone(void)
{
    graphActivate(fsGraph);
    for (i = 0; i < arrayMax(fsArr); i++) {
	arrp(fsArr, i, FEATURESERIES)->on = 0;
	graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
    }
    fsSelFinish();
}
static void fsToggleAnnBottom(void)
{
    fsAnnBottomOn = !fsAnnBottomOn;
    selectFeatures();
    fsSelFinish();
}
static void fsToggleAnnRight(void)
{
    fsAnnRightOn = !fsAnnRightOn;
    selectFeatures();
    fsSelFinish();
}
static void fsSel(int box)
{
    int *on;

    if (box-fsBoxStart < 0 || box-fsBoxStart > arrayMax(fsArr))
	return;
    

    on = &arrp(fsArr, box-fsBoxStart, FEATURESERIES)->on;

    graphActivate(fsGraph);


    if (*on) {
	*on = 0;
	graphBoxDraw(box, BLACK, WHITE);
    }
    else {
	*on = 1;
	graphBoxDraw(box, WHITE, BLACK);
    }

    fsSelFinish();
}
    

static void selectFeatures(void)
{
    int i, box;

    if (!graphActivate(fsGraph))
    {
	fsGraph = graphCreate (TEXT_SCROLL, "Dotter - Feature Series Selection Tool", 0, 0, 0.4, 0.4);
    }
    graphPop();
    graphRegister(PICK, fsSel);

    graphText("Pick buttons to select/unselect series", 1, 1);

    graphButton("All ", fsSelAll, 1, 3);

    box = graphButton("Show bottom series annotation", fsToggleAnnBottom, 8, 3);
    if (!fsAnnBottomOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);

    box = graphButton("Show right series annotation ", fsToggleAnnRight, 8, 4.5);
    if (!fsAnnRightOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);

    fsBoxStart = 1+graphButton("None", fsSelNone, 1, 4.5);

    if (fsArr) {
	for (i = 0; i < arrayMax(fsArr); i++) {
	    float 
		y=7+i*1.5,
		margin = 0.1;
	    box = graphBoxStart();
	    graphText(arrp(fsArr, i, FEATURESERIES)->name, 1, y);
	    graphRectangle(1-margin, y-margin, 
			   1+margin+strlen(arrp(fsArr, i, FEATURESERIES)->name), 
			   y+1+margin);
	    graphBoxEnd();
	    if (!arrp(fsArr, i, FEATURESERIES)->on) {
		graphBoxDraw(box, BLACK, WHITE);
	    }
	    else {
		graphBoxDraw(box, WHITE, BLACK);
	    }
	}
    }

    graphRedraw();
}


static void initAlignment(void)
{
    float x, y;
    int frame, height = 0 ;

    {
	/* static int warned=0;
	if (!warned) messout("The residue colours of the Aligment tool have been corrupted by Jean."
			     " Send complaints to mieg@kaa.crbm.cnrs-mop.fr");
	warned=1;*/
    }

    if (!CrosshairON) {
	messout("Turn on the crosshair !");
	return;
    }
    
    graphActivate(dotGraph);
    graphBoxDim (vLineBox, &x, 0, 0, 0);
    graphBoxDim (hLineBox, 0, &y, 0, 0);

    if (!graphActivate(alnGraph))
    {
	alnGraph = graphCreate (TEXT_SCROLL, "Dotter - Alignment Tool", 0, 0, 1.25, 0.25);
	graphMenu(alnmenu);
	graphRegister(KEYBOARD, keyboard);
    }

    graphPop();
    graphClear();
    
    if (blastp)
	height = 7;
    else if (blastn) {
	if (watsonOnly) height = 7;
	else height = 14;
    }
    else if (blastx) 
	height = 9;

    graphTextBounds(MAXALIGNLEN, height);
    graphColor(backgColor); 
    graphRectangle(0, 0, 1000, 1000);
    graphColor(BLACK);
	
    for (i = 0; i < MAXALIGNLEN; i++)
	qseqDisp[i] = 	sseqDisp[i] =
	    xqseqDisp[0][i] = xqseqDisp[1][i] = xqseqDisp[2][i] =
		qseqDispCrick[i] = sseqDispCrick[i] = 0;

    /* Maybe do this instead * /
    qseqDisp = handleAlloc(0, handle, MAXALIGNLEN+1);
    xqseqDisp[0] = handleAlloc(0, handle, MAXALIGNLEN+1);
    xqseqDisp[2] = handleAlloc(0, handle, MAXALIGNLEN+1);
    xqseqDisp[3] = handleAlloc(0, handle, MAXALIGNLEN+1);
    sseqDisp = handleAlloc(0, handle, MAXALIGNLEN+1);
    qseqDispCrick = handleAlloc(0, handle, MAXALIGNLEN+1);
    sseqDispCrick = handleAlloc(0, handle, MAXALIGNLEN+1);
    */

    if (blastx)
    {
	/* Draw static features */
	graphText(messprintf("%s 1:", qname), 0.5, 3);
	graphText(messprintf("%s 2:", qname), 0.5, 4);
	graphText(messprintf("%s 3:", qname), 0.5, 5);
	graphText(messprintf("%s:", sname), 0.5, 6);
	graphLine(NAMESIZE+4 + (float)ALIGNLEN/2, 2, NAMESIZE+4 + (float)ALIGNLEN/2, 8);
	
	/* alignment boxes */
	for (frame = 0; frame < 3; frame++)
	{
	    xqseqbox[frame] = graphBoxStart();
	    graphColorSquares (xqcolors[frame], NAMESIZE+4, 3+frame, ALIGNLEN, 1, tints);
	    graphTextPtr     (xqseqDisp[frame], NAMESIZE+4, 3+frame, ALIGNLEN);
	    graphBoxEnd ();
	}

	sseqbox = graphBoxStart();
	graphColorSquares (scolors, NAMESIZE+4, 6, ALIGNLEN, 1, tints);
	graphTextPtr     (sseqDisp, NAMESIZE+4, 6, ALIGNLEN);
	graphBoxEnd ();

	/* coordinate boxes */
	qposbox = graphBoxStart();
	graphTextPtr (qpos, NAMESIZE + ALIGNLEN/2+2, 1, 10);
	graphBoxEnd ();
	sposbox = graphBoxStart();
	graphTextPtr (spos, NAMESIZE + ALIGNLEN/2+2, 8, 10);
	graphBoxEnd ();
    }

    if (blastp || (blastn && !crickOnly)) 
    {
	/* Draw static features */
	graphText(messprintf("%s:", qname), 0.5, 3);
	graphText(messprintf("%s:", sname), 0.5, 4);
	graphLine(NAMESIZE+2 + (float)ALIGNLEN/2, 2, NAMESIZE+2 + (float)ALIGNLEN/2, 6);
	
	/* alignment boxes */
	qseqbox = graphBoxStart();
	graphColorSquares (qcolors, NAMESIZE+2, 3, ALIGNLEN, 1, tints);
	graphTextPtr     (qseqDisp, NAMESIZE+2, 3, ALIGNLEN);
	graphBoxEnd ();
	sseqbox = graphBoxStart();
	graphColorSquares (scolors, NAMESIZE+2, 4, ALIGNLEN, 1, tints);
	graphTextPtr     (sseqDisp, NAMESIZE+2, 4, ALIGNLEN);
	graphBoxEnd ();
	    
	/* coordinate boxes */
	qposbox = graphBoxStart();
	graphTextPtr (qpos, NAMESIZE + ALIGNLEN/2+2, 1, 10);
	graphBoxEnd ();
	sposbox = graphBoxStart();
	graphTextPtr (spos, NAMESIZE + ALIGNLEN/2+2, 6, 10);
	graphBoxEnd ();
    }

    if (blastn && !watsonOnly)
    {
	/* Draw static features */
	graphText("RevComp:", 0.5, 10);
	graphText(messprintf("%s:", sname), 0.5, 11);
	graphLine(NAMESIZE+2 + (float)ALIGNLEN/2, 9, NAMESIZE+2 + (float)ALIGNLEN/2, 13);
	
	/* alignment boxes */
	qseqboxCrick = graphBoxStart();
	graphColorSquares (qcolorsCrick, NAMESIZE+2, 10, ALIGNLEN, 1, tints);
	graphTextPtr     (qseqDispCrick, NAMESIZE+2, 10, ALIGNLEN);
	graphBoxEnd ();
	sseqboxCrick = graphBoxStart();
	graphColorSquares (scolorsCrick, NAMESIZE+2, 11, ALIGNLEN, 1, tints);
	graphTextPtr     (sseqDispCrick, NAMESIZE+2, 11, ALIGNLEN);
	graphBoxEnd ();
	
	/* coordinate boxes */
	qposboxCrick = graphBoxStart();
	graphTextPtr (qpos, NAMESIZE + ALIGNLEN/2+2, 8, 10);
	graphBoxEnd ();
	sposboxCrick = graphBoxStart();
	graphTextPtr (spos, NAMESIZE + ALIGNLEN/2+2, 13, 10);
	graphBoxEnd ();
    }    
    
    graphBoxDraw(0, backgColor, backgColor);
    graphRedraw();
    drawAlignment((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
}


static void toggleCrosshair(void)
{
    if (CrosshairON) {
	CrosshairON  = 0;
    }
    else {
	CrosshairON  = 1;
    }

    dotterRedraw();
}

static void toggleCrosshairPos(void)
{
    CrosshairPosON = (CrosshairPosON ? 0 : 1);
    dotterRedraw();
}
static void toggleCrosshairFullscreen(void)
{
    CrosshairFullscreenON = !CrosshairFullscreenON;
    dotterRedraw();
}


static void drawCrosshair(float x, float y)
{
    if (CrosshairON) {
	if (CrosshairFullscreenON) {
	    graphBoxShift (hLineBox, 0.0, y);
	    graphBoxShift (vLineBox, x, 0.0);
	}
	else {
	    graphBoxShift (hLineBox, LeftBorder, y);
	    graphBoxShift (vLineBox, x, TopBorder);
	}

	if (CrosshairPosON) {
	    sprintf(CrosshairPosText, "%s, %s", qpos, spos);
	    graphBoxShift (CrosshairPosBox, x+10, y+10);
	}
    }

    crossx = x;
    crossy = y;
}


static void boundaries(double *x, double *y)
{
    if (*x < LeftBorder) *x = LeftBorder;
    else if (*x > LeftBorder + qlen4) *x = LeftBorder + qlen4-1;
    if (*y < TopBorder) *y = TopBorder;
    else if (*y > TopBorder + slen4) *y = TopBorder + slen4-1;
}


static void LeftDrag (double x, double y) 
{
    boundaries(&x, &y);
    setQposSpos((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
    drawCrosshair(x, y);
    drawAlignment((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
}


static void LeftDown (double x, double y) 
{ 
    if (!CrosshairON) return;
    
    boundaries(&x, &y);
    setQposSpos((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
    drawCrosshair(x, y);
    drawAlignment((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);

    graphRegister (LEFT_DRAG, LeftDrag);
}


/* KEYBOARD handles x, y coords in sequence units! (Mouse does it in screen units)
*/
static void keyboard (int key)
{
    int x, y;

    switch (key) {
    case UP_KEY:   x = oldx;    y = oldy -1;   break;
    case DOWN_KEY: x = oldx;	y = oldy+1;    break;
    case LEFT_KEY: x = oldx-1;	y = oldy;      break;
    case RIGHT_KEY: x = oldx+1; y = oldy;      break;
    case '>':	x = oldx+1 ;	y = oldy+1 ;   break ;
    case '<':	x = oldx-1 ;	y = oldy-1 ;   break ;
    case '}':	x = oldx+1 ;	y = oldy-1 ;   break ;
    case '{':	x = oldx-1 ;	y = oldy+1 ;   break ;
    default: return;
    }

    if (x < 0) x = 0;
    else if (x > qlen/resfac-1) x = qlen/resfac-1;
    if (y < 0) y = 0;
    else if (y > slen-1) y = slen-1;

    graphActivate(dotGraph);
    setQposSpos(x, y);
    drawCrosshair((float)x/zoom+LeftBorder, (float)y/zoom+TopBorder);
    drawAlignment(x, y);

    oldx = x;
    oldy = y;
}


/* Find an executable and return its complete pathname.
 */
int findCommand (char *command, char **retp)
{
#if !defined(NO_POPEN)
    
    static char retstr[1025] ;
    char *path, file[1024], retval ;
    int found=0;

    /* Don't use csh - fails if the path is not set in .cshrc * /
       FILE *pipe;
       static char *cp, csh[]="/bin/csh";
    if (access(csh, X_OK)) {
	messout("Could not find %s", csh);
	return 0;
    }
    if (!(pipe = (FILE *)popen(messprintf("%s -cf \"which %s\"", csh, command), "r"))) {
	return 0;
    }

    while (!feof(pipe))
	fgets(retval, 1024, pipe);
    retval[1023] = 0;
    pclose(pipe);

    if ((cp = strchr(retval, '\n'))) *cp = 0;
    if (retp) *retp = retval;

    / * Check if whatever "which" returned is an existing and executable file * /
    if (!access(retval, F_OK) && !access(retval, X_OK))
        return 1;
    else
	return 0;
    */
    
    path = getenv("PATH");
    
    path = strtok(path, ":");
    while (path) {
	    
	strcpy(file, path);
	strcat(file,"/");
	strcat(file, command);
	if (!access(file, F_OK) && !access(file, X_OK)) {
	    found = 1;
	    break;
	}

	path = strtok(0, ":");
    }
    
    if (found) {
	strcpy(retstr, file);
	retval = 1;
    }
    else {
	strcpy(retstr, "Can't find Dotter in path");
	retval = 0;
    }

    if (retp) *retp = retstr;
    return retval ;

#endif
}


static void callDotter(int dotterZoom, int xstart, int ystart, int xend, int yend)
{
#if !defined(NO_POPEN)
    FILE *pipe;
    static char *noname = "NoName", *qnam, *snam;
    MSP *msp;

    /* Need a name for correct number of arguments */
    qnam = (*qname ? qname : noname);
    snam = (*sname ? sname : noname);

    if (xstart < 1)  xstart = 1;
    if (xend > qlen) xend = qlen;
    if (ystart < 1)  ystart = 1;
    if (yend > slen) yend = slen;


    /* Open pipe to new dotterBinary */
    if (!dotterBinary) { 
	printf("Looking for Dotter ...\n");
	if (!findCommand("dotter", &(dotterBinary))) {
	    messout("Failed to zoom in - %s.  "
		    "($PATH=%s)", dotterBinary, getenv("PATH"));
	    dotterBinary = 0;
	    return;
	}
    }
    printf("Calling %s with region: %d,%d - %d,%d\n", 
	   dotterBinary, xstart, ystart, xend, yend);
    fflush(stdout);

    pipe = (FILE *)popen(messprintf("/bin/csh -cf \"%s -z %d -q %d -s %d -S %s %d %s %d %s %s\"", 
				    dotterBinary, 
				    dotterZoom, 
				    xstart-1+qoffset, 
				    ystart-1+soffset, 
				    qnam, 
				    xend-xstart+1, 
				    snam, 
				    yend-ystart+1,
				    dotterBinary,
				    (Xoptions ? Xoptions : "")), 
			 "w");

    fwrite(qseq+xstart-1, 1, xend-xstart+1, pipe);
    fwrite(sseq+ystart-1, 1, yend-ystart+1, pipe);

    /* Pass on features */
    for (msp = MSPlist; msp; msp = msp->next) {
	fprintf(pipe, "%d %d %s %d %d %s %d %d\n", 
		msp->score, 
		msp->id, 
		msp->frame,
		msp->qstart +MSPoffset,
		msp->qend   +MSPoffset,
		msp->sname, 
		msp->sstart,
		msp->send);
    }
    fprintf(pipe, "%c\n", EOF);
    fflush(pipe);
#endif
}


static void callDotterParams(void)
{
    static int dotterZoom=0, 
    xstart=0, 
    ystart=0, 
    xend=0, 
    yend=0;

    if (!graphPrompt ("Dotter parameters: zoom (compression) factor, "
		      "xstart, ystart, xend, yend", 
		      messprintf("%d %d %d %d %d", 
				 dotterZoom, xstart, ystart, xend, yend),
		      "iiiii"))

	return;

    freeint(&dotterZoom);
    freeint(&xstart);
    freeint(&ystart);
    freeint(&xend);
    freeint(&yend);
    callDotter(dotterZoom, xstart, ystart, xend, yend);
}


static void MiddleUp (double x, double y) 
{
    int  xstart, ystart,
         xend, yend, t;

    boundaries(&x, &y);

    if (oldrectx) {
	graphXorLine(rectx, recty, rectx, oldrecty);
	graphXorLine(rectx, recty, oldrectx, recty);
	graphXorLine(oldrectx, recty, oldrectx, oldrecty);
	graphXorLine(rectx, oldrecty, oldrectx, oldrecty);
    }

    if (x < rectx) {
	t = rectx;
	rectx = x;
	x = t;
    }

    if (y < recty) {
	t = recty;
	recty = y;
	y = t;
    }

    xstart = (rectx-LeftBorder)*zoom*resfac;
    ystart = (recty-TopBorder)*zoom;

    xend = (x-LeftBorder)*zoom*resfac;
    yend = (y-TopBorder)*zoom;

    if (xend-xstart < 10 || yend-ystart < 10) return;	/* Accidental clicks */

    callDotter(0, xstart, ystart, xend, yend);
}


static void MiddleDrag (double x, double y) 
{
    boundaries(&x, &y);

    if (oldrectx) {
	graphXorLine(rectx, recty, rectx, oldrecty);
	graphXorLine(rectx, recty, oldrectx, recty);
	graphXorLine(oldrectx, recty, oldrectx, oldrecty);
	graphXorLine(rectx, oldrecty, oldrectx, oldrecty);
    }
    
    graphXorLine(rectx, recty, rectx, y);
    graphXorLine(rectx, recty, x, recty);
    graphXorLine(x, recty, x, y);
    graphXorLine(rectx, y, x, y);

    oldrectx = x;
    oldrecty = y;
}

static void MiddleDown (double x, double y) 
{ 
    boundaries(&x, &y);

    rectx = x;
    recty = y;
    oldrectx = oldrecty = 0;

    graphRegister (MIDDLE_DRAG, MiddleDrag);
    graphRegister (MIDDLE_UP, MiddleUp);
}


static void initWindow(char *winsize)
{
    double 
	exp1, exp2, exp3;
    int
	win1, win2, win3;

    /* Call winsizeFromlambdak even if we don't want to set the window
       size in order to get the other parameters (exp_res_score) */

    if (blastx) {
	win1 = winsizeFromlambdak(MATRIX, atob_0, abetsize, pepqseq[0], sseq, &exp1, &Lambda);
	win2 = winsizeFromlambdak(MATRIX, atob_0, abetsize, pepqseq[1], sseq, &exp2, &Lambda);
	win3 = winsizeFromlambdak(MATRIX, atob_0, abetsize, pepqseq[2], sseq, &exp3, &Lambda);
	exp_res_score = (exp1 + exp2 + exp3)/3.0;
	win = (win1 + win2 + win3)/3.0;
    }
    else if (blastn)
	win = winsizeFromlambdak(MATRIX, ntob, abetsize, qseq, sseq, &exp_res_score, &Lambda);
    else
	win = winsizeFromlambdak(MATRIX, atob_0, abetsize, qseq, sseq, &exp_res_score, &Lambda);
    


    if (!winsize || ace_upper(*winsize) == 'K') {
	if (win < 3) {
	    messout("Karlin/Altschul estimate of window size = %d ignored. Using 10 instead.\n", win);
	    win = 10;
	}
	if (win > 50) {
	    messout("Karlin/Altschul estimate of window size = %d ignored. Using 50 instead.\n", win);
	    win = 50;
	}
	return;
    }

    if (!atoi(winsize))
	fatal("Bad window size specification: %s", winsize);
    win = atoi(winsize);
}


static void setWindow(void)
{
    if (!graphPrompt ("Give window size", 
		      messprintf("%d", win), 
		      "iz"))
	return;

    freeint(&win) ;
    calcWindow();
    pixels_real2col(data);
    dotterRedraw();
    
    /*    
       sprintf(banner, "%s (horizontal) vs %s (vertical).  Window = %d, Pixel values = %d x score/residue, Matrix = %s", 
       qname, sname, win, pixelFac, MATRIX_NAME);
       graphBoxDraw(bannerbox, BLACK, WHITE);
       */  
}


static void dotterDestroy(void)
{

    /* Free Dotter stuff */
    handleDestroy(handle);

    if (blastx) {
	for (i = 0; i < 3; i++) {
	    messfree(pepqseq[i]);
	}
    }


    /* Free stuff messalloc'ed in calling routine (usually blixem or dotterMain) */
    if (qseq == sseq) /*mhmp 10.12.98 */
      messfree(qseq);
    else
      {
	messfree(qseq);
	messfree(sseq);
      }
    qseq = sseq = 0 ;

    /* Don't free MSP's since that will screw blixem up !!! */

    /* Never ever free pixel arrays here - done by X !!!!! */

    if (graphActivate(alnGraph)) graphDestroy();
    if (graphActivate(fsGraph)) graphDestroy();
    /* if (graphActivate(rampGraph)) graphDestroy();*/
}


static void initCrosshair()
{
    if (!CrosshairON) {
	menuSetFlags(menuItem(dotterMenu, toggleCrosshairPosStr), MENUFLAG_HIDE);
	menuSetFlags(menuItem(dotterMenu, toggleCrosshairFullscreenStr), MENUFLAG_HIDE);
	return;
    }
    else {
	menuUnsetFlags(menuItem(dotterMenu, toggleCrosshairPosStr), MENUFLAG_HIDE);
	menuUnsetFlags(menuItem(dotterMenu, toggleCrosshairFullscreenStr), MENUFLAG_HIDE);
    }

    /* Set up crosshair boxes 
     ************************/
    graphColor(BLUE);
	
    if (!CrosshairFullscreenON) {
	hLineBox = graphBoxStart();
	graphLine(LeftBorder, crossy, LeftBorder-2+qlen4, crossy);
	graphBoxEnd();
	
	vLineBox = graphBoxStart();
	graphLine(crossx, TopBorder, crossx, TopBorder-2+slen4);
	graphBoxEnd();
    }
    else {
	int nx, ny;
	graphFitBounds (&nx, &ny);
	hLineBox = graphBoxStart();
	graphLine(0.0, crossy, (float)nx+1000, crossy);
	graphBoxEnd();
	
	vLineBox = graphBoxStart();
	graphLine(crossx, 0.0, crossx, (float)ny+1000);
	graphBoxEnd();
    }

    if (CrosshairPosON) {

	setQposSpos((int)(crossx-LeftBorder)*zoom, (int)(crossy-TopBorder)*zoom);
	sprintf(CrosshairPosText, "%s, %s", qpos, spos);
	CrosshairPosBox = graphBoxStart();
	graphTextPtr(CrosshairPosText, crossx+10, crossy+10, 25);
	graphBoxEnd();
	graphBoxDraw(CrosshairPosBox, BLUE, TRANSPARENT);
    }
	
    graphColor(BLACK);
}


static void togglePixelmap(void)
{
    if (!PixelmapON) 
    {
	if (!pixelmap_done) {
	    calcWindow();
	    pixels_real2col(data);
	}
	PixelmapON = 1;
    }
    else 
	PixelmapON = 0;

    dotterRedraw();
}


static void toggleGrid(void)
{
    gridOn = (gridOn ? 0 : 1);
    dotterRedraw();
}


static void initGreyramp(void)
{
/*
    float 
	minScore,		/ * score/residue * /
	maxScore;		/ * score/residue * /

    float
	offset;			/ * Disused idea for compensating for raised 
				   noise levels in compressed images * /
*/
    int 
	min = 40, 
	max = 100;

    
    graphGreyRamp(max, min);
    return;


    /* The pixelFac is set so that a pixel value of 50 = exp_res_score.
       Therefore, the stuff below is obsolete.  * /
    
    minScore = 0.8 * exp_res_score;
    maxScore = 2.0 * exp_res_score;

    offset = log((double)zoom*zoom)/(Lambda*win);

    printf("exp_res_score=%.2f => minScore=%.2f, maxScore=%.2f; offset=%.2f\n", 
	   exp_res_score, minScore, maxScore, offset);

    min = (maxScore  +offset)*(float)pixelFac;
    max = (minScore  +offset)*(float)pixelFac;

    if (min < 0) min = 0;
    if (min > 255) min = 255;

    if (max < 0) max = 0;
    if (max > 255) max = 255;

    graphGreyRamp(max, min);*/
}


static void dotterRedraw(void)
{
    static int pw, ph, bannerPos;
    static float fw, fh;
    int fsCount=0;

    if (!graphActivate(dotGraph))
    {
	float w, h;
	int fsCount=0;

	/* (rbrusk): should initialize these for all new dotGraphs...*/
	vLineBox = hLineBox = 0 ;

	graphScreenSize(&w, &h, &fw, &fh, &pw, &ph);
	
	if (fsArr) fsCount = arrayMax(fsArr);

	dotGraph = graphCreate (PIXEL_SCROLL, 
				messprintf("Dotter %s vs. %s", qname, sname), 
				0, 0, 
				(LeftBorder+qlen4+60 + 2*fsCount*ph/fh + (MSPlist ? 40:0))*w/pw, 
				(TopBorder +slen4+60 + 2*fsCount*ph/fh + (MSPlist ? 40:0))*h/ph);
	graphNewMenu (dotterMenu);
	graphRegister(MIDDLE_DOWN, MiddleDown) ;
	graphRegister(KEYBOARD, keyboard);
	graphRegister(DESTROY, dotterDestroy ); 
	graphRegister(LEFT_DOWN, LeftDown);

	initGreyramp();
	graphRawMaps (rawMap, rawMapRev);
	graphRampTool();

	graphActivate(dotGraph);
    }
    graphClear();
    graphPop();
    
/*    sprintf(banner, "%s (horizontal) vs. %s (vertical).  Window = %d, Pixel values = %d x score/residue, Matrix = %s",
	    qname, sname, win, pixelFac, MATRIX_NAME);
    graphTextPtr(banner, 5, 5, strlen(banner)+5);
*/
    
    sprintf(banner, "%s (horizontal) vs. %s (vertical)", qname, sname);
    bannerPos = (strlen(banner)*pw/fw < qlen ? LeftBorder : 10);
    graphText(banner, bannerPos, 5);
    if (!printMode) graphButton("About", Help, bannerPos, 20);
    
    if (fsArr) fsCount = arrayMax(fsArr);
    graphPixelBounds(LeftBorder+qlen4+60 + 2*fsCount*ph/fh + (MSPlist ? 40:0), 
		     TopBorder +slen4+60 + 2*fsCount*ph/fh + (MSPlist ? 40:0));

    if (PixelmapON && !(BlastHSPsOn && BlastMode == BLASTGREY)) {
	graphPixelsRaw ((char *)data, qlen4, slen4, qlen4, LeftBorder, TopBorder);
    }
	
    if (BlastHSPsOn) {

	static int gapwarned = 0;

	if (BlastMode == BLASTGREY) {
	    if (!HSPpixels) { 
		HSPpixels = (UCHAR *)malloc(datalen);
		for (i=0; i < datalen; i++) HSPpixels[i] = 0;
	    }
	    graphPixelsRaw ((char *)HSPpixels, qlen4, slen4, qlen4, LeftBorder, TopBorder);
	}

	drawBlastHSPs();

	if (HSPgaps && !gapwarned) {
	    graphRedraw();
	    messout("Note: gapped HSPs are shown ungapped in Dotter.");
	    gapwarned = 1;
	}
    }

    initCrosshair();
    drawScale();
    drawGenes(MSPlist);

    setMenuCheckmarks();
    graphRedraw();

    if (graphActivate(alnGraph)) 
	graphPop();
    else {
	if (!alignmentInitialized  && !getenv("ACEDB_PROJECT")) {
	    initAlignment();
	    alignmentInitialized = 1;
	}
    }

    if (graphActivate(rampGraph)) 
	graphPop();

    graphActivate(dotGraph);
}


Graph dotter (char  type,
	      char *opts,
	      char *queryname,
	      char *queryseq,
	      int   qoff,
	      char *subjectname,
	      char *subjectseq,
	      int   soff,
	      int   qcenter,
	      int   scenter,
	      char *savefile,
	      char *loadfile,
	      char *mtxfile,
	      char *featurefile,
	      float memoryLimit,
	      int   zoomFac,
	      MSP  *MSPs,
	      int   MSPoff,
	      char *winsize,
	      int   pixelFacset)
{
    int  i;

    /* Reset global statics */
    resfac = PixelmapON = 1;
    blastp = blastn = blastx = selfcomp =
	Display_selfcomp = watsonOnly = crickOnly = reversedScale = pixelmap_done = 0;
    BlastHSPsOn = 0;
    saveFil = 0;

    switch(type) {
    case 'P':  
	blastp = 1; 
	abetsize = 20;  break;
    case 'N':  
	blastn = 1; 
	abetsize = 4;   break;
    case 'X':  
	blastx = 1; 
	resfac = 3; 
	abetsize = 20;  break;
    default: fatal("Invalid sequence type passed to Dotter: %c", type);
    }
  
    /* Option parsing */
    if (opts)
      while (*opts) {
	switch (*opts) {
	case 'D': Display_selfcomp = 1; break;
	case 'R': reversedScale = 1;    break;
	case 'H': 
	    BlastHSPsOn = 1;         
	    PixelmapON = 0;             break;
	case 'W': watsonOnly = 1;       break;
	case 'C': crickOnly = 1;        break;
	case 'G': HSPgaps = 1;          break;
	}
	opts++;
      }
    plusmin = ( reversedScale ? -1 : 1 );
    
    if (graphActivate(dotGraph)) {
	dotterDestroy();
	/* Don't free data or HSPpixels, since every time rawImage 
	   is called with new data, a new XImage struct is added to gSubDev->images.
	   These are later free'd by XDestroyImage in graphSubDevDestroy.

	   However, memory allocated this way never seems to be entirely free'd
	   by graphDestroy, in the sense that it's not always reused. 
	   Maybe an illusion? */
	   
	/* free(data);
	   if (HSPpixels) free(HSPpixels);*/
    }
    HSPpixels = 0;
    handle = handleCreate();
    banner = (char *)handleAlloc(0, handle, 1000);

    strncpy(qname, queryname, NAMESIZE); qname[NAMESIZE] = 0;
    qseq = queryseq;
    qoffset = qoff;

    strncpy(sname, subjectname, NAMESIZE); sname[NAMESIZE] = 0;
    sseq = subjectseq;
    soffset = soff;
    MSPlist = MSPs;
    MSPoffset = MSPoff;
  
    if (!(qlen = strlen(qseq))) fatal("queryseq is empty");
    if (!(slen = strlen(sseq))) fatal("subjectseq is empty");

    for (i = 0; i < qlen; i++) qseq[i] = ace_upper(qseq[i]);
    for (i = 0; i < slen; i++) sseq[i] = ace_upper(sseq[i]);

    if (!memoryLimit) memoryLimit = 0.5; /* Mb */

    if (!strcmp(qseq, sseq)) selfcomp = 1;

    if (blastn && !watsonOnly) {
	/* Reverse and complement qseq for drawAlignment */
	qrevcompl = handleAlloc(0, handle, qlen+1);
	revcomp(qrevcompl, qseq);
	/*for (i = 0; i < qlen; i++) qrevcompl[i] = bton[ntob_compl[qseq[qlen-i-1]]];
	qrevcompl[qlen] = 0;*/
    }
    
    if (blastx) {
	for (i = 0; i < 3; i++)
	    pepqseq[i] = translate(qseq+i, stdcode1);
    }


    /* Get score matrix */
    if (mtxfile) {
	readmtx(MATRIX, mtxfile);
	strncpy(MATRIX_NAME, mtxfile, 80);
    }
    else {
	if (blastn) {
	    DNAmatrix(MATRIX);
	    strcpy(MATRIX_NAME, "DNA+5/-4");
	}
	else {
	    mtxcpy(MATRIX, BLOSUM62);
	    strcpy(MATRIX_NAME, "BLOSUM62");
	}
    }


    /* Don't do batch processing if output file can't be opened */
    if (savefile) 
    {
	if (!(saveFil = fopen (savefile, "w")))
	    fatal("Failed to open %s", savefile);
    }
	

    if (loadfile)
    {
	/* Load existing dotplot file */
	loadPlot(loadfile);
    }
    else 
    {
	initWindow(winsize);

	/* Set pixelFac so that exp_res_score is at 1/5 of the range. 
	 * This positions exp_res_score at 51.2
	 */
	if (pixelFacset) pixelFac = pixelFacset;
	else pixelFac = 0.2*256/exp_res_score;
	
	if (!zoomFac)
	    zoom = (int)sqrt((qlen/resfac/1e6*slen - 1e-6)/memoryLimit) + 1;
	else
	    zoom = zoomFac;

	qlen4 = (int)ceil((double)qlen/resfac/zoom);
	if (qlen4 % 4) qlen4 += 4-(qlen4 % 4);
	if (qlen/resfac > qlen4*zoom) fatal("qlen/resfac > qlen4*zoom (%d > %d (%d*%d))",
				     qlen/resfac, qlen4*zoom, qlen4, zoom);

	slen4 = (int)ceil((double)slen/zoom);
	if (slen4 % 4) slen4 += 4-(slen4 % 4);
	if (slen > slen4*zoom) fatal("slen > slen4*zoom (%d > %d (%d*%d))",
				     slen, slen4, zoom, slen4*zoom);

	datalen = slen4*qlen4;
	data = (UCHAR *)malloc(datalen);
    }

    if (savefile) {
	calcWindow();
	savePlot();
    }
    else
    {
	dotterMenu = menuInitialise ("dotter", (MENUSPEC*)mainMenu) ;
	if (!MSPlist)
	    menuSuppress (dotterMenu, drawBlastHSPgrayStr);
	menuSetFlags(menuItem(dotterMenu, selectFeaturesStr), MENUFLAG_HIDE);

	if (qcenter < 1 || qcenter > qlen) {
	    oldx = qlen/resfac/2;
	    crossx = oldx/zoom + LeftBorder;
	}
	else {
	    crossx = LeftBorder + qcenter/resfac/zoom;
	    oldx = qcenter;
	}

	if (scenter < 1 || scenter > slen) {
	    oldy = slen/2;
	    crossy = oldy/zoom + TopBorder;
	}
	else {
	    crossy = TopBorder + scenter/zoom;
	    oldy = scenter;
	}
	
	RightBorder = LeftBorder-1 + qlen4 - qlen % 4;
	
	if (!fsArr) fsArr = arrayCreate(50, FEATURESERIES);
	if (featurefile) {
	    FILE *file;
	    if (!(file = fopen(featurefile, "r"))) messcrash("Failed to open %s", featurefile);
	    loadFeatures(file);
	}

	dotterRedraw();

	if (!loadfile && !BlastHSPsOn) calcWindow();

	pixels_real2col(data);

    }
    return FALSE;
}


static void readmtx(int MATRIX[24][24], char *mtxfile)
{
    FILE *fil;
    int row, col;
    char line[1025] = "#", *p;
    
    if (!(fil = fopen(mtxfile, "r")) &&
	!(fil = fopen(messprintf("%s/%s", getenv("BLASTMAT"), mtxfile), "r")))
	fatal("Failed to open score matrix file %s - not found in ./ or $BLASTMAT/.", mtxfile);
    
    /* Ignore header ... */
    while (!feof(fil) && *line == '#') 
	fgets(line, 1024, fil);

    /* Read in the pairwise scores */
    for (row = 0; row < 24; row++)
    {
	if (!(fgets(line, 1024, fil)))
	    fatal("Wrong number of rows in matrix file: %d (should be 24).", row);

	p = strtok(line, " \t\n");
	for (col = 0; col < 24; col++) 
	{
	    while (*p == '*' || isalpha((int) *p))
		p = strtok(NULL, " \t\n");
	    if (!p) fatal("Error on row %d in matrix file.", row);

	    MATRIX[row][col] = atoi(p);

	    p = strtok(NULL, " \t\n");
	}
    }

    printf("I read your score matrix %s.\n", mtxfile);
    fclose(fil);
}

static void mtxcpy(int dest[24][24], int src[24][24])
{
    int i, j;

    for (i = 0 ; i < 24 ; i++)
	for (j = 0 ; j < 24 ; j++)
	    dest[i][j] = src[i][j];
}


static void DNAmatrix(int mtx[24][24])
{
    int i, j;

    for (i = 0 ; i < 6 ; i++)
	for (j = 0 ; j < 6 ; j++) {
	    if ( i < 4 && j < 4) 
		mtx[i][j] = (i == j ? 5 : -4);
	    else 
		mtx[i][j] = -4;
	}
}


/* Use this routine to add -install in programs that use Dotter */
void argvAdd(int *argc, char ***argv, char *s)
{
    char **v;
    int i;

    v = (char **)malloc(sizeof(char *)*(*argc+2));

    /* Copy argv */
    for (i = 0; i < (*argc); i++)
	v[i] = (*argv)[i];

    /* Add s */
    v[*argc] = (char *)malloc(strlen(s)+1);
    strcpy(v[*argc], s);

    (*argc)++;

    /* Terminator - ANSI C standard */
    v[*argc] = 0;

    /* free(*argv); */   /* Too risky */
    
    *argv = v;
}

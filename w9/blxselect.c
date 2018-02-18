/*  Last edited: Sep  4 18:01 1996 (esr) */

/* $Id: blxselect.c,v 1.2 2003/03/17 05:04:41 mieg Exp $ */

/* BLXSELECT - select seqbl/exblx files for blixem in a user-friendly way
 *
 * Erik Sonnhammer, 940329
 *

   Date   Modification
--------  ---------------------------------------------------
95-05-31  Got rid of popen reading of ls output (Mac unfriendly)
95-06-21  Got rid of fixed nr of boxes; changed to names array
          Added keyboard control.
96-08-29  [1.2] Added automatic dottering of first match.
96-09-01  [1.2] Added selection of queries.

*/

#include "regular.h"
#include "graph.h"
#include "blxview.h"
#include "../wh/menu.h"
#include "key.h"

#define MAXLINE  512
#define MAXLENGTH 1000
#define boxColor LIGHTGRAY


static char   qname[MAXLINE+1], 
              EXT[128]=".seqbl", 
              SEQEXT[128]=".seq", 
              list[128], **gargv,
              blixelectVersion[] = "1.2",
              zeroString[] = "with at least one match ";
static int    zeroOK=1, lastbox=0, items=0, itemsPerLine, itemsPerCol, *gargc,
              backgColor = LIGHTGRAY,
              dotter_first = 1;
static Graph  g;
static Array  names, selected;
static void   Help(void);
static void   doselect(void);
static void   unselect(void);
static void   printSelected(void);
static void   exchangeSelected(void);

static MENU Menu ;
static MENUOPT menu[] = {  
    graphDestroy,     "Quit" ,
    doselect,           "Select this query",
    unselect,         "Unselect this query",
    printSelected,    "Print selected queries",
    exchangeSelected, "Exchange selected set",
    Help,             "Help",
    graphPrint,       "Print",
    graphCleanUp,     "Clean up",
    0, 0 } ;


static void Help(void)
{
    if (zeroOK) *zeroString = 0; 

    graphMessage (messprintf("\
Blixelect - the Blixem file chooser.\n\
\n\
Version %s\n\
Copyright (C) Erik Sonnhammer, 1995\n\
\n\
\n\
Each sequence %sis listed with the number of matches to it.\n\
\n\
LEFT MOUSE BUTTON:\n\
  Double click on sequence to view in Blixem.\n\
\n\
MIDDLE MOUSE BUTTON:\n\
  Use to drag scrollbars.\n\
\n\
RIGHT MOUSE BUTTON:\n\
  Menu.\n\
\n\n\
KEYSTROKES:\n\
  Arrow keys: Go to next sequence to view in Blixem.", blixelectVersion, zeroString));
}


static void callBlixem(box)
{
    char *name, seqfilename[MAXLINE+1], HSPfilename[MAXLINE+1];
    FILE         *seqfile, *HSPfile;

    if (box > arrayMax(names)) return;

    box--;

    name = arr(names, box, char*); 

    strncpy(seqfilename, name, MAXLINE);
    strcat(seqfilename, SEQEXT);

    strncpy(HSPfilename, name, MAXLINE);
    strcat(HSPfilename, EXT);

    strncpy(qname, name, MAXLINE);
	    
    if (!(HSPfile = fopen(HSPfilename, "r"))) {
	messdump("Cannot open HSP file %s\n", HSPfilename); 
	return;
    }
    if (!(seqfile = fopen(seqfilename, "r"))) {
	messdump("Cannot open sequence file %s\n", seqfilename); 
	return;
    }
    
    blxreadhsp(seqfile, HSPfile, qname, 1, 0, 
	       (dotter_first ? " M     d" : " M     ")
	       , gargc, gargv);

    fclose(seqfile);
    fclose(HSPfile);
}


static void boxPick (int box)
{
    if (!box) return;

    if (box == lastbox)
	callBlixem(box);
    else
    {
	/* Turn last box off */
	if (lastbox) {
	    if (array(selected, lastbox, int))
		graphBoxDraw(lastbox, BLACK, RED);
	    else
		graphBoxDraw(lastbox, BLACK, boxColor);
	}

	/* Turn this box on */
	graphBoxDraw(box, WHITE, BLACK);
	lastbox = box;
    }
}


static void keypressed (int key)
{
    int box;

    if (!lastbox) return;

    switch (key) {
    case UP_KEY:    box = lastbox-1;    break;
    case DOWN_KEY:  box = lastbox+1;    break;
    case LEFT_KEY:  box = lastbox-itemsPerCol; break;
    case RIGHT_KEY: box = lastbox+itemsPerCol; break;
    default: return;
    }

    if (box < 1 || box > arrayMax(names)) return;

    /* Turn last box off */
    if (array(selected, lastbox, int))
	graphBoxDraw(lastbox, BLACK, RED);
    else
	graphBoxDraw(lastbox, BLACK, boxColor);

    /* Turn this box on */
    graphBoxDraw(box, WHITE, BLACK);
    lastbox = box;

    callBlixem(box);
}


int countHSP(char *filename)
{
    FILE *HSPfile;
    int n=0;
    char line[MAXLENGTH+1], realname[MAXLINE+1];

    strcpy(realname, filename);
    strcat(realname, EXT);

    if (!(HSPfile = fopen(realname, "r"))) {
	fprintf(stderr, "Cannot open HSP file %s\n", realname); 
	return 0;
    }
    
    while (!feof(HSPfile)) {
	if (!fgets(line, MAXLENGTH, HSPfile)) break;
	if (*line != '#') n++;
    }

    fclose(HSPfile);

    return n;
}


/* SEQBLMENU generates a menu of seqbl files */
static void seqblmenu(FILE *file)
{
    char   text[MAXLINE+1], *cp, *c;
    int    i, x, y, len, maxLen=0, n, box;
    Array  counts;
    float  nx, ny;

    names = arrayCreate(10, char*);
    counts = arrayCreate(10, int);

    while (!feof (file))
    { 
	if (!fgets (text, MAXLINE, file)) break;
	if ((c = (char *)strchr(text, '\n'))) *c = 0;
	len = strlen(text);
	n = countHSP(text);
	if (len && (n || zeroOK)) {
	    if (len > maxLen) maxLen = len;
	    array(names, items, char*) = messalloc(len+1);
	    strcpy(array(names, items, char*), text); 
	    array(counts, items, int) = n;
	    items++;
	}
    }
    fclose (file);

    g = graphCreate (TEXT_FULL_SCROLL, messprintf("Blixelect - the Blixem file chooser (File: %s, %d seqs)", list, items), 0, 0, 0.7, 0.5);
    graphRegister (PICK, boxPick);
    graphRegister (KEYBOARD, keypressed);
    Menu = menuInitialise ("blixelect", (MENUSPEC*)menu) ;
    graphNewMenu (Menu);

    /* Make square shaped chooser. 1.5 typical character height/width ratio */
    itemsPerCol  = (int)(sqrt((double)(items*(maxLen+6)))/1.5);
    nx = (maxLen+6)*ceil((float)items/itemsPerCol)+2;
    ny = itemsPerCol+2;
    graphTextBounds (nx, ny);

    graphClear();
    graphBoxDraw(0, backgColor, backgColor);
    graphColor(backgColor); graphRectangle(0, 0, nx+100, ny+100); graphColor(BLACK);

    x = y = 0;
    for (i = 0; i < items; i++)
    {
	if (y && !(y % itemsPerCol)) {
	    /* Go to next column */
	    x++;
	    y = 0;

	    /* Print column separator */
	    graphColor(DARKGRAY);
	    graphLine(x*(maxLen+6), 1.5, x*(maxLen+6), 1.5+itemsPerCol);
	    graphColor(BLACK);
	}
	
	box = graphBoxStart();
	graphText (messprintf("%3d %s", arr(counts, i, int), arr(names, i, char*)),
		   1+x*(maxLen+6), 1.5+y);
	graphBoxEnd(); 
	graphBoxDraw(box, BLACK, boxColor);
	y++;
    }

    graphButton("Help", Help, 0.5, 0.1);

    selected = arrayCreate(items+1, int);
    
    graphRedraw() ;
}


void main(int argc, char **argv)
{
    FILE	*file;
    int          optc;
    extern int   optind;
    extern char *optarg;
    char        *optstring="q:s:zD";

    static char *cc_date =
#if defined(__DATE__)
	__DATE__
#else
	    ""
#endif
		;

    char *usage;
    static char usageText[] = "\
\n\
 Blixelect - select files to view in Blixem.\n\
\n\
 Reference:  Sonnhammer ELL & Durbin R (1994). A workbench for Large Scale\n\
 Sequence Homology Analysis. Comput. Applic. Biosci. 10:301-307.\n\
\n\
 Usage: blixelect [options] <file_of_sequence-filenames>\n\
\n\
 Options:\n\
 -s <ext>   sequence filename extension (default .seq).\n\
 -q <ext>   seqbl filename extension (default .seqbl).\n\
 -z         Don't display sequences with 0 matches.\n\
 -D         Don't automatically dotter first sequence.\n\
\n\
 o Each sequence-filename must have a corresponding seqbl file.\n\
 o To make the seqbl files from blast output, run MSPcrunch with option -q.\n\
 o The file_of_sequence-filenames should not contain the .seq extensions.\n\
\n\
 by Erik Sonnhammer (esr@sanger.ac.uk)\n\
 Version ";

    usage = messalloc(strlen(usageText) + strlen(blixelectVersion) + strlen(cc_date) + 20);
    sprintf(usage, "%s%s, compiled %s\n", usageText, blixelectVersion, cc_date);
    
    while ((optc = getopt(argc, argv, optstring)) != -1)
	switch (optc) {
	case 'q': strcpy(EXT, optarg);          break;
	case 's': strcpy(SEQEXT, optarg);       break;
	case 'z': zeroOK = 0;                   break;
	case 'D': dotter_first = 0;             break;
	}

    if (argc - optind == 1) {
	if (!(file = fopen(argv[argc-1], "r")))
	    fatal("Cannot open exbl file %s\n", argv[argc-1]);
	strcpy(list, argv[argc-1]);
    }
    else {
	fprintf(stderr, "%s\n", usage); exit(1);
    }
    
    gargc = &argc;
    gargv = argv;

    /* Add -install for private colormaps */
    argvAdd(&argc, &argv, "-install");

    graphInit (&argc, argv) ;

    seqblmenu(file);
    
    graphLoop(FALSE);
    graphFinish ();
}


void doselect(void) {
    if (!lastbox) return;
    array(selected, lastbox, int) = 1;
    graphBoxDraw(lastbox, BLACK, RED);
}
void unselect(void) {
    if (!lastbox) return;
    array(selected, lastbox, int) = 0;
    graphBoxDraw(lastbox, BLACK, boxColor);
}
void printSelected(void) {
    int i;

    for (i = 1; i < arrayMax(selected); i++)
	if (arr(selected, i, int)) printf("%s\n", arr(names, i-1, char*));
}
void exchangeSelected(void) {
    int i;

    for (i = 1; i <= items; i++) {
	if (arr(selected, i, int)) {
	    arr(selected, i, int) = 0;
	    graphBoxDraw(i, BLACK, boxColor);
	}
	else {
	    arr(selected, i, int) = 1;
	    graphBoxDraw(i, BLACK, RED);
	}
    }
}

/*  File: graphhelp.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: %W% %G%
 * Description: graphical agent of HTML help package.
 * Exported functions:
 * HISTORY:
 * Last edited: May 12 13:50 1999 (edgrif)
 * * May 12 13:50 1999 (edgrif): Add new graphSetBlockMode call.
 * * Aug 18 15:54 1998 (fw): created to resolve 
                             ifdef NON-GRAPHICs in w1/help.c
 * Created: Tue Aug 18 15:54:15 1998 (fw)
 *-------------------------------------------------------------------
 */


#include "graph.h"		/* graphical functions */
#include "../wgd/gd.h"		/* for GIF handling */
#include "key.h"		/* for keyboard control */
#include "dict.h"
#include "keyset.h"		/* for removal - ACEDB specific */
#include "help_.h"		/* for helpPackage internal functions */
#include "../wh/menu.h"

/************************************************************/
static void helpMenuQuit (MENUITEM);
static void helpMenuPrint (MENUITEM);
static void helpMenuIndex (MENUITEM);
static void helpMenuOpenHelpDoc (MENUITEM);
static void helpMenuOpenInBrowser (MENUITEM m);

static void helpButtonQuit (void);
static void helpButtonGoBack (void);
static void helpButtonGoForward (void);
static void helpButtonIndex (void);
static void helpButtonPrint (void);
static void helpButtonOpenHelpDoc (void);

static void helpResize (void) ;
static void helpDraw (void) ;
static void helpDestroy (void) ;

static void helpPick (int box);
static void helpOpenInBrowser (int box);
/************************************************************/




/************************************************************/
/******************                   ***********************/
/******************  public functions ***********************/
/******************                   ***********************/
/************************************************************/


/************************************************************/
GRAPH_FUNC_DEF void help (void)
{
  /* display help on subject that is registered 
     with currently active window */
  /* NOTE: This function has the correct prototype for a MENUOPT */
  char *subject = graphHelp(0);

  if (!subject)
    {
      messout ("Sorry, no help subject registered for this window !") ;
      return ;
    }

  helpOn (subject);

  return;
} /* help */
/************************************************************/



/************************************************************/
/******************                   ***********************/
/****************** private functions ***********************/
/******************                   ***********************/
/************************************************************/



/************************************************************/
static void helpMenuSpacer (MENUITEM m) {;} /* dummy */


static MENUSPEC helpControlMenu [] = {
  { helpMenuQuit,		"Quit"},
  { helpMenuPrint,		"Print"},
  { helpMenuSpacer,		""},
  { helpMenuIndex,		"Subject Index"},
  { helpMenuOpenHelpDoc,	"Open help document"},
  { helpMenuOpenInBrowser,	"Open current in browser"},
  { 0, 0}
};
static MENU controlMenu;	/* constructed from helpControlMenu[] */


static MENUOPT linkMenu[] = {
  { helpButtonGoBack,		"Back"},
  { helpButtonGoForward,	"Forward"},
  { menuSpacer,			""},
  { (GraphFunc)helpPick,		"Open link"},
  { (GraphFunc)helpOpenInBrowser,	"Open in browser"},
  { 0, 0}
};


static Graph     helpGraph = 0,
                 helpSubGraph = 0 ;
static Array     box2href = 0,
                 pageArray = 0,
                 pageList = 0;
static HtmlPage *currPage ;
static char      currHelpFilename[MAXPATHLEN];

static Array imageArray ;
static DICT *imageDict ;
static DICT *nameRefDict = 0 ;
static KEYSET nameRefLine = 0 ;
static int WINX ;
static float xPos ;
static int indent ;

static DICT *pageDict = 0 ;
static float yPos ;
static float chWidth, chHeight, lineHeight ;
static int WINMAXX=0, currPageListIndex = -1, currPageArrayIndex = -1 ;
static int goForwardButton, goBackButton;
/************************************************************/

GRAPH_FUNC_DEF BOOL graphHtmlViewer (const char *helpFilename)
{
  Graph oldGraph ;
  int pageDictIndex ;

  if (!helpFilename)
    {
      messerror ("Sorry, no help available !");
      return FALSE;
    }

  if (currPageListIndex == -1)
    {
      /**** init everything for new window ****/

      /* keep a list of pointers for each page */
      pageArray = arrayCreate (10, HtmlPage*) ;

      /* dictionary of pages, indexed by their filename */
      pageDict = dictCreate (10) ;

      /* keeps a list of indexes to pages already visited */
      pageList = arrayCreate (10, int) ;

      imageDict = dictCreate (2) ;
      imageArray = arrayCreate (2, gdImagePtr) ;
    }

  /* make a copy, because helpFilename, is a static within
     getSubjectHelpfilename, and may change be changed without
     out control here */
  strncpy (currHelpFilename, helpFilename, MAXPATHLEN);

  /* he we got this file already loaded ? */
  if (dictAdd(pageDict, helpFilename, &pageDictIndex))
    {
      /* filename is new */
      currPage = htmlPageCreate(helpFilename);

      array (pageArray, pageDictIndex, HtmlPage*) = currPage;
    }
  else
    {
      /* we've had this filename before */
      currPage = arr(pageArray, pageDictIndex, HtmlPage*);

      if (!currPage)		/* we couldn't get it before */
	{
	   /* try again - to give user a chance to correct
	      the filesystem (move directories, etc..)
	      i.e. perform an action such that we'll find it this time */
	  currPage = htmlPageCreate(helpFilename);
	}
    }
  currPageArrayIndex = pageDictIndex ;


  currPageListIndex += 1;
      
  /* break the list, if we are taking a different route
     through the pagelist */
  if (currPageListIndex < arrayMax(pageList)-1 &&
      pageDictIndex != arr(pageList, currPageListIndex, int))
    {
      arrayMax(pageList) = currPageListIndex+1 ;
    }

  array (pageList, currPageListIndex, int) = pageDictIndex ;


  /*********** open window or redraw with new currPage ************/

  if ((oldGraph = graphActivate (helpGraph)))
    {
      /*** we already have a helpGraph ****/
      graphPop () ;		/* bring it to the forground */
    }
  else
    {
      WINX = 80*8 ;

      if (!(helpGraph = graphCreate
	    (TEXT_FIT, "On-Line Help", 0, 0,
	     (WINX+5*8)/900.0, 524/900.0)))
	return FALSE;
      graphRegister(DESTROY, helpDestroy) ;
      graphRegister(RESIZE, helpResize) ;

      /* the help-graph will still stay active, when a graph blocks
	 others (like the filechooser) */
      graphSetBlockMode(GRAPH_NON_BLOCKABLE) ;

      controlMenu = menuInitialise("Help menu", helpControlMenu);
    }

  helpResize () ;		/* force redraw */
	
  graphActivate (oldGraph) ;

  return TRUE;
} /* graphHtmlViewer */
/************************************************************/

static void helpMenuQuit (MENUITEM m) { graphDestroy(); }
static void helpButtonQuit (void) { graphDestroy(); }

/************************************************************/

static void helpButtonGoBack (void)
/* can't be called, if no previous level exists
   anyway, because the button doesn't appear then */
{
  int idx ;

  if (currPageListIndex > 0)
    {
      idx = arr (pageList, --currPageListIndex, int) ;

      currPage = arr (pageArray, idx, HtmlPage*) ;
      
      helpDraw () ;
    }
  else
    graphBoxDraw (goBackButton, GRAY, WHITE);

  return ;
} /* helpGoBack */
/************************************************************/

static void helpButtonGoForward (void)
/* can't be called, if no further level exists
   anyway, because the button doesn't appear then */
{
  int idx ;

  if (currPageListIndex < arrayMax(pageList)-1)
    {
      idx = arr (pageList, ++currPageListIndex, int) ;
      
      currPage = arr (pageArray, idx, HtmlPage*) ;
      
      helpDraw () ;
    }
  else
    graphBoxDraw (goForwardButton, GRAY, WHITE);

  return ;
} /* helpGoForward */
/************************************************************/

static void helpMenuIndex (MENUITEM m) { helpButtonIndex() ; }
static void helpButtonIndex (void)
{
  graphHtmlViewer ("?") ;
  /* the '?' character as a helpFilename is a special subject/filename
     to the helpPackage - to display a dynamically created marked up 
     directory listing of the help-directory */

  return;
} /* helpMenuIndex */
/************************************************************/

static void helpMenuPrint (MENUITEM m) { helpButtonPrint() ; }
static void helpButtonPrint (void)
{
  graphActivateChild() ;  /* because page is displayed in a subgraph */
  graphPrint () ;

  return;
} /* helpMenuPrint */
/************************************************************/

static void helpMenuOpenHelpDoc (MENUITEM m) { helpButtonOpenHelpDoc(); }
static void helpButtonOpenHelpDoc (void)
{
  static char fname[FIL_BUFFER_SIZE], dname[DIR_BUFFER_SIZE];
  static char helpFilename[DIR_BUFFER_SIZE+FIL_BUFFER_SIZE];
  FILE *fil;

  strcpy (dname, helpGetDir());
  strcpy (fname, "");

  if (!(fil = filqueryopen (dname, fname, HELP_FILE_EXTENSION, "r",
			    "Select help document")))
    return ;
  filclose (fil);

  helpSetDir (dname);

  sprintf(helpFilename, "%s%s%s",
	  dname, SUBDIR_DELIMITER_STR, fname);

  graphHtmlViewer (helpFilename);

  return;
} /* helpMenuOpenHelpDoc */

/************************************************************/

static void helpMenuOpenInBrowser (MENUITEM m)
{
  /* open current file in browser */
  graphWebBrowser (currHelpFilename);
}


static void helpOpenInBrowser (int box)
{
  /* open the link (that the mouse is hovering over) 
     in an external browser */
  char *link;
  if (!box) return;

  if (!(link = arr(box2href, box, char*))) return;

  graphWebBrowser (link) ;
 
  return;
} /* helpOpenInBrowser */

/************************************************************/

static void helpPick (int box)
/* references the link belonging to the box and displays it */
{
  char *link ; 

  if (!box) return ;
  if (!(link = arr(box2href, box, char*))) return;

  if (link[0] == '#')
    /* jump to specific line in help file */
    {
      int kk = 0 ;
      if (nameRefDict && dictFind (nameRefDict, link+1, &kk) &&
	  kk >= 0 && keySetExists(nameRefLine) &&
	  kk < keySetMax(nameRefLine))
	graphGoto (1,  keySet(nameRefLine, kk)) ;
    }
  else if (strncmp(link, "http:", 5) == 0 || 
	   strncmp(link, "ftp:", 4) == 0)
    /* open web document in external browser */
    {
      graphWebBrowser (link);
    }
  else if ((strcasecmp (filGetExtension(link), "gif") == 0) ||
	   (strcasecmp (filGetExtension(link), HELP_FILE_EXTENSION) == 0))
    {
      if (! graphHtmlViewer (helpLinkGetFilename (link)) )
	messerror("Can't open link 'HREF=%s' (%s)",
		  link, messSysErrorText());
    }
  else
    messerror ("Can't display HREF=%s in help browser !", link);

  return;
} /* helpPick */
/************************************************************/

static void helpMiddleDown (double x, double y)
/* opens the link of the selection in a browser
   (just like the middle button over a link in Netscape) */
{
  char *link;
  int box=0;
  float xmin, ymin, xmax, ymax;

  /* find out which box we clicked on */
  for (box = 1; box < arrayMax(box2href); ++box)
    {
      /* try every box, and look if its area covers the mouseclick */
      graphBoxDim(box, &xmin, &ymin, &xmax, &ymax);
      if (xmin <= x && x <= xmax &&
	  ymin <= y && y <= ymax)
	{
	  break;
	}
    }
  if (!box) return ;

  if (!(link = arr(box2href, box, char*))) return;

  if ((strncmp(link, "http:", 5) == 0 || 
       strncmp(link, "ftp:", 4) == 0))
    {
      graphWebBrowser (link) ;
    }
  else if (link[0] == '#')
    {
      graphWebBrowser (messprintf ("%s%s", currHelpFilename, link)) ;
    }
  else
    {
      /* the browser will handle errors/warnings if the link is bad */
      graphWebBrowser (helpLinkGetFilename (link));
    }

  return;
} /* helpMiddleDown */
/************************************************************/

static void helpKeyboard (KEY key)
{
  float x1, y1, x2, y2 ;
  graphActivate (helpSubGraph) ;
  graphWhere (&x1, &y1, &x2, &y2) ;

  switch (key)
    {
    case UP_KEY:
      graphGoto ((x2-x1), (y2-y1)/2 + y1 - chHeight) ;
      break ;
    case DOWN_KEY:
      graphGoto ((x2-x1), (y2-y1)/2 + y1 + chHeight) ;
      break ;
    case ' ':
      graphGoto ((x2-x1), (y2-y1)/2 + y1 + (y2-y1) - chHeight*2) ;
      break ;
    }

  return;
} /* helpKeyBoard */
/************************************************************/

static void helpResize (void)
{
#define nint(x) (((x) > 0) ? ((int)(x + .5)) : -((int)(.5 - x)) )
  float x1, y1, x2, y2 ;

  graphActivate (helpGraph) ;
  graphWhere (&x1, &y1, &x2, &y2) ;
  WINX = nint((x2 - 5)*8) ;      /* leave two char space at left
				    and right border */
  WINMAXX = 0 ;

  if (graphActivate (helpSubGraph))
    {
      graphDestroy () ;
    }
  helpSubGraph = graphSubCreate (PIXEL_SCROLL, 0, 3,
				 x2-(4/8.0), y2-3-(3/13.0)) ;

  graphRegister(PICK, helpPick) ;
  graphRegister(MIDDLE_DOWN, helpMiddleDown);
  graphRegister(KEYBOARD, helpKeyboard) ;

  graphSetBlockMode(GRAPH_NON_BLOCKABLE) ;

  if (WINX > chWidth)
    helpDraw () ;

  return ;
} /* helpResize */
/************************************************************/

static void newLine (void)
{
  if (xPos > WINMAXX)
    WINMAXX = xPos ;

  if (xPos != indent)
    {
      xPos = indent ;
      yPos += lineHeight ;
      lineHeight = chHeight+2 ;
    }

  return;
} /* newLine */
/************************************************************/

static void blankLine (void)
{
  newLine () ;
  yPos += lineHeight ;

  return;
} /* newLine */
/************************************************************/

static void printSection (HtmlNode *node)
{
  int i, len, dx, dy ;
  char *cp, *start ;
  float hrefXPos = -1 ;
  static BOOL 
    MODE_PREFORMAT=FALSE, 
    MODE_HREF=FALSE, 
    /*     MODE_HEADER=FALSE,  */
    FOUND_NOBULLET_IN_LIST_NOINDENT=FALSE ;
  static int box = -1, itemNumber ;
  static char *currentLink ;
  static char buf[10000];	/* for word wrapping ops */

  graphTextInfo (&dx, &dy, &chWidth, &chHeight) ;
  if ((chHeight+2) > lineHeight)
    lineHeight = chHeight+2 ;

  switch (node->type)
    {
    case HTML_SECTION:
      printSection (node->left) ;
      if (node->right) printSection (node->right) ;
      break ;
    case HTML_COMMENT:
      /* do nothing */
      break ;
    case HTML_DOC:
    case HTML_HEAD:
    case HTML_BODY:
      if (node->left) printSection (node->left) ;
      break ;
    case HTML_TITLE:
      graphRetitle (node->text) ;
      break ;
    case HTML_HEADER:
      {
	/* vars local to this bit here to avoid confusion */
	int oldFont=0, oldFontHeight=0, oldColor ;

	oldColor = graphColor (RED) ;
	if (node->hlevel == 1)
	  oldFont = graphTextFormat (BOLD) ;
	if (node->hlevel < 3)
	  oldFontHeight = graphTextHeight (20) ;
	
	/* MODE_HEADER = TRUE ; */

	blankLine () ;

	indent = node->hlevel*2 * chWidth ;
	xPos = indent ;
	
	/* check, in case some bozo has done a thing like <H1></H1> */
	if (node->left) printSection (node->left) ; 

	blankLine () ;

	graphColor (oldColor) ;
	if (node->hlevel == 1)
	  graphTextFormat (oldFont) ;
	if (node->hlevel < 3)
	  graphTextHeight (oldFontHeight) ;
	
	/* MODE_HEADER = FALSE ; */
      }
      break ;

    case HTML_LIST:
      newLine () ;
      if (node->lstyle == HTML_LIST_BULLET || 
	  node->lstyle == HTML_LIST_NUMBER)
	indent += 4*chWidth ;

      itemNumber = 0 ;

      /* a list might not have a leftnode (a list item) */
      if (node->left) printSection (node->left) ;

      if (node->lstyle == HTML_LIST_BULLET ||
	  node->lstyle == HTML_LIST_NUMBER)
	indent -= 4*chWidth ;
      if (node->lstyle == HTML_LIST_NOINDENT &&
	  FOUND_NOBULLET_IN_LIST_NOINDENT)
	{
	  indent -= 4*chWidth ;
	  FOUND_NOBULLET_IN_LIST_NOINDENT = FALSE ;
	}

      blankLine () ;
      break ;

    case HTML_LISTITEM:
      ++itemNumber ;
      if (node->left)
	{
	  newLine () ;
	  if (node->lstyle == HTML_LIST_BULLET ||
	      node->lstyle == HTML_LIST_NOINDENT)
	    {
	      /* draw a round bullet */
	      graphFillArc (xPos - 1.5*chWidth, yPos + 0.5*chHeight,
			    0.5*chWidth, 0, 360) ;
	    }
	  else if (node->lstyle == HTML_LIST_NUMBER)
	    {
	      /* draw the number of the listitem */
	      graphText (messprintf ("%d.", itemNumber),
			 xPos - 3*chWidth, yPos) ;
	      indent += (strlen(messprintf ("%d", itemNumber))-1)*chWidth ;
	      xPos = indent ;
	    }
	  else if (node->lstyle == HTML_LIST_NOBULLET)
	    {
	      /* part of a <DL> noindented list, but a <DD>
		 item becomes indented, but no bullet */
	      /* if we come across the first NO_BULLET item, in
		 a HTML_LIST_NOINDENT, the LIST becomes indented */
	      if (!FOUND_NOBULLET_IN_LIST_NOINDENT)
		{
		  indent += 4*chWidth ;
		  xPos = indent ;
		  FOUND_NOBULLET_IN_LIST_NOINDENT = TRUE ;
		}
	    }
	  else if (node->lstyle == HTML_LIST_NOINDENT_NOBULLET)
	    {
	      /* if we are in a <DL> list and went to indentation
		 because of a <DD> item, a <DT> item brings back
		 the old indent-level (noindent for <DL>'s) */
	      if (FOUND_NOBULLET_IN_LIST_NOINDENT)
		{
		  indent -= 4*chWidth ;
		  xPos = indent ;
		  FOUND_NOBULLET_IN_LIST_NOINDENT = FALSE ;
		}
	    }
	  printSection (node->left) ;
	}
      if (node->lstyle == HTML_LIST_NUMBER)
	{
	  indent -= (strlen(messprintf ("%d", itemNumber))-1)*chWidth ;
	}
      else if (node->lstyle == HTML_LIST_NOBULLET)
	{
	  if (!FOUND_NOBULLET_IN_LIST_NOINDENT)
	    indent -= 4*chWidth ;
	}
      
      if (node->right)
	{
	  printSection (node->right) ;
	}
      break ;

    case HTML_HREF:
      if (node->link && node->isNameRef)
	{
	  int kk, line = yPos ;
	  if (!nameRefDict)
	    { nameRefDict = dictCreate(20) ; nameRefLine = keySetCreate() ; }
	  dictAdd (nameRefDict, node->link, &kk) ;
	  keySet(nameRefLine, kk) = line ;
	}

      if (node->link && !node->isNameRef)
	{
	  MODE_HREF = TRUE ;
	  currentLink = node->link ;
	}
      /* we have to check for leftnode, in case we have a thing
	 like <A HREF=...></A>. The HREF-node doesn't have a TEXT
	 node attached, and it would crash otherwise */
      if (node->left) printSection (node->left) ; 

      if (node->link && !node->isNameRef)
	{
	  MODE_HREF = FALSE ;
	  currentLink = 0 ;
	}
      break ;
    case HTML_TEXT:
      if (MODE_HREF)
	{
	  /* xPos now stands at the end of the last word,
	     remember that position so we know from where to start
	     the underlining from */
	  hrefXPos = xPos ;
	  
	  /* in TEXT-mode, and when we're not at the beginning of
	     a line, we add 1 to account for the spaces between words */
	  if (!MODE_PREFORMAT && xPos != indent)
	    hrefXPos += chWidth ;

	  /* hrefCp = node->text ; */

	  box = graphBoxStart () ; /* draw the HREF box */
	  array (box2href, box, char*) = currentLink ;
	}

      cp = node->text ;
      if (!MODE_PREFORMAT)
	stripSpaces (node->text) ;
      /* for MODE_PREFORMAT keeps all controls chars */

      while (*cp)
	{
	  len = 0 ;
	  start = cp ;
	  
	  if (!MODE_PREFORMAT)
	    {
	      while (*cp && !isspace((int)*cp)) { ++(cp) ; ++len ; }
	      if (*cp) ++cp ;	/* skip whitespace */
	    }
	  else
	    {
	      while (*cp && *cp != '\n') { ++(cp) ; ++len ; }
	      if (*cp) 
		{
		  ++cp ;	/* skip RETURN */
		  ++len ;	/* so we copy the RETURN into buf */
		}
	    }

	  memset (buf, 0, 10000) ;
	  strncpy (buf, start, len) ;
	  buf[len] = 0 ;

	  /* linewrapping of words/lines longer than WINX */
	  if (strlen(buf) > WINX)
	    {
	      cp = start + (int)(WINX/chWidth) ;
	      buf[(int)(WINX/chWidth)] = 0 ;
	      len = (int)(WINX/chWidth) ;
	    }
	  
	  /* word wrapping if not in preformatting mode */
	  if (!MODE_PREFORMAT)
	    {
	      if (xPos != indent)	/* not at start of line ... */
		xPos += chWidth ;	/* ... one space before the word */
	      if (xPos + len*chWidth > WINX)
		{
		  if (MODE_HREF)
		    {
		      graphLine (hrefXPos, yPos+chHeight+1, xPos, yPos+chHeight+1) ;
		      hrefXPos = indent ; /* hrefCp = cp ; */
		      
		      /* word wrap the HREF-box */
		      graphBoxEnd () ;
		      graphBoxDraw (box, BLUE, TRANSPARENT) ;
		      graphBoxMenu (box, linkMenu);

		      box = graphBoxStart () ;
		      array (box2href, box, char*) = currentLink ;
		    }
		  newLine () ;
		}
	      graphText (buf, xPos, yPos) ;
	      xPos += strlen(buf)*chWidth ; /* place xPos at the end of word */
	    }
	  else if (MODE_PREFORMAT)
	    {
	      int oldpos, stringpos, screenpos, ii ;
	      i = 0 ;

	      /* replace TABs with appropriate number of spaces */
	      while (buf[i])
		{
		  if (buf[i] == '\t')
		    {
		      /* oldpos is the position, that this TAB char
			 would go on the screen without TABifying
			 NOTE: xPos is always at least "indent" 
			 (to leave a left margin) */
		      oldpos = (xPos - indent)/chWidth + i ;

		       /* screenpos is the position of the TAB char 
			  after inserting spaces
			  NOTE: the TAB itself will be
			  overwritten by one space */
		      screenpos = (((oldpos/8)+1)*8) - 1 ;

		       /* stringpos is where the TAB should go
			  in the string, where it'll turn into a space
			  at that position */
		      stringpos = screenpos - (xPos-indent)/chWidth ;

		      /* shift all text from current position "i"
			 onwards */
		      for (ii = strlen(buf)-1; ii >= i ; --ii)
			buf[ii+(stringpos-i)] = buf[ii] ;

		      /* fill gap with spaces and also overwrite TAB
			 with a space */
		      for (ii = i; ii <= stringpos; ++ii)
			buf[ii] = ' ' ;

		      i = stringpos ;
		    }
		  ++i ;
		}

      /* don't use len, it might have changed when inserting spaces */
	      if (buf[strlen(buf)-1] == '\n')
		{
		  buf[strlen(buf)-1] = 0 ;
		  graphText (buf, xPos, yPos) ;
		  xPos += strlen(buf)*chWidth ;
		  newLine ();	/* for the '\n' */
		}
	      else
		{
		  graphText (buf, xPos, yPos) ;
		  xPos += strlen(buf)*chWidth ;
		}
	    }
	}
      if (MODE_HREF)
	{
	  /* underline HREF's */
	  graphLine (hrefXPos, yPos+chHeight+1, xPos, yPos+chHeight+1) ;

	  graphBoxEnd () ;
	  graphBoxDraw (box, BLUE, TRANSPARENT) ;
	  graphBoxMenu (box, linkMenu);
	}
      break ;


    case HTML_GIFIMAGE:
      {
	int imbox ;
	int imageDictIndex;
	unsigned char red[128], green[128], blue[128], myPixels[256] ;
	int i, j, ncol, x, y ;
	unsigned char *pixels, *p ;
	gdImagePtr theImage ;
	FILE *fil;

	/* do we need to load the image ? */
	if (dictAdd (imageDict, node->link, &imageDictIndex))
	  {
	    char *filename = helpLinkGetFilename (node->link);

	    theImage = 0;

	    /* load the image */
	    fil = fopen (filename, "r");

	    if (fil)
	      theImage = gdImageCreateFromGif (fil);

	    if (!theImage)
	      {
		printf ("Warning : Image %s couldn't be opened !!\n",
			filename) ;
		
		node->type = HTML_NOIMAGE;
		printSection (node);
		messfree (filename);
		return;
	      }

	    array (imageArray, imageDictIndex, gdImagePtr) = theImage ;
	    
	    fclose (fil) ;

	    if ((ncol = gdImageColorsTotal(theImage)) > 128)
	      messout ("Warning : Image has more colors (%d) "
		       "than I can manage (128)", ncol) ;
	  }
	else
	  {
	    theImage = arr(imageArray, imageDictIndex, gdImagePtr) ;
	  }

	xPos += chWidth ;

	ncol = gdImageColorsTotal(theImage);
	if (ncol > 128)
	  ncol = 128 ;

	if (graphWritableColors())
	  {
	    for (i = 0 ; i < ncol ; ++i)
	      { red[i] = gdImageRed(theImage, i) ;
		green[i] = gdImageGreen(theImage, i) ;
		blue[i] = gdImageBlue(theImage, i) ;
	      }
	    graphSetColorMaps (red, green, blue) ;
	  }
	graphRawMaps (myPixels, 0) ;
	for (i = 0 ; i < 128 ; ++i)
	  myPixels[i] = myPixels[2*i] ; /* undo pixel-sharing */
	
	x = gdImageSX(theImage) ; while (x % 8) ++x ;
	y = gdImageSY(theImage) ;
	pixels = (unsigned char*) malloc (x * y * sizeof(unsigned char)) ;

	for (j = 0 ; j < y ; ++j)
	  for (p = pixels + j*x, i = 0 ; i < gdImageSX(theImage) ; ++i)
	    *p++ = (unsigned char) myPixels[(int)theImage->pixels[i][j]] ;
	
	imbox = graphBoxStart () ;
	graphPixelsRaw ((char*)pixels, gdImageSX(theImage), y,
			x, xPos, yPos) ;
	graphBoxEnd () ;
	graphBoxDraw (imbox, -1, -1) ;

	if (MODE_HREF)
	  {
	    int oldColor ;

	    oldColor = graphColor (BLUE) ;
	    graphRectangle (xPos-1, yPos-1,
			    xPos+gdImageSX(theImage), yPos+y) ;
	    array (box2href, imbox, char*) = currentLink ;
	    graphColor (oldColor) ;
	  }

	xPos += gdImageSX(theImage) ;
	if (gdImageSY(theImage)+2 > lineHeight)
	  lineHeight = gdImageSY(theImage)+2 ;
      }
      break ;
    case HTML_NOIMAGE:
      {
	int imbox ;

	xPos += chWidth ;

	imbox = graphBoxStart () ;
	graphRectangle (xPos, yPos, xPos+20, yPos+20) ;
	graphText ("?", xPos+6, yPos+3) ;
	graphBoxEnd () ;
	graphBoxDraw (imbox, -1, -1) ;

	if (MODE_HREF)
	  {
	    int oldColor ;

	    oldColor = graphColor (BLUE) ;
	    graphRectangle (xPos-1, yPos-1,
			    xPos+21, yPos+21) ;
	    array (box2href, imbox, char*) = currentLink ;
	    graphColor (oldColor) ;
	  }

	xPos += 20 ;
	if (20+2 > lineHeight)
	  lineHeight = 20+2 ;
      }
      break ;
    case HTML_RULER:
      {
	int oldColor ;

	newLine () ;
	yPos += 0.5*lineHeight ;
	oldColor = graphColor (GRAY) ;
	graphLine (indent + 1, yPos - 1,
		   (WINX-chWidth) - 1, yPos - 1) ;
	graphColor (BLACK) ;
	graphLine (indent, yPos,
		   (WINX-chWidth), yPos) ;
	graphColor (GRAY) ;
	graphLine (indent + 1, yPos + 1,
		   (WINX-chWidth) - 1, yPos + 1) ;
	graphColor (BLACK) ;
	graphColor (oldColor) ;
	yPos += 0.5*lineHeight ;
      }
      break ;
    case HTML_PARAGRAPH:
      blankLine () ;
      break ;
    case HTML_LINEBREAK:
      newLine () ;
      break ;
    case HTML_BOLD_STYLE:
    case HTML_STRONG_STYLE:
      {
	/* vars local to this bit here to avoid confusion */
	int oldFont ;

	oldFont = graphTextFormat (BOLD) ;
	if (node->left)	printSection (node->left) ;
	graphTextFormat (oldFont) ;
      }
      break ;
    case HTML_ITALIC_STYLE:
      {
	/* vars local to this bit here to avoid confusion */
	int oldFont ;

	oldFont = graphTextFormat (ITALIC) ;
	if (node->left) printSection (node->left) ;
	graphTextFormat (oldFont) ;
      }
      break ;
    case HTML_CODE_STYLE:
      {
	/* vars local to this bit here to avoid confusion */
	int oldFont ;

	oldFont = graphTextFormat (PLAIN_FORMAT) ;
	if (node->left)	printSection (node->left) ;
	graphTextFormat (oldFont) ;
      }
      break ;
    case HTML_STARTBLOCKQUOTE:
      newLine () ;
      indent += 3*chWidth ; xPos = indent ;
      break ;
    case HTML_ENDBLOCKQUOTE:
      indent -= 3*chWidth ;
      blankLine () ;
      break ;
    case HTML_STARTPREFORMAT:
      MODE_PREFORMAT = TRUE ;
      newLine () ;
      break ;
    case HTML_ENDPREFORMAT:
      MODE_PREFORMAT = FALSE ;
      break ;
    case HTML_UNKNOWN:
      break;			/* compiler happiness */
    }

  return;
} /* printSection */
/************************************************************/

static void helpDraw (void)
{
  graphActivate (helpGraph) ;
  graphClear () ;
  graphNewMenu (controlMenu);
  
  dictDestroy (nameRefDict) ;
  keySetDestroy (nameRefLine) ;

  graphRetitle ("Online Help") ;

  graphButton ("Quit", helpButtonQuit, 2, 1) ;

  graphButton ("Print", helpButtonPrint, 9, 1) ; /* mhmp 07.10.97 */

  goBackButton = graphButton (" << ", helpButtonGoBack, 17, 1) ;

  if (currPageListIndex == 0)
    graphBoxDraw (goBackButton, GRAY, WHITE);

  goForwardButton = graphButton (" >> ", helpButtonGoForward, 23, 1) ;

  if (currPageListIndex == arrayMax(pageList)-1)
    graphBoxDraw (goForwardButton, GRAY, WHITE);

  graphText ("middle button click to follow link in Netscape", 29, 1);

  graphButton ("Index", helpButtonIndex, WINX/8.0 - 3, 1) ;
  graphLine (0, 3, WINX/8.0 + 5, 3) ;

  graphRedraw() ;

  /********/

  graphActivate (helpSubGraph) ;
  graphClear () ;

  graphTextHeight (0) ;		/* default size for plain text */
  graphColor (BLACK) ;		/* color for plain text */
  graphTextFormat (PLAIN_FORMAT) ; /* format for plain text */

  indent = 8*2 ;
  xPos = indent; yPos = 2 ;
  lineHeight = 13+2 ;

  box2href = arrayReCreate (box2href, 10, char*) ;

  if (currPage)
    printSection (currPage->root) ;
  else
    {
      /* no page - all attempts to display the current file failed */
      messout ("Help unavailable");
    }
  
  /* an extra newline, if cursor isn't at start of line */
  blankLine () ;
  graphPixelBounds (WINMAXX + (2*8), yPos) ;
     
  graphRedraw() ;

  return ;
} /* helpDraw */
/************************************************************/

static void helpDestroy (void)
{
  int k ;

  for (k = 0; k < arrayMax(pageArray); ++k)
    htmlPageDestroy (arr(pageArray, k, HtmlPage*)) ;

  arrayDestroy (pageList) ;
  arrayDestroy (pageArray) ;
  dictDestroy (pageDict) ;
  dictDestroy (nameRefDict) ;
  keySetDestroy (nameRefLine) ;

  arrayDestroy (imageArray) ;
  dictDestroy (imageDict) ;

  arrayDestroy (box2href) ;

  menuDestroy (controlMenu);

  currPageListIndex = -1 ;
  currPageArrayIndex = -1 ;

  helpGraph = -1 ;
  helpSubGraph = -1 ;

  return ;
} /* helpDestroy */
/************************************************************/


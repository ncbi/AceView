/*  Last edited: Mar  5 01:57 1996 (rd) */

/* $Id: graphtest.c,v 1.1.1.1 2002/07/19 20:22:55 sienkiew Exp $ */
#include "regular.h"
#include "graph.h"

/***************************************************************/

void squeal (void)
{ printf ("Ahh! I die!!!!!!!!!\n") ; }

static Array newboxes = 0 ;

void point (int box, double x, double y)
{ 
  static int lastBox = 0 ;

  if (!newboxes)
    newboxes = arrayCreate (4, float) ;

  if (box && box == lastBox)
    printf ("Double click box %d\n", box) ;
  else
    printf ("Point at %f, %f in box %d\n", x, y, box) ;

  if (lastBox)
    { graphBoxDraw (lastBox, BLACK, WHITE) ;
      lastBox = 0 ;
    }

  if (box)
    { graphBoxDraw (box, WHITE, BLUE) ;
      lastBox = box ;
    }
  else
    { int i = graphBoxStart() ;
      graphRectangle (x-1, y-1, x+1, y+1) ;
      graphBoxEnd () ;
      graphBoxDraw (i, BLACK, WHITE) ;
      array(newboxes, i, float) = x ;
    }
}

void boxInfo (void)
{
  int i ;

  if (!newboxes) return ;
  for (i = 0 ; i < arrayMax(newboxes) ; ++i)
    printf ("Rectangle %d is at %.2f\n", i, arr(newboxes, i, float)) ;
}

void where (void)
{ 
  float x1, y1, x2, y2 ;
  
  printf ("selected at %f %f", graphEventX, graphEventY) ;
  graphWhere (&x1, &y1, &x2, &y2) ;
  printf ("in viewport %f, %f to %f, %f\n", x1, y1, x2, y2) ;
}

static float oldx,oldy ;

void start (double x, double y)
{ graphXorLine (0.5,0.5,x,y) ; oldx = x ; oldy = y ; }

void stop (double x, double y)
{ graphXorLine (0.5,0.5,oldx,oldy) ; }

void drag (double x, double y)
{ stop (x,y) ; start (x,y) ;}

static FREEOPT options[] = { 4,"Options",
			     'd',"Quit",
			     'f',"Finish",
			     'p',"Pop",
			     't',"Retitle",
			   } ;
static char countchar = 0;
static Graph blockGraph ;
void keyboard (KEY k)
{
  static char prompt[128] = "Initial test prompt" ;

  switch (k)
    { 
    case 'a':
      printf("simulating a RIGHT_DOWN\n");
#if !MACINTOSH
      graphEvent(RIGHT_DOWN, 0.0, 0.0);
#endif
      break;
    case 'b': graphColor (BLUE) ; break ;
    case 'd': graphDestroy() ; break ;
    case 't': graphRetitle("Brand new title") ; break ;
    case 's':
      printf ("I will simulate an 'e'\n") ; 
#if !MACINTOSH
      graphEvent ('e',0.,0.) ; 
#endif
      break ;
    case 'f':
      graphFinish() ;
      break; 
    case 'r': graphColor (RED) ; break ;
    case 'p': graphPop () ; break ;
    case '0': case '1': case '2': case '3': case '4': 
      graphGoto (k-'0',0) ; 
      break ;
    case 'm': graphMessage ("Testing the message system") ;
    case 'u': graphUnMessage () ; break ;
    case 'o': graphOut (messprintf ("%s%s%s",
	       "Here is a fairly long test message to see how ",
	       "the mechanism works.  Can it do better than on ",
	       "the Sun?")) ; break ;
    case 'q': 
      if (graphQuery (prompt))
	printf ("got TRUE\n") ;
      else
	printf ("got FALSE\n") ;
      break ;
    case 'x':
      if (graphPrompt ("New prompt for query:",prompt,"t"))
	strcpy (prompt,freeword()) ;
      break ;
    case 'y':
      if (graphSelect (&k, options))
	printf ("selected option with key '%c'\n",k) ;
      break ;
    case 'z': /* test all ascii events */
      {
	printf("simulating an %c\n", (char)countchar);
#if !MACINTOSH
    graphEvent(countchar++, 0.0, 0.0) ;
#endif
	break ;
      }
    case 'k':
      if (blockGraph)
	printf ("Answer from blocked graph was %d\n",graphBlock()) ;
      break ;
    case '8': case '9':
      graphUnBlock (k) ;
      break ;
	
    default:
      printf ("Got a key press %d = '%c'\n",k,k) ;
    }
}

static MENUOPT menu[] = {
  graphDestroy, "Quit",
  graphPrint, "Print Screen",
  boxInfo, "Box info", 
  where, "Where",
  0, 0 } ;

void main (int argc, char **argv)
{
  Graph g1, g2, g3, g4 ;

  graphInit (&argc, argv) ;

  g1 = graphCreate (PLAIN,"plain",0,0,0.4,0.3) ;
  graphPointsize (0.05) ;
  graphRegister (KEYBOARD, keyboard) ;
  graphRegister (MIDDLE_DOWN, start) ;
  graphRegister (MIDDLE_DRAG, drag) ;
  graphRegister (MIDDLE_UP, stop) ;
  graphText ("Are the circles round?",0.3,0.7) ;
  graphCircle (0.5,0.5,0.1) ;
  graphRedraw () ;

  g2 = graphCreate (TEXT_FIT,"text_fit",0,0,0.5 ,0.5) ;
  graphRectangle (0.5, 0.5, 75, 40) ;
  graphRegister (DESTROY, squeal) ;
  graphRegister (KEYBOARD, keyboard) ;
  graphRegister (PICK, point) ;
  graphMenu (menu) ;

  graphBoxStart () ;
  graphLine (10, 10, 20, 10) ;
  graphText ("label", 12, 8.5) ;
  graphBoxEnd () ;

  graphRedraw () ;
/*   graphGIFname ("graph2") ; */

  g3 = graphCreate (TEXT_SCROLL,"text_scroll",0,0,1,0.5) ;
  graphRectangle (0.5, 0.5, 75, 40) ;
  graphTextBounds (80,25) ;
  graphRegister (KEYBOARD, keyboard) ;
  graphButton ("Squeal", squeal, 1, 1) ;
  graphRegister (PICK, point) ;
  graphMenu (menu) ;

  g4 = graphCreate (TEXT_FULL_SCROLL,"text_full_scroll",0,0,0.4,0.6) ;
  graphRectangle (0.5, 0.5, 75, 40) ;
  graphTextBounds (120, 60) ;
  graphRegister (KEYBOARD, keyboard) ;
  graphRegister (PICK, point) ;
  graphMenu (menu) ;

  graphRampTool () ;

  graphStart (g2) ;

  graphFinish () ;
}
 

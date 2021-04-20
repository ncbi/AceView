/*  File: metab.c
 *  Author: Stan Letovsky (letovsky@gdb.org)
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:		Fully contributed by Stan Letovsky
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 14:05 1998 (fw)
 * * Jun 4 15:28 1996 (rbrusk):
 *	-	change fn() => fn(void) forward function declarations to please VC++ compiler
 * * Jun  3 15:51 1996 (rd)
 * * Jun  3 15:31 1996 (rd): changed assert to StanAssert because assert was abused
 * * Jun  3 15:31 1996 (rd): included metab.h directly
 * Created: Sun Nov  6 18:53:55 1994 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: metab.c,v 1.5 2017/02/15 20:36:40 mieg Exp $ */

/* Notes on coordinates:

diag* and metab* and diagram C-structures use homogenous coords -- 1
unit of x = 1 unit of y -- while graph* routines use character-based
coords -- 1 unit of x = width of a character, 1 unit of y = height of
a character. Thus coords or distances must be converted whenever going
back and forth between the diag* or metab* routines and the graph*
routines. The conversion is applied to y-values or distances only; the
homogenous units are = width of 1 character, so xh = xc.  The function
diagTextAspect() returns the aspect ratio of the font, ie the
height/width.  
    yc = yh/diagTextAspect(); yh = yc*diagTextAspect();


Note also that because y increases downwards, normal analytic geometry
must be corrected. In particular, either:
       angles should be negated, especially atan2 output
or
       sin(angle) and dy = y2-y1 inputs to atan2 should be negated.
In the code the latter approach is used: angles are still positive counterclockwise
with 0 to the right, but sins and atan dy's are negated.

Note that the external functions sinc, cos and atan2 use radians,
while graphArc uses degrees. Internally, angles stored in arc and
tangent structures are in degrees. Most internal functions
operate on radians. 
*/

#include "acedb.h"

#include "display.h"
#include "key.h"
#include "lex.h"
#include "bs.h"
#include "systags.h"
#include "classes.h"
#include "sysclass.h"
#include "tags.h"
#include "query.h"
#include "graphAcedbInterface.h"
#include <ctype.h>
#include <math.h>

#define NAME 60
#define ISIZE sizeof(int)

#define NODE_RADIUS 1
#define NODE_LABEL_DX 2.0
#define NODE_LABEL_DY -.5
#define DEBUG 0

/* these dims in character widths */
#define ARROW_HEAD_LENGTH 1.5
#define ARROW_HEAD_HALF_WIDTH 0.5

/* RD commented these out since now all static in one file

extern float X, Y, Pi, DegreesPerRadian,  RadiansPerDegree;
extern struct element *Elements;
extern int Made, Expand;
extern KEY NodeColor, ArcColor, CofColor;

*/

enum element_type { DIAG_TEXT, DIAG_NODE, DIAG_ARC, DIAG_TANGENT };
enum node_type { NORMAL_NODE, POINT_NODE, TANARC_NODE };

struct diagram
  { char *name; 
    KEY key;
    struct element *elements, *selection;
    float xmax, ymax ; /* of selected object within its box */
    Graph graph;
    int modified;
    AC_HANDLE handle ;
  };

struct element
  { int box, mark, highlight;
    float x, y, x0, y0, width, height; /* redundant representation of boxes
				 used by layout for fast virtual relocation */
    KEY key, color;
    enum element_type type; 
    void *item;
    struct element *next;
  };
  
struct node 
  {
    float x, y;
    struct element *label;
    enum node_type type;
    int id, depth, i, indeg, outdeg, mark, cycle;
  };

struct arc
  {
  float radius, xc, yc, x1, y1, start_angle, arc_angle;
  int id, arrow, /* 0 - no arrowhead; 1 - arrowhead on to-node */
     curvature; /* 0 = straight, 1 = curved, -1 = oppositely curved */
  struct element *label, *to, *from; 
  struct element *tangent, *tangent_to;
  };

struct tangent
{  float radius, xc, yc, xt, yt, start_angle, arc_angle /* angles in degrees */;
   int curvature;
   struct element *tangent_to;
   struct text *from, *to;
 };

struct text
  {
  float x, y;
  int id, outline;
  char *name;
  };

/* function prototypes */
static void diagQuit(void);
static void diagPick(double x, double y);
static void diagSelect(double x, double y);

static void diagMoveElement(double x, double y);
static void RecomputeArcParams(struct arc *a);
static struct element *ID2Node(int);
static struct element *ID2Arc(int);
static struct element *MakeText(int id, char *name, float x, float y, int outline, KEY color);
static struct element *MakeNode(int id, struct element *label, float x, float y);
static struct element *MakeArc(int id, struct element *from_node, struct element *to_node,
			struct element *label, int curvature);
static struct diagram *MakeDiagram(KEY);
static struct diagram *GetDiagram(KEY);

static void metabPathway(void);
static void metabRecompute(void);
static struct element *Pos2Element(struct diagram *d, float x, float y);
static void DrawDiagram(struct diagram *d);
static void metabDrawText( struct element *e );
static void DrawNode( struct element *e);
static void DrawArc( struct element *r);
static void DrawCurvedArc( struct element *r);
static void DrawStraightArc( struct element *r );
static void diagArrowHead( float x, float y, float dir);
static void diagMoveNode( struct element *e, float x, float y);
static void diagMoveText( struct element *e, float x, float y);
static void diagRedrawArc( struct element *e);
static void diagMoveArc( struct element *e, float x, float y);

static void diagRedraw(void);
static void RecomputeCurvedArcParams(struct arc *a);
static float diagTextAspect(void);
static void diagStraightArc( float xfrom, float yfrom, float xto, float yto );
static void LayoutDiagram(struct diagram *d);
static void SaveDiagram(struct diagram *d);
static void msg(char *msg);

static void *notNULL (void *x);
static BOOL diagBoxOverlap(struct element *a, struct element *b);
static BOOL diagFixOverlap(struct element *a, struct element *b);
static BOOL IntervalsOverlap(float min1, float max1, float min2, float max2);
static BOOL InInterval(float val, float min, float max);
static BOOL LineSegsIntersect(float x11, float y11, float x12, float y12, 
			     float x21, float y21, float x22, float y22, 
			     float *xi, float *yi);

static void metabExpand(void);
static struct diagram *ActiveDiagram(void);
static void ExpandDiagram(KEY key);
static void metabArcLabels(void);

static void metabMakeDiagram0(KEY key, struct diagram *d);
static BOOL Movable(struct element *e);

/************* tags and classes needed for this module ****************/

static KEY _Components ;
static KEY _Diagram ;
static KEY _Label ;
static KEY _Pathway ;
static KEY _Node ;
static KEY _Arc ;
static KEY _TangentArc ;
static KEY _Major_Reactant ;
static KEY _Major_Product ;
static KEY _Cofactor_Reactant ;
static KEY _Cofactor_Product ;
static KEY _Additional_Cofactor_Reactant ;
static KEY _Additional_Cofactor_Product ;

static int _VPathwayDiagram ;

static void initialiseTags (void)
{
  KEY classe ;

  lexaddkey ("Components", &_Components, 0) ;
  lexaddkey ("Diagram", &_Diagram, 0) ;
  lexaddkey ("Label", &_Label, 0) ;
  lexaddkey ("Pathway", &_Pathway, 0) ;
  lexaddkey ("Node", &_Node, 0) ;
  lexaddkey ("Arc", &_Arc, 0) ;
  lexaddkey ("TangentArc", &_TangentArc, 0) ;
  lexaddkey ("Major_Reactant", &_Major_Reactant, 0) ;
  lexaddkey ("Major_Product", &_Major_Product, 0) ;
  lexaddkey ("Cofactor_Reactant", &_Cofactor_Reactant, 0) ;
  lexaddkey ("Cofactor_Product", &_Cofactor_Product, 0) ;
  lexaddkey ("Additional_Cofactor_Reactant", 
	     &_Additional_Cofactor_Reactant, 0) ;
  lexaddkey ("Additional_Cofactor_Product", 
	     &_Additional_Cofactor_Product, 0) ;

  lexaddkey ("PathwayDiagram", &classe, _VMainClasses) ; 
  _VPathwayDiagram = KEYKEY(classe) ;
}

/*****************/

extern BOOL treeDisplay (KEY key, KEY from, BOOL isOldGraph);

static char DIAGRAM[]  = "DIAGRAM"; 
/*
static char ElementTypeName[][10] =  { "Text","Node","Arc", "Tangent" };
static char NodeTypeName[][10] =  { "Normal","Point","Tanarc" };
*/

static float X, Y, Pi = 3.141592, DegreesPerRadian = 57.29578,
  RadiansPerDegree = 0.017453293 ;
static struct element *Elements;
static KEY NodeColor = RED, ArcColor = BLUE, CofColor = DARKGREEN;
static int Expand = 0, ArcLabels = 1, FudgeLabels = 1, Made = 0, 
		Debug = 0, DebugBox = -1, Uncaching = 0;

static AC_HANDLE handle ;
#define metabAlloc(x)	halloc(x,handle)

/* define StanAssert rather than use regular assert because this code
   relies on the asserted expression being evaluated for correct execution,
   which is not ANSI, and is not true under MFC in WIN32
*/
#define StanAssert(x) if (x) {;} else messcrash ("StanAssert failed")

static void msg(char *msg) { if(Debug) printf("%s\n",msg); }

static float pythag(float a, float b) { return sqrt(a*a+b*b); }
static void *notNULL (void *x) { StanAssert(x!=NULL); return x; }

/* Count number of rows stored under tag. 
   Currently assumes first element of row is a key or tag. */
static int bsNRows(OBJ obj, KEY tag)
{ KEY  item; int n = 1;
  if(!bsGetKeyTags(obj,tag,&item)) return 0;
  while(bsGetKeyTags(obj,_bsDown,&item)) n++;
  return n; 
}

static MENUOPT DiagramMenus[2][2][10] = 
{
  {
    { 
      { diagQuit,      "Quit" },
      { help,          "Help" },
      { graphPrint,    "Print Window" },
      { metabPathway,  "Info"},
      { metabRecompute,"Recompute"},
      { metabExpand,   "Expand"},
      { metabArcLabels,"Show Arc Labels"},
      { diagRedraw,    "Redraw" },
      { displayPreserve, "Preserve" },
      { NULL,          NULL }
    },
    { 
      { diagQuit,      "Quit" },
      { help,          "Help" },
      { graphPrint,    "Print Window" },
      { metabPathway,  "Info"},
      { metabRecompute,"Recompute"},
      { metabExpand,   "Expand"},
      { metabArcLabels,"Hide Arc Labels"},
      { diagRedraw,    "Redraw" },
      { displayPreserve, "Preserve" },
      { NULL,          NULL}
    }
  },
  {
    {
      { diagQuit,      "Quit" },
      { help,          "Help" },
      { graphPrint,    "Print Window" },
      { metabPathway,  "Info"},
      { metabRecompute,"Recompute"},
      { metabExpand,   "No Expand"},
      { metabArcLabels,"Show Arc Labels"},
      { diagRedraw,    "Redraw" },
      { displayPreserve, "Preserve" },
      { NULL,          NULL}
    },
    {
      { diagQuit,      "Quit" },
      { help,          "Help" },
      { graphPrint,    "Print Window" },
      { metabPathway,  "Info"},
      { metabRecompute,"Recompute"},
      { metabExpand,   "No Expand"},
      { metabArcLabels,"Hide Arc Labels"},
      { diagRedraw,    "Redraw" },
      { displayPreserve, "Preserve" },
      { NULL,          NULL}
    }
  }
};


static char *ElementLabel(struct element *e)
{ struct element *t = NULL;
  if(e==NULL) return "NULL";
  switch(e->type){
  case DIAG_NODE: t = ((struct node *)e->item)->label; break;
  case DIAG_ARC: t = ((struct arc *)e->item)->label; break;
  case DIAG_TEXT: t = e; break;
  default: return "Unknown Elt Type";
 }
 return t==NULL?"No label":((struct text *) t->item)->name;
}

static struct diagram *ActiveDiagram(void)
{ struct diagram *d; 
  StanAssert (graphAssFind(DIAGRAM, &d));
  handle = d->handle ;
  return d; 
}

static void graphBoxRedraw(int b){ graphBoxDraw(b, -1, -1); } 

static void metabExpand(void)
{  Expand = ! Expand; 
  graphMenu(DiagramMenus[Expand][ArcLabels]);
/*
  graphSetMenuOptText(DiagramMenu,metabExpand,Expand?"No Expand":"Expand");
*/
}



static void diagRedraw(void)
{
  graphWhiteOut(); 
  graphRedraw();
}

/* not used
static void diagResize(void)
{ struct element *e = ActiveDiagram()->elements;
return;
  graphWhiteOut();
  while(e!=NULL)
    { if(e->type==DIAG_TEXT && ((struct text *)e->item)->outline)
	  { graphBoxClear(e->box);
	    metabDrawText(e);
	  }
	e = e->next;
    }
  graphRedraw();
}
*/

static float diagLineDir(float xfrom, float yfrom, float xto, float yto)
/* return direction in radians of a directed line segment */
{ return(float) (xto-xfrom)?atan2((double) -(yto-yfrom), 
				  (double) xto-xfrom ):
			  (yto > yfrom ? -Pi/2 : Pi/2); 
} 


static void metabPathway(void)
{ treeDisplay(ActiveDiagram()->key,0,0); }

/*
static int WarnSave = 0;
static void metabSave()
{ if(!WarnSave)
    { graphOut("Warning: You must also save changes from main menu.");
      WarnSave = 1; }
  SaveDiagram(ActiveDiagram()); }
*/

static void metabRecompute(void)
{ struct diagram *new = MakeDiagram(ActiveDiagram()->key);
  graphClear();			/* RD was graphDestroy */
  DrawDiagram(new);
}

static void diagQuit(void)
{ if (ActiveDiagram()->modified && graphQuery("Save layout?"))
   SaveDiagram(ActiveDiagram());
  handleDestroy(handle) ;
  graphDestroy(); 
}
    
/*
static void metabOptimize(void)
{  OptimizeLayout(ActiveDiagram()); 
}
*/

static void metabArcLabels(void)
{ struct diagram *d = ActiveDiagram(); 
  struct element *e = d->elements;
  struct arc *a;
  while(e!=NULL)
    { if(e->type == DIAG_ARC)
	{ a = e->item;
	  if(a->label!=NULL) graphBoxClear(a->label->box); }
      e = e->next; }
  ArcLabels = !ArcLabels;

  graphMenu(DiagramMenus[Expand][ArcLabels]);

/*
  graphSetMenuOptText(DiagramMenu,metabArcLabels,
		      ArcLabels?"Hide Arc Labels":"Show Arc Labels");
*/
}


static void diagDrag(float *x, float *y, BOOL up) 
{ 
/*  printf("Drag (%6.2f, %6.2f)\n",*x,*y); */
  if(up) diagMoveElement(*x,*y);
}

static void diagSelect(double x, double y)
{
  struct diagram *d = ActiveDiagram();
  struct element *e;
msg("+diagSelect");
  d->selection = e = Pos2Element(d, x, y);
  if(e != NULL)
    { X = x; Y = y;
      graphBoxDraw(e->box,GREEN,-1);
      if(e->type!=DIAG_ARC) graphBoxDrag(e->box, diagDrag); 
    }
msg("-diagSelect");
}


static void diagMoveElement(double x, double y)
{ 
  struct diagram *d = ActiveDiagram();
  struct element *e = d->selection;

  if(e==NULL) return;
  d->modified = 1;
/*
  printf("Moving box origin of %s to (%6.2f, %6.2f)\n", ElementLabel(e),x, y); 
*/
  y *= diagTextAspect();
  graphBoxDraw(e->box,BLACK,-1);
  switch(e->type){
  case DIAG_TEXT: diagMoveText(e,x,y); break;
  case DIAG_NODE: diagMoveNode(e,x,y); break;
  case DIAG_ARC: diagMoveArc(e,x,y); break;
  default: StanAssert(0);
  }
  
  diagRedraw();
  d->selection = NULL;
}

static struct element *ID2Node(int id)
{ struct element *l = Elements;
  while (l != NULL)
    { if(l->type == DIAG_NODE)
	  { if(((struct node *)l->item)->id == id) return l; }
	l = l->next; }
  printf("Can't find node %d\n",id);
  StanAssert(0); 
  return 0 ; /* i suppose ? */
}

static struct element *ID2Arc(int id)
{ struct element *l = Elements;
  while (l != NULL)
    { if(l->type == DIAG_ARC)
	  { if(((struct arc *)l->item)->id == id) return l; }
	l = l->next; }
  printf("Can't find arc %d\n",id);
  StanAssert(0); 
  return 0 ; /* i suppose ? */
}

static struct element *Key2Node(KEY key)
{ struct element *l = Elements;
  while (l != NULL)
    { if(l->type == DIAG_NODE && ((struct node *)l->item)->type == NORMAL_NODE)
	  { if(l->key == key) return l; }
	l = l->next; }
  return NULL;
}

static struct element *ID2Label(int id)
{ struct element *l = Elements;
  if(id==-1) return NULL;
  while (l != NULL)
    { if(l->type == DIAG_TEXT)
	  { if(((struct text *)l->item)->id == id) return l; }
	l = l->next; }
  return NULL;
}

static struct element *MakeElement(enum element_type type, void *item, KEY color)
{ struct element *e = (struct element *) metabAlloc (sizeof(struct element));
  e->key = 0;
  e->box = -1;
  e->item = item;
  e->type = type;
  e->next = Elements;
  e->color = color;
  e->highlight = 0;
  Elements = e;
  return e; 
}

static struct element *MakeText(int id, char *name, float x, float y, int outline, KEY color)
{ if(name==NULL) return NULL;
  else { struct text *t = (struct text *) metabAlloc(sizeof(struct text));
	 struct element *e;
	 t->id = id;
	 t->x = x;
	 t->y = y;
	 t->outline = outline;
	 t->name = name;
	 e = MakeElement(DIAG_TEXT,t,color);
	 e->width = (float)strlen(name);
	 e->height = diagTextAspect();
	 return e;
       }
}

static struct element *MakeNode(int id, struct element *label, float x, float y)
{  struct node *n = (struct node *) metabAlloc(sizeof(struct node));
   n->id = id;
   n->x = x;
   n->y = y;
   n->type = NORMAL_NODE;
   n->label = label;
   return MakeElement(DIAG_NODE,n,NodeColor);
 }

static struct element *MakePoint(int id)
{  struct element *e = MakeNode(id,NULL,0,0);
   struct node *n = e->item;
   n->type = POINT_NODE;
   e->color = ArcColor;
   return e; 
 }

static struct element *MakeArc(int id, struct element *from_node, struct element *to_node, 
				struct element *label,	int curvature)
{ struct arc *a = (struct arc *) metabAlloc(sizeof(struct arc));
  StanAssert(from_node!=to_node);
  StanAssert(from_node!=NULL&&to_node!=NULL);
  a->id = id;
  a->radius = 0;
  a->xc = 0;
  a->yc = 0;
  a->arrow = 1;
  a->start_angle = 0;
  a->arc_angle = 0;
  a->curvature = curvature; 
  a->from = from_node;
  a->to = to_node;
  a->label = label;
  a->tangent_to = a->tangent = NULL;
  RecomputeArcParams(a);
  return MakeElement(DIAG_ARC,a,ArcColor);
}

static void RecomputeTangentArcParams(struct arc *tan)
{ /* Recompute tangent params */
  struct arc *a = tan->tangent_to->item;
  struct node *from = tan->from->item, *to = tan->to->item;
  struct text *flab = from->label->item, *tlab = to->label->item;
  float x1 = ((struct node *) a->from->item)->x,
        x2 = ((struct node *) a->to->item)->x,
        y1 = ((struct node *) a->from->item)->y,
        y2 = ((struct node *) a->to->item)->y,
        xm = (x1+x2)/2, ym = (y1+y2)/2;
  float xn, yn, xt, yt, n;
  float th;
/*  msg("+RAP TAN"); */
  if(a->curvature)
    { float theta = RadiansPerDegree*(a->start_angle + 0.5*a->arc_angle), 
	  c = cos(theta), s = -sin(theta); 
	/* [c,s] is unit normal from parent arc (a) center to tangent point */
/*	msg("Curved"); */
	th = theta;
	/* find tangent point t */
	xt = a->xc + a->radius*c; yt = a->yc + a->radius*s;
	/* find tangent arc center */
	tan->xc = xt - a->curvature*tan->curvature*tan->radius*c;
	tan->yc = yt - a->curvature*tan->curvature*tan->radius*s; }
  else { /* Straight Arc */
/*    msg("Straight"); */
    xt = xm; yt = ym;
    if(ym==y1) { yn = 1; xn = 0; }
    else { xn = 1; yn = -(xm-x1)/(ym-y1); }
    if(y2<y1){ xn = -xn;  yn = -yn; } /* why is this needed? who cares? it is. */
    n = tan->curvature*tan->radius/pythag(xn,yn); 
    xn *= n; yn *=n;
    tan->xc = xm + xn;  tan->yc = ym + yn; 
  }
/*  msg("Common");	 */
  /* vector n is parallel to line(from,to) and has size = radius of tangent arc */
  xn = x2-x1; yn = y2-y1;
  n = tan->radius/pythag(xn,yn); xn *= n; yn *=n;
  /* recompute positions of tangent arc nodes = center +/- n */
  from->x = tan->xc - xn; from->y = tan->yc - yn;
  to->x = tan->xc + xn; to->y = tan->yc + yn;
  tan->x1 = to->x; tan->y1 = to->y;
/*  printf("Tan arc from (%6.2f,%6.2f) to (%6.2f,%6.2f)\n",
	   from->x,from->y,to->x,to->y); */
  /* Recompute tanarc node label positions.
     Want each end of arc to point to center of label.
     Let vector u = [xu,yu] = unit c - t (tanarc center - tangent point)
     = vector directed from center towards tanarc center, hence parallel
     to tangents to arc at from and to points, directed towards them.
     
     Position the label so that its center is at the origin.
     The vector k*u (k large) will intersect the bounding box of the label at p.
     If -p is translated to the node position, the origin will be translated
     to the desired label center. The label coords are upper left from there.
     Letting n be the node pos, n + p + l is desired position of label,
     where l goes from label center to label origin (upper left).
     */
  if(! Uncaching)
    { float xu = tan->xc - xt, yu = tan->yc - yt, k = 100000;
      float fudge = 1.05, h = fudge*0.5/diagTextAspect(), 
      w = fudge*(float)strlen(flab->name)/2, xp, yp;
      xu *= k; yu *= k; 
      if(!(LineSegsIntersect(0,0,xu,yu,-w,-h,w,-h,&xp,&yp) /* upper */
	   ||LineSegsIntersect(0,0,xu,yu,w,-h,w,h,&xp,&yp) /* right */
	   ||LineSegsIntersect(0,0,xu,yu,w,h,-w,h,&xp,&yp) /* lower */
	   ||LineSegsIntersect(0,0,xu,yu,-w,h,-w,-h,&xp,&yp))) /* left */
	{ xp = 0; yp = 0; }
      flab->x = from->x + xp - w; flab->y = from->y + yp - h;
      w = (float)strlen(tlab->name)/2;
      if(!(LineSegsIntersect(0,0,xu,yu,-w,-h,w,-h,&xp,&yp) /* upper */
	   ||LineSegsIntersect(0,0,xu,yu,w,-h,w,h,&xp,&yp) /* right */
	   ||LineSegsIntersect(0,0,xu,yu,w,h,-w,h,&xp,&yp) /* lower */
	   ||LineSegsIntersect(0,0,xu,yu,-w,h,-w,-h,&xp,&yp))) /* left */
	{ xp = 0; yp = 0; }
      tlab->x = to->x + xp - w; tlab->y = to->y + yp - h;
    }
  RecomputeCurvedArcParams(tan); 
  if(0&&Made)
    { float AR = diagTextAspect();
      if(DebugBox>0) graphBoxClear(DebugBox);
      DebugBox = graphBoxStart();
      diagStraightArc(a->xc,a->yc,a->xc+3*cos(th),a->yc-3*sin(th));
      graphCircle(a->xc,a->yc/AR,0.5); /* small - arc center */
      graphCircle(tan->xc,tan->yc/AR,1.0); /* med - tan center */
      graphCircle(xt,yt/AR,1.5); /* large - tan pt */
      graphBoxEnd();
      graphBoxDraw(DebugBox,CYAN,TRANSPARENT);
    }
  
  /*  msg("-RAP TAN"); */
}

static void RecomputeArcParams(struct arc *a)
{ struct text *label = (a->label!=NULL)?a->label->item:NULL;
  float x1 = ((struct node *)(a->from->item))->x,
  y1 = ((struct node *)(a->from->item))->y,
  x2 = ((struct node *)(a->to->item))->x,
  y2 = ((struct node *)(a->to->item))->y,
  xm = (x1+x2)/2, ym = (y1+y2)/2;
  if(a->tangent_to!=NULL) { RecomputeTangentArcParams(a); return; }
  /*  msg("+RAP"); 
      msg(label!=NULL?label->name:"Unlabelled"); */
  switch(a->curvature){
  case 0:  
    { /* msg("Straight"); */
      a->xc = x1; a->yc = y1;
      a->x1 = x2; a->y1 = y2; 
      a->radius = pythag(x1-xm,y1-ym);
      if(label!=NULL && !Uncaching)
	{ if(1) label->x = xm - (float)strlen(label->name)/2;
	else if(x2<x1&&y2!=y1) /* Position labels of left-sloping arcs on left side */
	  { label->x = xm-1.5*NODE_LABEL_DX-(float)strlen(label->name); 
	    if(label->x<0.5) label->x = 0.5; }
	else label->x = xm+1.5*NODE_LABEL_DX; 
	  label->y = ym-diagTextAspect(); }}
    break;
  case 1:
  case -1:
    /* Compute arc params (center, start_angle, delta_angle)
       given node coords and preferred radius. */
    { float r = a->radius, rmin = pythag(x2-x1,y2-y1)/2, xc, yc, xc1, yc1, 
	  theta ;
	float A, B, C, D; /* quadratic equation params */
	float c, d; /* slope and intercept of perpendicular to line(x1,xto) */
/*	msg("Curved"); */
	if(fabs(y1-y2)>0.0001)
	  { /* msg("y1!=y2"); */
	    c = (x1-x2)/(y2-y1); 
	    d = ym - c*xm;
	    A = 1 + c*c;
	    B = 2*c*d - 2*x1 - 2*y1*c;
	    C = d*d - r*r + x1*x1 + y1*y1 - 2*y1*d;
	    D = B*B-4*A*C;
	    if(D>0)
		{ float yr;
/*		  msg("D>0"); */
		  /* 2 solutions, pick one */
		  D = sqrt(D);
		  /* (xc,yc) and (xc1,yc1) are the two solutions */
		  xc = (-B + D)/(2*A);
		  yc = c*xc + d;
		  xc1 = (-B - D)/(2*A);
		  yc1 = c*xc1 + d;
		  /* Translate (xc,yc) to the origin, rotate by - the slope of line(from,to). */
		  theta = diagLineDir(x1,y1,x2,y2); 
		  yr = (xc-xm)*sin(theta) + (yc-ym)*cos(theta);
		  /* Should now be (0,h); if sign(h) disagrees w/ curvature, use (xc1,yc1) */
		  if(!((a->curvature==1&&yr<0)||(a->curvature==-1&&yr>0)))
		    { xc = xc1; yc = yc1; } 
		}
	    else 
		{ /* 0 or 1 solutions: use midpoint of line between nodes; semi-circular arc */
/*		  msg("0/1 solns"); */
		  r = a->radius = rmin;
		  xc = xm; yc = ym;
		}
	  }
	else 
	  { /* line(from,to) is horizontal. 
	       perpendicular is vertical; centers will be on it. */
/*	    msg("horizontal"); */
	    xc = xm;
	    D = r*r - (xm-x1)*(xm-x1);
	    if(D>0)
		{ /* 2 solutions (xc,y1+D) and (xc,y1-D), pick one */
/*		  msg("2 solutions"); */
		  D = sqrt(D);
		  yc = y1 + D*((x2-x1>0)?-1:1)*a->curvature; }
	    else 
		{ /* 0 or 1 solution; use midpoint of line between nodes; 
		     semi-circular arc */
/*		  msg("0/1 solution"); */
		  r = a->radius = rmin;
		  xc = xm; yc = ym;
		}
	  }
/*	msg("Join"); */
	/* Store center coords */
	a->xc = xc; a->yc = yc; 
	a->x1 = x2; a->y1 = y2; 
/*	msg("Call RCAP"); */
	RecomputeCurvedArcParams(a); }
    break;
  default:
    printf("Invalid curvature %d\n",a->curvature);
  }
  if(a->tangent!=NULL) 
    { struct arc *tan = a->tangent->item;
      if(a->curvature)
	{ tan->radius = (a->radius * a->arc_angle * RadiansPerDegree)/4.0;
	  if(tan->radius < 2 ) tan->radius = 2; }
      else tan->radius = a->radius/4;
      if(tan->radius>6) tan->radius = 4;
      RecomputeTangentArcParams(tan);
    }
/*  msg("-RAP"); */
}

static void RecomputeCurvedArcParams(struct arc *a)
{/* Determine start and arc angles; label position, from radius, curvature & endpoints */
  float x1 = ((struct node *)(a->from->item))->x,
        y1 = ((struct node *)(a->from->item))->y,
        x2 = ((struct node *)(a->to->item))->x,
	  y2 = ((struct node *)(a->to->item))->y,
        xc = a->xc, yc = a->yc, r = a->radius;
  float a1, a2, da;
  struct text *label = (a->label!=NULL && !Uncaching)?a->label->item:NULL;
  a1  = diagLineDir(xc,yc,x1,y1);
  a2 = diagLineDir(xc,yc,x2,y2);
  if(a->curvature == -1){ float t = a1; a1 = a2; a2 = t; } 
  a->start_angle = a1*DegreesPerRadian;
  da = a2 - a1;
  if(da<0) da+=2*Pi;
  a->arc_angle = da*DegreesPerRadian; 
  if(label!=NULL)
    { /* Position label */
	a1 = -(a1+0.5*da);
	label->x = xc + (r+NODE_LABEL_DX)*cos(a1);
	label->y = yc + (r+NODE_LABEL_DX)*sin(a1);  }
}


static BOOL DirWithinArc(float dir, float a0, float da)
/* return true if dir is within arc from a0 to a0 + da (in radians) */
{ while(dir<a0) dir += 2*Pi;
  while(dir>a0+2*Pi) dir -= 2*Pi;
  return dir <= a0 + da; }

static BOOL PointWithinBox(float x, float y, int box)
{ float xmin, ymin, xmax, ymax;
  graphBoxDim(box, &xmin, &ymin, &xmax, &ymax);
  return (xmin <= x && xmax >= x && ymin <= y && ymax >= y); }

static float ElementScore(struct element *e, float x, float y)
{ switch(e->type){
 case DIAG_NODE:
 case DIAG_TEXT:
  if(PointWithinBox(x,y,e->box))
    { if(e->type==DIAG_NODE) return 100; else return 50; }
  else return 0; 
 case DIAG_ARC: 
  { float score = 0, distance = 1000; struct arc *a = e->item; 
    y = y*diagTextAspect();
    if(a->curvature)
      { if(DirWithinArc(diagLineDir(a->xc,a->yc,x,y),
			  a->start_angle*RadiansPerDegree,
			  a->arc_angle*RadiansPerDegree))
	    distance = fabs(pythag(y-a->yc,x-a->xc)-a->radius); }
    else /* Check distance from cursor to straight arc <= 2 */
	{ float x1 = a->xc, y1 = a->yc, x2 = a->x1, y2 = a->y1;
	  float dx = x2 - x1, dy = y2 - y1;
	  if(dx==0) { if( ((y1<=y)&&(y<=y2)) || ((y2<=y)&&(y<=y1))) distance = fabs(x-x1); }
	  else if(dy==0) { if( ((x1<=x)&&(x<=x2))||((x<=x)&&(x<=x1))) distance = fabs(y-y1); }
	  else { float a, b, c, d, xm, ym, dxm, dym, dot1212, dot121m, epsilon = 0.000001;
		   StanAssert(fabs(dx)>epsilon&&fabs(dy)>epsilon);
		   a = dy/dx; b = y1 - a * x1; /* y=ax+b is line(from,to) */
		   c = -1/a;  d = y - c*x;  /* y = cx+d is normal to that thru cursor pos */
		   xm = (d-b)/(a-c); ym = a*xm + b; /* (xm,ym) is intersection of the two lines */
		   dxm = xm - x1; dym = ym - y1; /* vector(dxm,dym) goes from x1 to xm */
		   /* if dot-product 121m of that vector with the from-to vector is positive
			and <= dod-product 1212 of from-to with itself, then m is within the from-to
			interval */
		   dot1212 = dx*dx + dy*dy; dot121m = dxm*dx + dym*dy;
		   if(dot121m>=0&&dot121m<=dot1212) distance = pythag(x-xm,y-ym); }}
    if(distance<=2) score = 1/(.02+distance); 
    if(a->tangent_to!=NULL) score = 0.8*score;
/*    printf("Score for %s %s = %6.3f\n",ElementTypeName[e->type],ElementLabel(e),score); */
    return score; }
 default:
  return 0 ;
}}

static struct element *Pos2Element(struct diagram *d, float x, float y)
{ struct element *e = d->elements, *best = NULL;
  float score, best_score = 0;
  while(e != NULL)
    { score = ElementScore(e,x,y);
	if(score>best_score) { best_score = score; best = e; }
      e = e->next;
    }
  return best;
}

static void diagMaxDims(struct diagram *d, struct element *e)
{ /* update diag xmax, ymax to reflect position of element. */
  float x,y;
  switch(e->type)
    {
    case DIAG_NODE: 
      { 
	struct node *n = e->item;
	x = n->x + NODE_RADIUS;
	y = n->y + NODE_RADIUS; 
	y *= 1/diagTextAspect();
      }
      break;

    case DIAG_ARC: 
      return; /* endpoint nodes will suffice */
      
    case DIAG_TEXT: 
      { 
	float x0, y0;
	graphBoxDim(e->box,&x0,&y0,&x,&y); 
      }
    default: break;
    }
  if(d->xmax < x) d->xmax = x;
  if(d->ymax < y) d->ymax = y; 
}

static void DrawDiagram(struct diagram *d)
{ struct element *e = d->elements;
  msg("+DrawDiagram");
  d->graph = graphActive();
  graphRetitle(d->name);
  graphTextHeight(0); 
  graphAssociate(DIAGRAM, d);
  graphPointsize(0.5); 
/*  graphRegister (LEFT_DOWN, diagSelect) ; */
/*  graphRegister (LEFT_UP, diagUnselect) ; */
  graphRegister (LEFT_DOWN, diagPick) ;
  graphRegister (MIDDLE_DOWN, diagSelect) ;
  graphRegister (MIDDLE_UP, diagMoveElement) ;
  graphMenu(DiagramMenus[Expand][ArcLabels]);
  d->xmax = d->ymax = 0;
  e = d->elements;
  while( e != NULL)
    {	/* printf("%s %s\n",ElementTypeName[e->type], ElementLabel(e)); */
	switch(e->type){
	case DIAG_NODE: DrawNode( e ); break;
	case DIAG_ARC: DrawArc( e ); break;
	case DIAG_TEXT: metabDrawText (e ); break;
	default: printf("something other than TEXT, NODE id found\n");
	  StanAssert(0); } 
	diagMaxDims(d,e);
	e = e->next;
    }
  d->ymax *= 1.1;
  graphTextBounds(d->xmax, d->ymax);
  graphRedraw ();
}

static void metabDrawText( struct element *e)
{ struct text *t = e->item;
  float x = t->x, y = t->y/diagTextAspect();
  graphColor(e->color);
  e->box = graphBoxStart();
  graphBoxDraw(e->box,-1,TRANSPARENT);
  graphText(t->name, x, y);
  if(t->outline) 
    graphRectangle(x-0.3,y-0.1, x+((float)strlen(t->name)+0.4),y+1.1);
  graphBoxEnd();
  if (e->key)
    graphBoxInfo (e->box, e->key, name(e->key)) ;
}

static void DrawNode( struct element *e)
{  struct node *n = e->item;
   float x = n->x, y = n->y/diagTextAspect();
   graphColor(e->color);
   e->box = graphBoxStart();
   graphBoxDraw(e->box,-1,TRANSPARENT);
   switch(n->type)
     {
     case NORMAL_NODE: 
       graphCircle (x, y, NODE_RADIUS); 
       break;

     case POINT_NODE:  
       graphPoint(x,y);
       break;

     default: break;
     }
   graphBoxEnd();
}

static float diagTextAspect(void)
{ /* returns the ratio of x/y; multiply y-distances times this to get homogenous coords. */
  int dx = 0, dy;
  float width, height, aspect = 1.625;
  if(graphActive()) 
    { graphTextInfo(&dx, &dy, &width, &height); 
	if(dx) aspect =  ((float) dy)/dx; }
  return aspect;
}

static void DrawArc( struct element *e )
{ struct arc *a = e->item;
  if(a->curvature) DrawCurvedArc( e ); else DrawStraightArc( e ); 
}
/*
static float NormalizeRadians(float r)
{ r = fmod(r,(Pi*2));
  if(r<0) r += Pi*2;
  return r; 
}
*/
static void DrawCurvedArc( struct element *e ) 
{ struct arc *a = e->item;
  float AR = diagTextAspect(), xc = a->xc, yc = a->yc, xto = a->x1, yto = a->y1;
  graphColor(e->color);
  e->box = graphBoxStart(); 
  graphBoxDraw(e->box,-1,TRANSPARENT);
  { struct node *from = a->from->item, *to = a->to->item; 
    float x1, y1, x2, y2;
    x1 = from->x; y1 = from->y; x2 = to->x; y2 = to->y;
    if(Debug) printf("Drawing %s Arc from %s at (%7.3f,%7.3f) to %s \n\
       at (%7.3f,%7.3f) c=(%7.3f,%7.3f), r=%7.3f, a0=%7.3f da=%7.3f\n",
	   a->curvature>0?"+":"-",
	   ElementLabel(a->from),x1,y1,ElementLabel(a->to),x2,y2, xc,  yc,
	   a->radius,  a->start_angle,  a->arc_angle);
  }
#ifdef WWW_ACEDB
#define sign(x) (((x)>0)?1:-1)
  graphFillArc(xc,yc/AR,a->radius,sign(a->curvature)*a->start_angle,
               sign(a->curvature)*a->arc_angle);
/*
  graphFillArc(xc,yc/AR,a->radius,a->start_angle,a->arc_angle); 
*/
#else
  graphArc(xc,yc/AR,a->radius,a->start_angle,a->arc_angle); 
#endif /* WWW_ACEDB */
  if(a->arrow) 
    { /* used to use tangent at to-point as arc direction, but this looks bad
	   for small arcs. Better would be direction of chord of length of arrowhead
	   ending at to. Approximate this by length of arc, which is angle in radians
	   times radius. p is at other end of this arc. */
	double ato, atop, ap ; float xp, yp ;
	/* ato = direction from center to "to" point. */
	ato = diagLineDir(xc,yc,xto,yto);
	/* atop = direction from "to" to p = point 1 arrowhead of arc away */
	atop = ARROW_HEAD_LENGTH / a->radius;
	/* ap = direction from p to to = arrowhead direction = ato +/- ap */
	ap = ato - a->curvature * atop;
	/* p is point on circumference at angle atop away from to; p-to
	   arc length = arrowhead length */
	xp = xc + a->radius*cos(ap); yp = yc - a->radius*sin(ap);
	/* arrowhead direction is from p to to */
	diagArrowHead( xto, yto, diagLineDir(xp,yp,xto,yto));
    }
  graphBoxEnd();
}

static void DrawStraightArc( struct element *e )
{ struct arc *a = e->item;
  float AR = diagTextAspect(), xfrom = a->xc, yfrom = a->yc, xto = a->x1, yto = a->y1;
  graphColor(e->color);
  e->box = graphBoxStart();
  graphBoxDraw(e->box,-1,TRANSPARENT);
  graphLine(xfrom, yfrom/AR, xto, yto/AR);
  if(a->arrow) diagArrowHead( xto, yto, diagLineDir(xfrom,yfrom,xto,yto));
  graphBoxEnd();
}

/* used for debugging only at present */
static void diagLine( float xfrom, float yfrom, float xto, float yto )
{ float AR = diagTextAspect();
  graphLine(xfrom, yfrom/AR, xto, yto/AR);
}

/* used for debugging only at present */
static void diagStraightArc( float xfrom, float yfrom, float xto, float yto )
{ msg("+diagStraightArc"); 
  diagLine(xfrom, yfrom, xto, yto);
  diagArrowHead( xto, yto, diagLineDir(xfrom,yfrom,xto,yto));
}


/* Draw an arrow head with point at (x,y) in direction dir (in radians).
   Template position of arrowhead is facing right, point at 0,0. */
static void diagArrowHead( float x, float y, float dir)
{ 
  int i,j;
  float tx[3] = { -ARROW_HEAD_LENGTH, 0, -ARROW_HEAD_LENGTH},  
  	ty[3] = { ARROW_HEAD_HALF_WIDTH, 0, -ARROW_HEAD_HALF_WIDTH },
        cosdir = cos(dir), sindir = -sin(dir) ;
  Array ArrowHead = arrayCreate(8,float);
  for(i=0;i<4;i++) 
    { j = (i<3)?i:0;
	/* Rotation:
	      x' = x*cos(d) - y*sin(d);
		y' = x*sin(d) + y*cos(d) */
	array(ArrowHead,2*i,float) =   (float) x + tx[j]*cosdir - ty[j]*sindir;  
	array(ArrowHead,2*i+1,float) = (float) (y + tx[j]*sindir + ty[j]*cosdir)/diagTextAspect();  
    }
  graphPolygon(ArrowHead);
}

/* Move node and label only */
static void diagMoveNode0(struct element *e, float x, float y)
{ struct node *n = e->item;
  struct element  *tr = n->label;

  float AR = diagTextAspect(),  xmin, ymin, xmax, ymax;
  graphBoxDim(e->box, &xmin, &ymin, &xmax, &ymax);
  graphBoxShift(e->box, x, y/AR);
  n->x += x-xmin;  n->y += y-ymin*AR;
  if(tr!=NULL) diagMoveText(tr, x + NODE_LABEL_DX, y + NODE_LABEL_DY);
}

/* Move node, label and arcs */
static void diagMoveNode(struct element *e, float x, float y)
{ struct element *e1;
/*
   printf("Moving %s\n",ElementLabel(e));
*/
  diagMoveNode0(e,x,y);
  e1 = ActiveDiagram()->elements;
  while(e1!=NULL)
    { if(e1->type == DIAG_ARC &&
	   ( ((struct arc *)e1->item)->to == e || ((struct arc *)e1->item)->from == e))
	  { RecomputeArcParams(e1->item);
	    diagRedrawArc(e1); }
	e1 = e1->next;
    }
}

static void diagMoveText( struct element *e, float x, float y)
{ struct text *t = e->item; 
  graphBoxShift(e->box, x, y/diagTextAspect() );
  t->x = x; t->y = y;
}

static void diagMoveArc( struct element *e, float x, float y)
{ struct arc *a = e->item;  
  struct node *from = a->from->item, *to = a->to->item;
  float AR = diagTextAspect();
  float x0 = X, y0 = Y*AR, 
        x1 = from->x, y1 = from->y, x2 = to->x, y2 = to->y, a1, b1, a2, b2;
  float d, xm = (x1+x2)/2, ym = (y1+y2)/2, xc = a->xc, yc = a->yc,
        xc0 = a->curvature?a->xc:xm, yc0 = a->curvature?a->yc:ym, 
        r = a->radius, rmin =  pythag(x2-x1,y2-y1)/2;
  float cm, xcm, ycm, xf, yf;
  if(a->tangent_to!=NULL) 
    { struct arc *parent = a->tangent_to->item;
	struct node *pfrom = parent->from->item, *pto = parent->to->item;
	struct element *tfrom = ((struct node *)a->from->item)->label, 
	               *tto = ((struct node *)a->to->item)->label;
	float xt = (pfrom->x+pto->x)/2.0, yt=(pfrom->y+pto->y)/2.0;
	a->radius = pythag(x-xt,y-yt);
	if((x-xt)*(a->xc-xt)+(y-yt)*(a->yc-yt)<0) a->curvature *= -1;
	RecomputeArcParams(a);
	diagRedrawArc(e);
	diagMoveText(tfrom,((struct text *)tfrom->item)->x,
			 ((struct text *)tfrom->item)->y);
	diagMoveText(tto,((struct text *)tto->item)->x,((struct text *)tto->item)->y);
	return;
    }
  /* Form equations of lines perpendicular to line(from,p) and line(to,p) */
  a1 = -(x1-x)/(y1-y); b1 = (y+y1)/2-a1*(x+x1)/2;
  a2 = -(x2-x)/(y2-y); b2 = (y+y2)/2-a2*(x+x2)/2;
  StanAssert(a1!=a2); /* p not colinear with line(from,to) */
  /* Center (xc,yc) is intersection of those lines. */
  xc = (b2-b1)/(a1-a2); yc = a1*xc + b1;
  /* cm is normalized vector from c to m */
  xcm = xm - xc; ycm = ym - yc;
  cm = pythag(xcm, ycm);
  cm = 1/cm; xcm *= cm; ycm *= cm;
  d = xcm*(x-xm) + ycm*(y-ym); /* d = scalar_product(cm,mc) */
  /* f is as distant from c as p, but lies on normal(line(from,to)) */
  xf = xm + xcm*d; yf = ym + ycm*d; 
  if(pythag(xf-xm,yf-ym)<=0.8)
    { /* Become a straight arc */
	a->curvature = 0;
	xc = x1; yc = y1; a->x1 = x2; a->y1 = y2;  }
  else { float theta, yr;
	   r = pythag(xc-x,yc-y);
	   if(r>rmin&&xcm*(xf-xm)+ycm*(yf-ym)<0) r = rmin;
	   /* Translate m to the origin, rotate f by - the slope of line(from,to). */
	   theta = diagLineDir(x1,y1,x2,y2); 
	   yr = (xc-xf)*sin(theta) + (yc-yf)*cos(theta);
	   if(yr<0) a->curvature = 1; else a->curvature = -1;
	 }

  if(0)
    { if(DebugBox>0) graphBoxClear(DebugBox);
	DebugBox = graphBoxStart();
	graphCircle(xc0,yc0/AR,0.5); /* OLD CENTER IS SMALL CIRCLE */
	graphCircle(xc,yc/AR,1);
	graphRectangle(xf-.5,(yf-.5)/AR,xf+.5,(yf+.5)/AR);
	graphLine(x1,y1/AR,x,y/AR); graphLine(x2,y2/AR,x,y/AR);
	graphLine((x1+x)/2,(y1+y)/(2*AR),xc,yc/AR); graphLine((x2+x)/2,(y2+y)/(2*AR),xc,yc/AR);
	diagStraightArc(xc0,yc0,xc,yc);
	diagStraightArc(x0,y0,x,y);
	graphBoxEnd();
	graphBoxDraw(DebugBox,CYAN,TRANSPARENT);
	printf("MoveArc setting r=%6.2f (rmin=%6.2f) \n",r,rmin);
    }
  a->xc = xc; a->yc = yc; a->radius = r;
  RecomputeArcParams(a);
  diagRedrawArc(e);
}
  
static void diagRedrawArc(struct element *e)
{  struct arc *a = e->item;
   struct text *t = a->label!=NULL?a->label->item:NULL;

   graphBoxClear(e->box);
   DrawArc(e);
   graphBoxRedraw(e->box);
   if(a->tangent!=NULL)
     { struct arc *tan = a->tangent->item;
	 struct node *from = tan->from->item, *to = tan->to->item;
	 diagRedrawArc(a->tangent);
	 diagMoveText(from->label,((struct text *)from->label->item)->x, 
			  ((struct text *)from->label->item)->y);
	 diagMoveText(to->label,((struct text *)to->label->item)->x, 
			  ((struct text *)to->label->item)->y);
     }
   if(t!=NULL) diagMoveText(a->label,t->x, t->y);
}

static BOOL metabInDiagram(KEY rxn, struct diagram *d)
{  struct element *e = d->elements; 
   while(e!=NULL) { if(e->key == rxn) return TRUE; e = e->next; }
msg("Couldn't find it!");
   return FALSE; }

BOOL metabDisplay (KEY key, KEY from, BOOL isOldGraph)
{ 
  OBJ obj = bsCreate(key); 
  KEY diagkey;
  struct diagram *d;
  static BOOL isFirst = TRUE ;

  if (isFirst)
    { initialiseTags () ;
      isFirst = FALSE ;
    }

  if(bsFindTag(obj,_Components)) /* if is pathway... */
    { 
      if (isOldGraph && !ActiveDiagram()->modified)
	{ handleDestroy(handle); /* re-use graph, free old data structures */
	  graphClear();
	}
      else
	displayCreate(METAB);
      if(Expand && graphAssFind(DIAGRAM, &d) && metabInDiagram(key,d))
	ExpandDiagram(key);
      else DrawDiagram(bsGetKey(obj,_Diagram,&diagkey)
		       ?GetDiagram(diagkey)
		       :MakeDiagram(key)); 
      bsDestroy (obj) ;
      return TRUE ;
    }
  treeDisplay (key, from, FALSE) ;
  bsDestroy (obj);
  return FALSE ;
}

static int EltVal(struct element *a)
{ switch(a->type){
 case DIAG_NODE: return ((struct node *)a->item)->id;
 case DIAG_ARC: return ((struct arc *)a->item)->id;
 case DIAG_TEXT: return ((struct text *)a->item)->id;
 default: StanAssert(0); return 0 ;
 }
}

static BOOL EltLT (struct element *a, struct element *b)
{ if(a->type == b->type) return EltVal(a) < EltVal(b);
  else return a->type < b->type;
}

/* probably not necessary -- trying to get at caching bug */
static void SortDiagram(struct diagram *d)
{ int done = 0;
  struct element **e, *tmp;
  if(d->elements==NULL) return;
  while(!done)
    { done = 1;
      e = &d->elements;
      while((*e)->next!=NULL)
	{ if(EltLT(*e,(*e)->next))
	    { tmp = *e;
	      (*e) = (*e)->next;
	      tmp->next = (*e)->next;
	      (*e)->next = tmp;
	      done = 0; }
	  e = &(*e)->next; }
    }
/*  { struct element *e = d->elements;
    while(e!=NULL)
      { printf("%s %s %d\n",ElementTypeName[e->type],ElementLabel(e),EltVal(e));
	e = e->next; }
  } */

}

/* Load a cached diagram structure from the DB;
   Create C-structures for diagram from stored B-structure. */
static struct diagram *GetDiagram(KEY key)
{  KEY nextKey, k; OBJ obj = bsCreate(key);  BSMARK pos = NULL;
   struct diagram *d ;
   struct element *e;
   int id, lab, done, curvature; 

   Made = 0;
   handle = handleCreate() ;
   d = (struct diagram *) metabAlloc(sizeof(struct diagram));
   d->handle = handle ;
   d->modified = 0;
   d->selection = NULL;
msg("+GetDiagram");
   Uncaching = 1;
   StanAssert(bsGetKey(obj,_Pathway,&nextKey));
   d->name = name(nextKey);
   d->key = nextKey;
   Elements = NULL;

   /* make Labels */
   StanAssert(bsGetData(obj,_Label,_Int,&id));
   done = 0;
   while(!done)
     { char *text, *buf; float x, y ;
       pos = bsMark(obj, pos);  
       StanAssert(bsGetData(obj,_bsRight,_Float,&x));
       StanAssert(bsGetData(obj,_bsRight,_Float,&y));
       StanAssert(bsGetData(obj,_bsRight,_Text,&text));
       buf = (char*) metabAlloc((strlen(text)+1)*sizeof(char));
       strcpy(buf,text);
       MakeText(id, buf, x, y, 0, NodeColor);
       bsGoto(obj, pos);
       done = ! bsGetData(obj,_bsDown,_Int,&id);
     }
   /* make Nodes */
   StanAssert(bsFindTag(obj,_Node));
   StanAssert(bsGetKeyTags(obj,_bsRight,&nextKey));
   done = 0;
   while(!done)
     { float x, y ; 
       pos = bsMark(obj, pos);	
       StanAssert(bsGetData(obj,_bsHere,_Int,&id));
       StanAssert(bsGetData(obj,_bsRight,_Float,&x));
       StanAssert(bsGetData(obj,_bsRight,_Float,&y));
       StanAssert(bsGetData(obj,_bsRight,_Int,&lab));
       e = MakeNode(id, ID2Label(lab), x, y);  
       if(bsGetKey(obj,_bsRight,&k)) e->key = k;
       bsGoto(obj, pos);
       done = ! bsGetKeyTags(obj,_bsDown,&nextKey);
     }
   /* make Arcs */
   if(bsFindTag(obj,_Arc) &&bsGetKeyTags(obj,_bsRight,&nextKey))
     { done = 0;
       while(!done)
	 { int from, to;
	   float radius; struct element *e; struct arc *a;
	   pos = bsMark(obj, pos);	
	   StanAssert(bsGetData(obj,_bsHere,_Int,&id));
	   StanAssert(bsGetData(obj,_bsRight,_Int,&lab));
	   StanAssert(bsGetData(obj,_bsRight,_Int,&from));
	   StanAssert(bsGetData(obj,_bsRight,_Int,&to));
	   StanAssert(bsGetData(obj,_bsRight,_Int,&curvature));
	   StanAssert(bsGetData(obj,_bsRight,_Float,&radius));
	   e = MakeArc(id,ID2Node(from),ID2Node(to),ID2Label(lab),curvature);
	   if(bsGetKey(obj,_bsRight,&k)) e->key = k;
	   a = e->item; a->radius = radius; 
	   bsGoto(obj, pos);
	   done = ! bsGetKeyTags(obj,_bsDown,&nextKey);
	 }}

   /* make Tangents */
   if(bsFindTag(obj,_TangentArc))
      { StanAssert(bsGetKeyTags(obj,_bsRight,&nextKey));
	done = 0;
	while(!done)
	  { int arcid, tanid; struct element *arcelt, *tanelt;
	    struct arc *arc, *tan;
	    pos = bsMark(obj, pos);	
	    StanAssert(bsGetData(obj,_bsHere,_Int,&arcid));
	    StanAssert(bsGetData(obj,_bsRight,_Int,&tanid));
	    arcelt = ID2Arc(arcid); tanelt = ID2Arc(tanid);
	    arc = arcelt->item; tan = tanelt->item;
	    arc->tangent = tanelt; tan->tangent_to = arcelt;
	    ((struct node *) tan->from->item)->type = TANARC_NODE;
	    ((struct node *) tan->to->item)->type = TANARC_NODE;
	    bsGoto(obj, pos);
	    done = ! bsGetKeyTags(obj,_bsDown,&nextKey);
	  }}

   /* set implied attributes */
   e = Elements;
   while(e!=NULL)
     { struct element *e1; struct text *label;
       if(e->type==DIAG_NODE&&(e1=((struct node *)e->item)->label)!=NULL)
	 { e1->key = e->key; 
	   label = e1->item; 
	   e1->color = ((struct node *) e->item)->type == TANARC_NODE?CofColor:NodeColor;
	 }
       else if(e->type == DIAG_NODE)
	 ((struct node *) e->item)->type = POINT_NODE;
       else if(e->type==DIAG_ARC&&(e1=((struct arc *)e->item)->label)!=NULL)
	 { e1->key = e->key;
	   label = e1->item; 
	   e1->color = ArcColor; 
	   label->outline = 1; }
       if(e->type == DIAG_ARC
	  && ((struct node *)((struct arc *)e->item)->to->item)->label==NULL)
	 ((struct arc *) e->item)->arrow = 0;
       if(e->type == DIAG_ARC) RecomputeArcParams(e->item);
       e = e->next;
     }

   d->elements = Elements;
msg("-GetDiagram");
   bsDestroy(obj);
   Made = 1;
/*   SortDiagram(d); */
   Uncaching = 0;
   return d;
}

void bsDump (OBJ obj);

static void SaveDiagram(struct diagram *d)
{ KEY pwkey = d->key, diagkey, addpos;
  OBJ pw, diag; char *dname;
  BSMARK pos = NULL; int anytans = 0;
  struct element *e;
  msg("+SaveDiagram");
  pw = notNULL(bsUpdate(pwkey));
  msg(name(pwkey));
  dname = (char*) metabAlloc(strlen(name(pwkey))+20);
  sprintf(dname,"%s Diagram Info",name(pwkey));
  lexaddkey(dname,&diagkey,_VPathwayDiagram);  
  StanAssert(diagkey!=0);
  diag = bsUpdate(diagkey);  
  StanAssert(diag!=NULL);

  if(!bsFindKey (pw,_Diagram,diagkey)) 
    { StanAssert(bsAddKey(pw,_Diagram,diagkey)); }

  if(bsFindTag(diag,_Pathway))
    { bsRemove(diag); }
  StanAssert(bsAddTag(diag,_Pathway));
  bsAddKey(diag,_Pathway,pwkey); 

  /* Save node info */
  if(bsFindTag(diag,_Node))
    { bsRemove(diag); }
  bsAddTag(diag,_Node);
  e = d->elements; addpos = _bsRight;
  while(e!=NULL)
    { if(e->type==DIAG_NODE)
	{ struct node *n = e->item; int label;
	  StanAssert(bsAddData(diag,addpos,_Int,&n->id));
	  pos = bsMark(diag, pos);
	  StanAssert(bsAddData(diag,_bsRight,_Float,&n->x));
	  StanAssert(bsAddData(diag,_bsRight,_Float,&n->y));
	  if(n->label==NULL) label = -1;
	  else label = ((struct text *)n->label->item)->id;
	  StanAssert(bsAddData(diag,_bsRight,_Int,&label));
	  if(e->key!=0 && n->type != POINT_NODE) 
	    { StanAssert(bsAddKey(diag,_bsRight,e->key)); }
	                     /* above fails sometimes ***** */
/*	  if(Debug) printf("%d %5.2f %5.2f %d %s\n",n->id,n->x,n->y,label,
			   e->key!=NULL?name(e->key):""); */
	  addpos = _bsDown;
	  bsGoto(diag, pos);
	}
      e = e->next;
    }
  /* Save arc info */
  if(bsFindTag(diag,_Arc)) bsRemove(diag);
  bsAddTag(diag,_Arc);
  e = d->elements; addpos = _bsRight;
  while(e!=NULL)
    { if(e->type==DIAG_ARC)
	{ struct arc *a = e->item; int label;
	  StanAssert(bsAddData(diag,addpos,_Int,&a->id));
	  pos = bsMark(diag, pos);
	  if(a->label==NULL) label = -1;
	  else label = ((struct text *)a->label->item)->id;
	  StanAssert(bsAddData(diag,_bsRight,_Int,&label));
	  StanAssert(bsAddData(diag,_bsRight,_Int,&((struct node *)a->from->item)->id));
	  StanAssert(bsAddData(diag,_bsRight,_Int,&((struct node *)a->to->item)->id));
	  StanAssert(bsAddData(diag,_bsRight,_Int,&a->curvature));
	  StanAssert(bsAddData(diag,_bsRight,_Float,&a->radius));
	  if(a->tangent!=NULL) anytans = 1;
	  if(e->key!=0) { StanAssert(bsAddKey(diag,_bsRight,e->key)); }
/*	  if(Debug) printf("%d %d %d %d %5.2f %s\n",label,
			   ((struct node *)a->from->item)->id,
			   ((struct node *)a->to->item)->id,
			   a->curvature,a->radius,
			   e->key!=NULL?name(e->key):""); */
	  addpos = _bsDown;
	  bsGoto(diag, pos);
	}
      e = e->next;
    }

  if(anytans)
    { /* Save tangent info */
      if(bsFindTag(diag,_TangentArc)) bsRemove(diag);
      bsAddTag(diag,_TangentArc);
      e = d->elements; addpos = _bsRight;
      while(e!=NULL)
	{ if(e->type==DIAG_ARC && ((struct arc *)e->item)->tangent!=NULL)
	    { struct arc *a = e->item;
	      StanAssert(bsAddData(diag,addpos,_Int,&a->id));
	      pos = bsMark(diag, pos);
	      StanAssert(bsAddData(diag,_bsRight,_Int,&((struct arc *)a->tangent->item)->id));
/*	      if(Debug) printf("%d %d\n",a->id,((struct arc *)a->tangent->item)->id); */
	      addpos = _bsDown;
	      bsGoto(diag, pos);
	    }
	  e = e->next;
	}
    }

  /* Save label info */
  if(bsFindTag(diag,_Label)) bsRemove(diag);
  bsAddTag(diag,_Label);
  e = d->elements; addpos = _bsRight;
  while(e!=NULL)
    { if(e->type==DIAG_TEXT)
	{ struct text *t = e->item; 
	  StanAssert(bsAddData(diag,addpos,_Int,&t->id));
	  pos = bsMark(diag, pos);
	  StanAssert(bsAddData(diag,_bsRight,_Float,&t->x));
	  StanAssert(bsAddData(diag,_bsRight,_Float,&t->y));
	  StanAssert(bsAddData(diag,_bsRight,_Text,t->name));
/*	  if(Debug) printf("%d %5.2f %5.2f %s\n",t->id,t->x,t->y,t->name); */
	  addpos = _bsDown;
	  bsGoto(diag, pos);
	}
      e = e->next;
    }
     
  bsSave(diag);
  bsSave(pw);
/* diag = bsCreate(diagkey);  bsDump(diag);  bsDestroy(diag); */
msg("-SaveDiagram");
}


/* MakeDiagram: Given key and obj of a Pathway object, construct C diagram structure.
   The reactions in the Pathway are listed under the Components tag, followed
   optionally by an Enzyme tag to dismbiguate the eznyme for reactions which can
   be catalyzed by more than one enzyme. Otherwise the enzyme is gotten from the
   Enzyme tag of the Reaction object. (Reactions and Pathways are currently the same
   object, Reaction_or_Pathway). The nodes of the pathway are defined by the
   (major) reactants and products of the component Reaction_or_Pathways. For every
   reaction an arc is constructed from its reactant to its product, labelled with
   the Enzyme name or else the Reaction_or_Pathway name. If the reaction has
   multiple reactants or products, auxiliary (point) nodes are created.
   A default layout is generated for the graph.
*/
static struct diagram *MakeDiagram(KEY key)
{  OBJ pw = bsCreate(key); 
   struct diagram *d ;

   Made = 0;
   handle = handleCreate() ;
   d = (struct diagram *) metabAlloc(sizeof(struct diagram));
   d->handle = handle ;
   d->modified = 0;		/* RD was 1 - no need if remade */
msg("+MakeDiagram");
   d->key = key;
   d->name = name(key);
   d->selection = NULL;
   Elements = NULL;
   metabMakeDiagram0(key,d);
   LayoutDiagram(d);
   Made = 1;
   bsDestroy(pw);
   SortDiagram(d);
   msg("-MakeDiagram");
   return d;
 }

static BOOL IntervalsOverlap(float min1, float max1, float min2, float max2)
{ float tmp; BOOL ans;
  if(min1>max1) { tmp = min1; min1 = max1; max1 = tmp; }
  if(min2>max2) { tmp = min2; min2 = max2; max2 = tmp; }
  if(max1<min2||max2<min1) ans = FALSE;
  else ans = TRUE; 
  return ans;
}

/* fuzzy to avoid roundoff error bugs */
static BOOL InInterval(float val, float min, float max) 
{ float tmp, epsilon = 0.0001;
  if(min>max) { tmp = min; min = max; max = tmp; }
  return (val+epsilon>=min&&val-epsilon<=max); }
  
static BOOL LineSegsIntersect(float x11, float y11, float x12, float y12, 
			     float x21, float y21, float x22, float y22, 
			     float *xi, float *yi)
{ float x, y, dx1 = x12-x11, dy1 = y12-y11, dx2 = x22-x21, dy2 = y22-y21;

  if(! (IntervalsOverlap(x11,x12,x21,x22) && IntervalsOverlap(y11,y12,y21,y22)))
    return FALSE; 
  if(dx1!=0&&dx2!=0)
    { float a1, b1, a2, b2;
	a1 = dy1/dx1; b1 = y11 - a1*x11;
	a2 = dy2/dx2; b2 = y21 - a2*x21;
	if(a1==a2) 
	  { if(b1==b2) /* Parallel line segments */
		{ if(InInterval(x11,x21,x22)) *xi = x11; 
		else if(InInterval(x12,x21,x22)) *xi = x12;
		else if(InInterval(x21,x11,x12)) *xi = x21;
		else *xi = x22;
		  if(InInterval(y11,y21,y22)) *yi = y11; 
		  else if(InInterval(y12,y21,y22)) *yi = y12;
		  else if(InInterval(y21,y11,y12)) *yi = y21;
		  else *yi = y22;
		  return TRUE; }
	    else return FALSE; }
	x = (b2-b1)/(a1-a2); y = a1*x + b1;
    }
  else if(dx1!=0)
    { float a1 = dy1/dx1, b1 = y11 - a1*x11;
	x = x21; y = a1*x + b1; }
  else if(dx2!=0)
    { float a2 = dy2/dx2, b2 = y21 - a2*x21;
	x = x11; y = a2*x + b2; }
  else { /* Overlapping verticals */
    *xi = x11; 
    if(InInterval(y11,y21,y22)) *yi = y11; 
    else if(InInterval(y12,y21,y22)) *yi = y12;
    else if(InInterval(y21,y11,y12)) *yi = y21;
    else *yi = y22;
    return TRUE; }
  /* Clipping: is (x,y) withing segs */
  if(InInterval(x,x11,x12) && InInterval(x,x21,x22)
     && InInterval(y,y11,y12) && InInterval(y,y21,y22))
    { *xi = x; *yi = y; return TRUE; }
  else return FALSE;
}

/*
static void spaces(int i){ int j;   for(j=0;j<i;j++) printf(" "); }
*/

static void metabMakeDiagram0(KEY key, struct diagram *d)
{  KEY rxnkey, mbkey; 
   OBJ rxn, pw = bsCreate(key);  
   BSMARK pos = NULL, pos1 = NULL;
   struct element  *e;
   int i, id = 1, j, done, done1 ;

   /* Make nodes */
   StanAssert(bsGetKey(pw,_Components,&rxnkey));
   done = 0;
   while(!done)
     { pos = bsMark(pw, pos);  
       StanAssert(bsGetKey(pw,_bsHere,&rxnkey));
       rxn = bsCreate(rxnkey);
       if(rxn!=NULL)
	 for(j=1;j<=2;j++)
	   { KEY tag = _Major_Reactant;
	     if(j==2) tag = _Major_Product;
	     if(bsGetKey(rxn,tag,&mbkey))
	       { done1 = 0;
		 while(! done1)
		   { pos1 = bsMark(rxn,pos1);
		     if(Key2Node(mbkey)==NULL)
		       { struct node *n;
		       e = MakeNode(id,NULL,0,0); id++ ;
			 e->key = mbkey;
			 n = e->item;
			 n->label = MakeText(id++,name(mbkey),0,0,0,NodeColor);
			 n->label->key = mbkey; }
		     /* bsGoto(rxn,pos1) -- not needed yet, but will be */
		     done1 = ! bsGetKeyTags(rxn,_bsDown,&mbkey); }}}
       bsDestroy(rxn);
       /* bsGoto(pw,pos) -- not needed yet, but will be */	 
       done = ! bsGetKeyTags(pw,_bsDown,&rxnkey); 
     }
   /* Make Arcs */
   StanAssert(bsGetKey(pw,_Components,&rxnkey));
   done = 0;
   while(!done)
     { pos = bsMark(pw, pos);  
       StanAssert(bsGetKey(pw,_bsHere,&rxnkey));
       { /* Skip reactions for which arcs already exist --- this
	    can happen during expansion. */
	 int exists = 0;
	 e = Elements; 
	 while(e!=NULL)
	   { if(e->type==DIAG_ARC && e->key == rxnkey) { exists = 1; break; }
	     e = e->next; }
	 if (exists) 
	   { done = ! bsGetKeyTags(pw,_bsDown,&rxnkey);
	     continue; }}
	 rxn = bsCreate(rxnkey);
	 if(rxn!=NULL) /* Foreach reaction in components of pathway: */
	   { int nr = bsNRows(rxn,_Major_Reactant), np = bsNRows(rxn,_Major_Product);
	     KEY r, p, enz; 
	     struct element *e, *from, *to, *lab, *arc; char *label;
	     if(!(nr&&np))
		 { char buf[300];
		   sprintf(buf,"No reactants/products listed for %s;\n\
It will not be included in the diagram.",name(rxnkey));
		   msg(buf);
		   graphOut(buf);
		   done = ! bsGetKeyTags(pw,_bsDown,&rxnkey); 
		   continue; }
	     if(nr==1)
		 { StanAssert(bsGetKey(rxn,_Major_Reactant,&r));
		   from = notNULL(Key2Node(r)); }
	     else { from = MakePoint(id++); from->key = rxnkey; }
	     if(np==1)
		 { StanAssert(bsGetKey(rxn,_Major_Product,&p));
		   to = notNULL(Key2Node(p)); }
	     else { to = MakePoint(id++); to->key = rxnkey; }
	     /* Get arc label: Enzyme name or Rxn/PW name.
		  should check for enzyme(s) at level of Component;
		  also handle multiple enzymes at either place. */
	     if(bsGetKey(rxn,_Enzyme,&enz)) label = name(enz);
	     else label = name(rxnkey);
	     
	     arc = MakeArc(id,from,to,lab = MakeText(id+1,label,0,0,1,ArcColor),0);
	     id += 2;
	     arc->key = rxnkey; lab->key = rxnkey;

	     /* Generate tangent arc for cofactor reactions */
	     for(i=1;i<=2;i++)
	       { KEY reactant = _Cofactor_Reactant, product = _Cofactor_Product;
		 int curvature = -1;
		 if(i==2)
		   { reactant = _Additional_Cofactor_Reactant;
		     product = _Additional_Cofactor_Product;
		     curvature = 1; }
		 if(bsFindTag(rxn,reactant))
		   { /* Assumes only 1 reactant and product per tangent arc */
		     KEY fromkey, tokey; struct element *tan; struct arc *tanarc;
		     struct element *from, *to; struct node *fromn, *ton;
		     fromkey = tokey = 0;
		     bsGetKey(rxn,_bsRight,&fromkey);
		     bsGetKey(rxn,product,&tokey);
		     from = MakeNode(id,MakeText(id+1,fromkey==0?"":name(fromkey),
						   0,0,0,CofColor),0,0);
		     id += 2;
		     to = MakeNode(id,MakeText(id+1,tokey==0?"":name(tokey),
						 0,0,0,CofColor),
				   0,0);
		     id += 2;
		     tan = MakeArc(id++,from,to,NULL,curvature);
		     tan->color = from->color = to->color = CofColor;
		     tanarc = tan->item; fromn = from->item; ton = to->item;
		     fromn->type = ton->type = TANARC_NODE;
		     ((struct arc *) arc->item)->tangent = tan;
		     tanarc->tangent_to = arc;
		     tanarc->radius = 4;
		     tanarc->from->key = fromkey;
		     tanarc->to->key = tokey;
		     fromn->label->key = fromkey; ton->label->key = tokey;
		   }}

/*	     if(np>1) ((struct arc *) arc->item)->arrow = 0; */
	     
	     if(nr>1) /* Create half-arcs for multiple reactants */
		 { int done = 0;
		   BSMARK pos = NULL;

		   StanAssert(bsGetKey(rxn,_Major_Reactant,&r));
		   while(! done)
		     { struct element *ea; 
		       pos = bsMark(rxn, pos);  
		       StanAssert(bsGetKey(rxn,_bsHere,&r));
		       e = notNULL(Key2Node(r));
		       ea =  MakeArc(id++,e,from,NULL,0);
		       /* 
			  ea->item->arrow = 0;
			  */
		       ea->key = rxnkey;
		       done = !bsGetKey(rxn,_bsDown,&r); }}
	     
	     if(np>1) /* Create half-arcs for multiple products */
	       { int done = 0;
		 BSMARK pos = NULL;

		 StanAssert(bsGetKey(rxn,_Major_Product,&p));
		 while(! done)
		   { struct element *ea;
		     pos = bsMark(rxn, pos);  
		     StanAssert(bsGetKey(rxn,_bsHere,&p));
		     e = notNULL(Key2Node(p));
		     ea = MakeArc(id++,to,e,NULL,0);
		     ea->key = rxnkey;
		     done = !bsGetKey(rxn,_bsDown,&p); }}
	     bsDestroy(rxn);
	   }
	 /* bsGoto(pw,pos) -- not needed yet, but will be */	 
	 done = ! bsGetKeyTags(pw,_bsDown,&rxnkey); }
   
   d->elements = Elements;
 }

/* ExpandDiagram: Given key of a Subpathway, restructure the Active Diagram to
   expand the subpwathway.
*/
static void ExpandDiagram(KEY key)
{  OBJ  pw = bsCreate(key); 
   struct diagram *d = ActiveDiagram();
   struct element **p, *e;
   int    done ;

   d->modified = 0;		/* RD - was 1 */
   printf("+ExpandDiagram expanding %s in %s\n",name(key),d->name);

   /* Verify that exansion is appropriate! */
   done = 0;
   e = d->elements; 
   while(e!=NULL) { if(e->key == key) { done = 1; break; } 
		    e = e->next; }
   if(! done) return;
   
   /* Eliminate unexpanded elements from diagram */
   e = d->elements; while(e!=NULL) { e->mark = 0; e = e->next; }
   e = d->elements; 
   done = 1;
   /* Mark expanded arc, its point-nodes if any, its tangent and its nodes if any */
   while(e!=NULL) 
     { if(e->type == DIAG_ARC && e->key == key)
	 { struct arc *a = e->item;
	   struct node *from, *to;
	   e->mark = 1; 
	   if(a->label!=NULL) a->label->mark = 1; 
	   from = a->from->item; to = a->to->item;
	   if(from->type == POINT_NODE) { a->from->mark = 1; done = 0; }
	   if(to->type == POINT_NODE) {  a->to->mark = 1; done = 0; }
	   if(a->tangent!=NULL)
	     { struct arc *tan = a->tangent->item;
	       a->tangent->mark = 1;
	       tan->from->mark = 1;
	       from = tan->from->item; 
	       if(from==NULL) msg("null from");
	       else if(from->label!=NULL) from->label->mark = 1;
	       tan->to->mark = 1;
	       to = tan->to->item; if(to->label!=NULL) to->label->mark = 1; }
	   break; }
       e = e->next; }
   /* Remove dangling arcs */
   if(!done)
     { e = d->elements;
       while(e!=NULL)
	 { if( (!e->mark) && e->type == DIAG_ARC )
	     { struct arc *a = e->item;
	       if(a->from->mark || a->to->mark)
		 { e->mark = 1;
		   if(a->label!=NULL) a->label->mark = 1; }}
	   e = e->next; } }
       
   /* Splice out marked diagram elements */
   p = &d->elements;
   while(*p!=NULL)
     { e = *p;
       if(e->mark) *p = e->next;
       else p =  &e->next; }
   Elements = d->elements;
   metabMakeDiagram0(key,d);
   LayoutDiagram(d);
   Made = 1;
   bsDestroy(pw);
   graphClear();		/* RD - was graphDestroy */
   DrawDiagram(d);
   msg("-ExpandDiagram");
 }

/********************* GRAPH LAYOUT SECTION *******************/
/* RD removed embedded function definitions to the head, because some
   compilers (e.g. OSF cc) don't accept them
*/

struct nlist { struct node *node; struct nlist *next; } ;

struct cycle { int id, length, mark;
	       enum cycle_type { MAIN, SHUNT, TANGENT } type;
	       struct nlist *nodes;
	       float phase;
	       struct cycle *class, *main, *next; } *cycles ;

static int cycid, ncycles ;
static int s, n, *g, *hord ;
static struct node **node, **stack ;

static void cycle_find (int i)
  /* 
	initialize marks to unseen
	call cycle_break from each root, or else from node[0] if no roots.

	Cycle-break does a depth-first traversal of g starting from i.
	Nodes are marked unseen (0), in-progress (1) or seen (2).
	A seen node is not part of any cycle.
	A cycle occurs when an in-progress node is encountered at
	the end of an arc. The arc is then eliminated from g.
	An adjacent cycle may also occur if a seen node is encountered that
	is part of a cycle any of whose elements are in progress.

	Cycle-break returns 1 if the node it was called on was not in-progress. 

	Keep track of cycles. Create a new cycle structure for each.
*/
{ int j,k;
  
  if(!node[i]->mark)
    { node[i]->mark = 1;	/* in progress */
      stack[s++] = node[i];
      k=i*n;
      for(j=0;j<n;j++) { if(g[k]) cycle_find(j); k++; }
      node[i]->mark = 2;	/* seen */
      s--; }
  else if(node[i]->mark == 1)	/* cycle detected, create structure */
    { struct cycle *new = (struct cycle *) metabAlloc(sizeof(struct cycle));
      struct nlist **p = &new->nodes;
      int k, in = 0, l = 0;
      new->id = cycid++; ++ncycles; new->class = new; new->type = MAIN;
      if(Debug) printf("Cycle#%d: ",new->id);
      for(k=0;k<s;k++)
	{ if(stack[k]==node[i] ) in = 1;
	  if(in)
	    { l++; 
	      stack[k]->cycle = 1;
	      (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
	      (*p)->node = stack[k];
	      (*p)->next = NULL;
	      if(Debug)printf("%s ",ElementLabel(stack[k]->label));
	      p = &(*p)->next; } }
      if(Debug)printf("\n");
      (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
      (*p)->node = node[i];
      (*p)->next = NULL;
      new->next = cycles;
      new->length = l; 
      cycles = new; }
  else if(node[i]->cycle) 
    {			/* Possible adjacent cycle:
			   find in_progress cycle that node[i] is on */
      struct cycle *c = cycles;
      int inprog = 0, containsi = 0, inprognode = 0;
      while(c!=NULL)
	{ struct nlist *l = c->nodes;
	  inprog = containsi = 0;
	  while(l!=NULL)
	    { if(l->node->i == i) containsi = 1;
	      if((!inprog) && l->node->mark == 1)
		{ inprog = 1; 
		  inprognode = l->node->i; }
	      if(containsi  && inprog) break;
	      l = l->next; }
	  if(containsi  && inprog) break;
	  c = c->next;
	}
      if(c!=NULL)
	{ struct cycle *new = (struct cycle *) metabAlloc(sizeof(struct cycle));
	  struct nlist **p = &new->nodes, *l;
	  int k, in = 0, len = 0;
	  if(Debug)printf("Adjacent cycle involving %s and %s in cycle#%d\n",
			  ElementLabel(node[i]->label),
			  ElementLabel(node[inprognode]->label),
			  c->id);
	  new->id = cycid++; ++ncycles; new->class = new; new->type = MAIN;
	  if(Debug)printf("Cycle#%d: ",new->id);
	  for(k=0;k<s;k++) 
	    { if(stack[k]->i == inprognode ) in = 1;
	      if(in)
		{ len++; 
		  stack[k]->cycle = 1;
		  (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
		  (*p)->node = stack[k];
		  (*p)->next = NULL;
		  if(Debug)printf("%s ",ElementLabel(stack[k]->label));
		  p = &(*p)->next; } }
	  l = c->nodes;
	  in = 0;
	  if(Debug)printf(" / ");
	  while(1)
	    { if(l->node->i == i) in = 1;
	      if(in)
		{ len++; 
		  (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
		  (*p)->node = l->node;
		  (*p)->next = NULL;
		  if(Debug)printf("%s ",ElementLabel(l->node->label));
		  p = &(*p)->next; 
		  if(l->node->i == inprognode) break; }
	      l = l->next;
	      if(l==NULL) l = c->nodes; }
	  if(Debug)printf("\n");
	  new->next = cycles;
	  new->length = len; 
	  cycles = new; 
	}
    }
}

static void phase(struct cycle *c)
      /* Phase = 0 implies first node in list is at top of cycle.
	 Compute optimal shift as the average of +shift to input nodes
	 and -shift to output nodes. */
{ float p = Pi/2, a = 0, y, da = 2*Pi/c->length;
  struct node *prev;
  int nl = 0;
  struct nlist *l = c->nodes;
  msg("Cycle phasing");
  while(l->next!=NULL)
    { struct node *n = l->node;
      p += (n->indeg-1)*a; 
      p += (n->outdeg-1)*(a+Pi);
      nl += n->indeg + n->outdeg - 2;
      a += da; 
      l = l->next; }
  if(nl) p = p/nl;
  c->phase = p;
  /* Invert upward-pointing arcs. */
  y = -sin(p);
  l = c->nodes->next;
  prev = c->nodes->node;
  while(l!=NULL)
    { struct node *cur = l->node;
      float y1;
	    p += da;
	    y1 = -sin(p);
	    if(y1<y) /* Invert arc */
	      { int k = prev->i*n+cur->i;
		if(g[k])
		  { g[k] = 0; g[cur->i*n+prev->i] = 1; 
		    if(Debug)printf("Inverting %s -> %s\n",
				    ElementLabel(prev->label),
				    ElementLabel(cur->label));
		  }
	      }
      y = y1; 
      prev = cur;
      l = l->next;
    }
}

 /* Cycle-breaking:
    Double-check to make sure arc-inversion has not introduced any new cycles. 
    If so, warn, break, and sauve qui peut.
    */
static int cycle_break (int i)
{ int j,k;
  if(!node[i]->mark)
    { node[i]->mark = 1; /* in progress */
      k=i*n;
      for(j=0;j<n;j++) 
	{ if(g[k]) 
	    { g[k] = cycle_break(j); 
	      if(0&&!g[k]) 
		graphMessage("Arc inversion introduced new cycles.  "
			     "Show this diagram to Stan."); 
	    }
	  k++; 
	}
      node[i]->mark = 2; /* seen */
      return 1; 
    }
  else if(node[i]->mark == 1) /* cycle detected */
    return 0;
  else return 1; 
}

static int maxd ;

static void set_depth (int i)
{ int j,k = i*n;
  for(j=0;j<n;j++)
    { if(g[k] && node[j]->depth < node[i]->depth+1)
	{ node[j]->depth = node[i]->depth+1;
	  if(node[j]->depth > maxd) maxd = node[j]->depth;
	  set_depth(j); 
	}
      k++; 
    }
}

static struct cycle *ID2cycle(int i)
{ struct cycle *c = cycles;
  while(c!=NULL)
    { if(c->id == i) return c;
      c = c->next; 
    }
  StanAssert(0); 
  return c ;
}

static int max_cycle_length (int i, int class)
{ int maxlen = 0,  j,k;
  /* spaces(s); printf("MCL %s\n",ElementLabel(node[i]->label)); */
  if(!node[i]->mark)
    { node[i]->mark = 1; /* in progress */
      stack[s++] = node[i];
      k=i*n;
      for(j=0;j<n;j++)
	{ if(hord[k]==class) 
	    { int len = max_cycle_length(j,class); 
	      if(len>maxlen) maxlen = len; }
	  k++;}
      node[i]->mark = 0; 
      s--; 
      return maxlen;
    }
  else if(node[i]->mark) /* cycle detected, return length */
    { int len = 0, in = 0;
      for(j=0;j<s;j++)
	{ if(stack[j]==node[i]) in = 1;
	  if(in) len++; }
      return len; 
    }
  return 0 ;
}

static struct cycle *longest_cycle(int i, int class, int len)
{ int j,k;
  struct cycle *c = NULL;
  /* spaces(s); printf("LC %s\n",ElementLabel(node[i]->label)); */
  if(!node[i]->mark)
    { node[i]->mark = 1; /* in progress */
      stack[s++] = node[i];
      k=i*n;
      for(j=0;j<n;j++)
	{ if(hord[k]==class) 
	    { c = longest_cycle(j,class,len);
	      if(c!=NULL) break; }
	  k++;}
      node[i]->mark = 0; 
      s--; 
      return c; 
    }
  else if(node[i]->mark == 1) 
    { int l = 0, in = 0;
      for(j=0;j<s;j++)
	{ if(stack[j]==node[i]) in = 1;
	  if(in) l++; 
	}
      if(l==len)
	{ struct cycle *new = (struct cycle *) metabAlloc(sizeof(struct cycle));
	  struct nlist **p = &new->nodes;
	  int k, in = 0;
	  new->id = cycid++; ++ncycles; new->class = new; 
	  new->type = MAIN;
	  if(Debug)printf("Cycle#%d: ",new->id);
	  for(k=0;k<s;k++)
	    { if(stack[k]==node[i] ) in = 1;
	      if(in)
		{ stack[k]->cycle = 1;
		  (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
		  (*p)->node = stack[k];
		  (*p)->next = NULL;
		  if(Debug)printf("%s ",ElementLabel(stack[k]->label));
		  p = &(*p)->next; 
		} 
	    }
	  if(Debug)printf("\n");
	  (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
	  (*p)->node = node[i];
	  (*p)->next = NULL;
	  new->next = cycles;
	  new->length = len; 
	  cycles = new; 
	  return new; 
	}
      return NULL; 
    }
  return NULL ;
}

static void make_shunt(int j, int k, struct cycle *maxcycle, struct cycle *eqclass)
{ struct cycle *new = (struct cycle *) metabAlloc(sizeof(struct cycle));
  struct nlist **p = &new->nodes;
  int j1 = j, k1 = k, len = 0;
  if(Debug)printf("Shunt: ");
  new->id = cycid++; 
  ++ncycles; new->class = new; new->type = SHUNT; new->main = maxcycle;
  while(1)
    { (*p) = (struct nlist *) metabAlloc(sizeof(struct nlist));
      (*p)->node = node[j1];
      (*p)->next = NULL; len++;
      if(! node[j1]->mark) node[j1]->mark = 2; /* mark for recheck */
      hord[j1*n+k1] *= -1;
      if(Debug)printf("%s ",ElementLabel(node[j1]->label));
      p = &(*p)->next; 
      if(k1==-1) break;
      j1 = k1;
      if(node[k1]->mark) k1 = -1;
      else { for(k1=0;k1<n;k1++)
	       if(hord[j1*n+k1]==eqclass->id) break;
	     if(k1==n) k1  = -1; }}
  if(Debug)printf("\n");
  /* Add to *end* of cycle list -- ensures main cycle is
     processed first, so on-cycle nodes have coords
     before shunt coords need to be computed. */
  { struct cycle **c = &cycles;
    while((*c)!=NULL) c = &(*c)->next;
    (*c) = new;
    new->length = len; 
    new->next = NULL; 
  }
}

static float chord_angle(struct node *a, struct node *b, struct cycle *c)
{ int nonchord = 0, len = 0; int in = 0;
  struct nlist *l = c->nodes;
  while(!(in && l->node == b))
    { if(Debug)printf("%s%s ",in?"+":"-",ElementLabel(l->node->label));
      if((!in) && l->node == a) in = 1;
      if(in) len++;
      if(len>c->length) { nonchord = 1; break; }
      l = l->next;
      if(l->node==c->nodes->node) 
	{ if(in) l = c->nodes;
	  else { nonchord = 1; break; }
	}
    }
  if(nonchord) return 0;
  return len; 
}

/*
   Graph layout alg 
   
   Uses a hierachical layout algorithm which performs reasonably on DAGS,
   augmented with cycle handling features. Attempts to minimize arc crossings,
   lay out cyclic pathways as circles, and run metabolic flow from top to bottom.
   
   The hierarchical algorithm works as follows:
   *break cycles by deleting or reversing arcs
   *compute max depth of each node
   *for each depth, order nodes consistent with constraints
   inherited from previous rows. For every ordering decision
   i < j, require that children of i < children j (except for
   shared children). These become constraints on ordering
   the children.

A1: Assume graph G is a planarizable DAG.

	Find ordering on children of each node that
	achieves near-planarity.

T1: (non)Theorem:
	let B & C be children of A.
	B < C is consistent with planarity iff
		tc(B) - tc(C)  < tc(B) intersect tc(C) < tc(C) - tc(B)
		is consistent w/ planarity

Counterexample: if A & B are at depth i, and C & D are at i+1,
	and A and B each have arcs to C & D, then some arcs will cross.
	However, it does a good job much of the time.

Representation is a consistent horizontal partial order which
is incrementally constrained. (DAG arrow direction (= metabolic flow)
is layed out from top to bottom, so horizontal po means ordering at 
right angles to DAG.)

Depth(0): order arbitrarily
Depth(D+1): assuming po at D was consistent and A and B are unordered,
then if A<B is consistent by T1, then do it, else do B<A.

ALG:
find cycles, choose phase, reverse upward arcs
find max depth of all nodes
for D=0 to max max depth
	while exists unordered A,B at depth D do
		assume order A<B
		tentative transitive closure (erasable -- use a different marker)
		if contradiction do B<A
	end
now have consistent horizontal partial order, which is total within each depth row.
Assign coords consistent w/ both partial orders.

Horizontal assignment now works by equispacing widest row, propagating
x-vals up and down to relatives. (Makes above pointless?)

*/

static void LayoutDiagram(struct diagram *d)
{ struct element *e; 
  int nsq, done; 
  /* hord[i,j] = hord[i*n+j] = 1 iff node[i] < node[j], 
     -1 iff node[j] < node[i],
     0 if unordered
     vord is same, describes transitive closure of flow order (downstreamness)
     ie vord[i,j] implies an arc or chain of arcs from i to j.
     hord describes left-right ordering 
     node[i] is the node with ->i = i
     nnodes[depth] is # of nodes at that depth
     the row[i] are the nodes at the current depth, in no particular order.
     */
  int *vord, *nnodes, i, j, depth ;
  float rowheight;
  struct node **row;
  struct nlist **layout;

  cycid = 1 ; ncycles = 0 ; cycles = NULL ; 
  n = 0 ; maxd = 0 ;

  msg("+Layout");
  /* n = # nodes */
  e = d->elements; 
  while(e!=NULL)
    { if(e->type==DIAG_NODE && ((struct node *)e->item)->type !=TANARC_NODE)  n++;
      e = e->next; }
  nsq = n*n;

  /* Allocate structures */
  hord = (int*) metabAlloc(nsq*sizeof(int));
  vord = (int*) metabAlloc(nsq*sizeof(int));
  g = (int*) metabAlloc(nsq*sizeof(int));
  node = (struct node **) metabAlloc(n*sizeof(struct node *));
  row = (struct node **) metabAlloc(n*sizeof(struct node *));
  stack = (struct node **) metabAlloc(n*sizeof(struct node *));

  /* Initialize order arrays to unordered. */
  for(i=0;i<nsq;i++) g[i] = vord[i] = hord[i] = 0;
  
  /* Set up array of nodes. node[i]; node[i]->i = i */
  e = d->elements; i = 0;
  while(e!=NULL)
    { if(e->type==DIAG_NODE&&((struct node *)e->item)->type!=TANARC_NODE) 
	{ struct node *nd = e->item;
	  nd->i = i;
	  node[i++] = nd; }
      e = e->next; }
  StanAssert(i==n);

  /* Initialize graph (g) based on flow DAG */
  e = d->elements; 
  while(e!=NULL)
    { if(e->type==DIAG_ARC&&((struct arc *)e->item)->tangent_to==NULL) 
	{ struct arc *a = e->item;
	  struct node *from = a->from->item, *to = a->to->item;
	  int i = from->i, j = to->i;
 	  g[i*n+j] = 1; }
      e = e->next; }

  /* Compute node degrees */
  for(i=0;i<n;i++){ node[i]->depth = node[i]->indeg = node[i]->outdeg = 0; }
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if(g[i*n+j]){ node[i]->outdeg++; node[j]->indeg++; }

  msg("Cycle finding");
  for(i=0;i<n;i++) { node[i]->cycle = 0; node[i]->mark = 0; /* unseen */ }
  s = 0 ;
  for(i=0;i<n;i++) if(! node[i]->mark) cycle_find (i); 

  if(ncycles>1)
    { /* Group cycles sharing arcs into equivalence classes.
	 Restructure to pull out longest cycle and shunts. 
	 Reuse hord array for each equivalence class.
	 */
      struct cycle *c = cycles;
      struct cycle *eqclass[10];
      int nclasses = 0;
      msg("=Class detection");
      while(c!=NULL)
	{ struct nlist *l = c->nodes->next;
	  struct node *prev = c->nodes->node;
	  if(Debug)printf("Cycle#%d: %s",c->id,ElementLabel(prev->label));
	  while(l!=NULL)
	    { int k = prev->i*n+l->node->i;
	      if(Debug)printf(" -> %s",ElementLabel(l->node->label));
	      if(!hord[k]) hord[k] = c->class->id;
	      else { struct cycle *a = ID2cycle(hord[k])->class, *b = c->class, *c;
		     if(a!=b)
		       { /* Merge classes */
			 if(a->id > b->id) { c = a; a = b; b = c; }
			 if(Debug)printf("Merging class %d into class %d.\n",
				b->id, a->id);
			 c = cycles;
			 while(c!=NULL)
			   { if(c->class == b) c->class = a;
			     c = c->next; }}}
	      prev = l->node;
	      l = l->next; }
	  if(Debug)printf("\n");
	  c = c->next; }
      c = cycles; while(c!=NULL) { c->mark = 0; c = c->next; }
      c = cycles; while(c!=NULL) { if(c->class!=c) c->class->mark = 1; c = c->next; }
      nclasses = 0;
      c = cycles; while(c!=NULL) { if(c->mark) eqclass[nclasses++] = c; c = c->next; }
      if(Debug)printf("Found %d =classes\n",nclasses);
      if(nclasses)
	{ /* canonicalize hord */
	  c = cycles;
	  while(c!=NULL)
	    { if(c->class != c)
		{ struct nlist *l = c->nodes->next;
		  struct node *prev = c->nodes->node;
		  while(l!=NULL)
		    { int k = prev->i*n+l->node->i;
		      hord[k] = c->class->id;
		      prev = l->node;
		      l = l->next; }}
	      c = c->next; }
	  /* find longest cycles & shunts for merged classes */
	  for(i=0;i<nclasses;i++)
	    { int k, maxlen;
	      struct cycle *maxcycle;
	      /* Delete cycles in eqclass */
	      { struct cycle **c = &cycles;
		while((*c)!=NULL)
		  { if((*c)->class == eqclass[i]) (*c) = (*c)->next;
		    else c = &(*c)->next; }}
	      /* Find length of longest cycle made of of class arcs */
	      for(k=0;k<n;k++) 
		node[k]->mark = 0; 
	      s = 0;
	      maxlen = max_cycle_length(eqclass[i]->nodes->node->i, 
					eqclass[i]->class->id);
	      if(Debug)printf("Longest cycle is %d\n",maxlen);
	      /* Construct a cycle of that length. */
	      for(k=0;k<n;k++) 
		node[k]->mark = 0;
	      s = 0;
	      maxcycle = longest_cycle(eqclass[i]->nodes->node->i,
				       eqclass[i]->class->id, maxlen);
	      /* Create shunts. */
	      { struct nlist *l = maxcycle->nodes->next;
		struct node *prev = maxcycle->nodes->node;
		/* Mark nodes in longest cycle. */
		for(k=0;k<n;k++) node[k]->mark = 0;
		/* Mark arcs in longest cycle by negating */
		while(l!=NULL) 
		  { hord[prev->i*n+l->node->i] *= -1;
		    prev = l->node;
		    l->node->mark = 1; 
		    l = l->next; }
		/* Find arcs in eqclass emanating from longest cycle but not on it. */
		for(j=0;j<n;j++)
		  if(node[j]->mark) /* if node on longest cycle or already seen */
		    for(k=0;k<n;k++)
		      if(hord[j*n+k] == eqclass[i]->id) /* arc not on cycle */
			{ node[j]->mark = 1; make_shunt(j,k,maxcycle,eqclass[i]); }
		{ int done = 0;
		  while(!done)
		    { done = 1;
		      for(j=0;j<n;j++)
			if(node[j]->mark == 2) /* if node on shunt and not yet checked */
			  { for(k=0;k<n;k++)
			      if(hord[j*n+k] == eqclass[i]->id) /* arc not on cycle */
				{ node[j]->mark = 1; make_shunt(j,k,maxcycle,eqclass[i]); done = 0;  }
			  }} }
		    } 
	    } /* fi nclasses */
	} /* fi canonicalize hord */
    } /* fi ncycles > 1 */

      /* Choose cycle phases, invert arcs. */
  { struct cycle *c = cycles; 
    while(c!=NULL) { if(c->type == MAIN) phase(c); c = c->next; } 
  }

  msg("Cycle Breaking");
  { for(i=0;i<n;i++) node[i]->mark = 0; /* unseen */
    for(i=0;i<n;i++) if(!node[i]->mark) cycle_break(i);  
  }

  /* Reompute node degrees */
  for(i=0;i<n;i++){ node[i]->depth = node[i]->indeg = node[i]->outdeg = 0; }
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if(g[i*n+j]){ node[i]->outdeg++; node[j]->indeg++; }

  /* Set depth of each node to its max depth -- safe now because graph is acyclic. 
     maxd = maximum depth. */
  { int anyroots = 0;
    for(i=0;i<n;i++) if(!node[i]->indeg) { anyroots = 1; set_depth(i); }
    if(! anyroots)
      { struct cycle *c = cycles;
	for(i=0;i<n;i++) node[i]->depth = -1;
	while(c!=NULL) 
	  { if(c->nodes->node->depth!=-1)
	      set_depth(c->nodes->node->i); 
	    c = c->next;
	  }}
  }

  for(i=0;i<nsq;i++) hord[i] = 0;   /* Reinit hord */

  for(i=0;i<n;i++)
    /* Correct root depths to be min(child depth) - 1 */
    if(!node[i]->indeg)
      { int k=i*n, minchild = maxd+1;
	for(j=0;j<n;j++)
	  if(g[k++] && node[j]->depth < minchild)
	    minchild = node[j]->depth;
	node[i]->depth = minchild-1; }
    else if(!node[i]->outdeg) /* align sinks at maxd */
      node[i]->depth = maxd; 
  layout = (struct nlist **) metabAlloc((maxd+1)*sizeof(struct nlist *));
  nnodes = (int*) metabAlloc((maxd+1)*sizeof(int));

  /* copy g to vord; allow both negative ane positive links */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if(g[i*n+j]) { vord[i*n+j] = 1; vord[j*n+i] = -1; }

  /* vord = Transitive closure of g */
  done = 0;
  while(! done)
    { int k,l,m;
      done = 1;
      for(k=0;k<n;k++)
	for(l=0;l<n;l++)
	  if(k!=l)
	    {
	      for(m=0;m<n;m++)
		if(m!=l&&m!=k&&vord[k*n+l]*vord[l*n+m]>0)
		  {
		    if(!vord[k*n+m]) 
		      {
			int mark = vord[k*n+l];
			done = 0;
			vord[k*n+m] = mark; vord[m*n+k] = -mark;
		      }
		    else if(vord[k*n+m]!=vord[k*n+l]) 
		      { StanAssert(0); }
		  }
	    }
    }
  
  /* Horizontally order nodes at each depth */
msg("for depth");
  for(depth=0;depth<=maxd;depth++)
    { nnodes[depth] = 0;
      /* row = list of nodes at this depth, # in nnodes[depth] */
      for(i=0;i<n;i++) if(node[i]->depth == depth) row[nnodes[depth]++] = node[i];
      if(Debug)printf("%d nodes at depth %d\n",nnodes[depth],depth);

      /* For each pair of nodes i and j */
      for(i=0;i<nnodes[depth]-1;i++)
	for(j=i+1;j<nnodes[depth];j++)
	  if(!hord[row[i]->i*n+row[j]->i])
	    { int k,l,m;
	      /* Tentatively assume order i < j -- use +/-2 as tentative mark */
	      hord[row[i]->i*n+row[j]->i] = 2; hord[row[j]->i*n+row[i]->i] = -2;
	      if(Debug)printf("Assuming %s < %s\n",ElementLabel(row[i]->label),
		     ElementLabel(row[j]->label));
msg("order kids");
	      /* order vertical children consistently with horizontal order */
	      done = 0;
	      for(k=0;k<n;k++)
		for(l=0;l<n;l++)
		  if(k!=l&& vord[row[i]->i*n+k]&&vord[row[j]->i*n+l]
		     && (! vord[row[j]->i*n+k]) && (!vord[row[i]->i*n+l]))
		    { int h = k*n+l;
		      if(hord[h]<0) done = -1;
		      if(! hord[h]) 
			{ hord[h] = 2; hord[l*n+k] = -2; 
			  if(Debug)printf("Implies %s < %s\n",
				 node[k]->label==NULL?"No label":
				 ((struct text *) node[k]->label->item)->name,
				 node[l]->label==NULL?"No label":
				 ((struct text *) node[l]->label->item)->name); }}
msg("tc ord");
	      /* Transitive closure of horizontal order */
	      while(! done)
		{ done = 1;
		  for(k=0;k<n;k++)
		    for(l=0;l<n;l++)
		      if(k!=l)
			for(m=0;m<n;m++)
			  if(m!=l&&m!=k&&hord[k*n+l]*hord[l*n+m]>0)
			    { if(!hord[k*n+m]) 
				{ done = 0;
				  hord[k*n+m] = hord[k*n+l]>0?2:-2;
				  hord[m*n+k] = -hord[k*n+m]; }
			    else if(hord[k*n+l]*hord[k*n+m]<0)
			      { done = -1;
				break; break; break; }}}

	      if(done>0) /* Accept order i < j */
		{ for(k=0;k<nsq;k++)
		    if(hord[k]==2) hord[k] = 1;
		    else if(hord[k]==-2) hord[k] = -1; }
	      else 
		{ /* Use order j < i */
		  msg("Rejecting order");
		  for(k=0;k<nsq;k++) if(hord[k]==2||hord[k]==-2) hord[k] = 0;
		  hord[row[j]->i*n+row[i]->i] = 1; hord[row[i]->i*n+row[j]->i] = -1;

		  /* order vertical children consistently with horizontal order */
		  for(k=0;k<n;k++)
		    for(l=0;l<n;l++)
		      if(vord[row[j]->i*n+k]&&vord[row[i]->i*n+l]
			 && (! vord[row[i]->i*n+k]) && (!vord[row[j]->i*n+l]))
			{ int h = k*n+l;
			  if(!hord[h]) { hord[h] = 1; hord[l*n+k] = -1; }
			  else if(hord[h]==-1) 
			    if(Debug)printf("B Contradiction placing %s < %s\n",
				   ElementLabel(node[k]->label),
				   ElementLabel(node[l]->label)); }
msg("tc hord1");
		  /* Transitive closure of horizontal order */
		  done = 0;
		  while(! done)
		    { done = 1;
		      for(k=0;k<n;k++)
			for(l=0;l<n;l++)
			  if(k!=l)
			    for(m=0;m<n;m++)
			      if(m!=l&&m!=k&&hord[k*n+l]*hord[l*n+m]>0)
				{ if(!hord[k*n+m]) 
				    { done = 0;
				      hord[k*n+m] = hord[k*n+l]>0?1:-1;
				      hord[m*n+k] = -hord[k*n+m]; }
				else if(hord[k*n+m]*hord[k*n+l]<0)
				  msg("C Contradiction");  }}}
	    } /* end foreach pair i,j in row */
msg("End ij");

      /* Sort row by hord */
      { int done = 0;
	while(! done)
	  { done = 1;
	    for(i=0;i<nnodes[depth]-1;i++) 
	      if(hord[row[i+1]->i*n+row[i]->i] == 1)
		{ struct node *tmp = row[i];
		  row[i] = row[i+1];
		  row[i+1] = tmp; 
		  done = 0; }}}

      /* Save sorted row in layout */
      { struct nlist **p = &layout[depth];
	layout[depth] = NULL;
	for(i=0;i<nnodes[depth];i++) 
	  { *p = (struct nlist*) metabAlloc(sizeof(struct nlist));
	    (*p)->node = row[i];
	    (*p)->next = NULL;
	    p = &(*p)->next; } }
      msg("end for depth in"); 
    } /* end foreach depth */
msg("End for depth after");

  /* Assign coords */	  
  msg("Layout");
  for(i=0;i<=maxd;i++)
    { struct nlist *l = layout[i];
      if(Debug)printf("%d:",i);
      while(l!=NULL)
	{ if(Debug)printf(" %s",ElementLabel(l->node->label));
	  l = l->next; }
      if(Debug)printf("\n"); }
  
  rowheight = 15;
  if(rowheight*maxd > 90) rowheight = 90/maxd;
  if(rowheight<4) rowheight = 4;

  msg("Assign coords");
  { int maxelts = 0, done, depth, i, j;
    /* Assign vertical coords */
    for(i=0;i<n;i++) 
      { node[i]->x = -1;
	node[i]->y = node[i]->depth*rowheight + 3;  }
    if(1) /* 1 means equispace each depth; 0 is experimental */
      { for(depth=0;depth<=maxd;depth++) 
	  { float h = 1.0, xscale = 70.0/(nnodes[depth]+1);  
	    struct nlist *l = layout[depth];
	    while(l!=NULL)
	      { i = l->node->i;
		node[i]->x = h++*xscale;
		l = l->next; }
	  } }
    else 
      {  /* Find widest row */
	msg("Widest row");
	for(depth=0;depth<=maxd;depth++) 
	  if(nnodes[depth]>maxelts)
	    { maxelts = nnodes[depth]; }
	/* Equispace nodes in widest row(s) */
	msg("Equispace widest");
	for(depth=0;depth<=maxd;depth++)
	  if(nnodes[depth]==maxelts)
	    { struct nlist *l = layout[depth];
	      while(l!=NULL)
		{ int h = 1; i = l->node->i;
		  for(j=0;j<n;j++) 
		    if(node[j]->depth==depth && hord[j*n+i]==1) h++;
		  node[i]->x = ((float) h)/(nnodes[depth]+1)*70;  
		  l = l->next; }
	    } 
	/* Propagate horizontal coords */
	msg("Propagate");
	done = 0; 
	while(! done)
	  { done = 1;
	    for(i=0;i<n;i++)
	      if(node[i]->x == -1)
		{ int inknown = 1, outknown = 1;
		  float inx = 0, outx = 0;

		  for(j=0;j<n;j++)
		    { if(g[i*n+j]) 
			{ if(node[j]->x == -1) outknown = 0;
			else outx += node[j]->x; }
		      if(g[j*n+i]) 
			{ if(node[j]->x == -1) inknown = 0; 
			else inx += node[j]->x; } }
		  if(inknown && node[i]->indeg)
		    { done = 0;  node[i]->x = inx / node[i]->indeg; }
		  else if(outknown && node[i]->outdeg)
		    { done = 0;  node[i]->x = outx / node[i]->outdeg; }
		}}
	
	/* Equispace any unassigned nodes */
	msg("Equispace rest");
	for(i=0;i<n;i++) 
	  if(node[i]->x == -1)
	    { int h = 1;
	      for(j=0;j<n;j++)
		if(node[j]->depth==node[i]->depth&&hord[j*n+i]==1)  h++;
	      node[i]->x = ((float) h)/(nnodes[depth]+1)*70; } 
      }
    /* Adjust any overlaps */
    for(depth=0;depth<maxd;depth++)
      { struct nlist *l = layout[depth]; float low = 5;
	while(l->next!=NULL)
	  { struct node *a = l->node, *b = l->next->node;
	    if(a->x == b->x)
	      { float high; int n = 2;
		{ struct nlist *m = l->next;
		  while(m!=NULL && m->node->x == a->x){ n++;  m = m->next; }
		  if(m==NULL) high = 70;
		  else high = m->node->x; }
		{ struct nlist *m = l; float x0 = a->x, x = low, delta = (high-low)/n;
		  while(m!=NULL && m->node->x == x0)
		    { x += delta; m->node->x = x;  m = m->next; }}}
	    low = a->x;
	    l = l->next; }}
  }
  /* Puff cycles */
  { struct cycle *c = cycles, *c1;
    for(i=0;i<nsq;i++) hord[i] = 0;
    while(c!=NULL)
      { if(Debug)printf("Puffing cycle#%d\n",c->id);
	if(c->type == MAIN)
	  { float xc, yc, a, da, r, r0;
	    struct nlist *l = c->nodes;
	    struct node *prev;
	    /* Figure cycle center, radius */
	    xc = yc = 0;
	    while(l->next!=NULL)
	      { struct node *n = l->node;
		xc += n->x; yc += n->y;
		l = l->next; }
	    xc = xc/c->length;
	    yc = yc/c->length;
	    da = 2*Pi/c->length;
	    /* First estimate of cycle radius based on spacing perimeter nodes
	       equally rowheight apart. */
	    r0 = r = rowheight*c->length/(2*Pi);
	    /* Expand radius if needed so that shunts are not too tightly
	       packed inside circle. */
	    c1 = cycles;
	    while(c1!=NULL)
	      { if(c1->type == SHUNT && c1->main == c )
		  { float ca, r1; struct nlist *last;
		    last = c1->nodes;
		    while(last->next!=NULL) last = last->next;
		    ca = da*chord_angle(c1->nodes->node,last->node,c);
		    if(ca)
		      { if(Debug)printf("rh=%6.2f l=%d ca=%6.2f sin=%6.2f\n",
			       rowheight,c1->length,ca,2*sin(0.5*ca));
			r1 = rowheight*c1->length/(2*sin(0.5*ca));
			if(r1>r)
			  { r = r1;
			    if(r>r0*3) r = r0*3; /* don't go crazy...*/
			  }}}
		c1 = c1->next; }
	    /* Assign node positions */
	    a = c->phase;
	    l = c->nodes;
	    while(l->next!=NULL)
	      { struct node *n = l->node;
		n->x = xc + r*cos(a);
		n->y = yc - r*sin(a); 
		a += da;
		l = l->next; }
	    /* Assign arc params -- use hord to keep track of on-cycle arcs */
	    l = c->nodes->next; prev = c->nodes->node;
	    while(l!=NULL)
	      { hord[prev->i*n+l->node->i] = c->id;
		prev = l->node;
		l = l->next; }
	    e = d->elements;
	    while(e!=NULL)
	      { if(e->type == DIAG_ARC)
		  { struct arc *a = e->item;
		    struct node *from = a->from->item, *to = a->to->item;
		    if(a->tangent_to == NULL && hord[from->i*n+to->i])
		      { a->curvature = 1;
			a->radius = r; }}
		e = e->next; }}
	else if(c->type==SHUNT)
	  /* Lay shunts out linearly. Assumes x & y of shunt endpoints
	     are already assigned. */
	  { float x1, y1, x2, y2, dx, dy;
	    struct node *a = c->nodes->node, *b;
	    struct nlist *l = c->nodes;
	    while(l->next!=NULL) l = l->next;
	    b = l->node;
	    x1 = a->x; y1 = a->y;
	    x2 = b->x; y2 = b->y;
	    StanAssert(x1 != 0 && y1 != 0 && x2 != 0 && y2 != 0);
	    if(Debug)printf("Shunt from %s (%6.2f, %6.2f) to %s (%6.2f, %6.2f)\n",
		   ElementLabel(a->label),x1,y1,ElementLabel(b->label),x2,y2);
	    dx = (x2-x1)/(c->length-1);
	    dy = (y2-y1)/(c->length-1);
	    l = c->nodes->next;
	    while(l->next!=NULL) 
	      { x1 += dx; y1 += dy;
		if(Debug)printf("Placing %s at (%6.2f,%6.2f)\n",
		       ElementLabel(l->node->label),x1,y1);
		l->node->x = x1; l->node->y = y1;
		l = l->next; }
	  }
	c = c->next; }}
	
  msg("Assignment complete");

  /* Position node labels */
  e = d->elements;
  while(e!=NULL)
    { if(e->type==DIAG_NODE && ((struct node *)e->item)->label!=NULL)
	{ struct node *n = e->item; 
	  struct text *t = n->label->item;
	  t->x = n->x + NODE_LABEL_DX;
	  t->y = n->y + NODE_LABEL_DY;  }
      e = e->next; }
  
  /* Recompute arcs */
  e = d->elements;
  while(e!=NULL)
    { if(e->type==DIAG_ARC && ((struct arc *)e->item)->tangent_to == NULL) 
	RecomputeArcParams(e->item);
      e = e->next; }
  
  /* Fudge label positions */
  if(FudgeLabels)
    { int done = 0;
      msg("Fudging Labels");
      e = d->elements;
      while(e!=NULL)
	{ if(e->type==DIAG_TEXT)
	    { struct text *t = e->item;
	      e->x = e->x0 = t->x; 
	      e->y = e->y0 = t->y;  }
	  e = e->next; }
      while(! done)
	{ done = 1;
	  e = d->elements;
	  while(e!=NULL)
	    { if(e->type==DIAG_TEXT && Movable(e))
		{ struct element *e1 = d->elements;
		  while(e1!=NULL)
		    { if(e1->type==DIAG_TEXT
			 &&((struct text *)e->item)->id < ((struct text *)e1->item)->id
			 &&(Movable(e) || Movable(e1))
			 &&diagBoxOverlap(e,e1))
			{ if(diagFixOverlap(e,e1))
			    { done = 0; 
			      break; } }
		      e1 = e1->next; }}
	      if(done==0) break; 
	      e = e->next; }}}
  msg("-Layout");
}

static BOOL diagBoxOverlap(struct element *a, struct element *b)
{ if(a->x + a->width < b->x
     || b->x + b->width < a->x
     || a->y + a->height < b->y
     || b->y + b->height < a->y)
    return FALSE;
  else return TRUE; }

static float sq(float x) { return x*x; }

static float badness(struct element *e)
{ float dx = e->x0 - e->x, dy = e->y0 - e->y;
  return dx*dx + dy*dy; 
}

/* Fix label overlap by moving one or both labels to eliminate
overlap, while trying to minimize total badness = distance
of label from desired position. Note that this fixes only
the pairwise overlap; it may introduce new overlaps with
other labels. Infinite loops of adjustment are prevented
by requiring that the badness of individual labels never
decreases, although this may prevent exploitation of
serendipitous opportunities for improved positions resulting
from vacancies left by other fixes. */

static BOOL diagFixOverlap(struct element *a, struct element *b)
{ float xa1, ya1, xb1, yb1, bada1, badb1, 
    bestbadness, bestax, bestay, bestbx, bestby,
    epsilon = 1.0,
    infinity = 10000000000000.0;

#define swap	{ struct element *tmp = a; a = b; b = tmp; }

#define save_best1	{ float bada = badness(a); \
			  if(bada > bada1 && \
			     (bada < bestbadness || \
			      (bada == bestbadness && \
			       randfloat() <= 0.5))) \
			    { bestax = a->x; bestay = a->y; \
				bestbadness = bada; }}

#define save_best2(text) { float bada = badness(a); \
			   float badb = badness(b); \
			   float bad = bada + badb; \
			   if(bada >= bada1 && badb >= badb1 && \
			      (bad < bestbadness || \
			       (bad == bestbadness && \
				randfloat() <= 0.5))) \
			     { msg(text); \
			       bestax = a->x; bestay = a->y; \
			       bestbx = b->x; bestby = b->y; \
			       bestbadness = bad; }}

  StanAssert(diagBoxOverlap(a,b));
  if(Movable(a))
    { if(Movable(b) && badness(b)<badness(a)) swap ;}
  else swap ;
  bestax = xa1 = a->x; bestay = ya1 = a->y; bada1 = badness(a); 
  bestbx = xb1 = b->x; bestby = yb1 = b->y; badb1 = badness(b);
  bestbadness = infinity;
  if(Movable(a)&&Movable(b))
    { StanAssert(bada1<=badb1); }
  /* Above: */
  a->y = b->y - (a->height + epsilon); save_best1;
  /* Below: */
  a->y = b->y + b->height + epsilon; save_best1;
  /* Left: */
  a->y = ya1;  a->x = b->x - (a->width + epsilon); save_best1;
  /* Left: */
  a->x = b->x + b->width + epsilon; save_best1;
  if (Movable(a) && Movable(b) && bestbadness > badb1)
    { /* Try solutions involving moving both:
	 there are 4, as above: A on left, right, top or bottom
	 but now displacement (in one axis only) is shared between
	 a and b in such a way that the badness of both ends up equal,
	 which means that the combined badness is minimal, and so is the
	 maximum badness. */
      bestbadness += badb1;
      /* A Left */
      a->y = ya1; b->y = yb1;
      a->x = 0.5*(sq(a->width + epsilon - b->x0) + sq(b->y - b->y0)
		  - sq(a->x0) - sq(a->y - a->y0))
	/(a->width - b->x0 - a->x0);
      b->x = a->x + a->width + epsilon;
      save_best2("A Left");
      /* B Left */
      a->y = ya1; b->y = yb1;
      b->x = 0.5*(sq(b->width + epsilon - a->x0) + sq(a->y - a->y0)
		  - sq(b->x0) - sq(b->y - b->y0))
	/(b->width - a->x0 - b->x0);
      a->x = b->x + b->width + epsilon;
      save_best2("A Right");
      /* A Below */
      a->x = xa1; b->x = xb1;
      a->y = 0.5*(sq(a->height + epsilon - b->y0) + sq(b->x - b->x0)
		  - sq(a->y0) - sq(a->x - a->x0))
	/(a->height - b->y0 - a->y0);
      b->y = a->y + a->height + epsilon;
      save_best2("A Below");
      /* B Below */
      a->x = xa1; b->x = xb1;
      b->y = 0.5*(sq(b->height + epsilon - a->y0) + sq(a->x - a->x0)
		  - sq(b->y0) - sq(b->x - b->x0))
	/(b->height - a->y0 - b->y0);
      a->y = b->y + b->height + epsilon;
      save_best2("A Above");
    }
  if(bestbadness<infinity)
    { struct text *t = a->item;
      t->x = a->x = bestax; t->y = a->y = bestay;
      t = b->item; t->x = b->x = bestbx; t->y = b->y = bestby;
      if(Debug) printf("Adjusting %s from (%6.2f,%6.2f) to (%6.2f,%6.2f) badness=%10.3f \n",
	     ElementLabel(a),xa1,ya1,a->x,a->y,bestbadness);
/*      StanAssert(! diagBoxOverlap(a,b));      */
      return TRUE; }
  else return FALSE;
}

static BOOL Movable(struct element *e) { return TRUE; }
  
/* not used 
static void diagPrintElement(struct element *e)
{ switch(e->type){
 case DIAG_ARC:
  { struct arc *a = e->item;
    printf("%s Arc %s from %s to %s\n",
	   (a->tangent_to!=NULL)?"Tangent":
	   a->curvature?((a->curvature==1)?"+Curved":"-Curved"):"Straight",
	   (a->tangent_to!=NULL)?ElementLabel(a->tangent_to):ElementLabel(e),
	   ElementLabel(a->from),ElementLabel(a->to)); }
  break;
 case DIAG_NODE:
  { struct node *n = e->item;
    printf("%s node %s\n",NodeTypeName[n->type],ElementLabel(e)); }
  break;
 default: printf("%s %s\n",ElementTypeName[e->type],ElementLabel(e));}}


static void diagPrint(struct diagram *d)
{ struct element *e = d->elements;
  printf("Diagram %s\n",d->name);
  while(e!=NULL)
    { diagPrintElement(e);
      e = e->next; }
  msg("EOD");
}
*/

static int HighlightColor = GREEN;

/*
static struct element *diagLabelOf(struct element *Label)
{ struct element *e = ActiveDiagram()->elements;
  while(e!=NULL)
    { switch(e->type){
    case DIAG_NODE: if(((struct node *) e->item)->label == Label) return e; break;
    case DIAG_ARC: if(((struct arc *) e->item)->label == Label) return e; }
      e = e->next; }
  return NULL; }
*/
static void diagHighlight0(struct element *e, int s)
{ if(e==NULL) return;
  if(s) graphBoxDraw(e->box,HighlightColor,TRANSPARENT);
  else  graphBoxDraw(e->box,e->color,TRANSPARENT); 
}

static void diagHighlight(struct element *e) 
{ int s; struct element *e1 = ActiveDiagram()->elements;
  if(e == NULL) return;
  s = e->highlight = ! e->highlight;
  if(Debug) printf("%sHighlight %s\n",s?"+":"-",ElementLabel(e));
  while(e1!=NULL)
    { if(e1->key == e->key) diagHighlight0(e1,s);
      e1 = e1->next; }}
   
static void diagPick(double x, double y)
{ /* Note: this is not an ACEDB pick handler, because
     can`'t rely on them to identify picked object.
     Hence this is a LEFT_DOWN handler, does its own
     highlighting. */
  struct diagram *d = ActiveDiagram(); 
  struct element *e, *oldsel = d->selection;
msg("+diagPick");
  e = Pos2Element(d, x, y);
  if(e!=NULL && e->key == 0) e = NULL;
  d->selection = e;
  diagHighlight(oldsel); 
  if(e == NULL) return;
  if(oldsel != NULL && e->key == oldsel->key)
    { d->selection = NULL;
      if(Debug) printf("Picking %s\n",ElementLabel(e)); 
      display(e->key,0,0); }
  else { if(Debug) printf("Pre-picking %s\n",ElementLabel(e)); 
	 diagHighlight(e); }
msg("-diagPick");
}
 

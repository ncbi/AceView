/*  dendrogram.c : implementation file
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Defines the class model behaviors for new ?Tree class
 *		which can draw dendrograms of phylogenetic or taxonomic trees.
 *
 * Exported:
 *	These routines read in a New Hampshire formatted tree:
 *	void readTreeFile(int treeType) ; // see dendrogram.h
 *	void readTaxonomicTree(void)
 *	void readDNATree(void)
 *	void readProteinTree(void)
 *
 *	No function currently provided for "Cell lineage" type trees -- put in explicitly
 *
 * Assumed Class Models:
 *
 * ?Tree    Description	?Text
 *          Type UNIQUE Taxonomy                        // One can add to these but 
 *                      DNA                             // only these two are "in the code"
 *                      Protein
 *                      Cell_Lineage
 *          Root UNIQUE	?TreeNode			// Unique tree node of "root" subtree
 *          Tree_Node   ?TreeNode			// hooks to ?TreeNode's of other ?Tree's embedding this ?Tree
 *          Display     No_Header                       // Suppresses display "header"
 *                      Descriptive_Labels              // Show descriptive labels
 *			Colour                          // Taxon and line colouring, if indicated
 *                      Normalization    UNIQUE Float   // Normalization factor for display (defaults to 1.0)
 *                      Bootstrap_Factor UNIQUE Float   // Normalization factor for bootstrap values (defaults to 1.0)
 *                      Hide_Bootstraps                 // Global visibility of tree bootstrap values
 *                      Alignment     UNIQUE Top        // How the dendrogram is to be drawn
 *                                           Middle
 *                                           Bottom
 *                                           Unrooted   Angle	Float // Orientation of unrooted tree in degrees
 *                                                      Spread	Float // -360..+360 degrees tree spread
 *          Reference	?Paper
 * 
 * ?TreeNode	Label           ?Text                    // Tree vertex label, e.g. sequence name
 *              Id              Int                      // Optional node number
 *              Description     ?Text
 *              Type UNIQUE     Root                     // Optional node labels
 *                              Interior
 *                              Leaf
 *              Distance  UNIQUE Float                   // branch length == from parent node
 *              Bootstrap UNIQUE Float                   // bootstrap value for node subtree
 *              Tree   UNIQUE   ?Tree
 *              Parent UNIQUE   ?TreeNode XREF  Child
 *              Child           ?TreeNode XREF  Parent
 *              Display  Hide                            // Hide the subtree (children) of this node
 *                       Colour UNIQUE  #Colour          // Fixes the colour of the subtree; 
 *                       Hide_Bootstraps                 // Subtree visibility of tree bootstrap values
 *                                                       // overridden by child node settings
 *              // "Contains" is now fully tag2 for associated data display
 *		// The tags listed below, however, have special status in the code at present
 *              Contains Embedded_Tree UNIQUE ?Tree XREF Tree_Node
 *			 Taxon	  UNIQUE  ?Taxon    XREF Tree_Node
 *                       DNA      UNIQUE  ?Sequence XREF Tree_Node
 *                       Protein  UNIQUE  ?Protein  XREF Tree_Node
 *                       Cell     UNIQUE  ?Cell     XREF Tree_Node  // Cell lineages...
 *                       URL	  UNIQUE  ?Url      XREF Tree_Node  // Associated data
 *                       Pick_me_to_call UNIQUE Text Text           // Can show a gif?
 *              Function Map ?Map XREF     Tree_Node #Map_Position  // Leaf node ordinate
 *                       Gene_Function     ?Gene_Function XREF PosTreeNode
 *                       Not_Gene_Function ?Gene_Function XREF NegTreeNode
 *
 * ?Taxon   Common_name UNIQUE ?Text
 *          Description  ?Text
 *          Rank  UNIQUE Superkingdom
 *                       Kingdom
 *                       Phylum
 *                       Subphylum
 *                       Superclass
 *                       Class
 *                       Subclass
 *                       Superorder
 *                       Order
 *                       Suborder
 *                       Superfamily
 *                       Family
 *                       Subfamily
 *                       Genus
 *                       Species ?Species XREF Taxon
 *                       No_Rank
 *          Other_names  Text
 *          Taxonomy     UNIQUE ?TreeNode       // Unique Taxonomy tree
 *          Tree_Node    ?TreeNode              // Other trees dereferencing this taxon?
 *          // Optional colour coding (i.e. on Phylogenetic trees)
 *          Display      Foreground_Colour   #Colour  
 *	                 Background_Colour   #Colour 
 *
 * ?Species Taxon UNIQUE ?Taxon XREF Species
 *          Common_name  UNIQUE ?Text
 *          Description  ?Text
 *          Loci         ?Locus
 *          Sequences    ?Sequence
 *          Proteins     ?Protein
 *          Reference    ?Paper
 *
 * ?Gene_Function Description Text
 *                Map ?Map XREF Gene_Function #Map_Position // For a feature map of genes
 *                Positive PosTreeNode ?TreeNode XREF Gene_Function
 *                Negative NegTreeNode ?TreeNode XREF Not_Gene_Function
 *                Reference ?Paper
 *
 * HISTORY:
 * Last edited: Jun 25 11:38 1999 (fw)
 * * Aug 31 09:38 1998 (rbrusk): several enhancements
 * *  Aug 22 01:04 1998 (rd)
 * Created: May 31 10:50 1998 (rbrusk)
 *-------------------------------------------------------------------
 */

#define DEFINE_OBJ
typedef struct sobj *OBJ ;
#include "acedb.h"

#include "bs_.h"
#include "dendrogram.h"
#include "graphAcedbInterface.h"			    /* graph -> ace header. */
#include "biblio.h"
#include "dump.h"
#include "dna.h"
#include "peptide.h"
#include "tree.h"		 /* for treeUpdate() */
#include "help.h"		 /* for helpOn() */
#include "graphAcedbInterface.h" /* for controlGetColour() */

/* Tree display: alignment */
#define TOP		1
#define MIDDLE		2
#define BOTTOM		4
#define UNROOTED	8	/* "star" unrooted phylogenetic tree */
#define pi		3.141592653

#define BUF_WIDTH	100

typedef struct neighbour_struct
{  int previous, next ; } NEIGHBOURS ;

/* dgram_box_struct flags... */
#define LABELBOX 0x1
#define NODEBOX  0x2
#define HIDDEN   0x4

typedef struct dgram_box_struct 
{
	char  flags ;
	OBJ   obj ;
	short leafPos ;
	short nbIdx ;
	char  colour ;
} DGRAM_BOX ;

/* Should find a labelled ?TreeNode object? */
#define GETOBJ(l,b,d) \
	if(!( (b) && \
	     ((d) = arrp((l)->boxes,(b),DGRAM_BOX)) && \
	     ((d)->flags & LABELBOX) \
	     ) ) return ; 

typedef struct LOOKSTUFF
  { int   magic ;        /* == MAGIC */
    KEY   key ;
    char  title[BUF_WIDTH],
          desc[BUF_WIDTH],
          common[BUF_WIDTH],
          taxon[BUF_WIDTH] ;
    OBJ   objet ;
    Graph graph ;
    int   activebox,  /* always a label box -- corresponding nodeBox is activebox - 1 */
	  rootbox,
	  firstbox,
	  lastbox ;
    KEY   activeKey ;
    AC_HANDLE handle ;
    Associator key2obj ;
    Associator key2hide ;
    Array boxes ;
    Associator obj2box ; 
    Array nodeNeighbours ;
    short nbSeq ;
    Array leaves ;
    short leafPos ;
    float xmax, ymax ;
    char  treeType ;	 /* TAXONOMIC_TREE, DNA_TREE, etc.*/
    float normalFactor,  /* Maximum Branch length display normalized size, in user coordinates */
          distFactor,	 /* Computed Branch distance display normalization factor */
          bootstrapFactor,  /* Bootstrap values normalization factor */
	  angle,	 /* orientation of an unrooted tree */
	  spread,	 /* arc "spread" of an unrooted tree */
	  phi ;		 /* angular decrement = f(spread, nLeaves) */
    char  alignment ;	 /* TOP, MIDDLE, BOTTOM */
    BOOL  bHeader,	 /* flag for presence/absence of graph header */
          bDescriptiveLabels ,	 /* flag for displaying object descriptions versus ?TreeNode names only */
          bColour ,	 /* flag to turn taxon colouring "highlight" on/off */
          bHideBootstraps ; /* current status of bootstrap value display */
    char  activeColour ; /* current subtree colouring (dynamic during drawSubTree */
    char  messageText[BUF_WIDTH] ;
    int   messageBox ;
  } *LOOK ;

/* When set to a look, the user can force
   reuse of the same display over and over again? */
static LOOK  fixedLook = 0 ;
static Graph fixedGraph = 0 ;  

static int MAGIC = 666999 ;

static void 
	dgramHelp(void),
	dgramHeader(void),
	dgramLabels(void),
	dgramColour(void),
	dgramHideBootstraps(void),
	dgramRoot(void),
	dgramHideSubtree(void),
	dgramShowText(void),
	dgramUpdateNode(void),
	dgramFirstLeaf(void),
	dgramLastLeaf(void),
	dgramBiblio(void),
	dgramAceDump(void),
	dgramLeafDump(void),
	dgramFastaDump(void),
	dgramAlignTop(void),
	dgramAlignMiddle(void),
	dgramAlignBottom(void),
#if defined(WIN32_TREE)
	dgramUnRooted(void),
#endif
	dgramFixedDisplay(void) ,
	dgramRescale(void) ;

static MENU dgramMenu ; 
static MENUOPT dgMenuItems[] = {
	{graphDestroy,		"Quit"},
	{dgramHelp,		"Help"},
	{graphPrint,		"Print"},
#if defined(WIN32)
	{graphCopy,		"Copy"}, 
#endif
	{menuSpacer, ""},
	{dgramRoot,		"Select Root Node"},
	{dgramFirstLeaf,	"Goto First Leaf   (Home key)"},
	{dgramLastLeaf,		"Goto Last Leaf    (End key)"},
	{dgramHideSubtree,	"Hide/Show Subtree (Delete key)"},
	{dgramShowText,		"Show Tree Node    (Enter)"},
	{dgramUpdateNode,	"Update Tree Node  (Insert key)"},
	{dgramBiblio,		"Biblio"},
#if defined(WIN32)
	{ 0, "Export Tree Data..."},
#endif
	{dgramAceDump,		"Dump tree data as .ace"},
	{dgramLeafDump,		"Dump tree leaf data"},
	{dgramFastaDump,	"Dump tree sequences as Fasta"},
#if defined(WIN32)
	{ 0, 0 }, 
#endif
	{menuSpacer, ""},
#if defined(WIN32)
	{ 0, "Configure Display..."},
#endif
	{dgramHeader,		"Toggle Header"},
	{dgramLabels,		"Toggle Labels (Tab key)"},
	{dgramColour,		"Toggle Colour (Esc key)"},
	{dgramHideBootstraps,	"Toggle Bootstrap Values"},
	{dgramAlignTop,		"Align Tree to Top"},
	{dgramAlignMiddle,	"Align Tree to Middle"},
	{dgramAlignBottom,	"Align Tree to Bottom"},
#if defined(WIN32_TREE)
	{dgramUnRooted,		"Draw as Unrooted Tree"},
#endif
	{dgramRescale,		"Rescale Tree"},
	{displayPreserve,	"Preserve"},
	{dgramFixedDisplay,	"Anchor Tree Display"},
#if defined(WIN32)
	{ 0, 0 }, 
#endif
	{ 0, 0 }
} ;
static void  dgramClear (LOOK look) ;
static void  dgramDestroy (void) ;
static void  dgramRedraw  (void) ;
static void  dgramPick    (int k) ;
static void  dgramKbd     (int k) ;
static void  dgramScroll  (double x, double y) ;
static OBJ   dgramGet     (LOOK look, KEY key) ;
static void  dgramDraw    (LOOK look);

#define LOOKGET(name)     LOOK look ; \
     			  if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("%s can't find graph",name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer",name) ;\
			  graphActivate(look->graph) 

#define NORMAL_MODE TRUE
#define SELECT_MODE FALSE
#define STD_FCOL BLUE
#define STD_BCOL PALEYELLOW
#define HL_BCOL GREEN

/* Utility routines for converting one class of 
   object reference into a ?Tree object reference */
static OBJ treeNode2Tree(KEY *activeKey, OBJ treeNode)
{
	KEY treeKey ;
	OBJ treeObj = 0 ;

	/* Retrieve the associated ?Tree, if available...*/
	if( strcasecmp(className(treeNode->root->key),"TreeNode") ||
	    !bsGetKey(treeNode,str2tag("Tree"),&treeKey) ||
	    strcasecmp(className(treeKey),"Tree") || 
	    pickType (treeKey) != 'B' || 
	    !(treeObj = bsCreate(treeKey)) )
	{
		*activeKey = 0 ;
		return 0 ;
	}
	else
	{
		/* Set the TreeNode KEY as the active key */
		*activeKey = treeNode->root->key ;

		/* and return the ?Tree object */
		return treeObj ;
	}
}

static OBJ containsTreeNode2Tree(KEY *activeKey, OBJ containsObj)
{
	KEY treeNodeKey ;
	OBJ treeNodeObj = 0, treeObj = 0 ;

	/* Retrieve the associated ?TreeNode, if available...*/
	if( strcasecmp(className(containsObj->root->key),"Cell") ||
	    (!bsGetKey(containsObj,str2tag("Tree_Node"),&treeNodeKey) ||
	    strcasecmp(className(treeNodeKey),"TreeNode") || 
	    pickType (treeNodeKey) != 'B' || 
	    !(treeNodeObj = bsCreate(treeNodeKey)) ||
	    !(treeObj = treeNode2Tree(activeKey,treeNodeObj))) )
	{
		if(treeNodeObj) bsDestroy(treeNodeObj) ; 
		return 0 ;
	}
	else
	{
		bsDestroy(treeNodeObj) ; 
		return treeObj ;
	}
}

static OBJ taxon2Tree(KEY *activeKey, OBJ taxonObj)
{
	KEY treeNodeKey ;
	OBJ treeNodeObj = 0, treeObj = 0 ;

	/* Retrieve the associated ?TreeNode, if available...*/
	if( (strcasecmp(className(taxonObj->root->key),"Taxon")) ||
	    (!bsGetKey(taxonObj,str2tag("Taxonomy"),&treeNodeKey) ||
	    strcasecmp(className(treeNodeKey),"TreeNode") || 
	    pickType (treeNodeKey) != 'B' || 
	    !(treeNodeObj = bsCreate(treeNodeKey)) ||
	    !(treeObj = treeNode2Tree(activeKey,treeNodeObj))) )
	{
		if(treeNodeObj) bsDestroy(treeNodeObj) ; 
		return 0 ;
	}
	else
	{
		bsDestroy(treeNodeObj) ; 
		return treeObj ;
	}
}

static OBJ species2Tree(KEY *activeKey, OBJ species)
{
	KEY taxonKey ;
	OBJ taxonObj = 0, treeObj = 0 ;

	/* Retrieve the associated ?Taxon, if available...*/
	if( strcasecmp(className(species->root->key),"Species") ||
	    !bsGetKey(species,str2tag("Taxon"),&taxonKey) ||
	    strcasecmp(className(taxonKey),"Taxon") || 
	    pickType (taxonKey) != 'B' || 
	    !(taxonObj = bsCreate(taxonKey)) ||
	    !(treeObj = taxon2Tree(activeKey,taxonObj)) )
	{
		if(taxonObj) bsDestroy(taxonObj) ; 
		return 0 ;
	}
	else
	{
		bsDestroy(taxonObj) ; 
		return treeObj ;
	}
}

/************************************************************/
/* The main display function...*/
BOOL dendrogramDisplay (KEY key, KEY from, BOOL isOldGraph)
{
	OBJ obj = 0, dgObj = 0 ;
	KEY activeKey = 0, dgKey = key ;
	float nfac ;
	KEY desckey ;
	LOOK look ;
	
	if ( pickType (key) != 'B' || !(obj = bsCreateCopy(key)) )
		goto abort1 ;

	if(!strcasecmp(className(key),"Tree"))
	{
		activeKey = dgKey = key ;
		dgObj = obj ;
	}
	else if((dgObj = treeNode2Tree(&activeKey,obj)))
	{
		bsDestroy(obj) ;
		dgKey = dgObj->root->key ;
	}
	else if((dgObj = containsTreeNode2Tree(&activeKey,obj)))
	{
		bsDestroy(obj) ;
		dgKey = dgObj->root->key ;
	}
	else if((dgObj = species2Tree(&activeKey,obj)))
	{
		bsDestroy(obj) ;
		dgKey = dgObj->root->key ;
	}
	else if((dgObj = taxon2Tree(&activeKey,obj)))
	{
		bsDestroy(obj) ;
		dgKey = dgObj->root->key ;
	}
	else goto abort2 ; /*  object type unknown for DtDendrogram? */

	/* In with the old...maybe */
	if (isOldGraph)
	{
		if (!graphAssFind (&MAGIC, &look)) \
			messcrash ("dendrogramDisplay can't find graph") ; \
		if (look->magic != MAGIC) \
			messcrash ("dendrogramDisplay received a wrong pointer") ;\
		dgramClear (look) ; 
	} 
	else if(fixedLook)
	{
		if( graphActivate(fixedGraph))
		{
			look = fixedLook ;
			dgramClear(look) ;
			isOldGraph = TRUE ;
		}
	}
	if(!isOldGraph)  /* otherwise, allocate a new one...*/
		look =(LOOK)messalloc(sizeof(struct LOOKSTUFF));

	look->magic = MAGIC;
	look->activeKey = activeKey ;
	look->key   = dgKey ;
	look->objet = dgObj ;
	*look->desc = 0 ;
        *look->common = 0 ;
        *look->taxon = 0 ;

	if (isOldGraph)
	{ 
		graphClear () ;
		graphGoto (0,0) ;
	}
	else
	{ 
		if (!displayCreate (DtDendrogram))
		{
			messout("Can't create dendrogram display?") ;
			goto abort3 ;
		}
		graphRegister (DESTROY, dgramDestroy) ;
		graphRegister (PICK, dgramPick) ;
		graphRegister (KEYBOARD, dgramKbd) ;
		graphRegister (MIDDLE_DOWN, dgramScroll) ;
		graphRegister (RESIZE, dgramRedraw) ;

		dgramMenu = menuInitialise ("Dendrogram Menu", (MENUSPEC*)dgMenuItems) ;
		graphNewMenu (dgramMenu) ;

		graphHelp("Dendrogram");
	}

	/* Determine type of tree */
	if( bsFindTag(look->objet, str2tag ("Taxonomy")) )
	{
		strcpy(look->title,"Taxonomy:") ;
		look->treeType = TAXONOMIC_TREE ;
	}
	else if( bsFindTag(look->objet, str2tag ("DNA")) )
	{
		strcpy(look->title,"DNA Sequence Phylogeny:") ;
		look->treeType = DNA_TREE ;
	}
	else if( bsFindTag(look->objet, str2tag ("Protein")) )
	{
		strcpy(look->title,"Protein Sequence Phylogeny:") ;
		look->treeType = PROTEIN_TREE ;
	}
	else if( bsFindTag(look->objet, str2tag ("Cell_Lineage")) )
	{
		strcpy(look->title,"Cell Lineage:") ;
		look->treeType = CELL_LINEAGE ;
	}
	else 
	{
		strcpy(look->title,"Dendrogram:") ;
		look->treeType = UNKNOWN_TREE ;
	}

	if(bsGetKey(look->objet,str2tag("Description"),&desckey))
		strncat(look->title,messprintf("\"%s\"",name(desckey)),BUF_WIDTH-strlen(look->title)) ;
	else
		strncat(look->title,messprintf("\"%s:%s\"", className(key),name (key)),BUF_WIDTH-strlen(look->title)) ;

	graphRetitle (look->title) ;
	look->graph	= graphActive() ;
	graphAssociate (&MAGIC, look) ;
	look->handle	= handleCreate() ;
	look->key2obj	= assCreate () ;
	look->key2hide	= assCreate () ;
	look->obj2box	= assCreate () ;
	look->boxes	= arrayReCreate(look->boxes,128,DGRAM_BOX) ;
	look->nodeNeighbours = arrayReCreate(look->nodeNeighbours, 128,NEIGHBOURS) ;

	look->leaves    = arrayCreate(64,int) ;
	look->leafPos   = -1 ;
	look->activebox = look->firstbox = look->lastbox = 0 ;

	/* Branch length display normalization factor */
	if ( bsGetData (look->objet, str2tag ("Normalization"), _Float, &nfac) )
	{
		if(nfac < 0.0)
		{
		      messout("Tree Normalization <= 0.0? Setting it to 1.0!") ;
		      nfac = 1.0 ;
		}
		look->normalFactor = nfac ; 
	}
	else 
		look->normalFactor = 1.0 ; /* if no distance found, assume 1.0 */

	/* Bootstrap normalization factor */
	if ( bsGetData (look->objet, str2tag ("Bootstrap_Factor"), _Float, &nfac) )
	{
		if(nfac < 0.0)
		{
		      messout("Bootstrap Normalization Factor <= 0.0? Setting it to 1.0!") ;
		      nfac = 1.0 ;
		}
		look->bootstrapFactor = nfac ; 
	}
	else 
		look->bootstrapFactor = 1.0 ; /* if no bootstrap normalization factor found, assume 1.0 */

	look->distFactor = 0.0  ;
	look->angle	 = 0.0  ;
	look->spread	 = 360.0 ;
	look->phi	 = 1.0  ;

	if( bsFindTag(look->objet, str2tag ("Top")) )
		look->alignment = TOP ;
	else if( bsFindTag(look->objet, str2tag ("Middle")) )
		look->alignment = MIDDLE ;
	else if( bsFindTag(look->objet, str2tag ("Bottom")) )
		look->alignment = BOTTOM ;
	else if( bsFindTag(look->objet, str2tag ("Unrooted")) )
	{
		float ang ;
		look->alignment = UNROOTED ;

		if(bsGetData (look->objet, str2tag ("Angle"), _Float, &ang))
		{
			/* coerce to within 0.0..360.0 range */
			if(ang < 0.0 || ang > 360.0)
			{
			      messout("Unrooted tree orientation angle %5.1f is outside range [0.0 .. 360.0]? Setting it to 0.0!",ang) ;
			      look->angle = 0.0 ;
			}
			else look->angle = ang ;
		}
		else 
			look->angle = 0.0  ; /* parallel w/ X axis, pointing right */

		if(bsGetData (look->objet, str2tag ("Spread"), _Float, &ang))
		{
			/* coerce to within -360.0..+360.0 range; sign determines curl of plot? */
			if(ang < -360.0 || ang > 360.0)
			{
			      messout("Unrooted tree spread angle %5.1f is outside range [-360.0 .. +360.0]? Setting it to 360.0!",ang) ;
			      look->angle = 360.0 ;
			}
			else look->angle = ang ;
		}
		else
			look->spread = 360.0  ; /* full circle, counterclockwise */
	}
	else if( bsFindTag(look->objet, str2tag ("Protein")) ||
		 bsFindTag(look->objet, str2tag ("DNA"))     ||
		 bsFindTag(look->objet, str2tag ("Cell_Lineage")))
		look->alignment = MIDDLE ;
	else /* "Taxonomy" and all other trees? ) */
		look->alignment = TOP ;

	if( bsFindTag(look->objet, str2tag ("No_Header")) )
		look->bHeader = FALSE ;
	else
		look->bHeader = TRUE ;

	if( bsFindTag(look->objet, str2tag ("Descriptive_Labels")) )
		look->bDescriptiveLabels = TRUE ;
	else
		look->bDescriptiveLabels = FALSE ;

	if( bsFindTag(look->objet, str2tag ("Colour")) )
		look->bColour = TRUE ;
	else
		look->bColour = FALSE ;

	look->activeColour = BLACK ;

	if( bsFindTag(look->objet, str2tag ("Hide_Bootstraps")) )
		look->bHideBootstraps = TRUE ;
	else
		look->bHideBootstraps = FALSE ;

	*look->messageText = 0;
	look->messageBox = -1;

	dgramDraw (look) ;
	return TRUE ;

abort3:	/* (mieg): should be localdestroy to kill look->obj and anything else */
	/* (rbrusk): only look->obj, assuming the look->objet==bsCreateCopy() needs to be destroyed? */
	obj = look->objet ;
	if( fixedLook == look )
	{
		fixedLook = 0 ;
		fixedGraph = 0 ;
	}
	messfree (look) ;
	
abort2:
	bsDestroy(obj) ;
abort1:
	/* fail gracefully: attempt to display object as an ACEDB tree...*/
	display (key,from, TREE) ;

	/* but return false? */
	return FALSE ;
}

static OBJ dgramGet(LOOK look, KEY key)
{
	OBJ obj = 0 ;
	if( !assFind(look->key2obj,assVoid(key),&obj))
	{
		obj = bsCreateCopy(key) ;
		assInsert(look->key2obj,assVoid(key),obj) ;
		if(!strcasecmp(className(key),"TreeNode"))
		  {
		    if( bsFindTag(obj,str2tag("Hide")))
		      assInsert(look->key2hide,assVoid(key),assVoid(TRUE)) ;
		    else
		      assInsert(look->key2hide,assVoid(key),assVoid(FALSE)) ;
		  }
	}
	return obj ;
}

static void dgramClear (LOOK look)
{
  OBJ obj ;
  void *key = 0; /* really a KEY, but... assNext prefers assVoid's */

  /* Need to bsDestroy all my cached dendrogram objects */
  while( assNext(look->key2obj,&key,&obj))
	  bsDestroy(obj) ;
  assDestroy (look->key2obj) ;
  assDestroy (look->key2hide) ;
  arrayDestroy (look->boxes) ;
  assDestroy (look->obj2box) ;
  arrayDestroy(look->nodeNeighbours) ;
  arrayDestroy(look->leaves) ;
  bsDestroy(look->objet) ;
  handleDestroy(look->handle) ;
}

static void dgramDestroy (void)
{
  LOOKGET("dgramDestroy") ;

  dgramClear(look) ;

  look->magic = 0 ;
  if( fixedLook == look )
  {
	fixedLook = 0 ;
	fixedGraph = 0 ;
  }
  messfree (look) ;

  graphAssRemove (&MAGIC) ;
}

/* This function sets ?TreeNode label text box highlighting 
   either to standard values or ?Taxon specified values */
static void labelColour(LOOK look, int box, BOOL bNormal )
{
	DGRAM_BOX  *lpDBox ;
	OBJ obj, keyObj ;
	int fcol, bcol, nodefcol, nodebcol ;
	KEY key = 0 ;

	GETOBJ(look,box,lpDBox)
	obj = lpDBox->obj ;

	/* Default "B & W" colouring */
	bcol = STD_BCOL ;
	nodefcol = fcol = STD_FCOL ;
	nodebcol = WHITE ;

	/* Check for overriding dendrogram colour mode */
	if(look->bColour)
	{
		/* if so, default to subtree colouring...*/
		if(lpDBox->colour != WHITE &&
		   lpDBox->colour != BLACK )
		{
			nodefcol = bcol = lpDBox->colour ;
			nodebcol = fcol = WHITE ;
		}
		
		/* ...but override with Taxon colour coding, if specified */
		if(!( ( 
			bsGetKey(obj, str2tag ("DNA"),&key) ||
	                bsGetKey(obj, str2tag ("Protein"), &key)
		       ) &&
		      (keyObj = dgramGet(look, key)) /* sequence or protein object? */
		         && 
		       bsGetKey(keyObj, str2tag ("Species"),&key)  /* associated species? */
		         &&
		      (keyObj = dgramGet(look, key)) /* ... should have a Taxon tag? */
		     )
		   )
		   /* no associated species found? try the ?TreeNode... */
			keyObj = obj ;

		if( bsGetKey(keyObj, str2tag ("Taxon"),&key) &&
		   (keyObj = dgramGet(look, key))) /* ?Taxon object w/colours? */
		{
			if(bsFindTag(keyObj, str2tag ("Foreground_Colour")))
				fcol = controlGetColour(keyObj);

			if(bsFindTag(keyObj, str2tag ("Background_Colour")))
				bcol = controlGetColour(keyObj);
		}
	}

	/* Select mode? Swap colours */
	if(!bNormal)
	{
		int col = fcol ;
		fcol = bcol ;
		bcol = col ;
		nodebcol = HL_BCOL ;
	}

	graphBoxDraw(box, fcol, bcol) ;

	/* Assume that a nodeBox exists at box - 1 and needs highlighting? */
	graphBoxDraw(box-1, nodefcol, nodebcol) ;
}

/************************************************************/

/* This function returns the description, if any, 
   associated with the specified object argument, of "treeType";
   The tag value is assumed to be a ?Text key in the model.
   Note: subsequent invocations overwrite the old value, so... */
static char *description(LOOK look, OBJ obj, int treeType)
{	
	KEY tkey ;

	*look->desc = 0 ;
	switch(treeType)
	{
		case DNA_TREE:
			if(bsGetKey(obj,str2tag("Brief_identification"),&tkey))
				strncpy((char*)look->desc,name(tkey),
					BUF_WIDTH-1) ;
			break ;
		case PROTEIN_TREE:
			if(bsGetKey(obj,str2tag("Title"),&tkey))
				strncpy((char*)look->desc,name(tkey),
					BUF_WIDTH-1) ;
			break ;
		case CELL_LINEAGE:
			if(bsGetKey(obj,str2tag("Remark"),&tkey))
				strncpy((char*)look->desc,name(tkey),
					BUF_WIDTH-1) ;
			break ;
		case TAXONOMIC_TREE:
		default: ; 
	}
	if(!*look->desc)
	{
		if(bsGetKey(obj,str2tag("Description"),&tkey))
			strncpy((char*)look->desc,name(tkey),
				BUF_WIDTH-1) ;
	}
	return look->desc ;
}

/* This function returns the common or descriptive name, if any, 
   associated with the specified object argument, of "treeType";
   The tag values is assumed to be a "Common_Name" ?Text key in the Taxon model,
   and "Title" ?Text keys in other models (Sequence, Protein).
   Note: subsequent invocations overwrite the old value, so... */
static char *commonName(LOOK look, OBJ obj, int treeType)
{	
	KEY tkey ;

	*look->common = 0 ;
	switch(treeType)
	{
		case TAXONOMIC_TREE:
			if( bsGetKey(obj,str2tag("Common_name"),&tkey) &&
		            class(tkey)==_VText)
				strncpy((char*)look->common,name(tkey),BUF_WIDTH-1) ;
			break ;
		case DNA_TREE:
		case PROTEIN_TREE:
			if(bsGetKey(obj,str2tag("Title"),&tkey) &&
		           class(tkey)==_VText)
				strncpy((char*)look->common,name(tkey),BUF_WIDTH-1) ;
			break ;
		case CELL_LINEAGE:
			if(bsGetKey(obj,str2tag("Fate"),&tkey) &&
		           class(tkey)==_VText)
				strncpy((char*)look->common,name(tkey),BUF_WIDTH-1) ;
			break ;
		default:  ;
	}
	/* If no common name yet found, look for a "Description ?Text" tag-value */
	if(!*look->common)
	{
		/* get first "Contains" tag2 object
		   Description ?Text tag value as name */ ;
		if( bsGetKey(obj,str2tag("Description"),&tkey) &&
		    class(tkey)==_VText )
			strncpy((char*)look->common,name(tkey),BUF_WIDTH-1) ;
	}
	return (char*)look->common ;
}

/* This function expects a ?TreeNode (or ?Species, or other (?)) object
   "obj" containing a Taxon tag with ?Taxon value; 
   returns a full description of the taxon;
    */
static char *taxon(LOOK look, OBJ obj)
{
	KEY key, tkey, rank ;
	char *cp ;

	*look->taxon = 0 ;
	if(bsGetKey(obj,str2tag("Taxon"),&key))
	{
		OBJ taxonObj ;
		if ((taxonObj = dgramGet(look, key)))
		{
			if(bsGetKeyTags(taxonObj,str2tag("Rank"),&rank))
			{
				strncpy(look->taxon, messprintf("%s: ",name(rank)),
						  BUF_WIDTH-1) ;
				strncat(look->taxon, messprintf("%s",name(key)),
					BUF_WIDTH-strlen(look->taxon)) ;
			}
			if(bsGetKey(taxonObj,str2tag("Common_name"),&tkey))
				strncat(look->taxon, messprintf(" (%s)",name(tkey)),
					BUF_WIDTH-strlen(look->taxon)) ;
			if(*(cp = description(look,taxonObj,TAXONOMIC_TREE)))
				strncat(look->taxon, messprintf(" - %s",cp),
					BUF_WIDTH-strlen(look->taxon)) ;
		}
		else strncpy(look->taxon,messprintf("Taxon %s",name(key)), BUF_WIDTH-1) ;
	}
	return (char*)look->taxon ;
}

/* obj is assumed to have a Species ?Species 
   or a Taxon ?Taxon value pair (which should be of rank "species"
   although this is not checked(?)) */
static char *species(LOOK look, OBJ obj)
{
	KEY key ;

	*look->taxon = 0 ;
	if(bsGetKey(obj,str2tag("Species"),&key))
		strncpy(look->taxon,messprintf("%s",name(key)), BUF_WIDTH-1) ;
	else if(bsGetKey(obj,str2tag("Taxon"),&key))
		strncpy(look->taxon,messprintf("%s",name(key)), BUF_WIDTH-1) ;
	return (char*)look->taxon ;
}

/* obj is assumed to be the current ?TreeNode object */
static void messageBox (LOOK look, OBJ obj)
{
	KEY key ;
	char* cp = 0 ;

	if(look->messageBox<0) return ;

	/* Clear the box */
	*look->messageText = 0 ;

	if(!obj) return ;

	switch(look->treeType)
	{
		case TAXONOMIC_TREE:
			if(*(cp = taxon(look,obj)))
				strncpy((char*)look->messageText,cp,
					BUF_WIDTH-1) ;
			break ;
		case DNA_TREE:
			if(bsGetKey(obj,str2tag("DNA"),&key))
			{
				OBJ SequenceObj = 0 ;
				strncpy((char*)look->messageText,
					messprintf("DNA:%s",name(key)),
					BUF_WIDTH-1) ;
				if ((SequenceObj = dgramGet(look, key)))
				{
					if(*(cp = description(look,SequenceObj,DNA_TREE)))
						strncat((char*)look->messageText,
							messprintf(" - %s",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
					else
						strncat((char*)look->messageText,
							messprintf(" - %s",name(SequenceObj->root->key)),
							BUF_WIDTH-strlen(look->messageText)) ;
					/* specify species if known */
					if(*(cp = species(look,SequenceObj)))
						strncat((char*)look->messageText,
							messprintf(" {%s}",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
					else if(*(cp = species(look,obj)))
						strncat((char*)look->messageText,
							messprintf(" {%s}",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
				}
				else if(*(cp = species(look,obj)))
					strncat((char*)look->messageText,
						messprintf(" {%s}",cp),
						BUF_WIDTH-strlen(look->messageText)) ;
			}
			break ;
		case PROTEIN_TREE:
			if(bsGetKey(obj,str2tag("Protein"),&key))
			{
				OBJ proteinObj = 0 ;
				strncpy((char*)look->messageText,
					messprintf("Protein:%s",name(key)),
					BUF_WIDTH-1) ;
				if ((proteinObj = dgramGet(look, key)))
				{
					if(*(cp = description(look,proteinObj,PROTEIN_TREE)))
						strncat((char*)look->messageText,
							messprintf(" - %s",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
					else
						strncat((char*)look->messageText,
							messprintf(" - %s",name(proteinObj->root->key)),
							BUF_WIDTH-strlen(look->messageText)) ;
					/* specify species if known */
					if(*(cp = species(look,proteinObj)))
						strncat((char*)look->messageText,
							messprintf(" {%s}",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
					else if(*(cp = species(look,obj)))
						strncat((char*)look->messageText,
							messprintf(" {%s}",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
				}
				else if(*(cp = species(look,obj)))
					strncat((char*)look->messageText,
						messprintf(" {%s}",cp),
						BUF_WIDTH-strlen(look->messageText)) ;
			}
			break ;
		case CELL_LINEAGE:
			if(bsGetKey(obj,str2tag("Cell"),&key))
			{
				OBJ cellObj = 0 ;
				strncpy((char*)look->messageText,
					messprintf("Cell:%s",name(key)),
					BUF_WIDTH-1) ;
				if ((cellObj = dgramGet(look, key)))
				{
					if(*(cp = description(look,cellObj,CELL_LINEAGE)))
						strncat((char*)look->messageText,
							messprintf(" - %s",cp),
							BUF_WIDTH-strlen(look->messageText)) ;
				}
			}
			break ;
		default: ;
	}
	if( !*look->messageText      && /* is the message empty? */
	    bsFindTag(obj,_Contains) && /* but the ?TreeNode "Contains" something? */
	    bsFindTag(obj,_bsRight)  && /* get to the tag2 tag */
	    bsFindTag(obj,_bsRight)  && /* get to the tag2 object */
	    bsGetKey(obj,_bsHere,&key))
	{
		OBJ containsObj ;
		if ((containsObj = dgramGet(look, key)))
		{
			strncpy((char*)look->messageText,
				messprintf("%s:%s",className(key),name(key)),
				BUF_WIDTH-1) ;
			if(*(cp = description(look,containsObj,UNKNOWN_TREE)))
				strncat((char*)look->messageText,
					messprintf(" - %s",cp),
					BUF_WIDTH-strlen(look->messageText)) ;
		}
	}

	look->messageText[BUF_WIDTH-1] = 0 ;
	graphBoxDraw (look->messageBox, BLACK, PALECYAN) ;
}

/************************************************************/

/* Sets the active ?TreeNode box, with the expected side effects;
   "box" is assumed to be the labelBox of the node */
static void selectActive(LOOK look, int box)
{
	float x,y,x1,y1,x2,y2 ;
	DGRAM_BOX  *lpDBox ;

	if(box)
	{
		GETOBJ(look,box,lpDBox)
		if(lpDBox->obj)
		{
			look->activeKey = lpDBox->obj->root->key ;
			if(look->activebox != box)
			{
				labelColour(look, look->activebox, NORMAL_MODE) ;
				look->activebox = box ;
			}
			labelColour(look, look->activebox, SELECT_MODE) ;

			/* Details about hit in messageText? */
			messageBox(look, lpDBox->obj) ;

			graphWhere(&x1,&y1,&x2,&y2) ;
			graphBoxDim(look->activebox,&x,&y,0,0) ;
			if(x<x1 || x+10.0>x2 || y<y1 || y+1.0>y2)
				graphGoto(x,y) ;

			return ; /* success! */
		}
	}
	if(look->activebox != 0)
		labelColour(look, look->activebox, NORMAL_MODE) ;
	look->activeKey = 0 ;
	look->activebox = 0 ;
	messageBox(look, 0) ; /* clears the message box */
}

/************************************************************/

/* returns TRUE if the graph was redrawn */
#define TOGGLE_NODE 0
#define SHOW_NODE   1
#define HIDE_NODE   2
static void dgramHide(int bShow)
{
	DGRAM_BOX  *lpDBox ;
	OBJ obj = 0 ;
	BS bs ;
	BOOL bHidden = TRUE, redraw = FALSE ;
	int box ;

	LOOKGET("dgramHide") ;

	if(!look->activebox)
	{
		messout("Pick a tree node first!") ;
		return ;
	}

	if(!((lpDBox = arrp(look->boxes,look->activebox,DGRAM_BOX)) && 
	     (lpDBox->flags & LABELBOX) 
	     ) ) return ;
	
	if(!(obj = lpDBox->obj)) return ;

	bs = obj->root ; 

	if (!assFind (look->key2hide,assVoid(bs->key), &bHidden))
	{
		/* default to Hide if flag not already set? */
		assInsert(look->key2hide,assVoid(bs->key),assVoid(TRUE)) ;
		redraw = TRUE ;
	}
	/* bHidden assumed to have a value by this point...*/
	else if (bShow & SHOW_NODE)
	{
		if(bHidden) /* only do something if the active box is hidden */
		{
			assRemove(look->key2hide,assVoid(bs->key)) ;
			assInsert(look->key2hide,assVoid(bs->key),assVoid(FALSE)) ;
			redraw = TRUE ;
		}
	}
	else if (bShow & HIDE_NODE)
	{
		if(!bHidden) /* only do something if the active box is hidden */
		{
			assRemove(look->key2hide,assVoid(bs->key)) ;
			assInsert(look->key2hide,assVoid(bs->key),assVoid(TRUE)) ;
			redraw = TRUE ;
		}
	}
	else 
	{
		assRemove(look->key2hide,assVoid(bs->key)) ;
		assInsert(look->key2hide,assVoid(bs->key),assVoid(!bHidden)) ;
		redraw = TRUE ;
	}

	if(redraw)
	{
		dgramDraw(look) ;

		/* mapping of active box, but not 
		   the object itself, will have changed? Reset back... */
		if(assFind(look->obj2box,obj,&box))
			selectActive(look, box) ;
	}
}

/************************************************************/

/* The object argument to this function is assumed to be a ?TreeNode...*/
#define TEXT_ONLY    FALSE
#define FULL_DISPLAY TRUE
#define DGRAM_TAG2
static void  dgramDisplayNode(LOOK look, OBJ obj, BOOL bFullDisplay)
{
	KEY key = 0, displayKey = 0, speciesKey = 0 ;
	BS bs = obj->root ;
	OBJ containedObj ;
	BOOL bTaxonSeen = FALSE, bDisplayBTree = FALSE ;
	Array containsKeys = 0 ;
	int i ;

	/* preserve the original dendrogram display */
	displayPreserve() ;

	/* Display Tree Node object if node 
	   'Contains' no other information...*/
	if( !bsFindTag (obj, _Contains) )
	     display (bs->key,look->key, "TREE") ;

        else { /* bsFindTag (obj, _Contains) == TRUE */
	   /*   Otherwise, display associated data. 
		Note: 'Contains' is now "tag2" */
	   containsKeys = arrayCreate (32, BSunit) ;
	   if (bsFlatten (obj, 2, containsKeys))
		for (i = 0 ; i < arrayMax(containsKeys) ; i += 2)
		{ 
			if (class(arr(containsKeys, i+1, BSunit).k)) 
			{
				key = arr(containsKeys, i+1, BSunit).k ;
				if(!strncasecmp(className(key),"Taxon",5)) {
					/* display ?Taxon as ACEDB text-only */ 
					bTaxonSeen = TRUE ;
				        bDisplayBTree = TRUE ;
				} else if(   /* display "Contains ?Cell's ...  */
					    !strncasecmp(className(key),"Cell",4) ||

					    /* ...?Species's ...  */
					    !strncasecmp(className(key),"Species",7) ||

					    /* ...and ?TreeNode's ...  */
					    !strncasecmp(className(key),"TreeNode",8) 
					) {
				        bDisplayBTree = TRUE  ; /* as ACEDB text-only in this context */

				} else if( /* otherwise, full display (maybe?) */
					    bFullDisplay &&
					    /* default display should *not* be TREE? */
					   (displayKey = pickDisplayKey (key)) &&
					    strncasecmp(name(displayKey),"TREE",4) 
					) 
				{
					displayPreserve() ;
					/* Note: if the default graph display fails,
					   is is assumed to already display a B Tree instead
					   so I suppress duplicate display of the B Tree? */
					display (key,look->key, 0) ;
				}
				/* iff default graph display was successful,
				   then also show the associated B-tree data :-) */
				if(bDisplayBTree) { 
				  displayPreserve() ;
				  display (key,look->key, "TREE") ;
			        }
			}
		}  /* end for all containsKeys... */
	   arrayDestroy(containsKeys) ;
	}  /* endif !bsFindTag (obj, _Contains) */

	if( !bTaxonSeen && 
	    (containedObj = dgramGet(look, bs->key)) &&
	     bsGetKey(containedObj,str2tag("Species"),&speciesKey) )
	{
		displayPreserve() ;
		/* Force ?Species to display as TREE in this context */
		display (speciesKey,look->key, "TREE") ;
	}

	if (bFullDisplay && bsFindTag(obj,str2tag("Pick_me_to_call")))
		externalDisplay(bs->key) ;
}

/************************************************************/

static void dgramUpdate (LOOK look, OBJ obj)
{
	/* First, display...*/
	display (obj->root->key,look->key, TREE) ;
	/* then update it? */
	treeUpdate() ;
}

/********************* Menu Routines *******************************/
static void dgramHelp(void)
{
      helpOn("Dendrogram") ;
}

static void dgramAlignTop(void)
{
	LOOKGET("dgramAlignTop") ;
	look->alignment = TOP ;
	dgramDraw (look) ;
}

static void dgramAlignMiddle(void)
{
	LOOKGET("dgramAlignMiddle") ;
	look->alignment = MIDDLE ;
	dgramDraw (look) ;
}

static void dgramAlignBottom(void)
{
	LOOKGET("dgramAlignBottom") ;
	look->alignment = BOTTOM ;
	dgramDraw (look) ;
}

#if defined(WIN32_TREE)
/* Unrooted "star" phylogenetic tree */
static void dgramUnrooted(void)
{
      LOOKGET("dgramUnrooted") ;
      look->alignment = UNROOTED ;
      {
		char ang[16], arc[16] ;
		strncpy(ang,messprintf("%.1f",look->angle),15) ;
		ang[15] = '\0' ;
		strncpy(arc,messprintf("%.1f",look->spread),15) ;
		arc[15] = '\0' ;
	getAng:
		if( messPrompt ("Unrooted tree orientation angle [0.0 <= n <= 360.0]:",ang, "f"))
		{
		      float f ;
		      freefloat(&f) ;
		      if(f < 0.0 || f > 360.0)
		      {
			      messout("Orientation angle is out of range 0.0 <= n <= 360.0!") ;
			      goto getAng ;
		      }
		      look->angle = f ;
	      }
	getArc:
		if( messPrompt ("Unrooted tree spread [-360.0 <= n <= 360.0]:",arc, "f"))
		{
		      float f ;
		      freefloat(&f) ;
		      if(f < -360.0 || f > 360.0)
		      {
			      messout("Tree spread is out of range -360.0 <= n <= 360.0!") ;
			      goto getArc ;
		      }
		      look->spread = f ;
	      }
      }
      dgramDraw (look) ;
}
#endif

/* anchor on current look */
static void dgramFixedDisplay(void)
{
	LOOKGET("dgramFixedDisplay") ;
	if(fixedLook == look)
	{
		fixedLook = 0 ; /* turns it off */
		fixedGraph = 0 ;
	}
	else
	{
		fixedLook = look ;
		fixedGraph = graphActive() ;
	}
}

static void dgramRescale(void)
{
	char nfacbuf[16] ;
	LOOKGET("dgramRescale") ;

	strncpy(nfacbuf,messprintf("%.1f",look->normalFactor),15) ;
	nfacbuf[15] = '\0' ;
getAgain:
	if( messPrompt ("Normalization factor [0.0 < n <= 100.0]:",nfacbuf, "f"))
	{
	      float f ;
	      freefloat(&f) ;
	      if(f <= 0.0 || f > 100.0)
	      {
		      messout("Normalization is out of range 0.0 < n <= 100.0!") ;
		      goto getAgain ;
	      }
	      look->normalFactor = f ;
	      dgramDraw (look) ;
      }
}

/* Toggle Labels */
static void dgramLabels(void)
{
	LOOKGET("dgramLabels") ;
	look->bDescriptiveLabels = !look->bDescriptiveLabels ;
	dgramDraw (look) ;
}

/* Toggle Header */
static void dgramHeader(void)
{
	LOOKGET("dgramHeader") ;
	look->bHeader = !look->bHeader ;
	dgramDraw (look) ;
}

/* Toggle Colour */
static void dgramColour(void)
{
	LOOKGET("dgramColour") ;
	look->bColour = !look->bColour ;
	dgramDraw (look) ;
}

/* Toggle Bootstrap display */
static void dgramHideBootstraps(void)
{
	LOOKGET("dgramHideBootstraps") ;
	look->bHideBootstraps = !look->bHideBootstraps ;
	dgramDraw (look) ;
}


/* Select Root Node */
static void dgramRoot(void)
{
	LOOKGET("dgramRoot") ;
	selectActive(look, look->rootbox) ;
}

/* Hide/Show Subtree */
static void dgramHideSubtree(void)
{
	dgramHide(TOGGLE_NODE) ;
}

/* Goto First Leaf */
static void dgramFirstLeaf(void)
{
	LOOKGET("dgramFirstLeaf") ;
	selectActive(look, look->firstbox) ;
}

/* Goto Last Leaf */
static void dgramLastLeaf(void)
{
	LOOKGET("dgramLastLeaf") ;
	selectActive(look, look->lastbox) ;
}

/* Biblio */
static void dgramBiblio(void)
{
	LOOKGET("dgramBiblio") ;
	biblioKey (look->key) ;
}

/* Show as text... */
static void dgramShowText(void)
{
	DGRAM_BOX  *lpDBox ;
	BS bs ;
	OBJ obj ;
	LOOKGET("dgramShowText") ;

	if (!look->activebox)
	{ 
		messout ("Please first select a node") ;
		return ;
	}
        GETOBJ(look,look->activebox,lpDBox)
	obj = lpDBox->obj ;

	if(obj)
		bs = obj->root ;
	else
	{
		messout(messprintf("Sorry, dgramShowText() could not retrieve tree node for box #%d?",
			look->activebox)) ;
		return ;
	}
        if (iskey(bs->key) == 2 )
	{ 
	  if(pickType(bs->key) == 'B')
		dgramDisplayNode(look,obj,TEXT_ONLY) ;
	  else
		messout("Sorry, I cannot display %s as an ACEDB objet", name(bs->key)) ;
	}
      else
	messout("No data associated to %s", name(bs->key)) ;
}

/* Update TreeNode Data */
static void dgramUpdateNode(void)
{
	DGRAM_BOX  *lpDBox ;
	LOOKGET("dgramUpdateNode") ;
	if (!look->activebox)
	{ 
		messout ("Please first select a node") ;
		return ;
	}
	GETOBJ(look,look->activebox,lpDBox)
	dgramUpdate (look, lpDBox->obj) ;
}

/************************************************************/

static void  dgramPick (int box)
{
	DGRAM_BOX  *lpDBox ;
	OBJ obj=0;
	BOOL    bLabelObj = FALSE, bNodeObj  = FALSE ;
	LOOKGET("dgramPick") ;

	if(box == look->messageBox)
	{
		graphPostBuffer (look->messageText) ; 
		return ;
	}

	if (!box)
	{ 
		if ( look->activebox > 0)
		{
			GETOBJ(look,look->activebox,lpDBox)
			obj = lpDBox->obj ;

		        if(obj) labelColour(look,look->activebox, NORMAL_MODE) ;
		}
		look->activebox = 0 ;
		return ;
	}

	if((lpDBox = arrp(look->boxes,box,DGRAM_BOX)))
	{
		obj = lpDBox->obj ;
		if(lpDBox->flags & LABELBOX)
			bLabelObj = TRUE ;
		else if(lpDBox->flags & NODEBOX)
		{
			bNodeObj  = TRUE ;
			box++ ;
		}
		if(bLabelObj || bNodeObj)
		{
			KEY labelText ;
			if(bsGetKey (obj, str2tag ("Label"), &labelText))
				graphPostBuffer (name(labelText)) ; 
			else
				graphPostBuffer (name (obj->root->key)) ; 
		}
	}

	if (box == look->activebox)         /* a second hit - follow it */
	{  
		if(bLabelObj) 
			/* Display associated data */
		  	dgramDisplayNode(look, obj, FULL_DISPLAY) ;
		else if (bNodeObj)
			dgramHide(TOGGLE_NODE) ;
		else
			selectActive(look, 0) ;
	}
	else  /* box != look->activebox; change now -- maybe */
	{
		if ( look->activebox > 0)
		{
			GETOBJ(look,look->activebox,lpDBox)
		        if(lpDBox->obj) labelColour(look,look->activebox, NORMAL_MODE) ;
		}

		if(bLabelObj || bNodeObj)
			selectActive(look, box) ;
		else
			selectActive(look, 0) ;
	}
}
/************************************************************/

static void  dgramScroll  (double x, double y)
{
	graphGoto((float)x,(float)y) ;
}

/************************************************************/

static void  dgramKbd (int k)
{
	DGRAM_BOX  *lpDBox ;
	OBJ obj = 0, nextObj = 0 ;
	KEY nextNode = 0 ;
	LOOKGET("dgramKbd") ;

	if(k == ESCAPE_KEY)   /* change colour mode? */
	{
		look->bColour = !look->bColour ;
		dgramDraw (look) ;
		return ;
	}

	if(k == TAB_KEY)   /* change label mode? */
	{
		look->bDescriptiveLabels = !look->bDescriptiveLabels ;
		dgramDraw (look) ;
		return ;
	}

	if(k == HOME_KEY)   /* goto first leaf */
	{
		selectActive(look, look->firstbox) ;
		return ;
	}
	else if(k == END_KEY)    /* goto last leaf */
	{
		selectActive(look, look->lastbox) ;
		return ;
	}
	else if (!look->activebox)
		return ; /* All remaining operations assume an active node? */

	else /* look->activebox > 0? */
	{
		/* Current box should be a valid ?TreeNode object for these operations */
		GETOBJ(look,look->activebox,lpDBox)
		if(!(obj = lpDBox->obj))
		{
			selectActive(look, 0) ;
			return ;
		}

		/* change to new ?TreeNode object? */
		switch (k)
		{
			case RETURN_KEY:
				/* Display associated data */
				dgramDisplayNode(look, obj,FULL_DISPLAY) ;
				return ;
			case INSERT_KEY:
				dgramUpdate (look, obj) ;
				return ;
			case DELETE_KEY:
				dgramHide(TOGGLE_NODE) ;
				return ;
			case LEFT_KEY :  /* goto parent of node */
				if(!( bsGetKey(obj,str2tag("Parent"),&nextNode) &&
				     (nextObj = dgramGet(look, nextNode))))
					nextObj = obj ;
				break ;
			case RIGHT_KEY : /* goto first child of subtree, expand if hidden */
				if(!(bsGetKey(obj,str2tag("Child"),&nextNode) &&
				    (nextObj = dgramGet(look, nextNode))))
					nextObj = obj ;
				else dgramHide(SHOW_NODE) ;
				break ;
			case DOWN_KEY : /* goto next 1st level child of current subtree, expand if hidden */
				{
					NEIGHBOURS nb; 
					GETOBJ(look,look->activebox,lpDBox)
					if(lpDBox->nbIdx>=0)
					{
						nb = arr( look->nodeNeighbours,
							  lpDBox->nbIdx,
							  NEIGHBOURS) ;
						selectActive(look, 0) ;
						look->activebox = nb.next ; 
						GETOBJ(look,look->activebox,lpDBox)
						nextObj = lpDBox->obj ;
					}
					else nextObj = obj ; /* stay put?*/
				}
				break ;
			case UP_KEY :   /* goto previous 1st level child of current subtree, expand if hidden */
				{
					NEIGHBOURS nb; 
					GETOBJ(look,look->activebox,lpDBox)
					if(lpDBox->nbIdx>=0)
					{
						nb = arr( look->nodeNeighbours,
							  lpDBox->nbIdx,
							  NEIGHBOURS) ;
						selectActive(look, 0) ;
						look->activebox =  nb.previous ;
						GETOBJ(look,look->activebox,lpDBox)
						nextObj = lpDBox->obj ;
					}
					else nextObj = obj ; /* stay put?*/
				}
				break ;
			case SPACE_KEY :  /* goto next leaf at same level */
				GETOBJ(look,look->activebox,lpDBox)
				if(lpDBox->leafPos >= 0)
				{
					look->leafPos = lpDBox->leafPos ;
				        if( ++look->leafPos >= arrayMax(look->leaves))
						look->leafPos = 0 ;
					
					selectActive(look, 0) ;
					look->activebox = arr(look->leaves,look->leafPos,int) ;
					GETOBJ(look,look->activebox,lpDBox)
					nextObj = lpDBox->obj ;
				}
				else nextObj = obj ; /* stay put?*/
				break ;
			case BACKSPACE_KEY :    /* goto previous sib at same level */
				GETOBJ(look,look->activebox,lpDBox)
				if(lpDBox->leafPos >= 0)
				{
					look->leafPos = lpDBox->leafPos ;
				        if( --look->leafPos < 0 )
						look->leafPos = arrayMax(look->leaves)-1 ;
					
					selectActive(look, 0) ;
					look->activebox = arr(look->leaves,look->leafPos,int) ;
					GETOBJ(look,look->activebox,lpDBox)
					nextObj = lpDBox->obj ;
				}
				else nextObj = obj ; /* stay put?*/
				break ;
			default:
				return ; /* nop for non-operational key */
		}
	}
	if(nextObj)
	{
		char *vp ; 
		int box ;
		if(!assFind(look->obj2box,nextObj,&vp)) 
			return ;
		box = (int)(vp -(char *)0) ;

		selectActive(look, box) ;
	}
	else selectActive(look, 0) ;
}

/************************************************************/
/* Returns the simple node label for a given ?TreeNode  */
static char *leafName(OBJ obj, int treeType)
{
	KEY key = 0 ;
	char *label = 0 ;

	if(!obj) return 0 ;

	switch(treeType)
	{
		case TAXONOMIC_TREE:
			if(bsGetKey(obj,str2tag("Taxon"),&key))
				label = name(key) ; /* Taxon name? */
			break ;
		case DNA_TREE:
			if(bsGetKey(obj,str2tag("DNA"),&key))
				label = name(key) ; /* DNA sequence name? */
			break ;
		case PROTEIN_TREE:
			if(bsGetKey(obj,str2tag("Protein"),&key))
				label = name(key) ; /* Protein name? */
			break ;
		case CELL_LINEAGE:
			if(bsGetKey(obj,str2tag("Cell"),&key))
				label = name(key) ; /* Cell name? */
			break ;
		default:  
			break ;
	}
	if(label == 0)
	{
		if(bsFindTag(obj,_Contains) &&
		   bsFindTag(obj,_bsRight) /* to the tag2 tag */ &&
		   bsFindTag(obj,_bsRight) /* to the tag2 value */ &&
		   bsGetKey(obj,_bsHere,&key))
			label = name(key) ; /* first "tag2" Contains object name? */
		else if(obj->root)
			label = name(obj->root->key) ;
	}
	return label ;
}
/************************************************************/

/* Returns a descriptive label for a given ?TreeNode primary
   object (?Taxon, ?Sequence, ?Protein or ?Cell) if available;
   if not "long" label is found, the primary object's name is returned;
   if no primary object, return NULL */
static char *longLabel(LOOK look, OBJ obj,int treeType )
{
	KEY key ;
	char *cp ;
	static Stack label = 0 ;

	if (!label)
		label = stackCreate (BUF_WIDTH) ;
	stackClear (label) ;
	catText(label,"") ;

	switch(treeType)
	{
		case TAXONOMIC_TREE:
			if(bsGetKey(obj,str2tag("Taxon"),&key))
			{
				OBJ taxonObj ;
				if ((taxonObj = dgramGet(look, key)))
				{
					if(*(cp = commonName(look,taxonObj,TAXONOMIC_TREE)))
						catText(label,cp) ;
					else
					        /* default back to taxon key name */
						catText(label,name(key)) ;
				}
			}
			break ;
		case DNA_TREE:
			if(bsGetKey(obj,str2tag("DNA"),&key))
			{
				OBJ SequenceObj ;
				if ((SequenceObj = dgramGet(look,key)))
				{
					if (*(cp = commonName(look,SequenceObj,DNA_TREE)))
						catText(label,cp) ;
					else
					        /* default back to DNA sequence name */
						catText(label,name(key)) ;
					/* specify species if known */
					if(*(cp = species(look, SequenceObj)))
						catText(label,messprintf(" {%s}",cp)) ;
				}
			}
			break ;
		case PROTEIN_TREE:
			if(bsGetKey(obj,str2tag("Protein"),&key))
			{
				OBJ proteinObj ;
				if ((proteinObj = dgramGet(look, key)))
				{
					if (*(cp = commonName(look,proteinObj,PROTEIN_TREE)))
						catText(label,cp) ;
					else
					        /* default back to Protein sequence name */
						catText(label,name(key)) ;
					/* specify species if known */
					if(*(cp = species(look, proteinObj)))
						catText(label,messprintf(" {%s}",cp)) ;
				}
			}
			break ;
		case CELL_LINEAGE:
			if(bsGetKey(obj,str2tag("Cell"),&key))
			{
				OBJ cellObj ;
				if ((cellObj = dgramGet(look, key)))
				{
					if(*(cp = commonName(look,cellObj,CELL_LINEAGE)))
						catText(label,cp) ;
					else
					        /* default back to Cell key name */
						catText(label,name(key)) ;
				}
			}
			break ;
		default: ;
	}
	/* Look for the first "Contains" tag2 object,
	    if present, for a common name? */
	if( !*stackText(label,0)     && /* is the label empty? */
	    bsFindTag(obj,_Contains) && /* but the ?TreeNode "Contains" something? */
	    bsFindTag(obj,_bsRight)  && /* get to the tag2 tag */
	    bsFindTag(obj,_bsRight)  && /* get to the tag2 object */
	    bsGetKey(obj,_bsHere,&key))
	{
		OBJ containsObj ;
		if ((containsObj = dgramGet(look, key)))
		{
			if(*(cp = commonName(look,containsObj,UNKNOWN_TREE)))
				catText(label,cp) ;
			else
				/* default back to "Contains" Object key name */
				catText(label,name(key)) ;
		}
	}
	return strnew(stackText(label,0),look->handle) ; 
}
/************************************************************/

static void biblioButton (void) 
{
  LOOKGET ("biblioButton") ;
  biblioKey (look->key) ;
}

/* This routine finds the number of tree leafs and upper
   (non-zero) bound of tree branch lengths */
static BOOL measureTree(LOOK look, KEY treeNode, int *lpnLeaves, float *max)
{
	OBJ obj ;
	KEY childNode ;
	float distance = 0.0 ;

	if (!(obj = dgramGet(look, treeNode))) return FALSE ;

	/* If !RootNode, check the current node for a distance */
	if ( bsFindTag (obj, str2tag ("Parent")))
	{
		if (!bsGetData (obj, str2tag ("Distance"), _Float, &distance))
			distance = 0.0 ; /* if no distance found, assume zero */
		if(distance < 0.0) distance = 0.0 ; /* reject negative lengths */
	}
	/* Start with current node distances setting the range...*/
	*max = distance ;

	/* then, visit the children recursively */
	if(bsGetKey(obj,str2tag("Child"),&childNode))
	{
		do
		{
			float cmax = 0.0 ;
			int   nChildLeaves = 0 ;
			if( !measureTree(look, childNode, &nChildLeaves, &cmax))
				return FALSE ; 
			if(cmax > *max) *max = cmax ;
			*lpnLeaves += nChildLeaves ;
		}
		while(bsGetKey(obj,_bsDown,&childNode)) ;
	}
	else /* isa Leaf -- count it */
		*lpnLeaves = 1 ;

	return TRUE ;
}

/************************************************************/

/* display character sizes */
static int tx, ty ;
static float cw, ch ; 

/* treeNode may be a root, interior or leaf "?TreeNode" object  */
#define VOFFSET 2.0 /* double spacing...*/

/* drawSubTree(): a recursive routine */
static int drawSubTree (LOOK look, 
			KEY treeNode, 
			float *x, float *y, 
			int *labelBox,
			float theta )
/* theta, lower, upper ignored except for unrooted "star" trees */
{
	DGRAM_BOX  dBox, *lpDBox ;
	OBJ	obj ;
	KEY	childNode ;
	int	oldActiveColour = -1 , oldGraphColour = -1, 
		nChildren = 0, nodeBox = 0 ;
	float   distance,lineLength, 
		nextNodeX = *x, nextNodeY = *y, nextNodeR, branchEnd,
		yAlignTop = *y, yAlignBottom = *y, ddy=0;
	char   *label = 0 ;
	BOOL	bIsaRoot, 
		bIsaLeaf, 
		bResetBootstraps,
		bHidden  = FALSE, 
		bFirstChild = TRUE ;

	if (!(obj = dgramGet(look, treeNode)))
	{
		messerror("drawSubTree can't access tree node object?") ;
		return 0 ;
	}
	bIsaRoot = !bsFindTag (obj, str2tag ("Parent")) ;
	bIsaLeaf = !bsFindTag(obj,str2tag("Child")) ;

	/* One can only hide a whole subtree's bootstraps? */
	if((bResetBootstraps = !look->bHideBootstraps))
		look->bHideBootstraps = bsFindTag(obj,str2tag("Hide_Bootstraps")) ;

	/* The graphBox containing the entire subtree... */
	graphBoxStart() ;

	/* If "colour" mode set and subtree colour is specified... */
	if( look->bColour && bsFindTag(obj, str2tag ("Colour")))
	{
		oldActiveColour = look->activeColour ;
		look->activeColour = controlGetColour(obj) ;
		oldGraphColour = graphColor(look->activeColour) ;
	}

	if (bIsaRoot)
		nextNodeR = branchEnd = *x ;
	else
	{
		if(look->distFactor == 0.0) /* dimensionless tree? */ 
			lineLength = 2.0 ; /* esthetically pleasing but meaningless branch length --  */
		else if (bsGetData (obj, str2tag ("Distance"), _Float, &distance))
		{
			lineLength = distance*look->distFactor ;
		}
		else lineLength = 1.0 ; /* small branch length, for cosmetic purposes only */
		if(lineLength<1.0) lineLength = 1.0 ; /* impose lower limit for cosmetic purposes */

		nextNodeR = branchEnd = *x+lineLength ;
	}

	if ( !look->bDescriptiveLabels || 
	     !*(label = longLabel(look,obj,look->treeType)) )
	{
		KEY labelText ;
		if (bsGetKey (obj, str2tag ("Label"), &labelText))
			label = (char*)strnew(name(labelText),look->handle) ;
		else
		{
		    if( look->treeType == TAXONOMIC_TREE ||
			look->treeType == CELL_LINEAGE   ||
			bsFindTag(obj, str2tag ("Leaf")) ||
			bsFindTag(obj, _Contains) ) 
				label = leafName(obj,look->treeType) ;
		    else /* ignore interior node names in other tree types */
				label = 0 ;
		}
	}

	if(look->alignment == UNROOTED)
	{
		nextNodeX += nextNodeR*cos(theta) ;
		nextNodeY += nextNodeR*sin(theta);
	}
	else
	{
		if(label) 
			nextNodeR += (1.5+strlen(label)*cw) ;

		nextNodeX = nextNodeR ;
		ddy = nextNodeY ;
	}

	/* If the subtree is visible, draw it */
	if (assFind (look->key2hide,assVoid(treeNode), &bHidden) && !bHidden)
	{
		int i = 0 ,j, childBox;
		Array subTreeBoxes = arrayCreate(8,int) ; 

		/* draw children recursively */
		if(!bIsaLeaf && bsGetKey(obj,str2tag("Child"),&childNode))
		{
			int nSubChildren ; 
			do 
			{ 
				graphColor(look->activeColour) ;
				if(look->alignment != UNROOTED) 
					nextNodeY = ddy ;

				nSubChildren = 
					drawSubTree(look,
						    childNode,
						    &nextNodeX,&nextNodeY, 
						    &childBox,
						    theta ) ; /* need to do something different here...*/
				array(subTreeBoxes,i++,int) = childBox ;
				if(bFirstChild)
				{
					yAlignTop = nextNodeY ; bFirstChild = FALSE ;
				}
				nChildren += nSubChildren ;
				if(look->alignment == UNROOTED)
					/* a simplistic version M1 of unrooted trees: fixed angle phi rotation */
					theta += nSubChildren*look->phi ; 
				else 
					ddy += nSubChildren*VOFFSET ;
			} 
			while(bsGetKey(obj,_bsDown,&childNode)) ;
			yAlignBottom = nextNodeY ;

		}
		else nChildren = 1 ; 

		if(look->alignment != UNROOTED &&
		   yAlignTop != yAlignBottom ) 
			graphLine((float)nextNodeX,(float)yAlignTop+0.5*ch,
				  (float)nextNodeX,(float)yAlignBottom+0.5*ch) ;

		/* Store list of neighbours for each 
		   first degree children of interior nodes */
		if(!bIsaLeaf)
		{
			for(j = 0; j < i; ++j)
			{
				NEIGHBOURS nb ;
				nb.previous = arr(subTreeBoxes,(j?j-1:i-1),int) ;
				nb.next = arr(subTreeBoxes,((j==i-1)?0:j+1),int) ;
				lpDBox = arrp(look->boxes,arr(subTreeBoxes,j,int),DGRAM_BOX) ;
				lpDBox->nbIdx = look->nbSeq++ ;
				array( look->nodeNeighbours,lpDBox->nbIdx,NEIGHBOURS) = nb ;
			}
		}
		arrayDestroy(subTreeBoxes) ;
	}
	else  /* otherwise, measure #leaves in hidden subtree but treat as a leaf? */
	{
		float dummy = 0.0 ;
		measureTree(look, treeNode, &nChildren, &dummy) ;
	}

	if(look->alignment != UNROOTED)
	{
		if( yAlignTop != yAlignBottom )
			switch(look->alignment)
			{
			case MIDDLE:
				*y = 0.5*(yAlignTop+yAlignBottom) ;
				break;
			case BOTTOM:
				*y = yAlignBottom ;
				break;
			case TOP:
				*y = yAlignTop ;
			default: /* do nothing - y stays the same */
				break;
			}
		else *y = yAlignTop ; /* degenerate case? */
	}

	/* Create the current ?TreeNode graphBox */
	nodeBox = graphBoxStart() ;

	graphColor(look->activeColour) ;
	dBox.leafPos = dBox.nbIdx = -1 ;
	dBox.flags   = NODEBOX ;
	dBox.colour  = look->activeColour ; 
	dBox.obj     = obj ;
	array(look->boxes,nodeBox,DGRAM_BOX) = dBox ;

	/* Draw this ?TreeNode's branch "distance" line & node bootstrap value? */
	if (!bIsaRoot) 
	{
		if(look->alignment == UNROOTED)
		{
			float nfac ;
			/* Distance line... */
			graphLine(*x,*y,nextNodeX,nextNodeY) ;

			/* Bootstrap values, if displayable... */
			if( !bHidden && !look->bHideBootstraps &&
			    bsGetData (obj, str2tag ("Bootstrap"), _Float, &nfac))
			{
				float fBootValue = nfac / look->bootstrapFactor ;
				if(fBootValue<= 1.0)
					graphText(strnew(messprintf("%1.1f%%",fBootValue*100.0),look->handle),
						  (*x+nextNodeX)/2,(*y+nextNodeY)/2 ) ;
				else
					graphText(strnew(messprintf("%1.1f",fBootValue),look->handle),
						  (*x+nextNodeX)/2,(*y+nextNodeY)/2 ) ;
			}

		}
		else
		{
			float nfac ;
			/* Distance line... */
			graphLine(*x,*y+0.5*ch,branchEnd,*y+0.5*ch) ;

			/* Bootstrap values, if displayable... */
			if( !bHidden && 
			    bsGetData (obj, str2tag ("Bootstrap"), _Float, &nfac) &&
			    (!look->bHideBootstraps ||
			      bsFindTag(obj,str2tag("Show_Bootstrap"))) 
			   )
			{
				float fBootValue = nfac / look->bootstrapFactor ;
				if(fBootValue<= 1.0)
					graphText(strnew(messprintf("%1.1f%%",fBootValue*100.0),look->handle),
						   branchEnd+1.0*cw,*y ) ;
				else
					graphText(strnew(messprintf("%1.1f",fBootValue),look->handle),
						   branchEnd+1.0*cw,*y ) ;
			}
		}
	}

	if(look->alignment == UNROOTED)
		graphCircle(nextNodeX+0.5,nextNodeY,0.5) ;
	else
		graphCircle(branchEnd,*y+0.5,0.5) ;

	if(!bIsaLeaf && bHidden)
	{
		Array temp = arrayCreate(6, float) ;
		float px, py ;
		if(look->alignment == UNROOTED)
		{
			px = nextNodeX ; py = nextNodeY ;
		}
		else
		{
			px = nextNodeX+1.0 ; py = *y+0.5*ch ;
		}
		graphLine(px,py,px+2.0,py) ;
		px += 2.0 ;
		array(temp, 0, float) = px ;
		array(temp, 1, float) = py-0.35;
		array(temp, 2, float) = px+0.7;
		array(temp, 3, float) = py;
		array(temp, 4, float) = px;
		array(temp, 5, float) = py+0.35;
		graphPolygon(temp);
		px += 2.0 ;
		graphText(strnew(messprintf("%d",nChildren),
				 look->handle),
			  px,py-0.35) ;
	}
	graphBoxEnd () ; /* nodeBox */

	/* Create the current ?TreeNode label Box */
	*labelBox = graphBoxStart() ;
	graphBoxInfo (*labelBox, treeNode, label); /* WWW */

	if (bIsaRoot)
	{
		NEIGHBOURS nb ;
		nb.previous =  nb.next =  *labelBox ;
		dBox.nbIdx = look->nbSeq++ ;
		array(look->nodeNeighbours,dBox.nbIdx,NEIGHBOURS) = nb ;
		look->rootbox = *labelBox ;
	}
	else dBox.nbIdx = -1 ;

	if( bIsaLeaf || bHidden) /* pay special attention to leaves/hidden nodes */
	{
		if(!look->firstbox)
			look->firstbox = look->lastbox = *labelBox ;
		else
			look->lastbox = *labelBox ;

		array (look->leaves, ++(look->leafPos), int) = *labelBox ;
		dBox.leafPos = look->leafPos ;
	}
	else dBox.leafPos = -1 ;

	/* Keep track of ?TreeNode/graph *labelBox one-to-one mapping */
	dBox.flags   = LABELBOX ;
	if(look->activeColour) 
		dBox.colour = look->activeColour ;
	else 
		dBox.colour = 0 ; 
	dBox.obj     = obj ;
	array(look->boxes,*labelBox,DGRAM_BOX) = dBox ;

	assRemove (look->obj2box, obj) ;
	assInsert (look->obj2box, obj, assVoid(*labelBox)) ; 

	/* Draw the node label, if present */
	if(look->alignment == UNROOTED)
	{
#if defined(WIN32_TREE) /* only WinAce has graphTextAlign() for the moment? */
		/* Cosmetic label alignment...*/
		if(theta >= 0.0 && theta <= 0.5*pi)
			graphTextAlign(LEFT|BOTTOM) ;
		else if(theta > 0.5*pi && theta <= pi)
			graphTextAlign(RIGHT|BOTTOM) ;
		else if(theta > pi && theta <= 1.5*pi)
			graphTextAlign(RIGHT|TOP) ;
		else graphTextAlign(LEFT|TOP) ;
#endif
		if(label)
			graphText (label,nextNodeX+0.5,nextNodeY) ;
		else
			graphText ("",nextNodeX+0.5,nextNodeY) ;
	}
	else
	{
		if(label)
			graphText (label,branchEnd+1.0,*y) ;
		else
			graphText ("",branchEnd+1.0,*y) ;
	}

	graphBoxEnd () ; /* *labelBox */
	labelColour(look, *labelBox, NORMAL_MODE) ;

	if(oldActiveColour>=0) look->activeColour = oldActiveColour;
	if(oldGraphColour>=0)   graphColor(oldGraphColour) ;

	graphBoxEnd() ; /* subtree graph box */

	if(look->xmax < nextNodeX) look->xmax = nextNodeX ;
	if(look->alignment == UNROOTED)
	{
		if(look->ymax < nextNodeY) look->ymax = nextNodeY ;
	}
	else
	{
		if(look->ymax < ddy) look->ymax = ddy ;
	}
	if(bResetBootstraps)
		look->bHideBootstraps = FALSE ;

	if(bIsaLeaf || bHidden)
		return 1 ;
	else 
		return nChildren ;
}

static void  dgramDraw (LOOK look)
{
	int nLeaves = 0 ;
	float mby, tdmax, x, y ;
	KEY root ;
	int box ;
	OBJ obj ;
	char brMax[30], brNormal[30], numLeaves[16], scale[24]  ;

	graphClear () ;
	graphTextFormat(FIXED_WIDTH) ;

	look->boxes = arrayReCreate(look->boxes,128,DGRAM_BOX) ;
	look->rootbox = look->activebox = look->firstbox = look->lastbox = 0 ;
	look->nodeNeighbours = arrayReCreate(look->nodeNeighbours, 64, NEIGHBOURS) ;
	look->nbSeq = 0 ;

	assClear (look->obj2box) ;
        arrayMax(look->leaves) = 0 ;

	look->leafPos   = -1 ;
	look->firstbox = look->lastbox = 0 ;

	if( look->treeType == TAXONOMIC_TREE || 
	    look->treeType == CELL_LINEAGE )
		mby = 3.0 ;
	else
		mby = 5.0 ;

	/* Draw the dendrogram here? */
	if(look->bHeader)
	{
		look->xmax = x = 2.0 ; 
		look->ymax = y = mby+2.0 ;
	}
	else
	{
		look->xmax = x = 2.0 ; 
		look->ymax = y = mby ;
	}
	/* Root Tree Node? */
	if( bsGetKey(look->objet, str2tag("Root"),&root) )
	{ 
		float theta, arc ;

		/* Measure maximum branch distance in the tree */
		if( measureTree(look, root,&nLeaves,&tdmax))
		{
			/* valid tdmax -- may be 0.0 
			   if tree is dimensionless (i.e. taxonomic) */
			if(tdmax > 0.0)
				look->distFactor = look->normalFactor/tdmax ;
			else
				look->distFactor = 0.0 ;
		}
		else
		{
			messerror("dgramDraw(): measureTree() failed?") ;
			return ;
		}

		graphTextInfo (&tx,&ty,&cw,&ch) ; /* to calibrate text window sizes...*/

		/* Calculate some parameters for unrooted trees,
		   converting degrees to radians */
		arc = pi*(look->spread/180.0) ;
		theta = pi*(look->angle/180.0) - arc/2.0 ;  
		look->phi = arc/nLeaves ;

		/* The main event! */
		drawSubTree (look, root, &x, &y, &box, theta ) ;

#if defined(WIN32_TREE)
		graphTextAlign(LEFT|TOP) ;
#endif
		/* Activate first box as the active box -- maybe */
		look->leafPos = 0 ;
	}
	else
	{
		messout("dgramDraw(): no \"Root ?TreeNode\" object found in ?Tree object? Can't draw tree!") ;
		return ;
	}

	if (look->xmax < BUF_WIDTH) look->xmax = BUF_WIDTH ;

	if(look->bHeader)
	{
		graphTextFormat(PLAIN_FORMAT) ;
		graphButton("Quit", graphDestroy, BUF_WIDTH - 4.5, 1.0) ;
		graphButton("Help", dgramHelp, BUF_WIDTH - 10.0, 1.0) ;

		if (biblioKeyPossible(look->key))
			graphButton("Biblio", biblioButton, BUF_WIDTH - 17.5, 1.0) ;

		if( !(	look->treeType == TAXONOMIC_TREE  || 
			look->treeType == CELL_LINEAGE ) )
		{
			strncpy(numLeaves,messprintf("# Leaves: %d",nLeaves),15) ;
			numLeaves[15] = '\0' ;
			graphText(numLeaves,1.0,1.4) ;
			if(tdmax > 0.0)
			{
				strncpy(scale,
					messprintf("Scale: %.1f = %.3f",
					           look->normalFactor/10.0,
						   tdmax/10.0), 23) ;
				scale[23] = '\0' ;
				graphText(scale,BUF_WIDTH - 19.0,2.8) ;
				graphLine(BUF_WIDTH - 15.0,4.2,BUF_WIDTH - 15.0,4.6) ;
				graphLine(BUF_WIDTH - 15.0,4.4,BUF_WIDTH - 5.0,4.4) ;
				graphLine(BUF_WIDTH - 5.0,4.2,BUF_WIDTH - 5.0,4.6) ;
				strncpy(brMax,messprintf   ("Max. Branch Length: %.3f",tdmax),29) ;
				brMax[29] = '\0' ;
				graphText(brMax,16.0,1.4) ;
				strncpy(brNormal,messprintf("Display Normalization: %.1f",look->normalFactor),29) ;
				brNormal[29] = '\0' ;
				graphText(brNormal,16.0,2.8) ;
			}
		}
		look->messageBox = graphBoxStart () ;
		graphTextPtr (look->messageText, 1.0, mby, BUF_WIDTH-1) ;
		graphBoxEnd () ;
		graphBoxDraw (look->messageBox, BLACK, PALECYAN) ;
	}
	else
	{
		look->messageBox = -1 ;
		graphTextFormat(PLAIN_FORMAT) ;
		if( !(	look->treeType == TAXONOMIC_TREE || 
			look->treeType == CELL_LINEAGE ) && 
			tdmax > 0.0)
		{
			strncpy(scale,
				messprintf("%.3f", tdmax/10.0), 23) ;
			scale[23] = '\0' ;
			graphText(scale,9.0,1.2) ;
			graphLine(6.0,2.4,6.0,2.8) ;
			graphLine(6.0,2.6,16.0,2.6) ;
			graphLine(16.0,2.4,16.0,2.8) ;
		}
	}
	graphTextBounds ((int)look->xmax+10,(int)look->ymax+5) ;
	graphRedraw() ;

	/* Activate the first object */
	if( look->activeKey && 
	   (obj = dgramGet(look, look->activeKey)) && 
	    assFind(look->obj2box,obj,&box))
		selectActive(look, box) ;

	else
	{
		DGRAM_BOX  *lpDBox ;
		GETOBJ(look,look->firstbox,lpDBox)
		if(lpDBox->obj)
			selectActive(look, look->firstbox) ;
		else    
			selectActive(look, 0) ;
	}
}

static void dgramRedraw (void)
{
	LOOKGET("dgramRedraw") ;
	dgramDraw(look) ;
}

/********************************/
/* Dendrogram Tree dumping code */
/********************************/
static char	dirName[DIR_BUFFER_SIZE] = "", 
		fileName[MAXPATHLEN] = "";

typedef void  (*DgramTreeDumpFn)(LOOK, OBJ, FILE *) ;

static void  dgramTreeAceDump(LOOK look, OBJ treeNodeObj, FILE *fil)
{
	KEY key ;

	/* dump the current ?TreeNode record...*/
	dumpKey (treeNodeObj->root->key, fil, 0) ;

	/* I do a shallow dump of these keys - users should 
	   dump these "contained" classes separately if required */
	if(bsGetKey(treeNodeObj,str2tag("Taxon"),&key))
	{
		dumpKey (key, fil, 0) ;
	}
	if(bsGetKey(treeNodeObj,str2tag("DNA"),&key))
	{
		dumpKey (key, fil, 0) ;
	}
	if(bsGetKey(treeNodeObj,str2tag("Protein"),&key))
	{
		dumpKey (key, fil, 0) ;
	}
	/* Dump child subtrees recursively */
	if(bsGetKey(treeNodeObj,str2tag("Child"),&key))
	{
		BSMARK childMark = 0 ;
		do 
		{ 
			OBJ subTreeObj ;
			childMark = bsMark(treeNodeObj,childMark) ;
			if ((subTreeObj = dgramGet(look, key)))
				dgramTreeAceDump(look, subTreeObj, fil) ;
			else 
				messerror("dgramTreeAceDump() can't access root tree node object?") ;
			bsGoto(treeNodeObj,childMark) ;
		} 
		while(bsGetKey(treeNodeObj,_bsDown,&key)) ;
		if(!bsMarkFree(childMark)) 
			messcrash("dgramTreeAceDump(): childMark deletion error?") ;
	}
}


static void dgramTreeLeafDump(LOOK look, OBJ treeNodeObj, FILE *fil)
{
	KEY key ;
	/* Dump leaves of child subtrees recursively */
	if(bsGetKey(treeNodeObj,str2tag("Child"),&key))
	{
		BSMARK childMark = 0 ;
		do 
		{ 
			OBJ subTreeObj ;
			childMark = bsMark(treeNodeObj,childMark) ;
			if ((subTreeObj = dgramGet(look, key)))
				dgramTreeLeafDump(look, subTreeObj, fil) ;
			else 
				messerror("dgramTreeLeafDump() can't access child tree node object?") ;
			bsGoto(treeNodeObj,childMark) ;
		} 
		while(bsGetKey(treeNodeObj,_bsDown,&key)) ;
		if(!bsMarkFree(childMark)) 
			messcrash("dgramTreeLeafDump(): childMark deletion error?") ;
	}
	else  /* found a leaf? dump it? */
	{
		char    *shortLabel = 0, 
			*descLabel = longLabel(look, treeNodeObj,look->treeType) ;
		KEY labelText ;

		if (bsGetKey (treeNodeObj, str2tag ("Label"), &labelText))
			shortLabel = (char*)strnew(name(labelText),look->handle) ;
		else
		{
		    if( look->treeType == TAXONOMIC_TREE ||
			look->treeType == CELL_LINEAGE   ||
			bsFindTag(treeNodeObj, str2tag ("Leaf")) ||
			bsFindTag(treeNodeObj, _Contains) ) 
				shortLabel = leafName(treeNodeObj,look->treeType) ;

		    else	/* Just take the ?TreeNode name */
				shortLabel = name(treeNodeObj->root->key) ; 
		}
		if(!*descLabel) descLabel = shortLabel ;

		fprintf(fil,
			"\"%s\",\"%s\"\n", 
			shortLabel, 
			descLabel) ;
	}
}

static void dgramTreeDump(DgramTreeDumpFn lpDgmDF, 
			  char *fext, 
			  char *header, 
			  BOOL bAceDump)
{
	FILE *fil = NULL ;
	DGRAM_BOX  *lpDBox ;
	OBJ obj = 0 ;
	KEY key ;
	char *cp ;
	BOOL bDumpTreeRecord = FALSE ;
	LOOKGET("dgramDump") ;

	/* First, try to retrieve obj of an active node? */
	if (  look->activebox && 
	     (lpDBox = arrp(look->boxes,look->activebox,DGRAM_BOX))
	    ) obj = lpDBox->obj ;

	/* If no active node, then attempt to  
	   retrieve the main ?Tree Root ?TreeNode obj */
	if(!obj && (obj = dgramGet(look, look->key))) 
	{
		if(!messQuery("Dumping the entire dendrogram ?Tree! Continue?"))
			goto abort ;
		
		/* ace dump the ?Tree record when indicated*/
		if(bAceDump) bDumpTreeRecord = TRUE ;

		/* ... then, retrieve the Root ?TreeNode */
		if( bsGetKey(obj, str2tag ("Root"), &key) )
			obj = dgramGet(look, key) ;
		else 
		{
			messerror("dgramDump() can't access ?Tree Root ?TreeNode ?") ;
			goto abort ;
		}
	}

	/* if object is still not retrieved at this point, then this is an error... */
	if(!obj)
	{
		messerror("dgramDump() can't access tree data?") ;
		goto abort ;
	}

	/* ...otherwise, by this point, obj should be set to the 
	   Root ?TreeNode from ?Tree or to the currently active ?TreeNode;
	   Open the output file using the root object name as the default file name */

	strncpy(fileName,name(obj->root->key),FIL_BUFFER_SIZE-1) ;
	fileName[FIL_BUFFER_SIZE-1] = '\0' ;
	for(cp=fileName;*cp;cp++) /* filter out forbidden characters */
		if(strchr( " \t\n\r\\/?\"\'|<>*:;.", *cp ))
			*cp = '_' ;

	if(!*dirName)
	  {
	    if (getenv ("ACEDB_DATA"))
	      strcpy (dirName, getenv("ACEDB_DATA")) ;
	    else if (filName ("rawdata", 0, "r"))
	      strcpy (dirName, filName ("rawdata", 0, "r")) ;
	  }

	if (!(fil = filqueryopen (dirName, fileName, fext, "w",
			    "Where do you wish to dump this tree data?")))
	{
		messout ("dgramDump() failed to open dump file") ;
		goto abort ;
	}

	/* Print the file output header */
	if(header && *header) fprintf(fil, "%s", header) ;

	/* dump the ?Tree record now! */
	if(bDumpTreeRecord) dumpKey (look->key, fil, 0) ;

	/* then call the caller defined dump function on the subtree */
	(lpDgmDF)(look, obj, fil) ;
	messout ("Tree data for %s written to %s\n", 
		 name(obj->root->key), fileName) ;
	goto finished ;

abort:
	if(fil) fprintf (fil,"\n// Dendrogram tree dump failed/aborted!\n") ;

finished:
	filclose (fil);
}

static void dgramAceDump(void)
{
	dgramTreeDump(dgramTreeAceDump,"ace",
		      "// Tree data dumped from dendrogram display\n\n",
		      TRUE /* ISA acedump */) ;
}

static void dgramLeafDump(void)
{
	dgramTreeDump(dgramTreeLeafDump,"txt",
		      "Format: Short label, Descriptive label {Species}\n\n",
		      FALSE /* NOT an acedump */) ;
}


static int nofseq ;
static void  dgramTreeFastaDump(LOOK look, OBJ treeNodeObj, FILE *fil)
{
	KEY key ;
	
	/* Dump any associated sequence? */
	switch(look->treeType)
	{
		case DNA_TREE:
			if(bsGetKey(treeNodeObj,str2tag("DNA"),&key))
			{
				nofseq++ ;
				dnaDumpFastAKey (key, fil, 0, 0) ;
			}
			break ;
		case PROTEIN_TREE:
			if(bsGetKey(treeNodeObj,str2tag("Protein"),&key))
			{
				nofseq++ ;
				pepDumpFastAKey (key, fil, 0, 0) ;
			}
			break ;
		default:  /* Otherwise, do nothing? */
			break ;
	}

	/* Dump child subtrees recursively */
	if(bsGetKey(treeNodeObj,str2tag("Child"),&key))
	{
		BSMARK childMark = 0 ;
		do 
		{ 
			OBJ subTreeObj ;
			childMark = bsMark(treeNodeObj,childMark) ;
			if ((subTreeObj = dgramGet(look, key)))
				dgramTreeFastaDump(look, subTreeObj, fil) ;
			else 
				messerror("dgramTreeFastaDump can't access child tree node object?") ;
			bsGoto(treeNodeObj,childMark) ;
		} 
		while(bsGetKey(treeNodeObj,_bsDown,&key)) ;
		if(!bsMarkFree(childMark)) 
			messcrash("dgramTreeFastaDump(): childMark deletion error?") ;
	}
}

static void dgramFastaDump(void)
{
	FILE *fil = 0 ;
	DGRAM_BOX  *lpDBox ;
	OBJ obj = 0 ;
	KEY key ;
	LOOKGET("dgramFastaDump") ;

	/* First, try to retrieve obj of an active node? */
	if (  look->activebox && 
	     (lpDBox = arrp(look->boxes,look->activebox,DGRAM_BOX))
	    ) obj = lpDBox->obj ;

	/* If no active node, then attempt to dump the main ?Tree
	   and retrieve Root ?TreeNode obj */
	if(!obj && (obj = dgramGet(look, look->key))) 
	{
		if(!messQuery("Dumping sequences for the entire dendrogram ?Tree! Continue?"))
			return ;

		/* Retrieve the Root ?TreeNode */
		if( bsGetKey(obj, str2tag ("Root"), &key) )
			obj = dgramGet(look, key) ;
		else 
		{
			messerror("dgramFastaDump() can't access ?Tree Root ?TreeNode ?") ;
			return ;
		}
	}
	
	/* By this point, obj should be set to the Root ?TreeNode from ?Tree 
	   or to the currently active ?TreeNode */
	if(obj)
	{
		switch(look->treeType)
		{
			case DNA_TREE:
				fil = dnaFileOpen () ; break ;
			case PROTEIN_TREE:
				fil = pepFileOpen () ; break ;
			default: break ;  /* Otherwise, do nothing? */
		}
		if (!fil)
		{
			messout ("dgramFastaDump() failed to open dump file") ;
			return ;
		}

		nofseq = 0 ;
		dgramTreeFastaDump(look, obj, fil) ;

		filclose (fil) ;

		messout ("%d FASTA sequence(s) written out for %s", 
			 nofseq, name(obj->root->key)) ;
	}
	else
		messerror("dgramFastaDump() can't access tree data?") ;
}

/************************************************************************ 
 * readTreeFile(): reads in New Hampshire phylogenetic trees into ACEDB
 ************************************************************************/

#define BUFSZ 64
static char	treeObjName[BUFSZ]  = "";
static FILE *ttfil = 0 ;	/* temporary .ace file for the tree */

static Stack parentNodes = 0 ;	/* Tree parsing stack, of current Parent nodes */

static int nextId = 0 ;
/* Types of TREENODES */
#define ROOT	  0
#define INTERIOR  1
#define LEAF	  2

typedef struct treenode {
	int id ;
	char name[BUFSZ+1] ;
} TREENODE, *LPTREENODE ;
static LPTREENODE currentParent ;
static TREENODE currentNode ;

static int parseTreeType = UNKNOWN_TREE ;

static BOOL skipSpaces (char **cp)
{
	char *start = *cp ;
	while (**cp && isspace((int)**cp)) { ++(*cp) ; }
	return (*cp == start ) ;
} /* skipSpaces */

static BOOL parseList(char **lpTreeList, char *name) ;
static int  leafNo = 0 ;

static void blankAtom()
{
	/* create a "blank" atom? */
	currentNode.id = ++nextId ;
	currentNode.name[0] = '\0' ;

	/* Create a tree node? */
	fprintf(ttfil,"TreeNode : \"%s%d\"\n",treeObjName,currentNode.id) ;
	fprintf(ttfil,"Id %d\n",currentNode.id) ;
	fprintf(ttfil,"Leaf\n") ;
	fprintf(ttfil,"Tree \"%s\"\n",treeObjName) ;
	if(currentParent->id>=0)
		fprintf(ttfil,"Parent \"%s%d\"\n",treeObjName,currentParent->id) ;
	fprintf(ttfil,"Map \"%s\" Position %d\n",treeObjName,++leafNo) ;
}

static BOOL isatomic(char c)
{
  if(isspace((int)c)) return FALSE ;
  switch(c)
    {
    case ',':
    case ':':
    case ';':
    case '(':
    case ')':
    case '[':
    case ']':
    case '{':
    case '}':
      return FALSE ;
    default:
      return TRUE ;	
    }
}

static char *parseName(char **cp)
{
	int pos = 0 ;
	static char buf[BUFSZ+1] ;

	skipSpaces(cp) ;
	while(pos < BUFSZ && **cp && ((!pos && isalpha((int)**cp))|| isatomic(**cp)))
	{ 
		if(**cp == '_')
			buf[pos++] = ' ' ; /* Substitute space for '_'*/
		else if(**cp == '\\')      /* escaped character follows */
			buf[pos++] = *(++(*cp)) ;
		else
			buf[pos++] = **cp ;
		++(*cp) ;
	}
	if(pos == BUFSZ ) /* atomic string overflow? skip extra's*/
		while(**cp && isatomic(**cp) ) ++(*cp) ;
	buf[pos] = '\0' ;
	return buf ; /* non-empty pointer but may be a blank string? */
}

static char rankName[BUFSZ+1] ;
static char taxonName[BUFSZ+1] ;
static  BOOL parseTaxon(char **lpTreeList)
{
	char *cp ;
		
	skipSpaces(lpTreeList) ;

	if( **lpTreeList != '{')
		return FALSE ;

	(*lpTreeList)++ ;
	if(*(cp = parseName(lpTreeList)))
	{
		strncpy(rankName,cp,BUFSZ) ;
		rankName[BUFSZ] = '\0' ;
	}
	else 
		rankName[0] = '\0' ;

	skipSpaces(lpTreeList) ;

	/* parse for species name? */
	if( **lpTreeList == ':')
	{
		(*lpTreeList)++ ; 
		skipSpaces(lpTreeList) ;
		if(*(cp = parseName(lpTreeList)))
		{
			strncpy(taxonName,cp,BUFSZ) ;
			taxonName[BUFSZ] = '\0' ;
		}
		else 
			taxonName[0] = '\0' ;
	}
	else taxonName[0] = '\0' ;

	skipSpaces(lpTreeList) ;
	if( **lpTreeList != '}')
	{
		messout("Unmatched \"}\" while parsing a taxon label?") ;
		return FALSE ;
	}
	(*lpTreeList)++ ;
	return TRUE ;
}

static BOOL printAtom(char *name)
{
	/* These atoms must have names...*/
	if( !name || !*name ) return FALSE ;

	currentNode.id = ++nextId ;
	strncpy(currentNode.name,name,BUFSZ) ;
	currentNode.name[BUFSZ] = '\0' ;

	/* Note that I label nodes anonymously, with sequence numbers, to ensure 
	   that they are always uniquely associated with a given "treeObjName" tree.
	   Because of UNIQUE Tree and Parent tag values in the models, failure to
	   observe this practice of anonymous ?TreeNode labelling may result in 
	   overwritten subtrees whenever new data is input with duplicate ?TreeNode names.
	   If this effect is desired, the user should be directed to edit the nodes manually */
	if(parseTreeType == TAXONOMIC_TREE || parseTreeType == CELL_LINEAGE)
		fprintf(ttfil,"TreeNode : \"%s\"\n",currentNode.name) ;
	else
		fprintf(ttfil,"TreeNode : \"%s%d\"\n",treeObjName, currentNode.id) ;

	fprintf(ttfil,"Id %d\n",currentNode.id) ;
	fprintf(ttfil,"Label \"%s\"\n",currentNode.name) ;
	fprintf(ttfil,"Leaf\n") ;
	fprintf(ttfil,"Tree \"%s\"\n",treeObjName) ;
	if( *currentParent->name && 
	   (parseTreeType == TAXONOMIC_TREE || parseTreeType == CELL_LINEAGE))
		fprintf(ttfil,"Parent \"%s\"\n",currentParent->name) ;
	else
		fprintf(ttfil,"Parent \"%s%d\"\n",treeObjName,currentParent->id) ;
	fprintf(ttfil,"Map \"%s\" Position %d\n",treeObjName,++leafNo) ;

	switch(parseTreeType)
	{
		case TAXONOMIC_TREE:
			/* do nothing: ?Taxon inserted elsewhere! */
			break ;
		case DNA_TREE:
			fprintf(ttfil,"Sequence \"%s\"\n",name) ;
			break ;
		case PROTEIN_TREE:
			fprintf(ttfil,"Protein \"%s\"\n",name) ;
			break ;
		case CELL_LINEAGE:
			fprintf(ttfil,"Cell \"%s\"\n",name) ;
			break ;
		default:  /* Otherwise, unknown tree type -- do nothing? */
			break ;
	}
	return TRUE ;
}

static char *parseNumber(char **cp)
{
	static char nbuf[BUFSZ+1] ;
	int pos = 0 ;

	skipSpaces(cp) ;

	if( **cp == '-' )
	{
		messout("Warning: negative branch length or bootstrap value encountered in tree file? Sign ignored!") ;
		++(*cp) ;
	}
	while(pos < BUFSZ && **cp && isdigit((int)**cp))
	{ 
		nbuf[pos++] = **cp ; ++(*cp) ; 
	}

	/* I'm going to force numbers to have at least one digit
	   before a decimal point, even if this is only a zero...*/
	if(pos && pos < BUFSZ && **cp == '.')
	{ 
		nbuf[pos++] = **cp ; ++(*cp) ; 
		while(pos < BUFSZ && **cp && isdigit((int)**cp))
		{ nbuf[pos++] = **cp ; ++(*cp) ; }
	}
	if(pos == BUFSZ ) /* number string overflow? skip extra's*/
		while(**cp && (isdigit((int)**cp)||**cp == '.') ) ;
	nbuf[pos] = '\0' ;

	return nbuf ;
}

/* I assume that I directly follow a current tree 
   member object treenode description in the .ace file stream */
static BOOL parseDistance(char **cp)
{
	char *nbuf = parseNumber(cp) ;

	if(*nbuf)
	{
		fprintf(ttfil, "Distance %s\n", nbuf) ;
		return TRUE ;
	}
	else return FALSE ;
}

/* I assume that I directly follow a current tree 
   member object treenode description in the .ace file stream */
static BOOL parseBootstrap(char **cp)
{
	char *nbuf ;

	skipSpaces(cp) ;

	if( **cp != '[' ) return FALSE ;

	++(*cp) ;

	nbuf = parseNumber(cp) ;

	skipSpaces(cp) ;

	if( **cp != ']' )
	{
		messout("Unmatched \"]\" while parsing a bootstrap value?") ;
		return FALSE ;
	}


	++(*cp) ;

	/* I assume that I directly follow a current ?TreeNode 
	   object description in the .ace file stream */
	if(*nbuf)
	{
		fprintf(ttfil, "Bootstrap %s\n", nbuf) ;
		return TRUE ;
	}
	else return FALSE ;
}

/* I assume that I directly follow a current tree 
   member object treenode description in the .ace file stream */
static void printTaxon(char *rank, char *taxon)
{
	if(*rank)
	{
		fprintf(ttfil, "Taxon \"%s\"\n\n",taxon) ;
		/* Taxon object too! */
		fprintf(ttfil, "Taxon : \"%s\"\n",taxon) ;
		fprintf(ttfil, "Taxonomy \"%s\"\n",taxon) ; /* treenode?*/
		fprintf(ttfil, "%s",rank) ;
		if(!strcasecmp(rank,"Species"))
		{
			fprintf(ttfil, " \"%s\"",taxon) ;
			fprintf(ttfil, "\n\n") ;
			fprintf(ttfil, "Species : \"%s\"\n",taxon) ;
			fprintf(ttfil, "Taxon \"%s\"",taxon) ;
		}
		fprintf(ttfil, "\n\n") ;

		/* Relate to the "standard" taxonomy tree?*/
		fprintf(ttfil, "TreeNode : \"%s\"\n",taxon) ;
		fprintf(ttfil, "Tree \"%s\"\n","NCBI Taxonomy") ; 
	}
}

static BOOL parseMember(char **cp)
{
	char   *name, label[BUFSZ+1], 
		rank[BUFSZ+1], taxon[BUFSZ+1] ;

	skipSpaces(cp) ;

	if( !**cp || **cp == ';')
	{
		messout("Premature tree termination (null or semi-colon) while parsing list member?") ;
		return FALSE ;
	}


	if( **cp == ',' )
	{
		blankAtom() ;
		/* Substitute blank to ensure that subsequent empty names 
		    preceeding the next ',' or ')' are scanned */
		**cp = ' ' ; 
		return TRUE ;
	}

	if( **cp == ')')
	{
		blankAtom() ;
		return TRUE ;
	}

	if( **cp == ':' )
	{
		++(*cp) ;
		blankAtom() ;
		if( !parseDistance(cp) )
		{
			messout("Expecting distance value after colon?") ;
			return FALSE ;
		}

		parseBootstrap(cp) ;

		skipSpaces(cp) ;
		if( **cp == ',' )
		{
			/* Substitute blank to ensure that subsequent empty names 
			    preceeding the next ',' or ')' are scanned */
			**cp = ' ' ; 
			return TRUE ;
		}
		if( **cp == ')') return TRUE ;
		messout("Missing comma or unmatched \")\" while parsing a distance value?") ;
		return FALSE ;
	}

	/* Parse for a {rank:taxon} label */
	if(parseTaxon(cp) && *taxonName)
	{
		strncpy(rank,rankName,BUFSZ) ;
		rank[BUFSZ] = '\0' ;
		strncpy(taxon,taxonName,BUFSZ) ;
		taxon[BUFSZ] = '\0' ;

		/* Check for a subtree here... */
		if( !((parseTreeType == TAXONOMIC_TREE &&
		    (parseList(cp,taxon)|| printAtom(taxon))) ||
		    parseList(cp,0)))
		  /* then */ blankAtom() ;

		skipSpaces(cp) ;

		/* Parse an branch "distance" if present,
		   for a blank atomic node, not a list...*/
		if( **cp == ':')
		{
			++(*cp) ;
			if( !parseDistance(cp) )
			{
				messout("Expecting distance value after colon following taxon label?") ;
				return FALSE ;
			}
		}

		parseBootstrap(cp) ;

		skipSpaces(cp) ;

		printTaxon(rank,taxon) ;

		if( **cp == ',' )
		{
			/* Substitute blank to ensure that subsequent empty names 
			    preceeding the next ',' or ')' are scanned */
			**cp = ' ' ; 
			return TRUE ;
		}
		if( **cp == ')') return TRUE ;
		messout("Missing comma or unmatched \")\" after parsing a taxon label?") ;
		return FALSE ;
	}

	/* Otherwise, parse for a full node */
	if(*(name = parseName(cp)))
	{
		strncpy(label,name,BUFSZ) ;
		label[BUFSZ] = '\0' ;
	}
	else label[0] = '\0' ; /* will force failure iff this null label hits printAtom()*/

	if(parseTaxon(cp))
	{
		strncpy(rank,rankName,BUFSZ) ;
		rank[BUFSZ] = '\0' ;
		if(*taxonName)
		{
			strncpy(taxon,taxonName,BUFSZ) ;
			taxon[BUFSZ] = '\0' ;
		}
		else if(parseTreeType == TAXONOMIC_TREE && *label)
		{
			strncpy(taxon,label,BUFSZ) ;
			taxon[BUFSZ] = '\0' ;
		}
		/* ignore the rank entirely, since we need a taxon name too? */
		else rank[0] = '\0' ; 
	}
	else rank[0] = '\0' ;

	if( parseList(cp,label) ||
	    printAtom(label) )
	{
		skipSpaces(cp) ;

		/* Parse an branch "distance" if present */
		if( **cp == ':')
		{
			++(*cp) ;
			if( !parseDistance(cp) )
			{
				messout("Expecting distance value after colon after parsing non-blank label?") ;
				return FALSE ;
			}
		}
		parseBootstrap(cp) ;

		skipSpaces(cp) ;

		printTaxon(rank,taxon) ;

		if( **cp == ',' )
		{
			/* Substitute blank to ensure that subsequent empty names 
			    preceeding the next ',' or ')' are scanned */
			**cp = ' ' ; 
			return TRUE ;
		}
		if( **cp == ')') return TRUE ;

		messout("Missing comma or unmatched \")\" while parsing a tree list member?") ;
		return FALSE ;
	}
	messout("Expecting non-NULL tree list label?") ;
	return FALSE ;
}

static int level = 0 ; 

static BOOL parseList(char **lpTreeList,char *name)
{
	char c, *label, listname[BUFSZ+1],
		 rank[BUFSZ+1], taxon[BUFSZ+1] ;
	LPTREENODE parent = currentParent ;

	if( name && *name )
	{
		strncpy(listname,name,BUFSZ) ;
		listname[BUFSZ] = '\0' ;
	}
	else
		listname[0] = '\0' ;

	skipSpaces(lpTreeList) ;

	if( **lpTreeList != '(')
	{
		return FALSE ;
	}

	/* Recursing into a subtree... */
	push(parentNodes,currentParent,LPTREENODE ) ;
	
	/* New subtree root node... */
	currentParent = messalloc(sizeof(TREENODE)) ;
	currentParent->id = ++nextId ;
	if(*listname)
		strcpy(currentParent->name,listname) ;
	else
		strncpy(currentParent->name,
		        messprintf("%s%d",treeObjName,currentParent->id),
		        BUFSZ) ;

	(*lpTreeList)++ ;
	skipSpaces(lpTreeList) ;
	while( (c = **lpTreeList) && (c != ')') && (c != ';') )
	{
		if( !parseMember(lpTreeList) ) 
			return FALSE ;
		/* assume the end of member in the .ace file? */
		fprintf(ttfil,"\n") ;  
	}

	if( **lpTreeList == ')') /* a list must terminate with a ')'? */
	{
		(*lpTreeList)++ ;  /* move ahead one character */

		/* Create a subtree interior node? */
		if(currentParent->name[0] && (parseTreeType == TAXONOMIC_TREE || parseTreeType == CELL_LINEAGE))
			fprintf(ttfil,"TreeNode : \"%s\"\n",currentParent->name) ;
		else if(currentParent->id>0)
			fprintf(ttfil,"TreeNode : \"%s%d\"\n",treeObjName,currentParent->id) ;
		else if(currentParent->id==0)
			fprintf(ttfil,"TreeNode : \"%s\"\n",treeObjName) ;

		/* initialized to empty, in case not found */
		rank[0] = taxon[0] = '\0' ;

		/* Parse for optional bootstrap label value if present? */
		if(*(label = parseNumber(lpTreeList)))
			fprintf(ttfil, "Bootstrap %s\n", label) ;
		/* ... or a postfix specified Label for the node? */
		else if(*(label = parseName(lpTreeList)))
		{
			/* Replace prefix name with postfix name, for node label? */
			strncpy(listname,label,BUFSZ) ;
			listname[BUFSZ] = '\0' ;

			fprintf(ttfil, "Label \"%s\"\n", listname) ;

			/* Check for a {taxon} here too? */
			if(parseTaxon(lpTreeList))
			{
				strncpy(rank,rankName,BUFSZ) ;
				rank[BUFSZ] = '\0' ;
				if(*taxonName)
				{
					strncpy(taxon,taxonName,BUFSZ) ;
					taxon[BUFSZ] = '\0' ;
				}
				else if(parseTreeType == TAXONOMIC_TREE && listname[0])
				{
					strncpy(taxon,listname,BUFSZ) ;
					taxon[BUFSZ] = '\0' ;
				}
			}
		}

		fprintf(ttfil,"Id %d\n",currentParent->id) ;
		if(currentParent->id == 0)
			fprintf(ttfil,"Root\n") ;
		else
			fprintf(ttfil,"Interior\n") ;
		fprintf(ttfil,"Tree \"%s\"\n",treeObjName) ;
		if(parent)
		{
			if(parent->name[0])
				fprintf(ttfil,"Parent \"%s\"\n",parent->name) ;
			else
				fprintf(ttfil,"Parent \"%s%d\"\n",treeObjName,parent->id) ;
		}

		printTaxon(rank,taxon) ;

		/* discard the old current parent */
		messfree(currentParent) ;

		/* the restore the previous parent*/
		currentParent = pop(parentNodes, LPTREENODE ) ;

		return TRUE ;
	}
	messout("Unmatched \")\" while parsing a tree structure?") ;
	return FALSE ;
}

/* Do it the "easy" way: create a temporary .ace file
   then read the temporary file into acedb; this avoids
   the tedium of my dealing with low level bs routines? */
static BOOL parseTree(char *TreeList)
{
	currentParent = 0 ;
	level = 0 ;
	leafNo = 0 ;
	nextId = 0 ;

	if( parseList(&TreeList,treeObjName)
	    && ( skipSpaces(&TreeList), *TreeList == ';'))
	{
		/* Store the tree? */
		return TRUE ;
	}
	messout("Improperly formed tree structure or missing semi-colon?") ;
	return FALSE ;
}

/* This routine reads in a text file containing a "New Hampshire" 
   formatted tree file, into a ?Tree and ?TreeNode ACEDB representation
   of the tree;  The function is self contained: the user is prompted 
   for the file name and the name of the ACEDB tree */
void readTreeFile(int treeType)
{
	FILE *fil = NULL ;
	char *textbuf = 0, *lpszTree = 0, *treeName = 0 ;
	int size;
	char *nameptr ;

	if(!*dirName)
	  {
	    if (getenv ("ACEDB_DATA"))
	      strcpy (dirName, getenv("ACEDB_DATA")) ;
	    else if (filName ("rawdata", 0, "r"))
	      strcpy (dirName, filName ("rawdata", 0, "r")) ;
	  }

	if( !(fil = filqueryopen (dirName, fileName, "ph", "r", "New Hampshire Tree File")))
		return ;

	/* determine filesize */
	fseek (fil, 0, SEEK_END);
	size = ftell (fil);
	rewind (fil);

	/* Read the whole file into memory? */
	lpszTree = textbuf = messalloc (size*sizeof (char) + 16);	
	fread (textbuf, sizeof (char), size, fil);
	textbuf[size] = '\0';		/* add string terminator */

	fclose (fil);

	/* Prompt the user for the name of the tree object? */
	treeObjName[0] = '\0' ;

	if(*(treeName = parseName(&lpszTree)))
		strncpy(treeObjName,(const char *)treeName,BUFSZ-1) ;
	else
	{		
		if( messPrompt("Please give a name for your ACEDB tree object:",treeObjName,"w"))
		{
			strncpy(treeObjName,freeword(),BUFSZ-1) ;
			treeObjName[BUFSZ-1] = '\0' ;
		}
		else
			strcpy(treeObjName,"Dendrogram") ;
	}

	/* Remember parseTreeType */
	parseTreeType = treeType ;

	/* Parse the tree object into ACEDB */
#if defined(WIN32) && defined(_DEBUG)
	ttfil = fopen("dgram.trace","w+") ;
#else
	ttfil = filtmpopen (&nameptr, "w+") ; /* temporary .ace file for tree data ?*/
#endif
	
	/* Write the tree object */
	fprintf(ttfil,"Tree : \"%s\"\n",treeObjName) ;
	switch(parseTreeType)
	{
		case TAXONOMIC_TREE:
			fprintf(ttfil,"Type Taxonomy\n") ;
			fprintf(ttfil,"Top\n") ;  
			break ;
		case DNA_TREE:
			fprintf(ttfil,"Type DNA\n") ;
			fprintf(ttfil,"Middle\n") ; 
			fprintf(ttfil,"Bootstrap_Factor 1.0\n") ; 
			break ;
		case PROTEIN_TREE:
			fprintf(ttfil,"Type Protein\n") ;
			fprintf(ttfil,"Middle\n") ; 
			fprintf(ttfil,"Bootstrap_Factor 1.0\n") ; 
			break ;
		case CELL_LINEAGE:
			fprintf(ttfil,"Type Cell_Lineage\n") ;
			fprintf(ttfil,"Middle\n") ;  
			break ;
		default:  /* Otherwise, unknown tree type -- do nothing? */
			break ;
	}
	fprintf(ttfil,"Root \"%s\"\n\n",treeObjName) ;

	if( !parentNodes ) 
		parentNodes = stackCreate(128) ;
	else 
		stackClear (parentNodes) ;

	if( parseTree(lpszTree) )
	{
		KEY cKey, tKey ;
		rewind(ttfil) ;
		parseOneFile (ttfil, 0) ;
		if (  lexword2key("Tree", &cKey, _VClass) &&
		      lexword2key (treeObjName, &tKey, cKey))
			display (tKey,0, 0) ;
	}
	else 
		messout("Tree file \"%s\" is not a valid New Hampshire tree?\n", fileName) ;

	fclose(ttfil) ;

#if !(defined(WIN32) && defined(_DEBUG))
	/* Free the temporary file */
	filtmpremove (nameptr) ;
#endif

	/* Free the text buffer */
	messfree(textbuf) ;
}

void readTaxonomicTree(void)
{
	readTreeFile(TAXONOMIC_TREE) ;
}

void readDNATree(void)
{
	readTreeFile(DNA_TREE) ;
}

void readProteinTree(void)
{
	readTreeFile(PROTEIN_TREE) ;
}

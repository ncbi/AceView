/*  $Id: acesybase.c,v 1.2 2005/06/20 01:01:11 mieg Exp $  */
/*  File: acesybase.c
 *  Author: Detlef Wolf
 *  Copyright (C) J Thierry-Mieg and Detlef Wolf, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: translated between ACEDB tree and acesyb representation
 * Exported functions:
       sybaseInit
       sybaseKill
       sybaseGet
       sybaseCommit
       sybaseStore
 * HISTORY:
 * Last edited: Sep  1 14:20 1994 (srk)

 * Created: Fri Oct 18 20:25:36 1991 (mieg)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include "acedb.h"
#include "acesybase.h"
#include "pick.h"
#include "lex_bl_.h"
#include "bs_.h"
#include "chrono.h"
#include "acesyb.h"
#include "systags.h"

/* #define CHRONO */

/* --------------------- module memory ---------------------------- */
static int currentNodeNo;  /* used in tree2array */
static KEY currentKey;     /* global for efficiency */
static baseNode;           /* index of the the first node of the
                              current object */
static KEY lastTag;        /* number of last tag converted */
static int valPos;         /* position after a tag */
static Array nodeBuf;      /* I/O buffer */
static int acesybaseTrace=0;  /* trace level, 0=none */

/* --------------------- local functions ---------------------------- */

/* --------------------------------------------------------------- */
static void sybasePrintRec(BS bs, int offset) {
  /* recursively descend into the tree, right before down */

  int i;
  int l=0;  /* number of chars printed */

  printf(" %d/%d(%s)", class(bs->key), KEYKEY(bs->key), bsText(bs)); l+=10+1;

  if (bs->right != NULL)
    sybasePrintRec(bs->right, offset+l);
  if (bs->down != NULL) {
    printf("\n");
    for (i=1;i<=offset;i++) printf(" ");
    sybasePrintRec(bs->down, offset);
  };
}


/* notes on bs_.h, that declares the cell data structure

bshow.key -- tells the type of cell:
  0 .. _LastN      : a data cell (Int, Float, Text); n.x holds data
                     I accept only _Int, _Float and _Text
  _LastN .. 2^24-1 : a tag; n is not relevant
  2^24 .. 2^32-1   : a pointer; n.key is redunant but I fill it anyway

my comments on the code:
typedef union  { int i;
                 unsigned int u;   -- not used
                 long L; -- not used
                 unsigned long U; -- not used
                 float f;
		 mytime_t time ;
                 TEXTKEY text;  -- not used
                 KEY key ;    -- not used, but I fill it 
               } BSdata;
 
typedef struct bshow *BS;
typedef struct btext *BT;
 
        struct bshow {BS  up, down, right;
                      KEY key;  -- type/tagNo/Pointer
                      BSdata n;  -- data depending on key
                      BT bt ; -- textdata (should be in BSdata, strange design)
                              -- bt is not allocated by BSalloc
                      int size;  
                     } ;
 
        struct btext { char *cp ; 
		       BS   bsm ;  only used during treedisp update 
                     } ;


*/
/* --------------------------------------------------------------- */
static void bs2node (BS bs, ayNodeType *np) {
  /* data conversion from ACEDB tree cell into acesyb node 
     text data: since it do not know, whether the caller
     will still need them, I allocate new memory and make a copy
     precondition: valPos is set
     postcondition: lastTag is set
                    valPos is set
  */
  if (bsIsTag(bs)) {   /* just a tag, translated to a boolean */
    np->dataType=AY_BOOL;
    lastTag= (np->data.boolData=KEYKEY(bs->key));
    valPos=0;   /* new tag starts */
  } else {
  valPos++;
  if (class(bs->key)) {   /* a pointer */
    np->dataType=AY_KEY;
    np->data.keyData=bs->key;
  } else {  /* must be a simple data type */
    switch (KEYKEY(bs->key)) {
      case _Int: np->dataType=AY_INT; np->data.intData=bs->n.i; break;
      case _Float: np->dataType=AY_FLOAT; np->data.floatData=bs->n.f;break;
      case _Text: np->dataType=AY_TEXT; 
                  np->data.textData=messalloc(strlen(bs->bt->cp));
                  strcpy(np->data.textData, bs->bt->cp);
                  break;
      case _DateType: np->dataType=AY_DATE; np->data.dateData=bs->n.time; break;
      default: messcrash("bs2node: unknown data type %d", KEYKEY(bs->key));
    } /* switch */
  }; /* simple data type */
  } /* not a tag */
};

/* --------------------------------------------------------------- */
static void node2bs (ayNodeType *np, BS bs) {
  /* data conversion from acesyb node into ACEDB tree cell 
     text data: since I assume that the current contents of the
     node buffer is discarded anyway, I give the text memory of the
     node as a present to the tree and zero the pointer in the node.
  */
  switch (np->dataType) {
/*    case AY_EMPTY: printf("AY_EMPTY\n"); break; */
    case AY_INT  : bs->key=_Int; bs->n.i=np->data.intData; break;
    case AY_FLOAT: bs->key=_Float; bs->n.f=np->data.floatData; break;
    case AY_TEXT : bs->key=_Text;
 	           bs->bt = BTalloc();
                   bs->bt->cp=np->data.textData;  
		   messcrash ("acesybase: bs->bt->cp should be allocated on obj->x->handle, this code shoudl be upgraded, sorry") ;
                   np->data.textData=NULL;  /* steal the memory */
                   break;
    case AY_BOOL : bs->key=np->data.boolData; break;
    case AY_KEY  : bs->key=np->data.keyData; 
                   bs->n.key=np->data.keyData;  /* necessary? */
                   break;
    case AY_DATE : bs->key=_DateType; bs->n.time=np->data.dateData; break;
    default : messcrash("node2bs: invalid data type %d",np->dataType);
  };
};

/* --------------------------------------------------------------- */
static void tree2array (BS bs) {
  /* convert an ACEDB tree into an acesyb array
     precondition: currentNodeNo contains the index of the node
                   to be generated;
                   currentKey contains the key of the current ACEDB object
     postcondition: (after last return)
                    nodeBuf is filled with contents of bs
     algo: the father creates its children, so he knows them and
           can point to them
  */
  ayNodeType *node=arrayp(nodeBuf, arrayMax(nodeBuf), ayNodeType);
  node->key=currentKey;
  node->thisNode=currentNodeNo;
  node->rightNode=NULL;
  node->downNode=NULL;
  bs2node(bs, node);        /* convert simple data */
  node->classNo=class(currentKey);
  node->tagNo=lastTag;  /* should get no of last tag converted by bs2node */
  node->valPos=valPos;  /* set by bs2node */
  if (bs->right) {
    node->rightNode=++currentNodeNo;
    tree2array(bs->right);
  };
  if (bs->down) {
    node->downNode=++currentNodeNo;
    tree2array(bs->down);
  };
}

/* --------------------------------------------------------------- */
static void array2tree (int nno, BS bs) {
  /* convert an acesyb array into an ACEDB tree 
     input: nno  -- index in nodeBuf of node to convert,
            bs   -- pointer to receiving cell
     precondition: nodeBuf is corrently wired
                   baseNode contains the index of the first node of the
                   current object
     postcondition: (after last return)
                    the tree is attached to bs
     algo: the father creates its children, so he knows them and
           can point to them
  */
  ayNodeType *np=arrp(nodeBuf, baseNode+nno, ayNodeType); 
  bs->right=NULL;
  bs->down=NULL;
  bs->key=np->key;
  node2bs(np, bs);        /* convert simple data */
  if (np->rightNode) {
    BS son=BSalloc();     /* all elements are filled with zeros */
    son->up=bs;           /* make the son point to his father */
    bs->right=son;
    array2tree(np->rightNode, son);
  };
  if (np->downNode) {
    BS son=BSalloc();
    son->up=bs;           /* make the son point upward to his father */
    bs->down=son;
    array2tree(np->downNode, son);
  };
}

/* --------------------- exported functions ---------------------------- */

/* --------------------------------------------------------------- */
void sybaseInit (void)
{
   if (getenv("ACESYBASE_TRACE")) {
     sscanf(getenv("ACESYBASE_TRACE"),"%d",&acesybaseTrace);
     printf("sybaseInit: trace level set to %d\n", acesybaseTrace);
   };
   nodeBuf=arrayReCreate(nodeBuf, 128, ayNodeType);
   baseNode=0;    /* assume: only one object at a time */
   ayBegin();
}

/* --------------------------------------------------------------- */
void sybaseKill (KEY key)
{ 
  LEXP q;
  ayNodeType *np;  

  chrono("SybaseKill") ;
  
  q=KEY2LEX(key);
  if(q->diskaddr != 0)
    { q->diskaddr = 0 ;           /* reset faked disk address */
      lexmark(class(key));
    }

  if (acesybaseTrace)
    printf("sybaseKill: key=%d/%d", class(key), KEYKEY(key));

  ayFreeNodeBuf(nodeBuf);
  np=arrayp(nodeBuf, 0, ayNodeType);
  np->key=key;
  np->thisNode=0;
  if (ayDelete(nodeBuf))
    messcrash("sybaseKill: could not find key=%d/%d", class(key), KEYKEY(key));

  chronoReturn () ;
  return ;
}


/* --------------------------------------------------------------- */
BS sybaseGet (KEY key) { 
  ayQueryType query;
  BS bs = BSalloc();
  if (acesybaseTrace)
    printf("sybaseGet: key=%d/%d\n", class(key), KEYKEY(key));

  query.mask.key=key;
  ayRead(nodeBuf, &query);
  if (!arrayMax(nodeBuf))
    messcrash("sybaseGet: key=%d not found", key);
  bs->up = NULL;        /* mark this to be the root */
  array2tree (0, bs);
  if (acesybaseTrace>1) {
    printf("bsTreeGet: ");
    sybasePrint(bs);
  };

  if (bs->key != key) 
    messcrash("sybaseGet: requested key %d, but got %d", key, bs->key);
  return bs ;
}

/* --------------------------------------------------------------- */
void sybaseCommit() {
  ayEnd(1); 
}

/* --------------------------------------------------------------- */
void sybaseStore (BS  bs)
{
  KEY key = bs->key ;
  LEXP q;
  chrono("bsTreeSybaseStore") ;
  
  if (acesybaseTrace)
    printf("sybaseStore: key=%d/%d, %d\n", class(key), KEYKEY(key), key);
  if (acesybaseTrace>1) {
    printf("sybaseStore:");
    sybasePrint(bs);
  };
  

  q=KEY2LEX(key);
  if(q->diskaddr != 1)    /* fake a disk address */
    { q->diskaddr = 1;
      lexmark(class(key));
    }
  currentNodeNo=0;         /* initialize tree2array */
  lastTag=0;
  valPos=0;
  currentKey=key;          /* store here for efficient passing to tree2array */
  ayFreeNodeBuf(nodeBuf);  /* initialize receiving array */
  tree2array(bs);          /* result is in nodeBuf */
  ayWrite(nodeBuf);

  chronoReturn(); 
  return ;
}



/* --------------------------------------------------------------- */
void sybasePrint(BS bs) {
  printf("--------tree %d/%d\n", class(bs->key), KEYKEY(bs->key));
  if (bs->down) printf("PROBLEM: root cell has down pointer\n");
  if (bs->right) sybasePrintRec(bs->right,1);
  printf("\n");
}


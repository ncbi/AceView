/*  File: aqlcheck.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
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
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: aqlcheck.c,v 1.13 2015/01/31 01:27:38 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  2 15:51 2004 (edgrif)
 * * Sep  7 10:31 1999 (fw): fixed propagation of where-clause when preDecls
 *              expand the tree into list of FROM_LOCs (SANgc05275)
 * * Mar 19 14:52 1999 (fw): recursive function to count number of columns in a table
 *              to enable checking of nested and consecutive table defs/ops
 * * Aug 25 16:26 1998 (fw): nodes containing table don't carry the column-types 
				in the node->name anymore
 * * Aug 14 16:18 1998 (fw): checkField() - find column names by their last
                             tag names, "Email" in a->Address[Email] or 
			     "Full_name" in a->cast->Full_name
 * * Aug 14 14:58 1998 (fw): case insensitive matching of var/tag-names
                 'select m->Cast, m from m in class movie where exists m->cast'
 * * Aug 14 12:38 1998 (fw): finally fixed 'select all @t' --> 'select :1 from @t'
 * * Aug 11 15:22 1998 (fw): preTableProcess - accept 'select all @tvar'
 * * Aug  5 17:34 1998 (fw): introduced threadsafe AQL handle
 * * Aug  3 11:59 1998 (fw): table-values created for TABLE_SFW and TABLE_OP nodes
                             now use the resultHandle store.
 * * Jul 27 15:48 1998 (fw): checked error-codes/messages
 * * Jul 24 11:17 1998 (fw): tableFunc-field syntax for '$x := count @t:2'
 * * Jul 22 17:37 1998 (fw): check argument of tableFunc to return just one column
 * * Jul 22 16:06 1998 (fw): changed postTableProcess to use a switch-construct
 * * Jul 20 15:39 1998 (fw): catch for illegal outside declaration error 959, 
                                caused by i.e. '$x := yeardiff(p->born, now)'
 * * Jul 16 17:17 1998 (fw): AC_HANDLE made global and renamed globalStoreHandle
 * * Jul 16 10:47 1998 (fw): aqlPostTableCheck() doesn't have table-values yet,
                             so we can inly check columnNumber compatabilty for TABLE_OPs
 * * Jul 15 13:47 1998 (fw): fixed checkField, which doesn't need a table value, 
                  but just the TABLE_SFW node to check sort and column fields
 * * Jul 15 10:10 1998 (fw): fixed pre-processing of LOCAL_FIELD to be similar
                             to FOLLOW_TAG and LOCAL_POS
 * * Jul 14 14:45 1998 (fw): removed markChildren(), as row variables
                             are now recognised by through a special syntax.
 * * Jul 10 10:18 1998 (fw): order_by sort field criterion is NOT case-sensitive.
 * * Jul  8 18:21 1998 (fw): functionality kind of complete now.
 * Created: Fri Oct 25 16:04:41 1996 (rd)
 *-------------------------------------------------------------------
 */

#include <wh/acedb.h>
#include <wh/aceio.h>
#include <waql/aql_.h>

/* #define DEBUG 1 more intermediate parsetree output */

/******************************************************************/

static BOOL  preProcessLocPos (AqlNode *inNode, AQL aql, int level);

static BOOL  preLocCheck (AqlNode *inNode, AQL aql, int level);

static BOOL  preDecls (AqlNode *inNode, AQL aql, int level);
static void  postDecls (AqlNode *inNode, AQL aql, int level);

static BOOL  preScope (AqlNode *inNode, AQL aql, int level);
static void  postScope (AqlNode *inNode, AQL aql, int level);

static BOOL  preTableProcess (AqlNode *inNode, AQL aql, int level);
static void  postTableProcess (AqlNode *inNode, AQL aql, int level);

static int   numberOfTableColumns (AqlNode *inNode, AQL aql);

static AqlResultDisposition aqlNeedsSorting(AqlNode *inNode, AqlNode **psortNode);
static void aqlSetParents(AqlNode *inNode);

static void checkEnvFinalise (void *block);
/******************************************************************/


/***************************************************************/
/******* check query structures, prepare for evaluation ********/

/* strategy:

+  expand implicit declarations in locs etc.
   
+    build Dict of variable names
+      establish whether object or row variables
       assign possible types, doing context checking

   evaluate constant numerical expressions where possible

+  for table_expr:
    build table and assign types where possible
+      check compatibility of columns in table operations

   check/assign node types as far as possible
*/

/***********************************************************************/

AqlQuery *aqlQueryCreate (const char *queryText, AC_HANDLE objectHandle)
/*
   (1) Creates a query object on the handle given.
   (2) A store-handle is created for this object - all allocation for
       the duration of the query processing is done on that handle,
       which can be destroyed after processing (remember to rescue the result table).
   (3) the passed querytext is copied to this object using the internal handle.
   (4) All unneeded whitespaces are stripped from the query->text.
    mieg 2015_01_25: the code was badly written
      it copied the str as many times as there were double spaces
      and used  strcpy (cws, cp) with overlapping strings, whci rightly failed on some harwares
*/
{
  AqlQuery *q ;
  char *cp, *cq ;

  q = (AqlQuery*)halloc(sizeof(AqlQuery), objectHandle);

  /* all memory allocated during the processing of the query will
     be put on this store handle, and destroyed at the end of processing */
  q->handle = handleCreate () ;

  /* make a copy of the query text for this query object */
  q->text = strnew ((char*)queryText, q->handle); 


  /* now strip out all unwanted characters from the query text */

  cp = cq = q->text;
  
  /* jump initial spaces */
  while (*cp == ' ' || *cp == '\t' || *cp == '\n' || *cp == '\r' || *cp == '\f')
    cp++ ; 
  /* copy all characters but not repeated spaces */
  while (*cp)
    {
      if (*cp == ' ' || *cp == '\t' || *cp == '\n' || *cp == '\r' || *cp == '\f')
	{
	  while (*cp == ' ' || *cp == '\t' || *cp == '\n' || *cp == '\r' || *cp == '\f')
	    cp++ ;
	  if (*cp != ',' && *cp != ')' && *cp != ']' ) *cq++ = ' ' ;
	}
      if (*cp) *cq++ = *cp++ ;
    }
  /* zero terminate */
  *cq = 0 ;

  return q ;
} /* aqlQueryCreate */

/**************************************************************/
/**** scope
  at each SFW make a new scope
  build dicts from bottom up (post-order traversal), 
  	transferring up anything not defined
  in the same tree traversal: 
    convert SFW_ALL to standard SFW
    expand complex locators at the same time
****/

AqlScope* aqlScopeCreate (AC_HANDLE handle)	/* also used in aql.c */
{
  AqlScope *s = (AqlScope*)halloc (sizeof(AqlScope), handle) ;
  s->dict = dictHandleCreate (32, handle) ;
  /*
    no longer needed
    dictAdd (s->dict, "__JUNK", 0) ; fill first entry as index 0 is used to signify something else(?)
  */
  s->loc = arrayHandleCreate (32, AqlLoc*, handle) ;
  return s ;
} /* aqlScopeCreate */

/***********************************************************************/

/* First stage query checking - scopes of variables etc. */
void aqlCheck1 (AQL aql)
{
  if (!aql)
    aqlError (aql, 700, 0, "aqlCheck1() - called with NULL aql pointer") ;

  if (!aql->query)
    aqlError (aql, 701, 0, "aqlCheck1() - called with NULL aql->query") ;

  /* no point continuing if we already had a parse-error */
  if (aql->IsError)
    aqlError (aql, 702, 0, "aqlCheck1() - cannot continue because of earlier error.") ;
  

  /******************** query pre-processing *******************/

  aql->query->root->up = NULL;
  aqlSetParents(aql->query->root);

  /***************************************/
  /* map the TABLE_SFW_ALL to TABLE_SFW,
   * process the loc[n] statements to loc[1][1]..n-times */
  aqlTraverse (aql->query->root, aql, preProcessLocPos, 0, 0) ;

#ifdef DEBUG
  /* extra debug output not really needed */
  if (aql->parseDebugLevel >= 2)
    {
      aceOutf (aql->debug_out, "\n\n"); 
      aqlNodeOut (aql->debug_out, aql->query->root, 0, "processed after preProcessLocPos", FALSE); 
      aceOutf (aql->debug_out, "\n\n"); 
    }
#endif /* DEBUG */

  /* create a self-contained object, which holds data for the pre-processing
   * it can be destroyed when check1 is finished */
  aql->checkEnv = (AqlCheckEnv*)halloc (sizeof(AqlCheckEnv), aql->query->handle);
  blockSetFinalise (aql->checkEnv, checkEnvFinalise);

  aql->checkEnv->declCount = 0 ; /* used to create unique names (integers 1..n) for new declarations */
  aql->checkEnv->pCurrDecl = 0 ; 
  aql->checkEnv->pCurrWhere = 0 ; 
  aql->checkEnv->currFromLoc = 0 ;
  aql->checkEnv->declStack = stackCreate (64);
  aql->checkEnv->whereStack = stackCreate (64);
  aql->checkEnv->fromStack = stackCreate (64);

  /***************************************/
  /* set up declarations, create tmp-variables for FROM_LOCs etc.. */
  /* make implicit decls explicit */
  aqlTraverse (aql->query->root, aql, preDecls, postDecls, 0) ;

#ifdef DEBUG
  /* extra debug output not really needed */
  if (aql->parseDebugLevel >= 2)
    {
      aceOutf (aql->debug_out, "\n\n"); 
      aqlNodeOut (aql->debug_out, aql->query->root, 0, "processed after decls", FALSE);
      aceOutf (aql->debug_out, "\n\n");
    }
#endif /* DEBUG */

  /**************************************/
  /* build scope & assign Loc structures */
  /* the scope for locators is on a per table_expr basis */
  aql->query->locScope = aqlScopeCreate (aql->query->handle);
  aql->checkEnv->currentLocScope = aql->query->locScope;

  /**************************************/
  /* the scope for table variables and scalar variables is on a per whole-query basis */
  /* we make a current scope and the global scope for the whole query becomes its parent */

  aql->checkEnv->currentTableScope = aql->query->tableScope;

  aql->checkEnv->currentScalarScope = aqlScopeCreate(aql->query->handle) ;
  aql->checkEnv->currentScalarScope->parent = aql->query->scalarScope; /* this is the environment scope of the query or NULL */

  /*******************************************************************/
  /* set up scopes by assigning variable-names to locators,          *
   * find locator definitions (the nodes that define a locator) etc. */
  aqlTraverse (aql->query->root, aql, preScope, postScope, 0) ;

#ifdef DEBUG
  /* extra debug output not really needed */
  if (aql->parseDebugLevel >= 2)
    {
      aceOutf (aql->debug_out, "\n\n");
      aqlNodeOut (aql->debug_out, aql->query->root, 0, "processed after scoping", FALSE);
      aceOutf (aql->debug_out, "\n\n");
    }
#endif /* DEBUG */

  /************************************************************************/
  /* create variable copies, which duplicates locator definitions         *
   * this'll make sure every from_loc->left LOC_xxx-node has a ->nxt node */
  aqlTraverse (aql->query->root, aql, preLocCheck, 0, 0) ;

  /* re-scope the parsetree, which is necessary in case the
   * previous command copied branches of the tree, whose new nodes
   * won't have correct locators set up */
  aqlTraverse (aql->query->root, aql, preScope, postScope, 0) ;
  
  /***************************************/
  /* process tables, init table variables etc. */
  aqlTraverse (aql->query->root, aql, preTableProcess, postTableProcess, 0) ;

  messfree (aql->checkEnv);

  return ;
} /* aqlCheck1 */




/************************************************************************/
/*********************** private functions ******************************/
/************************************************************************/

static void checkEnvFinalise (void *block)
{
  AqlCheckEnv *checkEnv = (AqlCheckEnv*)block;

  stackDestroy (checkEnv->declStack);
  stackDestroy (checkEnv->whereStack);
  stackDestroy (checkEnv->fromStack);

  return;
} /* checkEnvFinalise */

/******************** a little utility to make nodes ******************/

static AqlNode *makeNewNode (AqlNodeType type, AC_HANDLE handle)
{ 
  AqlNode* newNode = (AqlNode*)halloc (sizeof(AqlNode), handle) ;

  newNode->type = type ;

  /* init the value-type character code,
   * the character '0' means not yet initialised */
  newNode->vtype = '0';

  return newNode ;
} /* makeNewNode */



static BOOL preProcessLocPos (AqlNode *inNode, AQL aql, int level)	
/* called by aqlTraverse as a preFunc() as part of aqlCheck1() 
   so it is executed when going DOWN the tree in the recursive traversal */
/**************/
/* NOTE by RD : must be separated from preDecls() */
/**************/
{
  /* Map TABLE_SFW_ALL to TABLE_SFW, before starting surgery on the tree
   * This is done by creating a SELECT_FIELD dangling off the right leaf
   *  of the TABLE_SFW_ALL-inNode and then changing its type to TABLE_SFW  */
  if (inNode->type == nTABLE_SFW_ALL)
    { 
      AqlNode *fwFieldNode = inNode->left;
      AqlNode **pRightNode = &(inNode->right) ;

      /* SFW_ALL nodes can have only one FROM_xxx node - the "_DEF" locator */
      /* we transform 'select all class movie' into 'select _DEF from _DEF in class movie'
       * by creating a SELECT_FIELD that selects _DEF */

      if (fwFieldNode->type == nFROM_LOC)
	{ 
	  /* attach a SELECT_FIELD-node to the TABLE_SFW inNode
	   * and create a generic LOC_VAR node to be its left-node */
	  AqlNode *newSelectFieldNode = makeNewNode(nSELECT_FIELD, aql->query->handle) ;
	  AqlNode *newLocVarNode = makeNewNode(nLOC_VAR, aql->query->handle) ; 
	  
	  newSelectFieldNode->left = newLocVarNode ; 
	  newLocVarNode->name = fwFieldNode->name ; 
	  *pRightNode = newSelectFieldNode ; 
	  pRightNode  = &newSelectFieldNode->nxt ;

	}
      else if (fwFieldNode->type == nFROM_TABLE)
	{
	  /* Ideally , the expression 'select all @table' should be expanded to
	   * 'select a:1, a:2, a:3 from @table' for a table with 3 columns.
	   * However, at this point in the processing, we haven't yet mapped
	   * names to locs and can therefor not count the number of select-fields
	   * of the table by the name of '@table'.
	   * It becomes horribly complicated to wait until the name-scoping is
	   * done and append extra row-var select-fields in preTableProcess.
	   * We therefor expand 'select all @table' to 'select :1 from @table'
	   * regardless of the number of columns in @table (at least 1 is guaranteed).
	   */

	  /* attach a SELECT_FIELD-node to the TABLE_SFW inNode
	     and create a generic LOC_VAR node to be its left-node */
	  AqlNode *newSelectFieldNode = makeNewNode(nSELECT_FIELD, aql->query->handle) ;
	  AqlNode *newLocVarNode = makeNewNode(nLOC_VAR, aql->query->handle) ; 
	  AqlNode *newLocFieldNode = makeNewNode(nLOC_LOCAL_FIELD, aql->query->handle) ; 
	  
	  newSelectFieldNode->left = newLocVarNode ; 
	  newLocVarNode->name = fwFieldNode->name ; 
	  newLocVarNode->nxt = newLocFieldNode;
	  newLocFieldNode->number = 1;
	  *pRightNode = newSelectFieldNode ; 
	  pRightNode  = &newSelectFieldNode->nxt ;
	}

      /* this node is now like a normal select-from-where table expression */
      inNode->type = nTABLE_SFW ;
    }

  /* change all EXISTS_TAG nodes into EXISTS and 
   * deal with defaults for ->tag nodes in this context */
  if (inNode->type == nBOOL_EXISTS_TAG && inNode->left->nxt)
    {
      AqlNode *checkExistNode;

      /* go to the last node in the linked list. this is the one tested for non-NULL
       * e.g. "exists_tag a->address[email]" will become "exists a->address[email][0]"
       * or   "exists_tag a->email" will become "exists a->email[0]" */
      for (checkExistNode = inNode->left->nxt; checkExistNode->nxt; checkExistNode = checkExistNode->nxt);

      if (checkExistNode->type == nLOC_LOCAL_POS && checkExistNode->number == 0) /* we allow exists_tag <tag>[0] */
	{
	  for (checkExistNode = inNode->left->nxt; checkExistNode->nxt && checkExistNode->nxt->nxt; checkExistNode = checkExistNode->nxt);
	}

      if (checkExistNode->type == nLOC_FOLLOW_TAG_NAME || 
	  checkExistNode->type == nLOC_LOCAL_TAG_NAME)
	{
	  inNode->type = nBOOL_EXISTS;

	  /* treat free-ending tag nodes as tag[0] in an EXISTS_TAG context */
	  /* explicitly adding the [0] mode here will cause the next if-block
	   * to ignore this tag and not default it to tag[1] (when the recursion gets there)
	   * but later on the [0] node will be removed again (as no further eval is needed) */
	  if (!checkExistNode->nxt || checkExistNode->nxt->type != nLOC_LOCAL_POS)
	    {
	      AqlNode *newNode = makeNewNode (nLOC_LOCAL_POS, aql->query->handle) ;
	      
	      newNode->nxt = checkExistNode->nxt ;
	      checkExistNode->nxt = newNode ;
	      newNode->number = 0 ; 
	    }
	}
      else
	aqlError (aql, 903, inNode->pos, "Illegal argument to boolean clause `exists_tag' (tags-by-name only)");
    }

  /* add [1] to free-ending tag nodes */
  /* this means that obj->tag is interpreted as obj->tag[1]
   *  i.e. it will return the data that is one to the right of the tag,
   *  whereas obj->tag[0] is left untouched and tag[0] returns the tag itself */
  if (inNode->type == nLOC_FOLLOW_TAG_NAME || inNode->type == nLOC_LOCAL_TAG_NAME)
    { 
      if (!inNode->nxt || inNode->nxt->type != nLOC_LOCAL_POS)
        { 
	  AqlNode *newNode = makeNewNode (nLOC_LOCAL_POS, aql->query->handle) ;

	  newNode->nxt = inNode->nxt ;
	  inNode->nxt = newNode ;
	  newNode->number = 1 ;
	}
    }

  /* eliminate [0] and ->0 nodes, 
   * they will just return the actual object itself and 
   * they do not need this node as no evaluation is needed */
  if (inNode->nxt && 
      (inNode->nxt->type == nLOC_LOCAL_POS || inNode->nxt->type == nLOC_FOLLOW_POS) && 
      inNode->nxt->number == 0)
    { 
      AqlNode *skipNode = inNode->nxt ;
      inNode->nxt = inNode->nxt->nxt ;
      messfree (skipNode) ;
    }

  /* convert [k], k > 1 to k [1] copies */
  /* obj->tag[3] becomes obj->tag[1][1][1] */
  if (inNode->type == nLOC_LOCAL_POS)
    while (inNode->number > 1)
      {
	AqlNode *newNode = makeNewNode (nLOC_LOCAL_POS, aql->query->handle) ;

        newNode->nxt = inNode->nxt ; 
	inNode->nxt = newNode ;
	newNode->number = 1 ;
	--inNode->number ;
      }

  /* convert ->k (FOLLOW_POS for k > 1) into ->1 and (k-1)-many ->1 copies */
  /* obj->3 becomes obj->1[1][1] */
  if (inNode->type == nLOC_FOLLOW_POS) /* added (fw-980827) */
    while (inNode->number > 1)
      {
	AqlNode *newNode = makeNewNode (nLOC_LOCAL_POS, aql->query->handle) ;

        newNode->nxt = inNode->nxt ; 
	inNode->nxt = newNode ;
	newNode->number = 1 ;
	--inNode->number ;
      }


  return FALSE ;		/* FALSE tells aqlTraverse to continue in the tree */
} /* preProcessLocPos */

/******************************************************************/

/* a little utility to create a new FROM_xxx declaration
 * we need this so every declaration only has one from-clause
 * The conversion is as follows :
 *      "p->author->email" 
 * will be referred to as "tmp_4" with these new FROM_clauses
 *      from tmp_1 in p->Author[0], 
 *           tmp_2 in tmp_1[1], 
 *           tmp_3 in tmp_2->Email[0], 
 *           tmp_4 in tmp_3[1]
 * the new artificial FROM_clauses are created and hooked into the tree
 * by this function.
 */
static AqlNode *makeNewHookedDecl (AqlNodeType type,
                                   AqlNode *currNode,
                                   char *label,
                                   AQL aql)
{
  AqlNode *newNode;

  if (!(type == nFROM_TABLE || type == nFROM_LOC))
    aqlError (aql, 600, currNode->pos,
	      "makeNewHookedDecl() - type has to be FROM_LOC or FROM_TABLE") ;

  if (aql->parseDebugLevel >= 3)
    {
      aceOutf(aql->debug_out, "make artificial FROM_clause declaration, Label %s No %d - ",
                  label,
                  aql->checkEnv->declCount+1);
      aqlNodeOut (aql->debug_out, currNode, 10, "currNode", FALSE);
      aceOutf (aql->debug_out, "\n");
    }

  newNode = makeNewNode (type, aql->query->handle);

  /* insert the new node into the decalarations list */
  newNode->nxt      = *(aql->checkEnv->pCurrDecl);
  *(aql->checkEnv->pCurrDecl) =  newNode ; 
  aql->checkEnv->pCurrDecl    = &newNode->nxt ; /* rehook */

  if (aql->checkEnv->pCurrWhere)
    /* propagate the where-clause and attach it onto the new declaration */
    { 
      newNode->right     = *(aql->checkEnv->pCurrWhere);
      *(aql->checkEnv->pCurrWhere) = 0;
      aql->checkEnv->pCurrWhere    = &newNode->right;
    }

  /* the variable names of the artifical variables start with '_',
   * This would not be allowed by the lexer, so these variables can never
   * be mistaken for user-defined variables names 
   */
  newNode->name = strnew (messprintf ("_%s_%d", label, ++aql->checkEnv->declCount), aql->query->handle) ;

  return newNode ;
} /* makeNewHookedecl */




static BOOL preDecls (AqlNode   *inNode, 
		              AQL        aql,
                      int        level)
{
  static AqlNode *parent ;	/* nasty STATIC - oh my god */
  AqlNode *newNode, *newDeclaredNode ;

  if (aql->parseDebugLevel >= 4) {
      aqlNodeOut(aql->debug_out, inNode, level, "preDecls entered", FALSE);
      aqlCheckEnvOut(aql->debug_out, "preDecls", aql->checkEnv, level+1);
  }

/* Now start on remapping the declarations.  First, maintain our pointers. */

  if (inNode->type == nTABLE_SFW)		
    { 
      /* move into a new scope */
      push (aql->checkEnv->declStack, aql->checkEnv->pCurrDecl, AqlNode**) ;
      push (aql->checkEnv->whereStack, aql->checkEnv->pCurrWhere, AqlNode**) ;
      push (aql->checkEnv->fromStack, aql->checkEnv->currFromLoc, AqlNode*) ;

      aql->checkEnv->pCurrDecl = &inNode->left ; /* TABLE_SFW->left is a FROM_LOC or FROM_TABLE */
      aql->checkEnv->pCurrWhere = 0 ;
      aql->checkEnv->currFromLoc = 0 ;
#ifdef DEBUG
      aceOutf(aql->debug_out, "Setting currdecl\n");
      aqlNodeOut(aql->debug_out, inNode->left, level, "new currdecl", FALSE);
#endif
    }

  
  /* recursion is slightly awkward when going down the decl list itself
     we need to put new implicit decls before self when 
     processing ->left (the decl itself), after on ->right (the where
     clause) and ->nxt
  */

  if (aql->checkEnv->currFromLoc && 
      (inNode == aql->checkEnv->currFromLoc->right
       || inNode == aql->checkEnv->currFromLoc->nxt))
    { 
      aql->checkEnv->pCurrDecl = &(aql->checkEnv->currFromLoc->nxt) ;
      aql->checkEnv->pCurrWhere = &(aql->checkEnv->currFromLoc->right) ;
      aql->checkEnv->currFromLoc = 0 ;
#ifdef DEBUG
      aceOutf(aql->debug_out, "Changing currdecl\n");
      if (aql->checkEnv->currFromLoc == NULL)
      {
                aceOutf(aql->debug_out, "but nothing there\n");
      } else
      {
          aqlNodeOut(aql->debug_out, aql->checkEnv->currFromLoc->nxt, level, "change of currdecl", FALSE);
      }
#endif
    }

  if (inNode->type == nFROM_LOC || 
      inNode->type == nFROM_TABLE)
    {
      aql->checkEnv->currFromLoc = inNode ;
      /* we're switching to a new FROM_LOC, so we have to prevent
       * the where-clause to be propagated into that part 
       * of the parsetree (see SANgc05275) */
      aql->checkEnv->pCurrWhere = 0;
    }

/* now do the work itself */

  /* transfer TABLE based LOCs to straight row declarations */
  /* construct LOC_LOCAL_FIELDs from LOC_TABLE_FIELDS,
     these can then be processed in postTableProcess and go through
     the checkField process */
  if (inNode->type == nLOC_TABLE_FIELD_NAME || 
      inNode->type == nLOC_TABLE_FIELD)
    {				/* inNode will become a LOC_VAR */
      if (!aql->checkEnv->pCurrDecl)
	{
	  /* things would go horribly wrong if it is NULL (in makeNewHookedDecl) */
	  /* because, if pCurrDecl hasn't been initialised (by the occurence
	     of a TABLE_SFW-node) - it means we can't have a proper declaration
	     here. caused by e.g. '$x:= @t:2', which is a currently illegal query anyway */
	  
	  aqlError (aql, 958, inNode->pos, "invalid declaration");
	}
      
      newDeclaredNode = makeNewHookedDecl (nFROM_TABLE, inNode,
                                           inNode->type == nLOC_TABLE_FIELD_NAME ?
                                           "local_field_name" : "local_table_field",
                                           aql) ;
      newNode = makeNewNode (0, aql->query->handle) ;		/* assign type later */

      newDeclaredNode->left = inNode->left ; 
      inNode->left = 0 ;	/* transfer the table */
      
      if (inNode->type == nLOC_TABLE_FIELD_NAME)
	{
	  newNode->type = nLOC_LOCAL_FIELD_NAME ;
          newNode->name = inNode->name ;
	}
      else /* nLOC_TABLE_FIELD */
	{ 
	  newNode->type = nLOC_LOCAL_FIELD ;
	  newNode->number = inNode->number ;
	}

      /* insert newNode in linked-list behind newNode */
      newNode->nxt = inNode->nxt ; 
      inNode->nxt = newNode ;

      inNode->type = nLOC_VAR ; 
      inNode->name = newDeclaredNode->name ;
    }


  /* also transfer OBJECT() to a solo declaration */
  if (inNode->type == nOBJECT)
    if (!aql->checkEnv->currFromLoc || aql->checkEnv->currFromLoc->left != inNode || inNode->nxt)
      {					/* inNode will become a LOC_VAR */
        aqlTraverse (inNode->left, aql, preDecls, postDecls, level+1) ; /* first do left and right */
        aqlTraverse (inNode->right, aql, preDecls, postDecls, level+1) ;

	newDeclaredNode = makeNewHookedDecl (nFROM_LOC, inNode, "object", aql);
	newNode = makeNewNode (nOBJECT, aql->query->handle);
	newDeclaredNode->left = newNode ;

	newNode->left = inNode->left ; 
	newNode->right = inNode->right ; 

	inNode->left = 0 ; 
	inNode->right = 0 ;
	inNode->type = nLOC_VAR ; 
	inNode->name = newDeclaredNode->name ;
      }


  /* the tags, relpos etc. are trickiest - we have to leave one per decl */
  if (inNode->type == nLOC_LOCAL_POS   || inNode->type == nLOC_LOCAL_TAG_NAME  ||
      inNode->type == nLOC_FOLLOW_POS  || inNode->type == nLOC_FOLLOW_TAG_NAME ||
      inNode->type == nLOC_LOCAL_FIELD || inNode->type == nLOC_LOCAL_FIELD_NAME || /* this line added (fw-980715) */
      inNode->type == nLOC_METHOD) /* added (fw-980817) */
    {
      if (aql->checkEnv->currFromLoc && aql->checkEnv->currFromLoc->left->nxt == parent)
				/* inside declaration - handle differently */
	{				/* transfer parent to newNode->nxt */
      char *label = (aql->parseDebugLevel == 0) ?
          "transferparent" :
          (hprintf(aql->query->handle, "transfer_%s", aqlNodeTypeName(inNode->type)));

      if (aql->parseDebugLevel >= 2)
      {
          aceOutf(aql->debug_out, "Making transferparent for:\n");
          aqlNodeOutPretty(aql, aql->debug_out, inNode, FALSE, 8);
#if 1
          aceOutf(aql->debug_out, "\nwith parent:\n");
          aqlNodeOut(aql->debug_out, parent, 4, "parent", FALSE);
#endif
          aceOutf(aql->debug_out, "\n");
      }


	  newDeclaredNode = makeNewHookedDecl (nFROM_LOC, inNode, label, aql);
	  newNode = makeNewNode (nLOC_VAR, aql->query->handle) ;

	  newDeclaredNode->left   = newNode ;
	  newNode->name           = aql->checkEnv->currFromLoc->left->name ; 
	  newNode->nxt            = parent ;
	  parent->nxt             = 0 ;
	  aql->checkEnv->currFromLoc->left->name = newDeclaredNode->name ;
	  aql->checkEnv->currFromLoc->left->nxt  = inNode ;
	}
      else if (!aql->checkEnv->currFromLoc || aql->checkEnv->currFromLoc->left != parent)
				/* outside declaration */
	{				/*  transfer inNode to newNode->nxt */
	  if (!aql->checkEnv->pCurrDecl)
	    {
	      /* things would go horribly wrong if it is NULL (in makeNewHookedDecl) */
	      /* because, if pCurrDecl hasn't been initialised (by the occurence
		 of a TABLE_SFW-node) - it means we can't have a proper declaration
		 here. caused by e.g. '$x := yeardiff(p->born, now)' */
	      aqlError (aql, 959, inNode->pos, "illegal outside declaration");
	    }
	  
	  newDeclaredNode = makeNewHookedDecl (nFROM_LOC, inNode, "transferdecl", aql);
	  newNode = makeNewNode (nLOC_VAR, aql->query->handle) ;
	  
	  newDeclaredNode->left = newNode ;
	  newNode->name  = parent->name ; 
	  newNode->nxt   = inNode ;
	  parent->name   = newDeclaredNode->name ; 
	  parent->nxt    = inNode->nxt ; 
	  inNode = parent ;		/* so below, parent is set correctly! */
	}
    }
  parent = inNode ;	/* OK since loc expressions pass linearly through ->nxt pointers */
  if (aql->parseDebugLevel >= 4) {
      aqlNodeOut(aql->debug_out, inNode, level, "preDecls exited", FALSE);
  }
  return FALSE ;
} /* preDecls */

/************************************************************/

static void postDecls (AqlNode   *inNode,
		               AQL        aql,
                       int        level)
/* executed when coming UP the tree in the recursive traversal */
{
  if (inNode->type == nTABLE_SFW)
    { 
      /* getting back the current scope values */
      aql->checkEnv->pCurrDecl = pop (aql->checkEnv->declStack, AqlNode**) ;
      aql->checkEnv->pCurrWhere = pop (aql->checkEnv->whereStack, AqlNode**) ;
      aql->checkEnv->currFromLoc = pop (aql->checkEnv->fromStack, AqlNode*) ;
    }

  if (inNode == aql->checkEnv->currFromLoc)	/* at bottom of decls */
    aql->checkEnv->pCurrDecl = &inNode->nxt ;

  if (inNode->type == nFROM_LOC)
    aql->checkEnv->currFromLoc = 0 ;

  if (inNode->nxt && 
      (inNode->type == nLOC_LOCAL_POS   || inNode->type == nLOC_LOCAL_TAG_NAME ||
       inNode->type == nLOC_FOLLOW_POS  || inNode->type == nLOC_FOLLOW_TAG_NAME ||
       inNode->type == nLOC_LOCAL_FIELD || inNode->type == nLOC_LOCAL_FIELD_NAME || /* this line added (fw-980715) */
       inNode->type == nOBJECT || inNode->type == nLOC_METHOD))
    inNode->nxt = 0 ;	/* no more multiple LOC derivatives - see above */
} /* postDecls */


static BOOL preScope (AqlNode   *inNode, 
		              AQL        aql,
                      int        level)
     /* used as preFunc parameter for aqlTraverse,
	where it is executed on a node, before traversing its branches 
	- as used in aqlTraverse, returning TRUE means that no further traversal
	  of the branches and no execution of the postFunc will be performed.
      */
{
  AqlScope *scope ;
  int varNameDictIndex=0 ;

  /* declaration of variables */

  switch (inNode->type)
    {
    case nFROM_LOC:
      /* FROM_LOC->left  == locator,class
       * FROM_LOC->right == bool condition
       * FROM_LOC->nxt is another FROM_LOC  */
    case nFROM_TABLE:
      /* FROM_TABLE->left  == table
       * FROM_TABLE->right == bool condition
       * FROM_TABLE->nxt is another FROM_TABLE  */
    {
        AqlNode *naming_outer_node = inNode->left;
        AqlNode *naming_inner_node = naming_outer_node->nxt;
#if 0
        AqlNode *reference_node = naming_inner_node;
#else
        AqlNode *reference_node = naming_outer_node;
#endif
        char *reference_name = reference_node ? reference_node->name : NULL;
        int reference_number = reference_node ? reference_node->number : 0;
        AqlNode *locatorDefinitionNode = NULL ;
/* shouldn't we start that at NULL, to be on the safe side?
   What if we never set it / never find one -- or are we
   sure to? */
        int scopeDictIndex ;
	
        if (!dictAdd (aql->checkEnv->currentLocScope->dict, inNode->name, &varNameDictIndex))
        {
            /* e.g. generated by  >select all from a in class movie, a in class person<
               or in a more hidden way by the query >select all a,b< , because two variables are being
               declared (i.e. they both get a '_DEF' LOC_VAR), which is then picked out here */
            if (strcmp(inNode->name, "_DEF")==0)
                aqlError (aql, 988, inNode->pos, "only one declaration allowed for default iterator");
            else
                aqlError (aql, 989, inNode->pos, "multiple declaration of variable %s in same scope", inNode->name) ;
        }
	
        /* do we have it already?  if so point to that, identifying prefixes */
	
        /* first process node to left, out of order, because we need that */
        preScope (naming_outer_node, aql, level+1) ;	/* NB - non-recursive */

        if (aql->parseDebugLevel >= 4) {
            aceOutf(aql->debug_out, "Looking in scopes for locator matching name \"%s\" to suit node:\n", reference_name);
            aqlNodeOut(aql->debug_out, inNode, level, "locandum", FALSE);
            aqlNodeOut(aql->debug_out, naming_outer_node, level, "naming_outer_node", FALSE);
            aqlNodeOut(aql->debug_out, naming_inner_node, level, "naming_inner_node", FALSE);
            aqlNodeOut(aql->debug_out, reference_node, level, "reference_node", FALSE);
            aceOutf(aql->debug_out, "The name from the inner node is \"%s\" and the name from the outer node is \"%s\" and the reference name is \"%s\"\n",
		    naming_inner_node && naming_inner_node->name ? naming_inner_node->name : "NULL", naming_outer_node && naming_outer_node->name ? naming_outer_node->name : "NULL" , reference_name ? reference_name : "NULL");

            {
                Stack decls = aql->checkEnv->declStack;
                int top_decl = stackMark(decls);
                int decl_i;

                AqlNode **pdecl;

                for (pdecl = stackNext (decls, AqlNode**), decl_i = 0;
                     !stackAtEnd(decls);
                     pdecl = stackNext (decls, AqlNode**), decl_i++)
                {
                    aceOutf(aql->debug_out, "decl %d:\n", decl_i);
                    if (pdecl && *pdecl)
                    {
                        aqlNodeOut(aql->debug_out, *pdecl, level, "stackeddecl", FALSE);
                    }
                }
                stackCursor (decls, top_decl);
/* int    stackMark (Stack s)
   void   stackCursor (Stack s, int mark)
   TYPE   stackNext (Stack s, TYPE)
   double stackDoubleNext(Stack s)
   char*  stackNextText (Stack s) 
   BOOL   stackAtEnd (Stack s)
*/
            }
        }

        for (scope = aql->checkEnv->currentLocScope ; scope ; scope = scope->parent)
        {
	  for (scopeDictIndex = 1 ; scopeDictIndex <= dictMax (scope->dict) ; scopeDictIndex++)
            {
	      if (array(scope->loc, scopeDictIndex, AqlLoc*)) /* miewg, oct 2004, array not arr */
                {
                    locatorDefinitionNode = arr(scope->loc, scopeDictIndex, AqlLoc*)->definition ;
		
                    if (!locatorDefinitionNode)
                        /* it is possible for this to be NULL, but I'm unsure of the
                         * exact implications of this, so we better check (fw) */
                        continue;

                    if (naming_outer_node->loc && naming_outer_node->loc == locatorDefinitionNode->left->loc)
                    { 
                        AqlNode *yy = locatorDefinitionNode->left->nxt ;
		    
                        if (!reference_node && !yy)
                        {
                            goto done ;
                        }
                        if (!reference_node || !yy) /* allows explicit redefinition */
                        {
                            continue ;
                        }

                        if (reference_node->type == yy->type)
                        {
                                /* nodes that have the same type 
                                 * and share either their name or the 
                                 * number can have the same locator */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   I think that is the false assumption

   So what information do we have, that would let us make a better decision?

   We should probably never re-use a locator at a definition point.
   Does this cover all cases correctly?
   I'm pretty sure it's the right starting point, at least.

*/
                            if (reference_name && strcasecmp (reference_name, yy->name) == 0)
                            {
                                goto done ;
                            }

                            if (reference_number && (reference_number == yy->number))
                            {
                                goto done ;
                            }
                        }
                    }
                }
            }
        }
        locatorDefinitionNode = NULL ;
	
    done:
        if (locatorDefinitionNode 
            /* although we've found an existing locator to assign to this node
             * we don't use it if it is already a row-variable. This will create
             * separate locators for column-variables that iterate over
             * different columns on the same table.
             * See SANgc05250 for an explanation of the problem. */
            && !(locatorDefinitionNode->loc->source == IS_ROW))
        {

            if (aql->parseDebugLevel >= 4) {
                aceOutf(aql->debug_out, "Found existing loc %lx\n", locatorDefinitionNode->loc);
                aceOutf(aql->debug_out, "Assigning the locator already defined at node:\n");
                aqlNodeOut(aql->debug_out, locatorDefinitionNode, level, "existingloc", TRUE);
                aceOutf(aql->debug_out, "To be also the locator of node:\n");
                aqlNodeOut(aql->debug_out, inNode, level, "inNode", FALSE);
            }

            inNode->loc = locatorDefinitionNode->loc ;

            array(aql->checkEnv->currentLocScope->loc, varNameDictIndex, AqlLoc*) = inNode->loc;
        }
        else
        { 
            /* make new locator for inNode */
            inNode->loc = (AqlLoc*)halloc (sizeof(AqlLoc), aql->query->handle) ;
            inNode->loc->definition = inNode ;
            if (inNode->type == nFROM_TABLE)
                inNode->loc->source = IS_ROW ;

            array(aql->checkEnv->currentLocScope->loc, varNameDictIndex, AqlLoc*) = inNode->loc;
            if (aql->parseDebugLevel >= 4) {
                aceOutf(aql->debug_out, "Did not find existing loc, making new one at %lx for node:\n", inNode->loc);
                aqlNodeOut(aql->debug_out, inNode, level, "newlylocated", FALSE);
            }
        }
    }
      break;



    case nTABLE_ASSIGN:
      /* assignment of table variables */

      if (dictFind (aql->checkEnv->currentTableScope->dict, inNode->name, &varNameDictIndex))
	/* make new scope only if new defn shadows a former defn */
        { 
	  scope = aqlScopeCreate (aql->query->handle) ;
	  scope->parent = aql->checkEnv->currentTableScope ;
	  aql->checkEnv->currentTableScope = scope ;
	  inNode->flags |= AQL_SWITCH_SCOPE ;
	}

      /* do not merge - prefix merging is not so meaningful here */
      inNode->loc = (AqlLoc*)halloc (sizeof(AqlLoc), aql->query->handle) ;
      inNode->loc->definition = inNode ;

      dictAdd (aql->checkEnv->currentTableScope->dict, inNode->name, &varNameDictIndex) ;
      array(aql->checkEnv->currentTableScope->loc, varNameDictIndex, AqlLoc*) = inNode->loc;
      break;




    case nVAR_ASSIGN:
      /* assignment of scalar variables */

      if (dictFind (aql->checkEnv->currentScalarScope->dict, inNode->name, &varNameDictIndex))
	/* make new scope only if new defn shadows a former defn */
        { 
	  scope = aqlScopeCreate (aql->query->handle) ;
	  scope->parent = aql->checkEnv->currentScalarScope ;
	  aql->checkEnv->currentScalarScope = scope ;
	  inNode->flags |= AQL_SWITCH_SCOPE ;
	}

      /* do not merge - prefix merging is not so meaningful here */
      inNode->loc = (AqlLoc*)halloc (sizeof(AqlLoc), aql->query->handle) ;
      inNode->loc->definition = inNode ;

      dictAdd (aql->checkEnv->currentScalarScope->dict, inNode->name, &varNameDictIndex) ;
      array(aql->checkEnv->currentScalarScope->loc, varNameDictIndex, AqlLoc*) = inNode->loc;
      break;



    case nLOC_VAR:
    case nTABLE_VAR:
    case nVAR:
      /* uses of variables */

      if (inNode->type == nLOC_VAR)
	scope = aql->checkEnv->currentLocScope;
      else if (inNode->type == nTABLE_VAR)
	scope = aql->checkEnv->currentTableScope ;
      else
	scope = aql->checkEnv->currentScalarScope ;


      /* find the name of this LOC_VAR node in this scope or its parent scopes */
      while (!dictFind (scope->dict, inNode->name, &varNameDictIndex) && scope->parent)
	scope = scope->parent ;

      if (scope->parent || 
	  (scope && varNameDictIndex > 0)) /* added by (fw-980813) to cope with context vars in outermost scope */
	{
	  /* name found -> map the named node to the locator */
	  inNode->loc = arr(scope->loc, varNameDictIndex, AqlLoc*) ;
	}
      else
	/* can't find inNode->name in existing scopes */
	{
	  /* make a new scope-entry with this name */
	  inNode->loc = (AqlLoc*)halloc (sizeof(AqlLoc), aql->query->handle) ;

	  dictAdd (scope->dict, inNode->name, &varNameDictIndex) ;
	  array(scope->loc, varNameDictIndex, AqlLoc*) = inNode->loc;
	}
      break;


      
    case nTABLE_SFW:
      /* loc scope boundary */
      scope = aqlScopeCreate (aql->query->handle);
      scope->parent = aql->checkEnv->currentLocScope;
      aql->checkEnv->currentLocScope = scope;
      break;

    default:;
    } /* end switch */

  return FALSE ;
} /* preScope */



static void postScope (AqlNode *inNode, 
		               AQL      aql,
                       int      level)	
/* this function is called when the recusive traversal is on its way back UP the tree */
/* pop scopes where necessary */
{
  switch (inNode->type)
    {
    case nTABLE_SFW:
      aql->checkEnv->currentLocScope = aql->checkEnv->currentLocScope->parent ;
      break;


    case nTABLE_ASSIGN:
      if (inNode->flags & AQL_SWITCH_SCOPE)
	aql->checkEnv->currentTableScope = aql->checkEnv->currentTableScope->parent ;
      break;

    case nVAR_ASSIGN:
      if (inNode->flags & AQL_SWITCH_SCOPE)
	aql->checkEnv->currentScalarScope = aql->checkEnv->currentScalarScope->parent ;
      break;

    default:;
    }
} /* postScope */


/******************************************************************/

static AqlNode *aqlCopyNode (AqlNode *srcNode, 
			     AQL      aql)
     /* recursive, it'll copy the whole tree below the
	initial node given. */
     /* this function is currently restricted to copy values
	that are of any use during check1. */
{
  AqlNode *newNode;

  if (!srcNode) return 0;

  newNode = (AqlNode*)halloc (sizeof(AqlNode),
			      aql->query->handle) ;
  
  newNode->type = srcNode->type;
  newNode->name = strnew(srcNode->name, aql->query->handle);
  newNode->op = srcNode->op;
  newNode->number = srcNode->number;
  newNode->flags = srcNode->flags;
  newNode->vtype = srcNode->vtype;
  newNode->value = srcNode->value;
  newNode->pos = srcNode->pos;

  newNode->left = aqlCopyNode(srcNode->left, aql);
  newNode->right = aqlCopyNode(srcNode->right, aql);
  newNode->nxt = aqlCopyNode(srcNode->nxt, aql);

  return newNode;
} /* aqlCopyNode */


/********************************************************************
 * check from-clauses and handle empty declarations
 ********************************************************************/

static BOOL preLocCheck (AqlNode *inNode, 
			 AQL      aql,
                         int      level)
{
  if (inNode->type == nFROM_LOC)
    {
      if (!inNode->left)
	aqlError (aql, 601, inNode->pos, "preLocCheck() - inNode->left == NULL in FROM_LOC-node") ;

      switch (inNode->left->type)
	{
	case nLOC_VAR:
	case nLOC_FOLLOW_TAG: case nLOC_FOLLOW_TAG_NAME:
	case nLOC_FOLLOW_POS:
	case nLOC_LOCAL_TAG: case nLOC_LOCAL_TAG_NAME:
	case nLOC_LOCAL_POS:
	case nLOC_LOCAL_FIELD: case nLOC_LOCAL_FIELD_NAME:
	case nLOC_METHOD:
	case nLOC_TABLE_FIELD: case nLOC_TABLE_FIELD_NAME:
	  if (inNode->left->nxt == NULL)
	    /* the locator-variable FROM_LOC declaration 
	     * has an empty definition, i.e. is based on another existing definition */
	    {
	      AqlNode *originNode = 0 ;
	      
	      if (!inNode->left->loc)
		aqlError (aql, 602, inNode->pos,
			  "preLocCheck() - a LOC_xxx node has inNode->loc == NULL") ;

	      /* the FROM_LOC is based upon a LOC_xxx,
	       * which has no definition, just a reference to
	       * another FROM_LOC. (this is the ->left) */
          
          {
              AqlNode *definition = inNode->left->loc->definition;

              /* this shouldn't really be NULL, but some got through the parser so let's play safe */
              if (definition == NULL) {
                  aqlError (aql, 961, inNode->pos,  "Malformed FROM in SelectFromWhere");
              } else {
                  originNode = definition->left;
              }
          }
	      
	      inNode->left = aqlCopyNode (originNode, aql);
	      /* NOTE: that we need to re-scope the entire tree after we've 
	       *    copied nodes, so the new nodes get locators set up etc. */
	    }
	  break;
	default:;
	} /* end-switch inNode->left->type */
    } /* end-if inNode->type == FROM_LOC */
  
  return FALSE;
} /* preLocCheck */


/********************************************************************/
/**** handle tables/row variables etc.                           ****/
/********************************************************************/

static BOOL preTableProcess (AqlNode   *inNode, 
			                 AQL        aql,
                             int        level)
{
  if (inNode->type == nLOC_VAR)
    {
      /* the LOC_VAR inNode may not have the ->loc field assigned
	 This may happen when queries are nested.
	 Remember: because they are parsed top-down rather than inside out. */
      if (inNode->loc && inNode->loc->source && inNode->loc->source == IS_ROW)
	{
	  AqlNode *followNode = inNode->nxt ;
	  
	  if (!followNode || 
	      !(followNode->type == nLOC_LOCAL_FIELD))
	    {
	      /* generated by e.g. >@t := select a,a->cast from a in class movie; select a from a in @t<
		 where a is a row variable, but select fields have to be scalars
		 so we should write >select a:1 from a in @t< to select the first column of a row a
		 NOTE: if @t was a single column table array, we might
		 want to allow the above form without the :1 */
	      if (followNode &&
		  (followNode->type == nLOC_FOLLOW_TAG_NAME ||
		   followNode->type == nLOC_FOLLOW_POS))
		{
		  /* select a->cast from a in @table --or--
		   * select a->2 from a in @table */
		  aqlError (aql, 933, followNode->pos,
			    "Cannot dereference a row-based variable directly\n"
			    "Declare a column-variable which can then be dereferenced :\n"
			    "Use something like %s:1->%s instead!", strcmp(inNode->name, "_DEF") == 0 ? "" : inNode->name,
			    followNode->type == nLOC_FOLLOW_TAG_NAME ? followNode->name : messprintf("%d", followNode->number)) ;
		}

	      aqlError (aql, 932, inNode->pos, "%s is row based and must specify the column \n"
			"that is to be extracted from the iterating rows :\n"
			"Use something like %s:1 to refer to the first column in the row-variable %s",
			strcmp(inNode->name, "_DEF") == 0 ? "default iterator" : messprintf("variable %s", inNode->name),
			inNode->name, inNode->name) ;
	    }
	  
	}
    }
  


  /* check that a scalar var has been defined */
  if (inNode->type == nVAR)
    {
      if (inNode->loc)
	if (inNode->loc->definition || 
	    inNode->loc->isContext) /* we'll get the definition later, so don't worry now */
	  return FALSE;

      /* no locator, or missing def-node */
      aqlError (aql, 990, inNode->pos+1, "scalar variable %s undefined", inNode->name);
    }
  
  /* check that a table var has been defined */
  if (inNode->type == nTABLE_VAR)
    {
      if (inNode->loc)
	if (inNode->loc->definition || 
	    inNode->loc->isContext) /* we'll get the definition later, so don't worry now */
	  return FALSE;

      /* no locator, or missing def-node */
      aqlError (aql, 991, inNode->pos+1, "table variable %s undefined", inNode->name);
    }

  /* check that a locator var has been defined */
  if (inNode->type == nLOC_VAR)
    {
      if (inNode->loc)
	if (inNode->loc->definition)
	  return FALSE;

      /* no locator, or missing def-node */
      if (strcmp(inNode->name, "_DEF")==0)
	{
	  if (inNode->nxt)
	    aqlError (aql, 993, inNode->nxt->pos, "no unnamed iterator declaration");
	  else
	    /* position of error impossible to determine */
	    aqlError (aql, 993, 0, "no unnamed iterator declaration");
      }

      aqlError (aql, 992, inNode->pos+1, "iterator variable %s undeclared", inNode->name);
    }
  

  return FALSE ;		/* tell aqlTraverse to continue down the tree */
} /* preTableProcess */
		


/* (1) Check a node against validity for a given table,
       such nodes are for instance sort fields in the order_by construct 
       or the column name/number specified in a row variable 'select t:2 from t in @table'
   (2) Transform all xxx_FIELD_NAME nodes into xxx_FIELD nodes by
       mapping the field-name to a column name in the table specified
   NOTE: all LOC_TABLE_FIELD and *_NAME have already been converted to
   their equivalent structure using the LOC_LOCAL_FIELD and row-variables
*/
static void checkField (AqlNode   *fieldNode, 
			AqlNode   *tableDefNode,
			AQL        aql)
{
  int      columnNumber;
  AqlNode *selectFieldNode;

  if (tableDefNode->type != nTABLE_SFW)
    aqlError (aql, 650, 0,
	      "checkField() - can only check fields against a TABLE_SFW-type node") ;

  columnNumber = numberOfTableColumns(tableDefNode, aql);

  /* field selection by number */
  if (fieldNode->type == nSORT_FIELD || 
      fieldNode->type == nLOC_LOCAL_FIELD ||
      fieldNode->type == nEXPR_TABLE_FUNC_FIELD)
    /* the node specifies the column-number, so we just check the number for validity */
    { 
      if (fieldNode->number < 1 || fieldNode->number > columnNumber)
        aqlError (aql, 922, fieldNode->pos+1, "field number %d out of range", fieldNode->number) ;
    }

  /* field selection by name -> transform to number selection fields */
  else if (fieldNode->type == nSORT_FIELD_NAME || 
	   fieldNode->type == nLOC_LOCAL_FIELD_NAME ||
	   fieldNode->type == nEXPR_TABLE_FUNC_FIELD_NAME)
    { 
      int columnCounter ; 
      char *columnName ;
      BOOL isFound = FALSE;

      /* if a field is given by its name we have to find this name
	 in the definition of this table */
      
      for (selectFieldNode = tableDefNode->right, columnCounter = 0; selectFieldNode; selectFieldNode = selectFieldNode->nxt, ++columnCounter)
	{
	  columnName = 0;

	  if (selectFieldNode->name == NULL && selectFieldNode->op == 0)
	    continue;

	   if (selectFieldNode->op)
	     {
	       if (selectFieldNode->name && selectFieldNode->name[0] != '_')
		 columnName = selectFieldNode->name; /* explicit rename of field overrides name of expr_op */
	       else
		 columnName = aqlOpTypeName(selectFieldNode->op);
	     }
	   else if (selectFieldNode->name[0] != '_') /* field-node has a proper name */
	    columnName = selectFieldNode->name;
	   else
	    {
	      /* field has an artificial name, e.g. _DEF, _1, _2 etc...
		 that means it is a more complex, in case of FOLLOW_TAG_NAME
		 and LOCAL_TAG_NAME, we want the last tag in the chain to be the columnName */
	      /*
		 @t := select ->Cast->Full_name from class movie; select all @t:Full_name
		 @t := select ->Address[Email] from class person; select all @t:Email
	      */
	      /*
		 Example parsetree, to show how we'd get to the named node :-
		  left   TABLE_SFW name=0
		  ...
		    next   FROM_LOC name=_3             (loc=140545e20)	<-- definition of selectField-def->left
		      left   LOC_VAR name=_2            (loc=140545da0)
		      next   LOC_FOLLOW_TAG_NAME name=Full_name		<-- the node whose name we want
		    next   FROM_LOC name=_4             (loc=140545ea0) <-- definition of the selectField
		      left   LOC_VAR name=_3            (loc=140545e20)
		      next   LOC_LOCAL_POS number=1
		    right  LOC_VAR name=_4              (loc=140545ea0)	<-- selectField
	      */
	      if (selectFieldNode->loc->definition)
		if (selectFieldNode->loc->definition->left)
		  if (selectFieldNode->loc->definition->left->loc)
		    if (selectFieldNode->loc->definition->left->loc->definition->left)
		      if (selectFieldNode->loc->definition->left->loc->definition->left->nxt)
			columnName = selectFieldNode->loc->definition->left->loc->definition->left->nxt->name;
	    }

	  if (columnName && strcasecmp(fieldNode->name, columnName)==0)
	    { 
	      /* set the column number in the way the user would have
		 specified it, i.e. leftmost column = 1, the mapping
		 back into C-numbers (starting from 0) is taken care of 
		 later in loopVariables(), when loc->column is assigned */
	      fieldNode->number = columnCounter + 1 ;
	      
	      if (fieldNode->type == nSORT_FIELD_NAME) 
		fieldNode->type = nSORT_FIELD ;
	      else if (fieldNode->type == nLOC_LOCAL_FIELD_NAME)
		fieldNode->type = nLOC_LOCAL_FIELD ;
	      else if (fieldNode->type == nEXPR_TABLE_FUNC_FIELD_NAME)
		fieldNode->type = nEXPR_TABLE_FUNC_FIELD ;
	      
	      isFound = TRUE;
	      
	      break ;		/* out of for-loop */
	    }
	} /* end for */

      if (!isFound)
	aqlError (aql, 923, fieldNode->pos+1, "field name %s not in table", fieldNode->name) ;
    }
  else
    aqlError (aql, 651, 0,
	      "checkField() - invalid type of fieldNode") ;

  return ;
} /* checkField */


static int numberOfTableColumns (AqlNode   *inNode,
				 AQL        aql) /* recursive */
{
  AqlNode *selectFieldNode;
  int colNum;

  switch (inNode->type)
    {
    case nTABLE_SFW:
      for (selectFieldNode = inNode->right, colNum = 0 ; selectFieldNode ; selectFieldNode = selectFieldNode->nxt, ++colNum);
      return colNum;

    case nTABLE_VAR:
      return numberOfTableColumns(inNode->loc->definition, aql);

    case nTABLE_ASSIGN:
    case nTABLE_ORDER:
      return numberOfTableColumns(inNode->left, aql);

    case nTABLE_OP:
      if (inNode->vtype != 'T')	/* set up by postTableProcess() - see below */
	aqlError (aql, 660, 0,
		  "numberOfTableColumns() - TABLE_OP node not yet initialised with correct vtype 'T'");

      return inNode->value.T->ncol;

    default:
	aqlError (aql, 661, 0, "numberOfTableColumns() - unknown tableNode type");

    }

  return 0;			/* return something for compiler happiness */
}


static AqlResultDisposition aqlNeedsSorting(AqlNode *inNode, AqlNode **psortNode)
{
#if 0
    AqlNode *n;
    for (n = inNode; n != NULL; n = n->up)
    {
        switch (n->type)
        {
        case nTABLE_ORDER:
            *psortNode = n;
            return dTableForSorting;
        case nEXPR_TABLE_FUNC: return dOutputAlmostImmediately;
        }
    }
    return dOutputAlmostImmediately;
#else
    return dTableForSorting;
#endif
}

static void postTableProcess (AqlNode   *inNode, 
			                  AQL        aql,
                              int        level)
/*
   process nodes of type
   TABLE_ORDER      - check for presence and correctness of sortfields
   TABLE_OP         - first check for combinable column-counts between the two tables
   EXPR_TABLE_FUNC  - the table that is operated upon must not have more than 1 column
   LOC_VAR (IS_ROW) - check select-columns and process NAMEs to the node->number's
   SELECT_FIELD     - absorb node and re-attach
   TABLE_SFW        - attach select-fields on the right, init table value & set columnNames
*/
{
  switch (inNode->type)
    {
    case nTABLE_ORDER:
      /* check that the specified sort field actually exist for that table */
      { 
	AqlNode *sortFieldNode;	/* these are the order_by :name,:2 fields */
	
	/* the sort fields are a linked-list starting with TABLE_ORDER->right */
	/* check the fields for validity against the given table and
	   transform SORT_FIELD_NAMEs into SORT_FIELD by finding their column number */
	for (sortFieldNode = inNode->right ; sortFieldNode ; sortFieldNode = sortFieldNode->nxt)
	  checkField (sortFieldNode, inNode->left, aql) ; /* inNode->left is the TABLE_SFW-node */

	inNode->vtype = 'T';
      }
      break;
      
  
    case nTABLE_OP:
      /* check that the two tables are "combinable" 
	 by checking that their number of columns is the same */
      { 
	/* TABLE_OP is one of UNION, INTERSECT or DIFF */
	char    *typeString;
	int      leftColNum, rightColNum, col;
	
	/* count number of columns in both tables */
	leftColNum = numberOfTableColumns(inNode->left, aql);
	rightColNum = numberOfTableColumns(inNode->right, aql);

	if (leftColNum != rightColNum)
	  aqlError (aql, 920, inNode->pos, "table %s operator - column number mismatch", aqlOpTypeName(inNode->op)) ;
	
	/* NOTE: we can't check for type compatability yet,
	   as we haven't eveluated any FROM_LOCs yet to know their types */

	/* now make the table and it's types as uninitialised */
	typeString = (char*)messalloc (leftColNum+1);
	
	/* set column-types to be yet unknown, i.e. '0' as recognised by the table-package */
	for (col = 0 ; col < leftColNum; ++col)
	  typeString[col] = '0';
	typeString[leftColNum] = '\0'; /* string terminator for column-types string */

	/* All tables are created on the query's own handle,
	 *   which dies when the query execution is finished. 
	 * Before the handle is destroyed, 
	 *  the final result table will be copied.
	 */
	inNode->value.T = tableHandleCreate (100, typeString, aql->query->handle) ;

	inNode->vtype = 'T';

	messfree(typeString);
	break;
      }


    case nEXPR_TABLE_FUNC:
    case nEXPR_TABLE_FUNC_FIELD:
    case nEXPR_TABLE_FUNC_FIELD_NAME:
      /* the table (inNode->left) must not have more than one column
       * the table is parsed as 'safe_table' and can therefore be a
       * TABLE_SFW-node or a a TABLE_VAR */
      /*
       * Examples of two erroneous queries that return error No 921:
       *  (1) $a := min(select p, p->height from p in class person)
       *  (2) @t := select p, p->height from p in class person; $a := min @t
       * whereas the form with column selection is correct
       *  (1) $a := min(select p, p->height[imperial] from p in class person):2
       *  (2) @t := select p, H::p->height[imperial] from p in class person; $a := min @t:H
       */
      /* EXAMPLE: In order to count the number of rows in a table
       * it is wise to restrict the number of columns to just 1.
       * This greatly reduces confusion; for example
       *       select m, count(select m, m->cast, m->director) from m in class movie
       * The count function returns a useless number, 
       * because the result table is a cartesian product of
       * the three columns. Look at the would-be result of the nested table:
       *	Heat    Al Pacino       Michael Mann
       *	Heat    Al Pacino       Alan Smithee
       *	Heat    Robert De Niro  Michael Mann
       *	Heat    Robert De Niro  Alan Smithee
       *	Heat    Val Kilmer      Michael Mann
       *	Heat    Val Kilmer      Alan Smithee
       *	Heat    John Voight     Michael Mann
       *	Heat    John Voight     Alan Smithee
       *	Heat    Tom Sizemore    Michael Mann
       *	Heat    Tom Sizemore    Alan Smithee
       * which gives the count-result of 10 for the movie 'Heat', i.e. rubbish.
       */
      {
	AqlNode *tableDefNode = 0 ;
	
	/* find the TABLE_SFW node that generates the table we operate on */

	if (inNode->left->type == nTABLE_SFW)
	  tableDefNode = inNode->left;

	else if (inNode->left->type == nTABLE_ORDER)
	  tableDefNode = inNode->left->left;

	else if (inNode->left->type == nTABLE_VAR)
	  /* the locator of a TABLE_VAR node is created by a TABLE_ASSIGN-node
	     whose ->left neighbor is the TABLE_SFW-node that generates the table */
	  /* NOTE: the definition node is guaranteed to exist, see aqlError No 991 */
	  {
	    tableDefNode = inNode->left->loc->definition->left;
	    
	    if (tableDefNode->type == nTABLE_ORDER)
	      tableDefNode = tableDefNode->left;
	  }
	else
	  aqlError (aql, 671, 0,
		    "postTableTypeCheck() - tableFunc must only operate on a 'safe_table'");


	/* The syntax allows the column of the 'safe_table' to be specified, so there are three cases :
	 *   (1) TABLE_FUNC              --> no column selection, check for only ONE column
	 *   (2) TABLE_FUNC_FIELD        --> column selection by number check for validity
	 *   (3) TABLE_FUNC_FIELD_NAME   --> column selection by name, check name
	 */

	if (inNode->type == nEXPR_TABLE_FUNC)
	  /* no column selection, check for one column only */
	  {
	    if (numberOfTableColumns(tableDefNode, aql) != 1)
	      aqlError (aql, 921, tableDefNode->pos, "table-argument for tableFunc expression must generate ONLY 1 column");
	  }
	else
	  /* we have a column selection field, check field selection */
	  checkField(inNode, tableDefNode, aql);
      }	/* end-case EXPR_TABLE_FUNC_xx */
      break;




    /* row variables */
    case nLOC_VAR:		
      if (inNode->loc && inNode->loc->source == IS_ROW)
	/* test row variables that specify the column selection 
	   by name (LOC_LOCAL_FIELD_NAME) or number (LOC_LOCAL_FIELD) */
	{ 
	  AqlNode *locatorDefinitionNode = inNode->loc->definition ; /* this is a FROM_TABLE */
	  AqlNode *tableDefNode;
	  
	  if (!locatorDefinitionNode)
	    aqlError (aql, 675, 0,
		      "postTableProcess()->case nLOC_VAR: locatorDefinitionNode is NULL, \n"
		      "inNode can't be a LOC_VAR, if it is defined - see aqlError No 992");

	  if (locatorDefinitionNode->left->loc->isContext)
	    /* the context var has no definition node, 
	     * only aqlCheck2() will get to know its value, don't worry for now */
	    break;


	  /* look up the TABLE_SFW node in which the table is generated */
	  /* NOTE, locatorDefinitionNode->left is a TABLE_VAR, */
	  /*       hence locatorDefinitionNode->left->loc->definition is a TABLE_ASSIGN,
		   i.e. where the table variable was assigned from  and its ->left is the TABLE_SFW */
	  /* NOTE the definition node (TABLE_ASSIGN) is guaranteed to exist by aqlError No 991 */

	  if (!locatorDefinitionNode->left->loc->definition) 
	    aqlError (aql, 676, 0,
		      "postTableProcess()->case nLOC_VAR: TABLE_VAR->definition is NULL, \n"
		      "can't get here if table-var variable is defined - see aqlError No 991");

	  tableDefNode = locatorDefinitionNode->left->loc->definition->left;

	  if (tableDefNode->type == nTABLE_ORDER)
	    tableDefNode = tableDefNode->left;
	  
	  checkField (inNode->nxt, tableDefNode, aql) ;
	}
      break;



    case nSELECT_FIELD:
      /* can absorb SELECT_FIELD once LOC_VAR have been processed 
       * and any LOC_VARs around don't have any ->nxt pointers
       * !! NB post-order important - has to be done on the way UP the tree !! */

      /* a select-field can have a specific name, 
       * if the 'select email::p->Author->Email from ...' alias construct is used 
       * is this case we assign that alias-name to the locator-variable attached to this select-field */
      if (inNode->name) 
	/* may use var name if no field name given and real, see below */
        inNode->left->name = inNode->name ;

      if (inNode->nxt)		
	/* chain the locator of the next select-field to the locator of this select-field
	 * This will remove the SELECT_FIELD from the linkage in the tree. */
	inNode->left->nxt = inNode->nxt->left ;

      break;



    case nVAR_ASSIGN:
      /* the scalar value needs to be put in a single-element table, which we create upon the resultHandle store  */
      {
	char typeString[2];	/* a copy will be allocated by tableHandleCreate */
	
	/* make the typeString for one column */
	typeString[0] = '0';	/* type yet unknown */
	typeString[1] = '\0';
	
	/* make a one-row table */

	/* we create all tables on the query's own handle,
	   which dies when the query is finished. One final result table
	   however will survive, as it will be moved onto a different
	   handle before this one is destroyed */
	inNode->value.T = tableHandleCreate (1, typeString, aql->query->handle);

	inNode->vtype = 'T';
      }
      break;
      


    case nTABLE_SFW:
      /* create a table to be attached to the inNode of type TABLE_SFW and set the column names */
    { 
        AqlNode *selectFieldNode;
        char    *typeString;
        int      colNum, col;
	
        /* ->right are the where clauses */
        inNode->right = inNode->right->left ; /* SELECT_FIELD absorbed - see above
                                                 (which is already done, of course!) */

        colNum = numberOfTableColumns(inNode, aql);

        if (colNum < 1)
            /* can't happen syntactically. The parser expects select fields
	     * as specified in the yacc-syntax */
	    aqlError (aql, 678, 0,
		      "postTableProcess() - impossible query, %d select-fields in TABLE_SFW-node",
		      colNum) ;

        /* make a name for the table that will be generated by the select-from-where 
           the characters of the table name denote the types of the values in each column
           at this stage the types are unknown and the types are '0' (see makeNewNode() and zMake()) */

        /* allocate the memory for a string of that length and a NULL-string-terminator */
        typeString = (char*)messalloc (colNum+1) ;
	
        /* set column-types to be yet unknown, i.e. '0' as recognised by the table-package */
        for (col = 0 ; col < colNum; ++col)
            typeString[col] = '0';
        typeString[colNum] = '\0'; /* string terminator for column-types string */
	
	
        /* initialise the table-value for the results */
        {
            AqlNode *sortNode = NULL;
            inNode->resultConsumer = (AqlResultConsumer*)halloc(sizeof(AqlResultConsumer), aql->query->handle);
            inNode->resultConsumer->sorting = (aqlNeedsSorting(inNode, &sortNode) == dTableForSorting);
            inNode->resultConsumer->resultTableDict = dictHandleCreate(100, aql->query->handle);
            inNode->resultConsumer->nRows = 0;
            inNode->resultConsumer->dictStringStack = stackHandleCreate(50, aql->query->handle);
            inNode->resultConsumer->sortNode = sortNode;
        }
        /* we create all tables on the query's own handle,
           which dies when the query is finished. One final result table
           however will survive, as it will be moved onto a different
           handle before this one is destroyed */
        inNode->value.T = tableHandleCreate (1024, typeString, aql->query->handle) ; 
	

        /* assign the names of the fieldNodes to the respective columns in the table */
        for (selectFieldNode = inNode->right, col = 0 ; selectFieldNode ; selectFieldNode = selectFieldNode->nxt, ++col)
            if (selectFieldNode->name)
            {
                if (selectFieldNode->name[0] == '_') /* temp. name */
                    tableSetColumnName (inNode->value.T, col, "_") ;
                else
                    tableSetColumnName (inNode->value.T, col, selectFieldNode->name) ;
            }

        inNode->vtype = 'T';

        messfree(typeString);
    }
      break;

    default:;
    } /* end switch(inNode->type) */
} /* postTableProcess */




/************************************************************************/
/************* Now the utilities: first the traversal routine ***********/
/***************************************** *******************************/

void aqlTraverse (AqlNode   *inNode,
                  AQL        aql, 
                  BOOL (*preFunc)(AqlNode*,AQL,int), 
                  void (*postFunc)(AqlNode*,AQL,int),
                  int level)
     /* traverse the parsetree
	preFunc  - a boolean function, that acts on a tree-node is 
	           called before traversing the tree further.
	           Returning TRUE will not traverse further.
	postFunc - a void function, that is called after the
	           tree past this node has been explored
     */
{
  BOOL isStopTraversal = FALSE;

  if (!inNode) return ;

  /**************************************/

  /* execute preFunc on the current node */
  if (preFunc) 
    isStopTraversal = (*preFunc)(inNode, aql, level);		/* allow truncation of search, if preFunc returns TRUE */

  /**************************************/

  if (isStopTraversal)
    return;

  /* traverse the tree further DOWN */
  if (inNode->left) 
    aqlTraverse (inNode->left, aql, preFunc, postFunc, level+1) ;

  if (inNode->right) 
    aqlTraverse (inNode->right, aql, preFunc, postFunc, level+1) ;

  if (inNode->nxt) 
    aqlTraverse (inNode->nxt, aql, preFunc, postFunc, level+1) ;

  /**************************************/

  /* execute postFunc on the current node */
  /* the postFunc is therefor executed when the recursive
     traversal is on the way UP the tree */
  if (postFunc)
    (*postFunc)(inNode, aql, level) ;

  return;
} /* aqlTraverse */


/* Set the "up" links in a query tree --
   do this here to avoid fiddling around
   in parser code where there might be lots
   of places in which to make these links.
 */

void aqlSetParents(AqlNode *inNode)
{

    if (!inNode) return ;

    if (inNode->left)
    {
        inNode->left->up = inNode;
        aqlSetParents(inNode->left);
    }

    if (inNode->right)
    {
        inNode->right->up = inNode;
        aqlSetParents(inNode->right);
    }

    if (inNode->nxt)
    {
        inNode->nxt->up = inNode;
        aqlSetParents(inNode->nxt);
    }
} 

/*********************** eof *********************************************/

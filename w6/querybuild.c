/*  File: querybuild.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: @(#)querybuild.c	1.9 6/3/96
 * Description:
 **  Constructs general queries.                              **
 ***    Gary Aochi modified version of querydisp.c
 * Exported functions:
 **  qbuildCreate
 * HISTORY:
 * Last edited: Dec 16 19:44 1998 (rd)
 * * Feb  7 12:33 1994 (mieg): moved downselectedTriangle etc to menu.c
 * * Feb 09 15:30 1993 (aochi): split off qbe to qbedisp.c and query_.h
 *                              and cleaned up some code in querydisp.c
 * * Jun 14 10:14:40 PDT 1993 (aochi) : split querydisp.c into 2 files:
 *                                      querydisp.c (ctrl), querybuild.c
 * * Jan 18 12:55 1993 (mieg): Changed MAGIC and added some static def
 * (original querydisp.c created Dec 30 19:05:47 1991 (mieg))
 * Created: Jun 14 10:14:40 PDT 1993 (aochi)
 *-------------------------------------------------------------------
 */

#include <ctype.h>
#include "acedb.h"
#include "lex.h"
#include "bs.h"         /* used for bsTagsInClass() procedure */
#include "query_.h"
#include "query.h"
#include "../wh/menu.h"
#include "systags.h"  /* need for _bsRight gha 6/9/93 */
#include "classes.h"  /* need for _bsRight gha 6/9/93 */
#include "tree.h"     /* for treeChooseTagFromModel() */

/* declaration because of GHA modifications to graphtext(Scroll)entry */
extern int graphTextScrollEntry (char* text, int len, int wlen, float x, float y, void (*fn)());
     

#define ENTRY_NUM   4    /* Number of fields (same as # of user-entry boxes) */
#define SELECT_OFF  3    /* offset from ATTR, COND, VALU, CONJ to triangle */
    /* Mapping of user-entry boxes/fields: */
#define ATTRIBUTE   0
#define CONDITION   1    /* Numbering sequence is important...do not change */
#define VALUE       2    /*  without carefully observing functionality. */
#define CONJUNCTION 3    /* Some code depends on CONJUNCTION being last. */
#define CLASSTAG   -1
#define PRECLASS   -2
#define ALL_STUFF  -3    /* used only for clearQbuilder */

#define RELATED       0x01 /* for display styles in query interface */
#define ALLTEXT       0x02
#define SIMPLE_UI     0x04
#define QUOTE_ON      0x08
#define TAG_SEL       0x10

#define SETFLAG(y)  qbuild->flagmask = (qbuild->flagmask | y)
#define TESTFLAG(y) (qbuild->flagmask & y)
#define RMFLAG(y)   qbuild->flagmask = (qbuild->flagmask & ~y)

typedef struct QBUILDSTUFF
{ int   magic;        /* == MAGIC*/
  Graph graph ;
  int class_box;           /* box containing preclass and classtag */
  char preclass_entry[BUFFER_SIZE]; /* is Find/>? or Follow/> or text srch */
  char preclass_syntax[BUFFER_SIZE];
  int pre_triBox;          /* for selector */
  char classtag_entry[BUFFER_SIZE];
  char classtag_syntax[BUFFER_SIZE];
  int class_triBox;        

  int user_box;        /* user-entry box area */
  Array entries;       /* text strings from user-entry boxes */
  Array syntax;        /* actual translations of entries */
  int entry_index;     /* index of current user-entry box */
  int entry_set;       /* user-entry box area */
  Array attrMenuArray; /* keep track of menu for each QB line */

  char flagmask;       /* control mask, 8 bits of togglable settings */
} *QBUILD ;

Graph qbuildGraph = 0 ;

char resbuffer[QBUFF_MULT*BUFFER_SIZE]; /* for forming query commands */

int qbuild_selected_class = -1 ; 
static int selected_classtag = -1; 

static int  MAGIC = 92984 ;

#include "display.h"	/* must come after QBUILD definition */
     
static void qbuildDisplay(void) ;
static void qbuildPick (int k) ;
void qbuildCreate (void) ;
     
#define QBUILDGET(name)  QBUILD qbuild ; \
                       \
                       if (!graphActivate (qbuildGraph)) \
                          messcrash ("%s can't find qbuildGraph", name); \
                       if (!graphAssFind (&MAGIC, &qbuild)) \
                           messcrash ("%s can't find graph", name) ; \
                       if (qbuild->magic != MAGIC) \
                           messcrash ("%s received a wrong pointer", name) ; \
                       displayPreserve()  
     
/* translate from index to a box # (index ~ 1:ENTRY_NUM box) */
static int INDEX2BOX(int box, int in) { return (ENTRY_NUM*in + box + 1); } 

/* mainly for qbuildPick(): translate box # to index.(box = ENTRY_NUM*index) */
static int BOX2INDEX(int box, int in) { return (in-box-1)/ENTRY_NUM; } 

static BOOL isSelector(int box, int in) 
{ 
  if (((in - box - 1) % ENTRY_NUM) == SELECT_OFF)
    return TRUE;
  else
    return FALSE;
}
    
/*********** stuff local to the qbuild routine and initialiser **********/
/*** key lists maintenance ***/

static BOOL attr_sub_model;  /* for # operator */
static KEYSET classtagkeyset;     /* keyset of meta-tags in class */
static KEYSET attrkeyset;    /* keyset of tags in current class */  
static Array classlist = 0;

static void qbuildMenuExec(void);

static int qbuildClassMax() {
  int i, k ;
  static int n ;
  
  if (arrayExists(classlist))
    return n;

  n = 0;
  classlist = arrayCreate(30, Array);   /* init classlist */
  
  for (i = 1 ; i <= 256 ; i++) {
    if (pickList[i].name
	&& ( 
	    (pickList[i].type == 'B' && iskey (pickList[i].model) == 2)
	    || (pickVisible(i) && pickList[i].type != 'B' )
	    ))
      {
      array(classlist, n, Array) = arrayCreate(BUFFER_SIZE, char);
      k = KEYMAKE(i,0) ; /* class() macro >> 24 for key to class translation */
      strcpy(ARR2STRING(classlist, n), className(k));
      n++;
    }
  }
  arraySort(classlist, textAlphaOrder);

  return n;
}

#define KEYLIST_NAME(offset, ks) name(arr(ks, offset, KEY)) /* get list_name */
     
/*  removes old keylist (a keyset) and finds a new list of
 *  keys for the class c...sets to global ***keyset.
 */
static void newKeyList(int c, KEYSET *k) 
{
  keySetDestroy(*k);
  *k = bsTagsInClass(c);
  if (*k)
    arraySort(*k, keySetAlphaOrder);   /* alphanum-order keyset */
}

/* returns max size of keyset */
static int keyListMax(KEYSET k) 
{
  if (k) 
    return keySetMax(k);
  else return 0;
}

/*****************************************/

/*
 * removes extra whitespace (spaces) from both ends of text string
 * and from the middle (between series of nonwhite chars)
 * string is modified in place
 * Tabs unlikely (parsed in widget for name completion)
 */
static void removeWhite(char *text) 
{
  char *cp = text - 1, *copy = text ;

  if (text) 
    while (TRUE)
      switch(*++cp)
	{
	case 0: /* done, char = NULL, but remove trailing spaces  */
	  if (copy > text && (*(copy - 1) == ' ')) copy-- ;
	  *copy = 0 ;
	  return ;
	case ' ': case '\t': /* spaces only one allowed */
	  if (copy == text)
	    break ; /* not at beginning */
	  if (*(copy - 1) == ' ')
	    break ; /* no repeats */
	  *copy++ = ' ' ; /* accept one */
	  break ;
	default:
	  *copy++ = *cp ;
	  break ;
	}
}

static BOOL isNumeric(char *text)
{
  char *cp = text;
  BOOL no_decimal_pt = TRUE;

  if (cp)
    while (*++cp) {
      if (!isdigit((int)*cp)) {
	if (*cp == '.' && no_decimal_pt)
	    no_decimal_pt = FALSE;
	  else
	    return FALSE;  /* only one decimal point is valid */
      }
    }
  return TRUE;
}
	  

/*
 * takes current contents of user-entry fields and constructs
 * query-expression in proper syntax.
 */
static char *makeQuery(char *result)
{
  int i;

  QBUILDGET("makeQuery");

  messStatus("Forming Query Command");
  strcpy(result, qbuild->preclass_syntax);
  strcat(result, qbuild->classtag_syntax);

  for (i = 0; i < arrayMax(qbuild->syntax); i++) { /* concatenate query */
    if ((i % ENTRY_NUM == CONJUNCTION)
	&&
	!strcasecmp(ARR2STRING(qbuild->entries, i), "END"))
      break;
    else {
    
      strcat(result, " ");  /* add a space between words */
      
      if ((i % ENTRY_NUM == VALUE)
	  && strcasecmp(ARR2STRING(qbuild->entries, i), "")) {
	
	if (strcmp(ARR2STRING(qbuild->entries, i - 2), "ANY TAG")
	    && strcmp(ARR2STRING(qbuild->entries, i - 2), "ITEM NAME")
	    && (!strcasecmp(ARR2STRING(qbuild->entries, i - 1), "contains")
		|| !strcasecmp(ARR2STRING(qbuild->entries, i - 1), "is equal to")
		|| !strcasecmp(ARR2STRING(qbuild->entries, i - 1), "=")
		|| !strcasecmp(ARR2STRING(qbuild->entries, i - 1), "begins with")
		|| !strcasecmp(ARR2STRING(qbuild->entries, i - 1), "ends with")))
	  strcat(result, " = ");
	
	if (TESTFLAG(QUOTE_ON))
	  strcat(result, "\"");
	if (!strcasecmp(ARR2STRING(qbuild->entries, i - 1), "contains")
	    || !strcasecmp(ARR2STRING(qbuild->entries, i - 1), "ends with"))
	  strcat(result, "*");
      }
      
      if (i % ENTRY_NUM == ATTRIBUTE) {
	if (!strcasecmp(ARR2STRING(qbuild->entries, i+1), "does not exist")) {
	  if (strcasecmp(ARR2STRING(qbuild->syntax, i), ""))
	    strcat(result, " ! ") ;
	  else
	    messout("Warning (line %d):\n Bad Attribute name specified for condition \"does not exist\"", qbuild->entry_set+1);
	}
	else if (!strncmp(ARR2STRING(qbuild->entries, i+1), "COUNT ", 6)) {
	  if (strcasecmp(ARR2STRING(qbuild->syntax, i), ""))
	    strcat(result, " COUNT ") ;
	  else
	    messout("Warning (line %d):\n Bad Attribute name specified for codition \"COUNT ...\"", qbuild->entry_set+1);
	}
      }

      strcat(result, ARR2STRING(qbuild->syntax, i));
      if (i % ENTRY_NUM == VALUE 
	  && strcasecmp(ARR2STRING(qbuild->entries, i), "")) {
	if (!strcasecmp(ARR2STRING(qbuild->entries, i - 1), "contains")
	    || !strcasecmp(ARR2STRING(qbuild->entries, i - 1), "begins with"))
	  strcat(result, "*");
	if (TESTFLAG(QUOTE_ON))
	  strcat(result, "\"");
      }
    }
  }
  
  
  removeWhite(result);
  return result;
}


static void preformQuery(int index);


/*
 * on calling from menu/button, will clear the 
 * user-entry boxes' contents.
 */
static void clearQbuilder(int from) 
{
  int i;
  
  QBUILDGET("clearQbuilder");

 switch (from) {
 case ALL_STUFF:
   qbuild->entry_set = 0;
   qbuild->entry_index = 0;

   for (i = 0; i < ENTRY_NUM; i++) {
     strcpy(ARR2STRING(qbuild->entries, i), "");
     strcpy(ARR2STRING(qbuild->syntax, i), "");
   }

 case PRECLASS:   /* deleting preclass_entry or preclass_syntax not done */
   strcpy(qbuild->classtag_entry, "");
   strcpy(qbuild->classtag_syntax, "");

 case CLASSTAG:
   for (i = 0; i < ENTRY_NUM; i++) {
     strcpy(ARR2STRING(qbuild->entries, i), "");
     strcpy(ARR2STRING(qbuild->syntax, i), "");
   }

 case ATTRIBUTE:
 case CONDITION:
 case VALUE:
 case CONJUNCTION:
 default:
   if (from < ENTRY_NUM)
     i = CONJUNCTION;
   else
     i = from - (from % ENTRY_NUM) + CONJUNCTION ;

   /* remove extra entries here */
   arrayMax(qbuild->entries) = i + 1;  /* ENTRY_NUM = CONJUNCTION + 1 */
   arrayMax(qbuild->syntax) = i + 1;
  
   strcpy(ARR2STRING(qbuild->entries, i), "END"); 
 }

}


/* various handler routines for each query builder box */
static void preclass_handler(char *entered_txt) 
{
  int i;
  void *look;
  KEYSET tempSet ;
  BOOL flag = FALSE;
  
  QBUILDGET("preclass_handler");
  
  if (entered_txt == NULL)
    messcrash("Error in preclass_handler(): NULL pointer");
  else {
    removeWhite(qbuild->preclass_entry);
    
    if (!strcmp(qbuild->preclass_entry, "ALL DATA")) {
      RMFLAG(RELATED);
      flag = TRUE;
      strcpy(qbuild->preclass_syntax, "Find ");
    }
    else if (!strcmp(qbuild->preclass_entry, "KEYSET")) {
      keySetActive(&tempSet, &look);
      if (tempSet) {
	/* reset Taglist to class in main window */
	qbuild_selected_class = class(keySet(tempSet, 0));
      }
      newKeyList(qbuild_selected_class, 
		 &classtagkeyset);  /* use previous class */
      SETFLAG(RELATED);
      RMFLAG(ALLTEXT);
      flag = TRUE;
      strcpy(qbuild->preclass_syntax, "Follow ");
    }
    else {
      for (i = 0; i < qbuildClassMax(); i++) {
	if (!strcasecmp(qbuild->preclass_entry, ARR2STRING(classlist, i)))  {
	  /* update class to user-entry */
	  qbuild_selected_class = pickWord2Class(qbuild->preclass_entry);
	  newKeyList(qbuild_selected_class, &classtagkeyset);
	  SETFLAG(RELATED);
	  RMFLAG(ALLTEXT);
	  flag = TRUE;
	  strcpy(qbuild->preclass_syntax, "Find ");
	  strcat(qbuild->preclass_syntax, qbuild->preclass_entry);
	  strcat(qbuild->preclass_syntax, " ; Follow ");
	  break;
	}
      }
    }
    
    if (!flag)  {
      strcpy(qbuild->preclass_syntax, "");
      messout("Invalid Type of Data (Preclass Box): Try Again");
    }
    else  {
      clearQbuilder(PRECLASS);
      preformQuery(PRECLASS); 
      qbuildPick(qbuild->class_box + 1); 
    /* classtag box is first in qbuild->class_box */
    }
  }
}



/* 
 * on carriage return, changes active box to attribute entry
 * and sets the entry field for the classtag
 *
 * special cases:  "KEYSET's" uses class From main window
 *                                  or previous query result
 *
 * Also, sets up qbuild_selected_class value and re-initializes Attrlist.
 * This is only called for the very first set of entry-boxes.
 */
static void classtag_handler(char *entered_txt) 
{
  int i;
  void *look;
  KEYSET tempSet ;
  BOOL no_return = TRUE;
  QBUILDGET("classtag_handler");
  
  if (entered_txt == NULL)
    messcrash("Error in classtag_handler(): NULL pointer");
  else if (!strcmp(qbuild->preclass_entry, "ALL DATA")) { /* Class name */
    if (!strcmp(qbuild->classtag_entry, "ANY CLASS")) {
      strcpy(qbuild->classtag_syntax, "Text");
      SETFLAG(ALLTEXT);
      no_return = FALSE;
    }
    else if (!strcmp(qbuild->classtag_entry, "KEYSET's")) {
      strcpy(qbuild->classtag_syntax,"");
      
      keySetActive(&tempSet, &look);
      if (tempSet)
	/* reset Attrlist to class in main window */
	selected_classtag = class(keySet(tempSet, 0));
      newKeyList(selected_classtag, &attrkeyset);   /* use previous class */
      RMFLAG(ALLTEXT);
      no_return = FALSE;
    }
    else {
      for (i = 0; i < qbuildClassMax(); i++) {
	if (!strcasecmp(qbuild->classtag_entry, ARR2STRING(classlist, i))) {
	  strcpy(qbuild->classtag_syntax, qbuild->classtag_entry); 
	  
	  /* update class to user-entry */
	  selected_classtag = pickWord2Class(qbuild->classtag_entry);
	  newKeyList(selected_classtag, &attrkeyset);
	  RMFLAG(ALLTEXT);
	  no_return = FALSE;
	  break;
	}
      } 
    }
  }
  else {  /* class name or KeySet for preclass pointer ">" operation */
    for (i = 0; i < keyListMax(classtagkeyset); i++) {
      if (!strcasecmp(qbuild->classtag_entry, 
		      KEYLIST_NAME(i, classtagkeyset))) {
	strcpy(qbuild->classtag_syntax, qbuild->classtag_entry);
	
	/* update class for attribute user entry */
	selected_classtag = pickWord2Class(qbuild->classtag_entry);
	newKeyList(selected_classtag, &attrkeyset);
	RMFLAG(ALLTEXT);
	no_return = FALSE;
	break;
      }
    }
  }
  
  if (no_return) {
    strcpy(qbuild->classtag_syntax,"");
    messout("Invalid Entry for Class/Tag (Classtag Box): Try Again.");
  }
  else {
    if (TESTFLAG(ALLTEXT)) {
      clearQbuilder(CLASSTAG);
      strcpy(qbuild->classtag_entry, "ANY CLASS");
      strcpy(qbuild->classtag_syntax,"Text");
      strcpy(ARR2STRING(qbuild->entries, ATTRIBUTE), "ANY TAG");
      strcpy(ARR2STRING(qbuild->entries, CONDITION), "contains");
      strcpy(ARR2STRING(qbuild->entries, CONJUNCTION), "END");
      
      /* preformQuery() calls qbuildDisplay(); */
      preformQuery(CLASSTAG);
      qbuildPick(INDEX2BOX( qbuild->user_box, VALUE ));
    }
    else {
      clearQbuilder(CLASSTAG);
      /* preformQuery() calls qbuildDisplay(); */
      preformQuery(CLASSTAG);
      qbuildPick(INDEX2BOX( qbuild->user_box, ATTRIBUTE ));
    }
    return;   
  }
}

/* see attrMenuSet, attrMenuResult */
static FREEOPT *attrMenu;  

/* 
 * on carriage return, changes active box to condition entry
 * and sets the syntax field for the attribute (tag)
 * MUST BE SET.  This entry-box MUST have a valid value.
 *   > Checks with tags in current class for validity. <
 */
static void attribute_handler(char *entered_txt)
{
  int i,  special_case = 0;
  BOOL found = FALSE;
  BOOL text_conj_bad = FALSE, text_other_bad = FALSE;
  char *entry_a=0, *syntax_a=0;
  
  QBUILDGET("attribute_handler");
  
  if (entered_txt == NULL)
    messcrash("Error in attribute_handler(): NULL pointer");
  else {
    entry_a = ARR2STRING(qbuild->entries, 
			 ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);
    syntax_a = ARR2STRING(qbuild->syntax, 
			  ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);
    removeWhite(entry_a);
    
    /* search to see if tagname exists */
    
    if (!strcmp(entry_a, "")) {
      special_case = 1;
      found = TRUE;
    }
    else if (!strcmp(entry_a, "ANY TAG")) {
      if (qbuild->entry_set != 0) 
	text_conj_bad = TRUE;
      else if (strcmp(qbuild->classtag_entry, "ANY CLASS"))
	text_other_bad = TRUE;
      else
	found = TRUE;
      special_case = 2;
    }
    else if (!strcmp(entry_a, "ITEM NAME")) { 
      special_case = 1;
      found = TRUE;
    }
    else {
      for (i = 1; i <= attrMenu[0].key ; i++) {
	if (!strcasecmp(entry_a, attrMenu[i].text)) {
	  found = TRUE;
	  break;
	}
      }
    }
  }
  
  
  if (found) {
    if (special_case) {
      if (strcmp(ARR2STRING(qbuild->syntax,
			    CONDITION + qbuild->entry_set*ENTRY_NUM),"")) 
	messout("Warning (line %d):\n (ITEM NAME or ANY TAG) cannot use current Operator.", qbuild->entry_set+1);

      strcpy(syntax_a, "");

      /* preformQuery() calls qbuildDisplay(); */
      preformQuery(ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);
      if (special_case == 2)
	qbuildPick(INDEX2BOX( qbuild->user_box, 
			     VALUE + qbuild->entry_set*ENTRY_NUM));
      else
	qbuildPick(INDEX2BOX( qbuild->user_box, 
			     CONDITION + qbuild->entry_set*ENTRY_NUM));
    }
    else {
      strcpy(syntax_a, entry_a);
      
      /* calls qbuildDisplay */
      preformQuery(ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);

      qbuildPick(INDEX2BOX( qbuild->user_box, 
			   CONDITION + qbuild->entry_set*ENTRY_NUM));
    }
  }  
  else if (text_conj_bad){
    graphTextScrollEntry(entry_a, 0, 0, 0, 0, 0);
    messout("Warning (line %d):\n No conjunctive queries allowed in a text search.", qbuild->entry_set+1);
  }
  else if (text_other_bad) {
    graphTextScrollEntry(entry_a, 0, 0, 0, 0, 0);
    messout("Warning (line %d):\n Bad syntax for query. Check Attribute.", 
	    qbuild->entry_set+1);
  }
  else {                       /* do not set str to null because */
    /* may correct by changing classname */
    strcpy(syntax_a, entry_a);
    graphTextScrollEntry(entry_a, 0, 0, 0, 0, 0);
    messout("Invalid Tag/Attribute Name (line %d): Try Again", 
	    qbuild->entry_set+1);
  }
}                                  /* end fxn */


/* 
 * on carriage return, changes active box to value entry
 * and sets the entry field for the condition
 *
 * checks current entered_txt for validity as condition
 *  returns TRUE if valid
 * Rough translations of valid conditions:
 *	is equal to 			= or NULL
 *	is not equal to 		!=
 *	is greater than 		>
 *	is greater than or equal to 	>=
 *	is less than 			<
 *	is less than or equal to	<=
 *	exists                          NULL string
 *      <other 'operators'>             NULL string
 */

static void condition_handler(char *entered_txt)
{
  BOOL okay_flag = TRUE;
  char *entry_c, *syntax_c, *entry_a;
  QBUILDGET("condition_handler");
  
  entry_a = ARR2STRING(qbuild->entries, 
		       ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);
  entry_c = ARR2STRING(qbuild->entries, 
		       CONDITION + qbuild->entry_set*ENTRY_NUM);
  if (entered_txt == NULL)
    messcrash("Error in condition_handler(): NULL pointer");
  else {
    removeWhite(entry_c);
    syntax_c = ARR2STRING(qbuild->syntax, 
			  CONDITION + qbuild->entry_set*ENTRY_NUM);
    
    if (TESTFLAG(ALLTEXT) || !strcasecmp(entry_a, "ITEM NAME")) {
      if ( strcasecmp(entry_c, "exists")
	  && strcasecmp(entry_c, "contains")
	  && strcasecmp(entry_c, "is equal to")
	  && strcmp(entry_c, "=")
	  && strcasecmp(entry_c, "ends with")
	  && strcasecmp(entry_c, "begins with")) {
	messout("Invalid condition (line %d): Only five operators (exists, is equal to, contains, begins with, ends with) allowed for Attributes \"ANY TAG\" and \"ITEM NAME\": Try Again", qbuild->entry_set+1);
	strcpy(syntax_c,"");
	graphTextScrollEntry(entry_c, 0, 0, 0, 0, 0);
	okay_flag = FALSE;
      }
      strcpy(syntax_c, "");
    }
    else {
      if (!strcmp(entry_c, "!=")
	  || !strcmp(entry_c, ">")
	  || !strcmp(entry_c, ">=")
	  || !strcmp(entry_c, "<")
	  || !strcmp(entry_c, "<="))
	strcpy(syntax_c, entry_c);
      else if (   !strcasecmp(entry_c, "exists")
	       || !strcasecmp(entry_c, "does not exist")
	       || !strcasecmp(entry_c, "is equal to")
	       || !strcmp(entry_c, "=")
  /*    || !strcasecmp(entry_c, "COUNT is equal to") mhmp 05.06.98 */
	       || !strcasecmp(entry_c, "ends with")
	       || !strcasecmp(entry_c, "begins with")
	       || !strcasecmp(entry_c, "contains"))
	strcpy(syntax_c, "");
      else if (!strcasecmp(entry_c, "COUNT is equal to"))
	strcpy(syntax_c, "=");	       
      else if (!strcasecmp(entry_c, "is not equal to")
	       || !strcasecmp(entry_c, "COUNT is not equal to"))
	strcpy(syntax_c, "!=");
      else if (!strcasecmp(entry_c, "is greater than")
	       || !strcasecmp(entry_c, "COUNT is greater than"))
	strcpy(syntax_c, ">");
      else if (!strcasecmp(entry_c, "is greater than or equal to")
	       || !strcasecmp(entry_c, "COUNT is greater than or equal to"))
	strcpy(syntax_c, ">=");
      else if (!strcasecmp(entry_c, "is less than") 
	       || !strcasecmp(entry_c, "COUNT is less than"))
	strcpy(syntax_c, "<");
      else if (!strcasecmp(entry_c, "is less than or equal to")
	       || !strcasecmp(entry_c, "COUNT is less than or equal to"))
	strcpy(syntax_c, "<=");
      else if (!strcmp(entry_c, "")) {
	strcpy(syntax_c, "");
	messout("Warning (line %d): No Conditional operator specified", 
		qbuild->entry_set+1);
      }
      else {
	strcpy(syntax_c, "");
	messout("Invalid entry for Conditional operator (line %d): Try Again", 
		qbuild->entry_set+1);
	okay_flag = FALSE;
      }
    }   /* end of else (not ALLTEXT) */

    if (okay_flag) {
      if (!strcmp(ARR2STRING(qbuild->syntax, 
			     ATTRIBUTE + qbuild->entry_set*ENTRY_NUM), "")
	  &&
	  !strcasecmp(syntax_c, "contains"))
	messout("Warning (line %d):\n No Attribute name", qbuild->entry_set+1);
      
      preformQuery(CONDITION + qbuild->entry_set*ENTRY_NUM); 
      
      if (!strcasecmp(entry_c, "exists") ||
	  !strcasecmp(entry_c, "does not exist"))
	qbuildPick(INDEX2BOX( qbuild->user_box, 
			     CONJUNCTION + qbuild->entry_set*ENTRY_NUM ));
      else
	qbuildPick(INDEX2BOX( qbuild->user_box, 
			     VALUE + qbuild->entry_set*ENTRY_NUM ));
    }
  }
}


static void value_handler(char *entered_txt)
{
  char *entry_a, *syntax_a, *entry_c, *syntax_c;
  char *entry_v, *syntax_v, *quote_chr, *qc;
  QBUILDGET("value_handler");
  
  if (entered_txt == NULL)
    messcrash("Error in value_handler(): NULL pointer");
  else {
    entry_a = ARR2STRING(qbuild->entries, 
			 ATTRIBUTE + qbuild->entry_set*ENTRY_NUM); 
    syntax_a = ARR2STRING(qbuild->syntax, 
			  ATTRIBUTE + qbuild->entry_set*ENTRY_NUM); 
    entry_c = ARR2STRING(qbuild->entries, 
			 CONDITION + qbuild->entry_set*ENTRY_NUM); 
    syntax_c = ARR2STRING(qbuild->syntax, 
			  CONDITION + qbuild->entry_set*ENTRY_NUM); 
    entry_v = ARR2STRING(qbuild->entries, 
			 VALUE + qbuild->entry_set*ENTRY_NUM); 
    syntax_v = ARR2STRING(qbuild->syntax, 
			  VALUE + qbuild->entry_set*ENTRY_NUM); 
    
    removeWhite(entry_v);
    
    if (!strcmp(syntax_a, "")   /* possible bad value */
	&&
	(strcmp(entry_a, "ANY TAG")  /* valid, null-syntax ATTRs */
	 && strcmp(entry_a, "ITEM NAME"))
	&& (strcasecmp(entry_c, "contains")  /* valid cond w/ no ATTR */
	    && strcasecmp(entry_c, "begins with")
	    && strcasecmp(entry_c, "exists")
	    && strcasecmp(entry_c, "is equal to")
	    && strcasecmp(entry_c, "ends with")))
      messout("Warning (line %d):\n No Attribute name", qbuild->entry_set+1);
    else if ((strcmp(syntax_c,"")
	     ||
	      (strcasecmp(entry_c, "exists")
	       && strcasecmp(entry_c, "does not exist")
	       && strcasecmp(entry_c, "")))
	     && !strcmp(entry_v, ""))   /* no VAL for COND that requires VAL*/
      messout("Warning (line %d):\n No value (right-hand argument) for Conditional operator which requires one.",
	      qbuild->entry_set+1);
    else if (((!strcasecmp(entry_c, "does not exist") ||
	      !strcasecmp(entry_c, "exists"))
	      &&
	      !strcmp(syntax_c, ""))
	     &&
	     strcmp(entry_v,""))
      messout("Warning (line %d):\n Existence operators do not use values (right-hand argument)", 
	      qbuild->entry_set+1);
    
    strcpy(syntax_v, entry_v);
    if (isNumeric(syntax_v))
      RMFLAG(QUOTE_ON);
    else {
      SETFLAG(QUOTE_ON);
      quote_chr = strchr(syntax_v, '"');
      while (quote_chr && *quote_chr) {
	if (*(quote_chr - 1) != (char) 92) {
	  qc = quote_chr+strlen(quote_chr);
	  *(qc+1) = 0 ;
	  while (qc != quote_chr)
	    { *(qc) = *(qc-1) ;
	      --qc ;
	    }
	  *quote_chr = (char) 92;  /* add a \ */
	  quote_chr += 3 ;  /* skip over the 2 characters: \" */
	}
	else {
	  quote_chr += 1;  /* skip over the 1 character: " */
	}
	if (quote_chr && *quote_chr)
	  quote_chr = strchr(quote_chr, '"');
      }
    }
    
    preformQuery(VALUE + qbuild->entry_set*ENTRY_NUM);  /* in value handler */
    qbuildPick(INDEX2BOX( qbuild->user_box, 
			 CONJUNCTION + qbuild->entry_set*ENTRY_NUM ));
  }
}    


/* check to see if hash exists in model for current Tag label */
static BOOL checkTagForHash(int label_box, int *classe) {
  KEY t, new;
  OBJ obj;
  QBUILDGET ("checkTagForHash");
  if (lexword2key(ARR2STRING(qbuild->entries, label_box), &t, 0)) {
    obj = bsCreate(KEYMAKE(selected_classtag, 0));   /* 0 for system-model */
    bsGetKeyTags(obj, t, &t);
    while (t) {
      new = bsType (obj,_bsRight) ;  /* move to the right : tag ?value #hash */
      if (t == new)
	break;
      else
	t = new;
      *classe = class(t) ;
      if (*classe && t == KEYMAKE(*classe, 1) && 
	  pickType(KEYMAKE(*classe, 1)) == 'B') { 
	bsDestroy(obj);
	return TRUE;
      }
    }
    bsDestroy(obj);
  }
  return FALSE;
}


/* Checks current "entered_txt" for validity as conjunction
 *  returns TRUE if one of valid conjunctions (and &, or |, xor ^)
 */
static void conjunction_handler(char *entered_txt)
{
  BOOL okay_entry = TRUE, terminal_entry = FALSE;
  int i, max, classe;
  char *entry_j, *syntax_j, *entry_a, *syntax_a;
  
  QBUILDGET("conjunction_handler");
  
  entry_a = ARR2STRING(qbuild->entries, ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);
  syntax_a = ARR2STRING(qbuild->syntax, ATTRIBUTE + qbuild->entry_set*ENTRY_NUM);
  entry_j = ARR2STRING(qbuild->entries, CONJUNCTION + qbuild->entry_set*ENTRY_NUM);
  syntax_j = ARR2STRING(qbuild->syntax, CONJUNCTION + qbuild->entry_set*ENTRY_NUM);
  
  if (entered_txt == NULL)
    messcrash("Error in conjunction_handler(): NULL pointer");
  else if (!strcmp(syntax_a, "") && !strcmp(entry_a, "")) {
    messout("Conjunction must be preceded by an attribute and a condition. (line %d)", qbuild->entry_set+1);
    strcpy(syntax_j, "");
    strcpy(entry_j, "END");
    preformQuery(CONJUNCTION + qbuild->entry_set*ENTRY_NUM);  
    qbuild->entry_index = ATTRIBUTE;
  }
  else {
    removeWhite(entry_j);

    if (!strcmp(entry_j, "#") || !strcasecmp(entry_j, "and subtype")) {
      if (checkTagForHash(ATTRIBUTE + qbuild->entry_set, &classe)) {
	newKeyList(classe, &attrkeyset);
	strcpy(syntax_j, "#");
	attr_sub_model = TRUE;
      }
      else {
	strcpy(entry_j, "and");
	strcpy(syntax_j, "&");
	newKeyList(selected_classtag, &attrkeyset);
	attr_sub_model = FALSE;
	messout(messprintf("Sorry, but no #subtype for the Attribute %s exists in the model %s", 
			   ARR2STRING(qbuild->entries, 
				      ATTRIBUTE + qbuild->entry_set),
			   qbuild->classtag_entry));
      }
    }
    else  {
      if (attr_sub_model) {
	newKeyList(selected_classtag, &attrkeyset);
	attr_sub_model = FALSE;
      }

      if (!strcmp(entry_j, "&")
	  || !strcmp(entry_j, "|")
	  || !strcmp(entry_j, "^"))
	strcpy(syntax_j, entry_j);
      else if (!strcasecmp(entry_j, "and"))
	strcpy(syntax_j, "&");
      else if (!strcasecmp(entry_j, "or"))
	strcpy(syntax_j, "|");
      else if (!strcasecmp(entry_j, "xor")
	       || !strcasecmp(entry_j, "exclusive or"))
	strcpy(syntax_j, "^");
      else if (!strcasecmp(entry_j, "END")) {
	strcpy(syntax_j, "");
	terminal_entry = TRUE;
      }
      else {
	strcpy(syntax_j, "");
	okay_entry = FALSE;
	messout("Invalid Conjunction Entry (line %d): Try Again",
		qbuild->entry_set+1);
      }
    }

    if (okay_entry) {
      /* where to pick to (next box)... */
      if (terminal_entry) {
	qbuild->entry_index = ATTRIBUTE;
	arrayMax(qbuild->syntax) = qbuild->entry_set*ENTRY_NUM + CONJUNCTION + 1;
	arrayMax(qbuild->entries) = qbuild->entry_set*ENTRY_NUM + CONJUNCTION +1;
	preformQuery(CONJUNCTION + qbuild->entry_set*ENTRY_NUM);  
	qbuildMenuExec();
      }
      else if (!TESTFLAG(ALLTEXT)
	       && strcasecmp(entry_j, "END")) {
	/* expand number of user-entry boxes */
	max = arrayMax(qbuild->syntax);
	
	/* the following "if" probably should be 
	 * qbuild->entry_index == CONJ + offset 
	 * rather than qbuild->entry_set ==  ...
	 */
	if (qbuild->entry_set == max/ENTRY_NUM - 1)  {
	  /* add more user boxes */
	  for (i = max; i < max + ENTRY_NUM; i++)  {
	    array(qbuild->syntax, i, Array) = arrayCreate(BUFFER_SIZE, char);
	    strcpy(ARR2STRING(qbuild->syntax, i), "");
	    array(qbuild->entries, i, Array) = arrayCreate(BUFFER_SIZE, char);
	    strcpy(ARR2STRING(qbuild->entries, i), "");
	    if (i % ENTRY_NUM == CONJUNCTION)
	      strcpy(ARR2STRING(qbuild->entries, i), "END");
	  }
	}
	preformQuery(CONJUNCTION + qbuild->entry_set*ENTRY_NUM);  
	if (qbuild->entry_index < arrayMax(qbuild->entries))
	  qbuild->entry_index = qbuild->entry_index + 1; 
	qbuild->entry_set = qbuild->entry_index/ENTRY_NUM;
	qbuildPick(INDEX2BOX( qbuild->user_box, 
			     ATTRIBUTE + qbuild->entry_set*ENTRY_NUM ));
      }
      else {
	messout("Warning (line %d):\n No conjunctive queries allowed in a text search.", qbuild->entry_set+1);
	strcpy(syntax_j, "");
	strcpy(entry_j, "END");
	preformQuery(CONJUNCTION + qbuild->entry_set*ENTRY_NUM);  
	qbuild->entry_index = ATTRIBUTE;
      }
    }
  }
}


static void preformQuery(int index)
{
  int i;
  char *i_entry = NULL;
  
  QBUILDGET("preformQuery");
  
  messStatus("Forming Query Command");
  i = index % ENTRY_NUM;
  if (i >= 0)
    i_entry = ARR2STRING(qbuild->entries, index);

  /* most syntax checking occurs in *_handler() routines, but here is
   * a few last checks */
  switch (i) {
    /* know that Order of TYPES is ATTR, COND, VALUE, CONJ */
  case ATTRIBUTE:
    if ((index >= ENTRY_NUM)
	&& !strcasecmp(ARR2STRING(qbuild->entries, index - 1), "END"))
    /* should not reach this point   gha 6/16/93 */
    messout("Warning (line %d):\n Entries after (END) are ignored and will be removed.", 
	      index/ENTRY_NUM + 1);
    break;
  case CONJUNCTION:
    /* clear all entries after "END" */
    if (!strcasecmp(i_entry, "END")
	&& (index > ENTRY_NUM))  {
      arrayMax(qbuild->entries) = index+1;
      arrayMax(qbuild->syntax) = index+1;
    }
    break;

  default:
    break;
  }
  
  if (querGraph != 0 && graphExists(querGraph)) {
    if(!graphActivate(querGraph))
      messcrash("Cannot reactivate querGraph in querybuild.c: preformQuery.\n");
    messStatus("Forming Query Command");
    queryWindowShowQuery(makeQuery(resbuffer));

    if(!graphActivate(qbuildGraph))
      messcrash("Cannot reactivate qbuildGraph in querybuild.c: preformQuery.\n");
  }

  qbuildDisplay();
}

/* Check for mistakes in the user-entry boxes.
 * iterates through each entry handler...
 * not optimally designed right now b/c each handler
 * needs to set up menus and other values, thus must redraw
 * which makes this entire process slow.
 */
static void checkQueryFormat(void)
{
  int offset, max;
  int old_entry_index, old_entry_set;

  QBUILDGET("checkQuery");
  
  old_entry_index = qbuild->entry_index;
  old_entry_set = qbuild->entry_set;
  max = arrayMax(qbuild->entries);
  
  for (offset = 0; offset < max; offset+=ENTRY_NUM) {
    qbuild->entry_set = offset/ENTRY_NUM;
    qbuild->entry_index = offset + ATTRIBUTE;
    attribute_handler(ARR2STRING(qbuild->entries, ATTRIBUTE + offset));
    condition_handler(ARR2STRING(qbuild->entries, CONDITION + offset));
    value_handler(ARR2STRING(qbuild->entries, VALUE + offset));
    if (!strcasecmp(ARR2STRING(qbuild->entries, CONJUNCTION + offset), "END")
	&&
	!strcmp(ARR2STRING(qbuild->syntax, CONJUNCTION + offset), ""))
      break;
    conjunction_handler(ARR2STRING(qbuild->entries, CONJUNCTION + offset));
  }
  qbuild->entry_index = old_entry_index;
  qbuild->entry_set = old_entry_set;
}
/************************************/

/* MENU STUFF */

typedef struct MenuListStruct
{ 
  FREEOPT *menup;
  int size;
  struct MenuListStruct *next;
} *ML;

/* preclassMenu gets set once at init.
 * attrMenu and classtagMenu may vary as user builds query
 * attrMenu declared above attribute_handler().
 */
static ML attrMenus, classtagMenus;
static FREEOPT *preclassMenu, *classtagMenu;

static FREEOPT conjMenu[] =
{
  { 5, "CONJUNCTIONS" },
  { 1, "END" },
  { 2, "and" },
  { 3, "and subtype" },
  { 4, "or" },
  { 5, "exclusive or" },
};


static FREEOPT valueMenu[] =
{
  { 3, "VALUE OPERATORS" },
  { 1, "*" },
  { 2, "?" },
  { 3, "\"" },
};

static FREEOPT condMenu[] = 
{
  { 19, "CONDITIONS" },
  { 1, "exists" },
  { 2, "does not exist" },
  { 3, "contains" },
  { 4, "begins with" },
  { 5, "ends with" },
  { 6, "" },
  { 7, "is equal to" },
  { 8, "is not equal to" },
  { 9, "is greater than" },
  { 10, "is greater than or equal to" },
  { 11, "is less than" },
  { 12, "is less than or equal to" },
  { 13, "" },
  { 14, "COUNT is equal to" },
  { 15, "COUNT is not equal to" },
  { 16, "COUNT is greater than" },
  { 17, "COUNT is greater than or equal to" },
  { 18, "COUNT is less than" },
  { 19, "COUNT is less than or equal to" },
};


void destroyMenu(FREEOPT *menup, int max);

static void destroyML(ML menulist)
{
  ML curr,next;
  
  curr = menulist;
  if (curr) {
    for (next = curr->next; next != NULL ; curr = next, next = next->next)  {
      destroyMenu(curr->menup, curr->size);
      messfree(curr);
      curr = NULL;
    }
    destroyMenu(curr->menup, curr->size);
    messfree(curr);
    curr = NULL;
  }
}

static int MENU_HANDLE ;

static void destroyPreclassMenu()
{ char *v0 = 0,  *v1 ;
  int max;
  if (graphAssFind(&MENU_HANDLE, &v1)) 
    { max = v1 - v0 ;
      destroyMenu(preclassMenu, max);
      graphAssRemove(&MENU_HANDLE);
    }
}

/*  *  *  *  *  *  *  *  *  *  */

static void menuResultHandler(KEY k, int type)
{
  char *m_entry;
  int restype, resclass;
  KEY reskey ;
  Stack resstack;
  QBUILDGET("menuResultHandler");

  m_entry = ARR2STRING(qbuild->entries, type+qbuild->entry_set*ENTRY_NUM);

  switch (type) {
  case ATTRIBUTE:
    if (attrMenu != arr(qbuild->attrMenuArray, qbuild->entry_set, FREEOPT *))
      attrMenu = arr(qbuild->attrMenuArray, qbuild->entry_set, FREEOPT *);
    if (!strcmp((attrMenu+k)->text, "USE TREE TAG CHOOSER")) {
      resstack = stackCreate(10);
      stackClear(resstack);
      if (pickList[selected_classtag].type == 'B') {
	if (treeChooseTagFromModel(&restype, &resclass, 
				   selected_classtag, 
				   &reskey, resstack, 0)) {
	  strcpy(m_entry, stackText(resstack, 0));
	  attribute_handler(m_entry);
	}
      }
      else
	messout("Sorry, but the class-tag has no model for the Tag Chooser to display.");
      stackDestroy(resstack);
    }
    else {
      strcpy(m_entry, (attrMenu + k)->text);
      attribute_handler(m_entry);
    }
    break;
  case CONDITION:
    strcpy(m_entry, (condMenu + k)->text);
    condition_handler(m_entry);
    break;
  case VALUE:
    strcpy(m_entry, (valueMenu + k)->text);
    value_handler(m_entry);
    break;
  case CONJUNCTION:
    strcpy(m_entry, (conjMenu + k)->text);
    conjunction_handler(m_entry);
    break;
  default:
    messerror("Bad type in menuResultHandler");
  }
}  

static void attrMenuResult(KEY k) 
{ menuResultHandler(k, ATTRIBUTE); }

static void condMenuResult(KEY k)
{ menuResultHandler(k, CONDITION); }

static void valueMenuResult(KEY k)
{ menuResultHandler(k, VALUE); }

static void conjMenuResult(KEY k) 
{ menuResultHandler(k, CONJUNCTION); }


static void classtagMenuResult(KEY k, int box)
{
  int restype, resclass;
  KEY reskey ;
  Stack resstack;
  QBUILDGET("classtagMenuResult");

  if (!strcmp((classtagMenu+k)->text, "USE TREE TAG CHOOSER")) {
    resstack = stackCreate(10);
    stackClear(resstack);
    if (pickList[qbuild_selected_class].type == 'B') {
      if (treeChooseTagFromModel(&restype, &resclass, 
				 qbuild_selected_class, 
				 &reskey, resstack, FALSE)) {
	strcpy(qbuild->classtag_entry, stackText(resstack, 0));
	classtag_handler(qbuild->classtag_entry);
      }
    }
    else
      messout("Sorry, but the class has no model for the Tag Chooser to display.");
    stackDestroy(resstack);
  }
  else {
    strcpy(qbuild->classtag_entry, (classtagMenu + k)->text);
    classtag_handler(qbuild->classtag_entry);
  }
}


static void preclassMenuResult(KEY k, int box) {
  QBUILDGET("preclassMenuResult");
  
  strcpy(qbuild->preclass_entry, (preclassMenu + k)->text);
  preclass_handler(qbuild->preclass_entry);
}     

/*  *  *  *  *  *  *  *  *  *  */

static void preclassMenuSet() 
{ char *v0 = 0,  *v1 ;
  int i, max;
  
  QBUILDGET("preclassMenuSet");
  
  max = qbuildClassMax() + 2;
  
  /* max should be one larger for header-entry */
  preclassMenu = (FREEOPT *)messalloc((max + 1) *sizeof(FREEOPT));
  preclassMenu->key = max;
  /* suz, changed from arrayCreate to messalloc */
  preclassMenu->text = (char *) messalloc(BUFFER_SIZE);
  strcpy(preclassMenu->text, "DATA TYPE");
  
  for (i = 1; i <= max; i++)  {
    (preclassMenu + i)->key = i;
    /* suz, changed from arrayCreate to messalloc */
    (preclassMenu + i)->text = (char *) messalloc (BUFFER_SIZE);
    if (i == 1)
      strcpy((preclassMenu + i)->text, "ALL DATA");
    else if (i == 2)
      strcpy((preclassMenu + i)->text, "KEYSET");
    else                       /* change this (- 3) if (+ 2) in max changes */
      strcpy((preclassMenu + i)->text, ARR2STRING(classlist, i - 3));
  }
  
  v1 = v0 + max ;
  graphAssociate(&MENU_HANDLE, v1);
}

static void classtagMenuSet(int box)
{
  int i, max;
  FREEOPT *newMenu;
  ML newML, curr, prev;
  
  QBUILDGET("classtagMenuSet");
  
  if (!box)
    return;
  
  if (strcmp(qbuild->preclass_entry, "ALL DATA"))  {  /* Tag */
    max = keyListMax(classtagkeyset) + 1;
    
    /* max should be one larger for header-entry */
    newMenu = (FREEOPT *)messalloc((max + 1) *sizeof(FREEOPT));
    newMenu->key = max;
    /* suz, changed from arrayCreate to messalloc */
    newMenu->text = (char *) messalloc(BUFFER_SIZE);

    strcpy(newMenu->text, "TAGS IN CLASS");
    
    for (i = 1; i <= max; i++) {   
      (newMenu + i)->key = i;
      /* suz, changed from arrayCreate to messalloc */
      (newMenu + i)->text = (char *) messalloc(BUFFER_SIZE);
      if (i == 1)
	strcpy((newMenu + i)->text, "USE TREE TAG CHOOSER");
      else
	strcpy((newMenu + i)->text, KEYLIST_NAME(i - 2, classtagkeyset));
    }
  }
  else  {                                /* Classes */
    max = qbuildClassMax() + 2;
    
    /* max should be one larger for header-entry */
    newMenu = (FREEOPT *)messalloc((max + 1) *sizeof(FREEOPT));
    newMenu->key = max;
    /* suz, changed from arrayCreate to messalloc */
    newMenu->text = (char *) messalloc (BUFFER_SIZE);
    
    strcpy(newMenu->text, "CLASSES");
    
    for (i = 1; i <= max; i++)  {
      (newMenu + i)->key = i;
      /* suz, changed from arrayCreate to messalloc */
      (newMenu + i)->text = (char *) messalloc(BUFFER_SIZE);
      if (i == 1)
	strcpy((newMenu + i)->text, "ANY CLASS");
      else if (i == 2)
	strcpy((newMenu + i)->text, "KEYSET's");
      else 
	strcpy((newMenu + i)->text, ARR2STRING(classlist, i - 3));
    }
  }
  
  prev = classtagMenus;
  
  if (prev)
    curr = prev->next;
  else
    curr = NULL;
  for (; curr ; prev = curr, curr = curr->next)  {
    if (curr->size == max)  {
      /* check for match...i > max means match */
      for (i = 1; (i <= max) 
	   && !strcmp( (newMenu +i)->text, ((curr->menup)+i)->text);
	   i++);
      if (i > max)
	break;   /* only breaks on match */
    }
  }
  
  if (curr)  {  /* match found */
    prev->next = curr->next;   /* keep sorted in Most Recently Used */
    curr->next = classtagMenus;
    classtagMenus = curr;
    destroyMenu(newMenu, max);
  }
  else {
    newML = (ML) messalloc(sizeof(struct MenuListStruct));
    newML->menup = newMenu;
    newML->size = max;
    newML->next = classtagMenus;
    classtagMenus = newML;        /* current menu is always at front */
  }
  
  graphBoxFreeMenu(graphTextScrollEntry(qbuild->classtag_entry, 0, 0, 0, 0, 0),
		   (FreeMenuFunction) classtagMenuResult, 
		   classtagMenus->menup);

  classtagMenu = classtagMenus->menup;
}


static void attrMenuSet(int box)
{
  int i, max;
  BOOL istext = FALSE;
  FREEOPT *newMenu;
  ML newML, curr, prev;
  QBUILDGET("attrMenuSet");
  
  if (!box)
    return;

  if (!strcmp(qbuild->classtag_entry,"ANY CLASS"))  {
    istext = TRUE;
    max = 1;
  }
  else 
    max = keyListMax(attrkeyset) + 2;
  
  /* add 1 to #tags because of header-entry for FREEOPT */
  newMenu = (FREEOPT *)messalloc((1 + max)*sizeof(FREEOPT));
  newMenu->key = max;
  /* suz, changed from arrayCreate to messalloc */
  newMenu->text = (char *) messalloc(BUFFER_SIZE);
  strcpy(newMenu->text, "ATTRIBUTES"); 
  
  if (istext) {
    (newMenu + 1)->key = 1;
    /* suz, changed from arrayCreate to messalloc */
    (newMenu + 1)->text = (char *) messalloc(BUFFER_SIZE);
    strcpy((newMenu + 1)->text, "ANY TAG");
  }
  else 
    for (i = 1; i  <= max; i++)  {
      (newMenu + i)->key = i;
      /* suz, changed from arrayCreate to messalloc */
      (newMenu + i)->text = (char *) messalloc(BUFFER_SIZE);
      if (i == 1)
	strcpy((newMenu + i)->text, "USE TREE TAG CHOOSER");
      else if (i == 2)
	strcpy((newMenu + i)->text, "ITEM NAME");
      else
	strcpy((newMenu + i)->text, KEYLIST_NAME(i - 3, attrkeyset));
    }
  
  prev = attrMenus;
  
  if (prev)
    curr = prev->next;
  else
    curr = NULL;
  for (; curr ; prev = curr, curr = curr->next) {
    if (curr->size == max)  {
      /* check for match...i > max means match */
      for (i = 1; (i <= max) 
	   && !strcmp( (newMenu +i)->text, ((curr->menup)+i)->text);
	   i++);
      if (i > max)
	break;   /* only breaks on match */
    }
  }
  
  if (curr)  {  /* match found */
    prev->next = curr->next;   /* keep sorted in Most Recently Used */
    curr->next = attrMenus;
    attrMenus = curr;
    destroyMenu(newMenu, max);
  }
  else  {
    newML = (ML) messalloc(sizeof(struct MenuListStruct));
    newML->menup = newMenu;
    newML->size = max;
    newML->next = attrMenus;
    attrMenus = newML;        /* current menu is always at front */
  }
  
  i = BOX2INDEX(qbuild->user_box, box);
  
  graphBoxFreeMenu(graphTextScrollEntry(ARR2STRING(qbuild->entries, i), 
					0, 0, 0, 0, 0),
		   (FreeMenuFunction) attrMenuResult, attrMenus->menup);
  attrMenu = attrMenus->menup;
  array(qbuild->attrMenuArray, i/ENTRY_NUM, FREEOPT *) = attrMenu;
}


/* Called by QUIT */
static void qbuildDestroy (void)
{    
  int i ;
  QBUILD qbuild;
  
  if ( graphAssFind (&MAGIC, &qbuild)
      && (qbuild->magic == MAGIC))  {

    qbuild->magic = 0 ;

    for (i = 0; i < arrayMax(qbuild->entries); i++)  {
      arrayDestroy(arr(qbuild->entries, i, Array));
      arrayDestroy(arr(qbuild->syntax, i, Array)); 
    }
    arrayDestroy(qbuild->entries);
    arrayDestroy(qbuild->syntax);
    
    if (!graphExists(querGraph) && !graphExists(qbeGraph))
      arrayDestroy(query_undo_set);

    for (i = 0; i < qbuildClassMax(); i++)
      arrayDestroy(arr(classlist, i, Array));
    arrayDestroy(classlist);
    
    keySetDestroy(classtagkeyset);
    
    destroyML(attrMenus);
    arrayDestroy(qbuild->attrMenuArray);
    destroyML(classtagMenus);
    destroyPreclassMenu();

    messfree (qbuild) ;
    
    if(graphActivate(qbuildGraph))
      graphDestroy () ;
    qbuildGraph = 0 ;
  }
}

/**********************************************************/

/*
 * modified gha 9/4/92
 * Mouse pick function for query graph window...based upon 
 * format of boxes as called in queryCreate()...caveat programmer.
 */
static void qbuildPick (int box)
{    
  int i ;
  KEY key;
  QBUILDGET("qbuildPick");
  
  if (box <= 0)
    return ;

  if (box > qbuild->user_box)  {        /* for user entry boxes */
    i = BOX2INDEX(qbuild->user_box, box);
    qbuild->entry_set = i/ENTRY_NUM;
    if (isSelector(qbuild->user_box, box)) {
      messStatus("Selector");
      graphTextScrollEntry(ARR2STRING(qbuild->entries, i), 0, 0, 0, 0, 0);
      switch(i % ENTRY_NUM) {
      case ATTRIBUTE:
	if (graphSelect(&key, attrMenu))
	  attrMenuResult(key);
	break;
      case CONDITION:
	if (graphSelect(&key, condMenu))
	  condMenuResult(key);
	break;
      case VALUE:
	if (graphSelect(&key, valueMenu))
	  valueMenuResult(key);
	break;
      case CONJUNCTION:
	if (graphSelect(&key, conjMenu))
	  conjMenuResult(key);
	break;
      default:
	messerror("Out of bounds index in qbuildPick()");
      }
    }
    else if (i < arrayMax(qbuild->entries)) {  /* gha 3/8/93 */
      qbuild->entry_index = i ;
      graphTextScrollEntry(ARR2STRING(qbuild->entries, i), 0, 0, 0, 0, 0);
    }
  }
  else if (box > qbuild->class_box)  {       /* for class-/-tag boxes */

    if (box == qbuild->class_triBox) {
      graphTextScrollEntry(qbuild->classtag_entry, 0, 0, 0, 0, 0);
      messStatus("Selector");
      if (graphSelect(&key, classtagMenu)) {
	classtagMenuResult(key, 0);
      }
      i = 0;
    }
    else if (box == qbuild->pre_triBox) {
      messStatus("Selector");
      graphTextScrollEntry(qbuild->preclass_entry, 0, 0, 0, 0, 0);
      if (graphSelect(&key, preclassMenu)) {
	preclassMenuResult(key, 0);
      }
      i = 1;
    }
    else
      i = ((box - qbuild->class_box) < 4) ? 0 : 1;
      
    qbuild->entry_set = 0;
    qbuild->entry_index = 0;
    if (!i)
      { if (TESTFLAG(SIMPLE_UI))
	  graphTextScrollEntry(qbuild->classtag_entry, 0, 0, 0, 0, 0);
	else
	  graphTextScrollEntry(qbuild->classtag_entry, 0, 0, 0, 0, 0);
      }	  
    else if (!TESTFLAG(SIMPLE_UI) && i == 1)
      graphTextScrollEntry(qbuild->preclass_entry, 0, 0, 0, 0, 0);
  }
  else
    return;

  graphActivate(qbuild->graph) ;
}

/****************************************************************/

static void qbeCallup(void) { 
  qbeCreate();
  qbuildDisplay();  /* gha added to unhighlight builder button 4/29/93 */
}

static void qcommandCallup(void) { 
  queryCreate();
  qbuildDisplay();  /* gha added to unhighlight builder button 4/29/93 */
}

static void qbuildButtonUndo(void) {
  queryCommandUndo();
  qbuildDisplay();  /* gha added to unhighlight undo button 4/29/93 */
}

static void menuClearQuery(void) {
  clearQbuilder(ALL_STUFF);
  qbuildDisplay();
}
  
static void lNewKeySet(void)
{
  newKeySet("Query Builder Answer") ;
}

static void qbuildMenuExec(void)
{
  KEYSET oldSet, newSet ;
  void *look ;
  QBUILDGET("qbuildMenuExec");

  /* form on query interface window */
  if (querGraph != 0 && graphExists(querGraph)) {
    messStatus("Forming Query Command");
    if(!graphActivate(querGraph))
      messcrash("Cannot reactivate querGraph in qbuildMenuExec.\n");
    queryWindowShowQuery(makeQuery(resbuffer)); /* calls QUERGET inside*/
    queryCommandEdit();
    if(!graphActivate(qbuildGraph))
      messcrash("Cannot reactivate qbuildGraph in qbuildMenuExec.");
    return;
  }

  /* Execute query */
  if(!keySetActive(&oldSet, &look))
    lNewKeySet() ;
  if(keySetActive(&oldSet, &look)) { 
    messStatus("Searching") ;
    newSet = query(oldSet, makeQuery(resbuffer)) ;
    if (newSet != oldSet )  {
      if(keySetExists(query_undo_set) &&
	 query_undo_set != oldSet )
	keySetDestroy(query_undo_set) ;
      query_undo_set = oldSet ;
    }
    keySetShow(newSet,look) ;
  }
  else
    messout("First select a keySet, thank you") ;
  qbuildDisplay();  /* to unhighlight the "Search Locally" button */
}

static MENUOPT qbuildButtons[] =
{
  { qbuildDestroy,    "Quit" },
  { help,             "Help" },
  { graphPrint,       "Print" },
  { menuClearQuery,   "Clear" },
  { lNewKeySet,       "New KeySet" },
  { qbeCallup,        "By Example" },
  { qcommandCallup,   "Command Edit" },
  { qbuildMenuExec,   "Search Locally" },
  { qbuildButtonUndo, "Undo" },
  { 0, 0 }
} ;


/* Look at following functions to make sure not to disturb [#] indexing */
/* qbuildMore and MENU_MORE */
static void qbuildMore(void);
#define MENU_MORE   6    
static MENUOPT qbuildMenu[] =
{
  { qbuildDestroy,    "Quit" },
  { help,             "Help" },
  { graphPrint,       "Print Window" },
  { menuSpacer,       "" },
  { qbuildDisplay,    "Refresh" },
  { menuClearQuery,   "Clear Entries" },
  { qbuildMore,       "Show More" },  /* MUST BE ENTRY #6, see qbuildMore(), MENU_MORE */
  { menuSpacer,       "" },
  { lNewKeySet,       "New KeySet" },
  { qbeCreate,        "Query By Example" },
  { queryCreate,      "Query Command Editor" },
  { menuSpacer,       "" },
  { checkQueryFormat, "Test Syntax" },
  { qbuildMenuExec,   "Search Locally" },
  { queryCommandUndo, "Undo Search" },
  { 0, 0 }
} ;

static void qbuildMore(void)
{
  QBUILDGET("qbuildMore");
  if (TESTFLAG(SIMPLE_UI)) {
    RMFLAG(SIMPLE_UI);
    qbuildMenu[MENU_MORE].text = "Show Less";
  }
  else {
    SETFLAG(SIMPLE_UI);
    qbuildMenu[MENU_MORE].text = "Show More";
  }
  graphMenu(qbuildMenu) ;  
  qbuildDisplay();
}


/* open/show query-interface window  */ 
#define QB_OFFSET 7.5  
#define QC_OFFSET 12
static void qbuildDisplay(void)
{
  int i = 0, j = 0, k=0;
  int triBox;
  
  QBUILDGET("qbuildDisplay");
  
  graphRegister (PICK,(GraphFunc)qbuildPick) ;
  graphActivate(qbuild->graph) ;
  graphClear() ;
  
  graphButtons (qbuildButtons, 1, 1, 75) ;  /* put buttons on screen */

  graphTextFormat(BOLD);
  graphText("Query Builder:", 2, QB_OFFSET - 2); 
  graphTextFormat(PLAIN);
  graphText("(Click the right mouse button on the shaded boxes for choices.)", 
	    17, QB_OFFSET - 1.7); 
  
  qbuild->class_box = graphBoxStart();

  if (TESTFLAG(RELATED))  {
    graphText ("Find and select the ", 2, QB_OFFSET); 
    graphTextFormat(BOLD);
    graphText ("related", 22.4, QB_OFFSET); 
    graphTextFormat(PLAIN);
    graphRectangle(29.8, QB_OFFSET-0.3, 52.4, QB_OFFSET+1.3);
    classtagMenuSet(graphTextScrollEntry(qbuild->classtag_entry, 
					 BUFFER_SIZE - 1, 20, 30.5, QB_OFFSET, 
					 (EntryFunc) classtag_handler));

    qbuild->class_triBox = graphMenuTriangle(TRUE, 50.9, QB_OFFSET) ;
    
    graphText("items", 53, QB_OFFSET);
  }
  else {
    graphText ("Find and select ", 2, QB_OFFSET); 
    graphRectangle(17.6, QB_OFFSET-0.3, 40.2, QB_OFFSET+1.3);
    classtagMenuSet(graphTextScrollEntry(qbuild->classtag_entry, 
					 BUFFER_SIZE - 1, 20, 18.3, QB_OFFSET, 
					 (EntryFunc) classtag_handler));
    qbuild->class_triBox = graphMenuTriangle(TRUE, 38.7, QB_OFFSET) ;

    graphText("items", 41, QB_OFFSET);
  }    
  graphBoxFreeMenu(qbuild->class_triBox, 
		   (FreeMenuFunction) classtagMenuResult, classtagMenu);

  if (!TESTFLAG(SIMPLE_UI)) {
    graphTextFormat(BOLD);
    graphText ("from", 59, QB_OFFSET); 
    graphTextFormat(PLAIN);
    
    graphRectangle(63.3, QB_OFFSET-0.3, 80.9, QB_OFFSET+1.3);
    graphBoxFreeMenu(graphTextScrollEntry(qbuild->preclass_entry,
					  BUFFER_SIZE-1, 15, 64, QB_OFFSET, 
					  (EntryFunc) preclass_handler),
		     (FreeMenuFunction) preclassMenuResult, preclassMenu);
    qbuild->pre_triBox = graphMenuTriangle(TRUE, 79.4, QB_OFFSET) ;
    graphBoxFreeMenu(qbuild->pre_triBox, 
		     (FreeMenuFunction) preclassMenuResult, preclassMenu);

  }
  graphBoxEnd();
  
  qbuild->user_box = graphBoxStart(); 
  
  for (j = 0; j < arrayMax(qbuild->entries)/ENTRY_NUM; j++) {  /* j= # lines */
    k = ENTRY_NUM*j;  /* k is offset for indexing each line */
    i = 2 * j;        /* i is physical display location (vertical) */
    
    graphText("where", 4, QB_OFFSET + 2 + i);
    graphRectangle(9.8, QB_OFFSET + 1.7 + i, 32.4, QB_OFFSET + 3.3 + i);
    attrMenuSet(graphTextScrollEntry(ARR2STRING(qbuild->entries, 
						ATTRIBUTE + k), 
				     BUFFER_SIZE - 1, 20, 10.5, 
				     QB_OFFSET + 2 + i, 
				     (EntryFunc) attribute_handler));
    triBox = graphBoxStart();
    graphMenuTriangle(TRUE, 30.9, QB_OFFSET + 2 + i) ;
    graphBoxEnd();
    graphBoxFreeMenu(triBox, 
		     (FreeMenuFunction) attrMenuResult, attrMenu);

    graphRectangle(33.8, QB_OFFSET + 1.7 + i, 48.4, QB_OFFSET + 3.3 + i);
    graphBoxFreeMenu(graphTextScrollEntry(ARR2STRING(qbuild->entries, 
						CONDITION + k), 
				     BUFFER_SIZE - 1, 12, 34.5, 
				     QB_OFFSET + 2 + i, 
				     (EntryFunc) condition_handler),
		     (FreeMenuFunction) condMenuResult, condMenu);
    triBox = graphBoxStart();
    graphMenuTriangle(TRUE, 46.9, QB_OFFSET + 2 + i) ;
    graphBoxEnd();
    graphBoxFreeMenu(triBox, 
		     (FreeMenuFunction) condMenuResult, condMenu);

    graphRectangle(49.8, QB_OFFSET + 1.7 + i, 72.4, QB_OFFSET + 3.3 + i);
    graphBoxFreeMenu(graphTextScrollEntry(ARR2STRING(qbuild->entries, 
						     VALUE + k), 
					  BUFFER_SIZE - 1, 20, 50.5,
					  QB_OFFSET + 2 + i,
					  (EntryFunc) value_handler),
		     (FreeMenuFunction) valueMenuResult, valueMenu);
    triBox = graphBoxStart();
    graphMenuTriangle(TRUE, 70.9, QB_OFFSET + 2 + i) ;
    graphBoxEnd();
    graphBoxFreeMenu(triBox, 
		     (FreeMenuFunction) valueMenuResult, valueMenu);

    graphRectangle(73.8, QB_OFFSET + 1.7 + i, 84.4, QB_OFFSET + 3.3 + i);
    graphBoxFreeMenu(graphTextScrollEntry(ARR2STRING(qbuild->entries, 
						     CONJUNCTION + k), 
					  BUFFER_SIZE - 1, 8, 74.5, 
					  QB_OFFSET + 2 + i, 
					  (EntryFunc) conjunction_handler),
		     (FreeMenuFunction) conjMenuResult, conjMenu);
    triBox = graphBoxStart();
    graphMenuTriangle(TRUE, 82.9, QB_OFFSET + 2 + i) ;
    graphBoxEnd();
    graphBoxFreeMenu(triBox, 
		     (FreeMenuFunction) conjMenuResult, conjMenu);
  }      
  graphBoxEnd(); 
  graphRectangle(0.5, QB_OFFSET - 2.5 , 86, QC_OFFSET + i);
  graphTextBounds (80, (int) (k + 1.5*i + 2)) ;        /* see qbuildCreate */
  
  graphRedraw() ;
}

/********************  public routines   *****************************/

/* ordering routine used to alphabetize array lists with arraySort */
int textAlphaOrder(const void *a, const void *b)
{
  return lexstrcmp((const char *) arrp(*(const Array *) a, 0, char),
		   (const char *) arrp(*(const Array *) b, 0, char));
}
     

/* change the qbuild_selected class to class */
void qbuildNewClass(int classe)
{
  qbuild_selected_class = classe;
  newKeyList(qbuild_selected_class, &classtagkeyset);
}

/* destroy a messalloced FREEOPT menu */
void destroyMenu(FREEOPT *menup, int max)
{
  int i;
  char *temp ;
  for (i = 0; i  <= max; i++) {
    temp = (menup+i)->text ;
    messfree (temp); 
  }
  messfree(menup);
  menup = NULL;
}

/* Initializes and creates all data structure necessary for
 * query-window interface.  Also, calls qbuildDisplay() to pop up image.
 */
void qbuildCreate (void)
{
  QBUILD qbuild ;
  int i; 
  KEYSET oldSet = 0 ;
  void   *look = 0 ;
  
  if(graphActivate(qbuildGraph))
    {
      graphPop() ;
      return ;
    }
  
  qbuildGraph = displayCreate(DtQueryBuilder) ;
  if (!qbuildGraph)
    return ;
  
  qbuild=(QBUILD)messalloc(sizeof(struct QBUILDSTUFF));
  qbuild->magic = MAGIC;
  qbuild->graph = qbuildGraph ; /* provision for multi windows */
  qbuild->entries = arrayCreate(ENTRY_NUM, Array) ;
  qbuild->syntax = arrayCreate(ENTRY_NUM, Array) ;
  qbuild->class_box = 0;
  
  for (i=0; i < ENTRY_NUM; i++) {
    /* use  supplied array subroutines instead of messalloc/malloc */
    array(qbuild->entries, i, Array) = 	arrayCreate(BUFFER_SIZE, char);
    array(qbuild->syntax, i, Array) = arrayCreate(BUFFER_SIZE, char);
  }
  
  qbuild_selected_class = 0;
  qbuildClassMax() ;  /* init classlist */
  classtagkeyset = NULL;
  attrkeyset = NULL;

  keySetActive(&oldSet,&look);
  if (oldSet) {      /* if keyset exists, initialiaze the taglist */
    qbuild_selected_class = class(keySet(oldSet, 0));
    newKeyList(qbuild_selected_class, &classtagkeyset);
  }
  
  strcpy(qbuild->preclass_entry, "ALL DATA");
  strcpy(qbuild->preclass_syntax, "Find ");
  qbuild->pre_triBox = 0;
  strcpy(qbuild->classtag_entry, "ANY CLASS");
  strcpy(qbuild->classtag_syntax, "Text");
  qbuild->class_triBox = 0;
  strcpy(ARR2STRING(qbuild->entries, ATTRIBUTE), "ANY TAG");
  strcpy(ARR2STRING(qbuild->entries, CONDITION), "contains");
  strcpy(ARR2STRING(qbuild->entries, CONJUNCTION), "END"); 
  qbuild->entry_set = 0;
  qbuild->user_box = 1;
  qbuild->entry_index = VALUE;       /* initialize PICK to this ebox */
  SETFLAG(ALLTEXT);
  
  graphTextBounds (80,40) ;   /* needed to for text box sizing */
  graphRegister (DESTROY, qbuildDestroy) ;
  graphRegister (PICK, qbuildPick) ;
  graphRetitle ("Query Builder"); /* acedb 1.9.1 prepends "DtQuery :" */

  SETFLAG(SIMPLE_UI);
  qbuildMenu[MENU_MORE].text = "Show More";

  graphMenu(qbuildMenu) ;  
  graphAssociate(&MAGIC, qbuild);
  graphHelp ("Query_Builder") ;
  
  attr_sub_model = FALSE;
  attrMenus = NULL;
  qbuild->attrMenuArray = arrayCreate(1, FREEOPT *);
  classtagMenus = NULL;
  preclassMenuSet();
  
  qbuildDisplay();
  qbuildPick(INDEX2BOX( qbuild->user_box, qbuild->entry_index ));
}

/******************************************************************************************/
/******************************************************************************************/

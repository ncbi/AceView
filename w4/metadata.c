/*  File: metadata.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
    First run the script wtools/parsecode
    Then read the file generated that way called code.ace
    Then run this program
    This results in cross referencing the includes, tags and 
    classes used in the code
 * Exported functions:
     void metaCheckCode(void)
 * HISTORY:
 * Last edited: Dec  9 15:03 1998 (fw)
 * Created: Wed Apr 21 20:22:39 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: metadata.c,v 1.1.1.1 2002/07/19 20:22:56 sienkiew Exp $ */

#include "acedb.h"

#include "lex.h"
#include "a.h"
#include "bs.h"
#include "systags.h"
#include "session.h"
#include "sysclass.h"
#include "pick.h"

static int   _VInclude, _VSourceCode, _VTag ;     /* metadata class, not used by the kernel */

/********************************************************/

static void metaCheckTagConsistency(void)
{ KEY tag = 0 , metaTag = KEYMAKE(_VTag,0) ;
  int n ;
  OBJ obj ; char *cp ;

  if (iskey(KEYMAKE(_VTag,0)) != 2)
    return ;  /* i.e. the Tag model does not yet exists */
 
     /* tag to metatag injection */
  while (lexNext(_VSystem, &tag))
    { lexaddkey(name(tag), &metaTag, _VTag) ;
      obj = bsUpdate(metaTag) ;
      if (bsGetData(obj, _Parent_tag, _Text, &cp))
	{ if (lexstrcmp(cp, name(tag)))
	    messerror ("Tag %s has somehow been renamed as %s",
		       cp, name(tag)) ;
	  if (bsGetData(obj, _bsRight, _Int, &n) && n!= 458 && n!= 478)
	    { if ( n != tag )
		messerror ("Tag %s has been renumbered, old value %d, new %d",
			   name(tag), n, tag);
	    }
	}
      bsAddData(obj, _Parent_tag, _Text, name(tag)) ;
      n = tag ; /* to cast to int */
      bsAddData(obj, _bsRight, _Int, &n) ;
      bsSave (obj) ;  /* Those operation are often null, hence quick and harmless */
    }

   /* is it a surjection */

  metaTag = 0 ;
  while (lexNext(_VTag, &metaTag))
    if ((obj = bsCreate(metaTag)))
      { 
	if (bsGetData(obj, _Parent_tag, _Text, &cp) &&
	    bsGetData(obj, _bsRight, _Int, &n) &&
	     strncmp(cp,"__sys", 5) &&
	    n != 458 && n != 478 ) /* A fix, Because 458 is Semi_dominant, was demi-dominant */
	{ if (!lexword2key(cp, &tag, _VSystem))
	    messerror ("Tag %s %d has somehow disappeared from wspec",
		       cp, n) ;
	  if (lexstrcmp(cp, name(tag)))
	      messerror ("Tag %d has been renamed, old name %s, new name %s",
			 n, cp, name(tag)) ;
	  if (n != tag )
		messerror ("Tag %s has been renumbered, old value %d, new %d",
			   name(metaTag), n, tag);
	}
	bsDestroy (obj) ;  /* Those operation are often null, hence quick and harmless */
      }
}

/************************************************************/

void metaCheckCode(void)
{
  KEY key, codeText, code , metaTag , metaClass, incl ;
  char *cp, *cq , cutter;
  Stack s ;
  int level , nt = 0, nc= 0, ni = 0 ;
  OBJ Code ;


  if (!sessionGainWriteAccess())
    { 
      messout("Sorry, you cannot gain write access") ;
      return ;
    }

  lexaddkey ("Include", &key, _VMainClasses) ; _VInclude = KEYKEY (key) ;
  lexaddkey ("SourceCode", &key, _VMainClasses) ; _VSourceCode = KEYKEY (key) ;
  lexaddkey ("Tag", &key, _VMainClasses) ; _VTag = KEYKEY (key) ;
  
  metaCheckTagConsistency() ;  
  code = 0 ;

  while (lexNext(_VSourceCode, &code))
    {
      if((Code = bsUpdate(code)))
	{ 
	  if(bsFindTag(Code, _Includes))
	    bsRemove(Code) ;
	  if(bsFindTag(Code, _Uses_tags))
	    bsRemove(Code) ;
	  if(bsFindTag(Code, _Uses_class))
	    bsRemove(Code) ;

	  if (bsGetKey(Code, _SourceFile, &codeText) &&
	      (s = stackGet(codeText) ))
	    { stackCursor(s, 0) ;
	      while ((cp = stackNextText(s)))
		{ level = freesettext(cp, "") ;
		  while (freecard(level))
		    if ((cq = freeword()) && !strcmp(cq, "#include"))
		      { 
			freenext() ;
			cq = freeword() ;
			if (lexword2key(cq, &incl, _VInclude))
			  { ni++ ;
			    bsAddKey(Code, _Includes, incl) ;
			  }
			else
			  printf("Unrecognised include: %s in %s \n", cp, name(code)) ;
		      }
		    else while ( freenext(), 
				(cq = freewordcut(" {}()[],.;:=!#>?\"", &cutter), cq || cutter))
		      if (cq && *cq == '_' && *(cq+1) == 'V')
			{ if (lexword2key(cq+2, &metaClass, _VClass))
			   { nc++ ; bsAddKey(Code, _Uses_class, metaClass) ; }
			  else
			    fprintf(stderr, "Unrecognised class %s in %s \n", cq ? cq : "Null", name(code)) ;
			}
		      else if (cq && *cq == '_')
			{ if (lexword2key(cq+1, &metaTag, _VTag))
			    { nt++; bsAddKey(Code, _Uses_tags, metaTag) ; }
			  else
			    fprintf(stderr, "Unrecognised tag %s in %s \n", cq ? cq : "Null", name(code)) ;
			}
		}
	      stackDestroy(s) ;
	    }
	  bsSave(Code) ;
      
	}
    }

  printf ("#*#*#*# Found %d tags, %d classes %d inclusions \n\n", nt, nc, ni) ;
}

/***********************************/
/***********************************/
 

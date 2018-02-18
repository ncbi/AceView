/*  File: acediffsorted.c
 *  Author: R. Durbin and J. Thierry-Mieg
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: makes a difference file from two sorted ace files
 * Exported functions: <standalone>
 * Usage: acediff f1 f2 > f3
 * HISTORY:
 * Last edited: Jul  9 11:59 2003 (rnc)
 * * Aug 26 16:46 1999 (fw): added this header
 * * Dec 15 13:19 1994 (dok256): 
   built from acediff by Detlef Wolf (dok256) and Michael Schultz (dok408) from DKFZ.

   differences to acediff:

   - acediffs requires two sorted (lexstrcmp) input-files.
     the classnames must be sorted, and the object names in the classes must be sorted.
     the tags inside an object need not be sorted.
     the .ace files must no contain doubled objects from aliasing.


   - acediffs no longer needs to create and sort long temporary files;
     therefore it is faster than acediff and requires no temporary disk space.

   - acediffs can read compressed (.Z|.gz) input-files.




   what was changed in the code?


   - I added the subroutine getline_from_buffer (modification of getline)
     This sub reads the current object into the buffer and sorts it if required. 
     (because the data-lines inside objects are expected to be sorted)
     It returns the next line of the input-file.

     
   - all calls of getline in processDiff were changed to calls of getline_from_buffer
     This has the advantage, that the old code can still be used.

    - DNA and LongText Objects are now treated as objects with one line of data
      (in earlier versions, they were not read to memory but compared on the disk)

   - removed sub diff2ace, because the difference-file is written directly

   - removed cmpFileParts, printFilePart, because DNA and LongText 
     are now compared in memory.


     IMPORTANT:

    ........
    - A-classes:
      when reading first time pointers to the text portions of DNA and LongText
      objects in the input files are remembered and used at compare
      and at output time 
    - protection against line overflow, raised limit to 10000000

    - remove intermediate files
    - added -d option
 * Created: Thu Aug 26 16:43:26 1999 (fw)
 * CVS info:   $Id: acediffsorted.c,v 1.9 2017/01/19 21:52:27 mieg Exp $
 *-------------------------------------------------------------------
 */

#include "regular.h"
#include <wh/mytime.h>
#include "keyset.h"
#include "chrono.h"
#include "call.h"
#include "utils.h"
#include "lex.h"

/* arrays are predefined to this size */
#define LINES_IN_OBJ 500
#define DATA_SIZE 50
/* stacks are preextended to hold DATA_SIZE bytes */
/* constant to separate lines from the two input files */
/* #define MAXLINES 10000000 */
/* no longer needed */

extern char* freepos(void);
static int isWhite[256];
FILE *inputA, *inputB;           /* input files; open during whole execution */
char *LongTextEnd= "\n***LongTextEnd***\n"; /* terminator for LongText class */
int termPosFinal;                           /* = strlen(LongTextEnd) */
BOOL debugFlag = FALSE;

#define Text(z) (stackText(z,0))

struct Buffer  {
  Stack Obj_name, quotedObj_name; 
  /* Ob_name: class objname; quotedObj_name: class "obj_name" */
  Array data; 
  int line, maxlines;
  
  FILE *in; 
  BOOL end_of_Buffer, isAObj, isSorted;
  long objCount;
  char *filename;
};



void initBuffer(struct Buffer *buf, FILE *in, char *name) {
  buf->in = in;
  buf->line = buf->maxlines = buf->objCount = 0;
  buf->data = arrayCreate(LINES_IN_OBJ, Stack);
  buf->Obj_name = stackCreate(128);
  buf->quotedObj_name = stackCreate(128);
  buf->isAObj = buf->end_of_Buffer = FALSE;
  buf->filename = name;

} 


int order(const void *a, const void *b)
{
  Stack x = *(const Stack *)a, y = *(const Stack *)b;
  return lexstrcmp(Text(x), Text(y));
}


void oneLine_next_Object(struct Buffer *buf) 
/* read next object from .ace file into buffer  buf
   while joining conitutation lines and
   labeling every line with line number and object name;
   replace A-class object data with pointer into input files */
{
  BOOL  isAObj, isLongText = 0 ;
  int   termPos;            /* current index in terminator string */
  unsigned int   ch, lastch; /* characters read in DNA scan */
  BOOL  inObject = FALSE;
  BOOL	escaped, outQuotes, inSpace;
  Stack current_obj, current_obj_q; 
  char  *word, *s;
  char tmp[2];


  tmp[1] = 0;
  buf->line = 0; 
  buf->maxlines = 0;
  buf->isSorted = FALSE;

  freespecial ("\n\"\\") ;
  while (freeread (buf->in)) { /* read only one object */ 
    if (inObject) {
      if (*freepos()) { 
	escaped = FALSE; outQuotes = TRUE; inSpace = FALSE; 
	
	array(buf->data, buf->line, Stack) = 
	  stackReCreate(array(buf->data, buf->line, Stack), DATA_SIZE);
	
	for (s = freepos(); *s; ++s) { 
	  if (escaped) escaped = FALSE;
	  else { 
	    if (outQuotes) {
	      if (isWhite[(int)*s])
		inSpace = TRUE;
	      else if (inSpace) { 
		catText(array(buf->data, buf->line, Stack), " ");  
		inSpace = FALSE;
	      }
	    }
	    if (*s == '"')
	      outQuotes ^= 0x01;
	    escaped = (*s == '\\');
	  }
	  if (!inSpace) {
	    /*fputc (*s, buf->out);  */
	    tmp[0] = s[0];
	    catText(array(buf->data, buf->line, Stack), tmp);  
	  }
	}
	
	buf->line++; /* goto next buffer line */
	buf->maxlines++;
      } else {
	inObject = FALSE;
	buf->line = 0; 
	if (buf->maxlines) return; /* if buffer contains an object -> return it
				      else go on reading. (no confusion with empty objs */
      }
    }
    else if ((word = freeword())) { 
      
      /* determine storage type from classname*/
      if (strcmp(word,"LongText")==0) {
	isAObj= TRUE; isLongText= TRUE; 
      } else if (strcmp(word,"DNA")==0) {
	isAObj= TRUE; isLongText= FALSE;
      } else isAObj= FALSE;
      
      /* construct object name */
      current_obj = stackCreate(128);  
      current_obj_q = stackCreate(128);  
      catText(current_obj, Text(buf->Obj_name));  
      catText(current_obj_q, Text(buf->quotedObj_name));  
      
      stackClear(buf->Obj_name);
      stackClear(buf->quotedObj_name);
      buf->objCount++;
      catText (buf->Obj_name, word); 
      catText (buf->quotedObj_name, word);

      freestep (':');
      if (!(word = freeword())) 
	messcrash ("Error with object name in line %d", freeline(buf->in));
      
      catText (buf->Obj_name, " "); 
      catText (buf->quotedObj_name, " "); 
/*    catText (buf->Obj_name, freeprotect(word));  */
      catText (buf->Obj_name, word); 
      /* changed because the quotes destroyed the order, but as outquotes are
	 needed I put the quoted name in buf->quotedObj_name. I know this is not
	 very elegant but it works. */
      
      catText (buf->quotedObj_name, freeprotect(word));  
      
      if (lexstrcmp(Text(current_obj), Text(buf->Obj_name)) != -1) {
	printf("----\n%s->%s\n", Text(current_obj_q), Text(buf->quotedObj_name));
	messcrash("Input files are not sorted before line %d. file: %s\n", freeline(buf->in), buf->filename);
      }
      
      stackDestroy(current_obj);
      stackDestroy(current_obj_q);

      inObject = TRUE;
      buf->isAObj = FALSE;
	
      /* process A-objects here */
      if (isAObj) {
	array(buf->data, 0, Stack) = 
	  stackReCreate(array(buf->data, buf->line, Stack), DATA_SIZE);
	
	/* freeread() reads line-wise */
	/* , so now at new line = start of data; */
	/* "steal" intervening A-class data */
	if (isLongText) {
	  termPos=1;
	  while ( (ch=getc(buf->in)) != (unsigned char)EOF) {
	    tmp[0] = ch; 
	    catText(array(buf->data, 0, Stack), tmp);  
	    if (ch==LongTextEnd[termPos]) {
	      if (++termPos == termPosFinal) break;
	    } else termPos = (ch==LongTextEnd[0]) ? 1 : 0;
	  }
	} else { /* DNA */
	  lastch='\n';
	  while ((ch=getc(buf->in)) != (unsigned char)EOF) {
	    /* DNAs are now stored in Buffer */
	    tmp[0] = ch; 
	    catText(array(buf->data, 0, Stack), tmp); 
	    if (ch=='\n' && lastch =='\n') break;
	    lastch=ch;
	  }
	}
	buf->isAObj = TRUE;
	
	buf->line++; 
	buf->maxlines++;
	
	return;
      } /* is A-Object */
    }
    
  }
  
/*  if (line > MAXLINES) messcrash ("too many lines in input file. sorry."); */


/*   fclose (in); -- needs to stay open until the end */
/*   fclose (out); */
}

BOOL getLine_from_Buffer(struct Buffer *buf, Stack object, Stack quoted_object, Stack data){

  stackClear(object);
  stackClear(quoted_object);
  stackClear(data);


  if (buf->line < buf->maxlines) {

    pushText(data, Text(arr(buf->data, buf->line, Stack))); 

    pushText(object, Text(buf->Obj_name));  
    pushText(quoted_object, Text(buf->quotedObj_name));
    buf->line++;
    return TRUE;

  } else { /* buf->line == buf->maxlines */
    if (feof(buf->in)) {
      buf->end_of_Buffer = TRUE;
      return FALSE;
    }

    oneLine_next_Object(buf);
    if (feof(buf->in) && !buf->maxlines) { /* if the last object has no data ... */
      buf->end_of_Buffer = TRUE;
      return FALSE;
    } 


    pushText(object, Text(buf->Obj_name));
    pushText(quoted_object, Text(buf->quotedObj_name));
    pushText(data, Text(arr(buf->data, 0, Stack))); 
    buf->line = 1;

    return TRUE;
  }
}


BOOL getLine (FILE *in, int *line, Stack obName, Stack data)
{
  char* cp;
  int len;

  stackClear (obName);	
  stackClear (data);	

  freespecial ("\n\"\\") ;
  if (!freeread (in)) return FALSE;
  freeint (line);
  catText (obName, freeword());
  catText (obName, " ");
  catText (obName, freeprotect (freeword()));
  cp = freepos();
  if ((len = strlen(cp) > 1) && cp[len-1] == '\"' && cp[len-2] != '\\')
    cp[len-1] = 0;
  catText (data, cp);
	
  return TRUE;
}

  
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void processDiff(struct Buffer *bufA, struct Buffer *bufB, FILE *out) 
{
  int	cmpClass, cmpData = 0 ;
  /* reicht int? */
  int   cmpAObj;   /* 0 if A-objects are equal, 1 else */
/*  long  posA, posB, lenA, lenB; */
  Stack dataA, dataB, objA, objB, objAq, objBq;
  Stack obCurrA ;


  obCurrA = stackCreate (128);

  objA = stackCreate(128);
  objAq = stackCreate(128);
  dataA = stackCreate(128);
  objB = stackCreate(128);
  objBq = stackCreate(128);
  dataB = stackCreate(128);

  getLine_from_Buffer(bufA, objA, objAq, dataA); 
  getLine_from_Buffer(bufB, objB, objBq, dataB);/*DATAB*/ 

  while (!bufA->end_of_Buffer || !bufB->end_of_Buffer) {  
  /* check for End-of-File and if both have data, compare object-id */

    if (bufB->end_of_Buffer) cmpClass = -1;
    else if (bufA->end_of_Buffer) cmpClass = 1;
    else cmpClass = lexstrcmp(Text(objA), Text(objB));


    /* if same object, then check for storage type */
    cmpAObj = 0;

    /* if same object, then compare data depending on storage type */
    if (cmpClass==0) {
      if (bufA->isAObj) {
	cmpAObj = lexstrcmp(Text(dataA), Text(dataB));
	
	cmpData = 0; 
      } else { /* objA == objB */
	cmpAObj = 0;
	cmpData = lexstrcmp(Text(dataA), Text(dataB)); 


	if (cmpData != 0 && bufA->isSorted == FALSE) {
	  /* if the data inside two objects is not equal, the rest of the
	     buffer has to be sorted */

	  bufA->line--; bufB->line--; 
	  arraySortPos(bufA->data, bufA->line, order);
	  arraySortPos(bufB->data, bufB->line, order);

	  bufA->isSorted = TRUE;
	  bufB->isSorted = TRUE;

	  getLine_from_Buffer(bufA, objA, objAq, dataA);
	  getLine_from_Buffer(bufB, objB, objBq, dataB);

	  cmpData = lexstrcmp(Text(dataA), Text(dataB)); 
	}
	  
      }
    }

    /* copy record to output, depending on previous comparison 
       and advance consumed input streams */


    if (cmpClass < 0 || (cmpClass == 0  &&  cmpData < 0)) { 
/*      fprintf (out, "%8d %s -D %s\n", lineA, Text(objA), Text(dataA)); */
      if (lexstrcmp(Text(objA), Text(obCurrA))) {
	stackClear(obCurrA);
	catText(obCurrA, Text(objA));
	if (!strncmp(Text(objA), "DNA", 3) || !strncmp(Text(objA), "LongText", 8)){ 

/* insert LongText */
	  printf ("\n-D %s\n\n", Text(objAq));
	} else 
	  printf ("\n%s\n-D %s\n", Text(objAq), Text(dataA));
      } else {
	printf("-D %s\n", Text(dataA));
      }
      getLine_from_Buffer (bufA, objA, objAq, dataA);
    }
    else if (cmpAObj || cmpClass > 0 || (cmpClass == 0  &&  cmpData > 0)) {

/*      fprintf(out, "%8d %s %s\n", MAXLINES+lineB, Text(objB), Text(dataB)); */
      if (lexstrcmp(Text(objB), Text(obCurrA))) {
	stackClear(obCurrA);
	catText(obCurrA, Text(objB));
	printf("\n%s\n%s\n", Text(objBq), Text(dataB));
      } else {
	printf("%s\n", Text(dataB));
      }
      getLine_from_Buffer(bufB, objB, objBq, dataB);
      if (cmpAObj) getLine_from_Buffer(bufA, objA, objAq, dataA);
    } else { 
      getLine_from_Buffer(bufA, objA, objAq, dataA);
      getLine_from_Buffer(bufB, objB, objBq, dataB); 
    }
  }

    

  stackDestroy (objA);
  stackDestroy (dataA);
  stackDestroy (objB);
  stackDestroy (dataB); 

  fclose (out);
}

  


int main(int argc, char **argv)
{
  FILE *tmpC, *temp;
  char *cq, *cp;
  struct Buffer BufA, BufB; 

	
  if (argc == 4) {
    if (strcmp(argv[1], "-d")==0) {
      argv[1]=argv[2]; /* shift down */
      argv[2]=argv[3];
      debugFlag=TRUE;
    }
    else argc=1;  /* force usage message */
  }

  if (argc != 3 && argc != 4)
    { fprintf (stderr, 
"Usage  : acediffs [-d] oldfile[.gz] newfile[.gz] > diff\nPurpose: loading file newfile into a given database has the same\n         effect as loading oldfile followed by loading diff.\nBoth files have to be sorted (like acedumps) and must not contain doubled objects.\n");
      return -1;
    }
/* option -d turns on debugFlag, but is currently not in use */

  termPosFinal=strlen(LongTextEnd);
  freeinit ();

  isWhite[' '] = isWhite['\t'] = 1;

  inputA = fopen (argv[1],"r");
  inputB = fopen (argv[2],"r"); 

  if (!inputA || !inputB) { 
    fprintf (stderr, "problems opening files - aborting\n");
    return 1;
  }

  cp = argv[1];

  cq = cp + strlen(cp) - 2; 
  if (cq > cp && !strcmp(cq, ".Z")) { /* a compressed file */
    temp = (FILE*) callScriptPipe ("zcat",cp) ;
  } else if (cq-- && cq > cp && !strcmp(cq, ".gz")) { 
    temp = (FILE*) callScriptPipe ("gunzip",messprintf(" -c %s", cp)) ;
  } else {
    temp = fopen (cp,"r");
  }

  fclose(inputA);
  inputA = temp;


  cp = argv[2];

  cq = cp + strlen(cp) - 2; 
  if (cq > cp && !strcmp(cq, ".Z")) { /* a compressed file */
    temp = (FILE*) callScriptPipe ("zcat",cp) ;
  } else if (cq-- && cq > cp && !strcmp(cq, ".gz")) { 
    temp = (FILE*) callScriptPipe ("gunzip",messprintf(" -c %s", cp)) ;
  } else {
    temp = fopen (cp,"r");
  }

  fclose(inputB);
  inputB = temp;

  fprintf (stderr,"performing diff from %s to %s...\n", argv[1], argv[2]);
  fflush (stderr);

  /* connect the buffers to the input-files */
  initBuffer(&BufA, inputA, argv[1]);
  initBuffer(&BufB, inputB, argv[2]);

  tmpC = fopen ("tempC1","w");
  fprintf (stdout,"// acediff difference from %s to %s\n", argv[1],argv[2]);

  /* processDiff now reads from buffers instead of files */
  processDiff (&BufA, &BufB, tmpC);  
  if (!debugFlag) {unlink("tempA2"); unlink("tempB2");};

  if (!debugFlag) unlink("tempC1"); /* needed? */

  tmpC = fopen ("tempC2","r");


  fprintf (stdout,"\n// end of file\n");

  fprintf(stderr, "%ld objects read from %s\n", BufA.objCount, argv[1]);
  fprintf(stderr, "%ld objects read from %s\n", BufB.objCount, argv[2]);

  return 0;
}




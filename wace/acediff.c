/*  File: acediff.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1998
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
 * Description: makes a difference file from two ace files
 *              Usage: acediff f1 f2 > f3
 *              
 *              
 * Exported functions: None, stand alone utility.
 * HISTORY:
 * Last edited: May  9 12:11 2002 (edgrif)
 * * Dec 10 11:03 1998 (edgrif): Make declaration of main correct, and add
 *              this formal header, sort out dangling else's.
 * * Dec 15 13:19 1994 (dok256):
 *              A-classes:
 *                when reading first time pointers to the text portions of DNA and LongText
 *                objects in the input files are remembered and used at compare
 *                and at output time 
 *                protection against line overflow, raised limit to 10000000
 *                remove intermediate files
 *                added -d option
 *              
 * Created: Thu Dec 10 11:01:22 1998 (edgrif)
 * CVS info:   $Id: acediff.c,v 1.9 2014/02/20 00:12:22 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wh/array.h>

/* constant to separate lines from the two input files */
#define MAXLINES1 10000000
#define MAXLINES2 20000000

static int isWhite[256] ;
FILE *inputA, *inputB;           /* input files; open during whole execution */
char *LongTextEnd= "\n***LongTextEnd***\n"; /* terminator for LongText class */
int termPosFinal;                           /* = strlen(LongTextEnd) */
BOOL debugFlag = FALSE ;


#define Text(z) (stackText(z,0))



static void aceDiffAbort(const char *reason)
{
  perror(reason);
  exit(1);
}


/* compare portions of files: return 0 if equal, 1 if not equal */
static int cmpFileParts (FILE *inA, long posA, long lenA, FILE *inB, long posB, long lenB) 
{
  int a, b, num_chars ;
  int result = 0 ;

  if (lenA != lenB)
    result = 1 ;
  else
    {
      if (fseek (inA, posA, SEEK_SET) != 0)
	aceDiffAbort("fseek error on A");
      if (fseek (inB, posB, SEEK_SET) != 0)
	aceDiffAbort("fseek error on B");

      for (num_chars = lenA ; num_chars && result == 0 ; num_chars--)
	{ 
	  a = getc(inA) ;
	  b = getc(inB) ;
	  
	  if ((a == EOF) || (b == EOF))
	    aceDiffAbort("Error reading data in cmpFileParts");
	  else if (a != b)
	    result = 1 ;
	}
    }

  return result;
}

static void printFilePart (FILE *f, long pos, long len, FILE *out)
{  /* print len bytes starting at pos from file f on out */
  if (fseek (f, pos, SEEK_SET) != 0)
    aceDiffAbort("fseek error in printFilePart");
  for (; len; len--)
    if (putc(getc(f),out) == EOF)
      aceDiffAbort("Error writing data in printFilePart");
}


/* .ace file in to out while joining conitutation lines and
 * labeling every line with line number and object name;
 * replace A-class object data with pointer into input files */
static void oneLine (FILE *in, FILE *out)
{
  BOOL  isAObj, isLongText = FALSE, isDNA = FALSE, isPeptide  = FALSE ;
  long  startPos, textLen;  /* startposition and length of A-type object */
  int   termPos;            /* current index in terminator string */
  unsigned int   ch, lastch; /* characters read in DNA scan */
  BOOL  inObject = FALSE ;
  BOOL	escaped, outQuotes, inSpace ;
  Stack obName ;
  char  *word, *s ;
  int   line = 0 ;

  obName = stackCreate(128) ;

  freespecial ("\n\"\\") ;
  while (freeread (in))
    {
      ++line ;

      if (inObject)
	{
	  if (*freepos())
	    {
	      fprintf (out,"%8d %s ", line, Text(obName)) ;

	      escaped = FALSE ; outQuotes = TRUE ; inSpace = FALSE ;

	      for (s = freepos() ; *s ; ++s)
		{
		  if (escaped)
		    {
		      escaped = FALSE ;

		      /* Ace files mark continued lines with "\n\" so if we see
		       * a '\' followed by 'n' we need to increment the line count
		       * as the freeread() call concatenates such lines. */
		      if (*s == 'n')
			line++ ;
		    }
		  else
		    {
		      if (outQuotes)
			{
			  if (isWhite[(int)*s])
			    inSpace = TRUE ;
			  else if (inSpace)
			    { 
			      if (putc(' ', out) == EOF)
				aceDiffAbort("Error writing output");
			      inSpace = FALSE ;
			    }
			}
		      if (*s == '"')
			outQuotes ^= 0x01 ;
		      escaped = (*s == '\\') ;
		    }
		  if (!inSpace)
		    if (fputc(*s, out) == EOF)
		      aceDiffAbort("Error writing output"); 
		}
	      if (fputc('\n', out) == EOF)
		aceDiffAbort("Error writing output"); 
	    }
	  else
	    inObject = FALSE ;
	}
      else if ((word = freeword()))
	{ 
	  stackClear (obName) ;
	
	  /* determine storage type from classname*/
	  if (strcmp(word,"LongText")==0)
	    {
	      isAObj= TRUE; isLongText = TRUE; isDNA = FALSE; isPeptide = FALSE ;
	    } 
	  else if (strcmp(word,"DNA")==0)
	    {
	      isAObj= TRUE; isLongText = FALSE; isDNA = TRUE; isPeptide = FALSE ;
	    }
	  else if (strcmp(word,"Peptide")==0)
	    {
	      isAObj= TRUE; isLongText = FALSE; isDNA = FALSE; isPeptide = TRUE ;
	    }
	  else 
	    isAObj= FALSE;

	  /* construct object name */
	  catText (obName, word) ;
	  freestep (':') ;
	  if (!(word = freeword()))
	    messcrash ("Error with object name line %d",line) ;
	  catText (obName," ") ;
	  catText (obName, freeprotect (word)) ;
	  inObject = TRUE ;

	  /* process A-objects here */
	  if (isAObj)
	    {
	      startPos=ftell(in); /* freeread() reads line-wise, so now at new line = start of data; "steal" intervening A-class data */
	      textLen=0;
	      if (isLongText)
		{
		  termPos=1;
		  while ( (ch=getc(in)) != EOF) 
		    {
		      textLen++;
		      if (ch=='\n')
			line++;
		      if (ch==LongTextEnd[termPos]) 
			{
			  if (++termPos == termPosFinal) break;
			}
		      else 
			termPos= (ch==LongTextEnd[0])? 1 : 0;
		    }
		  if (errno) aceDiffAbort("Unexpected EOF");
		} 
	      else if (isDNA || isPeptide) 
		{
		  lastch='\n';
		  while ( (ch=getc(in)) != EOF) 
		    {
		      if (ch=='\n')
			line++;
		      if (ch=='\n' && lastch =='\n')
			break;
		      lastch=ch;
		      textLen++;
		    }
		  if (errno) aceDiffAbort("Unexpected EOF"); 
		}
	      if (fprintf (out,"%8d %s %ld %ld\n", 
			   line, Text(obName), startPos, textLen) == 0)
		aceDiffAbort("IO error");
	      inObject = FALSE;
	    } /* is A-Object */

	}
    }

  fprintf (stderr, "%d lines read",line) ;
  if (line > MAXLINES1) messcrash ("too many lines in input file. sorry.") ;
  stackDestroy (obName) ;
  /*   fclose (in) ; -- needs to stay open until the end */
  if (fclose (out) == EOF)
    aceDiffAbort("IO error");

  return ;
}

static BOOL getLine (FILE *in, int *line, Stack obName, Stack data)
{
  char* cp ;
  int len ;

  stackClear (obName) ;
  stackClear (data) ;

  freespecial ("\n\"\\") ;
  if (!freeread (in))
    return FALSE ;
  freeint (line) ;
  catText (obName, freeword()) ;
  catText (obName, " ") ;
  catText (obName, freeprotect (freeword())) ;
  cp = freepos() ;
  if ((len = strlen(cp) > 1) 
      && cp[len-1] == '\"' && cp[len-2] != '\\')
    cp[len-1] = 0 ;
  catText (data, cp) ;
  return TRUE ;
}

static void processDiff (FILE *inA, FILE *inB, FILE *out)
{
  int	lineA, lineB, cmpClass, cmpData ;
  BOOL  isAClass ;  /* 0, if B-class object, 1 if A-class object */
  int   cmpAObj ;   /* 0 if A-objects are equal, 1 else */
  long  posA, posB, lenA, lenB;
  Stack	obNameA, obNameB, dataA, dataB ;

  obNameA = stackCreate (64) ;
  dataA = stackCreate (64) ;
  obNameB = stackCreate (64) ;
  dataB = stackCreate (64) ;

  getLine(inA, &lineA, obNameA, dataA) ;
  getLine(inB, &lineB, obNameB, dataB) ;

  while (!feof (inA) || !feof (inB))
    { 

      /* check for End-of-File and if both have data, compare object-id */
      if (feof(inB))
	cmpClass = -1 ;
      else if (feof(inA))
	cmpClass = 1 ;
      else
	{
	  cmpClass = strcasecmp (Text(obNameA), Text(obNameB)) ;
	  if (! cmpClass)
	    cmpClass = strcasecmp (Text(obNameA), Text(obNameB)) ;
	}

      /* if same object, then check for storage type */
      cmpAObj = 0;
      isAClass = FALSE;
      if (cmpClass==0) {
	if (strncasecmp(Text(obNameA),"LongText ",9)==0 ||
	    strncasecmp(Text(obNameA),"DNA ",4)==0 ||
	    strncasecmp(Text(obNameA),"Peptide ",8)==0) 
	  isAClass= TRUE;
      }

      /* if same object, then compare data depending on storage type */
      if (cmpClass==0) {
	if (isAClass) {
	  sscanf(Text(dataA),"%ld %ld", &posA, &lenA);
	  sscanf(Text(dataB),"%ld %ld", &posB, &lenB);
	  cmpAObj = cmpFileParts (inputA, posA, lenA, inputB, posB, lenB);
	  cmpData = 0;
	}
	else {
	  cmpAObj = 0 ;
	  cmpData = strcasecmp (Text(dataA), Text(dataB)) ;
	  if (! cmpData)
	    cmpData = strcmp (Text(dataA), Text(dataB)) ;
	}
      }

      /* copy record to output, depending on previous comparison 
       and advance consumed input streams */
      if (cmpClass < 0 || (cmpClass == 0  &&  cmpData < 0))
	{ 
	  if (fprintf (out, "%8d %s -D %s\n", 
		       MAXLINES1+lineA, Text(obNameA), Text(dataA)) == 0)
	    aceDiffAbort("IO error") ;
	  getLine (inA, &lineA, obNameA, dataA) ;
	}
      else if (cmpAObj || cmpClass > 0 || (cmpClass == 0  &&  cmpData > 0))
	{ 
	  if (fprintf (out, "%8d %s %s\n", 
		       MAXLINES2+lineB, Text(obNameB), Text(dataB)) == 0)
	    aceDiffAbort("IO error");
	  getLine (inB, &lineB, obNameB, dataB) ;
	  if (cmpAObj) getLine (inA, &lineA, obNameA, dataA);
	}
      else
	{ getLine (inA, &lineA, obNameA, dataA) ;
	getLine (inB, &lineB, obNameB, dataB) ;
	}
    }

  stackDestroy (obNameA) ;
  stackDestroy (dataA) ;
  stackDestroy (obNameB) ;
  stackDestroy (dataB) ;
  if (fclose (inA) == EOF)
    aceDiffAbort("IO error");
  if (fclose (inB) == EOF) 
    aceDiffAbort("IO error");
  if (fclose (out) == EOF)
    aceDiffAbort("IO error");

  return ;
}

static void diff2ace (FILE *in, FILE *out)
{
  int line  ;
  long pos, len;
  Stack obName, obCurr, data ;

  obName = stackCreate (64) ;
  obCurr = stackCreate (64) ;
  data = stackCreate (64) ;

  while (getLine (in, &line, obName, data)) { 

    if (strcmp (Text(obName),Text(obCurr))) { /* new object starts */

      stackClear (obCurr) ;
      catText (obCurr, Text(obName)) ;
       if (putc('\n',out) == EOF)
	 aceDiffAbort("Error writing data in diff2ace");
      if (strncmp(Text(obName),"LongText ",9)==0 ||
	  strncmp(Text(obName),"DNA ",4)==0 ||
	  strncmp(Text(obName),"Peptide ",8)==0) { /* storage type A */
	if (strncmp(Text(data), "-D",2)==0) 
	  {
	    if (fprintf (out,"-D %s\n", Text(obName)) == 0)
	      aceDiffAbort("IO error") ;
	  }
	else { /* print A-class data */
	  if (fprintf (out,"%s\n", Text(obName)) == 0)
	    aceDiffAbort("IO error") ;
	  sscanf(Text(data),"%ld %ld", &pos, &len);
	  printFilePart(inputB, pos, len, out);
	}
	/* if same A-object is twice in file */
	catText (obCurr, "0x01") ;
	continue;
      }

      if (fprintf (out,"%s\n", Text(obName)) == 0)
	aceDiffAbort("IO error");
    }
    if (fprintf (out,"%s\n", Text(data)) == 0)
      aceDiffAbort("IO error");
  }
  
  stackDestroy (obName) ;
  stackDestroy (obCurr) ;
  stackDestroy (data) ;
  if (fclose (in) == EOF)
    aceDiffAbort("IO error");
  /*  fclose (out) ; */
}

int main(int argc, char **argv)
{
  FILE *tmpA, *tmpB, *tmpC ;

  if (argc == 4)
    {
      if (strcmp(argv[1], "-d")==0)
	{
	  argv[1]=argv[2]; /* shift down */
	  argv[2]=argv[3];
	  debugFlag=TRUE;
	  --argc ;
	}
      else argc=1;  /* force usage message */
    }

  if (argc != 3)
    {
      fprintf (stderr, "Usage  : acediff [-d] oldfile newfile > diff\n"
      "Purpose: loading file newfile into a given database has the same\n"
      "effect as loading oldfile followed by loading diff.\n"
      "    -d (debug) keep intermediary files temp[ABC][12]\n") ;
      return(EXIT_FAILURE) ;
    }

  termPosFinal=strlen(LongTextEnd);
  freeinit () ;

  isWhite[' '] = isWhite['\t'] = 1 ;

  inputA = fopen (argv[1],"r") ;
  inputB = fopen (argv[2],"r") ;
  tmpA = fopen ("tempA1","w") ;
  tmpB = fopen ("tempB1","w") ;
  if (!inputA || !inputB || !tmpA || !tmpB)
    {
      fprintf (stderr, "problems opening files - aborting\n") ;
      return(EXIT_FAILURE) ;
    }

  fprintf (stderr,"reading file 1 - ") ;
  fflush (stderr) ;
  oneLine (inputA,tmpA) ;
  fprintf (stderr," - sorting file 1\n") ;
  system ("LC_ALL=C sort -b -f -u -T . -k 2 tempA1 -o tempA2") ;
  if (!debugFlag)
    unlink("tempA1");

  fprintf (stderr,"reading file 2 - ") ;
  fflush (stderr) ;
  oneLine (inputB,tmpB) ;
  fprintf (stderr," - sorting file 2\n") ;
  system ("LC_ALL=C sort -b -f -u -T . -k 2 tempB1 -o tempB2") ;
  if (!debugFlag)
    unlink("tempB1");

  fprintf (stderr,"performing diff - ") ;
  tmpA = fopen ("tempA2","r") ;
  tmpB = fopen ("tempB2","r") ;
  tmpC = fopen ("tempC1","w") ;
  if (!tmpA || !tmpB ||!tmpC)
    {
      fprintf (stderr, "problems opening intermediate files - aborting\n") ;
      return(EXIT_FAILURE) ;
    }

  processDiff (tmpA, tmpB, tmpC) ;
  if (!debugFlag)
    {
      unlink("tempA2");
      unlink("tempB2");
    }

  fprintf (stderr,"sorting diff\n") ;
  /* if (system ("LC_ALL=C sort -b -f -T . +1 -3 +0 -1 tempC1 -o tempC2") != 0) */
  if (system ("LC_ALL=C sort -b -f -T . -k 2 -k 4r  -k 1 -k 2r  tempC1 -o tempC2") != 0)
    aceDiffAbort("Trouble running sort");
  if (!debugFlag)
    unlink("tempC1");
  
  tmpC = fopen ("tempC2","r") ;
  if (!tmpC)
    aceDiffAbort("Cannot open tempC2");
  
  fprintf (stdout,"// acediff difference from %s to %s\n",
	   argv[1],argv[2]) ;
  
  diff2ace (tmpC, stdout) ;
  if (!debugFlag)
    unlink("tempC2");

  fprintf (stdout,"\n// end of file\n") ;

  return(EXIT_SUCCESS) ;
}

 


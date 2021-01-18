/*  Last edited: Aug 10 14:50 1996 (esr) */

/* $Id: blxparser.c,v 1.10 2020/05/30 16:50:34 mieg Exp $ */

/* blxparser - parses MSPcrunch output for blixem
 *
 * by Erik Sonnhammer
 *

  Date      Modification
--------  ---------------------------------------------------
93-05-17  Created
94-01-29  Added SEQBL format support.
94-01-30  Added proper option parsing and pipe'ability.
94-03-27  Added (fudged) Tblastn support.
95-02-06  Added (fudged) Tblastx support.
          Added autmatic mode detection.
95-06-21  Query sequence parsing that allows spaces and *
-------------------------------------------------------------
*/

#include "regular.h"
#include "graph.h"
#include "blxview.h"
#include <ctype.h>


#define MAXLINE 4096

enum { EXBLX, SEQBL };
int HSPgaps=0;

int isAlnRes (char ch)
{
    if (isalpha(ch) || ch == '*' || (HSPgaps && isdigit(ch)) )
	return 1;
    else
	return 0;
}

int isAlnNonres (char ch)
{
    if (ch == '.')
	return 1;
    else
	return 0;
}

/* BLXREADHSP parses an exblx/seqbl file and calls blixem */
Graph blxreadhsp(FILE *seqfile, FILE *exblxfile, char *qname, int dispstart, int qoffset, char *opts, int *argc, char **argv)
{
    MSP  *msp=0, *dbhits, DBhead;
    char  line[MAXLINE+1], sname[MAXLINE+1];
    int   format, i, qlen, slen;
    char *q;    /* The genomic sequence (query=the translated seq) */
    char *c, *cp, *cq, ch ;
    /* char *pipe1, *pipe2; */


    /* Read in query sequence */
    if (seqfile == stdin) {
	Array arr = arrayCreate(5000, char);
	int i=0;

	if (!fgets(line, MAXLINE, seqfile)) messcrash("Error reading seqFile");
	sscanf(line, "%s", qname);
	
	while ((ch = fgetc(seqfile)) != '\n') {
	    if (isAlnRes(ch) || isAlnNonres(ch))
		array(arr, i++, char) = ch;
	}
	q = messalloc(arrayMax(arr)+1);
	cq = q;
	for (i = 0; i < arrayMax(arr);) *cq++ = arr(arr, i++, char);
	arrayDestroy(arr);
    }
    else {
	fseek(seqfile, 0, SEEK_END);
	q = messalloc(ftell(seqfile)+1);
	cq = q;
	fseek(seqfile, 0, SEEK_SET);
	while (!feof(seqfile))
	{
	    if (!fgets(line, MAXLINE, seqfile)) break;
	    while (strchr(line, '>')) {
		strncpy(qname, line+1, 255); qname[255]=0;
		if ((c = (char *)strchr(qname, ' '))) *c = 0;
		if ((c = (char *)strchr(qname, '\n'))) *c = 0;
		if (!fgets(line, MAXLINE, seqfile)) break;
	    }
	    
	    for (cp = line; *cp; cp++) if (isAlnRes(*cp) || isAlnNonres(*cp)) *cq++ = *cp;
	    *cq = 0;
	}
    }

    /* Determine input format */
    format = EXBLX;
    while ( (ch = fgetc(exblxfile)) == '#') {
	if (!fgets(line, MAXLINE, exblxfile)) break;
	for (c = line; *c; c++) *c = ace_lower(*c);

	if (strstr(line, "seqbl")) format = SEQBL;
	else if (strstr(line, "blastp")) *opts = 'P';
	else if (strstr(line, "tblastn")) *opts = 'T';
	else if (strstr(line, "tblastx")) *opts = 'L';
	else if (strstr(line, "blastn")) *opts = 'N';
	else if (strstr(line, "blastx")) *opts = 'X';
	else if (strstr(line, "hspgaps")) {
	    HSPgaps = 1;
	    opts[7] = 'G';
	}
	else messout("Unrecognised keyword(s) in infile: \"%s\"", line);
    }
    if (!feof(exblxfile)) ungetc(ch, exblxfile);


    /* Read in MSP's */
    /* printf("Reading in the MSP's. Please wait...\n");*/
    DBhead.next = NULL;
    dbhits = &DBhead;
    while (!feof(exblxfile))
      { 
	char * dummy ; /* mieg 2012_01_18 dummy prevent a warning on ubuntu:: ignoring return value of 'fgets', declared with attribute warn_unused_result */

	/* Parse # lines */
	if ((ch = fgetc(exblxfile)) == '#') {
	    dummy = fgets(line, MAXLINE, exblxfile);
	    if (!strncasecmp(line+1, "DESC", 4) &&
		!strncmp(line+6, sname, strlen(sname))) {
		msp->desc = messalloc(256);
		strncpy(msp->desc, line+6+strlen(sname)+1, 255);
		if ((c = strchr(msp->desc, '\n'))) *c = 0;
	    }

	    /* Make sure entire line is read */
	    while (strlen(line) == MAXLINE) 
	      dummy = fgets(line, MAXLINE, exblxfile);
	    if (dummy){} ;
	    continue;
	}
	else if ((unsigned char)ch == (unsigned char)EOF )  break;
				/* EOF checking to make acedb calling work */
	else if (!feof(exblxfile)) ungetc(ch, exblxfile);
	
	msp = (MSP *)messalloc(sizeof(MSP));

	

	if (fscanf(exblxfile, "%d%s%d%d%d%d%s", 
		   &msp->score, msp->frame, &msp->qstart, &msp->qend, 
		   &msp->sstart, &msp->send, sname) 
	    != 7)
	    messcrash("Error parsing MSPs");
	
	if (*opts == 'T') strcpy(msp->frame, "(+1)"); /* MSPcrunch gives sframe for tblastn - restore qframe */

	strncpy(msp->sname, sname, FULLNAMESIZE); 
	msp->next = 0;
	msp->sseq = 0;
	if (msp->score < -3) opts[4] = 'S';

	/* Convert to upper case (necessary?) */
	for (i=0; msp->sname[i]; i++) msp->sname[i] = ace_upper(msp->sname[i]);

	/* Convert subject names to fetchable ones if from NCBI server 

	   Rule 1: If there is a gi, use that.
	   Rule 2: If no gi, use the first and last non-blank field as db:id.
	*/
	if (strchr(msp->sname, '|')) {
	    char *p, *src;

	    src = messalloc(strlen(msp->sname)+1);
	    strcpy(src, msp->sname);
	    
	    p = strtok(src, "|");

	    if (!strcasecmp(p, "GI")) {
		/* Use only GI number */
		p = strtok(0, "|");
		strcpy(msp->sname, "gi");
		strcat(msp->sname, ":");
		strcat(msp->sname, p);
	    }
	    else {
		/* Try to make a proper db:name.  Use last non-blank field */
		char *db=p, *last="";

		p = strtok(0, "|");
		while (p) {
		    if (*p && *p != ' ') last = p;
		    p = strtok(0, "|");
		}
		strcpy(msp->sname, db);
		strcat(msp->sname, ":");
		strcat(msp->sname, last);
	    }

	    messfree(src);
	}

	dbhits->next = msp;
	dbhits = msp;

	qlen = abs(msp->qend - msp->qstart)+1;
	slen = abs(msp->send - msp->sstart)+1;

	if (*opts == ' ') {
	    /* Guess mode from prefix or from coordinates */

	    if (qlen == slen) {
		/* Could be blastp, blastn or tblastx */
		if (strchr(msp->sname, '_')) *opts = 'P';
		else if (strstr(msp->sname, "PIR:")) *opts = 'P';
		else if (strstr(msp->sname, "EM:")) *opts = 'N';
		else *opts = 'N';  /* Could be P or L just as well */
	    }
	    else if (qlen > slen) *opts = 'X';
	    else if (qlen < slen) *opts = 'T';
	}

	if (format == EXBLX) {

	    opts[5] = ' '; /* Don't use full zoom default */

	    /* Skip Description - go to next line */
	    while (fgetc(exblxfile) != '\n');
	}
	else if (format == SEQBL) {
	    int 
		allocated, 
		nchar=0,
		nres=0;
	    
	    /* Copy sequence into MSP */

#ifdef JUNK
    Old_way:
            if (*opts == 'T') {
	      //TBLASTN - Problem: Subject can be reversed.  If reveresed, We can't
              // allocate from start of HSP since we don't know the full length of
              //   * sseq.  - Just pass HSP piece 
		
                
                msp->sseq = messalloc(slen+1);
                c = msp->sseq;
            }
            else if (*opts == 'L') {
                // TBLASTX - same problem as tblastn 
                
                slen = (abs(msp->send - msp->sstart) + 1)/3;

                msp->sseq = messalloc(slen+1);
                c = msp->sseq;
            }
            else {
                if (msp->sstart > msp->send) messcrash("Reversed subjects are not allowed in this mode");
                msp->sseq = messalloc(msp->send+1);
                memset(msp->sseq, '-', msp->send); // Fill up with dashes 
                c = msp->sseq + msp->sstart - 1;
            }
#endif


	    /* New way, can handle MSPs with gaps: */

	    if (*opts == 'L') {
		slen = (abs(msp->send - msp->sstart) + 1)/3;
	    }

	    
	    allocated = slen;
	    c = messalloc(allocated+1);
	    while ((ch = fgetc(exblxfile)) != '\n') {
		if (nchar == allocated) {
		    char *p;
		    allocated *= 2;
		    p = messalloc(allocated+1);
		    strcpy(p, c);
		    messfree(c);
		    c = p;
		}

		if (isAlnRes(ch) || isAlnNonres(ch)) {
		    c[nchar++] = ch;
		    if (isAlnRes(ch)) nres++;
		}
		if (nres > slen) messcrash("Too many residues in MSP \"%s/%d-%d\": expected %d, got %d", 
						msp->sname, msp->sstart, msp->send, slen, nres);
	    }

	    /* printf("%s/%d-%d %s\n", msp->sname, msp->sstart, msp->send, c); */

	    if (*opts == 'T' || *opts == 'L') {
		msp->sseq = c;
	    }
	    else {
		if (msp->sstart > msp->send) messcrash("Reversed subjects are not allowed in this mode");
		msp->sseq = messalloc(msp->sstart+nchar+1);
		memset(msp->sseq, '-', msp->sstart); /* Fill up with dashes */
		strcpy(msp->sseq + msp->sstart - 1, c);
		messfree(c);
	    }
	}

	/* Convert colour strings to ints */
	if (msp->score == -3) {
	    static char *colorNames[NUM_TRUECOLORS] = {
		"WHITE", 
		"BLACK", 
		"LIGHTGRAY", 
		"DARKGRAY",
		"RED", 
		"GREEN", 
		"BLUE",
		"YELLOW", 
		"CYAN", 
		"MAGENTA",
		"LIGHTRED", 
		"LIGHTGREEN", 
		"LIGHTBLUE",
		"DARKRED", 
		"DARKGREEN", 
		"DARKBLUE",
		"PALERED", 
		"PALEGREEN", 
		"PALEBLUE",
		"PALEYELLOW", 
		"PALECYAN", 
		"PALEMAGENTA",
		"BROWN", 
		"ORANGE", 
		"PALEORANGE",
		"PURPLE", 
		"VIOLET", 
		"PALEVIOLET",
		"GRAY", 
		"PALEGRAY",
		"CERISE", 
		"MIDBLUE"
	    };
	    int i, colornr;
	    
	    for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
		if (!strcasecmp(colorNames[i], msp->sname)) colornr = i;
	    if (colornr == -1) {
		messout("Unrecognized color: %s\n", msp->sname);
		colornr = 0;
	    }
	    msp->id = colornr;
	}

	/* Harmonize upper and lower case - 
	   In (t)blastx modes the query is always translated to upper case 
	   if (isupper(*q) || *opts == 'L' || *opts == 'X') {
	   for (c = msp->sseq; *c; c++) *c = ace_upper(*c);
	   }
	   else  {
	   for (c = msp->sseq; *c; c++) *c = ace_lower(*c);
	   }*/
    }

    /* Simulate exons * /
    msp = (MSP *)messalloc(sizeof(MSP));
    msp->sseq = messalloc(100);
    msp->score  = -1;
    msp->id     = 100;
    msp->sstart = 80;
    msp->send   = 159;
    msp->qstart = 80;
    msp->qend   = 159;
    strcpy(msp->frame, "(+1)");
    strcpy(msp->sseq, " exon 2:");
    msp->next = DBhead.next;
    DBhead.next = msp;

    msp = (MSP *)messalloc(sizeof(MSP));
    msp->sseq = messalloc(100);
    msp->score  = -1;
    msp->id     = 100;
    msp->sstart = 5;
    msp->send   = 59;
    msp->qstart = 5;
    msp->qend   = 59;
    strcpy(msp->frame, "(+1)");
    strcpy(msp->sseq, " exon 1:");
    msp->next = DBhead.next;
    DBhead.next = msp;
    */

    if (*opts == 'N') opts[3] = 'R';
    return(blxview (q, qname, dispstart, qoffset, DBhead.next, opts));
}




/* Call an external shell command and print output in a text_scroll window
   from Erik Sonnhammer, 92-01-12.  From w6/external.c
*/

#define MAXLENGTH 1024

void externalCommand (char *command)
{
#if !defined(MACINTOSH)
    FILE *pipe ;
    char text[MAXLENGTH+1], *cp ;
    int line=0, len, maxlen=0;
    static Stack stack ;
    Graph old = graphActive() ;

    stack = stackReCreate (stack, 50) ;

    pipe = popen (command, "r") ;
    while (!feof (pipe))
    { 
	if (!fgets (text, MAXLENGTH, pipe)) break;
	len = strlen (text) ;
	if (len)
	{ 
	    if (text[len-1] == '\n') 
		text[len-1] = '\0';
	    pushText (stack, text) ;
	    line++;
	    if (len > maxlen) maxlen = len;
	}
    }
    pclose (pipe);

    graphCreate (TEXT_SCROLL, command, 0, 0, 0.7, 0.5);
    graphTextFormat(FIXED_WIDTH);
    graphTextBounds (maxlen, line);

    stackCursor(stack, 0) ;
    line = 0 ;
    while ((cp = stackNextText (stack)))
	graphText (cp, 0, line++);
    graphRedraw() ;
    graphActivate (old) ;
#endif
}


/* This was once 'lost'.  Should be in messubs.c
void messStatus(char * text) { }
*/


/* match to template with wildcards.   Authorized wildchars are * ? #
     ? represents any single char
     * represents any set of chars
   case-insensitive.   Example: *Nc*DE# fits abcaNchjDE23

   returns 0 if not found
           1 + pos of first sigificant match (i.e. not a *) if found
*/
/* This was once 'lost'.  is in libace
int pickMatch (char *cp, char *tp)
*/

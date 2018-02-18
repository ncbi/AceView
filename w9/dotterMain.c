/*  Last edited: Sep  7 17:58 1996 (esr) */

/* $Id: dotterMain.c,v 1.3 2012/01/20 01:21:36 mieg Exp $ */

/*  dotterMain.c

*/

#include <stdio.h>
#include <ctype.h>
#include "regular.h"
#include "graph.h"
#include "dotter_.h"

#define MAXLINE 1024

char Seqtype(char *seq)
{
    char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
    char *primenuc    = "ACGTUN";
    char *protonly    = "EFIPQZ";

    /* Simplified version of Sean Eddy's */

    int  pos;
    char c;
    int  po = 0;			/* count of protein-only */
    int  nt = 0;			/* count of t's */
    int  nu = 0;			/* count of u's */
    int  na = 0;			/* count of nucleotides */
    int  aa = 0;			/* count of amino acids */
    int  no = 0;			/* count of others */
  
    /* Look at the first 300 characters
     */
    for (pos = 0; seq[pos] && pos < 300; pos++)
    {
	c = ace_upper(seq[pos]);

	if (strchr(protonly, c)) 
	    po++;
	else if (strchr(primenuc, c)) {
	    na++;
	    if (c == 'T') nt++;
	    else if (c == 'U') nu++;
	}
	else if (strchr(aminos, c)) 
	    aa++;
	else if (isalpha(c)) 
	    no++;
    }

    if (po > 0) return 'P';
    else if (na > aa) return 'N';
    else return 'P';
}


static void strNamecpy(char *dest, char *src, int n)
{
    char *cp;

    while (*src && *src == ' ') src++;
    if ((cp = (char *)strchr(src, ' ')))  *cp = 0;
    if ((cp = (char *)strchr(src, '\n'))) *cp = 0;
    strncpy(dest, src, n);
    dest[n] = 0;
}

int main(int argc, char **argv)
{
    int     
	l, qoffset=0, 
	soffset=0, 
	selfcall=0, 
	qlen = 0, 
	slen = 0,
	revcompq = 0,
	dotterZoom = 0,
	count,
	install = 1,
	pixelFacset = 0;
    char   
	*qseq, *sseq, line[MAXLINE+1],
	qname[FULLNAMESIZE+1], sname[FULLNAMESIZE+1], *cp, *cc, *cq, type = 0, 
	*savefile = 0,
	*loadfile = 0,
	*featurefile = 0,
	*mtxfile = 0,
	opts[] = "D  ",		/* -D (mirror); -w/-c (watson/crick); -H (hsponly) */
	*winsize = 0,
	text[MAXLINE+1];
    FILE   
	*qfile, *sfile;
    float   
	memoryLimit=0;
    MSP   
	*MSPlist=0, *msp = 0 ;

    
    int          optc;
    extern int   optind;
    extern char *optarg;
    char        *optstring="b:cDf:Hil:M:m:p:q:rs:SW:wz:";

    extern char
	*dotterVersion,
	*dotterBinary;

    static char *cc_date = 
#if defined(__DATE__)
    __DATE__
#else
    ""
#endif
    ;

    char *usage;
    static char usageText[] = "\
\n\
 Dotter - Sequence dotplots with image enhancement tools.\n\
\n\
 Reference: Sonnhammer ELL & Durbin R (1995). A dot-matrix program\n\
 with dynamic threshold control suited for genomic DNA and protein\n\
 sequence analysis. Gene 167(2):GC1-10.\n\
\n\
 Usage: dotter [options] <horizontal_sequence> <vertical_sequence>  [X options]\n\
\n\
 Allowed types:                Protein        -      Protein\n\
                               DNA            -      DNA\n\
                               DNA            -      Protein\n\
\n\
 Options:\n\
\n\
 -b <file>      Batch mode, write dotplot to <file>\n\
 -l <file>      Load dotplot from <file>\n\
 -m <float>     Memory usage limit in Mb (default 0.5)\n\
 -z <int>       Set zoom (compression) factor\n\
 -p <int>       Set pixel factor manually (ratio pixelvalue/score)\n\
 -W <int>       Set sliding window size. (K => Karlin/Altschul estimate)\n\
 -M <file>      Read in score matrix from <file> (Blast format; Default: Blosum62).\n\
 -f <file>      Read feature segments from <file>\n\
 -i             Do NOT use installed private colormap, but share with other apps.\n\
 -H             Do not calculate dotplot at startup.\n\
 -r             Reverse and complement horizontal_sequence (DNA vs Protein)\n\
 -D             Don't display mirror image in self comparisons\n\
 -w             For DNA: horizontal_sequence top strand only (Watson)\n\
 -c             For DNA: horizontal_sequence bottom strand only (Crick)\n\
 -q <int>       Horizontal_sequence offset\n\
 -s <int>       Vertical_sequence offset\n\
\n\
 Some X options:\n\
 -acefont <font> Main font.\n\
 -font    <font> Menu font.\n\
\n\
 See http://www.sanger.ac.uk/dotter.html for more info.\n\
\n\
 by Erik Sonnhammer (esr@sanger.ac.uk)\n\
 Version ";
    
    usage = messalloc(strlen(usageText) + strlen(dotterVersion) + strlen(cc_date) + 20);
    sprintf(usage, "%s%s, compiled %s\n", usageText, dotterVersion, cc_date);

    
    while ((optc = getopt(argc, argv, optstring)) != EOF)
	switch (optc) 
	{
	case 'b': 
	    savefile = messalloc(strlen(optarg)+1);
	    strcpy(savefile, optarg);         break;
	case 'c': opts[1] = 'C';              break;
	case 'D': opts[0] = ' ';              break;
	case 'f': 
	    featurefile = messalloc(strlen(optarg)+1);
	    strcpy(featurefile, optarg);      break;
	case 'H': opts[2] = 'H';              break;
	case 'i': install = 0;                break;
	case 'l': 
	    loadfile = messalloc(strlen(optarg)+1);
	    strcpy(loadfile, optarg);         break;
	case 'M': 
	    mtxfile = messalloc(strlen(optarg)+1);
	    strcpy(mtxfile, optarg);          break;
	case 'm': memoryLimit = atof(optarg); break;
	case 'p': pixelFacset = atoi(optarg); break;
	case 'q': qoffset = atoi(optarg);     break;
	case 'r': revcompq = 1;               break;
	case 's': soffset = atoi(optarg);     break;
	case 'S': 
	    selfcall = 1;
	    strncpy(qname, argv[optind], FULLNAMESIZE); qname[FULLNAMESIZE]=0;
	    qlen = atoi(argv[optind+1]);
	    strncpy(sname, argv[optind+2], FULLNAMESIZE); sname[FULLNAMESIZE]=0;
	    slen = atoi(argv[optind+3]);
	    dotterBinary = messalloc(strlen(argv[optind+4])+1);
	    strcpy(dotterBinary, argv[optind+4]);
	                                      break;
	case 'W': 
	    winsize = messalloc(strlen(optarg)+1);
	    strcpy(winsize, optarg);          break;
	case 'w': opts[1] = 'W';              break;
	case 'z': dotterZoom = atoi(optarg);  break;
	default : fatal("Illegal option");
	}

    if (selfcall) /* Dotter calling dotter */
    {
	qseq = (char *)messalloc(qlen+1);
	sseq = (char *)messalloc(slen+1);

	if ((l=fread(qseq, 1, qlen, stdin)) != qlen) 
	    fatal("Only read %d chars to qseq, expected %d", l, qlen);
	qseq[qlen] = 0;

	if ((l=fread(sseq, 1, slen, stdin)) != slen) 
	    fatal("Only read %d chars to sseq, expected %d", l, slen);
	sseq[slen] = 0;

	/* Read in MSPs */
	while (!feof (stdin))
	{ 
	    if (!fgets (text, MAXLINE, stdin) || (unsigned char)*text == (unsigned char)EOF ) break;
	    
	    if (!MSPlist) {
		MSPlist = (MSP *)messalloc(sizeof(MSP));
		msp = MSPlist;
	    }
	    else {
		msp->next = (MSP *)messalloc(sizeof(MSP));
		msp = msp->next;
	    }
	    
	    if (sscanf(text, "%d%d%s%d%d%s%d%d", 
		       &msp->score, 
		       &msp->id, 
		       msp->frame,
		       &msp->qstart,
		       &msp->qend,
		       msp->sname, 
		       &msp->sstart, 
		       &msp->send) != 8)
		messout("Failed to parse line %s\n", text);

	    /*printf("%d %d %s %d %d %s %d %d\n", 
		   msp->score, 
		   msp->id, 
		   msp->frame,
		   msp->qstart,
		   msp->qend,
		   msp->sname, 
		   msp->sstart, 
		   msp->send);*/
	}
	fclose(stdin);
    }
    else
    {
	if (argc - optind < 2) {
	    fprintf(stderr, "%s\n", usage); 
	    exit(1);
	}
	else if(!(qfile = fopen(argv[optind], "r"))) {
	    fprintf(stderr,"Cannot open %s\n", argv[optind]); exit(1);
	}
	fseek(qfile, 0, SEEK_END);
	qlen = ftell(qfile);
	qseq = (char *)messalloc(qlen+1);
	fseek(qfile, 0, SEEK_SET);
	if (!(cp = (char *)strrchr(argv[optind], '/'))) cp = argv[optind]-1;
	strncpy(qname, cp+1, FULLNAMESIZE);
	qname[FULLNAMESIZE] = 0;

	if (!(sfile = fopen(argv[optind+1], "r"))) {
	    fprintf(stderr,"Cannot open %s\n", argv[optind+1]); exit(1);
	}
	fseek(sfile, 0, SEEK_END);
	slen = ftell(sfile);
	sseq = (char *)messalloc(slen+1);
	fseek(sfile, 0, SEEK_SET);
	if (!(cp = strrchr(argv[optind+1], '/'))) cp = argv[optind+1]-1;
	strncpy(sname, cp+1, FULLNAMESIZE);
	qname[FULLNAMESIZE] = 0;
	

	/* Store X options for zooming in */
	{
	    int i, len=0;
	    extern char *Xoptions;

	    for (i = optind+2; i < argc; i++)
		len += strlen(argv[i])+1;

	    Xoptions = messalloc(len+1);

	    for (i = optind+2; i < argc; i++) {
		strcat(Xoptions, argv[i]);
		strcat(Xoptions, " ");
	    }
	}

	/* Read in the sequences */
	l = count = 0;
	cc = qseq;
	while (!feof(qfile))
	{
	    if (!fgets(line, MAXLINE, qfile)) break;

	    /* Name headers */
	    if ((cq = (char *)strchr(line, '>'))) {
		cq++;
		if (++l == 1) {
		    strNamecpy(qname, cq, FULLNAMESIZE);
		}
		else {
		    /* Multiple sequences - add break line */
		    if (!MSPlist) {
			MSPlist = (MSP *)messalloc(sizeof(MSP));
			msp = MSPlist;
		    }
		    else {
			msp = MSPlist;
			while(msp->next) msp = msp->next;
			msp->next = (MSP *)messalloc(sizeof(MSP));
			msp = msp->next;
		    }
		    strNamecpy(msp->sname, cq, FULLNAMESIZE);
		    msp->score = -3;
		    msp->qstart = msp->qend = count;
		    *msp->frame = '1';
		    msp->id = DARKGREEN;
		}
	    }
	    else {
		for (cq = line; *cq; cq++) if (isalpha(*cq)) {
		    *cc++ = *cq;
		    count++;
		}
	    }
	}
	*cc = 0;
	
	l = count = 0;
	cc = sseq;
	while (!feof(sfile))
	{
	    if (!fgets(line, MAXLINE, sfile)) break;

	    /* Name headers */
	    if ((cq = (char *)strchr(line, '>'))) {
		cq++;
		if (++l == 1) {
		    strNamecpy(sname, cq, FULLNAMESIZE);
		}
		else {
		    /* Multiple sequences - add break line */
		    if (!MSPlist) {
			MSPlist = (MSP *)messalloc(sizeof(MSP));
			msp = MSPlist;
		    }
		    else {
			msp = MSPlist;
			while(msp->next) msp = msp->next;
			msp->next = (MSP *)messalloc(sizeof(MSP));
			msp = msp->next;
		    }
		    strNamecpy(msp->sname, cq, FULLNAMESIZE);
		    msp->score = -3;
		    msp->qstart = msp->qend = count;
		    *msp->frame = '2';
		    msp->id = DARKGREEN;
		}
	    }
	    else {
		for (cq = line; *cq; cq++) if (isalpha(*cq)) {
		    *cc++ = *cq;
		    count++;
		}
	    }
	}
	*cc = 0;
    }

    /* Determine sequence types */
    if (Seqtype(qseq) == 'P' && Seqtype(sseq) == 'P') {
	printf("\nDetected sequence types: Protein vs. Protein\n");
	type = 'P';
    }
    else if (Seqtype(qseq) == 'N' && Seqtype(sseq) == 'N') {
	printf("\nDetected sequence types: DNA vs. DNA\n");
	type = 'N';
    }
    else if (Seqtype(qseq) == 'N' && Seqtype(sseq) == 'P') {
	printf("\nDetected sequence types: DNA vs. Protein\n");
	type = 'X';
    }
    else fatal("Illegal sequence types: Protein vs. DNA - turn arguments around!\n\n%s", usage);

    if (revcompq) {
	if (type != 'X')
	    fatal("Revcomp'ing horizontal_sequence only needed in DNA vs. Protein");
	else {
	    cc = messalloc(qlen+1);
	    revcomp(cc, qseq);
	    messfree(qseq);
	    qseq = cc;
	}
    }

    /* Add -install for private colormaps */
    if (install) argvAdd(&argc, &argv, "-install");

    if (!savefile) graphInit(&argc, argv);

    graphLoop(dotter(type, opts, qname, qseq, qoffset, sname, sseq, soffset, 
		     0, 0, savefile, loadfile, mtxfile, featurefile, memoryLimit, dotterZoom, MSPlist, 0, winsize, pixelFacset));
    graphFinish () ;
    return 0 ;
}

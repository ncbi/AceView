/*  Last edited: Aug  9 14:39 1996 (esr) */

/* $Id: blxmain.c,v 1.3 2015/09/25 17:28:30 mieg Exp $ */

/* BLXMAIN - Standalone calling program for blxview. 
 *
 * by Erik Sonnhammer
 *
 */

#include "regular.h"
#include "graph.h"
#include "blxview.h"
#include "dotter_.h"

Graph blxreadhsp(FILE *seqfilename, FILE *exblxfile, char *qname, int a, int b, char *opts, int *argc, char **argv);

int main(int argc, char **argv)
{
    FILE        *seqfile, *exblxfile;
    char         seqfilename[256];
    static char  opts[32]=" MBr Z  ", qname[256];
    extern char *blixemVersion;
    int          dispstart=0, install = 1;

#if !defined(MACINTOSH)
    int          optc;
    extern int   optind;
    extern char *optarg;
    char        *optstring="bhilnpS:s:t";

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
 Blixem - display Blast matches as a multiple alignment.\n\
\n\
 Reference:  Sonnhammer ELL & Durbin R (1994). A workbench for Large Scale\n\
 Sequence Homology Analysis. Comput. Applic. Biosci. 10:301-307.\n\
\n\
 Usage: blixem [options] <sequencefile> <seqblfile> [X options] \n\
\n\
 Both <sequencefile> and <seqblfile> can be substituted by \"-\"\n\
 for reading from stdin (pipe).  If <sequencefile> is piped, the first\n\
 line should contain the sequence name and the second the sequence itself.\n\
\n\
 Options:\n\
 -s <mode>  Sorting mode at startup.\n\
               s = by Score\n\
               i = by Identity\n\
               n = by Name\n\
               p = by Position\n\
 -b         Don't start with Big Picture.\n\
 -S <#>     Start display at position #.\n\
 -i         Do NOT use installed private colormap, but share with other apps\n\
 -h         Help and more options.\n\
\n\
 Some X options:\n\
 -acefont <font> Main font.\n\
 -font    <font> Menu font.\n\
\n\
 To make the seqblfile from blast output, run MSPcrunch with option -q.\n\n\
\n\
 by Erik Sonnhammer (esr@ncbi.nlm.nih.gov)\n\
 Version";
    
    usage = messalloc(strlen(usageText) + strlen(blixemVersion) + strlen(cc_date) + 10);
    sprintf(usage, "%s %s, %s\n", usageText, blixemVersion, cc_date);

    while ((optc = getopt(argc, argv, optstring)) != EOF )
	switch (optc)
	{
	case 'b': opts[2] = 'b';         break;
	case 'h': 
	    fprintf(stderr, "%s\
\n\
 o To pipe MSPcrunch output directly to Blixem, use \"-\"\n\
   as the second parameter ([seqblfile]).  Example:\n\
\n\
   MSPcrunch -q <my.blast_output> | blixem <my.seq> -\n\
\n\
\n\
 o The BLAST program (blastp, blastn, blastx, tblastn, tblastx)\n\
   is automatically detected from the Blast output header by MSPcrunch\n\
   and is passed on to Blixem in the seqbl format (-q).  If\n\
   your output lacks the header, set the program with these options:\n\
\n\
   -p    Use blastp output (blastx is default)\n\
   -n    Use blastn output\n\
   -t    Use tblastn output\n\
   -l    Use tblastx output\n\
\n", usage); 
	    exit(-1);
	                                 break;
	case 'i': install = 0;                break;
	case 'l': opts[0] = 'L';         break;
	case 'n': opts[0] = 'N';         break;
	case 'p': opts[0] = 'P';         break;
	case 'S': 
	    if (!(dispstart = atoi(optarg)))
		messcrash("Bad diplay start position: %s", optarg); 
	    opts[1] = ' ';
	                                 break;
	case 's': 
	    if (*optarg != 's' && *optarg != 'i' && *optarg != 'n' && *optarg != 'p') {
		fprintf(stderr,"Bad sorting mode: %s\n", optarg); 
		exit(1);
	    }
	    opts[6] = *optarg;
	                                 break;
	case 't': opts[0] = 'T';         break;
	default : messcrash ("Illegal option");
	}

    if (argc-optind < 2) {
	fprintf(stderr, "%s\n", usage); 
	exit(1);
    }

    if (!strcmp(argv[optind+1], "-")) {
	exblxfile = stdin;
    }
    else
    {
	if (!(exblxfile = fopen(argv[optind+1], "r"))) {
	    fprintf(stderr,"Cannot open seqbl file %s\n", argv[optind+1]); 
	    exit(1);
	}
    }

    strcpy(seqfilename, argv[optind]);
#else
    printf("Seq file: "); scanf("%s", seqfilename);
    printf("HSP file: "); scanf("%s", exblxfilename);
    if (!(exblxfile = fopen(exblxfilename, "r"))) {
	fprintf(stderr,"Cannot open seqbl file %s\n", exblxfilename); exit(1);
    }
#endif
 
    if (!strcmp(argv[optind], "-")) {
	seqfile = stdin;
    }
    else if(!(seqfile = fopen(seqfilename, "r"))) {
	fprintf(stderr,"Cannot open sequence file %s\n", seqfilename); exit(1);
    }


    /* Add -install for private colormaps */
    if (install) argvAdd(&argc, &argv, "-install");

    graphInit(&argc, argv) ;

    strncpy(qname, seqfilename, 255); qname[255]=0;
    graphLoop(blxreadhsp(seqfile, exblxfile, qname, dispstart, 0, opts, &argc, argv));
    return 0 ;
}

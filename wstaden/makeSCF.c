/*
 * Copyright (c) Medical Research Council 1994-1998. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */


/*
 * makeSCF v3.02, 31/03/98
 *
 * Derived from the older makeSCF; this one has been rewritten to use the new
 * IO libraries and no longer performs quality clipping itself. It also writes
 * in a new format that is more easily compressed.
 */

#include <stdio.h>
#include "Read.h"
#include "traceType.h"
#include "xalloc.h"
#include "compress.h"

/*
 * Add our comments.
 * 1. Work out the comment length - simply add 1K to the current length.
 * 2. Copy the old comments to our new.
 * 3. Replace '\n' with ' '.
 * 4. Add our comments and switch.
 */
void add_comments(Read *r, char *name, int format) {
    Comments *cc;
    int clen;
    char *cp, *src;
    char buf[1024];

    /* 1. */
    if (r->info) {
	clen = strlen(r->info) + 1024;
    } else {
	clen = 1024;
    }

    /* 2. */
    if (NULL == (cc = (char *)malloc(clen)))
	return;

    if (r->info)
	strcpy(cc, r->info);
    else
	*cc = 0;
    
    /* 3. */
    cp = cc;
/*
    while (*cp) {
	if (*cp == '\n')
	    *cp = ' ';
	
	cp++;
    }
    *cp++ = '\n';
*/
  
    /* 4. */
    /* MACH */
    switch (format) {
    case TT_ABI:
	src = "ABI 373A";
	break;
    case TT_ALF:
	src = "Pharmacia A.L.F.";
	break;
    case TT_SCF:
	src = "SCF";
	break;
    case TT_CTF:
	src = "CTF";
	break;
    case TT_EXP:
	src = "Experiment file";
	break;
    case TT_PLN:
	src = "Plain text";
	break;
    default:
	src = "Unknown";
	break;
    }

    sprintf(buf, "CONV=makeSCF V3.00\nMACH=%s\nDATF=%s\nDATN=%s\n",
	    src, trace_type_int2str(format), name);

    strcat(cp, buf);

    if (r->info)
	free(r->info);

    r->info = cc;
}

void scale_trace8(Read *r) {
    double s;
    int i;

    if (r->maxTraceVal <= 255)
	return;
    
    s = ((double)255)/r->maxTraceVal;

    for (i = 0; i < r->NPoints; i++) {
	r->traceA[i] *= s;
	r->traceC[i] *= s;
	r->traceG[i] *= s;
	r->traceT[i] *= s;
    }

    r->maxTraceVal = 255;
}

/*
 * Here we just take the minimum trace value and subtract this from all others.
 * The assumption is that the signal will always be 'base line' on at least
 * one of the four channels.
 */
void subtract_background(Read *r) {
    int i, min;
    for (i = 0; i < r->NPoints; i++) {
	min = 999999;
	if (r->traceA[i] < min) min = r->traceA[i];
	if (r->traceC[i] < min) min = r->traceC[i];
	if (r->traceG[i] < min) min = r->traceG[i];
	if (r->traceT[i] < min) min = r->traceT[i];
	r->traceA[i] -= min;
	r->traceC[i] -= min;
	r->traceG[i] -= min;
	r->traceT[i] -= min;
    }
}

/*
 * Find the average background level of a trace, and subtract this from the
 * peak heights.
 *
 * NB. That method is flawed. For now take the minimum instead of average, but
 * this also has horrid flaws. See above method.
 */
void subtract_background_old(Read *r) {
    int i, j, min, bg, max = 0;
    int win_len = 501, win_len2 = win_len/2;
    int *background;

    if (NULL == (background = (int *)xmalloc((r->NPoints + 2 * win_len)
					     * sizeof(*background))))
	return;

    if (r->NPoints < win_len)
	win_len = r->NPoints;

    /* Find minimum trace levels at each point */
    for (i = 0, j = win_len2; i < r->NPoints; i++, j++) {
	min = INT_MAX;
	if (r->traceA[i] < min)
	    min = r->traceA[i];
	if (r->traceC[i] < min)
	    min = r->traceC[i];
	if (r->traceG[i] < min)
	    min = r->traceG[i];
	if (r->traceT[i] < min)
	    min = r->traceT[i];
	background[j] = min;
    }
    for (i = 0; i < win_len2; i++) {
	background[i] = background[i + win_len2];
	background[i + r->NPoints] = background[i + r->NPoints - win_len2];
    }

    /* Take lowest background over win_len and subtract it */
    for (i = 0; i < r->NPoints; i++) {
	/* Could optimise this considerably */
	bg = INT_MAX;
	for (j = 0; j < win_len; j++) {
	    if (background[i + j] < bg)
		bg = background[i + j];
	}

	r->traceA[i] -= bg;
	r->traceC[i] -= bg;
	r->traceG[i] -= bg;
	r->traceT[i] -= bg;

	if (r->traceA[i] > max) max = r->traceA[i];
	if (r->traceC[i] > max) max = r->traceC[i];
	if (r->traceG[i] > max) max = r->traceG[i];
	if (r->traceT[i] > max) max = r->traceT[i];
    }
    
    r->maxTraceVal = max;

    xfree(background);
}

/*
 * Find the maximum height of traces at the called bases. Use this to clip any
 * other bases.
 */
void reset_max_called_height(Read *r) {
    int i, max = 0;

    /* Find max */
    for (i=0; i < r->NBases; i++) {
	switch(r->base[i]) {
	case 'a':
	case 'A':
	    if (r->traceA[r->basePos[i]] > max)
		max = r->traceA[r->basePos[i]];
	    break;

	case 'c':
	case 'C':
	    if (r->traceC[r->basePos[i]] > max)
		max = r->traceC[r->basePos[i]];
	    break;

	case 'g':
	case 'G':
	    if (r->traceG[r->basePos[i]] > max)
		max = r->traceG[r->basePos[i]];
	    break;

	case 't':
	case 'T':
	    if (r->traceT[r->basePos[i]] > max)
		max = r->traceT[r->basePos[i]];
	    break;
	}
    }

    /* Clip to max */
    for (i = 0; i < r->NPoints; i++) {
	if (r->traceA[i] > max)
	    r->traceA[i] = max;
	if (r->traceC[i] > max)
	    r->traceC[i] = max;
	if (r->traceG[i] > max)
	    r->traceG[i] = max;
	if (r->traceT[i] > max)
	    r->traceT[i] = max;
    }
    if (r->maxTraceVal > max)
	r->maxTraceVal = max;
}

void rescale_heights(Read *r) {
    int win_len = 1000;
    int total = 0;
    int max, max2;
    int i, j;

    if (r->NPoints < 2*win_len + 1)
	return;

    max2 = win_len * r->maxTraceVal;

    for (i = 0; i < win_len; i++) {
	max = 0;
	if (r->traceA[i] > max) max = r->traceA[i];
	if (r->traceC[i] > max) max = r->traceC[i];
	if (r->traceG[i] > max) max = r->traceG[i];
	if (r->traceT[i] > max) max = r->traceT[i];
	total += max;
    }

    for (j = 0; i < r->NPoints; i++, j++) {
	max = 0;
	if (r->traceA[j] > max) max = r->traceA[j];
	if (r->traceC[j] > max) max = r->traceC[j];
	if (r->traceG[j] > max) max = r->traceG[j];
	if (r->traceT[j] > max) max = r->traceT[j];
	total -= max;

	max = 0;
	if (r->traceA[i] > max) max = r->traceA[i];
	if (r->traceC[i] > max) max = r->traceC[i];
	if (r->traceG[i] > max) max = r->traceG[i];
	if (r->traceT[i] > max) max = r->traceT[i];
	total += max;

	r->traceA[j] *= (double)max2 / total;
	r->traceC[j] *= (double)max2 / total;
	r->traceG[j] *= (double)max2 / total;
	r->traceT[j] *= (double)max2 / total;
    }

    for (; j < r->NPoints; j++) {
	r->traceA[j] *= (double)max2 / total;
	r->traceC[j] *= (double)max2 / total;
	r->traceG[j] *= (double)max2 / total;
	r->traceT[j] *= (double)max2 / total;
    }
}

static int convert(char *in, FILE *ofp, char *out, int format, int prec, int comp,
	    int normalise, int outFormat) {
    Read *r;

    if (NULL == (r = read_reading(in, format))) {
	fprintf(stderr, "%s: failed to read\n", in);
	return 1;
    }

    if (normalise) {
	subtract_background(r);
	reset_max_called_height(r);
	rescale_heights(r);
    }

    add_comments(r, in, format);
    if (prec == 1)
	scale_trace8(r);

    if (comp != -1)
	set_compression_method(comp);
    if (0 != (fwrite_reading(ofp, r, outFormat))) {
	fprintf(stderr, "%s: failed to write\n", out);
	read_deallocate(r);
	return 1;
    }

    read_deallocate(r);
    return 0;
}


void usage() {
    fprintf(stderr,
	    "makeSCF [-8] [-2] [-3] [-s] [-compress mode] [-normalise]\n"
	    "       -(abi|alf|scf|ctf|any) input_name [-(output|ctfout) output_name]\n");

    exit(1);
}

int main(int argc, char **argv) {
    int format = TT_ANY, r, prec = 0, version = 3, silent = 0;
    int compress_mode = -1;
    char *inf = NULL;
    char *outf = NULL;
    FILE *ofp = stdout;
    int normalise = 0;
    int outFormat = TT_SCF ;  /* mieg */

    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (strcmp(*argv, "-8") == 0) {
	    prec = 1;
	} else if (strcmp(*argv, "-2") == 0) {
	    version = 2;
	} else if (strcmp(*argv, "-3") == 0) {
	    version = 3;
	} else if (strcmp(*argv, "-normalise") == 0) {
	    normalise = 1;
	} else if (strcmp(*argv, "-s") == 0) {
	    silent = 1;
	} else if (strcasecmp(*argv, "-abi") == 0) {
	    format = TT_ABI;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-alf") == 0) {
	    format = TT_ALF;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-scf") == 0) {
	    format = TT_SCF;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-ctf") == 0) {   /* mieg */
	    format = TT_CTF;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-any") == 0) {
	    format = TT_ANY;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-output") == 0) {
	    outf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-ctfout") == 0) {  /* mieg */
	    outf = *++argv; outFormat = TT_CTF ;
	    argc--;
	} else if (strcasecmp(*argv, "-compress") == 0) {
	    compress_mode = compress_str2int(*++argv);
	    argc--;
	} else {
	    usage();
	}
    }

    if (!inf)
	usage();

    if (!silent) {
	printf("makeSCF v3.02\n");
	printf("Copyright (c) MRC Laboratory of Molecular Biology, 1998. All rights reserved.\n");
    }

    if (outf) {
	ofp = fopen(outf, "wb+");
	if (NULL == ofp) {
	    perror(outf);
	    return 1;
	}
    }

    set_scf_version(version);

    r = convert(inf, ofp, outf, format, prec, compress_mode, normalise, outFormat);
    fclose(ofp);

    return r;
}


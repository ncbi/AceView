/* 
 * Title:       seqIOABI
 * 
 * File: 	seqIOABI.c
 * Purpose:	Reading (not writing) of ABI sequences
 * Last update: Fri Sep 02, 1994
 *
 * Change log:
 * 27/11/90 SD     writeSeqABI() outputs header to sequence file:
 * format: ;{noOfBases}{leftCutOff}{basesWritten}{type}{tracefile}
 * eg:     ;   867    45    383ABI a09b7.s1RES
 * 28.11.90 SD  put undesirables under STLOUIS compilation flag
 * 11.12.90 SD  new static function tail to find file name in path name
 * 02.01.91 SD  Merged with St.L version
 * 15.01.91 SD  New include added (opp.h)
 * 30.07.91 SD  Those ole FWO_ field blues
 * 17.09.91 LFW changed STLOUIS compilation flag to SAVE_EDITS
 *              and AUTO_CLIP
 * 25.10.91 SD  Machine independant I/O...removed BIGENDIAN flag
 * 21.07.92 LFW Added finding of primer position
 * 11.11.92 LFW added section to actually check that the trace it
 *              is trying to open is an ALF file using traceType sub
 * 10.11.92 SD  FWO_ and S/N% interpretation. Comments for information
 *              window.
 * 05-Jul-93 SD Added code to check base positions are in order and adjust
 *              them if they are not
 * 02.09.94 JKB Change to use Read instead of Seq library.
 */


/*
 * In the absense of any better format to store our ABI data in we use
 * the Read structure. Hence this module should be considered part of the
 * Read libary.
 *
 * This library also requires use of the mach-io code for the endian
 * independent machine IO.
 * 
 * The ABI results file is controlled by an index found towards
 * the end --- this is pointed to by a longword found at `IndexPO'.
 * The index consists of a number of entries, each of which is
 * four character label followed by 6 long words. The first of these
 * long words holds a simple count (starting at 1) for those cases
 * where there are multiple entries with the same label. Entries should
 * be found by label (and count), rather than their index position,
 * because entries can be ommited or new ones added. This happens when
 * ABI changes the version of their software and also depending
 * on whether the data was analysed or unalaysed. We do, however,
 * make assumptions about the relative order of entries.
 * 
 * Ideally we would have a separate module which provides a number
 * of functions to extract the data we are interested in, keeping
 * the ABI format well wrapped up and out of harms way.
 * 
 * Note that we are relying on the endian-ness of the machine being
 * appropriate so we can just read long words in as integers. This
 * should be recoded to deal with running on different endians.
 */




/* ---- Imports ---- */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "seqIOABI.h"
#include "Read.h"
#include "abi.h"
#include "fpoint.h"    /* IMPORT: int_to_float */
#include "mach-io.h"
#include "xalloc.h"

/* ---- Constants ---- */

#define BasesPerLine 50 /* For output formatting */

#define baseIndex(B) ((B)=='C'?0:(B)=='A'?1:(B)=='G'?2:3)

static int header_fudge = 0;

int dump_labels(FILE *fp, off_t indexO) {
    off_t entryNum = -1;
    uint_4 entryLabel, entryLw1;

    do {
	entryNum++;

	if (fseek(fp, header_fudge+indexO+(entryNum*IndexEntryLength), 0) != 0)
	    return 0;

	if (!be_read_int_4(fp, &entryLabel))
	    return 0;

	if (!be_read_int_4(fp, &entryLw1))
	    return 0;

	if (entryLabel) {
	    unsigned char c1, c2, c3, c4;

	    c1 = (entryLabel >> 24) & 0xff;
	    c2 = (entryLabel >> 16) & 0xff;
	    c3 = (entryLabel >>  8) & 0xff;
	    c4 = (entryLabel >>  0) & 0xff;

	    if (!isprint(c1))
		break;

	    printf("%c%c%c%c %d\n", c1, c2, c3, c4, entryLw1);
	}
    } while (entryLabel);

    return 0;
}

/*
 * From the ABI results file connected to `fp' whose index starts
 * at byte offset `indexO', return in `val' the `lw'th long word
 * from the `count'th entry labelled `label'.
 * The result is 0 for failure, or index offset for success.
 */
int getABIIndexEntryLW(FILE *fp, off_t indexO,
		       uint_4 label, uint_4 count, int lw,
		       uint_4 *val) {
    off_t entryNum=-1;
    int i;
    uint_4 entryLabel, entryLw1;
    
    do {
	entryNum++;

	if (fseek(fp, header_fudge+indexO+(entryNum*IndexEntryLength), 0) != 0)
	    return 0;

	if (!be_read_int_4(fp, &entryLabel))
	    return 0;

	if (!be_read_int_4(fp, &entryLw1))
	    return 0;
    } while (!(entryLabel == label && entryLw1 == count));
    
    for(i=2; i<=lw; i++)
	if (!be_read_int_4(fp, val))
	    return 0;
    
    return indexO+(entryNum*IndexEntryLength);
}


/*
 * Gets the offset of the ABI index.
 * Returns -1 for failure, 0 for success.
 */
int getABIIndexOffset(FILE *fp, uint_4 *indexO) {
    uint_4 magic;

    /*
     * Initialise header_fudge.
     *
     * This is usually zero, but maybe we've transfered a file in MacBinary
     * format in which case we'll have an extra 128 bytes to add to all
     * our fseeks.
     */
    rewind(fp);
    be_read_int_4(fp, &magic);
    header_fudge = (magic == ABI_MAGIC ? 0 : 128);

    if ((fseek(fp, header_fudge + IndexPO, 0) != 0) ||
	(!be_read_int_4(fp, indexO)))
	return -1;
    else
	return 0;
}

/*
 * Get an "ABI String". These strings are either pointed to by the index
 * offset, or held in the offset itself when the string is <= 4 characters.
 * The first byte of the string determines its length.
 * 'string' is a buffer 256 characters long.
 *
 * Returns -1 for failure, string length for success.
 */
int getABIString(FILE *fp, off_t indexO, uint_4 label, uint_4 count,
		 char *string) {
    uint_4 off;
    uint_4 len;

    if ((off = getABIIndexEntryLW(fp, indexO, label, count, 4, &len))) {
	uint_1 len2;

	if (!len)
	    return 0;

	/* Determine offset */
	if (len <= 4)
	    off += 20;
	else
	    getABIIndexEntryLW(fp, indexO, label, count, 5, &off);

	/* Read length byte */
	fseek(fp, header_fudge + off, 0);
	be_read_int_1(fp, &len2);

	/* Read data (max 255 bytes) */
	fread(string, len2, 1, fp);
	string[len2] = 0;

	return len2;
    } else
	return -1;
}


/*
 * Read the ABI format sequence from FILE *fp into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s. In this respect we
 * are adhering (more or less) to the CSET_DEFAULT uncertainty code set.
 * 
 * Returns:
 *   Read *	- Success, the Read structure read.
 *   NULLRead	- Failure.
 */
Read *fread_abi(FILE *fp) {
    Read *read = NULLRead;
    int i;
    float fspacing;		/* average base spacing */
    uint_4 numPoints, numBases;
    uint_4 signalO;

    uint_4 fwo_;     /* base -> lane mapping */
    uint_4 indexO;   /* File offset where the index is */
    uint_4 baseO;    /* File offset where the bases are stored */
    uint_4 basePosO; /* File offset where the base positions are stored */
    uint_4 dataCO;   /* File offset where the C trace is stored */
    uint_4 dataAO;   /* File offset where the A trace is stored */
    uint_4 dataGO;   /* File offset where the G trace is stored */
    uint_4 dataTO;   /* File offset where the T trace is stored */
    uint_4 MCHN_O;   /* File offset where the machine name is stored */
    uint_4 PDMF_O;   /* File offset where the dye primer guff is stored */
    uint_4 SMPL_O;   /* File offset where the sample name is stored */
    uint_4 SMPL_len; /* Length of sample name is stored */
    

    

    /* Get the index offset */
    if (-1 == getABIIndexOffset(fp, &indexO))
	goto bail_out;
    
    /* Get the number of points */
    if (!getABIIndexEntryLW(fp,(off_t)indexO,DataEntryLabel,9,3,&numPoints))
	goto bail_out;	
    
    /* Get the number of bases */
    if (!getABIIndexEntryLW(fp,(off_t)indexO,BaseEntryLabel,1,3,&numBases))
	goto bail_out;	

    
    /* Allocate the sequence */
    if (NULLRead == (read = read_allocate(numPoints, numBases)))
	goto bail_out;	
    
    /*
     * The order of the DATA fields is determined by the field FWO_
     * Juggle around with data pointers to get it right
     */
    {
	uint_4 *dataxO[4];
	
	dataxO[0] = &dataCO;
	dataxO[1] = &dataAO;
	dataxO[2] = &dataGO;
	dataxO[3] = &dataTO;
	
	/* Get the Freak World Out (FWO?) field ... */
	if (!getABIIndexEntryLW(fp,(off_t)indexO,FWO_Label,1,5,&fwo_))
	    goto bail_out;

	/*Get the positions of the four traces */
	if (!(getABIIndexEntryLW(fp, (off_t)indexO, DataEntryLabel,  9, 5,
				 dataxO[baseIndex((char)(fwo_>>24&255))]) &&
	      getABIIndexEntryLW(fp, (off_t)indexO, DataEntryLabel, 10, 5,
				 dataxO[baseIndex((char)(fwo_>>16&255))]) &&
	      getABIIndexEntryLW(fp, (off_t)indexO, DataEntryLabel, 11, 5,
				 dataxO[baseIndex((char)(fwo_>>8&255))]) &&
	      getABIIndexEntryLW(fp, (off_t)indexO, DataEntryLabel, 12, 5,
				 dataxO[baseIndex((char)(fwo_&255))]))) {
	    goto bail_out;
	}
    }
    
    
    /*************************************************************
     * Read the traces and bases information
     *************************************************************/

    /* Read in the C trace */
    if (fseek(fp, header_fudge + (off_t)dataCO, 0) == -1)
	goto bail_out;

    for (i = 0; i < (read->NPoints); i++) {
	if (!be_read_int_2(fp, &(read->traceC[i])))
	    goto bail_out;
    }
    
    /* Read in the A trace */
    if (fseek(fp, header_fudge + (off_t)dataAO, 0) == -1)
	goto bail_out;

    for (i = 0; i < (read->NPoints); i++) {
	if (!be_read_int_2(fp, &(read->traceA[i])))
	    goto bail_out;
    }
    
    /* Read in the G trace */
    if (fseek(fp, header_fudge + (off_t)dataGO, 0) == -1)
	goto bail_out;
    
    for (i = 0; i < (read->NPoints); i++) {
	if (!be_read_int_2(fp, &(read->traceG[i])))
	    goto bail_out;
    }
    
    /* Read in the T trace */
    if (fseek(fp, header_fudge + (off_t)dataTO, 0) == -1)
	goto bail_out;

    for (i = 0; i < (read->NPoints); i++) {
	if (!be_read_int_2(fp, &(read->traceT[i])))
	    goto bail_out;
    }
    
    
    /* Compute highest trace peak */
    for (i=0; i < read->NPoints; i++) {
	if (read->maxTraceVal < read->traceA[i])
	    read->maxTraceVal = read->traceA[i];
	if (read->maxTraceVal < read->traceC[i])
	    read->maxTraceVal = read->traceC[i];
	if (read->maxTraceVal < read->traceG[i])
	    read->maxTraceVal = read->traceG[i];
	if (read->maxTraceVal < read->traceT[i])
	    read->maxTraceVal = read->traceT[i];
    }

    
    /* Read in the bases */
    if (!(getABIIndexEntryLW(fp, (off_t)indexO, BaseEntryLabel, 1, 5, &baseO)
	  && (fseek(fp, header_fudge + (off_t)baseO, 0) == 0) ))
	goto bail_out;

    for (i = 0; i < (read->NBases); i++) {
	int ch;
	
	if ((ch = fgetc(fp)) == EOF)
	    goto bail_out;

	read->base[i] = (ch == 'N') ? '-' : (char)ch;
	read->prob_A[i] = 0;
	read->prob_C[i] = 0;
	read->prob_G[i] = 0;
	read->prob_T[i] = 0;
    }
    read->base[i] = 0;
    
    
    /* Read in the base positions */
    if (!(getABIIndexEntryLW(fp, (off_t)indexO, BasePosEntryLabel, 1, 5,
			     &basePosO) &&
          (fseek(fp, header_fudge + (off_t)basePosO, 0) == 0)))
	goto bail_out;
    
    for (i = 0; i < (read->NBases); i++) {
	if (!be_read_int_2(fp, (uint_2 *)&read->basePos[i]))
	    goto bail_out;
    }
    
    
    
    /*************************************************************
     * Gather useful information - the comments field
     *************************************************************/
    {
	char comment[8192];
	char line[128], *p;
	int_4 spacing;
	int_4 ppos;
	int_4 commlen;
	uint_4 commentO;
	
	*comment = '\0';
	
	/* The ABI comments */
	if ((commentO = getABIIndexEntryLW(fp, (off_t)indexO, CMNTLabel, 1, 4,
					  (uint_4 *)&commlen))) {
	    uint_1 clen;
	    char commstr[256], *commstrp;

	    if (commlen > 4) {
		getABIIndexEntryLW(fp, (off_t)indexO, CMNTLabel, 1, 5,
				   (uint_4 *)&commentO);
	    } else {
		commentO += 20;
	    }

	    fseek(fp, header_fudge + commentO, 0);
	    be_read_int_1(fp, &clen);
	    
	    fread(commstr, clen, 1, fp);
	    commstr[clen] = 0;
	    commstrp = commstr;

	    if (clen) {
		do {
		    char line[300];
		    
		    if ((p = strchr(commstrp, '\n')))
			*p++ = 0;
		    
		    sprintf(line, "COMM=%s\n", commstrp);
		    strcat(comment, line);
		} while((commstrp = p)) ;
	    }
	}
	
	/* Get Signal Strength Offset */
	if (getABIIndexEntryLW(fp, (off_t)indexO, SignalEntryLabel, 1, 5,
			       &signalO)) {
	    int_2 C,A,G,T;
	    int_2 *base[4];
	    base[0] = &C;
	    base[1] = &A;
	    base[2] = &G;
	    base[3] = &T;

	    if (fseek(fp, header_fudge + (off_t)signalO, 0) != -1 &&
		be_read_int_2(fp, (uint_2 *)
			      base[baseIndex((char)(fwo_>>24&255))]) &&
		be_read_int_2(fp, (uint_2 *)
			      base[baseIndex((char)(fwo_>>16&255))]) &&
		be_read_int_2(fp, (uint_2 *)
			      base[baseIndex((char)(fwo_>>8&255))]) &&
		be_read_int_2(fp, (uint_2 *)
			      base[baseIndex((char)(fwo_&255))])) {
		sprintf(line, "SIGN=A=%d,C=%d,G=%d,T=%d\n",
			A, C, G, T);
		strcat(comment, line);
	    }
	}

	/* Get the spacing.. it's a float but don't worry yet */
	if (getABIIndexEntryLW(fp, (off_t)indexO, SpacingEntryLabel, 1, 5,
			       (uint_4*)&spacing)) {
	    fspacing = int_to_float(spacing);
	    sprintf(line, "SPAC=%-6.2f\n", fspacing);
	    strcat(comment, line);
	} else {
	    fspacing = (float) read->NPoints / (float) read->NBases;
	}

	
	/* Get primer position */
	if (getABIIndexEntryLW(fp, (off_t)indexO, PPOSLabel, 1, 5,
			       (uint_4 *)&ppos)) {
	    /* ppos stores in MBShort of pointer */
	    sprintf(line, "PRIM=%d\n", (ppos>>16));
	    strcat(comment, line);
	}

	/* Get Machine Name Offset */
	if (getABIIndexEntryLW(fp, (off_t)indexO, MCHNLabel, 1, 5, &MCHN_O)) {
	    if (fseek(fp, header_fudge + (off_t)MCHN_O, 0) != -1) {
		unsigned char l;
		char buffer[257];

		/* first byte is a length */
		fread(&l, sizeof(char), 1, fp);
		fread(buffer, l, 1, fp);
		buffer[l] = 0;
		if (strchr(buffer, '\n'))
		    *strchr(buffer, '\n') = ' ';
		sprintf(line, "MACH=%.*s\n", l, buffer);
		strcat(comment, line);
	    } 
	}
	
	/* Get Dye Primer Offset */
	if (getABIIndexEntryLW(fp, (off_t)indexO, PDMFLabel, 1, 5, &PDMF_O)) {
	    if (fseek(fp, header_fudge + (off_t)PDMF_O, 0) != -1) {
		unsigned char l;
		char buffer[256];
		
		/* first byte is a length */
		fread(&l, sizeof(char), 1, fp);
		fread(buffer, l, 1, fp);
		sprintf(line, "DYEP=%.*s\n", l, buffer);
		strcat(comment, line);
	    }
	}

	/* Get Sample Name Offset */
	/*
         * When the sample name length is less than or equal to four, the
	 * string is packed into the SMPL_O field itself. Otherwise SMPL_O
	 * is as we expect (the offset of the real name).
	 */
	if ((SMPL_O = getABIIndexEntryLW(fp, (off_t)indexO, SMPLLabel, 1, 4,
					&SMPL_len))) {
	    if (SMPL_len > 4) {
		getABIIndexEntryLW(fp,(off_t)indexO, SMPLLabel, 1, 5, &SMPL_O);
	    } else {
		SMPL_O += 20;
	    }
		
	    if (fseek(fp, header_fudge + (off_t)SMPL_O, 0) != -1) {
		unsigned char l;
		char buffer[256];
		/* first byte is a length */

		fread(&l, sizeof(char), 1, fp);
		fread(buffer, l, 1, fp);
		sprintf(line, "NAME=%.*s\n", l, buffer);
		strcat(comment, line);
	    }
	}
	
	/* dumplicate string and set info */
	{
	    char *s = (char *)xmalloc(strlen(comment)+1);
	    strcpy(s,comment);
	    read->info = s;
	}
    }



    /*************************************************************
     * Check base positions are in order
     *************************************************************/
    {
	float pos;
	int start;

	for (i = 1; i < read->NBases; ) {
	    if (read->basePos[i] < read->basePos[i-1]) {
		fprintf(stderr,"ted: Base positions are not in order. Fixing\n");

		/* pass 1 - find end of region */
		start = i - 1;
		pos = (float) read->basePos[i-1] + fspacing;
		for(;i < read->NBases && (int)read->basePos[i] < pos;i++) {
		    pos += fspacing;
		}

		/* calculate average base spacing */
		if (i < read->NBases )
		    fspacing = ((float) read->basePos[i] -
				(float) read->basePos[start]) /
				    (float)(i - start);

		/* pass 2 - adjust */
		i = start + 1;
		pos = (float) read->basePos[i-1] + fspacing;
		for(;i < read->NBases && (int)read->basePos[i] < pos;i++) {
		    read->basePos[i] = (int) pos;
		    pos += fspacing;
		}

	    } else {
		i++;
	    }
	}
    }

    
    /* SUCCESS */

    read->format = TT_ABI;
    return(read);

    /* FAILURE */
 bail_out:
    if (read)
	read_deallocate(read);

    return NULLRead;
}

/*
 * Read the ABI format sequence from file 'fn' into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s. In this respect we
 * are adhering (more or less) to the CSET_DEFAULT uncertainty code set.
 * 
 * Returns:
 *   Read *	- Success, the Read structure read.
 *   NULLRead	- Failure.
 */
Read *read_abi(char *fn) {
    Read *read;
    FILE *fp;

    /* Open file */
    if ((fp = fopen(fn, "rb")) == NULL)
	return NULLRead;

    read = fread_abi(fp);
    fclose(fp);

    if (read && (read->trace_name = (char *)xmalloc(strlen(fn)+1)))
	strcpy(read->trace_name, fn);

    return read;
}
    
/*
 * Write to an ABI file - unsupported.
 */
/* ARGSUSED */
int write_abi(char *fn, Read *read) {
    fprintf(stderr, "ABI write support is unavailable\n");
    return -1;
}

/*
 * Write to an ABI file - unsupported.
 */
/* ARGSUSED */
int fwrite_abi(FILE *fp, Read *read) {
    fprintf(stderr, "ABI write support is unavailable\n");
    return -1;
}


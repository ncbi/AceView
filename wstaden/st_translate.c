/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
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
 * File:	translate.c
 * Purpose:	Translates between different reading structures.
 * Last update:	01/09/94
 */

#include <stdio.h>
#include "misc.h"
#include "scf.h"
#include "Read.h"
#include "expFileIO.h"
#include "traceType.h"
#include "stadentranslate.h"
#include "xalloc.h"

static int topos(Read *r, int ind);

/*
 * Translates an Scf structure into a Read structure.
 * The Scf structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Read structure upon success.
 *    NULLRead upon failure.
 */
Read *scf2read(Scf *scf) {
    Read *read;
    register int i, i_end;
    TRACE max_val = 0;

    /* allocate */
    read = read_allocate(scf->header.samples, scf->header.bases);
    if (NULLRead == read)
	return NULLRead;

    /* copy the samples */
    i_end = scf->header.samples;
    read->NPoints = i_end;

    if (scf->header.sample_size == 1) {
	for (i = 0; i < i_end; i++) {
	    read->traceA[i] = scf->samples.samples1[i].sample_A;
	    read->traceC[i] = scf->samples.samples1[i].sample_C;
	    read->traceG[i] = scf->samples.samples1[i].sample_G;
	    read->traceT[i] = scf->samples.samples1[i].sample_T;

	    if (read->traceA[i] > max_val) max_val = read->traceA[i];
	    if (read->traceC[i] > max_val) max_val = read->traceC[i];
	    if (read->traceG[i] > max_val) max_val = read->traceG[i];
	    if (read->traceT[i] > max_val) max_val = read->traceT[i];
	}
    } else { /* sample_size == 2 */
	for (i = 0; i < i_end; i++) {
	    read->traceA[i] = scf->samples.samples2[i].sample_A;
	    read->traceC[i] = scf->samples.samples2[i].sample_C;
	    read->traceG[i] = scf->samples.samples2[i].sample_G;
	    read->traceT[i] = scf->samples.samples2[i].sample_T;

	    if (read->traceA[i] > max_val) max_val = read->traceA[i];
	    if (read->traceC[i] > max_val) max_val = read->traceC[i];
	    if (read->traceG[i] > max_val) max_val = read->traceG[i];
	    if (read->traceT[i] > max_val) max_val = read->traceT[i];
	}
    }

    read->maxTraceVal = max_val;
    
    /* copy the bases */
    i_end = scf->header.bases;
    read->NBases = i_end;

    for (i = 0; i < i_end; i++) {
	read->basePos[i] = scf->bases[i].peak_index;
	read->prob_A[i]  = scf->bases[i].prob_A;
	read->prob_C[i]  = scf->bases[i].prob_C;
	read->prob_G[i]  = scf->bases[i].prob_G;
	read->prob_T[i]  = scf->bases[i].prob_T;
	read->base[i]    = scf->bases[i].base;
    }
    read->base[i] = 0;
    
    /* allocate and copy the comments */
    if (scf->header.comments_size > 0 && scf->comments) {
	read->info = (char *)xmalloc(scf->header.comments_size+1);
	if (NULL == read->info) {
	    read_deallocate(read);
	    return NULLRead;
	}

	memcpy(read->info, scf->comments, scf->header.comments_size);
	read->info[scf->header.comments_size] = '\0';
    }

    /* other bits and pieces */
    read->leftCutoff = scf->header.bases_left_clip;
    read->rightCutoff = read->NBases - scf->header.bases_right_clip + 1;
    read->format = TT_SCF;

    return read;
}

/*
 * Translates a Read structure into a Scf structure.
 * The Read structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Scf structure upon success.
 *    NULL upon failure.
 */
Scf *read2scf(Read *read) {
    Scf *scf;
    register int i, i_end;
    int sample_size;

    /* allocate */
    sample_size = read->maxTraceVal >= 0x100 ? 2 : 1;
    scf = scf_allocate(read->NPoints, sample_size, read->NBases, 0, 0);
    if (NULL == scf)
	return NULL;

    /* copy the samples */
    i_end = read->NPoints;
    scf->header.samples = i_end;

    if (sample_size == 1) {
	scf->header.sample_size = 1;
	for (i = 0; i < i_end; i++) {
	    scf->samples.samples1[i].sample_A = read->traceA[i];
	    scf->samples.samples1[i].sample_C = read->traceC[i];
	    scf->samples.samples1[i].sample_G = read->traceG[i];
	    scf->samples.samples1[i].sample_T = read->traceT[i];
	}
    } else {
	scf->header.sample_size = 2;
	for (i = 0; i < i_end; i++) {
	    scf->samples.samples2[i].sample_A = read->traceA[i];
	    scf->samples.samples2[i].sample_C = read->traceC[i];
	    scf->samples.samples2[i].sample_G = read->traceG[i];
	    scf->samples.samples2[i].sample_T = read->traceT[i];
	}
    }

    /* copy the bases */    
    i_end = read->NBases;
    scf->header.bases = i_end;
    
    for (i = 0; i < i_end; i++) {
	scf->bases[i].peak_index = read->basePos[i];
	scf->bases[i].prob_A     = read->prob_A[i];
	scf->bases[i].prob_C     = read->prob_C[i];
	scf->bases[i].prob_G     = read->prob_G[i];
	scf->bases[i].prob_T     = read->prob_T[i];
	scf->bases[i].base       = read->base[i];
    }

    /* allocate and copy the comments */
    if (read->info) {
	scf->header.comments_size = strlen(read->info) + 1;
	scf->comments = (char *)xmalloc(scf->header.comments_size);
	if (NULL == scf->comments) {
	    scf_deallocate(scf);
	    return NULL;
	}

	memcpy(scf->comments, read->info, scf->header.comments_size - 1);

	/* just to make sure */
	scf->comments[scf->header.comments_size-1] = '\0';
    }

    /* other bits and pieces */
    scf->header.bases_left_clip = read->leftCutoff;
    scf->header.bases_right_clip = read->NBases - read->rightCutoff + 1;
    scf->header.code_set = CSET_DEFAULT;
    memcpy(scf->header.version, scf_version_float2str(SCF_VERSION), 4);

    return scf;
}

#define extend(e, entry, len) \
do { \
    (void)ArrayRef(e->entries[entry],e->Nentries[entry]++); \
    if (NULL == (exp_get_entry(e, entry) = (char *)xmalloc(len))) \
	return NULL; \
} while (0)

/*
 * Translates a Read structure and an Experiment file.
 * The Read structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Exp_info structure upon success.
 *    NULL upon failure.
 */
Exp_info *read2exp(Read *read, char *EN) {
    Exp_info *e;
    char *t = trace_type_int2str(read->format), *p;
    int l = strlen(EN)+1;

    if (NULL == (e = exp_create_info()))
	return NULL;

    /* Copy original exp file if present */
    if (read->orig_trace && read->orig_trace_format == TT_EXP) {
	int i, j, k;
	Exp_info *re = (Exp_info *)read->orig_trace;

	for (i = 0; i < MAXIMUM_EFLTS; i++) {
	    if (EFLT_SQ == i ||
		EFLT_QL == i ||
		EFLT_QR == i)
		continue;

	    if (0 == (k = exp_Nentries(re, i)))
		continue;

	    e->Nentries[i] = k;	    
	    ArrayRef(e->entries[i], e->Nentries[i]);
	    for (j = 0; j < k; j++) {
		arr(char *, e->entries[i], j) =
		    strdup(arr(char *, re->entries[i], j));
	    }
	}

    /* Otherwise create our EN, ID, LN and LT lines */
    } else {
	/* Entry name and ID lines */
	if ((p = strrchr(EN, '/')))
	    EN = p+1;
	extend(e, EFLT_EN, l);
	sprintf(exp_get_entry(e, EFLT_EN), "%s", EN);
	extend(e, EFLT_ID, l);
	sprintf(exp_get_entry(e, EFLT_ID), "%s", EN);

	/* Trace file & type */
	if (read->trace_name) {
	    char *cp;
	    if ((cp = strrchr(read->trace_name, '/')))
		cp++;
	    else
		cp = read->trace_name;
	    extend(e, EFLT_LN, strlen(cp)+1);
	    strcpy(exp_get_entry(e, EFLT_LN), cp);
	}

	if (read->format != TT_ANY) {
	    extend(e, EFLT_LT, strlen(t)+1);
	    strcpy(exp_get_entry(e, EFLT_LT), t);
	}
    }

    /* Output SQ, QL and QR lines */

    /* Cutoffs */
    if (read->leftCutoff) {
	extend(e, EFLT_QL, 5);
	sprintf(exp_get_entry(e, EFLT_QL), "%d", read->leftCutoff);
    }

    if (read->rightCutoff) {
	extend(e, EFLT_QR, 5);
	sprintf(exp_get_entry(e, EFLT_QR), "%d", read->rightCutoff);
    }

    /* Bases */
    extend(e, EFLT_SQ, read->NBases+1);
    strncpy(exp_get_entry(e, EFLT_SQ), read->base, read->NBases);
    exp_get_entry(e, EFLT_SQ)[read->NBases] = 0;

    return e;
}


/*
 * Controls the use of the SQ and ON lines when loading an experiment file.
 * The default (value&1 == 1) is to load these into the Read structure.
 * With value&1 == 0 we load the sequence directly from the trace file
 * (LT line).
 * value&2 controls whether to use the SL/SR fields when setting the cutoff.
 * value&2 == 0 implies to do so, and value&2 == 2 implies to not.
 *
 * The default use is to use the SQ and ON lines. Returns the old value.
 */
static int use_experiment_sequence = 1;
int read_experiment_redirect(int value) {
    int old = use_experiment_sequence;

    use_experiment_sequence = value;
    return old;
}

/*
 * Using RAWDATA to find a trace file
 */
static int find_trace_file(char *trace, char *path)
{
    if ( ! file_exists(trace) ) {
        char *s;
        char *rawData;
        /* try in rawData */
        /*
         * get environment details
         */
        rawData = (char *)getenv ("RAWDATA");
        if (rawData == NULL) return 1;

        if ((s=findfile(trace,rawData))==NULL) return 1;
        strcpy (path,s);
    } else
        /* ok */
        strcpy (path,trace);
    return 0;
}

/*
 * Translates an experiment file to a Read structure.
 * The Exp_info structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Read structure upon success.
 *    NULLRead upon failure.
 */
Read *exp2read(Exp_info *e) {
    Read *r;
    int q, s, ttype, err = 0;
    char *str;
    int use_exp = use_experiment_sequence;
    char path[1024];

    if (!exp_Nentries(e, EFLT_LN)) {
	err = 1;
    } else {
	/* Read the trace component of the experiment file */
	ttype = exp_Nentries(e,EFLT_LT)
	    ? trace_type_str2int(exp_get_entry(e, EFLT_LT)) : TT_ANY;

	err = find_trace_file(exp_get_entry(e, EFLT_LN), path);

	if (0 == err && (NULLRead == (r = read_reading(path, ttype))))
	    err = 1;
    }

    if (err) {
	use_exp = 1;
	r = read_allocate(0, 1);
    }
    
    /* Set the left cutoff (QL / SL) */
    q=-1;
    if (exp_Nentries(e, EFLT_QL))
	q = atoi(exp_get_entry(e, EFLT_QL));
    if ((use_exp&2) != 2) {
	s=-1;
	if (exp_Nentries(e, EFLT_SL))
	    s = atoi(exp_get_entry(e, EFLT_SL));
	if (q != -1 || s != -1)
	    r->leftCutoff = MAX(q, s);
    } else {
	if (q != -1)
	    r->leftCutoff = q;
    }

    /* Set the right cutoff (QR / SR) */
    q = INT_MAX;
    if (exp_Nentries(e, EFLT_QR))
	q = atoi(exp_get_entry(e, EFLT_QR));
    if ((use_exp&2) != 2) {
	s = INT_MAX;
	if (exp_Nentries(e, EFLT_SR))
	    s = atoi(exp_get_entry(e, EFLT_SR));
	if (q != INT_MAX || s != INT_MAX)
	    r->rightCutoff = MIN(q, s);
    } else {
	r->rightCutoff = q;
    }

    /* Bases and base positions, if desired */
    if (use_exp&1) {
	if (exp_Nentries(e, EFLT_SQ) && (str = exp_get_entry(e, EFLT_SQ))) {
	    int slen = strlen(str);
	    
	    if (NULL == (r->base = (char *)xrealloc(r->base, slen+1)))
		return NULLRead;
	    
	    strcpy(r->base, str);
	    r->NBases = slen;

	    if (exp_Nentries(e, EFLT_ON) &&
		(str = exp_get_entry(e, EFLT_ON))) {
		int_2 *opos;
		
		if (NULL == (r->basePos = (uint_2 *)xrealloc(r->basePos,
							     slen*2)))
		    return NULLRead;
		
		if ((opos = (int_2 *)xcalloc(slen, 2))) {
		    str2opos(opos, slen, str);
		    
		    read_update_opos(r, opos, slen, r->basePos,
				     NULL, NULL, NULL, NULL);
		}
	    }
	}
    }

    r->format = TT_EXP;
    r->orig_trace = e;
    r->orig_trace_format = TT_EXP;

    return r;
}

/*
 * Given an original positions array, update the 'basePos' and confidence
 * arrays to reflect the new sample numbers (effectively the "X coordinates").
 * The arrays to update (basePos, prob_[ACGT]) are assumed to have already
 * been allocated to the size r->NBases. The opos array is assumed to be
 * of size 'slen'.
 *
 * Specifying any of the arrays to be updated as NULL simply ignores
 * any updates required.
 * 
 * The updates are computed using 'r' and 'opos' ('slen' bases) and written
 * to the update arrays. It's possible to send over the real arrays
 * (r->basePos).
 */
void read_update_opos(Read *r, int_2 *opos, int slen, uint_2 *basePos,
		      char *prob_A, char *prob_C, char *prob_G, char *prob_T) {
    uint_2 *bp = r->basePos;
    int i, diff;
    int olen = r->NBases;

    diff = 1;
    for (i=0; i < slen; i++) {
	int d;
	
	/* Insertion */
	if (opos[i] == 0) {
	    continue;
	}
	
	/* Deletion */
	if ((d = opos[i] - (i+diff)) && (olen - (i+d)) > 0) {
	    if (basePos) memmove(&basePos[i], &basePos[i+d], (olen-(i+d))*2);
	    if (prob_A)  memmove(&prob_A[i],  &prob_A[i+d],  olen-(i+d));
	    if (prob_C)  memmove(&prob_C[i],  &prob_C[i+d],  olen-(i+d));
	    if (prob_G)  memmove(&prob_G[i],  &prob_G[i+d],  olen-(i+d));
	    if (prob_T)  memmove(&prob_T[i],  &prob_T[i+d],  olen-(i+d));
	    diff += d;
	}
    }
    for (i = 0; i < slen; i++) {
	if (opos[i] == 0) {
	    if (basePos) basePos[i] = 0;
	    if (prob_A)  prob_A[i] = 0;
	    if (prob_C)  prob_C[i] = 0;
	    if (prob_G)  prob_G[i] = 0;
	    if (prob_T)  prob_T[i] = 0;
	}
    }

    r->basePos = basePos;
    for (i = 0; i < slen; i++) {
	if (basePos) basePos[i] = topos(r, i);
    }
    r->basePos = bp;
}


/*
 * Find the position of an edited base
 */
static int topos(Read *r, int ind) {
    int back, forw;
    int backp, forwp;

    if (ind < 0)
        return 0;
    if (ind >= r->NBases)
        ind = r->NBases-1;

    /* The simple case */
    if (r->basePos[ind])
        return r->basePos[ind];

    /* Scan in both directions for last non zero position */
    for (back = ind-1; back >= 0 && r->basePos[back] == 0; back--)
        ;
    if (back < 0)
        back = 0;

    for (forw = ind+1; forw < r->NBases && r->basePos[forw] == 0; forw++)
        ;

    if (r->basePos[forw])
        forwp = r->basePos[forw];
    else
        forwp = r->NPoints;

    if (r->basePos[back])
        backp = r->basePos[back];
    else
        backp = forwp / 4;

    /* Interpolate */
    return (forwp - backp) * (ind - back) / (forw - back) + backp;
}

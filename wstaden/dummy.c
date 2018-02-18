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

#include <stdio.h>
#include "Read.h"

/*
 * Keep our compilers happy for our unavailable library formats
 */

/* --------------------------------- ABI --------------------------------- */

#ifdef NO_ABI
Read *read_abi(char *fn) {
    fprintf(stderr, "ABI read support is unavailable\n");
    return NULLRead;
}

int write_abi(char *fn, Read *read) {
    fprintf(stderr, "ABI write support is unavailable\n");
    return -1;
}

Read *fread_abi(FILE *fp) {
    fprintf(stderr, "ABI read support is unavailable\n");
    return NULLRead;
}

int fwrite_abi(FILE *fp, Read *read) {
    fprintf(stderr, "ABI write support is unavailable\n");
    return -1;
}
#endif

/* --------------------------------- ALF --------------------------------- */
#ifdef NO_ALF
Read *read_alf(char *fn) {
    fprintf(stderr, "ALF read support is unavailable\n");
    return NULLRead;
}

int write_alf(char *fn, Read *read) {
    fprintf(stderr, "ALF write support is unavailable\n");
    return -1;
}

Read *fread_alf(FILE *fp) {
    fprintf(stderr, "ALF read support is unavailable\n");
    return NULLRead;
}

int fwrite_alf(char *fn, Read *read) {
    fprintf(stderr, "ALF write support is unavailable\n");
    return -1;
}

#endif

/* --------------------------------- PLN --------------------------------- */
#ifdef NO_PLN
Read *read_pln(char *fn) {
    fprintf(stderr, "PLN read support is unavailable\n");
    return NULLRead;
}

int write_pln(char *fn, Read *read) {
    fprintf(stderr, "PLN write support is unavailable\n");
    return -1;
}

Read *fread_pln(FILE *fp) {
    fprintf(stderr, "PLN read support is unavailable\n");
    return NULLRead;
}

int fwrite_pln(FILE *fp, Read *read) {
    fprintf(stderr, "PLN write support is unavailable\n");
    return -1;
}

#endif


/* --------------------------------- CTF --------------------------------- */
#ifdef NO_CTF
Read *read_ctf(char *fn) {
    fprintf(stderr, "CTF read support is unavailable\n");
    return NULLRead;
}

int write_ctf(char *fn, Read *read) {
    fprintf(stderr, "CTF write support is unavailable\n");
    return -1;
}

Read *fread_ctf(FILE *fp) {
    fprintf(stderr, "CTF read support is unavailable\n");
    return NULLRead;
}

int fwrite_ctf(FILE *fp, Read *read) {
    fprintf(stderr, "CTF write support is unavailable\n");
    return -1;
}

#endif



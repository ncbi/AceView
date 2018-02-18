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
 * Keep our compilers happy for our unavailable library formats
 */

/* --------------------------------- ABI --------------------------------- */

Read *read_abi(char *fn);
Read *fread_abi(FILE *fp);
int write_abi(char *fn, Read *read);
int fwrite_abi(FILE *fp, Read *read);

/* --------------------------------- ALF --------------------------------- */
int write_alf(char *fn, Read *read);
int fwrite_alf(FILE *fp, Read *read);
Read *read_alf(char *fn);
Read *fread_alf(FILE *fp);

/* --------------------------------- PLN --------------------------------- */
int write_pln(char *fn, Read *read);
int fwrite_pln(FILE *fp, Read *read);
Read *read_pln(char *fn);
Read *fread_pln(FILE *fp);

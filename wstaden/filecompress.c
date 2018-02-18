/*  Last edited: Jun 18 10:29 1999 (fw) */
/*
 * Handles compression and decompression.
 * Two functions are available. One compresses files, and the other opens
 * (read only) a compressed file and returns a FILE pointer.
 * Neither of these two are likely to work under Windows or MacOS.
 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include "filecompress.h"
#define BS 8192

/*
 * This contains the last used compression method.
 */
static int compression_used = 0;

typedef struct {
    char magic[2];
    long int offset;
    char *compress;
    char *uncompress;
    char *extension;
} Magics;

/*
 * The list of magic numbers. The attempted order for compression is the
 * order of entries in this file.
 *
 * NB: bzip gives very good (better than gzip) results, is sometimes faster for
 * compression, but unfortunately much slower (4x?) for decompression. Most
 * people won't have it anyway.
 */
static Magics magics[] = {
    {{'B',   'Z'},	0,	"bzip",		"bzip -d <",	".bz"},
    {{'\037','\213'},	0,	"gzip",		"gzip -d <",	".gz"},
    {{'\037','\235'},	0,	"compress",	"uncompress < ",".Z"},
    {{'\037','\036'},	0,	"pack",		"pcat",		".z"},
};

void set_compression_method(int method) {
    compression_used = method;
}

int get_compression_method(void) {
    return compression_used;
}

/*
 * Converts compress mode strings (eg "gzip") to numbers.
 */
int compress_str2int(char *mode) {
    if (strcmp(mode, "bzip") == 0)
	return 1;
    else if (strcmp(mode, "gzip") == 0)
	return 2;
    else if (strcmp(mode, "compress") == 0)
	return 3;
    else if (strcmp(mode, "pack") == 0)
	return 4;
    else return 0;
}

/*
 * Compress a file using the method set in the compression_used value
 * (set by set_compression_method and fopen_compressed).
 *
 * If compression succeeds, we rename the file back its original name.
 *
 * When compression_used is 0 no compression is done.
 */
int compress_file(char *file) {
    int ret;
    char buf[2048];
    struct stat statbuf;

    /* Do nothing unless requested */
    if (compression_used == 0)
	return 0;

    /* Execute the compression program */
    sprintf(buf, "%s %s 1>/dev/null 2>/dev/null",
	    magics[compression_used-1].compress, file);
    if ((ret = system(buf)) != 0) {
	if (ret == -1)
	    perror(buf);
	else
	    fprintf(stderr, "%s: compression failed\n", file);
	return 1;
    }

    /* Rename the file back */
    if (-1 == stat(file, &statbuf) && errno == ENOENT) {
	sprintf(buf, "%s%s", file, magics[compression_used-1].extension);
	rename(buf, file);
    }

    return 0;
}

/*
 * Compress a file using the method set in the compression_used value
 * (set by set_compression_method and fopen_compressed).
 *
 * If compression succeeds, we rename the file back its original name.
 *
 * When compression_used is 0 no compression is done.
 */
int fcompress_file(FILE *fp) {
    char buf[BS];
    char *fname;
    FILE *newfp;
    int len, fd;

    /* Do nothing unless requested */
    if (compression_used == 0)
	return 0;

    /* It's also impossible if it's stdout as we can't rewind */
    if (fp == stdout)
	return 0;

    fname=strdup("/tmp/acedb.XXXXXX");
    if (! fname)
	return 0;

    fd = mkstemp(fname);
    if (fd < 0)
	return 0;

    newfp = fdopen(fd,"wb+");
    if (!newfp)
	return 0;

    fflush(fp);
    rewind(fp);
    do {
	len = fread(buf, 1, BS, fp);
	if (len > 0)
	    fwrite(buf, 1, len, newfp);
    } while (!feof(fp));
    fflush(newfp);
    fclose(newfp);

    /* Compress it */
    compress_file(fname);

    /* Copy it back */
    if (NULL == (newfp = fopen(fname, "rb")))
	return 0;

    rewind(fp);
    do {
	int len = fread(buf, 1, BS, newfp);
	if (len > 0)
	    fwrite(buf, 1, len, fp);
    } while (!feof(newfp));
    ftruncate(fileno(fp), ftell(fp));
    remove(fname);
    
    return 0;
}


/*
 * Returns a file pointer of an uncompressed copy of 'file'.
 * 'file' need not exist if 'file'.ext (eg file.gz)
 * exists and can be uncompressed.
 *
 * If ofp is non NULL then the original file pointer will also be returned
 * (opened for update) to allow writing back to the original file. In cases
 * of uncompressed data this is the same as the returned file pointer.
 */
FILE *fopen_compressed(char *file, FILE **ofp) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    int i, ret, fd, try, do_del = 1;
    char buf[2048], fext[1024], mg[3];
    char *fname, *fptr=0;
    FILE *fp;

    if (NULL == (fname = tmpnam(NULL)))
	return NULL;

    /*
     * Try opening the file and reading the magic number.
     * If this doesn't work, then don't worry - the filename may be
     * the original name which has been renamed due to compression.
     * (eg file.gz).
     */
    try = 1;
    fd = open(file, O_RDONLY);
    if (fd != -1) {
	if (2 != read(fd, mg, 2)) {
	    close(fd);
	} else {
	    try = 0;
	    fptr = file;
	}
    }
    for (i = 0; i < num_magics; i++) {
	/* If necessary, try opening the file as 'file'.extension */
	if (try) {
	    sprintf(fext, "%s%s", file, magics[i].extension);
	    fptr = fext;
	    if (-1 == (fd = open(fext, O_RDONLY)))
		continue;
	    
	    if (2 != read(fd, mg, 2)) {
		close(fd);
		continue;
	    }
	}
	
	/* Check the magic numbers */
	if (mg[0] == magics[i].magic[0] && mg[1] == magics[i].magic[1]) {
	    sprintf(buf, "%s %s 1>%s 2>/dev/null",
		    magics[i].uncompress, fptr, fname);
	    if ((ret = system(buf)) == 0) {
		compression_used = i+1;
		break;
	    }
	}
	
	if (try && fd != -1)
	    close(fd);
    }
    
    if (fd != -1) close(fd);
    
    if (i == num_magics) {
	/*
	 * It's only an error if the file couldn't be found. If it
	 * exists then we'll assume that it doesn't need uncompressing.
	 */
	if (try) {
	    return NULL;
	} else {
	    do_del = 0;
	    compression_used = 0;
	    fname = file;
	}
    }

    /*
     * We've now got the temporary file. Open it and unlink it.
     * We can also keep the original fp open for those who need to
     * write back to it.
     */
    if (NULL == (fp = fopen(fname, "r+b")))
	if (NULL == (fp = fopen(fname, "rb")))
	    return NULL;
    if (ofp) {
	if (file != fname)
	    *ofp = fopen(try ? fext : file, "r+b");
	else
	    *ofp = fp;
    }

    if (do_del) remove(fname);

    return fp;
}

/*
 * Returns a file pointer of an uncompressed copy of 'fp'.
 *
 * If ofp is non NULL then the original file pointer will also be returned
 * (opened for update) to allow writing back to the original file. In cases
 * of uncompressed data this is the same as the returned file pointer.
 */
FILE *freopen_compressed(FILE *fp, FILE **ofp) {
    int num_magics = sizeof(magics) / sizeof(*magics);
    char fname[L_tmpnam];
    char buf[BS], mg[3];
    FILE *newfp=0;
    int i;

    /* Test that it's compressed with 1 byte of magic (can't ungetc more) */
    mg[0] = fgetc(fp);
    ungetc(mg[0], fp);
    for (i = 0; i < num_magics; i++) {
	if (mg[0] == magics[i].magic[0]) {
	    break;
	}
    }
    if (i == num_magics) {
	return fp;
    }

    tmpnam(fname);

    /* Copy the file to newfp */
    newfp = fopen(fname, "wb+");
    if (NULL == newfp)
	return fp;

    do {
	int len = fread(buf, 1, BS, fp);
	if (len > 0)
	    fwrite(buf, 1, len, newfp);
    } while (!feof(fp));
    fflush(newfp);
    rewind(newfp);
    rewind(fp);
    
    /* Test that it's compressed with full magic number */
    fread(mg, 1, 2, newfp);
    rewind(newfp);
    for (i = 0; i < num_magics; i++) {
	if (mg[0] == magics[i].magic[0] && mg[1] == magics[i].magic[1]) {
	    break;
	}
    }
    if (i == num_magics) {
	/* Uncompressed, just return fp */
	remove(fname);
	return newfp;
    }

    /* It's compressed, so we use fopen_compressed on fname */
    fp = fopen_compressed(fname, ofp);
    remove(fname);
    if (newfp != NULL) {
	return fp;
    } else {
	return newfp;
    }
}

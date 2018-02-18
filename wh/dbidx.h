/*  Last edited: Jun  4 22:54 1996 (rd) */
/* 
   DBIDX.H

   Header file for EMBL sequence database indexing system

   Erik Sonnhammer, 930209

*/

/* $Id: dbidx.h,v 1.5 2015/09/11 19:32:29 mieg Exp $ */

/* LMB-MRC:
#define DEFAULT_SWDIR     "/nfs/al/pubseq/pubseq/seqlibs/swiss/"
#define DEFAULT_PIRDIR    "/nfs/al/pubseq/pubseq/seqlibs/pir/"
#define DEFAULT_WORMDIR   "/nfs/al/pubseq/pubseq/seqlibs/wormpep/"
#define DEFAULT_EMBLDIR   "/nfs/al/pubseq/pubseq/seqlibs/embl/"
#define DEFAULT_GBDIR     "/nfs/al/pubseq/pubseq/seqlibs/genbank/"
*/

/* SANGER: */
#define DEFAULT_SWDIR      "/nfs/disk3/pubseq/swiss/"
#define DEFAULT_PIRDIR     "/nfs/disk3/pubseq/pir/"
#define DEFAULT_WORMDIR    "/nfs/disk3/pubseq/wormpep/"
#define DEFAULT_EMBLDIR    "/nfs/disk3/pubseq/embl/"
#define DEFAULT_GBDIR      "/nfs/disk3/pubseq/genbank/"
#define DEFAULT_PRODOMDIR  "/nfs/disk3/pubseq/prodom/"
#define DEFAULT_PROSITEDIR "/nfs/disk3/pubseq/prosite/"


#define DEFAULT_IDXFILE   "entrynam.idx"
#define DEFAULT_DIVFILE   "division.lkp"  
#define DEFAULT_DBFILE    "seq.dat"
#define DEFAULT_ACTRGFILE "acnum.trg"
#define DEFAULT_ACHITFILE "acnum.hit"

#define max(a,b) (a > b ? a : b)

#include <stdio.h>
#include <ctype.h> 

#ifdef ACEDB
#include <mystdlib.h>
#if (!defined(WIN32)&&!defined(HP)&&!defined(CYGWIN)
  extern char  toupper(char),
               tolower(char); /* Otherwise I have to include regular.h */
#endif
#else
#include <stdlib.h>
#ifdef SUN
  extern void  fclose(FILE*),
              *memcpy(),
              *memset(),
               rewind();

  extern int   strlen(),
               strcmp(),
               strncmp(),
               getopt(),
               fprintf(),
               fwrite(),
               printf(),
               sscanf(),
               fseek(),
               fread(),
               fputc(),
               fputchar();
               
  extern char *strcpy(char*, char*),
              *strchr(),
               toupper(char),
               tolower(char);

#endif
#endif

#if !defined(max) /* should already be in (my)stdlib.h */
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

/* DATE as described in 3.4.2 in the EMBL CD-ROM Indices documentation 1.4
*/
typedef struct _Date 
{
  char noll;
  char year;
  char month;
  char day;

} Date;


/* HEADER contains the fields of the index file header
   as described in 3.4.2 in the EMBL CD-ROM Indices documentation 1.4
*/
typedef struct _Header {

  unsigned long  file_size;
  unsigned long  nrecords;
  unsigned short record_size;
  unsigned short entry_name_size;
           char  db_name[21];
           char  db_relnum[11];
           Date  db_reldate[4];

} Header;


/* ENTRYNAM contains a record of the index file ENTRYNAM.IDX
*/
typedef struct _Entrynam {

  char entry_name[32];  /* The actual size is calculated from the header of 
			   the index file:  size = record_size - 10 */
  unsigned long  ann_offset;
  unsigned long  seq_offset;
  unsigned short div_code;

} Entrynam;


/* ACNUM contains a record of the index file ACNUM.TRG
*/
typedef struct _Acnum {

  unsigned long  nhits;
  unsigned long  hit_rec_num;
  char acnum[32];  /* The actual size is calculated from the header of 
			   the index file:  size = record_size - 10 */
} Acnum;


/* ACHIT contains a record of the index file ACNUM.HIT
*/
typedef struct _Achit {

  unsigned long  entry_rec_num;

} Achit;



/* DIVISION contains a record of the division lookup table division.lkp
*/
typedef struct _Division {

  unsigned short div_code;
  char file_name[128]; /* This may vary. EMBL uses 12 */

} Division;


extern int 
    Reversebytes, 
    search, 
    display_head, 
    fasta,
    seqonly,
    raw_seq, 
    PIR_seq,
    ACsearch,
    AC2entryname,
    Exact_match,
    Mixdb,
    Mini;

Header *readheader(FILE *idx);
Entrynam *readEntrynam(FILE *idx, Header *head, int pos);
Entrynam *readnextEntrynam(FILE *idx, Header *head);
Acnum *readActrg(FILE *idx, Header *head, int pos);
void printEntrynam(Entrynam *rec, Header *head);
char *getname(FILE *idx, Header *head, int pos);
char *getdbfile(int div_code, char *divfile);
int bin_search(char *query, FILE *idx, Header *head, char *(*getname)(FILE *idx, Header *head, int pos) );
char *getAcname(FILE *idx, Header *head, int pos);

void printActrg(),
     printPIRseq(),
     printGBseq(),
     printEMBLseq(),
     printAnn(),
     printheader();

char *printFastaseq();

int  getAcentry_rec_num();


/********************** end of file *********************/
 

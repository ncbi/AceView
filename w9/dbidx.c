/*  Last edited: Oct 16 15:52 1998 (fw) */
/* 
   DBIDX.C

   Routines for EMBL CD-ROM Database Indices

   Erik Sonnhammer

Version Date   Modification
------- ------ ----------------------
0       930112 Created
1.0     930826 Added Wormpep support
1.1     940824 Changed bin_search to return record if unique substr, list of names
               if ambigous, unless -x is used.
*/

/* $Id: dbidx.c,v 1.5 2016/05/06 21:29:05 mieg Exp $ */ 

#include "dbidx.h"
#include <string.h> 
/*#include <ctype.h>*/
#include <stdarg.h>


#define  MAXLINE 1024

/* Global variables */
char dbsource[32]="";

int Reversebytes = 1;
int search       = 0; 
int display_head = 0; 
int fasta        = 0;
int seqonly      = 0;
int raw_seq      = 0;
int PIR_seq      = 0;
int ACsearch     = 0,
    AC2entryname = 0,
    Exact_match  = 1,
    Mixdb        = 0,
    Mini         = 0;

void fatal(char *format, ...)
{
    va_list  ap;

    fprintf(stderr, "\nFATAL ERROR: ");

    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);

    fprintf(stderr, "\n"); 
    exit(-1);
}


/* strcpy_noCR  Gets rid of these damned CR's in the EMBL db's
*/
void strcpy_noCR(char *dest, char *src)
{
    int i,j;

    for ( i=0, j=0; i<= strlen(src); i++ ) 
      if ((int)src[i] != 13) dest[j++] = src[i];
}


/* REVERSEBYTES changes order of n bytes at location ptr
   max capacity 256 bytes
*/
void reversebytes(char *ptr, int n)
{ 
  static char copy[256];
  int  i;

  memcpy(copy, ptr, n);  /* Note: strcpy doesn't work - stops at \0 */
  
  for(i=0; i<n; i++) 
    {
      ptr[i] = copy[n-i-1];
    }
}
  

int isseq(int c)
{
    return (isalpha(c) || (c == '*'));
}

Header *readheader(FILE *idx)
{
  Header *head;

  if (!(head = (Header *)malloc(sizeof(Header)))) {
      fprintf(stderr, "\nOut of Memory error\n");
      exit(1);
  }
  
  rewind(idx);
  if (fread(&head->file_size, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&head->file_size, 4);

  if (fread(&head->nrecords, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&head->nrecords, 4);

  if (fread(&head->record_size, 2, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&head->record_size, 2);

  head->entry_name_size = head->record_size - 10;
  /* Each record has 10 bytes in use (ver. 1.4)*/

  if (fread(head->db_name, 20, 1, idx) != 1) return NULL;
  head->db_name[20] = '\0';

  if (fread(head->db_relnum, 10, 1, idx) != 1) return NULL;
  head->db_relnum[10] = '\0';

  if (fread(head->db_reldate, 4, 1, idx) != 1) return NULL;

  /* set pointer to first record */
  fseek(idx, 300, 0); 

  return head;
}

void printheader(Header *head)
{
  printf("Size of index file    = %lu\n", head->file_size);
  printf("Number of records     = %lu\n", head->nrecords);
  printf("Size of index records = %u\n", head->record_size);
/*  printf("Size of entry names   = %u\n", head->entry_name_size); */
  printf("Database name         = %s\n", head->db_name);
  printf("Database Release nr   = %s\n", head->db_relnum);
  printf("Database Release date = %u-%u-%u\n\n", head->db_reldate->year,
  head->db_reldate->month, head->db_reldate->day);
}

void printActrg(Acnum *rec)
{
  printf("Number of hits             = %lu\n", rec->nhits);
  printf("Record number of first hit = %lu\n", rec->hit_rec_num);
  printf("Accession number           = %s\n\n", rec->acnum);
}


Acnum *readActrg(FILE *idx, Header *head, int pos)
{
  static Acnum rec;
  int name_size = head->record_size - 10;

  fseek(idx, 300 + head->record_size*pos, 0);

  if (fread(&rec.nhits, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.nhits, 4);

  if (fread(&rec.hit_rec_num, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.hit_rec_num, 4);

  if (fread(rec.acnum, name_size, 1, idx) != 1) return NULL;
  rec.acnum[name_size] = '\0';

  return &rec;
}


Entrynam *readEntrynam(FILE *idx, Header *head, int pos)
{
  static Entrynam rec;
  int name_size = head->record_size - 10;
  char *space;

  fseek(idx, 300 + head->record_size*pos, 0);

  if (fread(rec.entry_name, name_size, 1, idx) != 1) return NULL;
  rec.entry_name[name_size] = '\0';
  if(( space = (char *)strchr(rec.entry_name, ' ') )) *space = '\0';


  if (fread(&rec.ann_offset, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.ann_offset, 4);

  if (fread(&rec.seq_offset, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.seq_offset, 4);

  if (fread(&rec.div_code, 2, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.div_code, 2);

  return &rec;
}


Entrynam *readnextEntrynam(FILE *idx, Header *head)
{
  static Entrynam rec;
  int name_size = head->record_size - 10;

  if (fread(rec.entry_name, name_size, 1, idx) != 1) return NULL;
  rec.entry_name[name_size] = '\0';

  if (fread(&rec.ann_offset, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.ann_offset, 4);

  if (fread(&rec.seq_offset, 4, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.seq_offset, 4);

  if (fread(&rec.div_code, 2, 1, idx) != 1) return NULL;
  if (Reversebytes) reversebytes((char *)&rec.div_code, 2);

  return &rec;
}


void printEntrynam(Entrynam *rec, Header *head)
{
  printf("Entry name        = %s\n", rec->entry_name);
  printf("Annotation offset = %lu\n", rec->ann_offset);
  printf("Sequence offset   = %lu\n", rec->seq_offset);
  printf("Division code     = %u\n\n", rec->div_code);
}


/* PRINTSEQ prints the acutal sequence of entries from all databases
   (SwissProt, EMBL, Genbank and Wormpep)
*/
void printSeq(FILE *db, Entrynam *rec, int Start, int End)
{
  static char buff[MAXLINE+1];
  int i, Seqlen=0, Seqpos=0, Line_has_seq;

  /* Find out the sequence length by counting letters */
  fseek(db, rec->seq_offset, 0);
  fgets(buff, MAXLINE, db);

  /* Rodger's ridiculous Sequence Offset for PIR ... */
  if (!strncmp(buff, "SEQUENCE", 8)) { fgets(buff, MAXLINE, db); fgets(buff, MAXLINE, db); }

  while( !feof(db) && strncmp(buff, "//", 2) && buff[0] != '>')
  {
      for(i=0; i<strlen(buff); i++) 
	if (raw_seq) Seqlen++; 
	else if ( isseq(buff[i]) ) Seqlen++;
      fgets(buff, MAXLINE, db);
  }

  if (Start == -1) Start = 1;
  if (End   == -1) End = Seqlen;

  if (Start <= 0 || Start > Seqlen || End <= 0 || End > Seqlen)
      fatal("Coordinates out of range(1-%d): %d-%d\n", Seqlen, Start, End);

  fseek(db, rec->seq_offset, 0);
  fgets(buff, MAXLINE, db);

  /* Rodger's ridiculous Sequence Offset for PIR ... */
  if (!strncmp(buff, "SEQUENCE", 8)) { fgets(buff, MAXLINE, db); fgets(buff, MAXLINE, db); }

  while( !feof(db) && strncmp(buff, "//", 2) && buff[0] != '>' )
  {
      Line_has_seq = 0;

      for(i=0; i<strlen(buff); i++)
	if (raw_seq) putchar(buff[i]);
	else if ( isseq(buff[i]) )
	{
	  Seqpos++;
	  if (Seqpos >= Start && Seqpos <= End ) 
	  {
	      putchar(buff[i]);
	      Line_has_seq = 1;
	  }
	}
      if (Line_has_seq && !seqonly) putchar('\n');
      *buff = 0;
      fgets(buff, MAXLINE, db);
  }
  if (seqonly) putchar('\n');
}


/* PRINTPIRSEQ prints the sequence of PIR entries in FASTA format
*/
void printPIRseq(FILE *db, Entrynam *rec, int Start, int End)
{
  static char buff[MAXLINE+1];
  int i, TITLE=0;
  char junk[32], ACnr[32], TITLEline[MAXLINE]; 
  char *slash;

  strcpy(TITLEline, " ");

  /* Get the ACCESSION NUMBER */
  fseek(db, rec->ann_offset, 0);
  fgets(buff, MAXLINE, db);
  while( !feof(db) && strncmp(buff, "//", 2) )
  {
      if (!strncmp(buff, "ACCESSION", 9) )
      {
	  sscanf(buff, "%s%s", junk, ACnr);
          break;
      }
      fgets(buff, MAXLINE, db);
  }
  /* Check if ACCESSION number was found - pir3 doesn't have it,
     so add name as accession number if there is no  
  */
  if ( !strncmp(buff, "//", 2) ) strcpy( ACnr, rec->entry_name );

  if(( slash = (char *)strchr(ACnr, '\\') )) *slash = '\0';
  if (!seqonly) printf("%s ", ACnr);

  if (Mini) printf("\n");

  /* Get the TITLE line */
  fseek(db, rec->ann_offset, 0);
  fgets(buff, MAXLINE, db);
  TITLE = 0;
  while( !feof(db) && strncmp(buff, "//", 2))
  {
      if( TITLE && buff[0] != ' ' ) break;

      if( !strncmp(buff, "TITLE", 5) ) TITLE = 1;

      if( TITLE )
      {
	  if (Mini) printf(" %s", buff+16);
	  strcpy( &TITLEline[strlen(TITLEline) - 1], &buff[16] );
      }
      fgets(buff, MAXLINE, db);
  }
  for ( i=0; i< strlen(TITLEline); i++ ) 
      TITLEline[i] = toupper(TITLEline[i]);
  if (!Mini && !seqonly) printf("%s", TITLEline);

  TITLE = 0;
  if (Mini) while( !feof(db) && strncmp(buff, "//", 2))
  {
      if( TITLE && buff[0] != ' ' ) break;
      if( !strncmp(buff, "KEYWORDS", 8) ) TITLE = 1;
      if (TITLE) printf(" %s", buff+16);

      fgets(buff, MAXLINE, db);
  }
  
  /* Print the SEQUENCE */
  printSeq(db, rec, Start, End);
}


/* PRINTGBSEQ prints the sequence in FASTA format
   between Start and End, specific for GenBank
   */
void printGBseq(FILE *db, Entrynam *rec, int Start, int End)
{
    static char buff[MAXLINE+1];
    int i, TITLE=0;
    char junk[32], ACnr[32], TITLEline[MAXLINE]; 
    char *slash;
    
    strcpy(TITLEline, " ");

    
/*** Get the ACCESSION NUMBER */

    fseek(db, rec->ann_offset, 0);
    fgets(buff, MAXLINE, db);
    while( !feof(db) && strncmp(buff, "//", 2) )
	{
	    if (!strncmp(buff, "ACCESSION", 9) )
		{
		    sscanf(buff, "%s%s", junk, ACnr);
		    break;
		}
	    fgets(buff, MAXLINE, db);
	}
    
    if(( slash = (char *)strchr(ACnr, '\\') )) *slash = '\0';
    if (!seqonly) printf("%s ", ACnr);
    
    
/*** Get the TITLE line */

    fseek(db, rec->ann_offset, 0);
    fgets(buff, MAXLINE, db);
    TITLE = 0;
    while( !feof(db) && strncmp(buff, "//", 2) )
	{
	    if( !strncmp(buff, "DEFINITION", 10) ) TITLE = 1;
	    
	    if( TITLE )
	    {
		if (!strncmp(buff, "ACCESSION", 9)) break;
		strcpy_noCR( TITLEline + strlen(TITLEline) - 1, &buff[12] );
	    }
	    fgets(buff, MAXLINE, db);
	}
    for (i=0; i< strlen(TITLEline); i++) TITLEline[i] = toupper(TITLEline[i]);
    if (!seqonly) printf("%s", TITLEline);
    

/*** Print the SEQUENCE */

    printSeq(db, rec, Start, End);
}


/* PRINTEMBLSEQ prints out the sequence in FASTA format
   between Start and End, specific for EMBL and SWISSPROT
*/
void printEMBLseq(FILE *db, Entrynam *rec, int Start, int End)
{
  static char buff[MAXLINE+1];
  int  TITLE=0;
  char TITLEline[MAXLINE];
  char *slash;

  strcpy(TITLEline, " ");

  /* Do the ACCESSION NUMBER and DE field */
  fseek(db, rec->ann_offset, 0);
  fgets(buff, MAXLINE, db);
  while( !feof(db) && strncmp(buff, "//", 2) )
  {
      if (!strncmp(buff, "AC   ", 5) )
      {
	  /* Accession number starts on pos 5 and ends with a ; */
	if(( slash = (char *)strchr(buff, ';') )) *slash = '\0';
	  if (!seqonly) printf("%s", &buff[5]);
	  break;
      }
      fgets(buff, MAXLINE, db);
  }
  if( feof(db) || !strncmp(buff, "//", 2) )
  {
      fprintf(stderr, "Cannot find accession number\n"); exit(1);
  }

  if (Mini) printf("\n");

  /* Get the De lines */
  fgets(buff, MAXLINE, db);
  TITLE = 0;
  while( !feof(db) && strncmp(buff, "//", 2))
  {
      if (TITLE && strncmp(buff, "DE   ", 5) ) break;

      if (!strncmp(buff, "DE   ", 5) ) TITLE = 1;

      if (TITLE)
      {
	  if (Mini) printf(" %s", buff+5);
	  strcpy_noCR( &TITLEline[strlen(TITLEline) - 1], &buff[4] );
      }
      fgets(buff, MAXLINE, db);
  }

  /* Get rid of these damned CR's * /
  strcpy(TITLE2, TITLEline);
  for ( i=0, j=0; i<= strlen(TITLE2); i++ ) 
      if ((int)TITLE2[i] != 13) TITLEline[j++] = TITLE2[i];
  for ( i=0; i< strlen(TITLEline); i++ ) TITLEline[i] = toupper(TITLEline[i]);
  */

  if (!Mini && !seqonly) printf("%s", TITLEline);

  if (Mini) 
      while( !feof(db) && strncmp(buff, "//", 2))
  {
      if (/*!strncmp(buff, "CC   ", 5)  ||*/ !strncmp(buff, "KW   ", 5) ) 
	  printf("%s", buff+4);
      fgets(buff, MAXLINE, db);
  }
  
  /* Print the SEQUENCE */
  printSeq(db, rec, Start, End);
}


/* PRINTFASTASEQ prints out the sequence in FASTA format (native format)
   between Start and End, specific for Fasta databases
*/
char *printFastaseq(FILE *db, Entrynam *rec, int Start, int End)
{
  static char buff[MAXLINE+1];

  /* Do the DE field */
  fseek(db, rec->ann_offset, 0);
  fgets(buff, MAXLINE, db);
  if (!seqonly) printf("%s", (char *)strchr(buff, ' ')+1);

  /* Print the SEQUENCE */
  printSeq(db, rec, Start, End);

  return buff;
}


/* PRINTANN prints out the annotation of an entry in the dbfile db
*/
void printAnn(FILE *db, Entrynam *rec)
{
  static char buff[MAXLINE+1], clean[MAXLINE], pads[12], ichar[12];
  int i, mixedSeqstarted, fastaStarted, PDstarted;

  fseek(db, rec->ann_offset, 0);

  fgets(buff, MAXLINE, db);
  mixedSeqstarted = fastaStarted = PDstarted = 0;
  while( !feof(db) && strncmp(buff, "//", 2))
  {
      if (buff[0] == '>' && fastaStarted) break;  /* End of Fasta type db */
      fastaStarted = 1;

      strcpy_noCR(clean, buff);
      printf("%s", clean);
      fgets(buff, MAXLINE, db);

      /* Scale for ProDom families */
      if (*dbsource == 'D' && !strchr(rec->entry_name, '_')) {
	  if (!PDstarted) {
	      printf("                     ");
	      memset(pads, '-', 10);
	      for (i=10 ; i< strlen(buff) - 22; i += 10)
	      {
		  sprintf(ichar, "%d", i);
		  pads[10 - strlen(ichar)] = '\0';
		  printf("%s%d", pads, i);
	      }
	      for (i=0; i < (strlen(buff) - 22) % 10; i++) printf("-");
	      printf("\n");
	      PDstarted = 1;
	  }
      }

      if (Mixdb && buff[0] != ' ' && !mixedSeqstarted)  {
	  printf("\n");  /* printf ("\nSequence:\n"); */
	  mixedSeqstarted = 1;
      }
  }

  /* Print record separator */
  strcpy_noCR(clean, buff);
  if (!feof(db) && !Mixdb && buff[0] != '>') printf("%s", clean);
}
  

/* GETNAME gets the name of the n'th sequence in ENTRYNAM.IDX
*/
char *getname(FILE *idx, Header *head, int pos)
{
  static char name[16];
  int name_size;
  char *space;

  name_size = head->record_size - 10;

  fseek(idx, 300 + head->record_size*pos, 0);

  if (fread(name, name_size, 1, idx) != 1) {
    fprintf(stderr, "\nError reading index file at record %d\n", pos);
    exit(1);
  }

  name[name_size] = '\0';
  if ((space = (char *)strchr(name, ' '))) *space = 0;
  
  return name;
}


/* GETDBFILE returns the division database filename 
   given a division code and lookup table file
   If file can't be opened or is corrupted, NULL is returned
*/
char *getdbfile(int div_code, char *divfile)
{
  FILE *div;
  Header *divhead;
  int name_size, i;
  unsigned short code;
  static char name[128];

  if ((div = fopen(divfile, "r")) == NULL) return NULL;

  divhead = readheader(div);
  if (display_head) printheader(divhead);
  
  name_size = divhead->record_size - 2;

  fseek(div, 300, 0);

  for (i=0; i < divhead->nrecords; i++)
  {
    if (fread(&code, 2, 1, div) != 1) return NULL;
    if (Reversebytes) reversebytes((char *)&code, 2);

    if (fread(name, name_size, 1, div) != 1) return NULL; 

    if (code == div_code) break;
    }

  fclose(div);
  free(divhead);

  /* Terminate string after last character */
  for (i=0; i<strlen(name); i++) if (isspace(name[i])) name[i] = '\0';
  name[name_size] = '\0';

  /* Convert to lower case since the EMBL CD-ROM division.lkp contains uppercase chars */
  for (i=0; i<strlen(name); i++) name[i] = tolower(name[i]);

  return name;
}


/* GETacNAME gets the name of the n'th sequence in ACNUM.TRG
*/
char *getAcname(FILE *idx, Header *head, int pos)
{
    Acnum *a;

    a = readActrg(idx, head, pos);
    return a->acnum;
}

/* BIN_SEARCH performs a binary search in a sorted table of records
   and return the position of the record, starting with pos 0.
   If not found, it returns -1.
   If ambiguous, lists names that are possible.
*/

int bin_search(char *query, FILE *idx, Header *head, char *(*getname)(FILE *idx, Header *head, int pos) )
{
     int recless, recmore, pos, pos1, pos2, balance, newpos, qlen, i;
     char *posname, *space;
 
     recless  = 0;
     recmore  = head->nrecords;
     pos      = (recless + recmore)/2;
     posname  = getname(idx, head, pos);
     if ((space = (char *)strchr(posname, ' '))) *space = 0;
     qlen     = strlen(query);
     /* int matchlen = ( Exact_match ? max(strlen(posname), qlen) : qlen ); */

     /* Do binary search */
     while ((balance = strncmp(posname, query, qlen)))
     {
	 if (search) printf("%d %s %lu %d\n", pos, posname, strlen(posname), qlen);

	 if (balance > 0)  /* query is left of pos */
	 {
	     recmore = pos;
	     newpos = (recless + recmore)/2;
	 }
	 else  /* query is right of pos */
	 {
	     recless = pos;
	     newpos = (recless + recmore)/2;
	 }
	 if (newpos == pos) return -1;

	 pos = newpos;
	 posname = getname(idx, head, pos);
     }
     if (search) printf("%d %s %ld %d\n", pos, posname, strlen(posname), qlen);
     if (!Exact_match) return pos;


     /***** Check if hit is unique *****/

     /* If query matches entire dbname it must be unique */
     if (!strcmp(posname, query)) return pos;

     posname = getname(idx, head, pos-1);
     while (!strncmp(posname, query, qlen) && pos > 0) {
	 pos--;
	 posname = getname(idx, head, pos-1);
     }
     pos1 = pos;

     posname = getname(idx, head, pos);
     while (!strncmp(posname, query, qlen)) {

	 /* If query matches entire dbname it must be unique */
	 if (!strcmp(posname, query)) return pos;

	 pos++;
	 if (pos == head->nrecords) break;
	 posname = getname(idx, head, pos);
     }
     pos--;
     pos2 = pos;

     if (pos1 == pos2) return pos;
     else 
     {
	 fprintf(stderr, "\n\"%s\" is ambiguous, choose between:\n\n", query);
	 for (pos = pos1, i=1; pos <= pos2; pos++, i++) {
	     fprintf(stderr, "%-10s ", getname(idx, head, pos));
	     if (!(i % 6)) fprintf(stderr, "\n");
	 }
	 fprintf(stderr, "\n");
	 exit(0);
     }
}


int getAcentry_rec_num(FILE *achit, Header *achithead, int offset)
{
	unsigned long entry_rec_num;
	int seekpos;

	seekpos = 300 + offset*achithead->record_size;

        if (fseek(achit, seekpos, 0))
	{
	    fprintf(stderr, "\nError: Can't go to record %d in ACNUM.HIT \n",
		    offset);
	    exit(1);
	}
        if (fread(&entry_rec_num, 4, 1, achit) != 1) {
	  fprintf(stderr, "\nError reading ACNUM.HIT at offset %d\n", offset);
	  exit(1);
        }
        if (Reversebytes) reversebytes((char *)&entry_rec_num, 4);
        if (display_head) printf("Entry record number: %lu\n", entry_rec_num);
	return (int)(entry_rec_num - 1);
}


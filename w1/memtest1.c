/*
* using memtest:
*
* In memsubs.c, enable the malloc_log version of myMalloc() then compile
* and run your program.  It records all malloc() and free() operations
* that pass through myMalloc(), halloc(), etc.
*
* When you run your program, you get a file "malloc_log".  Run this
* program (memtest1) to translate it to "malloc_log.cooked", then 
* run memtest2 to reproduce the same sequence of malloc/free operations.
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

int cmpfn( const void *a , const void *b )
{
	int *aa, *ab;
	aa = (int *)a;
	ab = (int *)b;
	return *aa - *ab;
}

main(int argc, char **argv)
{
	FILE *f, *fout;
	int *memarray, *p;
	int malloc_count, free_count, mc, addr, size, i, n, realloc_count, addr2;
	char b[100];
	char *fname;
	fname = "malloc_log";
	if (argv[1])
		fname = argv[1];
	f = fopen(fname,"r");
	if (!f)
		{ perror(fname); exit(1); }

	fout = fopen("malloc_log.cooked","w");
	if (!fout)
		{ perror("malloc_log.cooked"); exit(1); }

	malloc_count = 0;
	free_count = 0;
	realloc_count = 0;

	for (;;)
		{
		switch (getc(f))
			{
		case 'r':
			realloc_count++;
			getc(f); getc(f); getc(f); getc(f);
			getc(f); getc(f); getc(f); getc(f);
			getc(f); getc(f); getc(f); getc(f);
			break;
		case 'm':
			malloc_count++;
			getc(f); getc(f); getc(f); getc(f);
			getc(f); getc(f); getc(f); getc(f);
			break;
		case 'f':
			free_count++;
			getc(f); getc(f); getc(f); getc(f);
			break;
		case -1:
			goto counted;
		default:
			printf("bad record type at %d\n",ftell(f));
			exit(1);
			}
		}
counted:
	printf("mallocs %d  reallocs %d  frees %d\n",malloc_count, realloc_count, free_count);

	memarray = malloc((malloc_count + realloc_count ) * sizeof(int));

	fseek(f, 0L, 0);

	mc = 0;

	for (;;)
		{

		switch (getc(f))
			{
		case 'r':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			addr2=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			size=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			break;
		case 'm':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			size=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			memarray[mc] = addr;
			mc++;
			break;
		case 'f':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			break;
		case -1:
			goto arrayed;
		default:
			printf("bad record type at %d\n",ftell(f));
			exit(1);
			}
		}
arrayed:
	printf("arrayed\n");

	qsort( memarray, mc, sizeof(int *), cmpfn);

	printf("sorted\n");

	fseek(f, 0L, 0);

	for (;;)
		{
		switch (getc(f))
			{
		case 'r':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			addr2=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			size=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			break;
		case 'm':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			size=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			p = bsearch( &addr, memarray, mc, sizeof(int*), cmpfn);
			if (!p)
				{ printf("malloc not found in array! %x\n", addr); break; }
			putc( 'm', fout );
			i = p - memarray;
			putc( i, fout);
			putc( i >> 8, fout);
			putc( i >> 16, fout);
			putc( i >> 24, fout);

			putc( size, fout);
			putc( size >> 8, fout);
			putc( size >> 16, fout);
			putc( size >> 24, fout);

			break;
		case 'f':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			p = bsearch( &addr, memarray, mc, sizeof(int*), cmpfn);
			if (!p)
				{ printf("free not found in array! %x\n", addr); break; }
			putc( 'f', fout );
			i = p - memarray;
			putc( i, fout);
			putc( i >> 8, fout);
			putc( i >> 16, fout);
			putc( i >> 24, fout);
			break;
		case -1:
			goto cooked;
		default:
			printf("bad record type at %d\n",ftell(f));
			exit(1);
			}
		}
cooked:
	printf("cooked\n");

}

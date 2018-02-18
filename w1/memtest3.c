/*
* produce output for graphing malloc/free activity
*/

#include <stdio.h>
#include <malloc.h>

main()
{
	FILE *f;
	int *memarray;
	int malloc_count, free_count, mc, addr, size, i, n;
	int total;
	char b[100];
	int bucket_number, bucket_mallocs, mallocs_per_bucket;
 	int bucket_events, events_per_bucket, currently_allocated, 
		bucket_allocated, bucket_freed;

	f = fopen("malloc_log.cooked","r");
	if (!f)
		{ perror("malloc_log.cooked"); exit(1); }

	mallocs_per_bucket = 5000;

	bucket_number = 0;	/* which bucket are we in */
	bucket_mallocs = 0;	/* how many mallocs in this bucket */

	for (;;)
		{
		switch (getc(f))
			{
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
			printf("bad record type at 0x%x %d 0%o \n",ftell(f),
				ftell(f), ftell(f));
			exit(1);
			}
		}
counted:

	memarray = calloc(malloc_count * sizeof(int), 1);

	fseek(f, 0L, 0);

	mc = 0;

	total = 0;

	bucket_events = 0;
	malloc_count  = 0;
	free_count  = 0;
	currently_allocated  = 0;
	bucket_allocated  = 0;
	bucket_freed  = 0;
	events_per_bucket = 10000;
	
	for (;;)
		{
		bucket_events++;
		if (bucket_events >= events_per_bucket)
			{
			printf("%d ", bucket_number);
			printf("%d ", malloc_count);
			printf("%d ", free_count);
			printf("%d ", currently_allocated);
			printf("%d ", bucket_allocated);
			printf("%d ", bucket_freed);
			printf("%d ", bucket_allocated - bucket_freed);
			printf("\n");
			bucket_number++;
			bucket_allocated = 0;
			bucket_freed = 0;
			bucket_events = 0;
			}
		switch (getc(f))
			{
		case 'm':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			size=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			memarray[addr] = size;
			currently_allocated = currently_allocated + size;
			bucket_allocated = bucket_allocated + size;
			malloc_count++;
			break;
		case 'f':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			currently_allocated = currently_allocated - memarray[addr];
			bucket_freed = bucket_freed + memarray[addr];
			memarray[addr] = 0;
			free_count++;
			break;
		case -1:
			goto done;
		default:
			printf("bad record type at 0x%x %d 0%o \n",ftell(f),
				ftell(f), ftell(f));
			exit(1);
			}
		}
done: 	;
}

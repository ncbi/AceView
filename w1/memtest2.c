#include <stdio.h>
#include <malloc.h>

main()
{
	FILE *f;
	int **memarray;
	int malloc_count, free_count, mc, addr, size, i, n;
	int total, event_count;
	char b[100];
	f = fopen("malloc_log.cooked","r");
	if (!f)
		{ perror("malloc_log.cooked"); exit(1); }
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
	printf("mallocs %d  frees %d\n",malloc_count, free_count);

	memarray = calloc(malloc_count * sizeof(void *), 1);

	fseek(f, 0L, 0);

	mc = 0;

	printf("overhead done - press return\n");
	gets(b);

	total = 0;
	event_count = 0;

	for (;;)
		{
		if (event_count++ > 1000*5000)
			{
			printf("pause\n");
			gets(b);
			event_count = 0;
			}
		switch (getc(f))
			{
		case 'm':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			size=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			if (memarray[addr])
				printf("block %d allocated when already allocated\n",addr);
			memarray[addr] = malloc(size);
			*memarray[addr] = size;
			total = total + size;
			break;
		case 'f':
			addr=getc(f)|(getc(f)<<8)|(getc(f)<<16)|(getc(f)<<24);
			if (memarray[addr] == 0)
				printf("block %d freed when already free\n", addr);
			else
				{
				free(memarray[addr]);
				total = total - *memarray[addr];
				}
			memarray[addr] = 0;
			break;
		case -1:
			goto done;
		default:
			printf("bad record type at 0x%x %d 0%o \n",ftell(f),
				ftell(f), ftell(f));
			exit(1);
			}
		}
done:
		printf("done loop - memory still allocated %d\n", total);
		gets(b);
}

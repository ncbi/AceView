#include <stdio.h>
#include <string.h>
#include <ctype.h>


#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

/*
* On Linux, many of the rusage fields are always 0, including 
* whichever one you are most interested in.
*/
void printrusage( char * arg )
{
static struct rusage old_r = { 0 };
static struct rusage t_old_r;
struct rusage r;
if (getrusage(RUSAGE_SELF, &r) < 0)
	{ perror("getrusage"); return; }

t_old_r = r;

r.ru_utime.tv_sec -= old_r.ru_utime.tv_sec;
r.ru_utime.tv_usec -= old_r.ru_utime.tv_usec;
if (r.ru_utime.tv_usec < 0)
	{ r.ru_utime.tv_sec++; r.ru_utime.tv_usec += 1000000; }

r.ru_stime.tv_sec -= old_r.ru_stime.tv_sec;
r.ru_stime.tv_usec -= old_r.ru_stime.tv_usec;
if (r.ru_stime.tv_usec < 0)
	{ r.ru_stime.tv_sec++; r.ru_stime.tv_usec += 1000000; }

r.ru_maxrss -= old_r.ru_maxrss;
r.ru_minflt -= old_r.ru_minflt;
r.ru_majflt -= old_r.ru_majflt;

fprintf(stderr,
"%s\nu %5d.%06d s %5d.%06d maxrss %-5d minflt %-5d majflt %-5d\n",
arg,
r.ru_utime.tv_sec, r.ru_utime.tv_usec,
r.ru_stime.tv_sec, r.ru_stime.tv_usec,
r.ru_maxrss, r.ru_minflt, r.ru_majflt);

old_r = t_old_r;

}


int main(int argc, char **argv)
{
	Associator *a;

	a = assCreate();

	printf("foo\n");

}


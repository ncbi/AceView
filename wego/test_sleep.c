#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

#include "wego.h"

/*
* Each thing that we might do in parallel is implemented in a function.
* That function has a struct full of input parameters and a struct full
* of result values.
*/

/***********************************************************************
* Here is the first function, named "think"
*/

/*
* here are the parameters; you can put whatever you want here.
*/
struct think_parameters 
	{
	int num;
	int delay;
	int queue_time;
	int verbose;
	};

/*
* here is the result; you can put whatever you want here.  It
* can be useful to return some information that tells you what
* the results are about.  In this example, the "num" field is
* just copied from the parameters.
*/
struct think_result 
	{
	int num;
	int queue_time;
	int start;
	int end;
	};

/*
* here is the function to run.  It doesn't do anything useful
* except serve as an example.
*/
void think( struct wego_task *t, void *parameters, void *result )
{
	struct think_parameters *p = (struct think_parameters *)parameters;
	struct think_result *r = (struct think_result *)result;

	if (p->verbose) 
		printf("think %d start\n",p->num);

	/*
	* When the main program picks up the result, it does not know
	* what the parameters were.  Copy some identifying information
	* from input to output.
	*/
	r->num = p->num;
	r->queue_time = p->queue_time;

	/*
	* do something (nothing much interesting in this example)
	*/
	r->start = time(0);
	sleep(p->delay);
	r->end = time(0);

	if (p->verbose) 
		printf("think %d end\n",p->num);
}

/*
* end of the function "think"
***********************************************************************/


int main(int argc, char **argv, char **envp)
{
	int x;
	struct think_result *r;

	struct wego_group *g;
	struct wego_task *t;

	#define N 10
	struct wego_task *a[N];

	/*
	* Set the maximum number of threads that can run concurrently.
	*
	* There is no speed benefit to making this number larger than
	* the number of idle CPUs in your system.
	*
	* For example, on an idle dual-processor system, you might
	* as well use 2.  On an 8 processor system with 1 processor
	* already completely busy, you would use 7.
	*
	* Your process WILL have more threads than you specify here,
	* but not all will be running.
	*/
	if (argv[1])
		wego_max_threads (atoi(argv[1]));
	else
		wego_max_threads (5);

	/*
	* create a group of tasks
	*/
	g = wego_new_group();

	for (x=0; x<N; x++)
		{
		/*
		* fill in a parameter structure. the parameters are
		* copied; you do not need to save the struct after
		* creating the task.
		*/
		struct think_parameters s;
		s.num = x;
		s.delay= 1 ;
		s.queue_time = time(0);
		s.verbose = 0;
		/*
		* run the task as part of a group.  The NULL for the return value
		* means that the library allocates memory for it; you can also
		* give a pointer to your own memory, but do not try to read it
		* before the task is finished.
		*/
		printf("CREATE %d\n",x);
		a[x] = wego_run(g,	/* group to join (NULL for no group) 	*/
			think, 		/* function to call 			*/
			&s, sizeof(s), 	/* struct of function parameters 	*/
			NULL, sizeof(*r)/* struct of return values 		*/
			);

		}

	/*
	* You can ask for the result of any particular wego_task at any time.
	* It blocks if the result is not ready.
	*/
	printf("RESULT 9\n");
	r = wego_result(a[9]);
	printf("	%d: %d %d %d  %d %d\n",
		r->num,r->queue_time, r->start,r->end,
		r->end-r->start, r->end - r->queue_time);

	/*
	* Or you can ask for the result of any task that is finished.
	*/
	for (;;)
		{
		/*
		* wego_fin finds a wego_task that is:
		*	in the group
		*	finished
		* 	has not had results collected by a call to wego_result
		* It returns NULL if it is impossible to meet those conditions.
		* (That is, if every task in the group is finished and has had
		*   it's results examined, wego_fin() returns NULL.)
		*
		* If no task is finished yet, it will wait for a task to finish.
		*/
		printf("FIN\n");
		t = wego_fin(g);
		if (!t)
			break;

		/*
		* wego_fin returns a task; we want the result struct.  This
		* looks just like what we did above.
		*/
		printf("RESULT\n");
		r = wego_result(t);
		printf("	%d: %d %d %d  %d %d\n",
			r->num,r->queue_time, r->start,r->end,
			r->end-r->start, r->end - r->queue_time);
		}
	return 0;
}

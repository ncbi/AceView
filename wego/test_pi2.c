/*
* PI example, using tasks within tasks
*
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "wego.h"

double splits = 10;

struct parameters 
	{
	/*
	* area to cover
	*/
	double	min_x, max_x, min_y, max_y;
	/*
	* step increments
	*/
	double  step_x, step_y;

	int	count, count2;
	int	verbose;
	};

/***********************************************************************
* second level task - computes small squares
*/

struct l2_result 
	{
	/*
	* number of points sampled
	*/
	int num_points;
	/*
	* number of points inside circle
	*/
	int circle_area;
	};

void l2_func( struct wego_task *t, void *parameters, void *result )
{
	struct parameters *p = (struct parameters *)parameters;
	struct l2_result *r = (struct l2_result *)result;
	double x,y;
	int circle_area = 0, num_points = 0;

	if (p->verbose) 
		printf("l2 %d.%d\t %f %f - %f %f start\n",p->count,p->count2,p->min_x,p->max_x,
			p->min_y, p->max_y);

	for (x=p->min_x; x<p->max_x; x=x+p->step_x)
		for (y=p->min_y; y<p->max_y; y=y+p->step_y)
			{
			double d = sqrt(x*x+y*y);
			if (d <= 1.0)
				circle_area++;
			num_points++;
			}

	r->num_points = num_points;
	r->circle_area = circle_area;

	if (p->verbose) 
		printf("l2 %d.%d\t %f %f end %d %d\n",p->count,p->count2,p->min_x,p->min_y,r->num_points,r->circle_area);
}

/***********************************************************************
* first level task - subdivide across y
*/

struct l1_result 
	{
	/*
	* number of points sampled
	*/
	int num_points;
	/*
	* number of points inside circle
	*/
	int circle_area;
	};

/*
* here is the function to run.  It doesn't do anything useful
* except serve as an example.
*/
void l1_func( struct wego_task *t, void *parameters, void *result )
{
	struct parameters *p = (struct parameters *)parameters;
	struct l1_result *r = (struct l1_result *)result;
	struct parameters l2p;
	struct l2_result *l2r;
	double x;
	int circle_area = 0, num_points = 0;
	struct wego_group *g;
	struct wego_task *tt;
	int count2 = 0;

	if (p->verbose) 
		printf("l1 %d.\t %f %f - %f %f start\n",p->count,p->min_x,p->max_x,
			p->min_y, p->max_y);

	g = wego_new_group();

	for (x=0; x<splits; x++)
		{
		/*
		* fill in a parameter structure. the parameters are
		* copied; you do not need to save the struct after
		* creating the task.
		*/
		l2p = *p;
		l2p.min_x = x / splits;
		l2p.max_x = (x+1)/ splits;
		l2p.count2 = count2++;
		if (p->verbose)
			printf("L1 %d.\t run\n",p->count);
		wego_run(g,	           /* group to join (NULL for no group) */
			l2_func, 		  /* function to call 			*/
			&l2p, sizeof(l2p),/* struct of function parameters 	*/
			NULL, sizeof(*l2r)/* struct of return values 		*/
			);

		}

	num_points = 0;
	circle_area = 0;
	for (;;)
		{
		if (p->verbose)
			printf("L1 %d.\t fin\n",p->count);
		tt = wego_fin(g);
		if (!tt)
			break;

		if (p->verbose)
			printf("L1 %d.\t result\n",p->count);
		l2r = wego_result(tt);
		if (p->verbose)
			printf("L1 %d.\t result = %d %d\n",p->count,l2r->num_points,l2r->circle_area);
		num_points += l2r->num_points;
		circle_area += l2r->circle_area;
		}


	r->num_points = num_points;
	r->circle_area = circle_area;

	if (p->verbose) 
		printf("l1 %d.\t %f %f end\n",p->count,p->min_x,p->min_y);
}


/***********************************************************************
* main program
*/

int main(int argc, char **argv, char **envp)
{
	int y, count;
	struct l1_result *r;
	struct wego_group *g;
	struct wego_task *t;
	int num_points, circle_area;
	struct parameters p;
	int verbose =1;

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

	p.min_x = 0.0;
	p.max_x = 1.0;
	p.step_x = .00005;
	p.step_x = .0001;
	p.step_y = p.step_x;
	p.verbose=verbose;
	count =0;

	for (y=0; y<splits; y++)
		{
		/*
		* fill in a parameter structure. the parameters are
		* copied; you do not need to save the struct after
		* creating the task.
		*/
		p.min_y = y / splits;
		p.max_y = (y+1)/ splits;
		p.count = count++;
		if (verbose)
			printf("MAIN run\n");
		wego_run(g,	/* group to join (NULL for no group) 	*/
			l1_func, 		/* function to call 			*/
			&p, sizeof(p), 	/* struct of function parameters 	*/
			NULL, sizeof(*r)/* struct of return values 		*/
			);

		}

	num_points = 0;
	circle_area = 0;
	for (;;)
		{
		if (verbose)
			printf("MAIN fin\n");
		t = wego_fin(g);
		if (!t)
			break;
		if (verbose)	
			printf("MAIN result\n");
		r = wego_result(t);
		if (verbose)
			printf("MAIN\t result = %d %d\n",r->num_points,r->circle_area);
		num_points += r->num_points;
		circle_area += r->circle_area;
		}

	printf("num_points	%d\n",num_points);
	printf("circle_area	%d\n",circle_area);
	printf("PI		%f\n",4.0*(double)circle_area/(double)num_points);
	return 0;
}

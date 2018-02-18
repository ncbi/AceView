/*
* Example program using multiple threads to calculate the 
* value of PI.
*
* Consider this monte carlo method:
*
* The area of a circle is PI * R^2.  
* Assume R = 1.  A quarter circle is inscribed in a
* square that covers from 0,0 to 1,1.
*
* Select a random point in the square.
* Compute the distance from 0,0 to the point.
* If the distance is <= 1, the point is in the circle;
* Add 1 to circle_area.
*
* After numerous random points are selected, 
*	circle_area / number_of_points == PI/4
*
*
* Now say we want something deterministic:  Instead of
* randomly selecting points, systematically choose them
* from a regular grid.  This ensures a uniform distribution over
* the range from 0,0 to 1,1.
*
* You can divide the points up into multiple tasks,
* each covering a certain slice of the square.  Each task
* returns the total circle area in that slice of the square.
* Add the results from each task, and 
* PI = 4 * circle_area / num_points
*
* This is only intended to demonstrate parallel processing with
* the wego library.  I know this is a crappy way to compute PI.
*
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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
	/*
	* area to cover
	*/
	double	min_x, max_x, min_y, max_y;
	/*
	* step increments
	*/
	double  step_x, step_y;

	int	verbose;
	};

/*
* here is the result; you can put whatever you want here.  It
* can be useful to return some information that tells you what
* the results are about.  In this example, the "num" field is
* just copied from the parameters.
*/
struct think_result 
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
void think( struct wego_task *t, void *parameters, void *result )
{
	struct think_parameters *p = (struct think_parameters *)parameters;
	struct think_result *r = (struct think_result *)result;
	double x,y;
	int circle_area = 0, num_points = 0;

	if (p->verbose) 
		printf("think %f %f - %f %f start\n",p->min_x,p->max_x,
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
		printf("think %f %f end\n",p->min_x,p->min_y);
}

/*
* end of the function "think"
***********************************************************************/


int main(int argc, char **argv, char **envp)
{
	int y;
	struct think_result *r;
	struct wego_group *g;
	struct wego_task *t;
	int num_points, circle_area;
	struct think_parameters p;
	double splits;

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
	p.step_y = p.step_x;
	p.verbose=1;

	splits = 10;
	for (y=0; y<splits; y++)
		{
		/*
		* fill in a parameter structure. the parameters are
		* copied; you do not need to save the struct after
		* creating the task.
		*/
		p.min_y = y / splits;
		p.max_y = (y+1)/ splits;
		wego_run(g,	/* group to join (NULL for no group) 	*/
			think, 		/* function to call 			*/
			&p, sizeof(p), 	/* struct of function parameters 	*/
			NULL, sizeof(*r)/* struct of return values 		*/
			);

		}

	num_points = 0;
	circle_area = 0;
	for (;;)
		{
		t = wego_fin(g);
		if (!t)
			break;

		r = wego_result(t);
		num_points += r->num_points;
		circle_area += r->circle_area;
		}

	printf("num_points	%d\n",num_points);
	printf("circle_area	%d\n",circle_area);
	printf("PI		%f\n",4.0*(double)circle_area/(double)num_points);
	return 0;
}

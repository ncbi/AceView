
From: "Thierry-Mieg, Jean  (NIH/NLM/NCBI)" <mieg@ncbi.nlm.nih.gov>
To: mieg@ncbi.nlm.nih.gov
Subject: colors X11 bug
Date: Sun, 26 Sep 2004 23:45:43 -0400
MIME-Version: 1.0

File libx.c for developing an X11 graphics driver 
Previous <http://www.cs.colorado.edu/~mcbryan/5229.03/mail/35.htm>
Next <http://www.cs.colorado.edu/~mcbryan/5229.03/mail/37.htm>  


/***********************************************************************
******
** file: libx.c
** This version works also with TrueColor displays!
**
**
**  Compile this with the test drawing program as follows:
**  gcc test_x11.c libx.c -lm -lX11 -L/usr/X11R6/lib -o test_x11
************************************************************************
******/
/*
*				NOTES ON X11 LIBRARY
*
*
* The following notes need to be supplemented by reading some X source
such as:
*
* Oliver Jones: Introduction to the X Window System, Prentice-Hall
* and (for reference)
* O'Reilly & Associates; Xlib Reference Manual Volume 2, O'Reilly &
Associates
*
* I grappled with trying to write notes for all this - and realized that
I would
* be rewriting the above books.  It is very hard to do it in a few
pages.
*
* So instead I wrote a bunch of working routines that do what you need.
* You will need to modify them to fit your situation.   I have tested
them using
* them in a program that has the form:
*	prog.c lib2d.c libx.c
*
* What I have done is to document various items so that you can get a
feel for
* what needs to be done.
*
* WARNING: the drawing routines supplied here accept INTEGER
coordinates.
*
*
*
* The DISPLAY variable:
*
* The X library is set to work in a default mode where you run the
program
* on the console of the machine you are using.  This is a natural mode
of
* use.   However, often there is a reason to run the program on a remote
* machine, but with output on the screen of your local system.  For
example your
* local system might be a slow SUN computer, and you want to run the
program
* on a fast remote SUN computer for speed, while sitting at your slow
SUN
screen.
*
* In such cases you need to set the DISPLAY environment variable ON THE
REMOTE
* MACHINE (SUNFAST).  For example:
*
* on SUNFAST:machine:0.0"
*
* where machine is the name of the SUNSLOW.
*
* Also the LOCAL MACHINE (SUNSLOW) needs to authorize access to its
screen by
* remote machines, the easiest way being to type:
*
* on SUNSLOW:
* xhost +
*
* which gives access to ALL remote systems (see man xhost).
*
* You are then ready to run your program on the SUNFAST, and it will pop
up a
* window on the SUNSLOW
*
* This code may need some small modifications for specific compilers and
* operating systems.
*
*/



/*
*				libx.c:
*
*	Device level 2D graphics driver for X11.
*	The "screen" is an X11 window with real pixels
*	X11 measures coordinates in the y direction from the top right
corner
*	of a window, which is point 0,0 of the window.  This needs to be
*	remembered in using this package.  The routine _enquire_device
reports
*	this fact by returning ydir = -1.
*
*
*				Visible Routines:
*
*
*		Initialization Routines:
*
*	_init_graphics(),
*		Must be called first.
*		Creates an X Window of size X11_COLS by X11_ROWS with a
*		default set of X11_COLORS colors. The colors consist of
a
*		default table of 8 basic colors (colors 0-7):
*			BLACK,RED,GREEN,BLUE,CYAN,MAGENTA,YELLOW,WHITE
*		and a default palette of X11_COLORS-8 smoothly varying
shades.
*		Here X11_COLS, X11_ROWS, X11_COLORS are environment
variables
*		and default to 512, 512, 128 respectively if not set.
*
*	_enquire_device(int *smin,int *smax,int *tmin,int *tmax,int
*ydir),
*		Sets smin,tmin and smax,tmax to the screen coordinates
of the
*		window boundaries: (0,X11_ROWS-1) and (X11_COLS-1,0).
*		Sets ydir to -1, signifying that y increases downwards
on screen.
*		That is how y is defined in X11 unfortunately.
*
*
*
*		Routines for Drawing:
*
*	_move(int x,int y),
*		Moves current point to window location x,y
*
*	_lineto(int x,int y),
*		Draws a line from current point to window location x,y
*		moves current point to window location x,y
*
*	_polygon(int col,int numv,int *vx,int *vy),
*		Draws a polgon, filled with color col from the color
table.
*		the polygon has numv vertices, and vx,vy are arrays
giving
*		the vertices in clockwise order.  The last vertex is
*		understood to connect over to the first one.  The
polygon
*		is one-sided - only the side from which the vertices
appear
*		clockwise is drawn.  The col here is a color entry in
the
*		Global color table - see below.
*
*	_text(char *string),
*		Draws a string in a standard font/size at the current
point.
*		Leaves the current point unchanged.
*
*
*
*		Routines for Handling Color - Overview:
*
*	The workstation is presumed to have a large color table with
some
*	maximum number of colors. We will call this the Global color
table.
*	In X this size is dynamic because other windows may already have
grabbed
*	a bunch of allowed colors.  We provide a set of routines below
that
*	allow the Global color table to be managed as a set of color
tables or
*	color palettes.
*
*	We allow the creation of up to MAX_COLOR_TABLES color tables
numbered
*	from 0.  Each color table has a number num_colors of colors
indexed from
*	0 to num_colors-1. The routine create_color_table() is used to
create a
*	new table and specify its table number and its number of colors.
*		colors = create_color_table(table,num_colors);
*	A table can have at most PIXELS colors in it.
*	A new table has no colors set.  A table may be freed later using
*	free_color_table().
*
*	The set of colors contained in all created color tables, taking
the
*	tables in increasing order, and the colors in order in each
table,
*	defines the Global color table.  It is introduced only for
backward
*	compatibility to previous versions of this library.  Only one
*	routine - _polygon(col,..) uses a col from this Global table.
*
*	Colors in each table can be set by the user using either of the
*	routines setup_table_color() or setup_table_colors().
*
*	Each table can be used either as a table of num_colors discrete
colors
*	or as a palette of continuously varying colors indexed from 0.0
to 1.0.
*	In the latter case, the colors are indexed to the closet value
in the
*	corresponding table.  For example if a table has 5 colors, a
palette
*	value of .95 would use the 5th color.  Use the routines:
*		set_color_from_Table()
*		set_color_from_Palette()
*	to set the current drawing color.
*
*	On startup, _init_graphics() creates two default color tables,
numbered
*	0 and 1, and defines these to be the default color table and
color
*	palette.  If you want to create your own tables or palettes, it
is
*	advisable to at least free table 1, since it consumes so many
colors.
*	Simplified versions of some  routines work on the default
tables:
*		set_color_from_table(), set_color_from_palette().
*	You can redefine the default color table and palette using:
*		set_default_color_table(), set_default_color_palette()
*
*	Drawing (lines, text) is always done in the current color.  The
current
*	color is set to be color 1 from the current default table by
*	_init_graphics(), _end_frame(), and _erase().
*	Polygons however are filled using a color specified as an
argument -
*	generally needed because 3D polygons are typically not drawn at
the
*	time they are created, but much later after hidden surface
steps.
*	as an argument.
*
*
*		Routines for Creating, Destroying Color Tables:
*
*	int available_colors():
*		Returns the number of unused Global color table colors.
*		Useful for deciding how many colors you can play with
below.
*
*	status = create_color_table(int table,int num_colors):
*		Creates a new table, with num_colors colors
*		No colors are actually defined in the table.
*		If num_colors colrs are not available, it does the best
possible.
*		Use setup_table_color(s) to store colors into the table.
*		status >0 indicates the number of colors allocated
*		status  0 indicates there are no available colors
*		status -1 indicates there are already MAX_TABLES tables.
*		status -2 indicates a request for more than PIXELS
colors
*		status -3 indicates the requested table already exists.
*		status -4 indicates failure to allocate space for the
colors.
*
*	status = free_color_table(table):
*		Frees a previously created table.
*		status = 0 signals an unknown table
*
*	set_default_color_table(table),
*		defines the default color table to be table
*		_init_graphics initializes this to table 0
*
*	set_default_color_palette(table),
*		defines the default color palette to be taken from table
table
*		_init_graphics initializes this to table 1
*
*
*		Routines for Setting up Colors in a Table:
*
*	status = setup_table_color(int table, int col,
*				unsigned red,unsigned green,unsigned
blue),
*
*		Defines color color in table table to be the specified
*		RGB value.  Each of R, G, B should be in range [0,65535]
*		status = 1 on success, 0 on failure.
*
*
*	status = setup_table_colors(int table,int num_colors,int
start_col,
*				unsigned *red,unsigned *green, unsigned
*blue),
*
*		Defines a set of num_colors colors in the  table table.
The
*		colors start at color start_col and run to color
*		start_col+num_colors-1.  The three arrays specify the
RGB values
*		and each element should be in the range [0,65535].
*		status = 1 on success, 0 on failure.
*
*
*		Routines for Setting Current Color:
*
*	status = set_color_from_Table(int table, int col),
*		Sets the current drawing color to color col from the
color
*		table table.
*
*	status = set_color_from_table(int col),
*		Sets the current drawing color to color col from the
default
*		color table.
*
*	status = set_color_from_Palette(table,double pal),
*		sets the current drawing color to color pal from the
default
*		color palette.
*
*	status = set_color_from_palette(double pal),
*		Sets the current drawing color to color pal from the
default
*		color palette.
*
*
*
*		Routines for Conversion between Color Representations:
*
*	status = color_from_Palette(int table,double pal,int *col),
*		Sets *col to the color corresponding to color 0 <= pal
<= 1 in
*		table table.
*
*	status = Table_to_Global(int table,int col,int *gcol);
*		Given a table and color col, sets *gcol to the
corresponding
*		Global color.
*
*	status = Global_to_Table(int gcol,int *table,int *col);
*		Given a Global table color gcol, sets *table to the
table
*		containing gcol, and *col to the color in that table.
*
*	status = Palette_to_Global(int table,double pal,int *gcol);
*		Given a table and color col, sets *gcol to the
corresponding
*		Global color.
*
*	_enquire_colors(int* num_table,int *num_palette);
*		Returns the number of colors in the default color table
and
*		palette.
*
*
*
*		Routines for Intensity Control:
*
*	set_linear_intensity(Numi,Mini,Maxi)
*		Initializes a range of Numi intensities in [Mini,Maxi]
*		The intensities are linearly seperated.
*
*	int linear_intensity(i)
*		Returns the ith linear intensity, or -1 on error.
*
*	set_eye_intensity(Numi,Mini,Maxi)
*		Initializes a range of Numi intensities in [Mini,Maxi]
*		The intensities are multiplicatively seperated.
*
*	int eye_intensity(i)
*		Returns the ith multiplicative intensity, or -1 on
error.
*
*
*
*		Routines for End of Frame, Erase, Cleanup:
*
*	_end_frame(),
*		display the current view
*
*	_erase(),
*		erase the screen - not effective until end_frame called!
*
*	_close_graphics().
*/





/*
*				Reference Sources:
*
* The following notes need to be supplemented by reading some X source
such as:
*
* Oliver Jones: Introduction to the X Window System, Prentice-Hall
* and (for reference)
* O'Reilly & Associates; Xlib Reference Manual Volume 2, O'Reilly &
Associates
*
*	Oliver McBryan
*/

/*
*				X DISPLAY:
*
*	The name of the display opened is obtained from the environment
*	variable DISPLAY; if DISPLAY is not set then the string "" is
used.
*	DISPLAY may be used to run the program on one work-station and
have
*	the image appear elsewhere.
*/

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

/******** MM: all prototypes to avoid warnings.. ******/
void _init_graphics();
static void x_init_graphics(int nx,int ny,int num_colors);
void _enquire_device(int *smin,int *smax,int *tmin,int *tmax,int *ydir);
static void init_window(int Nx,int Ny, int Numcolors,int *depth);
static void init_pixmap(int nx,int ny,int depth);
void _move(int x, int y);
void _cont(int x,int y);
void _line(int x1,int y1,int x2,int y2);
void _polygon(int gcol,int nvert,int *vx,int *vy);
void _Polygon(int nvert,int *vx, int *vy);
void _text(char *str);
static void xlabel(char *str); /* MM: should this be '_xlabel' */
static void init_colors(Display *cur_display,int screen,int depth, ulong
myforeground, ulong mybackground);
int available_colors();
int create_color_table(int table,int new_colors);
int free_color_table(int table);
static int good_table(int table);
int setup_table_color(int table,int color, unsigned red,unsigned
green,unsigned
blue);
int setup_table_colors(int table,int ncolor, int start, unsigned
*red,unsigned
*green,unsigned *blue);
void set_default_color_table(int table);
void set_default_color_palette(int table);
void _enquire_colors(int *numt,int *nump);
int set_color_from_Table(int table,int col);
int set_color_from_table(int col);
int set_color_from_Palette(int table,double pal);
int set_color_from_palette(double pal);
int color_from_Palette(int table,double pal,int *col);
int Table_to_Global(int table,int color,int *gcolor);
int Global_to_Table(int gcolor,int *table,int *color);
int Palette_to_Global(int table,double pal,int *gcol);
void set_linear_intensity(int Numi,int Mini,int Maxi);
int linear_intensity(int i);
void set_eye_intensity(int Numi,int Mini,int Maxi);
int eye_intensity(int i);
void _end_frame();
static void do_plot();
static void check_pixmap(Display *display,Pixmap pixmap); 
void _erase();
void _close_graphics();
static void x_close_graphics();
void flush_graphics();

/*****************************************************/

#define MAX_COLOR_TABLES 10

static char hello[] = {"Welcome to X Visualization"};


	/* Basic X data-structures shared by several routines below: */
static Display *mydisplay;
static int myscreen;
static int DDepth;	/* The number of bits of color */
static Window mywindow;
static GC mygc,pixgc;


static int NX, NY;      /* Window Size */

static int debug = 0;



		/* I: INITIALIZATION ROUTINES: */

/*
*				_init_graphics():
*
*	Routine called by lib2d.c to get the window created.
*	To get a bigger window, increase 512 below.
*	Default window is 512x512 with 128 colors.
*/

void _init_graphics()
{

 	int nx=512,ny=512,total_colors=128;
	char *s,*getenv();

	if ((s=getenv("X11_COLS"))) nx = atoi(s);
	if ((s=getenv("X11_ROWS"))) ny = atoi(s);
	if ((s=getenv("X11_COLORS"))) total_colors = atoi(s);
	x_init_graphics(nx,ny,total_colors);
}






/*
*				x_init_graphics():
*
*	This is an internal routine only.
*	Create a window as a specified size in pixels, and a specified
*	number of colors.  The first 8 colors (numbers 0-7) will be
black,
*	red, blue, green, cyan, magenta, yellow, white.
*	The colors from 8 to num_colors-1 will form a palette of
smoothly
*	varying colors.
*/

static void x_init_graphics(nx,ny,num_colors)
	int nx,ny,num_colors;
{

	NX = nx;
	NY = ny;
	init_window(nx,ny,num_colors,&DDepth);
	init_pixmap(nx,ny,DDepth);
}





/*
*				_enquire_device():
*
*	lib2d can call this routine to learn the size of the window
*	Not needed unless you want to use it.
*	Somehow lib2d must be told the X window size.
*/


void _enquire_device(smin,smax,tmin,tmax,ydir)
int *smin,*smax,*tmin,*tmax,*ydir;
{
	*smin = 0;
	*smax = NX-1;
	*tmin = NY-1;
	*tmax = 0;
	*ydir = -1;
}

/*
*				init_window():
*
*	init_window does all of the dirty work X requires to create a
window
*	on a display.  The name of the display opened is obtained from
the
*	environment variable DISPLAY; if DISPLAY is not set then the
string
""
*	is used.
*/

static int nx,ny;
static unsigned long myforeground, mybackground;

static void init_window(Nx,Ny,Numcolors,depth)
	int Nx,Ny,Numcolors,*depth;
{
	static int first = 1;
	int i,npalette,nalloc;
	char *getenv(), *display_name;
	static XSizeHints myhint;
	static Visual *visual;

	if (first) first = 0;
	else return;

	nx = Nx;
	ny = Ny;

		/* Open the Display and a Screen on it: */

	if ((display_name=getenv("DISPLAY")) == NULL) display_name =
"";
	if ((mydisplay = XOpenDisplay(display_name)) == NULL) {
		fprintf(stderr,"XOpenDisplay(display=%s)
failed\n",display_name);
		exit(1);
	}
	if (debug) printf("mydisplay = %d\n",(int)mydisplay);

	myscreen = DefaultScreen(mydisplay);
        /* MM: Changed the following line...  */
	if (debug) printf("myscreen = %d\n",myscreen);

	*depth = XDefaultDepth(mydisplay,myscreen);
	if (debug) printf("depth = %d\n",*depth);

	visual = XDefaultVisual(mydisplay,myscreen);
	if (debug) {
		if (*depth==1) {
			printf("One-plane Monochrome\n");
		}
		else {
			switch (visual->class) {
			case PseudoColor: printf("Multiplane
PseudoColor\n"); break;
			case GrayScale: printf("Multiplane
GrayScale\n"); break;
			case DirectColor: printf("Direct Color\n");
break;
			case TrueColor: printf("True Color\n"); break;
			case StaticColor: printf("Static Color\n");
break;
			case StaticGray: printf("Static Gray \n");
break;
			default: printf("Unknown visual class
%d\n",visual->class);
			}
		}
	}

		/* Create the colors needed: */

	myforeground = WhitePixel(mydisplay,myscreen);
	if (debug) printf("myforeground = %ld\n",myforeground);
	mybackground = BlackPixel(mydisplay,myscreen);
	if (debug) printf("mybackground = %ld\n",mybackground);
	
init_colors(mydisplay,myscreen,*depth,myforeground,mybackground);


		/* Set the 8 table colors: */
	if ((nalloc=create_color_table(0,8))!=8) {
		fprintf(stderr,
		"Able to allocate only %d colors of requested %d\n",
		nalloc,Numcolors);
		exit(1);
	}
	set_default_color_table(0);
	setup_table_color(0,0,0,0,0);
	setup_table_color(0,1,65535,0,0);
	setup_table_color(0,2,0,65535,0);
	setup_table_color(0,3,0,0,65535);
	setup_table_color(0,4,0,65535,65535);
	setup_table_color(0,5,65535,0,65535);
	setup_table_color(0,6,65535,65535,0);
	setup_table_color(0,7,65535,65535,65535);

		/* Set the palette colors: 8 .. Numcolors-1 */
	npalette = Numcolors - 8;
	if ((nalloc=create_color_table(1,npalette))!=npalette) {
		fprintf(stderr,
		"Able to allocate only %d colors of requested %d\n",
		nalloc+8,Numcolors);
		exit(1);
	}
	set_default_color_palette(1);
	/*
	set_eye_intensity(npalette,1024,65535);
        for (i=0; i<npalette; i++)
	
setup_table_color(1,i,eye_intensity(i),eye_intensity(i)/2,0);
	*/
	set_linear_intensity(npalette,1024,65535);
        for (i=0; i<npalette; i++)
	
setup_table_color(1,i,linear_intensity(i),linear_intensity(i),0);

		/* Hint to X on pixel location of window corner: */
	myhint.x = 32; myhint.y = 32;
		/* Hint to X on on pixel size of window : */
	myhint.width = nx; myhint.height = ny;
	myhint.flags = PPosition | PSize;

		/* Now create the window! */
	if ((mywindow = XCreateSimpleWindow(mydisplay,
		DefaultRootWindow(mydisplay), myhint.x,myhint.y,

		myhint.width,myhint.height,
5,myforeground,mybackground)) == (Window) NULL) {
		fprintf(stderr,"XCreateSimpleWindow failed\n");
		exit(1);
	}
	if (debug) printf("mywindow=%d\n", mywindow);

	XSetStandardProperties(mydisplay,mywindow,hello,hello,
				None, NULL,1,&myhint);
				/* Should be argc,argv instead of NULL,1
*/

		/* Create a Graphics Context (GC) for the window */

	if ((mygc = XCreateGC(mydisplay,mywindow,0,0)) == NULL) {
		fprintf(stderr,"Cannot created GC for window\n");
		exit(1);
	}
	if (debug) printf("mygc=%d\n",mygc);

	XSetBackground(mydisplay,mygc,mybackground);
	XSetForeground(mydisplay,mygc,myforeground);
	XClearWindow(mydisplay,mywindow);

		/* Specify the input events to care about: */
	XSelectInput(mydisplay,mywindow,ButtonPressMask | KeyPressMask |
		ExposureMask);

		/* Make the window visible: */
	XMapRaised(mydisplay,mywindow);
}







/*
*	Create a Pixmap - frame buffer equal in size to the window.
*	We will write into the pixmap.  Then display it on the screen.
*/

static Pixmap pixmap;
#define SURFACE pixmap

static void init_pixmap(nx,ny,depth)
	int nx,ny,depth;
{
	if ((pixmap = XCreatePixmap(mydisplay,mywindow,nx,ny,depth)) ==
(Pixmap) NULL)
{
		fprintf(stderr,"Failed to create pixmap\n");
		exit(1);
	}
	else check_pixmap(mydisplay,pixmap);

	if (debug) printf("pixmap = %d\n",pixmap);

	if ((pixgc = XCreateGC(mydisplay,pixmap,0,0)) == NULL) {
		fprintf(stderr,"Failed to create pixgc\n");
		exit(1);
	}
	if (debug) printf("pixgc = %d\n",pixgc);

	XSetBackground(mydisplay,pixgc,mybackground);
	set_color_from_table(0);
	XFillRectangle(mydisplay,pixmap,pixgc,0,0,nx,ny);
	set_color_from_table(1);
}











/*
*			II: DRAWING ROUTINES:
*/



/*
*				_move():
*
*	 Move the current point to x,y.
*/

static int xc, yc;	/* The current Point */

void _move(x,y)
        int x,y;
{
        xc = x;
        yc = y;
}


/*
*                               _cont():
*
*       Draw a line from current point to point x,y.
*       To implement _line() call _move then _cont.
*/

void _cont(x,y)
        int x,y;
{
        XDrawLine(mydisplay,SURFACE,pixgc,xc,yc,x,y);
        _move(x,y);
}



void _line(x1,y1,x2,y2)
        int x1,y1,x2,y2;
{
        XDrawLine(mydisplay,SURFACE,pixgc,x1,y1,x2,y2);
        _move(x2,y2);
}



/*
*                               _polygon():
*
*       Draw a 2D polygon with a color gcol from the Global color table.
*       The nvert vertices are in integer coords in the arrays vx, vy.
*/

                /* Largest number of vertices on a polygon: */
#define MAX_POINTS 100


void _polygon(gcol,nvert,vx,vy)
        int gcol;
        int nvert;
        int *vx, *vy;
{
        int table,col;

        if (!Global_to_Table(gcol,&table,&col)) {
                fprintf(stderr,"_polygon called with invalid color
%d\n",gcol);
                exit(1);
        };

        set_color_from_Table(table,col);
        _Polygon(nvert,vx,vy);
}



/*
*				_Polygon():
*
*	Draw a 2D polygon filled with the current color.
*	The nvert vertices are in integer coords in the arrays vx, vy.
*/

		/* Largest number of vertices on a polygon: */
void _Polygon(nvert,vx,vy)
	int nvert;
	int *vx, *vy;
{
	static XPoint p[MAX_POINTS];	/* X needs Xpoint data
structures */
	int i;

	if (nvert > MAX_POINTS) {
		fprintf(stderr,"Not more than %d points in polygons\n",
	
MAX_POINTS);
		exit(1); /* Or exit */
	}
	for (i=0; i<nvert; i++) {
		p[i].x = (short) vx[i];
		p[i].y = (short) vy[i];
	}
	XSetFillRule(mydisplay,pixgc,EvenOddRule);
	
XFillPolygon(mydisplay,SURFACE,pixgc,p,nvert,Complex,CoordModeOrigin);
}






/*
*				_text()
*
*	Draw a text string at the current point.
*	This will be called by text() in lib2d.
*	No control is provided on the font or size.
*/

void _text(str)
  char *str;
{
	XDrawImageString(mydisplay,SURFACE,pixgc,xc,yc,str,strlen(str));
}





/*
*				_xlabel():
*
*	Draw a string directly on the window at position 10,50:
*	This is purely intended for internal use.
*/

static void xlabel(str)
char *str;
{
	XDrawImageString(mydisplay,mywindow,mygc,10,50,str,strlen(str));
}








/*
*			ROUTINES FOR COLOR CONTROL
*/


#define MAX_COLOR_TABLES 10
static int color_table_created[MAX_COLOR_TABLES];

#define PIXELS 256
static unsigned long bwpixels[2];
static unsigned long *pixels[MAX_COLOR_TABLES];
static int num_colors[MAX_COLOR_TABLES];


/*
*			Routines for Creating Color Tables:
*/

static Colormap cmap;
static Display *display;
static int cells; /* Number of color cells */
#define PLANES 8
#define MAX_COLORS 200

static void
init_colors(cur_display,screen,depth,myforeground,mybackground)
	Display *cur_display;
	int screen,depth;
	unsigned long myforeground,mybackground;

{

	int i;

	display = cur_display;
	if (depth==1) {
		bwpixels[0] = mybackground;
		bwpixels[1] = myforeground;
		return;
	}

	cmap = XDefaultColormap(display,screen);
	cells = XDisplayCells(display,screen);
	if (debug) printf("cmap = %d cells = %d\n",cmap,cells);
	for (i=0; i<MAX_COLOR_TABLES; i++) color_table_created[i] = 0;
}


static int tot_num_colors = 0;

static unsigned long plane_masks[MAX_COLOR_TABLES][PLANES];

int available_colors()
{
	static Bool contig = False;
	static unsigned int nplanes = 0;
	unsigned long *Pixels;
	int result;
	int new_colors = PIXELS;
	unsigned long plane_mask = 0;

	Pixels = (unsigned long *)malloc(new_colors*sizeof (unsigned
long));
	if (Pixels==NULL) {
		fprintf(stderr,"Unable to allocate pixels wth
malloc()\n");
		return -1;
	}

	while (new_colors>0 && (result =
XAllocColorCells(display,cmap,contig,plane_masks[0],
		nplanes,Pixels,new_colors)) == 0) {
		if (debug) printf("XAllocColorCells(ncolors=%d) failed,
decreasing
number\n",new_colors);
		new_colors--;
	}
	if (debug) printf("Global color table has room for %d
colors\n",new_colors);
	if (new_colors)
XFreeColors(display,cmap,Pixels,new_colors,plane_mask);
	free(Pixels);

	return new_colors;
}




int create_color_table(table,new_colors)
int table,new_colors;
{
    /* MM: Not needed anymore...
	static Bool contig = False;
	static unsigned int nplanes = 0;
	int result;
    */
	unsigned long *Pixels;


	if (table < 0 || table >= MAX_COLOR_TABLES) {
		fprintf(stderr,"Color Tables must be numbered from 0 to
%d\n",
				MAX_COLOR_TABLES-1);
		return -1;
	}
	if (new_colors>PIXELS) {
		fprintf(stderr,"Color tables cannot have more than %d
entries\n"
			,PIXELS);
		return -2;
	}

	if (color_table_created[table]) {
		fprintf(stderr,"Color table %d already exists\n",table);
		return -3;
	}
	tot_num_colors += new_colors;
	Pixels = (unsigned long *)malloc(new_colors*sizeof (unsigned
long));
	if (Pixels==NULL) {
		fprintf(stderr,"Unable to allocate pixels wth
malloc()\n");
		return -4;
	}

/******************************************************************
**  MM: HERE WAS THE BUG: (new_colors>0) must be in front. !! 
**  This whole thing doesn't work anymore if You have a
**  X11-Server in TrueColor-Mode (= read only colormap..)

	while (new_colors>0 && (result =
XAllocColorCells(display,cmap,contig,plane_masks[table],
		nplanes,Pixels,new_colors)) == 0) {
		if (debug) printf("XAllocColorCells(ncolors=%d) failed,
decreasing
number\n",new_colors);
		new_colors--;
	}
	if (debug) printf("Color table %d created with %d
colors\n",table,new_colors);
	if (new_colors==0) {
		free(Pixels);
		return 0;
	}
******************************************************************/
	pixels[table] = Pixels;


	color_table_created[table] = 1;
	num_colors[table] = new_colors;

	return new_colors;
}



int free_color_table(table)
int table;
{
	unsigned long plane_mask = 0;

	if (table < 0 || table >= MAX_COLOR_TABLES) {
		fprintf(stderr,"Color Tables must be numbered from 0 to
%d\n",
				MAX_COLOR_TABLES-1);
		return 0;
	}
	if (!color_table_created[table]) {
		fprintf(stderr,"Color Table %d does not exist\n",table);
		return 0;
	}
	
XFreeColors(display,cmap,pixels[table],num_colors[table],plane_mask);
	free(pixels[table]);
	color_table_created[table] = 0;
	tot_num_colors -= num_colors[table];
        return 1; /* MM: hope this is correct */
}



static int good_table(table)
	int table;
{
	if (table < 0) return 0;
	if (table >= MAX_COLOR_TABLES) return 0;
	if (!color_table_created[table]) return 0;
	return 1;
}




/*
*			Routines for Setting Up Colors:
*/



int setup_table_color(table,color,red,green,blue)
int table,color;
unsigned red,green,blue;
{
    /*	int i,c; */
	static int flags = DoRed | DoGreen | DoBlue;
	static XColor def;

	if (!good_table(table)) {
		fprintf(stderr,"No color table %d exists\n",table);
		return 0;
	}
	if (color<0 || color >= num_colors[table]) {
		fprintf(stderr,"Attempt to set color outside range
[0,%d]\n",num_colors[table]);
		return 0;
	}
	def.red = red;
	def.green = green;
	def.blue = blue;
	def.pixel = pixels[table][color];
	def.flags = flags;
/*
*     Now store the arrays we created above into the color map.
*/

/* MM - here is the read-only version which takes the closest color */

        if (XAllocColor(display,cmap,&def)==0)
             printf("Problem wth. XAllocColor\n");
	pixels[table][color] = def.pixel;

/********  MM - doesn't work:

	XStoreColor(display,cmap,&def);
**************************************/
	return 1;
}






int setup_table_colors(table,ncolors,start,red,green,blue)
int table,ncolors,start;
unsigned *red,*green,*blue;
{
	int i,c;
	static int flags = DoRed | DoGreen | DoBlue;
	static XColor defs[PIXELS];

	if (!good_table(table)) {
		fprintf(stderr,"No color table %d exists\n",table);
		return 0;
	}
	if (start<0 || start+ncolors-1 >= num_colors[table]) {
		fprintf(stderr,"Attempt to set color outside range
[0,%d]\n",num_colors[table]-1);
		return 0;
	}
	for (i=0; i<ncolors; i++) {
		c = i + start;
		defs[i].red = red[i];
		defs[i].green = green[i];
		defs[i].blue = blue[i];
		defs[i].pixel = pixels[table][c];
		defs[i].flags = flags;

/*** MM: new code ... */
                if (XAllocColor(display,cmap,&(defs[i]))==0)
                      printf("XAllocColor failed !!\n");
                      pixels[table][c] = defs[i].pixel;
/*************************************************/                     
	}
/*
*     Now store the arrays we created above into the color map.
*/
/**** MM doesn't work :
	XStoreColors(display,cmap,defs,ncolors);
*********************************/
	return 1;
}



/*
*			Routines for Using Colors:
*/

int default_table = 0;
int default_palette = 1;

void set_default_color_table(table)
int table;
{
	default_table = table;
}




void set_default_color_palette(table)
int table;
{
	default_palette = table;
}




void _enquire_colors(numt,nump)
	int *numt,*nump;
{
	*numt = num_colors[default_table];
	*nump = num_colors[default_palette];
}




int set_color_from_Table(table,col)
	int table,col;
{
	unsigned long pix;

	if (debug) printf("set_color_from_Table(%d,%d)\n",table,col);
	if (col<0 || col>=num_colors[table]) return 0;
	pix = pixels[table][col];

	XSetForeground(mydisplay,pixgc,pix);
	return 1;
}


int set_color_from_table(col)
	int col;
{
	return set_color_from_Table(default_table,col);
}




int set_color_from_Palette(table,pal)
	int table;
	double pal;
{
	int cur_color;

	if (pal<0.0 || pal>1.0) return 0;
	cur_color = pal*(num_colors[table]-.001);
	XSetForeground(mydisplay,pixgc,pixels[table][cur_color]);
	return 1;
}


int set_color_from_palette(pal)
	double pal;
{
	return set_color_from_Palette(default_palette,pal);
}



/*
*		Routines for converting between color representations:
*/

int color_from_Palette(table,pal,col)
int table,*col;
double pal;
{
	if (!good_table(table)) return 0;
	if (pal<0.0 || pal>1.0) return 0;
	*col = (num_colors[table]-.0001)*pal;
	return 1;
}




int Table_to_Global(table,color,gcolor)
int table,color,*gcolor;
{
	int i;

	if (!good_table(table) || color<0 || color>=num_colors[table])
return 0;
	*gcolor = 0;
	for (i=0; i<table; i++) {
		if (!color_table_created[i]) continue;
		*gcolor += num_colors[i];
	}
	*gcolor += color;
	return 1;
}


int Global_to_Table(gcolor,table,color)
int gcolor,*table,*color;
{
	int i;

	if (gcolor<0 || gcolor >= tot_num_colors) return 0;
	for (i=0; i<MAX_COLOR_TABLES; i++) {
		if (!color_table_created[i]) continue;
		if (gcolor < num_colors[i]) break;
		gcolor -= num_colors[i];
	}
	if (i==MAX_COLOR_TABLES || !color_table_created[i]) return 0;
	*table = i;
	*color = gcolor;
	return 1;
}

int Palette_to_Global(table,pal,gcol)
int table,*gcol;
double pal;
{
	int col;

	if (!color_from_Palette(table,pal,&col)) return 0;
	if (col<0) return 0;
	return Table_to_Global(table,col,gcol);
}





/*
*			Routines for Controlling Intensity:
*/

static int numi;
static int mini, maxi; /* X11 uses 65535 for maxi */

void set_linear_intensity(Numi,Mini,Maxi)
	int Numi,Mini,Maxi;
{
	numi = Numi;
	mini = Mini;
	maxi = Maxi;
}




	/* mini <= linear_intensity(i) <= maxi, as 0 <= i <= numi-1 */

int linear_intensity(i)
	int i;
{
	if (i<0 || i>=numi) return 0;
	return mini + (i*(maxi-mini))/(numi-1);
}



static double r_eye;
#include <math.h>

void set_eye_intensity(Numi,Mini,Maxi)
	int Numi,Mini,Maxi;
{
	numi = Numi;
	mini = Mini;
	maxi = Maxi;
		/* mini * r**(numi-1) = maxi */
	r_eye = pow( (double)maxi/(double)mini, 1.0/(double)(numi-1) );
}





	/* mini <= eye_intensity(i) <= maxi, as 0 <= i <= numi-1 */

int eye_intensity(i)
	int i;
{
	if (i<0 || i>=numi) return 0;
	return (int) (mini*pow(r_eye,(double)i) + .00001);
}



/*
*			END OF FRAME AND ERASING, CLEANUP
*/



/*
*
*				end_frame():
*
*	This actually draws the picture
*	Should be called by end_frame() in lib2d.
*/

void _end_frame()
{
	do_plot();
	if (debug) printf("do_plot done\n");
	flush_graphics(); /* Must flush to get the complete picture on
screen */
	if (debug) printf("flush_graphics done\n");
}





int wait_for_prompt = 1;

/*
*				do_plot():
*
*	This routine actually draws the plots.  All other routines
simply
*	draw into the Pixmap.  This routines watches for an appropriate
event
*	such has having on overlaid window removed (Expose) or a mouse
button
*	pressed, and then dumps the Pixmap to the window.
*/

static void do_plot()
{
	char *prompt =
	"Type c (or mouse button) to continue, a for animation or q to
quit:
";
	int i;
	static int movie = 0;
	static char text[100];
	static XEvent myevent;
	static KeySym mykey;
	int num_events;
	static int done = 0;



	if (wait_for_prompt) {
		set_color_from_table(1);
	
XDrawImageString(mydisplay,SURFACE,pixgc,50,ny-50,prompt,strlen(prompt))
;
		if (debug) printf("XDrawImageString done\n");
	}

		/* Here is where we actually draw on the screen: */
	start:
	if (SURFACE == pixmap)
	
XCopyArea(mydisplay,pixmap,mywindow,pixgc,0,0,nx,ny,0,0);
	if (debug) printf("XCopyArea done\n");
	XFlush(mydisplay);
	if (debug) printf("XFlush done\n");

	while (!done) {
		num_events =
XEventsQueued(mydisplay,QueuedAfterReading);
		if (debug && num_events) printf("XEventsQueued found %d
events\n",num_events);
		if (num_events==0) {
			XFlush(mydisplay);
			check_pixmap(mydisplay,pixmap);
			if (wait_for_prompt) continue;
			return;
		}

			/* Get the Event: */
		XNextEvent(mydisplay,&myevent);
		switch(myevent.type) {
		case Expose:
			if (myevent.xexpose.count==0) goto start;
			break;

		case MappingNotify:
		        XRefreshKeyboardMapping((XMappingEvent
*)&myevent); 
			break;                                

		case ButtonPress:
			_erase();
			return;
			break;

		case KeyPress:
			_erase();
			if (movie) return;
			i = XLookupString((XKeyEvent
*)&myevent,text,10,&mykey,0);
			if (i==1 && text[0]=='q') {
				x_close_graphics();
				exit(0);
			}
			if (i==1 && text[0]=='c') return;
			if (i==1 && text[0]=='a') {
				movie = 1;
				return;
			}
			break;

		default:
			/* printf("This is an unknown event type\n"); */
			break;
		}
		if (movie) return;
		continue;

	}
}






static void check_pixmap(display,pixmap)
Display *display;
Pixmap pixmap; 
{
	Window root;
	int x,y;
	unsigned int width, height,border_width;
	int depth;
	Status status;

	status = XGetGeometry(display,pixmap,&root,&x,&y,
			&width,&height,&border_width,&depth);
	if (status!=1) {
		fprintf(stderr,"status=0: Failed to create pixmap\n");
		exit(1);
	}
}




/*
*				_erase():
*
*	Erases the Pixnap or the screen
*/

void _erase()
{
	set_color_from_table(0);
	if (SURFACE == pixmap)
XFillRectangle(mydisplay,pixmap,pixgc,0,0,nx,ny);
	else XClearWindow(mydisplay,mywindow);
}







/*
*				_close_graphics():
*
*	Call _close_graphics() from lib2d close_graphics() to
*	cleanup and restore screen
*/

void _close_graphics()
{
	x_close_graphics();
}

static void x_close_graphics()
{
	XFreeGC(mydisplay,mygc);
	XDestroyWindow(mydisplay,mywindow);
	XCloseDisplay(mydisplay);
	exit(0);
}




/*
*				flush_graphics():
*
*	Flushes any graphics from buffers onto the screen.
*/

void flush_graphics()
{
	XFlush(mydisplay);
}
--

Previous <http://www.cs.colorado.edu/~mcbryan/5229.03/mail/35.htm>
Next <http://www.cs.colorado.edu/~mcbryan/5229.03/mail/37.htm> 

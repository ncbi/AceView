/*  Last edited: Dec 12 19:37 1994 (mieg) */

/* $Id: plot.h,v 1.2 2011/09/22 16:38:47 mieg Exp $ */

void plotHisto (char *title, Array a) ;
/* mhmp 10.03.99 float xmin*/
void plotShiftedHisto (char *title, char *subtitle, Array a, float xmin, float xmax, float scale, int xfac) ;
void plotHistoRemove (void) ;
/* mhmp 11.05.01 */
int regular (int p) ;
void plot2d (char *title, char *subtitleX, char *subtitleY, Array a) ;

/**************/


/**********************************************************/
/**************	  GraphImage Code CLH	*******************/
/****** As described by RD, CLH, SEL, FEH 3/29/94 *********/
/**********************************************************/

/* $Id: graphimage.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

#ifndef _GRAPH_IMAGE_H_
#define _GRAPH_IMAGE_H_

enum GraphImageType { Grey8, Colour8, Colour32 } ;
enum GraphImageFileType { kNoFile, kGraphPICT } ;

struct GraphImageStruct {
	void*		imagePtr ;
	int		imageType ;
	int		fileType ;
	int		width ;
	int 		height ;
	float		top;
	float		left;
	float		bottom;
	float		right;
	} ;
typedef struct GraphImageStruct GraphImageStruct,  *GraphImagePtr;

GraphImageStruct *graphImageRead (char *filename) ;		/* reads from file */
GraphImageStruct *graphImageCreate (int type, int width, int height) ;
void graphImageDestroy (GraphImageStruct *gim) ;

/*	Pushes the image on to the drawing stack	*/
void graphImageDraw (GraphImageStruct *gim, float x0, float y0, float x1, float y1);

GraphImageStruct *graphImageResize (GraphImageStruct *gim, float width, float height) ;

BOOL graphImageGetData (GraphImageStruct *gim, unsigned char *data, int type,
                        int x, int y, int w, int h, int lineLength) ;
				/* true color space - lineLength in pixels */
BOOL graphImageSetData (GraphImageStruct *gim,  unsigned char *data, int type,
                        int x, int y, int w, int h, int lineLength) ;

BOOL graphImageRamp (GraphImageStruct *gim, int min, int max) ;
				/* contrast, brightness - can fail, e.g. for Colour32 for now! */


#endif

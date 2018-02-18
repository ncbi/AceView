/*  Last edited: Oct 13 13:00 1993 (ulrich) */

  /* topology.h */

/* $Id: topology.h,v 1.2 2009/05/20 17:25:53 mieg Exp $ */

typedef struct { KEY a , b ; int group , type ; } LINK ;
#define linkFormat "kkii"

typedef struct { KEY a ; int group, pos ; } VERTEX ;
#define vertexFormat "kii"

int topoVertexOrder (const void *v1, const  void *v2) ;
int topoLinkOrder (const void *v1, const  void *v2) ;
int topoConnectedComponents(Array links, Array vx) ;

void chronologicalOrdering(Array mm, Array cocolO, Array lilinO, Stack names) ;



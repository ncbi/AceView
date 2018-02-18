/*  Last edited: Nov 27 19:01 1995 (mieg) */
/*********************************************************
  myNetwork.h
  --------------------------------------------------------
  generated at Tue Jul 25 14:26:05 1995
  by snns2c ( Bernward Kett 1995 ) 
*********************************************************/

/* $Id: myNetwork.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

extern int myNetwork(float *in, float *out, int init);
#ifdef JUNK
static struct {
  int NoOfInput;    /* Number of Input Units  */
  int NoOfOutput;   /* Number of Output Units */
  int(* propFunc)(float *, float*, int);
} myNetworkREC = {240,4,myNetwork};
#endif

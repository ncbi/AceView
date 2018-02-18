/*  Last edited: Oct 14 10:33 1995 (mieg) */
/*********************************************************
  myNetwork.c
  --------------------------------------------------------
  generated at Fri Aug 25 14:16:49 1995
  by snns2c ( Bernward Kett 1995 ) 
*********************************************************/
/* @(#)myNetwork.c	1.2    3/20/97 */

#include <math.h>
#include "myNetwork.h"

#define Act_Logistic(sum, bias)     ( 1.0/(1 + exp(-sum-bias) ) )

#ifdef JUNK
mieg: feb 98 what is this for ?
static struct {
  int NoOfInput;    /* Number of Input Units  */
  int NoOfOutput;   /* Number of Output Units */
  int(* propFunc)(float *, float*, int);
} myNetworkREC = {240,4,myNetwork};

#endif

int myNetwork(float *in, float *out, int init)
{
  int member, source;
  float sum;
  enum{OK, Error, Not_Valid};

typedef struct UT {
          float act;         /* Activation       */
          float Bias;        /* Bias of the Unit */
          int   NoOfSources; /* Number of predecessor units */
   struct UT   *sources[144]; /* predecessor units */
          float weights[144]; /* weights from predecessor units */
        } UnitType, *pUnit;

  pUnit unit;

  /* unit definition section (see also UnitType) */
  static UnitType Units[224] = 
  {
    { 0.0, 0.0, 0, {0 /* NO SOURCES */}, {0.0 /* NO MEMBERS*/} },
    { /* unit 1 (Old: 1) */
      0.0, 0.041410, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 2 (Old: 2) */
      0.0, -0.136190, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 3 (Old: 3) */
      0.0, 0.058670, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 4 (Old: 4) */
      0.0, 0.021370, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 5 (Old: 5) */
      0.0, 0.072520, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 6 (Old: 6) */
      0.0, -0.154570, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 7 (Old: 7) */
      0.0, 0.193770, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 8 (Old: 8) */
      0.0, -0.033640, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 9 (Old: 9) */
      0.0, 0.136250, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 10 (Old: 10) */
      0.0, 0.046510, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 11 (Old: 11) */
      0.0, -0.076400, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 12 (Old: 12) */
      0.0, 0.176460, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 13 (Old: 13) */
      0.0, -0.159940, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 14 (Old: 14) */
      0.0, 0.134580, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 15 (Old: 15) */
      0.0, 0.136370, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 16 (Old: 16) */
      0.0, -0.013230, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 17 (Old: 17) */
      0.0, -0.041660, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 18 (Old: 18) */
      0.0, -0.033080, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 19 (Old: 19) */
      0.0, 0.092010, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 20 (Old: 20) */
      0.0, 0.043810, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 21 (Old: 21) */
      0.0, 0.082640, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 22 (Old: 22) */
      0.0, -0.096950, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 23 (Old: 23) */
      0.0, 0.080590, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 24 (Old: 24) */
      0.0, 0.169780, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 25 (Old: 25) */
      0.0, 0.038010, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 26 (Old: 26) */
      0.0, -0.142950, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 27 (Old: 27) */
      0.0, -0.176790, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 28 (Old: 28) */
      0.0, -0.065130, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 29 (Old: 29) */
      0.0, -0.138590, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 30 (Old: 30) */
      0.0, 0.198900, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 31 (Old: 31) */
      0.0, 0.015050, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 32 (Old: 32) */
      0.0, -0.013040, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 33 (Old: 33) */
      0.0, -0.115150, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 34 (Old: 34) */
      0.0, 0.093760, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 35 (Old: 35) */
      0.0, -0.193100, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 36 (Old: 36) */
      0.0, 0.077290, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 37 (Old: 37) */
      0.0, -0.040340, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 38 (Old: 38) */
      0.0, -0.043490, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 39 (Old: 39) */
      0.0, 0.115020, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 40 (Old: 40) */
      0.0, -0.154360, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 41 (Old: 41) */
      0.0, 0.078140, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 42 (Old: 42) */
      0.0, 0.139260, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 43 (Old: 43) */
      0.0, 0.064940, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 44 (Old: 44) */
      0.0, 0.045010, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 45 (Old: 45) */
      0.0, -0.057440, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 46 (Old: 46) */
      0.0, -0.101420, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 47 (Old: 47) */
      0.0, -0.041450, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 48 (Old: 48) */
      0.0, -0.012650, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 49 (Old: 49) */
      0.0, 0.016260, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 50 (Old: 50) */
      0.0, -0.061000, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 51 (Old: 51) */
      0.0, 0.069130, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 52 (Old: 52) */
      0.0, -0.178550, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 53 (Old: 53) */
      0.0, 0.052780, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 54 (Old: 54) */
      0.0, -0.177560, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 55 (Old: 55) */
      0.0, 0.197030, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 56 (Old: 56) */
      0.0, -0.006640, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 57 (Old: 57) */
      0.0, 0.090890, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 58 (Old: 58) */
      0.0, 0.190280, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 59 (Old: 59) */
      0.0, -0.036780, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 60 (Old: 60) */
      0.0, 0.118010, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 61 (Old: 61) */
      0.0, -0.151080, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 62 (Old: 62) */
      0.0, 0.176250, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 63 (Old: 63) */
      0.0, 0.095660, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 64 (Old: 64) */
      0.0, 0.078870, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 65 (Old: 65) */
      0.0, -0.156600, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 66 (Old: 66) */
      0.0, 0.000720, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 67 (Old: 67) */
      0.0, -0.084010, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 68 (Old: 68) */
      0.0, 0.085430, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 69 (Old: 69) */
      0.0, -0.145970, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 70 (Old: 70) */
      0.0, -0.070220, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 71 (Old: 71) */
      0.0, 0.016600, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 72 (Old: 72) */
      0.0, -0.183850, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 73 (Old: 73) */
      0.0, -0.109870, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 74 (Old: 74) */
      0.0, 0.049380, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 75 (Old: 75) */
      0.0, 0.108540, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 76 (Old: 76) */
      0.0, 0.058190, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 77 (Old: 77) */
      0.0, 0.079870, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 78 (Old: 78) */
      0.0, -0.067910, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 79 (Old: 79) */
      0.0, -0.087590, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 80 (Old: 80) */
      0.0, -0.026380, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 81 (Old: 81) */
      0.0, -0.129790, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 82 (Old: 82) */
      0.0, 0.043760, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 83 (Old: 83) */
      0.0, -0.127510, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 84 (Old: 84) */
      0.0, -0.137600, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 85 (Old: 85) */
      0.0, 0.127810, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 86 (Old: 86) */
      0.0, -0.177360, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 87 (Old: 87) */
      0.0, 0.030050, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 88 (Old: 88) */
      0.0, -0.008270, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 89 (Old: 89) */
      0.0, 0.173740, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 90 (Old: 90) */
      0.0, -0.165400, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 91 (Old: 91) */
      0.0, -0.153030, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 92 (Old: 92) */
      0.0, -0.104550, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 93 (Old: 93) */
      0.0, 0.040430, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 94 (Old: 94) */
      0.0, -0.075300, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 95 (Old: 95) */
      0.0, -0.104620, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 96 (Old: 96) */
      0.0, 0.038000, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 97 (Old: 97) */
      0.0, 0.149900, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 98 (Old: 98) */
      0.0, 0.006150, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 99 (Old: 99) */
      0.0, 0.111010, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 100 (Old: 100) */
      0.0, -0.149250, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 101 (Old: 101) */
      0.0, -0.011530, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 102 (Old: 102) */
      0.0, 0.199430, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 103 (Old: 103) */
      0.0, -0.144210, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 104 (Old: 104) */
      0.0, 0.192520, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 105 (Old: 105) */
      0.0, -0.125960, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 106 (Old: 106) */
      0.0, 0.102850, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 107 (Old: 107) */
      0.0, 0.074170, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 108 (Old: 108) */
      0.0, -0.186290, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 109 (Old: 109) */
      0.0, -0.174220, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 110 (Old: 110) */
      0.0, -0.123770, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 111 (Old: 111) */
      0.0, 0.003160, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 112 (Old: 112) */
      0.0, 0.111950, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 113 (Old: 113) */
      0.0, -0.030540, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 114 (Old: 114) */
      0.0, 0.084390, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 115 (Old: 115) */
      0.0, 0.071570, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 116 (Old: 116) */
      0.0, 0.095470, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 117 (Old: 117) */
      0.0, 0.130400, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 118 (Old: 118) */
      0.0, 0.199270, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 119 (Old: 119) */
      0.0, 0.182090, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 120 (Old: 120) */
      0.0, 0.103530, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 121 (Old: 121) */
      0.0, 0.033820, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 122 (Old: 122) */
      0.0, -0.080650, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 123 (Old: 123) */
      0.0, 0.111260, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 124 (Old: 124) */
      0.0, -0.001560, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 125 (Old: 125) */
      0.0, 0.173190, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 126 (Old: 126) */
      0.0, 0.042770, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 127 (Old: 127) */
      0.0, 0.008210, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 128 (Old: 128) */
      0.0, 0.112740, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 129 (Old: 129) */
      0.0, 0.112200, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 130 (Old: 130) */
      0.0, -0.166480, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 131 (Old: 131) */
      0.0, 0.059910, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 132 (Old: 132) */
      0.0, 0.122920, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 133 (Old: 133) */
      0.0, 0.115510, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 134 (Old: 134) */
      0.0, -0.053470, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 135 (Old: 135) */
      0.0, 0.178570, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 136 (Old: 136) */
      0.0, -0.113360, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 137 (Old: 137) */
      0.0, 0.187730, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 138 (Old: 138) */
      0.0, 0.022360, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 139 (Old: 139) */
      0.0, 0.129430, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 140 (Old: 140) */
      0.0, -0.172870, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 141 (Old: 141) */
      0.0, -0.163860, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 142 (Old: 142) */
      0.0, 0.010860, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 143 (Old: 143) */
      0.0, -0.148680, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 144 (Old: 144) */
      0.0, -0.078230, 0,
      {0 /* no source units */},
      {0.0 /* NO MEMBERS */}
    },
    { /* unit 145 (Old: 145) */
      0.0, -0.564390, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.113570, 0.068300, -0.262320, 0.293560, -0.156540, 
-0.031620, -0.188520, 0.079240, 0.051260, 0.135310, -0.279760, 
-0.089350, 0.315580, 0.042600, -0.059070, 0.041160, 0.219180, 
-0.155250, 0.225690, -0.275600, 0.369160, -0.236390, 0.169150, 
-0.327450, 0.192260, -0.077390, 0.163310, -0.256870, 0.120640, 
-0.002370, 0.056500, -0.392640, 0.048820, 0.075880, -0.164810, 
-0.181660, -0.326570, 0.190090, -0.300270, -0.397840, -0.375540, 
0.130430, -0.067940, -0.423100, -0.119070, 0.215750, -0.339940, 
-0.360650, -0.147270, 0.102460, 0.016250, -0.325600, -0.165710, 
-0.018970, 0.022810, -0.241660, -0.021890, -0.086590, -0.133730, 
0.016960, -0.204780, 0.250280, 0.032620, 0.081770, 0.111790, 
0.036740, 0.163400, -0.101210, 0.095040, -0.062440, 0.124320, 
0.173730, 0.191510, 0.176380, 0.188800, 0.137050, 0.344990, 
-0.010090, 0.027790, 0.233100, 0.009490, 0.251230, 0.077120, 
-0.050950, 0.031670, -0.069330, -0.245970, 0.221950, 0.248100, 
0.108920, -0.094320, 0.110790, 0.258500, 0.031890, -0.283550, 
0.287450, 0.241380, -0.018730, -0.112330, 0.120730, 0.005690, 
-0.021730, 0.111940, 0.088240, -0.101080, -0.196600, 0.076340, 
-0.180370, 0.093710, -0.067070, -0.054020, -0.042360, -0.261660, 
-0.158720, 0.073980, -0.246600, -0.116250, -0.168210, 0.254750, 
0.063810, -0.146580, -0.130320, 0.318190, -0.058010, 0.022430, 
-0.130810, 0.200550, 0.185670, -0.067760, 0.088230, -0.107400, 
0.086240, -0.084450, 0.080140, -0.070200, 0.119520, -0.041810, 
-0.062820, -0.228540, 0.203920, 0.083820, 0.311750, -0.124490, 
0.271690}
    },
    { /* unit 146 (Old: 146) */
      0.0, -0.642670, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.214690, -0.476250, 0.000320, 0.050550, -0.280720, 
-0.356650, 0.086800, 0.059320, 0.008630, -0.406250, -0.091250, 
-0.097830, 0.037730, -0.262830, 0.206820, 0.218810, -0.009830, 
-0.282030, -0.086500, 0.174930, 0.042160, -0.166880, 0.181140, 
0.307580, -0.063290, -0.080730, 0.212710, 0.040130, -0.212350, 
0.099330, 0.010760, 0.197600, -0.158800, 0.262740, -0.078820, 
0.023990, -0.112690, -0.004710, -0.245370, -0.042400, -0.088380, 
0.111870, -0.294790, -0.098050, -0.117800, -0.133210, -0.059790, 
0.168860, -0.173390, 0.191530, -0.111310, 0.098840, 0.022260, 
0.070190, 0.024830, 0.075940, -0.015310, -0.138600, 0.156300, 
0.061460, 0.074270, 0.059720, 0.079160, 0.076590, 0.390300, 
-0.053130, 0.358550, 0.197390, 0.466600, 0.217350, 0.150540, 
0.230700, 0.428370, 0.110460, 0.451960, 0.316990, 0.661580, 0.490090, 
0.436050, 0.388000, 0.502570, 0.491010, 0.438760, 0.521450, 0.104480, 
0.214430, 0.266670, 0.339350, -0.087270, 0.147020, 0.192280, 
0.318610, 0.114040, 0.085990, 0.076390, 0.008310, 0.060790, 0.092360, 
0.085980, 0.104940, 0.018930, 0.188380, 0.122230, -0.295380, 
0.102770, 0.110790, 0.340730, -0.119810, -0.207400, 0.016690, 
0.145850, -0.224020, 0.107480, -0.278380, 0.390500, -0.244610, 
-0.203390, 0.009340, 0.497770, -0.180710, -0.227830, 0.134700, 
0.375290, -0.402080, 0.010030, -0.100310, 0.307680, -0.266340, 
-0.064990, 0.147600, 0.115550, -0.290680, 0.021460, 0.300440, 
0.238630, -0.068430, -0.231200, 0.263410, -0.073160, -0.227650, 
-0.183100, 0.345050, 0.079630, -0.037230}
    },
    { /* unit 147 (Old: 147) */
      0.0, -0.361710, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.112730, 0.076880, 0.260000, -0.069090, -0.099400, -0.108900, 
-0.112200, 0.175820, -0.174820, -0.077470, -0.225770, 0.084160, 
-0.280620, -0.077020, -0.264450, -0.037450, 0.023220, -0.256220, 
-0.424030, 0.310700, 0.155990, -0.061640, -0.431410, 0.225380, 
-0.082270, -0.303070, -0.087280, -0.024190, 0.273050, 0.049570, 
-0.305390, 0.060300, 0.176940, -0.166490, -0.251540, -0.011200, 
-0.105070, -0.015050, -0.260490, 0.115260, 0.048340, -0.001570, 
-0.179730, 0.106730, -0.201820, 0.186700, -0.207030, -0.113980, 
-0.095830, -0.011200, 0.187960, 0.077100, -0.093390, -0.116250, 
0.053030, -0.000810, -0.217080, 0.012800, 0.296390, 0.211480, 
-0.180150, -0.285540, 0.180420, -0.019550, -0.017440, -0.169340, 
0.086810, 0.279690, -0.268690, -0.160050, 0.248250, 0.101810, 
0.034170, -0.086970, 0.168400, 0.216980, -0.159510, -0.079510, 
-0.119240, 0.007590, -0.008500, -0.165020, -0.156670, 0.088070, 
0.099300, -0.327480, -0.196210, 0.227490, -0.090290, -0.191110, 
-0.058990, -0.006600, 0.156860, 0.021810, -0.241390, -0.075380, 
-0.075510, -0.201890, -0.266420, -0.114960, 0.011050, -0.150520, 
0.033370, -0.051690, 0.187010, 0.019570, -0.071250, 0.013110, 
-0.038010, -0.106700, -0.247510, -0.066960, -0.209220, 0.229510, 
-0.254680, -0.126540, -0.027010, 0.301390, -0.097520, -0.181040, 
-0.058670, 0.182940, -0.137910, 0.158760, -0.242030, 0.241920, 
-0.303070, -0.088650, -0.025260, -0.041950, -0.113320, 0.090450, 
0.157470, -0.058740, -0.139230, 0.114880, 0.034590, -0.163820, 
-0.058740, -0.058710, 0.303460, -0.019930, -0.300190, 0.156770}
    },
    { /* unit 148 (Old: 148) */
      0.0, -0.327750, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.738280, -0.608760, 0.200330, 0.433900, -0.830440, 
-0.411310, 0.278050, 0.025810, -0.565740, -0.215440, 0.055880, 
-0.027150, -0.333110, -0.373960, 0.363410, 0.086040, 0.074360, 
-0.647590, 0.078600, 0.187340, 0.291770, -0.651910, 0.165280, 
-0.233660, 0.161270, -0.604270, 0.188150, -0.247730, 0.125650, 
-0.498650, 0.059620, -0.244110, -0.140070, -0.307030, 0.195410, 
-0.346290, -0.362600, -0.590530, 0.203070, -0.436710, -0.245620, 
-0.498160, 0.212920, -0.194120, -0.308980, -0.303920, 0.180110, 
-0.063950, -0.535910, -0.420290, 0.234510, -0.097940, -0.582570, 
-0.546650, -0.011960, 0.147410, -0.511130, -0.350120, 0.385560, 
-0.099310, -0.150890, -0.183430, 0.378010, -0.007040, 0.237870, 
-0.391580, 0.327770, 0.018800, 0.419270, -0.302380, 0.170780, 
0.099740, 0.305580, -0.067820, 0.057850, -0.062550, 0.402890, 
0.214450, 0.303270, -0.048370, 0.434420, 0.023060, -0.049780, 
0.169110, 0.384760, 0.132370, -0.204980, -0.173590, 0.468250, 
-0.007180, -0.166050, 0.135470, 0.219920, -0.168260, -0.070130, 
0.088750, 0.127680, -0.182510, -0.218790, 0.237410, 0.220570, 
-0.171020, -0.251880, -0.044140, -0.087560, 0.047090, -0.211550, 
0.136780, -0.036280, -0.135760, -0.221330, 0.329970, -0.388910, 
0.015500, 0.030670, 0.250860, -0.279100, 0.047890, 0.273230, 
0.372960, -0.555750, 0.042820, 0.259760, 0.422690, -0.354010, 
0.334250, 0.351900, 0.517160, -0.643130, 0.173530, 0.357860, 
0.311280, -0.457930, 0.130170, 0.119700, 0.302330, -0.596970, 
0.285840, 0.018740, 0.006450, -0.636910, -0.103540, 0.058010, 
0.071200}
    },
    { /* unit 149 (Old: 149) */
      0.0, -0.277560, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.022150, -0.331970, -0.454630, 0.224870, 0.108900, -0.047940, 
-0.308710, 0.192030, 0.004060, -0.128990, -0.028880, 0.262910, 
0.288670, 0.069800, -0.090660, 0.303420, 0.213680, -0.149390, 
-0.013830, -0.111880, -0.113790, -0.114830, -0.277090, 0.001330, 
-0.204650, -0.293060, -0.049860, -0.135390, -0.283370, -0.420040, 
-0.434810, -0.214040, -0.066210, -0.202730, -0.284580, -0.057890, 
-0.317990, -0.416200, -0.206930, -0.059360, -0.229970, -0.206940, 
-0.003750, -0.012040, -0.162590, -0.291020, -0.271640, 0.161990, 
0.063120, -0.136620, -0.145990, -0.032350, -0.035460, 0.234250, 
-0.132170, -0.036980, 0.051230, 0.374260, -0.139900, 0.017150, 
-0.293300, 0.363300, -0.403890, -0.152660, -0.288570, 0.709540, 
-0.203890, -0.337220, -0.444920, 0.440580, -0.318950, -0.064620, 
-0.139100, 0.723960, 0.025450, -0.122430, -0.296720, 0.563330, 
-0.136580, -0.159180, -0.109120, 0.611900, 0.095300, -0.038900, 
-0.233540, 0.517120, 0.185150, -0.280750, -0.248230, 0.598590, 
0.374200, -0.109800, -0.549190, 0.400220, 0.561260, -0.211440, 
-0.397850, 0.573600, 0.324450, 0.002620, -0.349880, 0.323830, 
0.265630, -0.064830, -0.338450, 0.343180, -0.023330, 0.163710, 
-0.360340, 0.188860, 0.080190, -0.090340, -0.108020, -0.103680, 
-0.287690, -0.137390, -0.283380, -0.216710, -0.013580, -0.117400, 
-0.008010, -0.557130, -0.274940, 0.036320, 0.062540, -0.414480, 
-0.221910, 0.007050, -0.293020, -0.646070, -0.426700, -0.181520, 
-0.205630, -0.578850, -0.417700, -0.226430, -0.305460, -0.648140, 
-0.046450, -0.039690, -0.296130, -0.499390, 0.142610, 0.197460}
    },
    { /* unit 150 (Old: 150) */
      0.0, -0.078060, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.262120, 0.006320, -0.061130, 0.125490, -0.385820, 0.065610, 
0.126870, 0.108110, -0.460870, 0.191030, -0.065570, 0.168700, 
-0.260500, 0.068830, 0.026900, 0.053250, -0.216080, -0.315940, 
-0.250780, 0.271980, -0.287940, -0.278150, -0.276190, -0.055690, 
-0.146610, -0.275330, -0.266790, -0.070180, -0.262830, -0.406640, 
-0.366580, -0.176620, -0.314830, -0.566990, -0.139320, 0.002510, 
-0.317100, -0.361700, -0.237270, -0.158500, -0.179850, -0.313550, 
0.019150, -0.241000, -0.209490, -0.044640, 0.065320, -0.261760, 
-0.092320, 0.041460, -0.220370, 0.023540, 0.079590, -0.148430, 
-0.367130, -0.070650, 0.179550, -0.101260, -0.252400, 0.068110, 
0.358340, -0.180990, -0.406890, 0.214470, 0.254700, -0.445220, 
-0.243670, 0.140330, 0.254480, -0.200920, -0.212550, 0.480600, 
0.107980, -0.517740, 0.067980, 0.560960, 0.083560, -0.220540, 
-0.084260, 0.593680, 0.339370, -0.147700, -0.074490, 0.551510, 
0.378830, -0.139160, 0.129180, 0.474050, 0.292790, -0.227910, 
-0.246850, 0.328990, 0.610170, -0.004020, -0.292140, 0.148580, 
0.218400, -0.059880, -0.319980, -0.062160, 0.031130, 0.091200, 
-0.401400, -0.129110, -0.163150, 0.169430, -0.219710, -0.271750, 
-0.119000, 0.209550, -0.313690, -0.214350, -0.046930, 0.411030, 
-0.396750, -0.021930, -0.159250, 0.332010, -0.396520, -0.216980, 
-0.337040, 0.396260, -0.135260, -0.115850, -0.386470, 0.474190, 
-0.103910, -0.099850, -0.206130, 0.431130, -0.278050, 0.152850, 
-0.396870, 0.240970, -0.193600, 0.100710, -0.308080, 0.075960, 
-0.023760, -0.061310, -0.087280, 0.379730, 0.074140, -0.242760}
    },
    { /* unit 151 (Old: 151) */
      0.0, -0.467590, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.020630, -0.058790, 0.098150, -0.088000, 0.020880, -0.353200, 
0.107150, -0.162380, 0.061820, -0.045530, -0.240610, -0.150410, 
-0.196440, -0.051780, -0.136640, 0.024000, -0.015750, -0.301050, 
0.112360, -0.098050, -0.166360, 0.043940, 0.027160, -0.070980, 
-0.049430, -0.145260, -0.075120, 0.000600, 0.225100, -0.040070, 
-0.040970, 0.169980, 0.230090, 0.132780, 0.287280, -0.105280, 
-0.108480, 0.019760, 0.212800, 0.187240, 0.130340, -0.157550, 
0.276370, 0.236050, 0.047070, 0.025370, 0.194750, 0.043690, 
-0.182700, -0.143390, 0.260440, -0.193200, -0.174030, 0.074740, 
0.099030, 0.030950, 0.057430, -0.140890, 0.113420, 0.098520, 
-0.064910, 0.194200, 0.214870, 0.071250, -0.111560, 0.201530, 
0.111760, 0.136800, 0.231100, -0.046430, -0.204170, -0.001790, 
0.145280, 0.139950, 0.035570, -0.014850, -0.031740, -0.028450, 
-0.159170, 0.088870, -0.132810, 0.344880, -0.044090, 0.046440, 
-0.009500, 0.264500, 0.201520, -0.042480, -0.281510, 0.092450, 
0.101890, 0.137350, -0.281740, 0.146990, 0.048990, -0.183020, 
-0.110710, 0.084410, 0.023860, -0.019440, 0.013570, 0.206400, 
-0.054130, -0.227320, -0.168740, 0.035040, 0.123850, 0.045270, 
-0.146100, -0.147380, -0.047190, 0.029240, -0.072970, 0.039260, 
-0.035100, -0.153140, -0.237860, -0.079360, -0.001790, 0.049090, 
-0.032740, -0.072870, -0.018670, -0.190580, 0.001870, 0.201890, 
-0.116090, 0.012250, -0.183970, 0.114840, -0.077060, -0.252090, 
-0.164310, -0.133420, -0.303370, -0.229440, 0.180120, 0.031080, 
-0.331850, 0.159300, 0.246510, -0.086160, -0.147430, 0.164090}
    },
    { /* unit 152 (Old: 152) */
      0.0, -0.301280, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.147350, -0.258710, 0.045610, -0.165110, -0.204550, 
-0.196180, 0.052220, -0.101480, 0.030890, -0.026530, 0.179260, 
0.000510, -0.130310, 0.086450, 0.142150, 0.061780, -0.206950, 
0.048780, -0.145910, -0.051440, -0.131560, -0.125250, -0.010950, 
-0.127510, -0.250700, -0.137610, -0.150030, 0.048940, -0.084840, 
-0.190960, -0.134120, -0.144950, -0.076710, -0.008490, -0.204560, 
0.043120, -0.107160, -0.161000, -0.160940, -0.109570, -0.087870, 
-0.125820, -0.187350, 0.158410, 0.058730, -0.255520, -0.074500, 
0.024020, -0.167450, -0.239590, 0.160850, -0.020000, -0.158940, 
-0.071960, -0.071340, 0.123060, 0.208600, 0.115560, 0.081960, 
-0.126280, 0.088120, -0.114450, -0.120730, -0.180030, 0.031420, 
0.205660, -0.121020, 0.054760, -0.271810, -0.086540, 0.131610, 
0.037720, -0.172000, -0.059740, 0.065840, 0.183630, 0.007900, 
0.231700, 0.023930, 0.041740, 0.075040, 0.147440, 0.196880, 
-0.122290, 0.120270, 0.186500, -0.006680, -0.147740, -0.221260, 
0.149120, 0.039750, 0.038970, -0.094620, 0.006830, -0.077890, 
-0.036520, -0.224140, -0.030540, 0.189450, 0.120330, -0.062770, 
-0.126080, -0.061490, -0.203960, 0.030900, -0.016610, -0.192410, 
-0.066140, -0.238640, 0.103850, -0.006420, -0.028440, 0.044800, 
-0.107330, -0.086350, 0.114290, 0.116100, -0.186350, -0.030890, 
-0.014860, 0.007690, -0.082490, -0.180830, 0.000760, -0.213920, 
-0.289360, -0.026700, -0.082320, 0.131710, 0.026170, 0.066310, 
-0.151110, 0.047970, 0.051570, -0.217150, -0.268860, 0.078470, 
-0.143500, -0.158250, 0.051870, -0.118010, -0.228430, 0.106130, 
0.067180}
    },
    { /* unit 153 (Old: 153) */
      0.0, -0.087520, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.173040, -0.045520, 0.223400, -0.126060, 0.027660, -0.270140, 
0.196750, 0.204130, 0.023610, -0.218960, 0.343240, 0.208610, 
0.011700, 0.165970, 0.253170, 0.072290, -0.085780, -0.051040, 
0.112260, -0.093620, 0.114410, -0.070190, 0.075290, -0.182200, 
-0.037360, 0.071980, 0.153990, -0.118670, 0.172660, -0.067200, 
0.154600, -0.249370, 0.107840, -0.226150, 0.109050, 0.063430, 
0.012780, -0.320140, 0.057280, 0.107750, -0.056690, 0.011800, 
0.161240, 0.104070, -0.000020, -0.170890, -0.041740, -0.151570, 
-0.026290, -0.237930, 0.062990, -0.176640, -0.066860, -0.195560, 
0.168110, -0.101530, 0.134320, -0.018260, -0.026630, -0.257510, 
-0.158380, -0.000630, -0.326100, -0.235460, -0.237830, -0.039520, 
0.029080, -0.074620, 0.109180, 0.091430, -0.019690, -0.370310, 
0.000340, 0.244710, -0.199250, -0.132060, -0.161480, 0.020140, 
0.024300, 0.023780, 0.182600, 0.266040, -0.091300, -0.007230, 
-0.067990, 0.094840, -0.004120, -0.269830, 0.118440, 0.238540, 
0.238280, -0.066730, -0.010800, 0.196070, 0.360380, -0.124330, 
0.048890, 0.212660, 0.017300, 0.024900, 0.113910, 0.066030, 0.288970, 
-0.219070, -0.088050, 0.054870, -0.047840, -0.194580, 0.097260, 
0.007070, 0.053200, -0.106470, -0.111300, -0.106230, 0.023680, 
-0.070490, -0.141130, -0.404920, -0.302640, 0.124910, -0.315490, 
-0.280460, -0.306930, 0.260440, -0.210210, -0.119070, -0.479130, 
0.007980, -0.377190, -0.338970, -0.142780, 0.096490, -0.272660, 
-0.203790, -0.118740, -0.094670, -0.456970, -0.257170, -0.012050, 
-0.245880, -0.160940, -0.200010, -0.078350, -0.154090}
    },
    { /* unit 154 (Old: 154) */
      0.0, -1.195360, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.237340, 0.085900, -0.153180, -0.109570, 0.186030, -0.054490, 
-0.638830, -0.112320, -0.318320, -0.135010, -0.969910, -0.221750, 
-0.226270, -0.194850, -0.835490, -0.078830, -0.437920, -0.138960, 
-0.939010, -0.215540, -0.306000, 0.095700, -0.495200, 0.089540, 
-0.131590, -0.119960, -0.425160, 0.425510, 0.236370, 0.272850, 
-0.306800, 0.116410, 0.431420, 0.378530, -0.069660, -0.104890, 
0.037240, 0.744450, 0.132400, 0.154820, -0.173630, 0.667760, 
-0.184320, 0.104990, -0.028510, 0.562720, -0.138110, -0.028190, 
-0.099900, 0.324580, -0.079760, 0.301780, -0.279460, 0.254860, 
0.508710, 0.609040, -0.251030, -0.222780, 0.680850, 0.879630, 
-0.246120, -0.559640, 1.116620, 1.036150, -0.035330, -0.773360, 
0.821450, 0.925620, -0.177120, -0.947870, 1.161470, 0.777490, 
-0.286900, -1.270870, 1.038620, 0.824810, -0.713750, -1.232760, 
1.005870, 0.435250, -0.728800, -1.212540, 0.760550, 0.263260, 
-0.182980, -1.080030, 0.387040, 0.146130, 0.081690, -0.727250, 
-0.502660, 0.380540, 0.112580, -0.228190, -0.718700, 0.225700, 
0.445820, 0.004410, -0.392700, 0.314710, 0.277200, -0.363420, 
-0.506860, 0.163310, 0.111920, -0.229420, 0.076660, 0.354890, 
0.047510, 0.271620, 0.440330, 0.379130, -0.397550, 0.284320, 
0.668620, 0.221930, -0.122660, 0.409970, 0.697430, 0.207700, 
0.109840, 0.316260, 0.888020, 0.201000, 0.250330, 0.253560, 0.781640, 
0.219120, 0.386340, 0.071460, 0.612560, -0.090270, 0.709800, 
-0.067990, 0.018060, 0.221360, 0.849750, -0.121710, -0.186910, 
0.076400, 1.035870, 0.037380, -0.864600, -0.085250}
    },
    { /* unit 155 (Old: 155) */
      0.0, 0.024400, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.272760, 0.147300, -0.310530, -0.364510, -0.186720, 
0.361090, -0.283560, -0.400850, -0.108540, 0.063280, -0.280680, 
-0.253670, -0.323730, 0.178470, -0.339570, -0.180970, -0.197760, 
-0.107290, -0.534400, 0.082350, 0.081800, 0.016730, -0.344790, 
0.066990, -0.152680, -0.259880, -0.375240, -0.105890, -0.060010, 
-0.173940, -0.135410, -0.244830, -0.168940, -0.260720, -0.290340, 
-0.472430, 0.026520, -0.327440, -0.296200, -0.451500, -0.236940, 
-0.508450, 0.138690, -0.422770, -0.091340, -0.489560, 0.105140, 
-0.219140, -0.083010, -0.178920, -0.361990, -0.019480, 0.075130, 
-0.276800, -0.396730, 0.068010, 0.016840, -0.553740, -0.464930, 
0.058680, 0.209030, -0.246330, -0.332240, 0.259830, 0.016700, 
-0.356370, -0.284430, 0.401690, 0.028660, -0.595140, -0.195170, 
0.575990, 0.136640, -0.358750, -0.175500, 0.614450, 0.050630, 
-0.496070, -0.226680, 0.502830, -0.035120, -0.365580, -0.288770, 
0.576160, 0.308850, -0.022920, -0.416270, 0.559080, 0.271470, 
-0.321170, -0.288280, 0.113090, 0.199170, -0.124840, -0.441860, 
0.251350, 0.035210, 0.087250, -0.113950, 0.068100, 0.121220, 
-0.174460, -0.264810, -0.070010, 0.050090, -0.027520, -0.169530, 
-0.120730, -0.152960, 0.068450, 0.279600, -0.351460, -0.113120, 
0.103400, 0.284820, -0.069750, -0.166160, 0.298220, -0.078530, 
-0.327780, -0.001530, 0.296110, 0.125030, -0.227880, 0.157870, 
0.304080, -0.145900, -0.061900, 0.031590, 0.346270, 0.079950, 
-0.193480, -0.083340, 0.313390, -0.091630, -0.411740, 0.044960, 
0.116860, 0.228520, -0.646270, -0.065490, 0.002570, 0.071970, 
-0.319710}
    },
    { /* unit 156 (Old: 156) */
      0.0, 0.345120, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.270080, -0.877210, 0.150610, 0.673680, -0.312210, 
-0.845540, 0.610520, 0.497710, -0.032300, -0.890660, 0.462410, 
0.392600, -0.129840, -0.778020, 0.621640, 0.297650, -0.239980, 
-0.733250, 0.542260, 0.186570, -0.362670, -0.695010, 0.627860, 
-0.098360, -0.181730, -0.538860, 0.662960, -0.305250, -0.307680, 
-0.664440, 0.296070, -0.353230, -0.448410, -0.491410, 0.342570, 
-0.157180, -0.584970, -0.500130, 0.293390, -0.130240, -0.392410, 
-0.324730, 0.385770, -0.417920, -0.231080, -0.373420, 0.083560, 
-0.024590, -0.320110, -0.326050, 0.124270, -0.117170, -0.051550, 
-0.486250, 0.060430, 0.252890, -0.279850, -0.125840, -0.339570, 
0.106730, -0.029070, -0.292960, -0.128660, 0.320560, -0.006240, 
-0.358770, -0.263260, 0.476710, -0.437910, -0.150030, -0.445910, 
0.241360, -0.265100, 0.172800, -0.395390, 0.190790, -0.223410, 
0.197900, -0.409870, 0.328670, -0.064280, -0.006920, -0.136330, 
0.478090, 0.042360, -0.057340, 0.033440, 0.491280, 0.152500, 
-0.158050, 0.035440, 0.371350, 0.184780, -0.370800, -0.009600, 
0.387850, 0.094630, -0.356020, 0.198680, 0.304410, 0.231340, 
-0.562140, -0.029510, 0.305040, 0.036370, -0.512600, 0.105520, 
0.261730, 0.143810, -0.548260, -0.065220, 0.383450, -0.188500, 
-0.324400, 0.035370, 0.026590, -0.149070, -0.502700, -0.136460, 
0.291450, -0.379670, -0.108520, -0.131150, 0.280250, -0.833740, 
-0.252150, -0.000620, 0.111280, -0.970600, -0.122780, -0.045230, 
0.319730, -1.018430, -0.219690, -0.054460, 0.386420, -0.882620, 
-0.280600, 0.242480, 0.348580, -1.041760, -0.344670, 0.290160, 
0.048080}
    },
    { /* unit 157 (Old: 157) */
      0.0, -0.406550, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.122480, -0.163110, 0.230690, 0.255550, 0.160240, -0.117320, 
-0.148820, -0.016290, -0.161340, -0.145420, -0.020030, -0.132220, 
-0.401200, -0.456080, -0.020400, 0.084190, -0.360740, -0.439110, 
0.249320, 0.022590, -0.291020, -0.544300, 0.276250, 0.164380, 
-0.291670, -0.215560, 0.193930, 0.005050, -0.349710, -0.339580, 
0.066620, 0.086890, -0.096430, -0.035480, 0.068310, 0.207500, 
-0.125560, -0.208820, -0.020250, 0.164240, 0.019980, 0.095970, 
-0.011780, 0.019250, -0.240220, 0.215530, -0.315940, 0.200490, 
-0.283920, -0.089840, -0.105490, 0.106820, -0.226310, 0.085220, 
0.208030, 0.042790, -0.239090, -0.352780, 0.503010, 0.018310, 
-0.106060, -0.148410, 0.373930, -0.090310, -0.258990, -0.339810, 
0.376910, -0.309680, -0.265500, -0.327680, 0.559140, -0.376240, 
-0.204890, -0.643450, 0.454870, -0.575290, -0.127150, -0.511720, 
0.079860, -0.528650, 0.047170, -0.372440, 0.231230, -0.586110, 
0.144530, -0.266930, 0.180920, -0.392380, 0.076870, -0.427800, 
-0.088410, -0.398540, -0.085030, -0.214960, -0.296220, -0.385050, 
0.060850, -0.093680, 0.018590, -0.309710, 0.283490, -0.002630, 
-0.044770, -0.169370, 0.104370, -0.080820, -0.049570, 0.092030, 
-0.074310, -0.066170, -0.111640, -0.165750, 0.103990, -0.109230, 
-0.037540, -0.077160, -0.016020, -0.173190, 0.113750, -0.278470, 
-0.009840, 0.183500, 0.259320, -0.148140, -0.146250, 0.051240, 
0.359010, -0.420690, 0.115690, 0.025300, 0.203810, -0.306130, 
0.121000, 0.194500, 0.184480, 0.023250, 0.324830, 0.077390, 
-0.176020, -0.051430, 0.057860, 0.020560, 0.085700, -0.152440}
    },
    { /* unit 158 (Old: 158) */
      0.0, 0.073050, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.110520, 0.169450, 0.129950, -0.058210, 0.409180, -0.090590, 
-0.043270, 0.187660, 0.023180, -0.375330, 0.175130, -0.079150, 
0.114700, -0.223970, -0.272050, -0.134020, 0.314510, -0.248560, 
-0.388760, 0.071100, 0.380700, -0.136770, -0.145190, -0.087120, 
0.513590, -0.294560, -0.310290, 0.052030, 0.518320, -0.092920, 
-0.050070, 0.145380, 0.320330, -0.067710, -0.205200, 0.306480, 
0.523700, 0.004990, -0.083340, -0.068050, 0.621590, -0.290090, 
0.108550, 0.283210, 0.420070, -0.115830, 0.203220, 0.085910, 
0.455520, -0.188460, -0.066630, -0.010650, 0.210950, 0.127480, 
-0.194440, 0.169450, -0.106310, 0.182990, -0.307980, 0.280110, 
-0.267980, -0.195140, -0.563230, 0.157000, -0.459310, 0.035780, 
-0.657090, 0.070850, -0.721480, -0.141240, -0.632770, -0.188770, 
-0.882710, -0.474010, -0.626240, 0.028540, -1.038040, -0.585410, 
-0.489330, -0.407900, -1.104770, -0.304240, -0.455350, -0.304300, 
-0.792950, -0.294470, -0.134960, -0.004660, -0.855980, -0.555880, 
-0.149990, -0.079150, -0.383490, -0.233780, -0.089490, -0.208040, 
-0.261070, -0.309250, 0.142440, 0.013190, -0.331900, -0.445110, 
0.005450, -0.028640, -0.152590, -0.372300, 0.158340, -0.028410, 
-0.016910, -0.150240, 0.254350, -0.016510, -0.312210, -0.025540, 
0.416720, 0.075660, -0.100960, 0.090940, 0.307710, 0.087210, 
-0.247390, 0.381550, 0.452060, 0.124140, -0.543100, 0.126330, 
0.185260, 0.125730, -0.286260, 0.446330, 0.218100, 0.107490, 
-0.317000, 0.108930, 0.222490, 0.353350, -0.373760, 0.113660, 
0.215680, 0.182450, -0.165560, 0.157280, 0.186260, 0.151870}
    },
    { /* unit 159 (Old: 159) */
      0.0, -0.330180, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.110870, 0.826800, -0.367910, -0.693730, -0.395610, 
0.227380, -0.584950, -0.343070, -0.383940, 0.168930, -0.589490, 
-0.457590, -0.659440, -0.167370, -0.258770, -0.466340, -0.496740, 
-0.287970, -0.160890, -0.355320, -0.692310, -0.150780, 0.216520, 
-0.089340, -0.637460, 0.050330, 0.470690, -0.037880, -0.694230, 
-0.119350, 0.488100, -0.197140, -0.884640, 0.206320, 0.756120, 
-0.391180, -0.743290, 0.230530, 0.685650, -0.330080, -0.698910, 
0.210350, 0.364780, -0.523310, -0.939390, 0.162990, 0.222360, 
-0.589800, -1.001020, 0.104720, 0.303730, -0.626590, -0.937870, 
-0.191890, -0.040670, -0.388970, -0.668410, -0.262090, 0.311840, 
-0.230350, -0.710320, -0.198090, 0.354230, 0.131600, -0.498780, 
-0.472210, 0.241710, 0.366040, -0.564190, -0.288690, 0.379950, 
0.566990, -0.334300, -0.484990, 0.265650, 0.558640, -0.458240, 
-0.469480, 0.393210, 0.677290, -0.442500, -0.185000, 0.039560, 
0.485450, -0.504610, -0.216080, 0.149680, 0.378620, -0.312420, 
-0.053260, -0.064210, 0.625830, -0.350030, 0.288800, 0.184530, 
0.627100, -0.332480, 0.153610, 0.103940, 0.357310, -0.256070, 
0.183130, -0.088410, 0.092750, -0.013560, 0.012140, 0.021340, 
0.070760, -0.042940, 0.370310, -0.014830, -0.384550, -0.286800, 
0.250720, 0.289700, -0.549290, -0.024840, 0.533800, 0.183350, 
-0.352940, 0.339980, 0.101150, 0.285270, -0.723930, 0.466960, 
0.331830, 0.419350, -0.480260, 0.283810, -0.068900, 0.450090, 
-0.740480, 0.304460, 0.028280, 0.439370, -0.365200, 0.343200, 
-0.200660, 0.242360, -0.337830, 0.371660, -0.076860, -0.375190, 
0.015470}
    },
    { /* unit 160 (Old: 160) */
      0.0, -0.322370, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.140640, 0.406300, -0.127750, -0.895670, -0.198390, 
-0.112570, -0.261000, -0.756920, -0.072970, -0.185720, -0.380060, 
-0.328140, -0.415250, -0.133750, -0.217410, -0.295380, -0.221630, 
-0.036610, -0.163530, -0.084920, -0.543910, -0.002030, -0.323110, 
0.110430, -0.417410, -0.180750, -0.222810, -0.086430, -0.602580, 
-0.140050, 0.147750, -0.021460, -0.509630, -0.053250, 0.115300, 
-0.304100, -0.479090, -0.226250, 0.399370, -0.515470, -0.604130, 
-0.220910, 0.609230, -0.435180, -0.289940, -0.323780, 0.543190, 
-0.203560, -0.279590, -0.275310, 0.224400, -0.410900, -0.106830, 
-0.339820, 0.139710, -0.043650, -0.303100, -0.193050, -0.069760, 
0.256850, -0.333310, -0.667110, -0.107660, 0.770070, -0.077500, 
-0.471870, -0.436570, 0.795430, -0.249420, -0.561450, -0.435760, 
1.114870, -0.480760, -0.763340, -0.335860, 0.955230, -0.237990, 
-0.762920, -0.307110, 1.134000, -0.403550, -0.586440, -0.331920, 
0.937880, -0.256240, -0.420640, -0.311460, 0.578740, -0.232590, 
-0.345540, -0.312030, 0.502640, -0.203460, 0.194460, -0.269100, 
0.406540, -0.303960, 0.134930, -0.145740, 0.155950, -0.315260, 
-0.335110, -0.018600, 0.261540, 0.015300, -0.515860, -0.010390, 
-0.044450, -0.000920, 0.082300, -0.204750, -0.304640, 0.059930, 
0.282010, -0.187510, -0.485490, 0.224020, 0.099390, -0.287880, 
-0.505590, 0.660170, 0.281900, -0.475070, -0.496870, 0.568060, 
0.413030, -0.307760, -0.772010, 0.288310, 0.310330, -0.279350, 
-0.766590, 0.500230, 0.184130, -0.254940, -0.711050, 0.375750, 
0.415750, -0.324640, -0.633420, 0.254730, 0.007510, -0.095630, 
-0.151720}
    },
    { /* unit 161 (Old: 161) */
      0.0, -0.245840, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.188980, -0.341490, 0.088200, -0.162260, 0.163950, 
-0.144190, 0.109270, -0.133870, 0.202580, -0.210730, 0.035700, 
-0.113860, -0.059490, -0.228370, -0.116180, -0.065530, 0.026660, 
-0.144810, -0.006040, -0.167820, 0.149240, -0.130060, -0.234150, 
0.031960, 0.281520, -0.399650, -0.231170, 0.083990, 0.154400, 
-0.053610, -0.296320, -0.004300, 0.350120, -0.075180, 0.040270, 
-0.109540, 0.206090, -0.409620, -0.257770, 0.215030, 0.243350, 
-0.269880, 0.003450, 0.035340, 0.091400, -0.005100, 0.012800, 
0.158160, 0.212560, 0.013400, -0.068350, 0.043800, 0.272550, 
0.020050, -0.011680, 0.023630, -0.096880, -0.051000, 0.023710, 
0.214280, 0.128170, 0.045220, -0.147600, 0.123230, -0.252430, 
0.107210, -0.183820, -0.044850, -0.056900, -0.064550, 0.021180, 
-0.024210, -0.427440, -0.183670, -0.219070, -0.256280, -0.399360, 
-0.350640, -0.034200, -0.230120, -0.524800, -0.422060, -0.043620, 
-0.081950, -0.165950, -0.112720, -0.209870, -0.014280, -0.469300, 
-0.301090, 0.143230, 0.146330, -0.260670, -0.162790, -0.104960, 
-0.073730, -0.072310, 0.033280, 0.159940, 0.186520, -0.341920, 
-0.035090, -0.139660, 0.029470, -0.223380, -0.007300, 0.111610, 
-0.007840, 0.103940, 0.164860, -0.252700, 0.249550, 0.161810, 
0.036910, -0.112930, 0.050630, -0.069050, -0.067990, -0.183210, 
0.169660, -0.107830, -0.038890, -0.116820, 0.131170, -0.012330, 
-0.070770, -0.214900, 0.062150, 0.211700, -0.088550, -0.094090, 
0.046810, 0.132700, 0.020890, -0.234980, 0.248450, 0.071950, 
-0.137790, 0.046950, 0.085780, 0.097970, 0.150940, -0.197730, 
0.076110}
    },
    { /* unit 162 (Old: 162) */
      0.0, 0.035220, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.278970, 0.825480, -0.212810, -0.071810, -0.303230, 
0.590640, -0.134630, -0.310580, -0.235150, 0.742210, -0.052180, 
0.084200, -0.543640, 0.553490, -0.182760, 0.001540, -0.539210, 
0.585620, -0.076800, -0.155180, -0.498430, 0.618250, -0.219120, 
-0.127400, -0.640330, 0.511260, -0.085810, -0.008910, -0.551790, 
0.333000, -0.249820, -0.082000, -0.181380, 0.644090, -0.161270, 
0.209710, -0.253170, 0.202380, -0.233650, -0.014970, 0.073390, 
0.181730, -0.107990, 0.195450, 0.134410, 0.006100, -0.168300, 
-0.038590, 0.156200, 0.117720, 0.088360, 0.069790, 0.266620, 
0.207060, -0.337010, -0.293000, -0.005170, 0.254830, -0.397600, 
-0.263620, 0.085750, 0.256190, -0.416530, -0.188760, -0.094730, 
0.230630, -0.637260, -0.174550, -0.092350, 0.193590, -0.263380, 
0.160380, -0.224590, 0.153600, -0.454090, 0.067530, -0.107100, 
-0.079450, -0.493560, 0.259640, -0.327780, 0.157320, -0.351970, 
0.302860, -0.384100, -0.028890, -0.193700, 0.197830, -0.206370, 
0.143140, 0.013340, -0.161840, -0.046950, 0.261340, 0.159800, 
0.140560, -0.132120, 0.174180, 0.176110, -0.114030, -0.157570, 
-0.024750, -0.066810, 0.041240, -0.257340, 0.010180, -0.045990, 
0.076590, -0.174820, -0.151220, -0.132110, -0.017710, -0.123030, 
0.031730, -0.401030, 0.011590, 0.174390, -0.204740, -0.471030, 
0.388640, 0.107390, -0.122710, -0.674320, 0.305790, -0.212020, 
-0.069970, -0.496130, 0.147040, -0.260580, -0.256890, -0.448750, 
0.205090, -0.185600, -0.110460, -0.364520, 0.014640, -0.048100, 
-0.138000, -0.021950, -0.073690, 0.010270, -0.173000, 0.104300, 
0.008090}
    },
    { /* unit 163 (Old: 163) */
      0.0, -0.457520, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.362510, -0.018990, -0.053850, -0.573700, -0.287030, 
0.119220, -0.044480, -0.399470, -0.238880, -0.079380, -0.165470, 
-0.293590, -0.277040, 0.089800, -0.353650, -0.266450, -0.416560, 
0.146640, -0.371870, 0.071580, -0.413870, -0.147990, -0.161350, 
0.062820, -0.313750, -0.151320, -0.241820, -0.097860, -0.075490, 
0.105990, 0.156850, -0.100930, -0.258590, 0.165880, 0.213700, 
-0.121580, -0.284860, -0.015360, 0.242200, -0.216110, -0.367750, 
0.137670, 0.284020, -0.242360, -0.351380, 0.023370, 0.126240, 
-0.280130, -0.229060, -0.167290, -0.012640, 0.000260, -0.169660, 
0.052300, 0.052420, 0.156600, -0.166310, 0.076150, -0.256300, 
0.273210, -0.224880, 0.004910, -0.190180, 0.071000, -0.305740, 
-0.178090, 0.022740, 0.262640, -0.317080, 0.034410, -0.004090, 
0.509540, -0.208550, -0.149890, -0.082210, 0.455490, -0.089080, 
-0.381060, -0.035690, 0.263410, -0.217270, 0.015700, -0.223080, 
0.423620, -0.339830, 0.053490, -0.101000, 0.272240, -0.165890, 
-0.218770, -0.166160, 0.098970, -0.064160, 0.188070, 0.029680, 
0.193970, 0.098160, 0.053970, 0.097560, 0.106290, 0.023590, 
-0.170920, 0.300230, 0.050840, -0.208690, -0.235240, 0.281660, 
0.012600, -0.297260, -0.061790, 0.170050, -0.291030, -0.053640, 
0.048900, -0.055350, -0.252330, 0.077730, -0.003430, 0.177480, 
-0.349400, 0.124630, -0.023850, 0.082850, -0.206140, 0.059080, 
0.165560, 0.177110, -0.190630, -0.056350, 0.110980, 0.026170, 
-0.421950, -0.326010, -0.213680, 0.127780, -0.394890, -0.328290, 
0.119420, 0.168060, -0.232220, -0.034030, 0.105910, -0.112130, 
-0.366130}
    },
    { /* unit 164 (Old: 164) */
      0.0, -0.053220, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.141650, -0.074910, 0.232160, 0.136460, -0.133820, 
-0.022520, 0.188490, 0.179660, 0.060350, -0.260340, 0.026560, 
-0.045330, -0.088450, -0.306820, 0.067100, -0.054490, -0.374280, 
-0.096960, 0.126250, 0.021320, -0.236490, -0.054780, 0.021050, 
0.187670, -0.143760, -0.230510, 0.160110, -0.080370, -0.084120, 
-0.347310, 0.054770, -0.055180, -0.220060, -0.339010, 0.137330, 
0.065890, -0.091510, 0.032570, 0.070180, 0.038150, -0.137180, 
-0.185450, 0.254320, 0.025080, 0.106850, -0.187830, 0.232050, 
0.073350, -0.139580, 0.133530, 0.177170, 0.058450, 0.019610, 
-0.005170, 0.205390, -0.111920, -0.162200, -0.084430, 0.174510, 
0.020430, -0.190450, -0.002740, -0.145740, -0.018300, -0.047750, 
0.034270, 0.106150, -0.167360, -0.211670, -0.296120, -0.267070, 
0.049060, -0.477090, -0.017440, -0.136790, 0.019970, -0.355060, 
0.008620, -0.094380, -0.224780, -0.024050, -0.093860, -0.220000, 
-0.175890, -0.168560, 0.132000, -0.256640, -0.185000, -0.107530, 
-0.229890, 0.111110, -0.048730, 0.113270, 0.058840, -0.101970, 
-0.059800, -0.087620, -0.279640, -0.102910, -0.164000, 0.026720, 
-0.205120, 0.077350, -0.133890, 0.098830, 0.029040, -0.056850, 
0.245450, -0.130950, 0.046340, 0.117190, 0.009060, 0.062940, 
-0.172690, -0.127660, 0.177120, -0.147320, -0.010200, -0.135190, 
0.175570, 0.025240, -0.290880, 0.016400, 0.157560, 0.020530, 
-0.158690, -0.083800, 0.172360, -0.265290, 0.035930, -0.182290, 
0.078050, -0.179090, -0.181240, 0.036750, 0.178970, -0.264930, 
0.063650, -0.073970, 0.017210, -0.303480, -0.105320, 0.149620, 
0.039920}
    },
    { /* unit 165 (Old: 165) */
      0.0, -0.562940, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.229980, 0.000740, 0.072830, -0.182560, -0.050690, 0.091850, 
0.212670, 0.062060, 0.060200, -0.073360, 0.020620, -0.172030, 
-0.057840, -0.000250, 0.178810, 0.105450, 0.059240, 0.104850, 
0.113710, -0.140170, 0.256780, 0.395430, 0.226830, -0.011170, 
0.154620, 0.205630, 0.122400, -0.204260, 0.068810, 0.357380, 
0.072090, -0.059450, 0.127840, 0.130760, -0.082200, -0.248310, 
0.103230, 0.330940, 0.015430, 0.063540, 0.127250, 0.189610, 
-0.226020, 0.104660, -0.029310, 0.095570, -0.033410, 0.108740, 
-0.057120, -0.055600, 0.078770, -0.236090, -0.058410, -0.300330, 
-0.136460, -0.263900, 0.129740, -0.210610, -0.379560, -0.114030, 
0.284380, -0.282950, -0.160250, -0.012270, 0.412160, -0.251610, 
-0.162020, -0.234270, 0.348510, -0.296770, -0.278240, -0.057160, 
0.448050, -0.242550, -0.112680, 0.096110, 0.569070, -0.153790, 
-0.278180, -0.065760, 0.331580, 0.038700, -0.030100, -0.046700, 
0.234300, 0.010140, -0.244740, 0.059740, 0.182100, 0.025370, 
-0.034170, -0.015040, 0.091520, -0.119400, -0.321110, 0.062440, 
-0.096780, 0.144300, -0.049870, -0.064160, 0.114000, -0.205680, 
-0.020590, 0.086600, -0.104040, -0.031280, -0.346720, 0.056790, 
0.133410, -0.059790, -0.229010, 0.275430, 0.094340, -0.337450, 
-0.055750, -0.055710, 0.047660, -0.214520, 0.001970, 0.023380, 
-0.133760, -0.409380, -0.058520, 0.192370, -0.224950, -0.028790, 
-0.060440, 0.191440, -0.367510, -0.332020, -0.128480, -0.130440, 
-0.331650, 0.073540, -0.212370, -0.215280, -0.349410, 0.015430, 
0.076560, -0.100380, -0.337070, 0.023960, -0.065430, -0.099210}
    },
    { /* unit 166 (Old: 166) */
      0.0, -0.258900, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.088730, 0.340750, -0.091660, -0.083370, -0.194550, 0.135620, 
0.071530, -0.140950, -0.087570, 0.101250, -0.033930, -0.036410, 
-0.208620, -0.043280, -0.051240, -0.387910, -0.247360, -0.173220, 
-0.237260, -0.152020, -0.282150, -0.115880, 0.088000, 0.011200, 
-0.345330, -0.066770, 0.116720, 0.074230, -0.350930, 0.073120, 
0.033460, -0.050800, -0.344560, 0.136980, -0.078370, -0.257050, 
-0.344380, 0.160890, 0.234590, -0.087170, -0.338450, 0.156650, 
0.032080, -0.298670, -0.150340, -0.155240, 0.030730, 0.105290, 
-0.089440, -0.027870, -0.096100, -0.175080, -0.057550, 0.027550, 
0.002310, 0.142680, -0.404270, 0.156070, -0.221000, -0.098340, 
-0.370930, -0.106440, -0.103240, 0.322150, -0.345340, 0.021750, 
-0.251090, 0.107800, -0.228950, -0.097000, -0.185190, 0.276160, 
-0.309340, -0.088650, -0.053060, 0.185590, -0.519210, -0.350630, 
-0.006550, 0.121470, -0.154340, -0.075900, -0.150610, 0.329630, 
-0.305170, -0.088790, 0.110580, 0.045740, -0.086770, -0.061710, 
-0.090930, 0.211900, 0.020760, -0.249560, 0.136180, -0.055490, 
-0.043190, -0.112320, -0.131980, 0.031540, -0.050670, -0.123210, 
0.248990, 0.154060, -0.200770, -0.197100, 0.022560, -0.135500, 
-0.142980, -0.211990, -0.017060, -0.133680, -0.101190, 0.060460, 
0.192250, -0.150440, -0.003920, -0.103150, 0.333800, -0.258920, 
-0.174290, 0.012770, 0.166740, -0.421980, -0.221020, -0.148890, 
0.247310, -0.486900, 0.100290, 0.137760, 0.070110, -0.424050, 
-0.095290, -0.224680, 0.157660, -0.016570, -0.123440, -0.151590, 
0.163130, -0.304750, -0.150640, 0.028930, -0.000600, -0.297890}
    },
    { /* unit 167 (Old: 167) */
      0.0, -0.397840, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.512350, -0.032530, 0.296310, 0.208820, -0.588160, 0.176700, 
-0.010870, 0.117870, -0.508970, 0.223240, 0.026590, -0.059650, 
-0.170920, 0.074100, 0.271610, 0.088540, -0.118930, -0.048170, 
0.217920, 0.160620, 0.199320, 0.094260, -0.041640, -0.146520, 
0.069940, 0.016030, -0.004480, 0.036240, 0.079750, -0.091750, 
0.054700, -0.028080, 0.126730, -0.025470, -0.077730, -0.058820, 
-0.117090, 0.002750, -0.142040, 0.018460, -0.254810, 0.267180, 
-0.445930, -0.079690, -0.336560, 0.153060, -0.083500, -0.138010, 
-0.438230, 0.106810, -0.166150, 0.040180, -0.163170, -0.306200, 
0.226980, -0.062690, 0.016920, -0.583820, 0.336160, -0.392260, 
0.191880, -0.331160, 0.357100, -0.252140, 0.226800, -0.595660, 
0.208530, -0.547580, 0.520220, -0.602610, 0.302150, -0.235810, 
0.845990, -0.719670, 0.645790, -0.621520, 0.714180, -0.258130, 
0.287830, -0.327560, 0.688650, -0.308180, 0.206240, -0.405380, 
0.381890, -0.362270, -0.064460, -0.250070, 0.618960, -0.375430, 
-0.077870, -0.373800, 0.215000, -0.280010, -0.238870, -0.318000, 
0.340710, -0.110120, -0.266970, -0.125540, 0.362550, -0.083670, 
-0.591810, -0.004580, 0.215640, -0.103530, -0.457240, -0.052100, 
-0.043030, -0.192700, -0.331650, -0.159360, -0.223760, 0.056110, 
-0.216800, -0.157950, -0.115220, -0.332880, -0.053200, -0.017780, 
-0.287530, -0.108480, 0.129380, 0.151900, -0.058010, -0.100410, 
0.319900, 0.186590, -0.287180, -0.000250, 0.313370, 0.196240, 
0.085820, -0.113610, 0.277720, 0.342390, -0.172730, 0.293060, 
0.134670, 0.216050, -0.139040, 0.222100, 0.001470, 0.106150}
    },
    { /* unit 168 (Old: 168) */
      0.0, 0.090310, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.083710, 0.265930, -0.034830, -0.188210, 0.306430, 0.216680, 
-0.226990, -0.212870, -0.016620, 0.007360, -0.026420, 0.069570, 
-0.048010, 0.205080, -0.074650, -0.154980, -0.170690, -0.104500, 
0.078110, 0.025470, -0.272180, 0.231440, 0.564010, 0.157010, 
-0.498310, -0.165210, 0.326180, 0.518090, -0.488530, -0.215490, 
0.338260, 0.560750, -0.283410, 0.096620, -0.002830, 0.561850, 
-0.181390, 0.128660, -0.120030, 0.362450, -0.052660, 0.130300, 
-0.354520, 0.240350, -0.111570, 0.309970, -0.283500, 0.381900, 
0.077760, 0.225740, -0.173730, 0.211280, -0.076720, -0.030750, 
-0.001570, 0.012220, -0.175370, -0.166680, -0.147510, 0.119480, 
-0.200600, -0.216560, 0.227400, -0.404310, -0.598780, -0.223530, 
0.175850, -0.812410, -0.572100, -0.373160, -0.052680, -1.107090, 
-0.433290, -0.294860, -0.130260, -1.321640, -0.459910, -0.094700, 
0.032790, -1.360720, -0.643630, -0.225000, 0.026100, -1.326450, 
-0.528300, -0.296020, -0.313930, -1.246830, -0.094690, -0.284220, 
-0.153610, -0.824660, -0.212380, 0.062770, -0.255920, -0.521610, 
-0.062860, 0.244350, -0.317560, -0.570510, -0.005530, 0.388580, 
-0.545680, -0.336770, -0.184850, 0.442430, -0.247910, -0.555120, 
0.179100, 0.413740, -0.161680, -0.704450, 0.187650, 0.242450, 
-0.059910, -0.704920, 0.118750, 0.277070, 0.099190, -1.089700, 
-0.039850, 0.293340, 0.308820, -0.835700, -0.148300, 0.049890, 
0.194310, -0.902380, 0.197720, 0.277170, 0.137810, -0.973340, 
0.123990, -0.175940, -0.125490, -0.672570, 0.171530, -0.069740, 
0.022270, -0.332870, 0.103180, -0.277230, 0.152720, -0.332520}
    },
    { /* unit 169 (Old: 169) */
      0.0, -0.007000, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.838940, 0.284760, -0.222150, 0.667610, -0.566680, 0.371490, 
-0.226990, -0.092780, 0.001500, 0.373800, -0.479750, -0.392640, 
0.019970, 0.006450, -0.137070, -0.694580, 0.237970, -0.225690, 
-0.142990, -0.611880, 0.502230, -0.171320, 0.112220, -0.825840, 
0.560160, -0.432910, 0.006070, -0.703600, 0.309970, -0.341750, 
-0.181610, -0.454010, 0.122540, -0.193570, -0.376630, -0.837560, 
-0.271550, 0.481910, -0.652240, -0.711410, -0.407530, 0.484000, 
-0.440240, -0.525520, -0.337150, 0.844000, -0.566220, -0.699680, 
-0.578390, 0.642960, -0.067410, -0.447430, -0.472430, 0.047140, 
0.389890, -0.654740, 0.143650, 0.055470, 0.499920, -0.495280, 
0.132580, -0.164750, 1.236960, -0.649020, 0.803760, -0.280150, 
1.531900, -0.828940, 0.734880, -0.584420, 1.524940, -0.726520, 
1.136650, -0.482850, 1.546520, -0.823310, 1.167840, -0.228710, 
1.115120, -0.558770, 1.054140, -0.600370, 0.775920, -0.286990, 
0.924880, -0.604170, -0.037210, -0.264000, 1.113320, -0.536590, 
-0.157980, -0.309000, 0.654480, -0.430440, -0.660680, -0.157660, 
0.574890, -0.745990, -0.740250, -0.319430, 0.453070, -0.332610, 
-0.991560, -0.259920, 0.352910, 0.038760, -0.884350, -0.514940, 
0.053520, -0.181200, -0.495780, -0.324530, -0.142250, -0.572370, 
-0.188260, -0.433820, -0.033500, -0.445860, -0.174530, -0.509970, 
-0.140710, -0.154710, 0.193560, -0.262470, 0.151700, -0.181550, 
0.351850, 0.030070, 0.247440, -0.105230, 0.132180, 0.259440, 
0.428210, 0.057150, 0.033920, 0.451550, 0.142760, -0.235350, 
-0.071530, 0.490100, -0.278610, 0.045860, -0.128420, 0.007170}
    },
    { /* unit 170 (Old: 170) */
      0.0, -0.047660, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.507830, -0.765730, -0.574320, 0.393660, -0.103540, 
-0.783610, -0.120160, 0.219250, 0.043140, -0.533920, -0.054450, 
0.317320, 0.034360, -0.379920, -0.111250, 0.359890, 0.189190, 
-0.603640, 0.090860, 0.401160, 0.239690, -0.629780, 0.078490, 
-0.057340, 0.329640, -0.828390, -0.092100, 0.053900, -0.133170, 
-0.694370, -0.510740, -0.174860, -0.220470, -0.897890, -0.258820, 
0.064030, -0.402680, -0.660230, -0.332650, -0.019320, -0.387580, 
-0.859590, -0.540670, -0.412110, -0.566110, -0.727470, -0.400660, 
-0.290590, -0.591430, -0.612940, -0.318600, -0.344710, -0.400910, 
-0.396600, -0.596210, -0.326320, -0.172030, -0.522950, -0.400680, 
-0.353670, -0.184170, -0.538230, -0.620470, -0.214270, 0.189220, 
-0.249680, -0.383220, -0.404830, 0.168340, -0.569570, -0.222910, 
0.044580, 0.290510, -0.474190, -0.259420, -0.068710, 0.117070, 
-0.124860, 0.064980, 0.012710, 0.069460, -0.018790, -0.001540, 
-0.065660, 0.128010, -0.259540, 0.036180, 0.028420, -0.133830, 
0.085710, -0.047890, -0.267720, -0.337230, -0.011270, -0.074560, 
-0.060720, -0.315020, 0.016770, -0.097500, -0.306110, 0.059050, 
-0.011100, 0.192130, -0.399290, -0.141970, 0.184320, 0.028480, 
-0.079820, -0.063920, -0.117690, 0.286440, -0.068960, 0.065070, 
-0.202050, 0.285350, 0.150180, -0.174690, -0.202240, 0.335720, 
0.408910, -0.193600, -0.270760, 0.151990, 0.259410, -0.314220, 
0.047490, 0.323080, 0.226950, -0.479940, 0.158630, 0.264440, 
0.157670, -0.411560, 0.105850, -0.019040, -0.042800, -0.433650, 
0.058840, 0.051800, 0.177040, -0.152840, -0.140690, -0.011270, 
-0.000390}
    },
    { /* unit 171 (Old: 171) */
      0.0, 0.164200, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {1.160640, -1.122250, -0.185730, -1.126690, 0.886990, 
-1.339600, -0.197060, -0.606530, 0.462220, -1.446650, 0.059780, 
-0.338810, 0.350660, -0.810600, -0.519830, -0.404000, -0.262080, 
-0.373020, -0.389300, -0.549960, -0.532320, -0.022200, -0.358240, 
-0.359900, -0.495280, 0.215370, -0.386210, -0.246230, -0.362060, 
0.464330, -0.113410, -0.189180, -0.276030, 0.051800, 0.343840, 
-0.236980, -0.146630, 0.152610, 0.358640, 0.087660, 0.039940, 
-0.195100, 0.277180, -0.057560, 0.655400, -0.064860, 0.120870, 
0.063180, 0.941880, 0.245220, -0.325040, -0.047760, 0.684510, 
0.746320, -0.942770, 0.092880, 0.377630, 1.005970, -1.691580, 
0.275200, -0.110320, 0.920610, -2.124980, 0.654900, -0.510940, 
1.092560, -2.332580, 0.900680, -0.798200, 1.522010, -2.384210, 
0.990410, -0.643440, 1.222810, -2.202230, 1.583800, -0.868080, 
1.278960, -1.795720, 1.605700, -0.600040, 1.143640, -1.406790, 
1.332570, -0.374720, 1.306400, -1.043950, 1.126900, -0.256330, 
1.082170, -0.303260, 0.741060, -0.332940, 0.729610, 0.231430, 
0.178250, -0.136940, 0.646630, 0.791990, 0.094680, 0.080760, 
0.008930, 1.026380, -0.042730, 0.014210, -0.562990, 1.137660, 
0.148280, -0.010420, -0.193570, 0.765270, 0.296480, 0.212940, 
-0.434810, 0.207660, 0.210110, 0.140540, -0.129230, -0.190030, 
-0.009200, 0.081240, -0.098460, -0.540780, -0.481790, 0.091410, 
-0.131420, -0.604270, -0.355140, -0.291960, 0.090140, -0.598060, 
-0.716650, -0.176850, 0.323080, -0.306250, -0.892600, -0.101940, 
0.304960, 0.034440, -0.793000, 0.535590, 0.158400, 0.163340, 
-0.353600}
    },
    { /* unit 172 (Old: 172) */
      0.0, 0.298030, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.149730, -0.079770, -0.693240, -1.166250, -0.243760, 
-0.142110, -0.331300, -0.593800, 0.131830, -0.017060, -0.259850, 
-0.335140, -0.084760, 0.116170, -0.051060, -0.145490, -0.040780, 
0.371940, -0.235640, 0.107450, -0.252020, 0.230700, -0.167560, 
0.165730, -0.382000, 0.403740, -0.335640, -0.021900, -0.183910, 
0.448800, -0.037990, -0.088380, -0.214140, 0.385680, -0.063430, 
-0.214860, -0.236350, 0.186350, -0.089930, -0.488690, 0.070780, 
0.268860, 0.251710, -0.152650, -0.014020, 0.046200, -0.008940, 
-0.095010, 0.270350, -0.049720, -0.429710, -0.101070, 0.202850, 
-0.188430, -1.167250, -0.257120, 0.802750, -0.614860, -1.303010, 
0.160260, 0.829270, -0.580590, -1.673480, 0.128630, 0.862640, 
-0.356400, -1.260070, 0.438730, 0.637950, -0.203030, -1.289220, 
0.942440, 0.639330, -0.013910, -1.244250, 1.237650, 0.491040, 
0.037890, -1.086660, 1.535480, 0.633200, 0.073720, -1.026630, 
1.445560, 0.686050, 0.371650, -0.659890, 1.009270, 0.633780, 
-0.084370, -0.818360, 0.910290, 0.595830, -0.010130, -0.945630, 
0.407180, 0.136010, -0.007420, -1.119620, 0.075100, 0.101540, 
0.107100, -1.083460, 0.124570, 0.037620, 0.222420, -0.900820, 
0.113630, 0.251160, 0.170340, -0.492850, 0.104450, 0.566910, 
0.475310, -0.245310, 0.037400, 0.660560, 0.401650, -0.153100, 
0.383390, 0.462570, 0.335580, 0.051810, 0.215770, 0.339830, 0.047710, 
0.130090, 0.013040, -0.092830, 0.140930, 0.314180, -0.428900, 
-0.377090, 0.039540, 0.273840, -0.577550, -0.543200, -0.303160, 
0.548930, -0.755370, -0.759390, -0.354050, 0.619370, -0.239100}
    },
    { /* unit 173 (Old: 173) */
      0.0, -0.311260, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.055010, -0.502040, -0.329300, -0.285600, -0.211880, 
-0.320530, 0.107280, -0.201800, -0.174800, -0.338770, 0.076630, 
-0.186790, 0.013210, -0.182330, -0.097460, -0.046140, -0.063680, 
-0.217050, -0.010000, -0.024540, -0.270240, 0.038710, 0.035470, 
-0.083150, -0.245630, -0.057500, 0.095670, -0.102460, -0.211560, 
-0.328260, 0.093870, 0.023610, -0.201290, -0.423610, 0.233090, 
0.124330, -0.272210, -0.516340, 0.026760, 0.084710, -0.208230, 
-0.255880, -0.113650, -0.033560, -0.235630, -0.183110, -0.073720, 
0.064830, 0.083920, -0.040470, 0.103010, 0.017110, 0.021370, 
-0.055550, -0.102460, -0.170370, 0.005420, 0.171470, 0.124420, 
-0.139640, 0.079050, 0.075390, -0.191470, -0.188030, -0.128170, 
0.484990, -0.092450, -0.120360, 0.068590, 0.190860, 0.095330, 
-0.022180, 0.134070, 0.262000, -0.044950, 0.072750, 0.083350, 
0.288720, 0.018390, -0.145510, 0.099850, 0.135800, -0.112260, 
0.063060, -0.146160, 0.130000, -0.075230, -0.158970, -0.266550, 
0.292900, 0.145920, -0.176510, -0.122440, 0.328870, -0.050210, 
-0.299250, -0.311520, -0.036610, -0.082660, -0.050520, -0.160330, 
0.040510, 0.110850, -0.140230, 0.054050, -0.192880, -0.201720, 
0.091810, 0.170240, -0.063810, -0.009050, -0.271370, -0.020940, 
-0.372530, -0.103090, -0.323410, 0.111090, -0.501950, -0.240680, 
-0.385490, 0.126860, -0.635970, -0.418020, -0.325170, 0.337480, 
-0.526240, -0.435770, -0.475080, 0.281260, -0.564490, -0.296470, 
-0.100510, 0.322480, -0.358400, -0.378980, -0.253200, -0.064340, 
-0.457550, -0.150720, 0.050180, 0.159930, -0.376010, -0.371790, 
-0.079600}
    },
    { /* unit 174 (Old: 174) */
      0.0, -0.353810, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.279500, 0.110620, 0.097000, 0.173130, -0.222560, -0.079180, 
-0.166100, 0.231300, -0.292630, -0.099790, 0.143100, -0.027310, 
-0.074890, -0.238760, -0.103020, -0.081460, -0.127130, -0.302930, 
-0.123300, 0.000180, -0.087500, -0.221550, -0.231800, 0.154340, 
0.008860, -0.126920, -0.121230, -0.188550, -0.050050, -0.026300, 
-0.144200, -0.142030, -0.273390, 0.082710, -0.208060, -0.228420, 
-0.281420, 0.115870, -0.071340, -0.233440, -0.045620, 0.057550, 
0.003540, 0.071300, -0.023420, 0.183010, 0.078140, 0.080760, 
-0.100100, 0.003030, 0.035040, 0.085610, 0.104250, 0.030770, 
-0.204800, -0.125170, 0.089620, -0.102740, -0.183820, 0.078630, 
-0.123620, -0.124070, -0.126070, -0.205360, -0.234560, -0.287910, 
-0.119200, 0.021540, -0.042190, -0.249790, -0.182340, 0.109440, 
0.056360, -0.033360, 0.089210, 0.066650, 0.008410, -0.078770, 
-0.168450, -0.196650, -0.168250, -0.094820, -0.250950, -0.043260, 
0.098260, -0.156910, -0.265970, -0.040730, 0.065600, -0.043290, 
0.109100, -0.191400, 0.107700, -0.105140, -0.027800, -0.208560, 
-0.139880, -0.202660, 0.097520, 0.138980, 0.152340, 0.132510, 
-0.081920, 0.020140, 0.019440, 0.017650, -0.250540, 0.050800, 
-0.263730, 0.082720, -0.142010, -0.161720, -0.008930, 0.015150, 
0.034370, -0.180780, -0.233170, 0.007190, -0.241690, 0.179690, 
-0.230020, -0.324540, -0.137190, 0.075990, -0.167550, 0.003910, 
-0.112020, 0.018700, -0.020500, 0.022580, 0.018320, 0.126740, 
-0.096850, 0.101700, -0.239300, 0.041670, -0.090430, 0.046830, 
-0.193520, -0.169790, -0.109770, 0.208050, 0.101050, -0.157860}
    },
    { /* unit 175 (Old: 175) */
      0.0, -0.292140, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.173140, 0.063790, -0.092040, -0.329610, -0.089710, 
-0.215140, -0.323100, 0.398960, -0.484220, -0.129600, -0.341580, 
0.809860, -0.371060, 0.241120, -0.615180, 0.812730, -0.471430, 
0.276690, -0.475430, 0.927600, -0.270470, 0.458520, -0.265440, 
0.862280, -0.430120, 0.585190, -0.147900, 0.853010, -0.423670, 
0.482170, 0.113100, 0.847360, -0.563530, 0.206340, 0.409530, 
0.566390, -0.419250, 0.128930, 0.476860, 0.645060, -0.502150, 
-0.083620, 0.384850, 0.332580, -0.635500, -0.215550, 0.247900, 
0.160370, -0.355490, -0.317950, 0.220140, 0.111140, -0.417630, 
0.077520, 0.136170, 0.309240, -0.373770, -0.113550, -0.286950, 
0.482550, -0.618280, -0.139440, -0.509210, 0.395700, -0.572320, 
0.098690, -0.591420, 0.496420, -0.645210, 0.230170, -0.656330, 
0.590210, -0.575570, 0.144970, -0.544470, 0.815010, -0.746250, 
0.313350, -0.530140, 0.442250, -0.774100, 0.441270, -0.399920, 
0.651970, -0.650700, 0.425480, -0.189110, 0.521600, -0.792640, 
0.442150, 0.025440, 0.395030, -0.617610, 0.586560, 0.151180, 
0.188030, -0.302320, 0.293580, 0.229320, 0.299250, -0.384390, 
0.114540, 0.217120, -0.112510, -0.335470, -0.219380, -0.022220, 
-0.002880, -0.148900, 0.274640, 0.157810, 0.182450, 0.190320, 
0.302990, -0.218080, 0.007490, -0.008020, 0.241300, -0.453360, 
-0.050100, -0.044240, 0.279920, -0.294680, 0.033560, 0.142960, 
0.015060, -0.141200, -0.184070, -0.096950, 0.018420, -0.134100, 
-0.451380, -0.218120, -0.077700, 0.117810, -0.367660, -0.148120, 
0.113740, 0.165430, -0.173830, -0.168590, -0.086460, 0.249310, 
-0.239800}
    },
    { /* unit 176 (Old: 176) */
      0.0, -0.329190, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.282300, 0.161360, 0.359730, -0.170690, -0.410300, 0.000250, 
0.049810, -0.035550, -0.434540, -0.132370, -0.203630, -0.053570, 
-0.336330, -0.392820, 0.000390, -0.233930, -0.040450, -0.454940, 
-0.096310, -0.046340, -0.217960, -0.659570, -0.080150, -0.029700, 
0.031650, -0.692860, 0.196970, 0.383710, 0.181430, -0.638160, 
0.425650, 0.072930, 0.008070, -0.564580, 0.380790, 0.283060, 
-0.036650, -0.213370, 0.470480, 0.042510, -0.059450, -0.299690, 
0.355900, 0.255470, -0.119210, -0.042910, 0.259200, 0.133150, 
-0.018660, 0.133650, 0.202180, -0.010220, -0.219750, -0.012430, 
0.423450, 0.011200, 0.277670, 0.017230, 0.624280, 0.293930, 0.250670, 
0.203210, 0.599370, 0.266570, 0.433710, -0.108290, 0.346860, 
0.513790, 0.518730, 0.080730, 0.521790, 0.436840, 0.257780, 
-0.218980, 0.379770, 0.505260, 0.134250, -0.202260, 0.050240, 
0.237880, 0.353920, -0.266260, 0.238170, 0.291880, 0.150990, 
-0.047880, 0.178660, 0.260660, 0.224550, 0.102020, -0.078060, 
0.365080, 0.011880, -0.054090, -0.337340, 0.222610, 0.265730, 
-0.050600, -0.211510, 0.074860, 0.249570, -0.188390, -0.038560, 
0.115920, 0.092230, 0.044500, -0.192890, -0.278560, -0.031170, 
-0.086360, -0.148480, -0.293720, -0.040370, -0.113100, -0.162090, 
-0.462060, -0.319570, -0.053760, -0.315700, -0.328500, -0.010510, 
-0.255390, -0.281980, -0.572130, -0.016900, -0.084530, -0.197370, 
-0.250300, -0.033670, -0.205340, -0.412450, -0.111420, 0.101980, 
-0.055010, -0.459100, -0.285480, 0.246230, 0.017650, -0.328490, 
-0.126970, 0.040690, -0.171820, -0.504920, -0.011430}
    },
    { /* unit 177 (Old: 177) */
      0.0, 0.277070, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.116460, -0.106730, 0.284730, -0.051890, -0.231560, 
-0.263400, -0.216780, 0.023600, -0.339310, -0.122010, -0.043080, 
0.441670, -0.335180, -0.042700, -0.289880, 0.307570, -0.430610, 
-0.273340, 0.075430, 0.600080, -0.033990, -0.359240, 0.023680, 
0.708890, -0.042620, -0.133070, 0.402760, 0.413250, 0.127370, 
-0.139690, 0.530460, 0.317440, 0.228580, -0.018760, 0.614020, 
0.134790, 0.017410, -0.168630, 0.261640, 0.200630, -0.262600, 
0.042950, 0.451240, 0.273040, -0.376780, -0.173140, 0.063720, 
-0.112010, -0.257890, 0.180950, -0.014350, 0.044110, -0.331210, 
0.104750, 0.233740, 0.256660, -0.443600, -0.254620, 0.206920, 
0.152840, -0.324370, -0.322940, 0.065680, 0.259090, -0.238050, 
-0.232930, 0.454610, 0.308570, -0.369870, -0.589940, 0.357020, 
0.253250, -0.215290, -0.568860, 0.473440, 0.078000, -0.119490, 
-0.793520, 0.420090, -0.212670, -0.558560, -0.689930, 0.131340, 
-0.269860, -0.310260, -0.723150, 0.039090, -0.269080, -0.370630, 
-0.318760, -0.160770, -0.192850, -0.198400, -0.218170, -0.179870, 
-0.029890, -0.038150, 0.112950, -0.256250, -0.059970, -0.078560, 
0.163480, -0.476760, 0.201740, 0.276930, -0.147330, -0.558310, 
0.056940, 0.268670, 0.328170, -0.490300, 0.074420, 0.097970, 
0.430000, -0.595290, 0.106580, -0.084320, 0.543310, -0.488000, 
0.230620, 0.257000, 0.219660, -0.421130, 0.208600, 0.125670, 
0.121730, -0.698640, 0.366990, 0.330150, 0.140090, -0.689920, 
0.473360, 0.379970, -0.043270, -0.474900, 0.346200, 0.351500, 
0.130980, -0.376410, 0.162570, 0.372880, -0.297670, -0.549510, 
0.248590}
    },
    { /* unit 178 (Old: 178) */
      0.0, -0.086850, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.680110, -0.947990, -0.281890, 0.308520, -0.464350, 
-0.861330, -0.251340, 0.370150, -0.073720, -0.780320, 0.025100, 
0.059000, 0.180380, -0.893580, -0.082600, -0.038130, 0.485280, 
-0.784220, -0.000580, 0.044280, 0.062790, -1.036040, 0.075560, 
0.060120, 0.140540, -0.892200, -0.226340, -0.328280, 0.102660, 
-1.121620, -0.516190, -0.163660, -0.300130, -1.064220, -0.296990, 
-0.290010, -0.589940, -0.904850, -0.449670, -0.317850, -0.537190, 
-1.036560, -0.350850, -0.460950, -0.442310, -0.915120, -0.563850, 
-0.416220, -0.757670, -0.640110, -0.450120, -0.213360, -0.489850, 
-0.521830, -0.428010, -0.441730, -0.257030, -0.460270, -0.574480, 
-0.165760, 0.137780, -0.629720, -0.350670, -0.125680, 0.240690, 
-0.688600, -0.414190, -0.313670, 0.073570, -0.742500, -0.284370, 
-0.141430, 0.102650, -0.380820, -0.213050, 0.206020, 0.184710, 
-0.542350, -0.426570, 0.296440, 0.261880, -0.186520, -0.428490, 
0.292030, 0.125400, -0.334340, -0.304930, 0.217770, 0.027200, 
-0.278710, -0.140620, -0.040230, -0.043270, -0.271680, -0.251780, 
-0.034740, -0.205480, -0.295130, -0.099060, -0.266430, 0.037190, 
0.234500, -0.114710, -0.266960, -0.067310, 0.293680, -0.039290, 
-0.192070, -0.092580, -0.145180, 0.424200, -0.153870, -0.131920, 
-0.274980, 0.189810, 0.096460, -0.073670, -0.400770, 0.407830, 
0.202220, -0.428570, -0.144350, 0.536830, 0.339570, -0.256020, 
-0.256720, 0.413690, 0.274950, -0.644710, 0.001460, 0.343510, 
0.322280, -0.605660, 0.087140, 0.234730, 0.326820, -0.441440, 
0.146830, 0.207190, 0.204550, -0.654810, -0.032610, 0.143390, 
0.003820}
    },
    { /* unit 179 (Old: 179) */
      0.0, -0.215780, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.082130, -0.102800, -0.154940, -0.088250, -0.034370, 
0.093350, 0.068810, -0.032220, 0.059030, -0.002570, -0.215910, 
-0.048250, 0.044060, -0.177190, -0.188120, -0.164680, -0.062450, 
-0.074380, -0.048470, -0.199910, -0.097700, -0.147730, -0.162430, 
-0.028620, -0.069310, 0.147680, -0.207590, 0.007770, -0.080740, 
0.030070, -0.249610, -0.017820, -0.216630, 0.139690, -0.106220, 
-0.225510, 0.157250, 0.140270, -0.189430, -0.141500, -0.057380, 
-0.108100, 0.119000, -0.140320, 0.088340, 0.152100, -0.157590, 
-0.158890, -0.040400, 0.037190, 0.070560, 0.042880, -0.177670, 
0.195210, -0.106580, 0.025210, -0.279660, 0.054120, 0.242620, 
-0.191940, -0.111080, -0.095170, -0.091270, 0.055690, -0.265050, 
-0.157150, -0.041560, -0.025080, -0.235410, -0.021440, -0.004690, 
-0.255180, -0.192690, -0.074500, -0.124890, -0.275380, -0.056010, 
0.027990, -0.141140, -0.402220, -0.203920, -0.159170, -0.035090, 
-0.402550, 0.020620, -0.100430, -0.072230, -0.351380, -0.164890, 
-0.138060, -0.058660, -0.082970, -0.030270, -0.043830, -0.126760, 
0.065080, -0.087630, -0.072750, -0.001010, 0.025800, -0.142350, 
0.103020, -0.078350, -0.013960, 0.043050, 0.024280, -0.088940, 
-0.269480, -0.163290, 0.167510, -0.160430, 0.032850, 0.056230, 
-0.142890, -0.094640, -0.213180, -0.001240, 0.154540, 0.079040, 
-0.277750, -0.047140, -0.133770, 0.052330, 0.025590, -0.041920, 
-0.174010, -0.002220, -0.204110, 0.065670, -0.017100, 0.008250, 
0.026980, 0.191490, -0.064980, -0.196080, -0.173850, 0.106450, 
0.094640, -0.078820, -0.211000, 0.154770, 0.222320, -0.101250, 
-0.012100}
    },
    { /* unit 180 (Old: 180) */
      0.0, -0.710130, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.387110, -0.469790, 0.078910, -0.248990, -0.109160, 
-0.407130, 0.322260, 0.113610, -0.014400, -0.547560, 0.456260, 
-0.036750, 0.261000, -0.386540, 0.148480, -0.005960, 0.038640, 
-0.384200, 0.001990, -0.223510, -0.002300, -0.521580, 0.014890, 
-0.345050, -0.060190, -0.306110, -0.136930, -0.277260, 0.028860, 
-0.613220, -0.261620, -0.081090, 0.117520, -0.727020, -0.330280, 
0.091850, 0.124830, -0.859420, -0.126330, 0.361390, 0.284200, 
-0.552260, -0.200860, 0.227340, 0.213180, -0.346940, -0.229360, 
-0.064200, 0.452760, -0.102310, -0.136420, 0.004650, 0.350940, 
0.381440, -0.328590, 0.021610, 0.192080, 0.918450, -0.384280, 
-0.159310, -0.258110, 1.147340, -0.633950, -0.179710, -0.132560, 
1.415070, -0.346830, -0.257640, -0.192550, 1.393320, -0.598930, 
-0.010330, -0.411540, 1.310480, -0.770250, 0.101620, -0.006210, 
1.257810, -0.626390, 0.009940, -0.019480, 1.047620, -0.472450, 
-0.000950, -0.156400, 0.912370, -0.227810, -0.261360, -0.272310, 
0.531640, -0.104680, -0.373980, -0.234870, 0.484480, 0.064880, 
-0.180460, -0.691630, 0.376390, 0.001040, -0.127400, -0.801560, 
0.160000, -0.012260, 0.010460, -0.719120, 0.254600, -0.197560, 
-0.195480, -0.273020, -0.009080, -0.100320, -0.141870, -0.311070, 
-0.599290, -0.454200, 0.074300, 0.107280, -0.681290, -0.345340, 
-0.240000, 0.093940, -0.724840, -0.645470, -0.073250, 0.169460, 
-1.069800, -0.494730, -0.236230, -0.095420, -1.240780, -0.661490, 
-0.299840, 0.143580, -1.018610, -0.384940, -0.366590, -0.150640, 
-0.929610, -0.563780, -0.121890, -0.204720, -0.954690, -0.558800, 
0.212030}
    },
    { /* unit 181 (Old: 181) */
      0.0, -0.175480, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.096320, -0.117850, 0.429060, -0.202410, -0.044890, 
0.031520, 0.341210, 0.085800, -0.331530, 0.068340, 0.529370, 
-0.087840, -0.294370, 0.094230, 0.396190, 0.122830, -0.261670, 
-0.195010, 0.325160, 0.165060, 0.070920, 0.074290, 0.102250, 
-0.105370, -0.001570, -0.205220, 0.215160, -0.088650, -0.219140, 
0.148400, 0.041690, -0.316880, -0.049260, 0.084140, -0.112160, 
-0.265750, 0.009320, 0.227880, 0.136980, -0.266000, -0.026090, 
0.214350, 0.241280, -0.224110, 0.189670, 0.255520, 0.189340, 
-0.239030, -0.141120, -0.118230, 0.216630, 0.075050, -0.011860, 
-0.146620, -0.125140, -0.034620, -0.041620, -0.298350, -0.089000, 
-0.180740, -0.028080, -0.285400, -0.191760, -0.024970, 0.166490, 
-0.170060, 0.015830, -0.182680, 0.019640, -0.297170, -0.259940, 
-0.206350, 0.090210, -0.113510, -0.024690, -0.075540, 0.053050, 
-0.050170, -0.015920, -0.079350, 0.114840, 0.125260, -0.000380, 
0.077110, -0.057980, 0.045720, -0.089080, 0.197710, 0.085280, 
-0.232450, -0.125430, -0.003380, 0.369970, -0.249250, 0.101980, 
0.215720, 0.064180, -0.385640, 0.092710, 0.091450, 0.299320, 
-0.309690, 0.155690, -0.183770, 0.195800, -0.092030, -0.080800, 
0.008310, 0.171910, -0.147660, -0.057860, 0.315410, -0.317620, 
-0.263740, -0.205760, 0.299530, -0.410570, -0.384650, -0.229400, 
0.230920, -0.342660, -0.354870, -0.305030, 0.229270, -0.328230, 
-0.237170, -0.397850, 0.297160, -0.584870, -0.324770, -0.403500, 
0.075430, -0.404870, -0.011070, -0.022940, -0.140800, -0.340040, 
-0.274060, -0.053590, -0.158050, -0.392990, 0.071680, 0.299590, 
-0.312520}
    },
    { /* unit 182 (Old: 182) */
      0.0, -0.298740, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.106070, -0.081480, 0.200160, 0.055770, 0.152860, -0.061670, 
0.089230, -0.144160, -0.076870, -0.018110, 0.012990, -0.138030, 
-0.078780, -0.202380, 0.082810, 0.142100, -0.306740, -0.071850, 
-0.252480, 0.078260, -0.187660, -0.263530, -0.060880, 0.202210, 
-0.104670, -0.025940, -0.085660, 0.085050, -0.208590, -0.020700, 
-0.159590, 0.139680, 0.121590, -0.000610, -0.118600, -0.286290, 
0.046800, 0.016150, 0.246200, -0.205110, -0.121940, -0.040890, 
0.131140, -0.023590, -0.074110, -0.209970, -0.072890, -0.159590, 
0.082540, 0.143080, 0.199170, -0.055280, 0.074140, -0.216010, 
-0.032880, -0.040600, -0.008970, 0.099850, -0.016440, 0.026970, 
0.139280, -0.198580, -0.280420, 0.170710, 0.008190, -0.122520, 
-0.134160, 0.352720, 0.033000, -0.086080, -0.253650, 0.213300, 
-0.004450, -0.351710, -0.035180, 0.273070, 0.087260, -0.349310, 
-0.210240, 0.024250, -0.076740, -0.129310, -0.088410, -0.054980, 
0.052200, -0.211330, 0.107110, 0.066300, 0.179550, -0.274970, 
-0.000550, 0.022570, 0.291310, -0.082320, -0.252700, 0.168920, 
0.190700, -0.200630, -0.142240, 0.119300, 0.045910, 0.107260, 
-0.151740, 0.075630, 0.156430, 0.136330, -0.076740, -0.233830, 
-0.100350, 0.117010, 0.033590, 0.028090, -0.057250, 0.005990, 
-0.077430, 0.123040, -0.170990, -0.024270, -0.349880, -0.059440, 
-0.078390, -0.091960, -0.251330, -0.201310, -0.125740, -0.102990, 
-0.281920, -0.033090, -0.030110, -0.049140, -0.049740, -0.243180, 
0.002760, 0.160220, -0.046100, -0.310080, 0.075710, -0.007920, 
-0.130200, -0.238560, 0.136190, 0.047640, -0.058940, -0.137230}
    },
    { /* unit 183 (Old: 183) */
      0.0, -0.347100, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.419990, -0.103420, 0.104690, 0.262090, -0.166980, 0.041930, 
-0.045810, 0.182920, 0.025300, 0.157050, 0.037640, 0.240490, 
0.229560, 0.311850, -0.002730, 0.076500, 0.214750, -0.001300, 
0.155100, -0.336230, 0.343670, -0.112770, -0.095000, -0.365340, 
-0.132580, -0.308300, -0.090840, -0.454100, 0.125550, -0.497960, 
-0.437950, -0.114700, 0.057980, -0.854870, -0.227980, 0.283730, 
0.046760, -1.081970, -0.323750, 0.175180, -0.083440, -0.995170, 
-0.174380, 0.206820, 0.146800, -0.520760, -0.123830, -0.174740, 
-0.128850, -0.388700, -0.272220, -0.309710, -0.009490, 0.082760, 
-0.304310, -0.287070, -0.302270, 0.892630, -0.597970, -0.460360, 
-0.285220, 1.303290, -0.707280, -0.277440, -0.718570, 1.597410, 
-0.510350, -0.244700, -0.478990, 1.486570, -0.536770, -0.152980, 
-0.383330, 1.400230, -0.574530, 0.165120, -0.349570, 1.578660, 
-0.473060, 0.042260, 0.040980, 1.161790, -0.492760, 0.082520, 
-0.311390, 1.063280, -0.137730, -0.112330, -0.460820, 0.756400, 
0.132180, 0.140090, -0.680840, 0.723560, 0.069880, -0.053910, 
-0.443140, 0.603100, -0.096160, -0.081280, -0.549910, 0.288250, 
0.068050, -0.020620, -0.421850, 0.187220, -0.144250, 0.281680, 
-0.423910, -0.124470, -0.158150, -0.065660, -0.089300, -0.609670, 
-0.327160, -0.045310, -0.128600, -1.166590, -0.306640, -0.198000, 
-0.108840, -1.311120, -0.206970, -0.101070, 0.124570, -1.414200, 
-0.407200, -0.302910, 0.158580, -1.500530, -0.551050, -0.378660, 
-0.171840, -1.531150, -0.354440, -0.167650, -0.081010, -1.323550, 
-0.284340, -0.204850, -0.262550, -0.966720, -0.210370, 0.223840}
    },
    { /* unit 184 (Old: 184) */
      0.0, -0.267230, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.112910, 0.142310, 0.321640, 0.275980, 0.318390, -0.029960, 
0.071970, 0.195240, 0.362740, -0.142310, -0.133240, -0.182760, 
0.136010, 0.145560, -0.088140, -0.183800, 0.096840, 0.037290, 
0.036540, -0.085700, 0.108420, 0.040320, 0.009460, 0.024170, 
-0.023990, -0.069040, 0.017740, -0.036930, -0.223610, -0.296190, 
0.036010, -0.114130, -0.161090, -0.273460, -0.026030, -0.018780, 
0.087060, -0.053630, -0.004630, -0.018770, 0.056720, -0.098960, 
-0.131140, 0.023950, -0.200890, 0.077040, -0.081710, 0.298630, 
0.137750, 0.192670, -0.134470, 0.094160, -0.044410, 0.150080, 
-0.095990, -0.011670, -0.278080, 0.200100, -0.141970, 0.244360, 
-0.065240, -0.157700, 0.058390, 0.163210, -0.227510, -0.260500, 
-0.222570, -0.134030, -0.561190, -0.165630, -0.379870, -0.232440, 
-0.499350, -0.186940, -0.274430, -0.290820, -0.676990, -0.442060, 
-0.302980, -0.302920, -0.708660, -0.196830, -0.261200, -0.052310, 
-0.359200, -0.294340, -0.346760, 0.005680, -0.340030, -0.409440, 
-0.074890, 0.059750, -0.141410, -0.183540, -0.105140, -0.017290, 
-0.028960, -0.430740, -0.153930, 0.227180, 0.087340, -0.211170, 
0.087630, 0.174270, -0.139840, -0.388380, 0.012880, -0.094340, 
-0.161900, -0.246120, 0.169090, -0.076160, 0.059050, -0.191310, 
0.180920, 0.006560, -0.176990, -0.257690, 0.320590, -0.243880, 
-0.162790, -0.150800, 0.102990, -0.178760, 0.081850, -0.012490, 
0.104430, -0.153870, -0.142670, -0.135010, 0.378520, 0.132770, 
0.041340, -0.210770, 0.111910, 0.146250, 0.077930, 0.104330, 
0.204700, 0.261850, 0.235570, -0.076810, -0.058610, 0.162340}
    },
    { /* unit 185 (Old: 185) */
      0.0, -1.256630, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.450100, 0.770640, 0.028010, 0.686670, -0.154960, 0.590280, 
-0.167180, 0.364620, -0.201860, 0.504480, -0.490700, -0.044880, 
-0.538120, -0.235730, -0.745570, -0.371780, -0.600400, -0.278160, 
-0.795100, -0.144730, -0.367870, -0.836370, -0.839310, -0.003350, 
0.032740, -0.692370, -0.646430, 0.389480, 0.292360, -0.768180, 
-0.440050, 0.148140, 0.571740, -0.595450, -0.520720, 0.293600, 
0.742430, -0.461230, -0.370080, 0.614890, 0.727350, -0.157860, 
-0.022930, 0.449910, 0.589110, 0.012800, 0.101000, 0.449150, 
0.631900, 0.361950, 0.515160, 0.211750, 0.219950, 0.553180, 1.084620, 
-0.032040, 0.124040, 0.718490, 1.280270, 0.061900, -0.301230, 
0.930550, 1.566710, -0.470180, -0.734470, 0.769590, 1.158370, 
-0.341130, -1.050930, 0.766780, 1.164760, -0.809570, -1.129710, 
0.871650, 1.166670, -0.754980, -1.189680, 0.673690, 0.949160, 
-0.959810, -0.904900, 0.449580, 1.115580, -1.039350, -0.889810, 
0.388650, 1.098330, -0.722350, -0.532810, 0.569190, 0.817420, 
-0.352670, -0.766530, 0.743470, 0.841790, -0.046250, -0.393150, 
0.349450, 1.027920, -0.033670, 0.010730, 0.206230, 0.842540, 
0.327640, -0.069770, -0.267460, 0.461950, 0.079060, -0.042300, 
-0.521050, 0.240410, -0.120680, -0.518830, -0.670670, -0.174700, 
-0.011410, -0.220220, -0.908240, -0.409360, -0.304930, -0.319050, 
-0.842440, -0.191330, -0.197330, 0.192320, -0.541010, -0.457790, 
-0.057440, 0.174200, -0.116680, -0.192950, 0.203900, 0.530470, 
0.228020, -0.206170, 0.231320, 0.695680, 0.402230, -0.219240, 
0.791440, 0.628950, 0.849770, 0.002560, 0.510990}
    },
    { /* unit 186 (Old: 186) */
      0.0, -0.471820, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.027960, 0.178750, -0.114970, -0.302120, 0.139460, -0.057540, 
-0.202780, 0.066820, -0.034280, -0.265620, -0.040690, 0.042870, 
0.009240, -0.130240, -0.249550, 0.002670, -0.131270, -0.205200, 
-0.267810, 0.025090, -0.042390, -0.184270, 0.040690, -0.128490, 
-0.139420, 0.020070, 0.002960, -0.000060, -0.294390, -0.119070, 
-0.085510, -0.119500, 0.040850, -0.127750, 0.236980, -0.145180, 
-0.253020, 0.010570, 0.222860, 0.104970, 0.044880, 0.080570, 
0.187000, 0.099060, -0.035690, -0.167800, -0.023650, 0.018390, 
-0.115280, 0.176950, -0.104920, -0.034310, -0.194710, -0.193090, 
-0.051380, -0.052080, 0.024920, -0.116870, -0.267080, 0.217590, 
-0.218410, -0.032850, -0.240780, 0.008810, -0.142800, -0.144290, 
-0.200830, 0.035760, -0.009410, 0.086300, -0.163720, 0.347300, 
-0.332770, 0.083640, -0.302710, -0.027250, -0.262690, -0.203460, 
-0.293280, 0.253170, -0.021840, -0.055670, -0.141430, 0.243820, 
-0.243370, -0.023460, -0.053460, 0.042970, 0.042100, -0.045890, 
-0.179890, 0.051380, -0.018040, 0.167710, -0.012240, 0.073960, 
0.018560, 0.055020, -0.139150, 0.075480, -0.104550, 0.148510, 
-0.031890, 0.086730, -0.033380, -0.107360, 0.082400, -0.113540, 
-0.205200, -0.007550, -0.252870, 0.053030, 0.100350, -0.099970, 
0.009560, 0.036490, -0.184300, 0.191320, -0.185780, -0.283390, 
-0.132950, 0.080440, -0.360530, -0.043770, 0.185510, -0.117280, 
-0.188440, -0.304390, -0.089510, 0.186900, -0.338480, -0.155390, 
-0.010040, 0.094340, 0.021320, -0.307150, 0.166530, 0.021500, 
-0.105820, -0.284030, 0.040050, -0.154690, -0.295450, 0.010790}
    },
    { /* unit 187 (Old: 187) */
      0.0, 0.052940, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.008910, -0.127350, -0.069740, 0.190140, -0.351040, 
-0.052410, 0.161540, 0.111900, -0.212650, 0.045910, 0.010950, 
-0.049390, -0.278870, 0.028280, -0.153080, 0.049720, -0.106660, 
-0.304020, 0.220250, 0.239440, -0.466640, -0.299540, -0.153130, 
0.267140, -0.240950, -0.348940, 0.040340, 0.172040, -0.442580, 
-0.161670, -0.099930, 0.096950, -0.143180, 0.013960, 0.166740, 
0.191060, -0.302810, -0.126830, 0.024470, 0.235520, -0.350230, 
-0.214780, 0.168480, 0.200390, -0.170960, -0.062430, 0.183370, 
0.090450, -0.081880, -0.221700, 0.006730, 0.274290, -0.199060, 
-0.180550, -0.097700, 0.247830, -0.024920, -0.334350, 0.122330, 
0.126120, -0.233110, -0.218910, -0.203990, -0.163650, -0.101100, 
-0.327100, 0.115860, 0.108540, -0.017110, -0.159680, -0.180360, 
-0.305520, -0.151260, -0.345090, -0.395080, -0.385400, 0.061040, 
-0.242230, -0.388250, -0.164180, -0.195070, -0.293640, -0.316960, 
-0.314500, -0.108260, -0.123790, -0.470300, -0.355110, -0.059260, 
-0.104260, -0.282180, -0.157230, -0.060480, -0.427890, -0.173180, 
0.001000, -0.179910, -0.161090, -0.364880, 0.177290, -0.112120, 
-0.039750, -0.192980, 0.268250, -0.197500, -0.059770, -0.353950, 
0.393420, -0.160110, -0.282950, -0.205430, 0.040950, 0.044610, 
-0.158960, 0.056760, 0.031700, -0.225280, -0.179890, 0.104160, 
0.178330, -0.208590, -0.061320, 0.217770, 0.232390, -0.248280, 
0.059130, 0.110100, 0.018140, -0.274600, -0.034030, -0.169110, 
0.217870, -0.275010, -0.121470, -0.086880, 0.243850, -0.260640, 
-0.003900, 0.016620, -0.053720, -0.212830, -0.016450, 0.110160, 
-0.017800}
    },
    { /* unit 188 (Old: 188) */
      0.0, -0.355070, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.153910, -0.050950, -0.128840, 0.028790, 0.016830, 0.035740, 
-0.127500, 0.088480, 0.041380, -0.074250, 0.077800, 0.009470, 
-0.204840, -0.149740, 0.104000, 0.041930, -0.139490, -0.295970, 
0.096530, 0.096120, -0.102360, -0.246430, -0.248510, -0.214290, 
0.042540, -0.245650, -0.022310, -0.117920, 0.069120, -0.001680, 
-0.125840, -0.199710, 0.150900, -0.201810, -0.107040, 0.079430, 
0.057950, -0.047780, 0.108000, -0.052070, -0.016270, 0.043790, 
0.132410, -0.056440, 0.035460, 0.110440, -0.130850, -0.176840, 
-0.134370, 0.054750, -0.057250, 0.184930, 0.160460, -0.022510, 
-0.170470, -0.161640, 0.056170, -0.087990, 0.028850, -0.145080, 
-0.036220, 0.002990, -0.058470, 0.071540, -0.283110, -0.025840, 
0.012290, 0.013330, -0.305580, 0.046080, -0.077740, -0.106230, 
-0.308840, -0.047290, -0.118600, 0.009840, -0.282140, 0.265780, 
0.105530, -0.233430, -0.182830, -0.065250, -0.156060, -0.068220, 
-0.076070, -0.062830, 0.177560, -0.234800, -0.117000, 0.091630, 
0.130200, -0.068310, 0.012050, -0.107210, 0.076950, 0.082160, 
0.053930, 0.194420, -0.031610, 0.069370, -0.327720, 0.150110, 
-0.109680, -0.202460, -0.178530, 0.100690, 0.092960, -0.066650, 
0.024130, 0.079290, -0.100370, -0.195210, -0.136790, 0.074280, 
-0.059920, -0.026910, -0.127630, -0.195640, -0.104370, -0.178670, 
-0.230570, -0.089530, 0.022760, 0.078020, -0.195690, -0.272550, 
-0.127580, 0.121050, 0.032630, -0.143780, 0.040430, -0.118320, 
-0.016020, -0.269780, -0.183160, -0.089890, 0.014150, -0.212000, 
-0.231990, -0.122210, -0.116940, -0.027050, 0.068690, 0.021100}
    },
    { /* unit 189 (Old: 189) */
      0.0, -0.415920, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.157750, -0.033640, -0.231530, -0.137780, -0.166160, 
-0.101720, 0.069000, -0.101830, -0.141910, -0.075400, -0.064100, 
0.154210, -0.020280, 0.084570, 0.053430, -0.232090, -0.245590, 
-0.102480, 0.008730, -0.198470, -0.279090, 0.040840, -0.017910, 
-0.115860, -0.021500, -0.168030, -0.045540, -0.045960, -0.164100, 
0.158080, 0.015090, -0.035770, -0.009790, -0.018650, -0.090520, 
-0.192100, -0.013520, 0.122390, -0.292000, -0.057440, -0.128620, 
-0.046540, -0.187900, 0.084260, 0.094640, -0.063410, -0.015290, 
-0.137890, -0.192060, 0.161080, 0.019180, 0.181940, 0.157140, 
-0.184210, 0.121540, -0.164090, -0.058510, -0.245540, -0.005840, 
-0.042260, -0.063070, 0.102630, 0.045040, -0.247340, -0.190920, 
0.070460, 0.138920, 0.043060, -0.326780, -0.232680, -0.157970, 
-0.343200, -0.033190, -0.227860, -0.087510, -0.387370, 0.000300, 
-0.100910, -0.047670, -0.123100, 0.032890, 0.007160, -0.206020, 
-0.096350, 0.007300, 0.027500, -0.169410, -0.305500, -0.251740, 
-0.071750, -0.119960, -0.337760, -0.059270, 0.021400, -0.022890, 
-0.303660, 0.007650, 0.064900, -0.176610, -0.212020, -0.152580, 
0.018180, -0.038770, -0.083000, 0.116250, 0.022850, 0.116920, 
-0.086440, -0.102110, 0.010190, -0.176520, 0.096000, 0.035690, 
-0.156490, -0.086420, -0.082970, 0.083870, 0.155880, 0.072980, 
-0.098950, -0.139960, -0.144380, 0.119460, -0.133170, -0.046000, 
0.066790, 0.116260, -0.206020, 0.005450, 0.085780, -0.113980, 
-0.130320, -0.034560, -0.041710, -0.080100, -0.197070, 0.138540, 
-0.001850, -0.057190, 0.047970, 0.097610, 0.195030, -0.042150, 
0.016360}
    },
    { /* unit 190 (Old: 190) */
      0.0, -0.197810, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.094490, 0.007440, -0.108820, -0.055760, -0.283160, 
0.260680, 0.206140, 0.069780, -0.299810, 0.058560, 0.279750, 
0.039300, -0.086130, 0.031640, 0.000310, -0.006600, -0.155000, 
-0.030290, -0.014410, 0.247990, 0.075000, -0.013430, -0.231230, 
-0.092070, -0.123610, -0.297930, -0.146340, 0.101870, -0.055820, 
-0.042960, -0.189500, 0.042790, -0.181320, 0.221900, -0.103720, 
0.000260, -0.086240, 0.266070, 0.038730, -0.301940, -0.164340, 
0.370020, -0.244140, -0.088790, -0.098240, 0.258660, -0.151380, 
-0.062990, 0.005530, 0.416880, -0.280790, -0.046130, 0.173360, 
0.103420, 0.012930, 0.082860, 0.057440, 0.223920, -0.096830, 
-0.233090, -0.133040, -0.125450, -0.329380, 0.017750, -0.090180, 
-0.345560, -0.405130, -0.035600, -0.180830, -0.319000, -0.178460, 
-0.174620, 0.094500, -0.281440, -0.445790, -0.146560, 0.061010, 
-0.332820, -0.358100, -0.091820, 0.040760, -0.310030, -0.339060, 
0.165510, 0.243280, -0.213850, -0.305490, -0.115730, 0.183310, 
-0.150850, -0.298470, 0.216490, 0.299770, -0.081870, -0.413080, 
0.228640, 0.230680, -0.246560, -0.276830, 0.181390, 0.228790, 
0.200410, -0.205980, 0.028830, -0.204260, -0.027400, -0.430930, 
0.089540, -0.241090, -0.195930, -0.144950, 0.202280, -0.025970, 
-0.174140, -0.237620, 0.250390, -0.486000, -0.316210, -0.235710, 
0.299170, -0.180900, -0.123520, -0.312740, 0.179160, -0.474090, 
0.005200, -0.406190, 0.057220, -0.509390, 0.013690, -0.466170, 
0.173080, -0.397780, -0.086310, -0.239520, 0.172470, -0.274670, 
0.084140, -0.009230, 0.005500, -0.148570, 0.024780, 0.041200, 
-0.185440}
    },
    { /* unit 191 (Old: 191) */
      0.0, -0.730190, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.174100, -0.184650, -0.327700, 0.448460, 0.299740, -0.431840, 
-0.183070, 0.312390, 0.191750, -0.545160, 0.229890, 0.261520, 
0.396480, -0.382340, 0.090550, -0.195490, 0.342760, 0.065550, 
0.059430, -0.521600, 0.097780, 0.242980, 0.280710, -0.457240, 
0.227150, 0.410490, 0.631580, -0.467580, 0.223460, 0.308750, 
0.623420, -0.516240, 0.766090, 0.213350, 0.549740, -0.489520, 
0.736680, -0.042770, 0.173020, -0.701250, 1.002270, -0.128100, 
0.006880, -0.660160, 0.647220, -0.158580, 0.166800, -0.326730, 
0.571850, -0.299780, 0.082980, -0.363320, 0.430940, -0.224350, 
-0.156810, -0.275810, 0.108260, -0.080280, -0.089110, -0.187180, 
-0.373170, -0.209800, -0.261590, 0.121840, -0.566910, -0.322320, 
-0.464980, -0.114180, -0.715430, -0.300330, -0.584770, 0.128620, 
-0.785890, -0.607560, -0.701630, 0.018660, -0.984590, -0.678940, 
-0.301840, -0.111540, -1.381010, -0.456560, -0.510050, 0.127620, 
-1.512140, -0.527390, -0.452050, -0.100970, -1.427000, -0.269680, 
-0.013620, 0.175110, -1.223340, -0.263590, -0.251440, -0.012410, 
-0.697370, 0.212620, -0.056670, -0.005950, -0.693830, 0.137560, 
0.043950, 0.059020, -0.214600, 0.320400, -0.061440, 0.288810, 
-0.561840, 0.752470, 0.087550, 0.104320, -0.328080, 0.939630, 
-0.103020, 0.370670, -0.721120, 0.843450, -0.117120, 0.325260, 
-0.771790, 1.129950, 0.003070, 0.331130, -0.746130, 0.944520, 
-0.105790, 0.639990, -0.810640, 0.654470, -0.166820, 0.686830, 
-0.665480, 0.818820, -0.083420, 0.941560, -0.504560, 0.452350, 
0.068390, 1.089270, -0.267280, 0.463350, -0.224280, 1.147690}
    },
    { /* unit 192 (Old: 192) */
      0.0, -0.253670, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.014240, 0.064710, -0.142470, -0.049720, 0.108210, 
-0.231050, -0.052140, 0.142460, 0.093290, -0.239870, 0.026400, 
-0.069370, -0.088160, -0.292570, 0.087780, -0.209690, -0.138850, 
0.024150, -0.102280, -0.048810, -0.037740, -0.278810, 0.072500, 
-0.016890, -0.159930, -0.257600, 0.034430, 0.001730, -0.257500, 
-0.121460, -0.012670, -0.118700, -0.204210, -0.000730, -0.168070, 
-0.141710, -0.028500, -0.149150, -0.248570, 0.014050, -0.099170, 
0.137040, -0.055240, -0.069570, 0.053250, -0.152860, -0.030270, 
-0.143450, -0.066550, -0.007320, -0.001930, 0.017060, 0.065710, 
-0.164920, -0.083580, 0.073150, -0.160700, 0.129810, 0.173010, 
-0.079840, 0.065940, 0.161840, 0.203430, 0.257140, -0.097910, 
0.115990, -0.038660, -0.015910, 0.079040, 0.003800, 0.103710, 
-0.012740, -0.095660, 0.201220, 0.143330, -0.105280, -0.270800, 
-0.182430, -0.017960, 0.179200, -0.233230, -0.157620, 0.001090, 
0.001220, 0.093800, -0.103960, 0.149500, 0.259900, 0.104570, 
-0.258740, 0.037970, 0.199270, -0.195380, 0.018120, -0.069170, 
0.187050, -0.124770, -0.029760, -0.062240, 0.147110, -0.074530, 
-0.041420, 0.034500, 0.243460, 0.126480, 0.052890, 0.110630, 
-0.149370, -0.208330, 0.067100, -0.185940, -0.173210, 0.013920, 
-0.087920, -0.162670, -0.004700, -0.077760, -0.246660, -0.153380, 
-0.084110, -0.176290, 0.006180, 0.086560, 0.045120, 0.059790, 
-0.079850, 0.075720, 0.052390, -0.042290, -0.094690, 0.120300, 
0.162410, -0.051280, -0.030760, -0.171140, 0.140570, 0.144460, 
-0.060300, 0.037840, 0.305500, 0.161500, 0.093580, -0.150680, 
0.104600}
    },
    { /* unit 193 (Old: 193) */
      0.0, -0.528530, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.086920, 0.213380, -0.183500, -0.224010, -0.364060, 
0.194750, -0.337510, -0.136580, -0.488100, -0.008240, -0.592130, 
-0.131050, -0.438130, -0.151490, -0.683580, -0.063760, -0.278960, 
-0.271140, -0.446070, -0.049990, -0.273800, -0.017350, -0.235100, 
-0.085950, -0.090800, -0.144810, -0.230560, 0.148060, -0.230470, 
-0.213980, -0.097140, -0.015980, 0.010720, -0.212090, -0.041090, 
0.020550, 0.000090, -0.041720, -0.249920, 0.047900, -0.300890, 
-0.054870, 0.060450, -0.033680, -0.279990, -0.184550, -0.018230, 
-0.115240, -0.042930, -0.268540, 0.036850, 0.054000, -0.075450, 
-0.070500, 0.084430, -0.032040, -0.150370, -0.002470, 0.248650, 
-0.142620, 0.132470, 0.066180, 0.516320, 0.235470, -0.025540, 
-0.000140, 0.570320, -0.060730, 0.292350, 0.029690, 0.574280, 
0.191330, 0.030910, -0.126170, 0.447510, -0.048790, 0.356460, 
-0.321360, 0.468090, 0.124080, 0.138730, -0.116730, 0.554460, 
0.018870, 0.099790, -0.141780, 0.156870, 0.092030, 0.382660, 
-0.235040, 0.109800, 0.031930, 0.014580, -0.017530, 0.149320, 
-0.018660, 0.163180, -0.079050, -0.207090, 0.214160, 0.118540, 
-0.142400, -0.066910, -0.027940, 0.037510, -0.057700, 0.024600, 
0.019220, 0.063460, -0.113000, -0.154730, -0.155000, 0.096480, 
-0.191140, 0.046010, 0.083920, 0.186980, -0.398510, -0.213980, 
-0.065750, 0.095450, -0.144160, 0.020010, -0.179820, 0.162970, 
-0.380480, -0.069260, -0.116930, 0.452440, -0.454740, -0.012650, 
-0.012050, 0.305170, -0.344210, -0.192970, 0.021830, 0.050620, 
-0.228370, -0.284450, 0.045200, 0.133570, -0.331570, -0.387700, 
0.110560}
    },
    { /* unit 194 (Old: 194) */
      0.0, -0.667630, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.246200, 0.060170, -0.178750, -0.027620, -0.044470, 
0.016630, 0.077230, -0.219930, 0.028450, 0.078570, 0.132050, 
-0.498280, -0.041070, 0.253620, 0.344700, -0.405890, -0.033300, 
-0.177390, 0.118820, -0.189250, -0.280890, 0.149670, 0.194840, 
0.051320, 0.018680, 0.135740, -0.031800, -0.024270, -0.303210, 
0.052890, -0.208290, -0.110320, -0.324120, 0.439620, -0.161380, 
-0.067320, 0.009600, 0.392050, -0.631850, 0.068550, 0.084160, 
0.484390, -0.608580, 0.142280, 0.008750, 0.440130, -0.523880, 
0.007370, 0.129110, 0.169370, -0.241890, -0.347500, 0.340160, 
-0.167830, -0.158820, -0.134650, 0.364770, -0.453110, 0.407260, 
-0.635130, 0.182560, -0.166340, 0.512780, -0.644120, 0.253050, 
-0.377420, 0.672870, -0.644840, 0.330070, -0.271380, 0.807250, 
-0.766400, 0.391900, -0.384100, 0.747380, -0.669080, 0.682810, 
-0.128580, 0.657520, -0.644020, 0.651260, -0.271380, 0.364830, 
-0.643640, 0.504660, -0.231150, 0.350770, -0.560080, 0.376520, 
-0.238940, 0.207770, -0.665560, 0.391250, -0.281630, -0.149920, 
-0.790680, 0.199620, -0.334390, -0.403580, -0.830970, -0.007540, 
0.073460, -0.120720, -0.703700, 0.041600, 0.196420, -0.311410, 
-0.462790, -0.041410, 0.178930, -0.243580, -0.773000, 0.273660, 
-0.065980, 0.060750, -0.698720, 0.317070, 0.094530, -0.149990, 
-0.855890, 0.223320, 0.506500, 0.102380, -0.670640, -0.053640, 
0.381520, -0.072310, -0.342760, 0.224210, 0.688220, 0.138330, 
-0.460420, 0.241760, 0.469150, -0.117080, -0.178330, 0.008830, 
0.420850, -0.155870, -0.199940, -0.217680, 0.656290, 0.068220, 
0.094400}
    },
    { /* unit 195 (Old: 195) */
      0.0, -0.294020, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.005900, -0.753250, 0.182280, -0.010660, -0.063190, 
-0.254910, 0.035140, -0.017050, 0.309240, -0.331170, 0.177990, 
0.234670, 0.457350, -0.102180, 0.336760, 0.264940, 0.637570, 
-0.338770, 0.233280, 0.024400, 0.244470, -0.315160, 0.219060, 
0.091280, 0.449450, -0.437660, 0.122970, -0.206080, 0.101830, 
-0.267760, -0.101790, 0.040300, 0.153630, -0.319150, -0.099750, 
-0.021090, -0.164920, -0.688670, -0.292940, 0.166710, 0.009020, 
-0.380260, -0.129820, -0.047780, -0.370230, -0.296890, -0.239670, 
0.186610, -0.133980, -0.408390, -0.084200, -0.019160, -0.259360, 
-0.298630, -0.186710, -0.105250, -0.043840, -0.063660, -0.180080, 
-0.469870, 0.045570, 0.047450, -0.050480, -0.468730, 0.136680, 
0.536530, -0.158010, -0.533600, 0.040390, 0.624570, -0.087020, 
-0.415010, 0.251150, 0.779210, 0.025410, -0.176250, 0.513180, 
0.444140, -0.063110, -0.413210, 0.273750, 0.505800, 0.250240, 
-0.257020, 0.154310, 0.539790, 0.001090, -0.492490, -0.111460, 
0.254470, 0.315730, -0.237750, -0.135540, 0.160270, 0.146910, 
-0.317520, -0.309370, 0.095930, 0.187510, -0.446190, -0.208530, 
0.273560, -0.097270, -0.459900, -0.257660, 0.261010, -0.214600, 
-0.171840, -0.213260, 0.129350, -0.268010, -0.062530, 0.086010, 
0.036530, -0.153130, -0.187220, -0.269540, -0.181950, -0.141840, 
0.052160, -0.107640, 0.037690, -0.071890, 0.183470, -0.334470, 
-0.131990, -0.007330, -0.030170, -0.364480, -0.090850, 0.018260, 
-0.186350, -0.249780, 0.053180, 0.054210, -0.339620, -0.703170, 
0.129190, 0.141560, -0.207660, -0.796480, -0.182030, 0.357630, 
-0.273440}
    },
    { /* unit 196 (Old: 196) */
      0.0, 0.125050, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.198420, 0.411340, -0.192940, -0.216770, -0.150880, 0.336390, 
-0.401090, 0.170820, -0.183970, 0.283560, -0.319410, 0.024490, 
-0.175180, 0.072580, -0.460210, 0.238240, -0.157380, 0.204440, 
-0.306790, 0.086920, -0.026980, 0.063150, -0.497900, 0.392860, 
-0.287560, -0.320050, -0.227370, 0.205510, -0.354170, -0.426220, 
-0.220290, -0.014690, -0.077990, -0.573380, -0.232390, -0.245990, 
0.001850, -0.712440, 0.010050, -0.278380, -0.067080, -0.586330, 
-0.098100, -0.231810, -0.221870, -0.709820, -0.026150, -0.210210, 
0.150350, -0.438190, -0.037910, -0.220760, -0.159710, -0.377970, 
0.021820, 0.223330, -0.085470, -0.429870, -0.104740, 0.285050, 
-0.065310, -0.554780, -0.004590, 0.238340, -0.063880, -0.566020, 
-0.251950, 0.616630, 0.003280, -0.291220, -0.081870, 0.393260, 
0.061590, -0.228080, -0.089290, 0.531980, 0.100710, -0.332670, 
0.017190, 0.450270, -0.169020, -0.240680, -0.185540, 0.291840, 
-0.021640, 0.059220, -0.227650, 0.470740, 0.015940, -0.042530, 
-0.285720, 0.085030, 0.428940, 0.062580, -0.202610, 0.096800, 
0.335430, -0.022940, -0.039710, -0.054300, -0.088810, -0.087560, 
-0.018140, -0.042440, -0.044790, 0.211550, -0.129690, -0.321130, 
-0.078140, 0.135020, 0.076930, -0.346660, 0.060530, 0.161570, 
0.085370, -0.023090, -0.078960, 0.316990, -0.162050, -0.130600, 
-0.037690, 0.525430, -0.097930, -0.166720, 0.012080, 0.310510, 
0.062410, -0.162450, -0.075200, 0.490160, 0.096550, -0.370150, 
0.014180, 0.184340, -0.171360, -0.405010, -0.010130, 0.000460, 
-0.174030, -0.455330, 0.278150, -0.201570, -0.078180, -0.493390}
    },
    { /* unit 197 (Old: 197) */
      0.0, -0.635440, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.535030, -0.314580, -0.134890, 0.070090, -0.339530, 
-0.017140, -0.325740, -0.097650, -0.371610, -0.027290, -0.278510, 
0.045210, -0.462670, -0.296710, -0.468670, -0.204100, -0.342730, 
-0.015850, -0.478350, -0.195990, -0.111730, -0.168090, -0.386180, 
-0.149110, -0.013890, -0.053230, -0.379400, -0.062550, 0.167300, 
-0.249460, 0.025490, -0.217620, -0.003100, -0.045090, -0.032140, 
-0.280610, -0.064490, -0.159730, -0.351000, -0.037990, -0.267860, 
0.096450, -0.330000, 0.067610, -0.090540, 0.165210, 0.116740, 
0.015540, -0.475040, -0.203650, -0.092660, 0.037390, -0.435340, 
-0.002370, 0.087030, 0.118760, -0.173710, 0.059340, 0.370650, 
-0.100560, -0.406430, -0.050200, 0.216750, -0.300520, 0.080430, 
-0.255360, 0.351390, -0.179350, 0.001200, -0.150470, 0.357980, 
-0.293810, 0.077080, -0.277570, 0.394410, -0.409670, -0.008070, 
-0.279520, 0.262900, -0.308210, 0.131370, -0.001620, 0.149420, 
-0.270500, 0.210810, -0.026460, -0.036600, -0.324510, 0.182120, 
-0.160110, -0.071390, 0.059580, -0.202250, -0.205690, -0.348380, 
0.103220, 0.112770, -0.087280, -0.198040, 0.149130, -0.200660, 
0.017780, -0.144590, 0.000260, 0.082940, -0.085380, -0.133290, 
0.131360, -0.180750, 0.150030, -0.374200, 0.288700, -0.194930, 
0.043120, -0.284160, 0.174010, 0.010210, 0.163520, -0.091640, 
0.336350, 0.031770, 0.091350, -0.389850, 0.247940, 0.171900, 
0.208410, -0.205180, 0.216920, 0.015210, -0.196370, -0.049950, 
0.402140, -0.006900, 0.144270, -0.396210, 0.605400, 0.064500, 
-0.123450, -0.449050, 0.543820, -0.294290, -0.155360, -0.434580, 
0.078580}
    },
    { /* unit 198 (Old: 198) */
      0.0, -0.749830, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.260910, 0.476880, 0.051460, -0.049470, 0.257360, 0.054640, 
-0.217580, -0.060030, 0.034080, 0.134290, -0.358000, -0.277820, 
0.025570, -0.158370, -0.471170, -0.150380, -0.352450, -0.099550, 
-0.372310, 0.063500, 0.109930, -0.219960, -0.252180, 0.238010, 
-0.055080, -0.186780, -0.057670, 0.277100, 0.270250, -0.246930, 
0.043500, 0.111850, 0.005660, 0.005390, -0.015050, 0.050740, 
-0.020610, 0.012420, 0.103640, -0.300060, 0.080110, 0.143910, 
0.024160, -0.018890, 0.122130, 0.166330, -0.125620, -0.308070, 
-0.034350, 0.211050, 0.125380, -0.133300, -0.001150, 0.014580, 
0.416910, -0.081450, -0.004990, -0.166090, 0.559200, -0.077730, 
-0.282470, -0.288990, 0.390190, 0.191670, -0.125980, -0.362330, 
0.589540, -0.052190, -0.184720, -0.739920, 0.509390, 0.053530, 
-0.132390, -0.884430, 0.311020, 0.061440, -0.214230, -0.832140, 
0.541460, -0.037070, -0.261340, -0.511120, 0.343190, -0.140840, 
-0.165590, -0.535550, -0.101180, -0.041520, 0.113940, -0.411230, 
-0.240860, 0.143440, -0.105010, -0.359650, -0.047240, 0.047030, 
0.290150, -0.016740, 0.023710, -0.101980, 0.178050, -0.307920, 
-0.128040, 0.089230, 0.311020, -0.026870, 0.120220, -0.000420, 
0.152020, 0.105990, -0.083210, 0.191430, -0.105240, 0.119780, 
0.317810, 0.116850, -0.103000, 0.007090, 0.070420, 0.194830, 
0.013480, -0.221550, 0.060880, 0.159350, -0.095230, -0.132420, 
0.193150, 0.056430, 0.186220, -0.058370, -0.019470, 0.185060, 
0.315980, -0.250640, -0.062650, 0.366900, 0.163510, 0.036610, 
-0.064450, 0.243550, 0.324630, -0.027970, -0.051620, 0.131380}
    },
    { /* unit 199 (Old: 199) */
      0.0, 0.023940, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.037570, -0.647870, -0.860650, -0.238480, -0.454350, 
-0.394820, -0.426630, -0.337880, 0.000600, -0.506350, -0.265280, 
-0.225950, -0.099280, -0.438170, -0.157820, -0.238830, 0.066610, 
-0.295650, -0.375460, 0.159100, -0.135150, -0.330560, -0.044370, 
0.092530, 0.056360, -0.491920, -0.384430, -0.099520, -0.227070, 
-0.848580, -0.132080, -0.207710, -0.296530, -0.737690, -0.318760, 
-0.269140, -0.088030, -0.771860, -0.140150, -0.368300, -0.210770, 
-1.085590, -0.129470, -0.334730, -0.323830, -1.116470, -0.269970, 
-0.480360, -0.273890, -1.046580, -0.514030, -0.462070, -0.303520, 
-0.916660, -0.504970, -0.111000, 0.034920, -0.971040, -0.703280, 
0.159310, 0.752780, -0.967800, -0.528910, 0.556150, 0.550210, 
-0.779470, -0.387200, 0.568450, 0.583810, -0.725900, -0.523590, 
0.702670, 0.703120, -0.782140, -0.476320, 0.963270, 0.352990, 
-0.426630, -0.374960, 0.922680, 0.422450, -0.501930, -0.095280, 
0.887310, 0.228120, -0.218690, -0.212050, 0.421220, 0.192340, 
-0.405210, -0.493660, 0.116510, 0.472730, -0.599320, -0.554710, 
-0.037710, -0.217510, -0.526840, -0.510550, -0.091200, -0.196150, 
-0.078640, -0.428420, -0.610730, -0.171480, -0.034690, -0.169930, 
-0.375260, -0.210240, -0.062320, 0.074950, -0.280820, 0.181130, 
0.332380, 0.170250, -0.289770, 0.306380, 0.315670, 0.442020, 
-0.294030, 0.171450, 0.468040, 0.257670, 0.031870, 0.253060, 
0.387020, 0.422540, 0.190690, 0.141370, 0.283850, 0.428610, 0.134850, 
0.218970, 0.238810, 0.104770, 0.008890, -0.247230, -0.192780, 
0.125270, -0.054150, -0.202860, -0.267970, 0.033520, -0.056010}
    },
    { /* unit 200 (Old: 200) */
      0.0, -0.032820, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.036420, 0.025260, 0.001690, 0.029270, -0.068350, -0.121860, 
0.019990, -0.120880, 0.120280, 0.082640, -0.165760, -0.114100, 
0.119860, -0.005090, 0.075990, 0.193850, -0.093000, 0.165520, 
-0.176140, 0.175270, 0.048450, 0.087730, -0.011020, 0.028480, 
-0.186600, 0.040410, 0.006250, 0.001560, 0.021380, -0.237640, 
0.000490, 0.066440, -0.171130, 0.104290, -0.200530, 0.080440, 
-0.043100, -0.101670, -0.034770, 0.033020, -0.182200, -0.257950, 
0.096460, -0.101490, -0.179370, 0.018830, -0.239990, -0.001580, 
0.145780, 0.069770, -0.122510, -0.160660, -0.020300, 0.036950, 
0.111470, -0.138700, 0.032890, -0.188230, -0.105540, 0.074390, 
-0.255070, 0.108080, 0.157730, -0.312420, -0.032070, -0.185260, 
0.125220, -0.191480, -0.312500, -0.223060, -0.252000, -0.171980, 
-0.118630, -0.018500, -0.272000, -0.175870, -0.265200, -0.137200, 
-0.106290, -0.411200, -0.209390, -0.074010, -0.303030, -0.063340, 
-0.196160, -0.189130, -0.032360, -0.037480, 0.036720, -0.194830, 
-0.027890, 0.067730, -0.013800, -0.077660, 0.079490, -0.150820, 
-0.264470, -0.093960, -0.177470, -0.222620, 0.072290, -0.036890, 
-0.306020, -0.150540, -0.295760, 0.043640, 0.052120, -0.173310, 
-0.269050, -0.080180, -0.258070, -0.143850, -0.281170, 0.096190, 
-0.249690, 0.082160, -0.052250, -0.225620, 0.104030, -0.022140, 
-0.275580, -0.086270, -0.069020, 0.119920, -0.262190, -0.038790, 
-0.161980, 0.011340, -0.168370, -0.162650, -0.004110, -0.180360, 
-0.217000, -0.196650, -0.083240, -0.154490, -0.080420, -0.013820, 
-0.014080, 0.143950, -0.282560, -0.136590, -0.168560, -0.213790}
    },
    { /* unit 201 (Old: 201) */
      0.0, -0.041630, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.172670, 0.074030, 0.058420, 0.052310, -0.274020, -0.120650, 
-0.239680, 0.043740, 0.028530, 0.105860, -0.259950, -0.190810, 
-0.248150, -0.019440, -0.017110, -0.038380, -0.115520, 0.040760, 
0.004110, 0.180130, -0.281550, -0.006980, 0.152400, 0.235410, 
-0.240220, 0.171570, 0.059900, -0.035260, -0.399930, -0.100430, 
0.136150, 0.042850, -0.042270, 0.023230, -0.193030, 0.178200, 
-0.086500, 0.152720, 0.087370, -0.021240, -0.048460, -0.058010, 
0.059280, 0.145700, 0.022010, 0.029570, 0.058470, 0.067180, 0.089980, 
-0.108710, -0.169320, 0.169270, -0.036150, -0.047290, 0.034830, 
0.167710, 0.109450, -0.024680, 0.103510, -0.176220, -0.169920, 
-0.299750, -0.028770, -0.277760, -0.282930, -0.299000, -0.078350, 
-0.296780, -0.234560, -0.119110, -0.092300, -0.527220, -0.184470, 
-0.118300, 0.131310, -0.379640, -0.136030, -0.122680, -0.050890, 
-0.706490, -0.012960, -0.309970, 0.002820, -0.687100, 0.113300, 
-0.275180, 0.051780, -0.407190, 0.255390, -0.285160, -0.366620, 
-0.379300, 0.290270, -0.119900, -0.335940, -0.431990, -0.092040, 
-0.317890, -0.276460, 0.026460, 0.233920, -0.103350, -0.268270, 
-0.139620, -0.081780, 0.170390, -0.323750, -0.192100, 0.073090, 
0.228260, 0.000660, -0.077490, -0.258370, -0.075510, -0.083510, 
-0.003810, -0.119710, 0.193130, 0.005280, -0.183720, -0.131790, 
-0.001230, -0.024400, -0.318540, -0.195400, 0.062980, 0.062480, 
-0.226880, -0.171060, 0.079890, 0.032560, -0.084560, -0.280320, 
-0.140030, 0.219160, 0.020820, -0.007570, -0.037240, -0.047740, 
-0.095850, -0.209490, 0.132280, -0.029470, -0.174330}
    },
    { /* unit 202 (Old: 202) */
      0.0, -1.163030, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.522220, 0.216740, -0.555770, 0.187640, -0.059690, 0.234220, 
-0.204350, 0.010840, -0.241900, 0.084100, 0.101020, -0.039200, 
0.083320, 0.339220, 0.438830, -0.231190, 0.114710, 0.117040, 
0.374580, -0.510140, 0.137560, 0.107220, 0.462150, -0.558170, 
-0.010250, 0.306900, -0.062600, -0.269490, -0.176950, -0.091220, 
-0.146180, -0.123240, -0.076400, 0.191070, -0.329710, 0.065810, 
0.234090, -0.015630, -0.472200, 0.029330, 0.182220, -0.090370, 
-0.545840, 0.236500, 0.538390, 0.080460, -0.505530, -0.051990, 
0.334700, -0.091260, -0.345290, 0.006550, 0.598010, -0.002060, 
-0.092460, -0.552720, 0.447970, -0.049300, 0.190800, -0.943280, 
0.502630, 0.049070, 0.365440, -1.107740, 0.417440, 0.113740, 
0.271000, -1.004400, 0.313340, 0.241630, -0.006100, -0.926240, 
0.678160, 0.664860, 0.037260, -0.714990, 0.707950, 0.705710, 
0.155450, -0.366000, 0.773820, 0.587810, 0.062520, -0.551010, 
0.227300, 0.358610, 0.233100, -0.497680, 0.006880, -0.017660, 
0.272450, -0.651790, -0.420380, -0.042390, 0.564100, -0.680210, 
-0.565160, 0.003780, 0.673330, -0.516030, -0.425700, 0.321530, 
0.302950, -0.346000, -0.562710, 0.065220, 0.489260, -0.587470, 
-0.469200, 0.031300, 0.393760, -0.548730, -0.158040, -0.154230, 
-0.009790, -0.488580, 0.228520, 0.098290, -0.199930, -0.550820, 
0.407300, 0.229500, -0.055160, -0.431290, 0.421550, 0.534400, 
-0.372720, -0.308080, 0.358050, 0.520450, -0.474290, -0.412820, 
0.446590, 0.356990, -0.341650, -0.081410, 0.515040, 0.448360, 
-0.358730, -0.131050, 0.353810, 0.433530, -0.294660, 0.376060}
    },
    { /* unit 203 (Old: 203) */
      0.0, -0.473100, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.019810, -0.044660, 0.029660, 0.131780, -0.097240, -0.335910, 
-0.091280, 0.091520, 0.185110, -0.289120, 0.136420, 0.058670, 
0.065620, -0.144880, 0.240680, 0.134620, 0.021600, -0.057240, 
-0.015320, -0.008240, -0.180700, -0.074610, -0.145180, -0.186430, 
-0.082600, -0.094880, -0.231100, -0.251480, 0.029640, 0.026990, 
-0.022570, -0.141680, -0.103810, -0.053730, 0.001580, 0.105390, 
0.105430, -0.427950, -0.243870, 0.284830, -0.132750, -0.082250, 
-0.007830, 0.205110, 0.048560, -0.351490, -0.248630, -0.060490, 
0.204380, -0.090420, -0.094340, 0.012170, 0.046930, 0.041350, 
-0.013580, 0.159300, 0.018590, 0.014580, 0.058000, -0.088050, 
-0.308190, 0.268800, -0.055950, -0.298910, -0.155550, 0.378710, 
-0.288270, -0.244930, -0.136620, 0.654520, -0.441460, -0.368760, 
-0.236260, 0.263470, -0.197210, -0.307840, -0.287490, 0.215340, 
-0.289500, -0.323030, -0.291860, 0.346320, 0.001330, -0.331970, 
-0.227080, 0.416950, -0.068850, 0.017360, -0.346000, 0.417260, 
0.163130, -0.036770, -0.321200, 0.098730, 0.226750, -0.069220, 
-0.326120, 0.385920, 0.395530, -0.152260, -0.459820, 0.303950, 
0.006750, -0.000450, -0.153990, -0.106620, 0.205450, 0.088200, 
-0.313590, -0.073490, -0.037240, -0.093920, -0.032010, -0.176090, 
-0.035970, 0.106160, -0.011580, -0.254400, -0.291390, 0.094440, 
-0.236330, -0.101540, -0.101460, -0.263980, -0.036330, -0.178970, 
-0.006190, -0.203350, 0.030540, 0.057580, -0.037170, -0.171970, 
-0.349270, -0.235570, -0.276660, -0.207250, -0.131220, -0.094800, 
-0.029990, -0.051700, -0.345180, -0.171920, -0.100200, 0.145240}
    },
    { /* unit 204 (Old: 204) */
      0.0, 0.720830, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.651500, -0.470020, -0.404660, -0.627080, -0.548960, 
-0.248870, -0.064000, -0.072200, -0.431620, -0.384040, 0.320680, 
0.123710, -0.142680, -0.111660, 0.673950, 0.163320, 0.279070, 
0.242120, 0.776400, 0.358230, 0.156330, 0.259180, 0.563420, 0.448510, 
0.324650, 0.118780, 0.468430, 0.401710, 0.110960, 0.148880, 0.459440, 
0.488570, 0.162110, 0.149850, 0.354200, 0.449090, 0.148790, 
-0.119520, 0.497420, 0.354310, 0.380630, -0.097510, 0.283100, 
0.652090, 0.085640, -0.094220, 0.202640, 0.693010, -0.054550, 
-0.309740, -0.089410, 0.705620, 0.068930, -0.254080, -0.151970, 
0.505240, 0.270320, 0.176700, -0.302890, 0.394380, 0.456410, 
0.340900, -0.518130, 0.783300, 0.743780, 0.576920, -0.470810, 
0.889530, 0.690070, 0.637810, -0.487260, 1.138160, 1.304930, 
0.773690, -0.338460, 1.299520, 1.393080, 0.832320, -0.697750, 
1.548120, 1.522100, 1.065110, -0.868500, 1.241840, 1.007740, 
0.950120, -0.942820, 0.972480, 0.593780, 0.819440, -0.802790, 
0.704270, 0.618930, 0.481790, -0.578500, 0.543410, 0.119300, 
0.472030, -0.565980, -0.011200, -0.372240, 0.534670, -0.175910, 
-0.157420, -0.209690, 0.221670, -0.195180, 0.042300, -0.333980, 
-0.124010, -0.131860, -0.252690, 0.073250, -0.260620, -0.290480, 
0.019240, -0.193610, -0.474220, -0.183810, -0.388950, -0.339070, 
-0.129340, -0.012730, -0.244970, -0.340620, -0.520370, -0.401770, 
-0.449340, -0.189630, -0.001220, -0.072810, -0.560250, -0.499480, 
-0.105130, -0.024230, -0.756430, -0.633370, -0.125640, -0.148430, 
-0.731740, -0.223960, 0.225320, 0.004160, -0.485860}
    },
    { /* unit 205 (Old: 205) */
      0.0, -0.924620, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-1.195510, -0.452310, -0.520740, 0.330340, -0.825560, 
-0.102850, -0.743310, 0.188430, -0.816590, -0.206120, -0.827840, 
-0.311290, -0.786290, -0.151370, -0.638000, -0.414120, -0.394240, 
-0.060360, -0.885800, -0.415660, -0.148090, 0.219760, -0.790650, 
-0.495580, 0.034320, 0.293670, -0.488590, -0.085620, -0.092240, 
0.246870, -0.235270, -0.081760, -0.006640, -0.124110, -0.387990, 
0.059120, -0.051110, -0.333730, -0.083040, -0.115400, -0.149380, 
-0.416080, -0.349590, -0.164490, -0.522780, -0.190590, -0.093870, 
-0.247360, -0.774590, -0.262280, 0.157700, -0.437430, -0.572360, 
-0.050210, 0.733230, -0.679930, -0.776470, -0.015050, 1.056090, 
-0.924790, -0.573010, 0.541710, 1.316330, -0.948950, -0.261340, 
0.866960, 1.344170, -0.735540, -0.069080, 0.935020, 0.744430, 
-0.961520, 0.123900, 0.834780, 0.791650, -0.592240, 0.121360, 
0.579540, 0.688950, -0.523790, 0.266430, 0.587830, 0.685520, 
-0.573430, 0.125280, 0.364580, 0.583580, -0.413330, 0.224550, 
0.455180, 0.335280, -0.137640, 0.020770, 0.069520, -0.016420, 
-0.129110, -0.158010, 0.055990, 0.129230, 0.069800, -0.179330, 
-0.229880, 0.032850, 0.075130, 0.047420, -0.056170, -0.145020, 
0.545620, -0.030230, -0.271260, -0.378270, 0.595260, 0.261520, 
-0.365350, -0.436380, 0.232270, 0.421160, -0.667680, -0.778780, 
0.188680, 0.577410, -0.569220, -0.993180, 0.326860, 0.646650, 
-0.454280, -0.767860, 0.178020, 0.350380, -0.407270, -0.475550, 
0.247450, 0.286640, -0.193230, -0.586550, 0.274610, 0.003140, 
0.107020, -1.117400, 0.031910, -0.103190, 0.106540, -1.225960, 
0.040000}
    },
    { /* unit 206 (Old: 206) */
      0.0, -0.354630, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.225390, 0.305920, 0.303880, -0.235840, 0.335650, 0.233910, 
0.032210, -0.300580, 0.009650, -0.068990, -0.150570, -0.211030, 
0.290300, 0.017100, -0.288710, -0.111520, 0.021050, -0.053440, 
-0.120500, -0.039630, 0.079210, -0.163930, -0.179310, -0.205030, 
0.114930, 0.011180, -0.008870, -0.092740, -0.138690, 0.099150, 
0.049390, 0.008090, 0.037000, -0.014240, 0.131140, -0.117260, 
0.147760, -0.049080, 0.018990, 0.047560, 0.062290, -0.249730, 
0.078460, -0.161540, 0.031720, 0.014570, 0.089420, -0.063020, 
-0.158760, 0.069270, -0.190620, -0.084060, 0.051820, -0.003630, 
0.020350, 0.056230, 0.057940, 0.122020, 0.010970, -0.034950, 
-0.219170, -0.056320, -0.256820, 0.212230, -0.193000, -0.271670, 
-0.277330, 0.113180, -0.461420, -0.283160, -0.382340, 0.118710, 
-0.396760, -0.391510, -0.315280, 0.467290, -0.737730, -0.415570, 
-0.235030, 0.459700, -0.758830, -0.524130, -0.053770, 0.186810, 
-0.494320, -0.169040, -0.147710, 0.098750, -0.332320, -0.130090, 
0.102810, 0.151950, -0.217170, -0.047360, 0.079280, 0.057510, 
-0.179260, 0.074430, 0.112070, 0.161250, 0.161780, -0.113090, 
0.139530, 0.190440, -0.027350, -0.257700, 0.055950, 0.028100, 
0.111640, 0.025960, -0.017790, -0.170260, -0.087730, -0.072940, 
0.259450, -0.211660, -0.112260, 0.218690, 0.309790, -0.178450, 
-0.220190, -0.131990, 0.461280, -0.238680, -0.065710, -0.050470, 
0.276840, -0.350610, -0.247310, -0.142860, 0.165590, -0.270410, 
0.006880, 0.009030, 0.206460, -0.042470, 0.097980, -0.001810, 
0.184780, -0.039160, 0.168140, 0.076740, 0.007960, -0.136240}
    },
    { /* unit 207 (Old: 207) */
      0.0, -0.170900, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.202230, -0.590280, -0.066270, 0.402570, -0.006980, 
-0.060750, -0.285460, 0.147800, 0.277880, -0.082270, -0.258020, 
0.148150, 0.124870, -0.308830, 0.110900, 0.040370, 0.476470, 
-0.298450, 0.065150, -0.296420, 0.394450, -0.205810, 0.119580, 
-0.246150, 0.388550, -0.169810, 0.038650, -0.484730, 0.198310, 
-0.459550, -0.006820, -0.153580, 0.181970, -0.221950, -0.295790, 
-0.421140, -0.239320, -0.422160, -0.155820, -0.105710, -0.168500, 
-0.389230, -0.354980, -0.292240, -0.083170, -0.304530, -0.138010, 
0.066270, -0.421370, -0.088620, -0.091780, -0.079100, -0.390360, 
-0.040170, 0.073470, -0.170790, -0.330230, -0.177000, 0.187950, 
-0.059960, -0.010570, 0.006900, 0.062870, -0.224500, -0.023880, 
0.244140, 0.033960, 0.014910, 0.155950, 0.444160, 0.252860, 
-0.201160, -0.059840, 0.460450, -0.127450, -0.212850, -0.089990, 
0.461130, 0.059730, -0.229900, 0.211610, 0.443760, 0.065960, 
0.005060, 0.175260, 0.265140, -0.071930, 0.009970, -0.070070, 
0.087020, 0.177560, 0.184600, -0.094230, -0.028280, 0.149490, 
0.121640, -0.138620, -0.118280, 0.182930, 0.260300, -0.494680, 
-0.001760, -0.137650, 0.253790, -0.648880, -0.134530, 0.151930, 
0.426930, -0.639390, -0.111170, -0.154280, 0.242760, -0.388220, 
-0.482320, 0.043060, 0.483860, -0.292960, -0.254780, 0.109240, 
0.466750, -0.546390, -0.130260, 0.057330, 0.308490, -0.432540, 
-0.284120, -0.005330, 0.478410, -0.017510, -0.270090, 0.125770, 
0.488030, -0.034920, -0.006760, -0.028200, 0.387230, -0.357800, 
0.012500, 0.033220, 0.361180, -0.201220, 0.037890, 0.031450, 0.288420}
    },
    { /* unit 208 (Old: 208) */
      0.0, -0.066710, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.399640, 1.263780, -0.171530, -0.394580, -0.267650, 
0.945630, -0.682530, -0.472540, -0.160110, 0.708570, -0.838530, 
-0.582160, -0.613040, 0.146240, -0.829220, -0.644350, -0.461380, 
-0.059370, -0.667900, -0.253710, -0.778460, -0.195770, 0.123250, 
-0.011920, -0.812300, -0.250380, 0.156730, 0.129780, -0.811720, 
-0.042410, 0.525220, 0.174380, -0.459360, -0.258070, 0.429010, 
0.239310, -0.683310, 0.003670, 0.261700, 0.157940, -0.871540, 
0.142490, 0.377860, -0.131410, -0.806570, 0.136520, 0.299570, 
-0.208650, -0.778080, -0.165130, 0.488400, -0.191610, -0.690720, 
0.124360, 0.530630, -0.221160, -1.221620, -0.227040, 1.045350, 
0.122380, -1.057140, -0.155660, 0.829940, 0.259620, -0.532520, 
-0.524510, 1.092930, 0.074390, -0.378190, -0.895790, 1.144150, 
0.007710, -0.247640, -1.317030, 0.868810, 0.262710, -0.442160, 
-1.404670, 0.887870, -0.048300, -0.346960, -1.355710, 0.462160, 
0.218350, -0.282520, -1.020670, -0.214810, 0.560460, -0.196410, 
-0.328830, -0.265690, 0.725070, -0.097700, 0.264000, -0.069940, 
0.798200, 0.020090, 0.346720, -0.203170, 0.611940, -0.119790, 
0.352770, 0.227120, 0.644630, 0.162010, 0.380640, 0.137520, 0.055630, 
0.097320, 0.530110, 0.135220, -0.029300, -0.163410, 0.463870, 
-0.091700, -0.090890, -0.271700, 0.562580, 0.194570, -0.481550, 
-0.025090, -0.037430, 0.061260, -0.713950, 0.120480, 0.108620, 
0.526900, -0.570910, 0.184360, -0.297580, 0.350580, -0.130420, 
0.638420, -0.574140, 0.456280, 0.233940, 0.479480, -0.481030, 
0.171280, 0.334440, 0.420660, -0.357670, -0.455060, 0.489490}
    },
    { /* unit 209 (Old: 209) */
      0.0, -0.068080, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.050270, 0.101330, -0.219620, 0.003580, 0.035650, 0.011820, 
-0.124700, 0.367800, 0.145870, -0.105130, 0.073300, 0.212180, 
0.008510, -0.227030, -0.160330, 0.290950, 0.001700, -0.395800, 
-0.027610, -0.005060, -0.156710, -0.358700, -0.223710, 0.207080, 
-0.477960, -0.246220, -0.364120, 0.149520, -0.342110, -0.485750, 
-0.150180, 0.062110, -0.223570, -0.434840, -0.252170, -0.171350, 
-0.180190, -0.201650, -0.320480, -0.229150, -0.091940, -0.208900, 
-0.174110, -0.201610, 0.185250, -0.113930, -0.243340, -0.026030, 
-0.005500, -0.203220, 0.037620, -0.196270, -0.157500, 0.224880, 
0.049560, -0.283710, 0.159880, 0.252980, -0.075500, 0.023210, 
-0.006090, 0.398050, -0.106220, -0.080410, -0.322320, 0.534130, 
-0.317720, -0.364370, -0.298160, 0.433580, -0.017530, -0.118880, 
-0.224700, 0.540330, -0.118810, -0.222370, -0.409530, 0.510160, 
0.014460, -0.047590, -0.386040, 0.214450, 0.250330, -0.128460, 
-0.197520, 0.408380, 0.436430, -0.010460, -0.150100, 0.296920, 
0.431180, -0.175000, -0.079050, 0.401300, 0.342930, -0.170480, 
-0.256140, 0.374170, 0.131280, -0.082260, -0.300030, 0.459020, 
0.341190, -0.051350, -0.119990, 0.340890, -0.028900, -0.005300, 
-0.180390, -0.102460, 0.066680, 0.027590, -0.086150, -0.030740, 
-0.323680, 0.185020, -0.199300, -0.146860, -0.413020, 0.213450, 
0.038370, -0.435920, -0.484130, -0.125140, -0.376390, -0.371320, 
-0.235830, -0.043510, -0.272320, -0.269340, -0.337550, 0.101020, 
-0.274350, -0.136590, -0.072080, -0.281020, -0.322130, -0.313100, 
0.042350, 0.019720, 0.010280, -0.205700, 0.222430, -0.159700}
    },
    { /* unit 210 (Old: 210) */
      0.0, -0.584830, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.020600, -0.026830, 0.136020, -0.100600, 0.105750, 0.034290, 
-0.052650, 0.028100, -0.066960, 0.060900, 0.097860, -0.132300, 
-0.108520, -0.026210, -0.098560, 0.060140, -0.084570, -0.047210, 
0.003650, 0.013150, -0.220810, -0.241030, 0.131430, -0.097370, 
-0.085280, -0.083740, -0.123390, -0.145100, -0.208500, -0.209330, 
0.017280, -0.100180, -0.193610, -0.058180, -0.118330, 0.132440, 
0.067070, -0.048300, -0.117030, -0.259240, -0.128520, -0.044320, 
0.130540, -0.226270, -0.250540, 0.070790, 0.088880, -0.211610, 
0.107030, -0.003860, 0.201470, -0.194360, 0.036980, -0.164450, 
0.071400, 0.067020, -0.130970, 0.134810, -0.022880, 0.047980, 
-0.161120, -0.126080, 0.231400, 0.017170, -0.016010, 0.070690, 
-0.100840, 0.072810, 0.045430, -0.044330, 0.070070, 0.024790, 
-0.095590, -0.024010, -0.027440, 0.051450, -0.086970, -0.038300, 
-0.170950, 0.129590, -0.154360, -0.030860, 0.078910, 0.159080, 
0.031180, -0.220900, 0.025480, -0.139130, 0.046610, 0.029800, 
0.135610, 0.034140, -0.171140, -0.100930, 0.024400, -0.161710, 
0.108910, -0.193210, 0.155370, 0.073190, -0.164090, -0.074430, 
-0.024340, -0.202460, 0.034420, -0.268730, 0.149980, -0.006480, 
-0.084560, 0.087330, -0.011550, -0.205890, -0.093030, 0.100320, 
0.104060, -0.133160, -0.212440, 0.136280, 0.055380, -0.181500, 
0.033350, 0.009830, 0.115310, -0.010240, -0.023520, -0.015220, 
0.027490, -0.158640, -0.113380, 0.124230, 0.097790, -0.059240, 
-0.025530, 0.077700, 0.120280, 0.064830, 0.028380, 0.086680, 
-0.203980, 0.039550, 0.086030, -0.215540, -0.196220, 0.059740}
    },
    { /* unit 211 (Old: 211) */
      0.0, -0.094920, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.261840, -0.516180, 0.050900, 0.026290, -0.216040, 
-0.240150, 0.115640, 0.194090, 0.194810, -0.249280, 0.339280, 
0.022710, 0.147770, -0.263470, 0.427650, -0.003780, 0.070770, 
-0.128360, 0.330350, -0.078750, -0.002830, -0.125780, 0.107180, 
0.020400, 0.024480, -0.307960, 0.155360, 0.023230, -0.132360, 
-0.432440, -0.020460, -0.223510, -0.191570, -0.418190, 0.070360, 
-0.228670, 0.152080, -0.553310, 0.104760, -0.112880, 0.101540, 
-0.520440, 0.289320, -0.177030, 0.109030, -0.169670, 0.077280, 
0.119050, -0.109720, -0.224830, 0.015890, 0.022600, -0.095230, 
-0.136200, -0.088620, 0.034910, -0.178170, -0.275520, 0.166790, 
-0.266430, -0.074470, -0.337870, -0.186410, -0.437350, -0.177610, 
-0.248380, 0.056560, -0.405060, -0.163260, 0.053660, 0.015010, 
-0.319890, -0.183670, 0.117750, 0.006060, -0.219940, -0.153010, 
0.042600, -0.022830, -0.041670, 0.155840, 0.177510, -0.148210, 
-0.317270, 0.110810, 0.014440, -0.072710, -0.210600, -0.224180, 
-0.093970, 0.255680, -0.347490, -0.223500, 0.127760, 0.307770, 
-0.320130, -0.129920, 0.082690, 0.196820, -0.047750, -0.015470, 
0.081890, 0.267210, -0.315060, -0.015030, -0.105570, 0.064490, 
-0.115570, 0.039640, -0.197300, 0.204790, 0.071210, 0.118030, 
-0.200170, 0.068130, -0.090450, -0.257390, 0.049440, 0.026450, 
0.139920, -0.349920, -0.165450, 0.104240, -0.017460, -0.165400, 
-0.242710, 0.039950, 0.027200, -0.436450, -0.184300, 0.049550, 
0.066440, -0.167640, 0.135830, 0.183880, 0.118210, -0.365530, 
0.032940, 0.156320, 0.106040, -0.498030, -0.058230, 0.172700, 
0.192290}
    },
    { /* unit 212 (Old: 212) */
      0.0, -0.612570, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.125520, -0.443200, -0.110940, 0.229110, -0.257370, 
-0.270830, 0.044000, -0.035100, -0.077970, -0.078790, -0.132500, 
0.031020, 0.211140, -0.050310, 0.293430, 0.111440, -0.002750, 
0.115390, -0.068530, 0.114980, 0.199230, 0.245540, 0.154280, 
0.125590, 0.018470, 0.081850, -0.092930, -0.233380, 0.015670, 
0.078610, -0.028240, -0.020590, -0.111460, 0.183650, 0.031120, 
-0.034010, 0.139970, 0.119390, -0.207570, 0.098530, -0.172970, 
-0.002950, -0.222700, 0.264740, -0.167710, 0.062850, 0.114530, 
0.194200, 0.118310, -0.201630, 0.140450, -0.016360, -0.086010, 
-0.111250, 0.053090, -0.013700, 0.092930, -0.145950, 0.150020, 
-0.060060, 0.185120, -0.189700, -0.113180, -0.197550, 0.343710, 
0.025000, 0.084480, -0.352130, 0.351470, 0.266810, 0.114640, 
-0.016270, 0.323020, 0.356560, 0.115300, -0.195430, 0.411020, 
0.199550, 0.108430, -0.103140, 0.459200, 0.400920, 0.142870, 
0.020910, 0.199670, 0.151440, -0.069090, -0.027830, 0.099250, 
0.378270, 0.178600, -0.075860, 0.097880, 0.047040, -0.034120, 
-0.280340, -0.071040, 0.056110, 0.111300, -0.286220, -0.083610, 
-0.171300, 0.305470, -0.199000, -0.022090, 0.052630, 0.143960, 
-0.129230, -0.243380, 0.012790, 0.122190, -0.223240, 0.007670, 
-0.313340, -0.066170, 0.084730, -0.160250, -0.356470, 0.134240, 
-0.000990, -0.146060, -0.062790, -0.127800, -0.086060, -0.119010, 
0.143890, 0.013210, -0.021940, -0.340780, -0.117410, 0.133250, 
0.088870, -0.181190, -0.059350, -0.169140, 0.094820, -0.329300, 
0.040900, -0.039100, -0.138890, -0.237460, 0.070070, -0.007670, 
0.041000}
    },
    { /* unit 213 (Old: 213) */
      0.0, 0.071140, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.005960, -0.120600, 0.213040, -0.161410, 0.140350, 0.174290, 
0.204120, 0.008430, 0.261010, -0.153170, 0.030530, 0.013030, 
-0.021800, -0.070360, -0.176510, 0.100430, 0.076010, -0.030600, 
-0.123450, 0.126900, 0.261240, -0.028530, -0.076640, -0.047200, 
0.170570, 0.137930, 0.040500, -0.211380, 0.005540, 0.109980, 
-0.104110, 0.062720, 0.079240, -0.059210, -0.201420, -0.238270, 
0.101010, 0.067120, -0.010610, -0.123600, -0.125990, -0.075870, 
0.112970, -0.073630, -0.115040, 0.169440, 0.030360, -0.204580, 
-0.087000, 0.182750, -0.082560, -0.115930, -0.233440, 0.305620, 
-0.124980, 0.218370, -0.276110, 0.099880, -0.097920, 0.216110, 
-0.074600, 0.018230, -0.214410, 0.266150, -0.398170, -0.185900, 
-0.191750, 0.086610, -0.299340, 0.003420, -0.236520, 0.139870, 
-0.577770, -0.107380, -0.111120, 0.009870, -0.617520, -0.126590, 
0.021490, -0.032970, -0.438760, -0.425760, -0.341960, 0.143270, 
-0.304200, -0.352110, -0.123460, 0.203820, 0.027960, -0.419420, 
-0.379860, 0.220990, -0.158150, -0.350200, -0.290850, -0.102090, 
-0.090490, -0.302180, 0.010260, 0.136260, 0.163480, -0.211380, 
-0.005960, 0.039560, 0.216640, -0.096750, -0.002710, -0.142980, 
0.306880, -0.239680, 0.120910, 0.037260, 0.255880, 0.012570, 
-0.063110, -0.299680, 0.108150, -0.002670, -0.184260, -0.082760, 
-0.016050, 0.086600, 0.017340, -0.311050, -0.041920, 0.010090, 
0.129220, -0.019020, 0.104720, -0.261520, 0.007910, -0.081400, 
0.089580, -0.225170, -0.080160, 0.048160, 0.211770, -0.082170, 
0.160710, 0.091150, 0.174200, -0.037320, 0.098770, -0.084430}
    },
    { /* unit 214 (Old: 214) */
      0.0, -0.213080, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.205190, -0.054110, -0.111170, 0.152170, 0.408920, 0.013310, 
-0.209990, 0.067330, 0.136000, 0.004620, -0.180080, -0.051840, 
0.077790, 0.200650, -0.179350, -0.150160, -0.070310, 0.042340, 
0.049790, 0.174370, -0.357270, 0.110640, -0.063130, -0.059170, 
-0.227000, 0.034360, 0.042510, 0.054370, -0.499610, -0.081810, 
0.201190, 0.330460, -0.349890, 0.138480, -0.105880, 0.016370, 
-0.443870, 0.076140, -0.265860, 0.303010, -0.182650, 0.175310, 
-0.075860, -0.048330, -0.164750, 0.130890, -0.149790, 0.031800, 
0.002140, 0.227030, -0.072740, 0.025710, -0.130710, -0.128360, 
0.189880, 0.227260, -0.414420, 0.104740, 0.043660, 0.066120, 
-0.190470, -0.110910, 0.245820, -0.079830, -0.621060, -0.290460, 
-0.029140, -0.404230, -0.337250, -0.297890, -0.035240, -0.790750, 
-0.398290, -0.296650, 0.136470, -0.975490, -0.595540, -0.368180, 
0.238430, -1.274550, -0.705130, -0.288100, -0.109540, -1.058530, 
-0.220170, -0.350120, -0.150680, -0.903340, -0.247240, -0.336260, 
-0.114790, -0.827650, 0.068270, -0.051980, -0.344270, -0.635410, 
0.062380, -0.060040, -0.056460, -0.440360, 0.090610, -0.014220, 
-0.220140, -0.417650, 0.091460, 0.153740, -0.267550, -0.308440, 
0.171730, 0.074820, -0.180210, -0.614970, 0.100740, 0.365880, 
-0.004990, -0.685960, 0.070060, 0.279260, -0.056130, -0.721440, 
-0.079900, 0.375660, 0.248310, -0.850850, -0.204350, 0.014210, 
0.171080, -0.513070, 0.220700, 0.113760, 0.158850, -0.589220, 
0.324330, -0.257090, 0.001960, -0.434550, 0.298370, -0.188100, 
0.049540, -0.276390, 0.016140, -0.080350, 0.030310, -0.031640}
    },
    { /* unit 215 (Old: 215) */
      0.0, -0.259910, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.415150, 0.621070, 0.134590, -0.094720, 0.079070, 0.155400, 
0.025650, -0.178090, -0.227750, 0.169800, 0.171180, 0.099290, 
-0.330420, 0.082020, 0.067940, -0.025240, -0.287910, -0.332990, 
0.104060, 0.112690, -0.127450, -0.106660, -0.037490, -0.014960, 
-0.265110, -0.323730, -0.096260, 0.146530, -0.096910, -0.271200, 
-0.253580, -0.112870, -0.111080, 0.021260, -0.041480, -0.191070, 
-0.230030, 0.027050, 0.036650, -0.034780, -0.205590, 0.036540, 
0.147570, 0.007060, -0.097290, 0.084640, 0.086330, -0.193920, 
0.141950, 0.035960, -0.084700, 0.131640, -0.058110, 0.163050, 
-0.041710, 0.208810, 0.099050, 0.048870, 0.000680, 0.267470, 
-0.162610, -0.305000, -0.484230, 0.301730, -0.098130, -0.395030, 
-0.487140, 0.260230, 0.066930, -0.200990, -0.528360, 0.245030, 
-0.117160, -0.355390, -0.488820, 0.316680, -0.302080, -0.362500, 
-0.548300, 0.277290, -0.160650, -0.147290, -0.199850, 0.199470, 
-0.085580, -0.096520, -0.111940, -0.026830, 0.107880, -0.128360, 
-0.040010, 0.183140, 0.187380, 0.067490, -0.094610, -0.151090, 
-0.153310, -0.016750, -0.095690, -0.167290, -0.210020, 0.183120, 
-0.201620, 0.107800, -0.047320, -0.047310, -0.102240, -0.169870, 
-0.165940, -0.018350, -0.194650, 0.002640, -0.050720, 0.078400, 
-0.394150, -0.085240, 0.050090, 0.145330, -0.508230, 0.236090, 
0.047590, 0.335260, -0.485770, 0.240430, -0.079840, -0.051590, 
-0.629630, 0.111080, 0.058960, 0.290570, -0.508310, -0.209120, 
-0.247340, -0.029090, -0.322730, -0.148290, 0.091180, 0.193400, 
-0.353130, -0.388420, -0.011460, -0.038800, -0.062460, -0.371090}
    },
    { /* unit 216 (Old: 216) */
      0.0, -0.265240, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.210370, 0.232530, 0.351030, -0.009410, -0.173340, 0.317570, 
0.055580, -0.154990, -0.039530, 0.084310, -0.182080, 0.099150, 
-0.035940, -0.203860, -0.103990, 0.082040, -0.074680, -0.258240, 
0.117220, 0.151540, 0.268940, -0.235340, 0.053720, 0.415470, 
0.153610, -0.274670, 0.064640, 0.154090, 0.222670, -0.134420, 
0.059310, 0.186530, 0.080390, -0.272910, 0.111260, -0.458800, 
-0.310130, -0.012910, -0.060100, -0.552380, -0.375480, 0.023270, 
-0.078840, -0.665180, -0.443400, -0.011020, 0.097050, -0.287880, 
-0.170700, -0.086950, -0.260070, -0.342280, -0.261670, -0.189200, 
-0.040440, -0.173220, -0.083620, -0.437880, 0.254200, 0.100320, 
-0.003140, -0.382740, 0.019540, 0.344380, 0.093980, -0.609520, 
0.190830, 0.444430, 0.313990, -0.482340, 0.359690, 0.596910, 
0.202590, -0.268360, 0.221290, 0.443270, 0.192990, -0.253560, 
0.100450, 0.357870, -0.005030, -0.058250, -0.225070, 0.555150, 
0.130650, -0.062790, -0.082470, 0.549740, 0.292570, 0.037620, 
-0.203110, 0.495570, 0.382070, 0.098540, -0.215300, 0.252920, 
0.511470, -0.165530, -0.227600, 0.244050, 0.494290, 0.190920, 
0.021130, 0.170790, 0.397160, 0.171570, 0.101750, -0.177840, 
0.213850, 0.001350, 0.242160, -0.315650, -0.351990, 0.223440, 
-0.039230, -0.253660, -0.160500, 0.197990, 0.050550, -0.391380, 
-0.558560, 0.249780, 0.117270, -0.089930, -0.244490, 0.051080, 
0.367260, -0.055580, -0.082870, 0.220190, 0.162610, 0.138590, 
0.031720, 0.158500, -0.129450, -0.196500, 0.121260, 0.024410, 
-0.038060, -0.033880, 0.237100, 0.245690, -0.031670, -0.319220}
    },
    { /* unit 217 (Old: 217) */
      0.0, 0.102300, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {-0.218980, -0.131040, -0.092820, 0.240330, -0.217280, 
0.131790, -0.263290, -0.141900, -0.210820, 0.141340, -0.057500, 
0.027380, -0.274550, 0.240210, 0.011820, -0.019450, -0.090180, 
0.278920, -0.198200, 0.015460, -0.189760, 0.353110, -0.098490, 
0.050000, -0.324830, 0.009440, 0.015360, 0.080910, -0.303870, 
-0.007670, -0.062310, 0.161670, -0.066300, -0.204440, -0.199250, 
0.135410, 0.086250, -0.025290, -0.056790, 0.080710, 0.167810, 
-0.130830, -0.351410, -0.019040, 0.240850, -0.191460, -0.330510, 
0.261650, 0.122720, -0.163540, -0.059130, -0.003990, 0.195660, 
0.006980, -0.000330, -0.149410, -0.246780, -0.181830, -0.053820, 
-0.238130, -0.151790, -0.286040, 0.118780, -0.257880, -0.389800, 
-0.128880, -0.029650, -0.722060, -0.095480, -0.072510, 0.047430, 
-0.636940, -0.247020, -0.097880, 0.020760, -0.742930, -0.097380, 
-0.202620, -0.139330, -0.947730, -0.273810, -0.155700, -0.264780, 
-0.653520, -0.102570, -0.363710, -0.142910, -0.688970, 0.001360, 
-0.360510, -0.040020, -0.659500, 0.188070, -0.160230, -0.312480, 
-0.621150, -0.027690, 0.103670, -0.280520, -0.528250, -0.125490, 
0.309570, -0.135460, -0.159370, 0.051610, 0.583280, -0.051180, 
-0.146850, -0.000450, 0.314180, -0.220530, -0.229360, -0.279320, 
0.216880, 0.071810, -0.577270, -0.053010, 0.104880, -0.134520, 
-0.364730, -0.351320, 0.183380, 0.118930, -0.629890, -0.147100, 
-0.025100, -0.135110, -0.545380, -0.209360, -0.005360, -0.003160, 
-0.217950, -0.211550, -0.115200, -0.061060, -0.413340, 0.027150, 
-0.071660, -0.014330, -0.197650, -0.215940, -0.147480, -0.117650, 
-0.172240}
    },
    { /* unit 218 (Old: 218) */
      0.0, -0.325550, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.127650, 0.014620, 0.173060, -0.102460, -0.175250, 0.150340, 
-0.052680, -0.115980, -0.119600, -0.181050, -0.138680, 0.119310, 
-0.120830, 0.112520, -0.170930, -0.035200, -0.201210, -0.037520, 
-0.025630, -0.021330, -0.049540, -0.085120, -0.323150, -0.119980, 
-0.175790, 0.096790, -0.230040, -0.173150, -0.146810, 0.099970, 
-0.221740, -0.127550, -0.107350, -0.125840, 0.031700, 0.074970, 
-0.082360, 0.135060, -0.202130, -0.215720, -0.117540, 0.210080, 
0.081020, -0.092740, 0.013400, 0.235920, 0.061680, 0.203530, 
0.008250, 0.129670, -0.038230, -0.129130, -0.138400, 0.222590, 
0.068140, -0.086100, -0.007630, 0.080890, 0.062810, 0.055470, 
0.088880, 0.057680, -0.141760, -0.178370, -0.040210, -0.209750, 
0.125240, -0.026090, -0.195690, 0.120800, -0.086710, 0.070630, 
-0.230360, 0.097290, 0.034400, -0.164490, -0.198740, 0.130970, 
-0.135890, -0.166160, -0.124150, 0.104520, -0.157510, -0.154050, 
-0.042360, -0.221580, 0.189310, -0.043850, -0.011640, 0.082770, 
-0.140630, 0.015590, 0.114410, 0.054280, 0.092720, 0.062530, 
0.094470, 0.122270, 0.104660, 0.054250, 0.099110, 0.128610, 
-0.067870, 0.150030, -0.264720, 0.191180, -0.131930, -0.195810, 
-0.044450, -0.008460, -0.228300, 0.134730, -0.138380, 0.150890, 
-0.356380, -0.170460, -0.302100, -0.083300, -0.002130, 0.027830, 
-0.233990, -0.161180, -0.090850, -0.200130, -0.316870, -0.251280, 
-0.271100, -0.193930, -0.340400, -0.205790, -0.427830, -0.216610, 
-0.050190, -0.161550, -0.282450, -0.127160, -0.098630, 0.020280, 
-0.001130, -0.107940, 0.014320, 0.086490, -0.176220, -0.120830}
    },
    { /* unit 219 (Old: 219) */
      0.0, -0.900020, 144,
      {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 
6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 
12, Units + 13, Units + 14, Units + 15, Units + 16, Units + 17, Units 
+ 18, Units + 19, Units + 20, Units + 21, Units + 22, Units + 23, 
Units + 24, Units + 25, Units + 26, Units + 27, Units + 28, Units + 
29, Units + 30, Units + 31, Units + 32, Units + 33, Units + 34, Units 
+ 35, Units + 36, Units + 37, Units + 38, Units + 39, Units + 40, 
Units + 41, Units + 42, Units + 43, Units + 44, Units + 45, Units + 
46, Units + 47, Units + 48, Units + 49, Units + 50, Units + 51, Units 
+ 52, Units + 53, Units + 54, Units + 55, Units + 56, Units + 57, 
Units + 58, Units + 59, Units + 60, Units + 61, Units + 62, Units + 
63, Units + 64, Units + 65, Units + 66, Units + 67, Units + 68, Units 
+ 69, Units + 70, Units + 71, Units + 72, Units + 73, Units + 74, 
Units + 75, Units + 76, Units + 77, Units + 78, Units + 79, Units + 
80, Units + 81, Units + 82, Units + 83, Units + 84, Units + 85, Units 
+ 86, Units + 87, Units + 88, Units + 89, Units + 90, Units + 91, 
Units + 92, Units + 93, Units + 94, Units + 95, Units + 96, Units + 
97, Units + 98, Units + 99, Units + 100, Units + 101, Units + 102, 
Units + 103, Units + 104, Units + 105, Units + 106, Units + 107, 
Units + 108, Units + 109, Units + 110, Units + 111, Units + 112, 
Units + 113, Units + 114, Units + 115, Units + 116, Units + 117, 
Units + 118, Units + 119, Units + 120, Units + 121, Units + 122, 
Units + 123, Units + 124, Units + 125, Units + 126, Units + 127, 
Units + 128, Units + 129, Units + 130, Units + 131, Units + 132, 
Units + 133, Units + 134, Units + 135, Units + 136, Units + 137, 
Units + 138, Units + 139, Units + 140, Units + 141, Units + 142, 
Units + 143, Units + 144},
      {0.150850, -0.466000, 0.046440, 0.120870, 0.208400, -0.445190, 
0.200170, 0.020040, 0.232050, -0.536550, 0.374090, 0.253580, 
0.061860, -0.450900, 0.069570, -0.014920, -0.085620, -0.669930, 
0.080890, 0.088230, 0.069230, -0.761340, -0.160030, -0.252330, 
-0.027260, -0.521420, -0.210580, 0.066830, 0.319300, -0.710380, 
-0.207130, 0.197110, 0.433420, -0.627270, -0.213310, 0.138490, 
0.153290, -0.331830, 0.119560, 0.126340, 0.217240, -0.146550, 
0.030140, 0.268580, 0.128290, -0.224580, 0.087630, 0.211220, 
0.182610, 0.059170, -0.226710, 0.020680, 0.289990, 0.402510, 
-0.267010, 0.186590, -0.119100, 0.698840, -0.028930, -0.244390, 
-0.321800, 0.555690, -0.138800, -0.289230, -0.315560, 0.544930, 
-0.473740, -0.120270, -0.638010, 0.747580, -0.612300, -0.505500, 
-0.707250, 0.411170, -0.644950, -0.119350, -0.689130, 0.544080, 
-0.517460, -0.481950, -0.429210, 0.374640, 0.007010, -0.238190, 
-0.531830, 0.509680, 0.222770, -0.309590, -0.411830, 0.455630, 
0.155380, -0.086710, -0.689030, 0.561030, 0.333810, -0.159140, 
-0.516520, 0.580970, 0.099320, 0.071990, -0.583660, 0.350050, 
0.332130, -0.117100, -0.260000, 0.180260, 0.014160, -0.031060, 
-0.237990, 0.107670, -0.119740, -0.041060, -0.058200, -0.092810, 
-0.259040, -0.014090, -0.253560, -0.084720, -0.159270, -0.046370, 
0.092680, -0.196160, -0.558880, -0.298570, 0.059720, -0.166430, 
-0.179510, -0.253730, 0.056070, -0.335900, -0.225830, -0.329740, 
-0.346750, -0.150220, -0.189710, -0.154960, -0.388890, 0.169020, 
0.034810, 0.186410, -0.139290, -0.073410, 0.022820, 0.183440}
    },
    { /* unit 220 (Old: 220) */
      0.0, 0.206130, 75,
      {Units + 145, Units + 146, Units + 147, Units + 148, Units + 
149, Units + 150, Units + 151, Units + 152, Units + 153, Units + 154, 
Units + 155, Units + 156, Units + 157, Units + 158, Units + 159, 
Units + 160, Units + 161, Units + 162, Units + 163, Units + 164, 
Units + 165, Units + 166, Units + 167, Units + 168, Units + 169, 
Units + 170, Units + 171, Units + 172, Units + 173, Units + 174, 
Units + 175, Units + 176, Units + 177, Units + 178, Units + 179, 
Units + 180, Units + 181, Units + 182, Units + 183, Units + 184, 
Units + 185, Units + 186, Units + 187, Units + 188, Units + 189, 
Units + 190, Units + 191, Units + 192, Units + 193, Units + 194, 
Units + 195, Units + 196, Units + 197, Units + 198, Units + 199, 
Units + 200, Units + 201, Units + 202, Units + 203, Units + 204, 
Units + 205, Units + 206, Units + 207, Units + 208, Units + 209, 
Units + 210, Units + 211, Units + 212, Units + 213, Units + 214, 
Units + 215, Units + 216, Units + 217, Units + 218, Units + 219},
      {0.097170, -0.747730, 0.033560, 1.106090, 0.222110, 1.252980, 
-0.350270, 0.207070, 0.069240, -2.077320, 0.224230, 3.248130, 
-0.055120, -1.561320, -2.481230, -1.915180, -0.048150, 0.107220, 
-0.452680, 0.976190, -0.723540, -0.290320, -0.309650, 1.490060, 
1.961320, 1.454210, -4.160790, 1.894920, -0.973060, 0.641240, 
-2.301580, -1.144120, -0.977800, 2.791080, 0.601950, -2.508340, 
0.685300, 0.019240, -1.959790, 0.006080, -2.529160, -0.205420, 
1.355900, 0.333990, 0.392390, 1.149200, -3.862110, 0.443230, 
-0.816740, 1.426590, -0.069200, 0.561990, -0.695320, -1.158190, 
1.439960, 0.879510, 0.936870, 0.820590, -0.210080, 2.840200, 
-2.810350, -0.905150, 0.873520, -3.528950, 0.960640, 0.088580, 
0.932970, -0.518210, 0.366860, 0.658960, 0.232610, -0.000020, 
1.630140, 0.748500, -1.345950}
    },
    { /* unit 221 (Old: 221) */
      0.0, -1.027170, 75,
      {Units + 145, Units + 146, Units + 147, Units + 148, Units + 
149, Units + 150, Units + 151, Units + 152, Units + 153, Units + 154, 
Units + 155, Units + 156, Units + 157, Units + 158, Units + 159, 
Units + 160, Units + 161, Units + 162, Units + 163, Units + 164, 
Units + 165, Units + 166, Units + 167, Units + 168, Units + 169, 
Units + 170, Units + 171, Units + 172, Units + 173, Units + 174, 
Units + 175, Units + 176, Units + 177, Units + 178, Units + 179, 
Units + 180, Units + 181, Units + 182, Units + 183, Units + 184, 
Units + 185, Units + 186, Units + 187, Units + 188, Units + 189, 
Units + 190, Units + 191, Units + 192, Units + 193, Units + 194, 
Units + 195, Units + 196, Units + 197, Units + 198, Units + 199, 
Units + 200, Units + 201, Units + 202, Units + 203, Units + 204, 
Units + 205, Units + 206, Units + 207, Units + 208, Units + 209, 
Units + 210, Units + 211, Units + 212, Units + 213, Units + 214, 
Units + 215, Units + 216, Units + 217, Units + 218, Units + 219},
      {-1.025040, -0.715430, -0.998660, -1.959870, 1.680030, 
-1.113600, -0.207270, 0.387500, 0.195840, -2.792120, -1.179700, 
-0.288130, -0.548650, 0.994200, -0.578160, -0.516450, 0.406050, 
1.139910, -0.355620, 0.132520, -0.252290, -0.263490, -1.317970, 
0.493760, -3.211480, -0.014250, 1.999940, 0.625060, 1.079630, 
-0.284420, 0.633820, -0.733410, -1.748420, -0.509980, -0.092190, 
2.284570, -0.310250, -0.631940, 3.385140, -0.507260, 0.926100, 
-0.235980, -0.453790, 0.303220, 0.155650, -0.497530, 1.710930, 
-0.114880, 0.389560, -0.165880, 0.864610, -1.367050, -1.024740, 
-1.489280, -0.619900, 0.220810, -0.337900, 1.980530, 0.992080, 
0.327220, 0.980970, -0.500210, -0.025110, -0.919090, 1.133060, 
-0.302180, 0.413470, 0.116710, -0.898650, -0.119230, 0.038650, 
-2.254900, 0.183090, 0.196610, 1.498290}
    },
    { /* unit 222 (Old: 222) */
      0.0, -0.612210, 75,
      {Units + 145, Units + 146, Units + 147, Units + 148, Units + 
149, Units + 150, Units + 151, Units + 152, Units + 153, Units + 154, 
Units + 155, Units + 156, Units + 157, Units + 158, Units + 159, 
Units + 160, Units + 161, Units + 162, Units + 163, Units + 164, 
Units + 165, Units + 166, Units + 167, Units + 168, Units + 169, 
Units + 170, Units + 171, Units + 172, Units + 173, Units + 174, 
Units + 175, Units + 176, Units + 177, Units + 178, Units + 179, 
Units + 180, Units + 181, Units + 182, Units + 183, Units + 184, 
Units + 185, Units + 186, Units + 187, Units + 188, Units + 189, 
Units + 190, Units + 191, Units + 192, Units + 193, Units + 194, 
Units + 195, Units + 196, Units + 197, Units + 198, Units + 199, 
Units + 200, Units + 201, Units + 202, Units + 203, Units + 204, 
Units + 205, Units + 206, Units + 207, Units + 208, Units + 209, 
Units + 210, Units + 211, Units + 212, Units + 213, Units + 214, 
Units + 215, Units + 216, Units + 217, Units + 218, Units + 219},
      {0.572940, -0.408750, 0.391960, 0.109540, -1.511690, -0.855750, 
-0.057340, -0.119790, -1.014960, 1.937700, -1.335370, -0.606490, 
-0.108570, 1.835520, 0.841080, -1.798030, 0.825940, -2.198610, 
-0.069280, -0.051620, -0.916960, 0.629520, -0.785750, 2.459300, 
1.589310, -2.252010, -2.400460, -3.154280, -0.197730, -0.090180, 
-0.999050, 1.302850, 0.867780, -2.502490, 0.534750, -1.919330, 
-1.272930, -0.295760, -2.361280, 1.488770, 2.600070, -0.127780, 
0.704900, -0.052280, 0.236250, -0.861770, 3.073350, 0.399460, 
0.636280, -0.605150, -1.606860, -0.818940, 1.135510, -0.433530, 
-2.819970, 0.303900, 0.522310, -0.210470, -0.555960, -2.879570, 
3.051570, 0.650080, 0.072600, 1.705100, -1.388910, 0.109150, 
-0.889360, -0.755180, 0.863930, 2.187780, -1.201350, 0.326470, 
1.270490, -0.201230, -1.194260}
    },
    { /* unit 223 (Old: 223) */
      0.0, -1.045990, 75,
      {Units + 145, Units + 146, Units + 147, Units + 148, Units + 
149, Units + 150, Units + 151, Units + 152, Units + 153, Units + 154, 
Units + 155, Units + 156, Units + 157, Units + 158, Units + 159, 
Units + 160, Units + 161, Units + 162, Units + 163, Units + 164, 
Units + 165, Units + 166, Units + 167, Units + 168, Units + 169, 
Units + 170, Units + 171, Units + 172, Units + 173, Units + 174, 
Units + 175, Units + 176, Units + 177, Units + 178, Units + 179, 
Units + 180, Units + 181, Units + 182, Units + 183, Units + 184, 
Units + 185, Units + 186, Units + 187, Units + 188, Units + 189, 
Units + 190, Units + 191, Units + 192, Units + 193, Units + 194, 
Units + 195, Units + 196, Units + 197, Units + 198, Units + 199, 
Units + 200, Units + 201, Units + 202, Units + 203, Units + 204, 
Units + 205, Units + 206, Units + 207, Units + 208, Units + 209, 
Units + 210, Units + 211, Units + 212, Units + 213, Units + 214, 
Units + 215, Units + 216, Units + 217, Units + 218, Units + 219},
      {-1.024760, -1.866050, 0.213500, -1.549310, -0.967090, 
0.544680, -0.385380, -0.017120, -0.195420, 1.049290, 0.950250, 
-2.032770, -1.709170, 0.873250, 2.327450, 2.219200, 0.491670, 
-0.017270, 1.301560, -0.019340, -0.894820, 1.035350, -1.692170, 
-2.777180, -1.838120, -0.529940, 2.873730, 1.173710, 0.334920, 
0.149800, 1.684900, 0.281660, 1.681810, -1.045710, -0.097050, 
-1.463530, -0.264100, 0.421560, -2.526020, -0.030430, -1.488160, 
0.300280, -0.664240, -0.053830, -0.327620, 0.066530, 1.484590, 
-0.236660, 0.236740, -2.867200, -1.429240, 1.517850, -0.543760, 
-0.153160, 1.577660, -0.069700, -0.907120, -3.280270, -0.501320, 
1.635270, -2.217890, 0.665510, -2.028010, -0.173720, -0.301640, 
-0.099800, -0.520740, -1.083840, 0.659270, -1.898860, 1.108630, 
-0.244750, -1.210120, 0.036450, -0.860570}
    }

  };


  /* layer definition section (names & member units) */

  static pUnit Input[144] = {Units + 1, Units + 2, Units + 3, Units + 
4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, Units + 15, Units + 
16, Units + 17, Units + 18, Units + 19, Units + 20, Units + 21, Units 
+ 22, Units + 23, Units + 24, Units + 25, Units + 26, Units + 27, 
Units + 28, Units + 29, Units + 30, Units + 31, Units + 32, Units + 
33, Units + 34, Units + 35, Units + 36, Units + 37, Units + 38, Units 
+ 39, Units + 40, Units + 41, Units + 42, Units + 43, Units + 44, 
Units + 45, Units + 46, Units + 47, Units + 48, Units + 49, Units + 
50, Units + 51, Units + 52, Units + 53, Units + 54, Units + 55, Units 
+ 56, Units + 57, Units + 58, Units + 59, Units + 60, Units + 61, 
Units + 62, Units + 63, Units + 64, Units + 65, Units + 66, Units + 
67, Units + 68, Units + 69, Units + 70, Units + 71, Units + 72, Units 
+ 73, Units + 74, Units + 75, Units + 76, Units + 77, Units + 78, 
Units + 79, Units + 80, Units + 81, Units + 82, Units + 83, Units + 
84, Units + 85, Units + 86, Units + 87, Units + 88, Units + 89, Units 
+ 90, Units + 91, Units + 92, Units + 93, Units + 94, Units + 95, 
Units + 96, Units + 97, Units + 98, Units + 99, Units + 100, Units + 
101, Units + 102, Units + 103, Units + 104, Units + 105, Units + 106, 
Units + 107, Units + 108, Units + 109, Units + 110, Units + 111, 
Units + 112, Units + 113, Units + 114, Units + 115, Units + 116, 
Units + 117, Units + 118, Units + 119, Units + 120, Units + 121, 
Units + 122, Units + 123, Units + 124, Units + 125, Units + 126, 
Units + 127, Units + 128, Units + 129, Units + 130, Units + 131, 
Units + 132, Units + 133, Units + 134, Units + 135, Units + 136, 
Units + 137, Units + 138, Units + 139, Units + 140, Units + 141, 
Units + 142, Units + 143, Units + 144}; /* members */

  static pUnit Hidden1[75] = {Units + 145, Units + 146, Units + 147, 
Units + 148, Units + 149, Units + 150, Units + 151, Units + 152, 
Units + 153, Units + 154, Units + 155, Units + 156, Units + 157, 
Units + 158, Units + 159, Units + 160, Units + 161, Units + 162, 
Units + 163, Units + 164, Units + 165, Units + 166, Units + 167, 
Units + 168, Units + 169, Units + 170, Units + 171, Units + 172, 
Units + 173, Units + 174, Units + 175, Units + 176, Units + 177, 
Units + 178, Units + 179, Units + 180, Units + 181, Units + 182, 
Units + 183, Units + 184, Units + 185, Units + 186, Units + 187, 
Units + 188, Units + 189, Units + 190, Units + 191, Units + 192, 
Units + 193, Units + 194, Units + 195, Units + 196, Units + 197, 
Units + 198, Units + 199, Units + 200, Units + 201, Units + 202, 
Units + 203, Units + 204, Units + 205, Units + 206, Units + 207, 
Units + 208, Units + 209, Units + 210, Units + 211, Units + 212, 
Units + 213, Units + 214, Units + 215, Units + 216, Units + 217, 
Units + 218, Units + 219}; /* members */

  static pUnit Output1[4] = {Units + 220, Units + 221, Units + 222, 
Units + 223}; /* members */

  static int Output[4] = {220, 221, 222, 223};

  for(member = 0; member < 144; member++) {
    Input[member]->act = in[member];
  }

  for (member = 0; member < 75; member++) {
    unit = Hidden1[member];
    sum = 0.0;
    for (source = 0; source < unit->NoOfSources; source++) {
      sum += unit->sources[source]->act
             * unit->weights[source];
    }
    unit->act = Act_Logistic(sum, unit->Bias);
  };

  for (member = 0; member < 4; member++) {
    unit = Output1[member];
    sum = 0.0;
    for (source = 0; source < unit->NoOfSources; source++) {
      sum += unit->sources[source]->act
             * unit->weights[source];
    }
    unit->act = Act_Logistic(sum, unit->Bias);
  };

  for(member = 0; member < 4; member++) {
    out[member] = Units[Output[member]].act;
  }

  return(OK);
}

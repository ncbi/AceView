/*
Code from 1993,  by Hank Wallace, downloaded from the web 2008_01_25:
This source code accompanies the article, "Using The Golay Error Detection And
Correction Code", by Hank Wallace. This program demonstrates use of the Golay
code.

   Usage: G DATA Encode/Correct/Verify/Test

     where DATA is the data to be encoded, codeword to be corrected,
     or codeword to be checked for errors. DATA is hexadecimal.

   Examples:

     G 555 E      encodes information value 555 and prints a codeword
     G ABC123 C   corrects codeword ABC123
     G ABC123 V   checks codeword ABC123 for errors
     G ABC123 T   tests routines, ABC123 is a dummy parameter

This program may be freely incorporated into your programs as needed. It
compiles under Borland's Turbo C 2.0. No warranty of any kind is granted.

*/

#include "stdio.h"
#include "stdio.h"

#define POLY  0xAE3  /* or use the other polynomial, 0xC75 */

/* ====================================================== */

unsigned long golay(unsigned long cw)
/* This function calculates [23,12] Golay codewords.
   The format of the returned longint is
   [checkbits(11),data(12)]. */
{
  int i;
  unsigned long c;
  cw&=0xfffl;
  c=cw; /* save original codeword */
  for (i=1; i<=12; i++)  /* examine each data bit */
    {
      if (cw & 1)        /* test data bit */
        cw^=POLY;        /* XOR polynomial */
      cw>>=1;            /* shift intermediate result */
    }
  return((cw<<12)|c);    /* assemble codeword */
}

/* ====================================================== */

static int parity(unsigned long cw)
/* This function checks the overall parity of codeword cw.
   If parity is even, 0 is returned, else 1. */
{
  unsigned char p;

  /* XOR the bytes of the codeword */
  p=*(unsigned char*)&cw;
  p^=*((unsigned char*)&cw+1);
  p^=*((unsigned char*)&cw+2);

  /* XOR the halves of the intermediate result */
  p=p ^ (p>>4);
  p=p ^ (p>>2);
  p=p ^ (p>>1);

  /* return the parity result */
  return(p & 1);
}

/* ====================================================== */

static unsigned long syndrome(unsigned long cw)
/* This function calculates and returns the syndrome
   of a [23,12] Golay codeword. */
{
  int i;
  cw&=0x7fffffl;
  for (i=1; i<=12; i++)  /* examine each data bit */
    {
      if (cw & 1)        /* test data bit */
        cw^=POLY;        /* XOR polynomial */
      cw>>=1;            /* shift intermediate result */
    }
  return(cw<<12);        /* value pairs with upper bits of cw */
}

/* ====================================================== */

static int weight(unsigned long cw)
/* This function calculates the weight of
   23 bit codeword cw. */
{
  int bits,k;

  /* nibble weight table */
  const char wgt[16] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

  bits=0; /* bit counter */
  k=0;
  /* do all bits, six nibbles max */
  while ((k<6) && (cw))
    {
      bits=bits+wgt[cw & 0xf];
      cw>>=4;
      k++;
    }

  return(bits);
}

/* ====================================================== */

static unsigned long rotate_left(unsigned long cw, int n)
/* This function rotates 23 bit codeword cw left by n bits. */
{
  int i;

  if (n != 0)
    {
      for (i=1; i<=n; i++)
        {
          if ((cw & 0x400000l) != 0)
            cw=(cw << 1) | 1;
          else
            cw<<=1;
        }
    }

  return(cw & 0x7fffffl);
}

/* ====================================================== */

static unsigned long rotate_right(unsigned long cw, int n)
/* This function rotates 23 bit codeword cw right by n bits. */
{
  int i;

  if (n != 0)
    {
      for (i=1; i<=n; i++)
        {
          if ((cw & 1) != 0)
            cw=(cw >> 1) | 0x400000l;
          else
            cw>>=1;
        }
    }

  return(cw & 0x7fffffl);
}

/* ====================================================== */

static unsigned long golay_correct (unsigned long cw, int *errs)
/* This function corrects Golay [23,12] codeword cw, returning the
   corrected codeword. This function will produce the corrected codeword
   for three or fewer errors. It will produce some other valid Golay
   codeword for four or more errors, possibly not the intended
   one. *errs is set to the number of bit errors corrected. */
{
  unsigned char
    w;                /* current syndrome limit weight, 2 or 3 */
  unsigned long
    mask;             /* mask for bit flipping */
  int
    i,j;              /* index */
  unsigned long
    s,                /* calculated syndrome */
    cwsaver;          /* saves initial value of cw */

  cwsaver=cw;         /* save */
  *errs=0;
  w=3;                /* initial syndrome weight threshold */
  j=-1;               /* -1 = no trial bit flipping on first pass */
  mask=1;
  while (j<23) /* flip each trial bit */
    {
      if (j != -1) /* toggle a trial bit */
        {
          if (j>0) /* restore last trial bit */
            {
              cw=cwsaver ^ mask;
              mask+=mask; /* point to next bit */
            }
          cw=cwsaver ^ mask; /* flip next trial bit */
          w=2; /* lower the threshold while bit diddling */
        }

      s=syndrome(cw); /* look for errors */
      if (s) /* errors exist */
        {
          for (i=0; i<23; i++) /* check syndrome of each cyclic shift */
            {
              if ((*errs=weight(s)) <= w) /* syndrome matches error pattern */
                {
                  cw=cw ^ s;              /* remove errors */
                  cw=rotate_right(cw,i);  /* unrotate data */
                  return(s=cw);
                }
              else
                {
                  cw=rotate_left(cw,1);   /* rotate to next pattern */
                  s=syndrome(cw);         /* calc new syndrome */
                }
            }
          j++; /* toggle next trial bit */
        }
      else
        return(cw); /* return corrected codeword */
    }

  return(cwsaver); /* return original if no corrections */
} /* golay_correct */

/* ====================================================== */

int golay_decode(int correct_mode, int *errs, unsigned long *cw)
/* This function decodes codeword *cw in one of two modes. If correct_mode
   is nonzero, error correction is attempted, with *errs set to the number of
   bits corrected, and returning 0 if no errors exist, or 1 if parity errors
   exist. If correct_mode is zero, error detection is performed on *cw,
   returning 0 if no errors exist, 1 if an overall parity error exists, and
   2 if a codeword error exists. */
{
  unsigned long parity_bit;

  if (correct_mode)               /* correct errors */
    {
      parity_bit=*cw & 0x800000l; /* save parity bit */
      *cw&=~0x800000l;            /* remove parity bit for correction */

      *cw=golay_correct(*cw, errs);     /* correct up to three bits */
      *cw|=parity_bit;            /* restore parity bit */

      /* check for 4 bit errors */
      if (parity(*cw))            /* odd parity is an error */
        return(1);
      return(0); /* no errors */
    }
  else /* detect errors only */
    {
      *errs=0;
      if (parity(*cw)) /* odd parity is an error */
        {
          *errs=1;
          return(1);
        }
      if (syndrome(*cw))
        {
          *errs=1;
          return(2);
        }
      else
        return(0); /* no errors */
    }
} /* golay_decode */

/* ====================================================== */

static void golay_test(void)
/* This function tests the Golay routines for detection and correction
   of various patterns of error_limit bit errors. The error_mask cycles
   over all possible values, and error_limit selects the maximum number
   of induced errors. */
{
  unsigned long
    error_mask,         /* bitwise mask for inducing errors */
    trashed_codeword,   /* the codeword for trial correction */
    virgin_codeword;    /* the original codeword without errors */
  unsigned char
    pass=1,             /* assume test passes */
    error_limit=3;      /* select number of induced bit errors here */
  int
    error_count;        /* receives number of errors corrected */

  virgin_codeword=golay(0x555); /* make a test codeword */
  if (parity(virgin_codeword))
    virgin_codeword^=0x800000l;
  for (error_mask=0; error_mask<0x800000l; error_mask++)
    {
      /* filter the mask for the selected number of bit errors */
      if (weight(error_mask) <= error_limit) /* you can make this faster! */
        {
          trashed_codeword=virgin_codeword ^ error_mask; /* induce bit errors */

          golay_decode(1,&error_count,&trashed_codeword); /* try to correct bit errors */

          if (trashed_codeword ^ virgin_codeword)
            {
              printf("Unable to correct %d errors induced with error mask = 0x%lX\n",
                weight(error_mask),error_mask);
              pass=0;
            }

          if (1 /*kbhit()*/) /* look for user input */
            {
              if (getchar() == 27) return; /* escape exits */

              /* other key prints status */
              printf("Current test count = %ld of %ld\n",error_mask,0x800000l);
            }
        }
    }
  printf("Golay test %s!\n",pass?"PASSED":"FAILED");
}

/* ====================================================== */

int golay_main(int argument_count, char *argument[])
{
  int i,j;
  unsigned long l,g;
  const char *errmsg = "Usage: G DATA Encode/Correct/Verify/Test\n\n"
             "  where DATA is the data to be encoded, codeword to be corrected,\n"
             "  or codeword to be checked for errors. DATA is hexadecimal.\n\n"
             "Examples:\n\n"
             "  G 555 E      encodes information value 555 and prints a codeword\n"
             "  G ABC123 C   corrects codeword ABC123\n"
             "  G ABC123 V   checks codeword ABC123 for errors\n"
             "  G ABC123 T   tests routines, ABC123 is a dummy parameter\n\n";

  if (argument_count != 3)
    {
      printf(errmsg);
      return 0 ;
    }

  if (sscanf(argument[1],"%lx",&l) != 1)
    {
      printf(errmsg);
      return 0 ;
    }

  switch (toupper(*argument[2]))
    {
    case 'U': 
      {
	long int nn = 0, ii ;
	for (ii = 0 ; ii < 0x800000 ; ii++)
	  {
	    l = ii ;
	    if (!golay_decode(0,&i,&l)) nn++ ;
	  }
	printf ("Out of %ld words, %ld are Golay\n", ii, nn) ;
      }
      case 'E': /* encode */
        l&=0xfff;
        l=golay(l);
        if (parity(l)) l^=0x800000l;
        printf("Codeword = %lX\n",l);
        break;

      case 'V': /* verify */
        if (golay_decode(0,&i,&l))
          printf("Codeword %lX is not a Golay codeword.\n",l);
        else
          printf("Codeword %lX is a Golay codeword.\n",l);
        break;

      case 'C': /* correct */
        g=l; /* save initial codeword */
        j=golay_decode(1,&i,&l);
        if ((j) && (i))
          printf("Codeword %lX had %d bits corrected,\n"
                 "resulting in codeword %lX with a parity error.\n",g,i,l);
        else
          if ((j == 0) && (i))
            printf("Codeword %lX had %d bits corrected, resulting in codeword %lX.\n",g,i,l);
          else
            if ((j) && (i == 0))
              printf("Codeword %lX has a parity error. No bits were corrected.\n",g);
            else
              if ((j == 0) && (i == 0))
                printf("Codeword %lX does not require correction.\n",g);
        break;

      case 'T': /* test */
        printf("Press SPACE for status, ESC to exit test...\n");
        golay_test();
        break;

      default:
        printf(errmsg);
        return 0 ;
    }
  return 0 ;
}

/* end of G.C */

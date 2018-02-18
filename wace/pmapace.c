/*  File: pmapace.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: stand alone program
 	translates results of ascify.for on vax into ace file
 * HISTORY:
 * Last edited: Sep  8 17:00 1994 (srk)
 * * Jan 29 23:11 1992 (rd): reads vaxmap entries
 * Created: ages ago
 *-------------------------------------------------------------------
 */

/* $Id: pmapace.c,v 1.1.1.1 2002/07/19 20:23:11 sienkiew Exp $ */
 
      /***************************************************************/
      /***************************************************************/
      /**  File pmapace.c :                                         **/
      /**  standalone preprocessor for the ACeDB program.           **/
      /**  Reads Sulston's physical map data.                       **/
      /***************************************************************/
      /***************************************************************/
      /*                                                             */
      /*         R.Durbin & J.Thierry-Mieg.                          */
      /*         last modified  5/11/90  by RMD                      */
      /*         heavily changed to fit new ascify.for on the VAX    */
      /*                                                             */
      /***************************************************************/
 
#include "regular.h"	/* standard declarations */

#define MAXBUF 256
#define MAXNCLONE 25000
 
/***************************************************/

int   nline = 0 ;
char  buffer[MAXBUF] ;
FILE *in ;

void getRecord (void)
{
  char *cp ;

  for (cp = buffer ; ; ++cp)
    { *cp = fgetc(in) ;
      if (*cp == '\n' || *cp == EOF)
	break ;
      if (*cp == '\"')
	*cp = '\'' ;   /* RMD this fix looks OK; else \" sets counting wrong */
    }
  *cp = 0 ;
  ++nline ;
}

/****************************************************/

char* cloneName[MAXNCLONE] ;

void getCloneNames (void)	/* pass once over data to get names */
{
  char *cp ;
  int  kclone = 2 ;
 
  getRecord () ;		/* first line is comment */

  while (TRUE)
    { do getRecord () ; 
        while (*buffer && buffer[1] != 'C') ;
      if (!*buffer)
	break ;
      
      for (cp = &buffer[3+8] ; *--cp == ' ' ; ) ; 
      *++cp = 0 ;
      cloneName[kclone] = (char*) messalloc (strlen (&buffer[3]) + 1) ;
      strcpy (cloneName[kclone++],&buffer[3]) ;
    }

  printf ("%d clone names read in\n",kclone-2) ;
  rewind (in) ;
  nline = 0 ;
}

/****************************************************/

void main (int argc, char **argv)
{
  FILE  *out ;
static char  *chromName[] = {"","I","II","III","IV","V","X"} ;
  int   left_pos_in_contig, length_of_clone, k_canonical, flag ;
  int   nb_of_bands, first_band, gel_number, contig_number ;
  int   chrom, insitu_left, insitu_right ;
  char  clone[9],gene[9],remark[41],*acedb_data,*cp ;
  int   nclone = 0 ;
  float vaxmap ;
 
  if (!(acedb_data = getenv ("ACEDB_DATA")))
    messcrash ("Env. variable ACEDB_DATA must be set to raw data directory") ;
  if (argc != 2)
    messcrash (
      "Usage: pmapace project : raw file is $(ACEDB_DATA)/pmap/xxx.asc\n") ;
 
  in = fopen (messprintf ("%s/pmap/%s.asc",acedb_data,argv[1]),"r") ;
  out = fopen (messprintf ("%s/pmap/%s.ace",acedb_data),"w") ;
 
  if (!in || !out)
    goto close;

  getCloneNames () ;
 
       /* The first record is a comment containing the project number, etc */
  getRecord () ;
  fprintf (out,"%s\n",buffer) ;

  clone[8] = gene[8] = remark[40] = 0 ;

  getRecord () ;
  while (*buffer)
    {
      if (sscanf (buffer," C %8c %i %i %i %i %i %i %i %i",
	      clone, &nb_of_bands, &first_band, &gel_number, &contig_number,
	      &left_pos_in_contig, &length_of_clone, &k_canonical, &flag) != 9)
	messcrash ("Bad input line %d : %s",nline,buffer) ;

      for (cp = &clone[8] ; *--cp == ' ' && cp >= clone ; ) ;
      *++cp = 0 ;

            /*Show the result*/
      if (!(nclone%1000))
        printf ("\n%d : %s",nclone,clone) ;
      else if (!(nclone%20))
        { putchar ('.') ;
	  fflush (stdout) ;	/* to make it appear immediately */
	}
 
              /* Write to ace files */
      if  (*clone != '!')          /* deleted clone */
        { fprintf (out,"Clone \"%s\"\n",clone) ;
 
          if (nb_of_bands)
            { fprintf (out,"gel_number : %d\n",gel_number) ;
              fprintf (out,"bands : %d %d\n",first_band,nb_of_bands) ;
            }
 
          if (k_canonical > 10000000)
            fprintf (out,"funny_match_to : %s\n",
                     cloneName[k_canonical-10000000]) ;
          else if (k_canonical > 0)
            fprintf (out,"exact_match_to : %s\n",cloneName[k_canonical]) ;
          else if (k_canonical < -1)
            fprintf (out,"approximate_match_to : %s\n",
		     cloneName[-k_canonical]) ;
          else if (contig_number)
	    fprintf (out,"pmap : ctg%d %d %d\n",contig_number,
		     left_pos_in_contig,left_pos_in_contig+length_of_clone) ;
 
          if (flag)
            { if (flag & 1) fprintf (out,"autopos\n") ;
              if (flag & 16) fprintf (out,"cosmid_grid\n") ;
              if (flag & 32) fprintf (out,"canon_for_cosmid\n") ;
              if (flag & 64) fprintf (out,"YAC_polytene : 1\n") ;
              if (flag & 2 || flag & 4)
                fprintf (out,"flag : %d\n",flag) ;
            }
 
	  getRecord () ;
	  while (*buffer && buffer[1] == 'R')
	    {
	      if (sscanf (&buffer[3+8+40]," %i %i %i %f",
			  &chrom,&insitu_left,&insitu_right,&vaxmap) != 4)
		messcrash ("Bad input line %d : %s",nline,buffer) ;

	      strncpy (gene,&buffer[3],8) ;
	      for (cp = &gene[8] ; *--cp == ' ' && cp >= gene ; ) ;
	      *++cp = 0 ;
	      
	      strncpy (remark,&buffer[11],40) ;
	      for (cp = &remark[40] ; *--cp == ' ' && cp >= remark ; ) ;
	      *++cp = 0 ;

	      if (*gene)
                fprintf (out,"gene : \"%s\"\n",gene) ;

              if (*remark)
		if (strstr (remark, "PCR"))
		  fprintf (out,"PCR_remark : \"%s\"\n",remark) ;
		else if (remark[0] == 'Y')
		  fprintf (out,"Y_remark : \"%s\"\n",remark) ;
		else
		  fprintf (out,"General_remark : \"%s\"\n",remark) ;
 
              if (chrom)
                fprintf (out,"chromosome : %s\n",chromName[chrom]) ;
 
              if (insitu_right)
                fprintf (out,"in_situ : %d %d\n",insitu_left,insitu_right) ;

	      if (vaxmap)
		fprintf (out,"vaxmap : %6.2f\n",vaxmap) ;

	      getRecord () ;
            }
 

          fprintf (out,"\n") ;
	  ++nclone ;
        }
      else
	do getRecord () ; while (*buffer && buffer[1] == 'R') ;

    }
 
close :
  printf ("Processed %d active clones\n",nclone) ;
}

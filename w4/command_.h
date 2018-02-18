/*  File: command_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2000
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Internal header for ace commands.
 * HISTORY:
 * Last edited: Oct  5 12:40 2000 (edgrif)
 * Created: Thu May 11 19:43:40 2000 (edgrif)
 * CVS info:   $Id: command_.h,v 1.2 2011/02/24 00:42:26 mieg Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_COMMAND_PRIV_H
#define DEF_COMMAND_PRIV_H

#include <wh/command.h>					    /* command public header. */


/* Callers are returned an opaque reference to this struct, which contains   */
/* all the state associated with a command.                                  */
/*                                                                           */
typedef struct _AceCommandStruct
{ 
  magic_t *magic ;					    /* == &AceCommand_MAGIC */
  ACEOUT dump_out;					    /* used for show -f / aql -o etc 
							     * used to be outLevel/outfil */
  BOOL noWrite ;
  SPREAD spread ; 
  KEYSET ksNew, ksOld, kA ;
  Stack ksStack ;
  int nns ;						    /* number of entries on stack */
  unsigned int choix ;
  unsigned int perms;
  Array aMenu ;
  KEY lastCommand ;
  BOOL showStatus , beginTable ;
  BOOL showTableTitle;
  int nextObj, lastSpread ; 
  int cumulatedTableLength ;
  int encore ;
  int minN ;						    /* startup */
  int maxN ;						    /* -1 or, max number of desired output */
  int notListable ;					    /* empty and alias objects */
  BOOL dumpTimeStamps ;
  BOOL dumpComments ;
  BOOL noXref ;
  BOOL allowMismatches ;
  BOOL noClassName ;
  COND showCond ;
  char beauty ;						    /* ace, tace, perl */
  char fastaCommand ; /* for DNA, Peptide more */
} AceCommandStruct ;



/* Note the restricted possibilities here...                                 */
/* CMD_FILESPEC_NONE => cmd does not use files at all                        */
/* CMD_FILESPEC_OPTIONAL                                                     */
/*                   => if file specified, must be last word on line.        */
/* CMD_FILESPEC_OPTIONAL & CMD_FILESPEC_FLAG                                 */
/*                   => will be word following flag on line.                 */
/* CMD_FILESPEC_COMPULSORY                                                   */
/*                   => file must be specified, must be last word on line.   */
typedef enum _CmdFilespecType
{
  CMD_FILESPEC_NONE         = 0x0,			    /* cmd. does not use files. */
  CMD_FILESPEC_FLAG	    = 0x1,			    /* flag used to signal file. */
  CMD_FILESPEC_OPTIONAL     = 0x2,			    /* file is optional. */
  CMD_FILESPEC_COMPULSORY   = 0x4			    /* file is compulsory (must be last
							       word on command line). */
} CmdFilespecType ;


/* Describes how files are specified within acedb command line commands.     */
/* Notes:                                                                    */
/*                                                                           */
typedef struct _CmdFilespecStruct
{
  BOOL any_files ;					    /* Does cmd have any files ? */
  CmdFilespecType infile ;				    /* How is file specified ? */
  char *infile_flag ;					    /* flag, e.g. "-f", */
  CmdFilespecType outfile ;				    /* How is file specified ? */
  char *outfile_flag ;					    /* flag, e.g. "-f", */
} CmdFilespecStruct ;



/* Command struct, defines format and use of an ace command.                 */
/* choix = one of the following:
   CHOIX_UNIVERSAL: universal
   CHOIX_NONSERVER: non-server
      CHOIX_SERVER: server
         CHOIX_GIF: gif
*/
typedef struct
{
  unsigned int choix ;					    /* Which command "sets" this command
							       belongs to. */
  unsigned int perm ;					    /* What level of access (e.g. read) is */
							    /* required to run the command.*/
  KEY key ;						    /* id by which command is known. */
  char *text ;						    /* help text for command. */
  CmdFilespecStruct filespec ;				    /* How in/out files are specified. */
} CHOIX ;




/* Internal to the acecommand package.                                       */
/*                                                                           */
CHOIX *uAceCommandGetChoixMenu(void) ;			    /* return ptr to array of all possible */
							    /* commands. */

Array uAceCommandChoisirMenu(unsigned int choix, unsigned int perms) ;
							    /* Return an array of allowed commands */
							    /* for a given command subset and
							       permission level.*/



#endif /* DEF_COMMAND_PRIV_H */

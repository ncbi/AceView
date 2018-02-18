/*  File: command.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: header file for functions in command.c
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 24 22:22 2001 (edgrif)
 * * Nov 23 09:46 1999 (edgrif): Insert new choixmenu values for admin cmds.
 * Created: Thu Aug 26 17:06:53 1999 (fw)
 * CVS info:   $Id: command7_.h,v 1.3 2011/02/25 22:39:02 mieg Exp $
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_COMMAND__7H
#define ACEDB_COMMAND__7H

#include "acedb.h"
#include "query.h"		/* for COND struct */
#include "spread.h"		/* for SPREAD struct */
#include "dict.h"

#define COMMAND_MAGIC 248035

/* Commands in acedb are grouped, these groups can be summed to produce sets */
/* of commands for particular users/executables (n.b. this is a once only    */
/* operation, we don't turn these on/off like normal flags).                 */
/* (see w3/tacemain.c as an example.)                                        */
/*                                                                           */
/*                                                                           */
/* CHOIX_RPC is a temporary hack while the rpc server is still supported.    */
/* The RPC server only includes CHOIX_RPC, while the socket server includes  */
/* CHOIX_RPC and CHOIX_SERVER. Once the rpc server goes we can adjust the    */
/* code to only use CHOIX_SERVER.                                            */
/*                                                                           */
typedef enum CmdChoix_
{
  CHOIX_UNDEFINED =  0U,
  CHOIX_UNIVERSAL =  1U,
  CHOIX_NONSERVER =  2U,
  CHOIX_SERVER    =  4U,
  CHOIX_GIF       =  8U,
  CHOIX_RPC       = 16U					    /* Temporary hack while rpc server is
							       still used. */
} CmdChoix ;

typedef enum CmdPerm_
{
  PERM_UNDEFINED =  0U,
  PERM_READ      =  1U,
  PERM_WRITE     =  2U,
  PERM_ADMIN     =  4U
} CmdPerm ;


/* Each command has a key by which it is known, currently this is defined as */
/* a KEY, but this is not correct because a KEY is an unsigned int and the   */
/* code uses negative values. This will change eventually, for now we just   */
/* use negative values and it all comes out in the wash.                     */
/*                                                                           */
/* Special values:                                                           */
enum {CMD_NOTVALID = 0, CMD_FAILED = -99, CMD_TERMINATE = -1} ;


/* Set of symbols that are the KEYs for general admin commands.              */
/* These commands are in the 400 range, other commands should not use this   */
/* range.                                                                    */
enum {CMD_NEWLOG = 400} ;


/* Set of symbols that are the KEYs for server only commands.                */
/* These commands are in the 500 range, other commands should not use this   */
/* range.                                                                    */
enum {CMD_SHUTDOWN = 500, CMD_WHO, CMD_VERSION, CMD_WSPEC} ;

/* Set of symbols that are the KEYs for server admin commands.               */
/* These commands are in the 600 range, other commands should not use this   */
/* range.                                                                    */
enum {CMD_PASSWD = 600, CMD_USER, CMD_GLOBAL, CMD_DOMAIN,
      CMD_REMOTEPARSE, CMD_REMOTEPPARSE, CMD_REMOTEDUMP, CMD_SERVERLOG} ;


/* Return result from the file arg finding call.                             */
/* CMD_FILE_NONE =>  cmd does not use any files OR filename not specified.   */
/* CMD_FILE_FOUND => a filename was specified correctly (check returned args */
/*                   for input and/or output file).                          */
/* CMD_FILE_ERR   => filename specified but incorrectly.                     */
/*                                                                           */
typedef enum CmdFilespecArg_ {CMD_FILE_NONE, CMD_FILE_FOUND, CMD_FILE_ERR } CmdFilespecArg ;

typedef struct _AceCommandStruct
{ 
  int magic ;
  AC_HANDLE h ;
  BOOL noWrite ;
  SPREAD spread ; 
  KEYSET ksNew, ksOld, kA ;
  Stack ksStack ;
  int nns ; /* number of entries on stack */
  int choix ;
  Array aMenu ;
  KEY lastCommand ;
  BOOL isPerl, showStatus , beginTable ;
  int nextObj, lastSpread ; 
  int cumulatedTableLength ;
  int encore ;
  int minN ; /* startup */
  int maxN ; /* -1 or, max number of desired output */
  int notListable ; /* empty and alias objects */
  BOOL dumpTimeStamps ;
  BOOL dumpComments ;
  BOOL noXref ;
  BOOL allowMismatches ;
  BOOL noClassName ;
  COND showCond ;
  int outLevel ;
  FILE *outfile ;
  char beauty ; /* ace, tace, perl */
	/*
	* For beauty, possible values are:  (see w4/dump.c dumpKey2)
	*  p  perl
	*  j  java format
	*  J  another java format
	*  a  ace file format
	*  A  ace file format (same?)
	*  C  Vahan's acec library format
	*  H  ???
	*  h  ???
	*  m  minilib format - much like ace file but shows all tags
	*     and has special markings for type A objects
	*
	* Sanger also uses: kKxXyY
	* Sanger has an XML output format
	*  
	*/
  char fastaCommand ; /* for DNA, Peptide more */
  DICT *dict ;
  Array namedSets ;
} *COMMAND_LOOK ;

#endif /* ACEDB_COMMAND__7H */

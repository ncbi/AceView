/*  File: acelib.c
 *  Author: Jean Thierry-Mieg (mieg@crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: acelib.c,v 1.10 2008/11/18 02:53:54 mieg Exp $
 * Description: public programming interace to acedb
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  5 17:41 1999 (fw)
 * Created: Sun Nov  9 23:02:41 1997 (rd)
 *-------------------------------------------------------------------
 */

/* First try at actual implementation of the public
ace interface

in this implemantation i suppose that i will link this
code into xace and use it to write say the genetic map

so i do not yet care about the ACE context which is the default

In a separate implementation, i will put all this in a clint
code and then have to open the database etc

*/


/*******************************************************************/
/*******************  session managment ****************************/

/*******************************************************************/
/*******************************************************************/
/* aceQuit needs work, in client server mode, you would not
   exit(), just close connection.
   also if commit = true and this does not work, there should
   be an error message 
   */

/*********************************************************************/
/*********************************************************************/

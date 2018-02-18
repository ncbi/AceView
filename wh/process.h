/* File: process.h
 *  Authors: Petr Kocab    (p.kocab@dkfz-heidelberg.de)
 *           Martin Senger (m.senger@dkfz-heidelberg.de)
 *
 *-------------------------------------------------------------------
 * Description:
 *    Deals with child processes, and signal handlings.
 *
 * Exported functions:
 *    
 *
 * HISTORY:
 *    Created: Thr Apr 20 1995 (Senger, Kocab)
 *-------------------------------------------------------------------
 */
/* $Id: process.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#ifndef DEFINE_PROCESS_h
#define DEFINE_PROCESS_h
/******************************************************************************
*
*   --- include files --- 
*
******************************************************************************/
#include <sys/wait.h>
#include <signal.h>
#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include "mystdlib.h"
/*
#include <stdlib.h>
*/

/******************************************************************************
*
*   --- constants, literals --- 
*
******************************************************************************/

#define MAX_PATH         256

/* --- names of the environment variables --- */
#define MOSAIC_SLEEP    "ACEDB_SLEEP_MOSAIC"
#define HTML_BROWSER    "ACEDB_HTML_BROWSER"


#define SLEEP_MOSAIC     atoi(getenv(MOSAIC_SLEEP) ? getenv(MOSAIC_SLEEP) : "5")
#define DEFAULT_HTML_BROWSER (getenv(HTML_BROWSER) ? getenv (HTML_BROWSER) : "Mosaic")

/* --- patterns for file names --- */
#define MOSAIC_FILE     "/tmp/Mosaic.%d"                   /* %d <== Mosaic_pid */

/* --- special identifiers in table of processes --- */
#define TCP_ITEM_FREE   (-1)
#define TCP_ITEM_USED   0
#define TCP_ITEM_HELP   1
#define TCP_ITEM_ACTION 2


/* --- child process status --- */
#define CHILD_STOPPED       0     /* process is stopped */
#define CHILD_EXITED        1     /* process is exited */
#define CHILD_TERMINATED    2     /* process is terminated */
#define CHILD_C_TERMINATED  3     /* process is terminated with core */
#define CHILD_NOT_BORN      4     /* process is not existing */
#define CHILD_ACTIVE        5     /* process is active */

/******************************************************************************
*
*   --- structures, typedefs --- 
*
******************************************************************************/

/* --- each child process has this entry in the table of child processes --- */
typedef struct {
   int      id;            /* internal process description */
   int      pid;           /* system identification of the process */
   int      status;        /* process status */
   int      code;          /* process exit code or terminating signal */
   time_t   last_change;   /* time of the last status change */
   void     (* fce)(int);  /* user function */
   } Child, *PtrChild;

/******************************************************************************
*
*   --- function prototypes --- 
*
******************************************************************************/
int  init_ipc (void);
void register_child (int pid, int id, void (* fce)(int));
int  get_process_info (int pid, int id, PtrChild info);
int  exec_child(char **cmd_line, int id, void(*fce)(int), 
	       int in_fd, int out_fd);
void suspend_sigchld();
void resume_sigchld ();

#endif
 

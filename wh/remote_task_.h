
/*  File: remote_task_.h
 *  channels interface : derived from the GO language channels 
 *  Author: Danielle and Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov) Mark Sienkiewicz (sienkiew@stsci.edu)
 *  Copyright (C) J Thierry-Mieg, 2015
 * -------------------------------------------------------------------
 * Created: mars 2015, mieg
 *-------------------------------------------------------------------
 */

#ifndef REMOTE_TASK__H
#define REMOTE_TASK__H

struct taskStruct { 
  int *magic ; 
  AC_HANDLE h ; 
  const char *program_name ;
  TASK_CONFIG config ;
  BOOL clientSide ;
  BOOL serverSide ;
  const char *host ; 
  const char *command ;
  char *error ; 
  int port, timeout ; 
  int client_fd ;
  struct acetcp_descriptor *d ; 
  DICT *channelDict ;
  Array channels ;
  Array clientTasks ;
  CHAN *lifeChannel ;
  CHAN *writer ;
  char *acetcp_response ;
  int acetcp_response_length ;
  int acetcp_response_safe ;  /* length of the currently allocated acetcp_response */
  pthread_mutex_t mutex ;
  void *wegoListener ;
  int signature ;
  int lifeCycle ;
  void *group ; /* WEGO_GROUP */
} ;

#endif

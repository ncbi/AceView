/*  Last edited: May 15 15:03 1992 (mieg) */

/* $Id: igdevent.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */
/*
 created March 25, 1992
 by Petr Kocab      dok416@cvx12.dkfz-heidelberg.de



This is the header file for functions enabling communication between parent
and child processes over sockets. The parent process forks and establishes
a "permanent" communication channel with the execed child. The channel consists
of two socketpairs, one used for data transfer and the other for transfer control.

The parent can communicate with MAX_EVENT_CHANNEL = 5 children.
Both ends may send and receive events. An event has the header and the data (body).
The header contains information on the event type and the data transfer type. See
the fredscEvent structure.

*/


#ifndef _IGDEVENT_
#define _IGDEVENT_

#ifndef _IGDCALL_

#include <stdio.h>
#include <errno.h>
#ifndef DEF_MYSTDLIB_H
#include <stdlib.h>     /* J.Mieg has his own definition stdlib.h */
#endif
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>

#endif

#define MAX_COMMANDS        2

          /* to distinguish parent and child ends of socketpairs */
#define MAIN_PROC   0
#define CHILD_PROC  1

typedef int FredEventType;

          /* allowed values: */
#define ERROR_MESSAGE       ((FredEventType)(1 << 0 )) 
#define DISPLAY_MESSAGE     ((FredEventType)(1 << 1 ))
#define ACE_QUERY           ((FredEventType)(1 << 2 ))
#define ACE_GREP            ((FredEventType)(1 << 3 ))
#define ACE_DATA            ((FredEventType)(1 << 4 ))
#define ACE_DATA_REQUEST    ((FredEventType)(1 << 5 ))
#define ACE_MODELS          ((FredEventType)(1 << 6 ))
#define ACE_MODEL_REQUEST   ((FredEventType)(1 << 7 ))
#define STATUS_REQUEST      ((FredEventType)(1 << 8 ))
#define QUIT_CHILD          ((FredEventType)(1 << 9 ))
#define DATA_OK             ((FredEventType)(1 << 10))
#define DATA_NOT_OK         ((FredEventType)(1 << 11))

          /* system events      */
#define END_TRANSFER         0
#define ERROR_IN_TRANSFER   -1

          /* what will go thru the channel?  */
#define FILENAME            1
#define DATA                2

         /* child filedesc's to be used in this example  */
#define FD_CHILD_0          4
#define FD_CHILD_1          5

         /* error message                    */
#define C_EM_KILL_EVENT_TYPE     "event type <%il> was kill\n"

typedef int DataType;
          /* allowed values: */
#define NONE_               0
#define DATA_IN_FILE        1
#define DATA_IN_FILE_DESC   2
#define DATA_IN_MEM_BUF     3

#define MAX_EVENT_CHANNEL   5         /* max # of simultaneous channels allowed     */
#define LEN_MY_BUF          1024*5    /* length of the buffer for read(), write()   */
#define MAX_LEN_FILE_NAME   128
#define INIT_SIZE           1024      /* initial size for malloc()                  */



typedef struct {
    char      *bptr;
    long      blen;
    } MemBuff ;      /* memory buffer */

typedef union {
    char      *filename;
    int       filedesc;
    MemBuff   buff;
    }   Data;

typedef long FredId;

typedef struct {
    FredEventType     typeE;
    FredId            replyTo;
    DataType          typeD;
    Data              buffer;
    }  fredscEvent ;  /* event control structure */

typedef int Event_Handle;


/* prototypes of public functions  */

/*#*******************************************************************
**     Event_Handle  fredInit ( char *cmd_line[], int fd1, int fd2 )
**
**        *cmd_line[] is a NULL-terminated list of strings, this is 
**                    the command line to exec the child
**        fd1, fd2    are the child's file descriptors to be used for
**                    communication
**         
**      This function is called from the parent only, for each forked
**      child just once. The Event_Handle identifies communication channels.
**      
**       Returns -1 when error
**
*********************************************************************/


/*#*******************************************************************
**     int  fredGetEvent ( handle, fredscEvent *fe )
**
**          This function receives an event from the channel 'handle'.
**          fe is the control structure for events: 
**            fe->typeE  holds the FredEventType
**            fe->typeD  holds the DataType
**            fe->buffer holds the location of data
**          You can fill some of the fields before calling fredGetEvent(),
**          but the function may as well override this. E.g. when receiving
**          only filename for data, then any fe->typeD is overwritten to
**          DATA_IN_FILE. In all other cases received data are written either
**          to given filedesc or memory buffer.
**
**          Returns 0 if there is no event available, 1 if succeeds, and
**          -1 when error.
**          
*********************************************************************/


/*#*******************************************************************
**     FredId  fredSendEvent ( Event_Handle handle, fredscEvent *fe )
**
**          This function sends an event to the channel 'handle'.
**          fe is the control structure for events: 
**            fe->typeE  holds the FredEventType
**            fe->typeD  holds the DataType
**            fe->buffer holds the location of data
**          You must fill all the fields of fe before calling fredGetEvent(),
**
**          Returns 1 if succeeds, and
**          -1 when error.
**          
*********************************************************************/


/*#*******************************************************************
**     void  fredKillEvent ( fredscEvent *fe )
**
**        Frees the memory allocated for the event.
**
*********************************************************************/


/*#*******************************************************************
**     void  fredCloseEvent_handle ( Event_Handel ha )
**
**        Closes the communication channel ha.
**
*********************************************************************/
#ifndef PROTO
#define PROTO(x)        x 
#define VPROTO(x)        x 
#define ___MAZ___
#endif

int           fredInit PROTO(( char *cmd_line[])) ;
int           my_dup2  PROTO(( int fd1, int fd2 ));
int           read_sc  PROTO(( int sd, char *buff, int size_buff));
FredId        fredSendEvent PROTO((fredscEvent *fe));
int           fredGetEvent PROTO((fredscEvent *fe)); 
int           buff_to_mem PROTO(( char **memptr, long *curr_size, 
                                  long *size_all, char *buff, long length));
void          fredKillEvent PROTO((fredscEvent *fe));
int           wait_to_events PROTO(( fredscEvent *fe, 
                                     FredEventType  set_event));


extern int   close  PROTO(( int fd));
extern int   socketpair PROTO(( int d, int type, int protocol, int sv[2] ));
extern int   dup2  PROTO(( int fd1, int fd2 ));
extern int   fork  PROTO((void));
extern void  _exit PROTO(( int e_stst));
extern int   execvp PROTO(( char *path, char *argv[]));
extern int   select PROTO(( int width, fd_set *r_s, 
                            fd_set *w_s, fd_set *x_s, 
                            struct timeval *timout ));
extern void  bzero PROTO(( char *buff, int length));
#ifndef DEF_MYSTDLIB_H
extern int   read  PROTO(( int sd, char *buff, int size_buff));
extern int   write PROTO(( int sd, char *buff, int size_buff));
extern void  printf VPROTO((char  *fmt, ...));
/* extern void * memcpy   (void *dest, const void *src, size_t n); */
char    * strcpy   (char *dest, const char *src);
size_t  strlen   (const char *s);
#endif


#ifdef ___MAZ___ 
#undef PROTO
#undef VPROTO
#undef ___MAZ___
#endif

#endif

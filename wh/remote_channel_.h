/*  File: remote_channel_.h
 *  channels interface : derived from the GO language channels 
 *  Author: Danielle and Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov) Mark Sienkiewicz (sienkiew@stsci.edu)
 *  Copyright (C) J Thierry-Mieg, 2015
 * -------------------------------------------------------------------
 * Created: mars 2015, mieg
 *-------------------------------------------------------------------
 */

#ifndef REMOTE_CHANNEL__H
#define REMOTE_CHANNEL__H

CHAN *uTaskCreateReaderChannel (TASK *task, const char *channelName, int max, int size) ;
CHAN *uTaskCreateWriterChannel (TASK *task, const char *channelName, int max, int size) ;
typedef enum { DONE=1, TIMEOUT } TASK_LIFE ;

#define taskCreateReaderChannel(_t,_nam,_max,_type) uTaskCreateReaderChannel((_t),(_nam),(_max),sizeof(_type))
#define taskCreateWriterChannel(_t,_nam,_max,_type) uTaskCreateWriterChannel((_t),(_nam),(_max),sizeof(_type))


#endif



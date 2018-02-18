/*  Last edited: Dec  3 17:27 1998 (fw) */

/* $Id: disk.h,v 1.3 2015/09/18 22:13:50 mieg Exp $ */

                     /* disk.h                             */
                     /* public functions of disksubs.c     */
                     /* handling the disk.                 */

#ifndef DEFINE_DISK_H
#define DISK_H

#ifndef DEF_BP
#define DEF_BP
typedef void* BP ;
#endif

#ifndef DEF_DISK
#define DEF_DISK
#if defined(ACEDB4)
typedef KEY  DISK;             /*holds disk addresses as 4 bytes unsined int */
#define NULLDISK 0
#else
typedef long  DISK;             /*holds disk addresses as file offsets*/
#define NULLDISK 0L
#endif /* ACEDB4 */
#endif /* DEF_DISK */
 
 
void diskPrepare    (void) ;
                      /*To be called once in the disk lifetime*/
                      /*prepares DISKSIZE blocks of storage and the*/
                      /*corresponding BAT file */
void diskSave       (void) ; 
                      /* Rewrites the Block allocation table BAT*/
                      /* of the main database file */
void diskWriteSuperBlock (BP bp) ;
                      /* General database management */
void diskavail(unsigned long int *free, unsigned long int *used, unsigned long int *plus, unsigned long int *minus, unsigned long int *dr, unsigned long int* dw) ;
                      /*gives disk statistics*/
void diskalloc      (BP bp);
                      /*modifies the BAT and gives a free block address*/
                      /* in sequence */
void diskfree       (BP bp) ;
                      /* returns d to the BAT */
void diskblockread  (BP p,DISK d) ; /* reads a block from disk */
void diskblockwrite (BP p) ; /* writes a block to disk */
BOOL dataBaseCreate (void) ;
void dataBaseAccessInit (void) ;
void diskFlushLocalCache(void);
BOOL saveAll(void) ; 

void diskFuseOedipeBats(Array fPlusArray, Array fMinusArray)  ;
BOOL diskFuseBats(KEY fPlus, KEY fMinus, 
		  KEY sPlus, KEY sMinus,
		  int *nPlusp, int *nMinusp) ;

void disknewStatus (void) ;
void diskShutDown (void) ;
#endif


 

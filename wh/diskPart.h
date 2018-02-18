/*  Last edited: Dec  1 20:33 1995 (mieg) */
 
  /* description  des partitions composants l'espaces des blocs */

/* $Id: diskPart.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

/* partitions description must reside on a unique disk block : */
/*#define MAX_PARTITIONS (BLOC_SIZE - (sizeof(int)*2)) / sizeof(PARTITION)*/
#define MAX_PARTITIONS 20

/* ATTENTION :
 * only a  multiple of INITIALDISKSIZE blocks can be used in each partition
 * See wh/disk.h , wspec/cachesize.wrm (DISK option)
 */
/* default maimum DataBase size in blocks */
#define MAX_DATABASE_SIZE 100 * 1024 

#define UNIX_FILE_SYSTEM 1
#define RAW_FILE_SYSTEM  2
#define LG_NAME          128

typedef struct {
    int type;             /* unix file system, raw, distant, ...            */
    char hostname[LG_NAME]; /* machine on which is the partition            */
    char fileSystem[LG_NAME];/* partition file system or "ACEDN if relative */
    char fileName[LG_NAME];  /* file name or "NONE" if not a file system    */
    int  maxBlocks;       /* max. size in n x  blocks                       */
    int  currentBlocks;   /* current size in n x blocks                     */
    int offset;           /* number of first non reserved block             */
    } PARTITION;

typedef struct {
    int       lastPartition; /* (i.e. number of active partitions)          */
    int       nbPartitions;  /* partition number                            */
    PARTITION partitions[MAX_PARTITIONS];
    } BASE_DEF;  


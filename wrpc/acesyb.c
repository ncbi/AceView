/*  $Id: acesyb.c,v 1.1.1.1 2002/07/19 20:23:32 sienkiew Exp $  */
/*  File: acesyb.c
 #  Author: Detlef Wolf (D.Wolf@dkfz-heidelberg.de)
 #  Copyright (C) European Data Resource, 1994
 #-------------------------------------------------------------------
 # This file is part of the IGD project
 #
 # Description: sybase storage for ACEDB B-classes
 # Main data structure: array of nodes (representing an ACEDB tree object)
 #
 # HISTORY:
 # Last edited: Tue Jul  5 04:18:24 WET DST 1994  detlef
# ***** test version, has lots of checks, since I had some
        memory allocation anomalities
        checks can be removed if running stable  *******
Mon Jul 11 08:26:24 WET DST 1994  detlef: add reaper, limit cachesize
 * * Jun 14 16:52 1994 (dok256): still to do:
data type AY_COMMENT handling
 # Created: Mon Jun  6 13:42:51 1994 (dok256)
 #-------------------------------------------------------------------
 */

#include <stdio.h>

#include "sybfront.h"
/* original sybfront.h declared BOOL, had to be replaced */

#include <sybdb.h>

#include "acesyb.h"


/*------------------------ acesyb module memore / sybase specific --------- */
static RETCODE        result_code;
static DBPROCESS     *dbproc;       /* Our connection with SQL Server. */
static LOGINREC      *login;        /* Our login information. */

/*------------------------------- acesyb module memory / general --------- */
static Array ayNodeBuf=0;  /* sybase cache, Array of ayNodeType */
static int aySize=0;       /* size of ayNodeBuf in bytes */
static int ayReapPos=0;    /* start point for cache reaper search */

static Array ayDeleteSet=0;  /* allocated once globally to save time */
static Array ayUpdateSet=0;
static Array ayUpdateSet_buf=0;
static Array ayInsertSet=0;
#define ayCmpLen (sizeof(KEY)+sizeof(int))
static char *typeNames[]={"char","int","float","bool","key","date"};
#define      typeNameCnt 6

static char *debugbuf;
static int ayTrace=0;  /* Tracelevel, 0=no trace */
static int ayMaxSize=0;  /* Maximum permanent cache size,
                            however temporarily it can be more */
#define sqlMax 4096 
                     /* maximum length of an SQL command */
#define AY_MAXSIZE 65536
                     /* default size of sybase cache */
/*------------------------------------- local sybase functions --------- */
#define ERR_CH     stderr
#define OUT_CH     stdout

int ay_err_handler(dbproc, severity, dberr, oserr, dberrstr, oserrstr)
DBPROCESS       *dbproc;
int             severity;
int             dberr;
int             oserr;
char            *dberrstr;
char            *oserrstr;
{
	if ((dbproc == NULL) || (DBDEAD(dbproc)))
		return(INT_EXIT);
	else 
	{
		fprintf (ERR_CH, "DB-Library error:\n\t%s\n", dberrstr);

		if (oserr != DBNOERR)
		fprintf (ERR_CH, "Operating-system error:\n\t%s\n", oserrstr);

		return(INT_CANCEL);
	}
}

int ay_msg_handler(dbproc, msgno, msgstate, severity, msgtext, 
                srvname, procname, line)

DBPROCESS       *dbproc;
DBINT           msgno;
int             msgstate;
int             severity;
char            *msgtext;
char            *srvname;
char            *procname;
DBUSMALLINT     line;

{
	fprintf (ERR_CH, "Msg %ld, Level %d, State %d\n", 
	        msgno, severity, msgstate);

	if (strlen(srvname) > 0)
		fprintf (ERR_CH, "Server '%s', ", srvname);
	if (strlen(procname) > 0)
		fprintf (ERR_CH, "Procedure '%s', ", procname);
	if (line > 0)
		fprintf (ERR_CH, "Line %d", line);

	fprintf (ERR_CH, "\t%s\n", msgtext);
        if (severity>0) messcrash("SYBASE error %s",msgtext);
	return(0);
}

/*------------------------------------- debugging aids --------- */

void printNodeHeading() {
  printf("       key   nno right  down |class tag pos| type  data\n");
}

void printNode(ayNodeType *n) {
  printf("%c%9d %5d %5d %5d | %3d %3d %3d | %1d ",
     n->temperature, n->key, n->thisNode, n->rightNode, n->downNode,
     n->classNo, n->tagNo, n->valPos, n->dataType);
  switch (n->dataType) {
    case AY_EMPTY: printf("AY_EMPTY\n"); break;
    case AY_INT  : printf("%d\n",n->data.intData); break;
    case AY_FLOAT: printf("%f\n",n->data.floatData); break;
    case AY_TEXT : printf("%d %s\n",(int)(n->data.textData), n->data.textData); break;
    case AY_BOOL : printf("%d\n",n->data.boolData); break;
    case AY_KEY  : printf("%d\n",n->data.keyData); break;
    case AY_DATE : printf("%d\n",n->data.dateData); break;
    default : printf("invalid data type\n");
  };
}

void printNodeArray(char *title, Array a) {
  int i;
  if (a) {
    printf("Array %s of ayNodeType with %d entries:\n", title, arrayMax(a));
    printNodeHeading();
    for (i=0; i<arrayMax(a); ++i)
      printNode(arrp(a,i,ayNodeType));
  }
  else
    printf("Array %s of ayNodeType is NULL");
}

int nodeSize(ayNodeType *np);
void checkSize(char *title) {
  /* check consistency of aySize with the ayNodeBuf */
  register int i; 
  register int size=0;
  for (i=0; i<arrayMax(ayNodeBuf); i++) 
    size+=nodeSize(arrp(ayNodeBuf,i,ayNodeType));
  if (size!=aySize) {
    printf("%s: checkSize: mismatch aySize=%d, size=%d\n", title, aySize, size);
    printNodeArray("ayNodeBuf", ayNodeBuf);
    messcrash("checkSize");
    aySize=size;  /* correct it */
  }
  else
    if (ayTrace>1) printf("%s: checkSize: ok, aySize=%d\n", title, aySize);
}

/*------------------------------------- local functions --------- */

/*-------------------------------------------------------------------*/
int nodeSize(ayNodeType *np) {
  /* determine the size of a node */
  if (np->dataType==AY_TEXT)
    return(sizeof(ayNodeType)+strlen(np->data.textData));
  return(sizeof(ayNodeType));
};

/*-------------------------------------------------------------------*/
static void sqlExec() {
  /* execute current SQL buffer that does not return rows */
  if (ayTrace>1) {
    debugbuf=messalloc(dbstrlen(dbproc)+1);
    dbstrcpy(dbproc, 0, dbstrlen(dbproc), debugbuf);
    printf("dbcmd=\"%s\"\n",debugbuf);
    messfree(debugbuf);
  };

  dbsqlexec(dbproc);
  while ((result_code = dbresults(dbproc)) != NO_MORE_RESULTS) {
    if (result_code == SUCCEED) {
      if (ayTrace>1) 
        printf("sqlExec: dbcount=%d\n", DBCOUNT(dbproc));
    };
  };
}

/*-------------------------------------------------------------------*/
static void dbAddData(DBPROCESS *dbproc, ayNodeType *np) {
  /* add the data from node np to the SQL command buffer */
  if (np->dataType>=typeNameCnt)
    fprintf(stderr,"dbAddData: datatype %d not supported\n",np->dataType);
  switch (np->dataType) {
   case AY_INT  : dbfcmd(dbproc, "%d",np->data.intData); break;
   case AY_FLOAT: dbfcmd(dbproc, "%f",np->data.floatData); break;
   case AY_TEXT : dbfcmd(dbproc, "\"%.252s\"",np->data.textData); break;
   case AY_BOOL : dbfcmd(dbproc, "%d",np->data.boolData); break;
   case AY_KEY  : dbfcmd(dbproc, "%d",np->data.keyData); break;
   case AY_DATE : dbfcmd(dbproc, "%d",np->data.dateData); break;
   default : fprintf(stderr,"dbAddData: invalid data type %d\n",np->dataType);
   };
}

/*-------------------------------------------------------------------*/
static int isEqual(ayNodeType *a, ayNodeType *b) {
  /* return 0 if contents of both nodes a and b are equivalent */
  if (memcmp((char *)a, (char *)b, sizeof(ayNodeType))==0) return 0;
  if (a->dataType!=AY_TEXT || b->dataType!=AY_TEXT
   || !(a->data.textData) || !(b->data.textData)) return 1; 
  if (a->key == b->key &&
      a->thisNode==b->thisNode &&
      a->rightNode==b->rightNode &&
      a->downNode==b->downNode &&
      a->classNo==b->classNo &&
      a->tagNo==b->tagNo &&
      a->valPos==b->valPos &&
      (strcmp(a->data.textData, b->data.textData)==0))
    return 0;
  return 1;
}

/*-------------------------------------------------------------------*/
static int ayOrder(ayNodeType *a, ayNodeType *b) {
  /* compare two nodes in the key and thisNode fields,
     return 0, if equal; -1 if a<b; +1 if a>b
  */
  if (a->key==b->key)
    if (a->thisNode==b->thisNode) return 0;
    else return (a->thisNode<b->thisNode)? -1: 1;
  else return (a->key<b->key)? -1: 1;
}

/*-------------------------------------------------------------------*/
static int intOrder(int *a, int *b) {
  /* compare two integers,
     return 0, if equal; -1 if a<b; +1 if a>b
  */
  if (*a==*b) return 0;
  if (*a<*b) return -1;
  return 1;
}

static void printDeleteSet() {
  register int i;
  printf("deleteSet(): %d, deleteSet->ayNodeBuf:", arrayMax(ayNodeBuf));
  for (i=0; i<arrayMax(ayDeleteSet); i++) 
    printf(" %d",arr(ayDeleteSet,i,int));
  printf("\n");
}

/*-------------------------------------------------------------------*/
static void deleteSet() {
  /* deletes the nodes pointed to by ayDeleteSet in ayNodeBuf 
     decrements aySize
     precondition: ayDeleteSet is in ascending order
     algo: gather the pieces that remain after the deletion
           shift them left so that there are no gaps;
           thus an element is at most moved once
  */
  int i;
  int len;
  int topos;    /* next position to be filled */
  int delpos;   /* index of element to deleted */
  int nextpos;  /* index of next elment to be deleted */
  ayNodeType *np;
  if (!ayDeleteSet || arrayMax(ayDeleteSet)<=0) return;
  if (ayTrace>2) printNodeArray("before delete:ayNodeBuf",ayNodeBuf);
  if (ayTrace>2) printDeleteSet();

  /* check that ayDeleteSet is in ascending order and has no dups */
  delpos=-1;  /* last position seen, misuse delpos */
  for (i=0; i<arrayMax(ayDeleteSet); i++) {
    if (delpos>=arr(ayDeleteSet,i,int))
      messcrash("internal error: invalid ayDeleteSet");
    delpos=arr(ayDeleteSet,i,int);
  };
  if (arr(ayDeleteSet,arrayMax(ayDeleteSet)-1,int)>=arrayMax(ayNodeBuf))
    messcrash("deleteSet out of bounds");
  if (arrp(ayNodeBuf, arr(ayDeleteSet,0,int), ayNodeType)->thisNode)
    messcrash("deleteSet does not start with thisNode==0");

  checkSize("deleteSet start");

  topos=arrayMax(ayNodeBuf)-arrayMax(ayDeleteSet);
  for (i=arrayMax(ayDeleteSet)-1; i>=0; i--) {
    delpos=arr(ayDeleteSet,i,int);   /* get position in ayNodeBuf */
    np=arrp(ayNodeBuf,delpos,ayNodeType);
    /* if the element contained string data, delete it */
    aySize -= nodeSize(np); 
    if (np->dataType==AY_TEXT)
      messfree(np->data.textData);
    arrayRemove(ayNodeBuf,np,ayOrder);
  };
  if (arrayMax(ayNodeBuf)!=topos)
    messcrash("unexpected number of elements");

  if (ayTrace>2) printNodeArray("after delete:ayNodeBuf",ayNodeBuf);
  checkSize("deleteSet start");
}


/*-------------------------------------------------------------------*/
static void reaper() {
  /*  check, if buffer is over the memory limit and if yes, then
      bring it down to 80% full
  Precondition: ayReapPos -- next location to look at (memorized across calls)
                ayMaxSize -- ideal buffer size
                aySize    -- current Buffer size
                ayNodeBuf -- the cache
  Postcondition: ayNodeBuf's size is less than ayMaxSize
  */
  if (aySize < ayMaxSize) return;  /* nothing to do */ 
  {
  KEY curKey;
  ayNodeType *np;
  int delete;    /* flag for delete mode */
  int toReap=(aySize-ayMaxSize) + ((float) ayMaxSize * .2);  
                            /* number of bytest that should be freed */
  arrayMax(ayDeleteSet)=0; /* index into ayNodeBuf */

  if (ayTrace>2) printf("reaper wants %d bytes out of %d.\n", toReap, aySize);

  /* a few plausibility checks on that delicate structure: */
  if (!arrayMax(ayNodeBuf))
    messcrash("reaper finds aySize and ayNodeBuf inconsistent");
  checkSize("reaper");
  /* now begin real work */

  /* position on start of object */
    /* position on the first node of next object, i.e. a node
       where thisNode==0, and remember that node in np */
  if (ayReapPos>=arrayMax(ayNodeBuf)) ayReapPos=0;
  while (ayReapPos<arrayMax(ayNodeBuf)
     && (np=arrp(ayNodeBuf, ayReapPos, ayNodeType))->thisNode)
     ayReapPos++;
  if (ayReapPos>=arrayMax(ayNodeBuf)) ayReapPos=0;
  if (ayTrace>1) printf("reaper starts at node %d\n", ayReapPos);
  if (ayTrace>1) printNodeArray("reaper start ayNodeBuf", ayNodeBuf);

  if (arrp(ayNodeBuf, ayReapPos, ayNodeType)->thisNode)
    messcrash("ayReapPos not at start of object");

  while (toReap>0) {
    /* pointer is on start of object */
    np=arrp(ayNodeBuf, ayReapPos, ayNodeType);
    if (np->thisNode) 
      messcrash("reaper: thisNode!=0");

    /* hot nodes get cold, cold nodes get reaped 
       since the reaper gathers here and there it eventually
       visits a node twice. It must rembers it (via AY_DEL)
       in order not to reap a node twice */
    delete=0;
    if (np->temperature==AY_HOT)  
      np->temperature=AY_COLD;
    else if (np->temperature!=AY_DEL) 
      { np->temperature=AY_DEL;  delete=1; };

    /* skip over complete object and mark it on request */
    curKey=np->key;
    do {
        if (delete) {
/* printf("reap add ayReapPos %d\n", ayReapPos);  */
          toReap-=nodeSize(arrp(ayNodeBuf, ayReapPos, ayNodeType));
          array(ayDeleteSet, arrayMax(ayDeleteSet), int)=ayReapPos;
        };
        ayReapPos++;
      } while (ayReapPos<arrayMax(ayNodeBuf) &&
               arrp(ayNodeBuf, ayReapPos, ayNodeType)->key==curKey);

    if (ayReapPos>=arrayMax(ayNodeBuf)) ayReapPos=0;

  };  /* while toReap */
  if (ayTrace>1) printf("reaper ends at node %d\n", ayReapPos);

  /* deleteSet requires an ordered set */
  arraySort(ayDeleteSet, intOrder);
  deleteSet();    /* deletes nodes from ayDeleteSet in ayNodeBuf */
  if (ayTrace>2) printf("reaper left %d bytes alive\n", aySize);
  };
}  

/*-------------------------------------------------------------------*/
static void dupString(ayNodeType *from, ayNodeType *to) {
  /* duplicated the string part of an element, if it exists */
  if (from->dataType==AY_TEXT) {
    to->data.textData=  messalloc(strlen(from->data.textData));
    strcpy(to->data.textData, from->data.textData);
  };
};

/*-------------------------------------------------------------------*/
static void StringSteal(ayNodeType *from, ayNodeType *to) {
  /* moved the text of node from into node to and zeros the text in from,
     if it exists 
  */
  if (from->dataType==AY_TEXT) {
    to->data.textData=  from->data.textData;
    from->data.textData=NULL;
  };
};

/*-------------------------------------------------------------------*/
static int readSybase(KEY key) {
  /* copy object with key key from SYBASE database into ayNodeBuf, if it exists
     returns: 1 if object read, 0 if it does not exist
     postcondition: object is in ayNodeBuf
  */
  int nodeCount, i;
  ayNodeType rn;          /* result node */
  int pos;
  ayNodeType queryNode;   /* dummy to construct a search pattern */
  queryNode.key=key;
  queryNode.thisNode=0;
  if (arrayFind(ayNodeBuf, &queryNode, &pos, ayOrder)) {
    arrp(ayNodeBuf, pos, ayNodeType)->temperature=AY_HOT;
    if (ayTrace>1)
      printf("readSybase: buffer hit on key %d\n", key);
    return 1;
  };

  /* object was not in buffer, so try and get it from sybase */
  nodeCount=0;  /* numer of nodes found */
  /* retrieve for all data types */
  for (i=0; i<typeNameCnt; i++)
    dbfcmd(dbproc,"select key, thisNode, rightNode, downNode, classNo, tagNo, valPos, value from nodes_%s where key= %d",typeNames[i],key);
  dbsqlexec(dbproc);

  if (ayTrace>4) {
    debugbuf=messalloc(dbstrlen(dbproc)+1);
    dbstrcpy(dbproc, 0, dbstrlen(dbproc), debugbuf);
    printf("dbcmd=\"%s\"\n",debugbuf);
    messfree(debugbuf); 
  };

  while ((result_code = dbresults(dbproc)) != NO_MORE_RESULTS) {
    if (result_code == SUCCEED) {
      while (dbnextrow(dbproc) != NO_MORE_ROWS) {

  /*  1    2         3          4         5        6      7        8   */
  /* key, thisNode, rightNode, downNode, classNo, tagNo, valPos, value */
        /* transfer result from sybase buffer into a node structure */
        if (key!=(*((DBINT *)dbdata(dbproc, 1))))
          messcrash("ayRead: received wrong key %d\n",
                  (*((DBINT *)dbdata(dbproc, 1))));
        rn.key=(*((DBINT *)dbdata(dbproc, 1)));
        dbconvert(dbproc, SYBINT2, dbdata(dbproc,2), (DBINT)-1,
                          SYBINT4, &rn.thisNode    , (DBINT)-1);
        dbconvert(dbproc, SYBINT2, dbdata(dbproc,3), (DBINT)-1,
                          SYBINT4, &rn.rightNode   , (DBINT)-1);
        dbconvert(dbproc, SYBINT2, dbdata(dbproc,4), (DBINT)-1,
                          SYBINT4, &rn.downNode    , (DBINT)-1);
        dbconvert(dbproc, SYBINT1, dbdata(dbproc,5), (DBINT)-1,
                          SYBINT4, &rn.classNo     , (DBINT)-1);
        dbconvert(dbproc, SYBINT2, dbdata(dbproc,6), (DBINT)-1,
                          SYBINT4, &rn.tagNo       , (DBINT)-1);
        dbconvert(dbproc, SYBINT2, dbdata(dbproc,7), (DBINT)-1,
                          SYBINT4, &rn.valPos      , (DBINT)-1);

        switch (DBCURCMD(dbproc)) {
        /* transfer the value part, according to type */
        /* the case numbers correspond to the order in typeNames[] */
        /* "char","int","float","bool","key","date" and */
        case 1: rn.dataType=AY_TEXT;
                rn.data.textData= 
                         messalloc(dbdatlen(dbproc, 8)+1); /* for \0 */
                dbconvert(dbproc,   /* append \0 to terminate string */
                          SYBCHAR, dbdata(dbproc,8), dbdatlen(dbproc,8),
                          SYBCHAR, rn.data.textData, (DBINT)-1);
                break;
        case 2: rn.dataType=AY_INT;
                rn.data.intData=(*((DBINT *)dbdata(dbproc, 8)));
                break;
        case 3: rn.dataType=AY_FLOAT;
                /* hopefully the compiler is clever enough to
                   convert that double (DBFLT8) into float (4 Bytes) */
                rn.data.floatData=(*((DBFLT8 *)dbdata(dbproc, 8)));
                break; 
        case 4: rn.dataType=AY_BOOL;
/* turned into smallint                    dbconvert(dbproc,
                          SYBBIT,  dbdata(dbproc,8), (DBINT)-1, 
                          SYBINT4, rn.data.boolData, (DBINT)-1); */
                dbconvert(dbproc,
                          SYBINT2, dbdata(dbproc,8), (DBINT)-1,
                          SYBINT4, &rn.data.boolData, (DBINT)-1);
                break;

        case 5: rn.dataType=AY_KEY;
                rn.data.keyData=(*((DBINT *)dbdata(dbproc, 8)));
                break;
        case 6: rn.dataType=AY_DATE;
                /* should have a fancy conversion, but not now */
                rn.data.dateData=(*((DBINT *)dbdata(dbproc, 8)));
                break;
        default: 
            messcrash("ayRead: DBCURCMD(dbproc)=%d", DBCURCMD(dbproc));
	    }; /* switch */

        /* copy into buffer */
/* printf("readSybase: insert key %d thisNode %d\n", rn.key, rn.thisNode); */
        rn.temperature=AY_HOT;
        if (!arrayInsert(ayNodeBuf, &rn, ayOrder))
          fprintf(stderr,"ayRead: key=%d, thisNode=%d already in buffer
inconsistent database?\n", rn.key, rn.thisNode);
        aySize += nodeSize(&rn); 
        nodeCount++;
        }; /* loop over the rows */
      }; /* command succeeded */
    }; /* loop over the batch */
    if (arrp(ayNodeBuf,0,ayNodeType)->thisNode)
      messcrash("sybaseRead: (0)->thisNode exists");
    return (nodeCount>0) ? 1 : 0;
}

/*-------------------------------------------------------------------*/
static int copy1Object(Array from, Array to, KEY key) {
  /* copy nodes with key from from to to
     inputs: from is sorted by key
             key  -- key of object to copy
             to   -- must already exist as an array
     output: to   -- contains exactly the object identified by key
                     or is empty, if key not found
             returns number of nodes copied
  */
    int topos=0;
    int frompos;
    ayNodeType queryNode;   /* dummy to construct a search pattern */
    queryNode.key=key;
    queryNode.thisNode=0;

    if (!to)
      messcrash("copy1Object: to-Array must already exist");
    ayFreeNodeBuf(to);
    setCompareLength(from, ayCmpLen);  /* in arrayFind */
    if (arrayFind(from, &queryNode, &frompos, ayOrder))
      do {
        array(to, topos, ayNodeType) = 
                                    arr(from, frompos, ayNodeType);
        dupString(arrp(from,frompos,ayNodeType),
                  arrp(to, topos, ayNodeType));     /* copy also the string */
        topos++; frompos++;
      } while (frompos<arrayMax(from) && 
               (arrp(from, frompos, ayNodeType))->key==key);
    if (ayTrace>1) printf("copy1Object: %d nodes copied.\n",topos);
    return topos;
}

/* ----------------------------- exported functions ------------------ */

/*-------------------------------------------------------------------*/
void ayFreeNodeBuf(Array a) {
  /* clears array a with elements of type ayNodeType,
     freeing also memory allocated for char data
  */
  register int i;
  register ayNodeType *np;
  if (!a) return;  /* nothing to free */
  for (i=0; i<arrayMax(a); i++) {
    np=arrp(a,i,ayNodeType);
    if (np->dataType==AY_TEXT)  /* if AY_TEXT, the must exist */
      messfree(np->data.textData);
  };
  arrayMax(a)=0;
}

/*-------------------------------------------------------------------*/
void ayBegin() {
  /* initialize module and begin a sybase session
  precondition: environment variable AY_USER contains the sybase user name
                environment variable AY_PW contains the sybase password
                environment variable DSQUERY contains the sybase server name
                the default database of AY_USER is the database to operate on
                environment variable AY_TRACE may contain trace level
                                              (the higher the more output)
                environment variable AY_CACHE may contain chachesize in bytes
  postcondition: ayEnd(), ayRead(), ayWrite() can be called
  */
  char *AY_USER;
  char *AY_PW;

  /* tracing requested? */
  if (getenv("AY_TRACE")) {
    sscanf(getenv("AY_TRACE"),"%d",&ayTrace);
    printf("ayBegin: trace level set to %d\n", ayTrace);
  };

  /* cachesize given by environment? */
  if (getenv("AY_CACHE")) {
    sscanf(getenv("AY_CACHE"),"%d",&ayMaxSize);
    printf("ayBegin: cachesize set to %d\n", ayMaxSize);
  } else ayMaxSize=AY_MAXSIZE;

  /* initialize general part of module memory */
  ayNodeBuf=arrayReCreate(ayNodeBuf, 100, ayNodeType);
  setCompareLength(ayNodeBuf, ayCmpLen);  /* in arrayFind */

  /* initialize the sybase part of module memory */
  if (dbinit() == FAIL)  exit(ERREXIT);
		
  /* Install error-handling and message-handling routines */	
  dberrhandle(ay_err_handler);
  dbmsghandle(ay_msg_handler);

  /* setup login info */
  /* getenv is declared in mystdlib.h */
  if (!(AY_USER=getenv("AY_USER")))
    messcrash("please set environment variable AY_USER to the sybase user's name");
  if (!(AY_PW=getenv("AY_PW")))
    messcrash("please set environment variable AY_PW to the sybase user's password");
  login = dblogin();
  DBSETLUSER(login, AY_USER);
  DBSETLPWD(login, AY_PW);
  DBSETLAPP(login, "acesyb");

  /* login to the sybase server: open a database process */
  dbproc = dbopen(login, NULL);

  dbcmd(dbproc, "begin transaction acesyb");
  dbsqlexec(dbproc);

  /* set up these arrays (only once) */
  ayDeleteSet=arrayReCreate(ayDeleteSet,100,int); /* index into ayNodeBuf */
  ayUpdateSet=arrayReCreate(ayUpdateSet,100,int); /* index into wnodes */
  ayUpdateSet_buf=arrayReCreate(ayUpdateSet_buf,100,int); /* index into ayNodeBuf*/
  ayInsertSet=arrayReCreate(ayInsertSet,100,int); /* index into wnodes */

};  /* begin a session */

/*-------------------------------------------------------------------*/
void ayEnd(BOOL commit) {
  /* End a session;
     precondition: ayBegin has been called
     input: commit=TRUE    -- changes are made permanent
            commit==FALSE  -- changes are undone
  */
  if (commit)
    dbcmd(dbproc, "commit transaction acesyb");
  else
    dbcmd(dbproc, "rollback transaction acesyb");
  printf("ayEnd(commit=%d) sql_retcode=%d (ok=1)\n",commit,dbsqlexec(dbproc));
  ayFreeNodeBuf(ayNodeBuf);
  aySize=0;
};    /* end current session and store results */

/*-------------------------------------------------------------------*/
void ayRead(Array ayNodes, ayQueryType *query) {
   /* retrieve nodes given by query, return result in ayNodes,

output: ayNodes -- complete retrieved object(s)
        ayNodes has an independent copy of the nodes, the caller
        of ayRead can modify also the strings included in some
        nodes as s/he likes.
returns -- if call ok, else crash;
input: a query,
possible queries:
1. with exact key:  query->key!=NULL
   retrieve all nodes of the object with key
2. with query->key==NULL
   and  query->dataType==AY_EMTPY:
   search for whole classes and existing tags
   2.1 query->classNo==NULL and query->tagNo==NULL
   2.2 query->classNo==NULL and query->tagNo!=NULL
   2.3 query->classNo==NULL and query->tagNo==NULL
   2.4 query->classNo!=NULL and query->tagNo!=NULL
3. with query->key==NULL
   and  query->dataType!=AY_EMTPY
   and  query->tagNo!=NULL
   and  query->valPos!=NULL
   and  query->data != NULL
   search for a specific value in the given position
   (convention: tagNo=NULL, valPos=NULL searches the object names)

   precondition: ayBegin() was successfully called
   postcondition:
     ayNodeBuf contains the retrieved data
*/


if (!ayNodes)
  messcrash("ayRead: ayNodes must already exist");
ayFreeNodeBuf(ayNodes);  /* prepare for receiving */

if (query->mask.key!=NULL) {  /* an exact key query */

  if (ayTrace>1) printf("ayRead: key=%d\n", query->mask.key);
  if (ayTrace>3) printNodeArray("at start of ayRead, ayNodeBuf", ayNodeBuf);
    
  if (!readSybase(query->mask.key)) return;

  /* must be in buffer now, so copy1Object must succeed */
  if (!copy1Object(ayNodeBuf, ayNodes, query->mask.key))
    messcrash("ayRead: inconsistent Buffer, should contain object\n");

  reaper();  /* trim cache */
  if (ayTrace>1) printNodeArray("ayRead result ayNodes", ayNodes);
  if (ayTrace>3) printNodeArray("after ayRead: ayNodeBuf", ayNodeBuf);
  return;
} /* exact key query */
else {
  printf("non-exact queries not yet treated.\n");
};
messcrash("ayRead: problematic call");
} /* ayRead() */

/*-------------------------------------------------------------------*/
void ayWrite(Array wnodes) {
  /* write a set of Nodes 
  input: wnodes  -- set of nodes to write; contains either none or all 
                     nodes of an object; must be sorted by key and thisNode, 
                     i.e. all nodes of one key are contiguous
                 since I assume that the wnodes buffer is not needed
                 after the writes, I steal the text memory from it,
                 leaving zero text pointers in wnodes.
  precondition: ayNodeBuf may contain nodes read before --
                   if a node of an object is in ayNodeBuf, then
                   all nodes of that object are in ayNodeBuf
  postcondtion: the objects in wnodes are in ayNodeBuf and in the
                sybase database
  */
  
  /* algorithm:
  comparing the wnodes with ayNodeBuf yields three sets of nodes:
  nodes to delete/update/insert.
  these are translated into sybase commands.
  To save space, these sets are not implemented as Arrays of Nodes,
  but Arrays of integers which are interpreted as indices of nodes
  in ayNodeBuf.
  */

  KEY curKey;
  int i,j,p,q;
  int sqlWritten;   /* number bytes written into SQL command buffer */
  ayNodeType *npq;  /* node pointer in ayNodeBuf */
  ayNodeType *npp;  /* node pointer in wnodes*/

  if ( !wnodes || !arrayMax(wnodes) ) {
    fprintf(stderr,"ayWrite: nothing to do\n");
    return;
  };
  setCompareLength(wnodes, ayCmpLen);  /* for arrayFind */
  /* clear the indexing arrays */
  arrayMax(ayDeleteSet)=0; /* index into ayNodeBuf */
  arrayMax(ayUpdateSet)=0; /* index into wnodes */
  arrayMax(ayUpdateSet_buf)=0; /* index into ayNodeBuf */
  arrayMax(ayInsertSet)=0; /* index into wnodes */
  
  if (ayTrace>1) printNodeArray("ayWrite input: wnodes", wnodes);
  if (ayTrace>2) printNodeArray("ayWrite at start ayNodeBuf", ayNodeBuf);

  /* get current state from database into buffer */
  for (p=0; p<arrayMax(wnodes); p++) {
    npp=arrp(wnodes,p,ayNodeType);
    npp->temperature=AY_HOT;   /* mark all nodes from input as hot, so 
                                  when they get copied into
                           the cache, just these entries are considered hot */
    if (npp->thisNode) continue;
    readSybase(npp->key);
  };

  /* determine the delta sets,
     i.e. compare the wnodes and the ayNodeBuf to arrive at
     the minimal actions necessary to for SYBASE
     the algorithm relies on both arrays being sorted
     by (key, thisNode) */
  p=0; 
  while (p<arrayMax(wnodes)) {   /* loop over the input nodes */
    npp=arrp(wnodes,p,ayNodeType);
    curKey=npp->key;
    if (npp->thisNode !=0)
      messcrash("ayWrite: p=%d, first node of key %d is %d instead of 0",
                 p, npp->key, npp->thisNode);
    if (arrayFind(ayNodeBuf, npp, &q, ayOrder)) { /* buffered object */
      while (p<arrayMax(wnodes) && q<arrayMax(ayNodeBuf) 
             && (npp=arrp(wnodes,p,ayNodeType))->key==curKey
             && (npq=arrp(ayNodeBuf,q,ayNodeType))->key==curKey) {
        if (npp->thisNode == npq->thisNode) {
          if (isEqual(npp, npq)==0) {
            if (ayTrace>1) 
              printf("ayWrite: buffer hit on key=%d, thisNode=%d\n",
                   npp->key, npp->thisNode);
            /* identical data -> a buffer hit */
          }
          else { /* identical keys, but different data */
            /* remember both versions */
            array(ayUpdateSet, arrayMax(ayUpdateSet), int)=p;
            array(ayUpdateSet_buf, arrayMax(ayUpdateSet_buf), int)=q;
          };       
          p++; q++;
        }
        else if (npp->thisNode > npq->thisNode) {
          /* ayNodeBuf(npq) is not in wnodes */
          array(ayDeleteSet, arrayMax(ayDeleteSet), int)=q++;
        }
        else {
          /* wnodes(npp) is not in ayNodeBuf */
          array(ayInsertSet, arrayMax(ayInsertSet), int)=p++;
        };
      }; /* loop over the nodes of one object */
      /* might be that there are some more old lines in ayNodeBuf: */
      while ( q<arrayMax(ayNodeBuf) 
             && (arrp(ayNodeBuf,q,ayNodeType))->key==curKey) {
        array(ayDeleteSet, arrayMax(ayDeleteSet), int)=q++;
      };
      /* and remaining nodes in wnodes */
      while ( p<arrayMax(wnodes) 
             && (arrp(wnodes,p,ayNodeType))->key==curKey) {
        array(ayInsertSet, arrayMax(ayInsertSet), int)=p++;
      };
    } /* buffered object */
    else {  /* object from wnodes not in ayBufNode */
      /* mark all its nodes for insert */
      do 
        array(ayInsertSet, arrayMax(ayInsertSet), int)=p++;
      while (p<arrayMax(wnodes) &&
             (arrp(wnodes,p,ayNodeType))->key==curKey);
    };
  };   /* loop over the input nodes */

  if (ayTrace>1) {
    printf("updateSet->wnodes:");
    for (i=0; i<arrayMax(ayUpdateSet); i++) 
      printf(" %d",arr(ayUpdateSet,i,int));
    printf("\n");
    printf("updateSet_buf->ayNodeBuf:");
    for (i=0; i<arrayMax(ayUpdateSet_buf); i++) 
      printf(" %d",arr(ayUpdateSet_buf,i,int));
    printf("\n");
    printDeleteSet();
    printf("insertSet->wnodes:");
    for (i=0; i<arrayMax(ayInsertSet); i++) 
      printf(" %d",arr(ayInsertSet,i,int));
    printf("\n");
  };

  /* now that the sets are determined, convert them into SQL actions */
  /* the order updates-deletes-inserts is important because
     updates and deletes use indices on the ayNodeBuf that
     are invalidated by deletes */

  /* updates --------------------- */
  if (arrayMax(ayUpdateSet)) {
    int fromtype, totype;
    if (ayTrace>4) printNodeArray("ayWrite/before update:ayNodeBuf",ayNodeBuf);
    if (ayTrace>4) printNodeArray("ayWrite/before update:wnodes",wnodes);
    sqlWritten=0;       /* number bytes written into SQL command buffer */
    for (j=0; j<arrayMax(ayUpdateSet); j++) {

      /* get parameters of current node */
      p=arr(ayUpdateSet,j,int);              /* index in wnodes */
      q=arr(ayUpdateSet_buf,j,int);          /* index in ayNodeBuf */   
      npp=arrp(wnodes, p, ayNodeType);
      npq=arrp(ayNodeBuf, q, ayNodeType);
      /* check a basic assumption */
      if (npp->thisNode != npq->thisNode ||
          npp->key != npq->key)
        messcrash("ayWrite/update: inconsistent wnodes / ayNodeBuf");
      fromtype=npq->dataType;
      totype=npp->dataType;

      if (fromtype==totype) {
        /* table does not change, create an update statement */
        dbfcmd(dbproc,
          "Update nodes_%s set rightNode=%d, \
          downNode=%d, classNo=%d, tagNo=%d, valPos=%d, value=",
          typeNames[totype], npp->rightNode, npp->downNode,
          npp->classNo, npp->tagNo, npp->valPos);
        dbAddData(dbproc, npp);
        dbfcmd(dbproc," where key=%d and thisNode=%d ", npp->key, npp->thisNode);
      }  /* same type */
      else { /* different type  -> issue insert + delete */
        dbfcmd(dbproc,
          "Delete from nodes_%s where key=%d and thisNode=%d ", 
           typeNames[fromtype], npq->key, npq->thisNode);
        dbfcmd(dbproc,
          "insert into nodes_%s values (%d, %d, %d, %d, %d, %d, %d,",
          typeNames[totype], npp->key, npp->thisNode, npp->rightNode, 
          npp->downNode, npp->classNo, npp->tagNo, npp->valPos);
        dbAddData(dbproc, npp);
        dbfcmd(dbproc, ")");
      };

      /* copy the updates from wnodes to ayNodeBuf */
      aySize = aySize -nodeSize(npq) +nodeSize(npp);
      arr(ayNodeBuf,q,ayNodeType)=arr(wnodes,p,ayNodeType);
      StringSteal(npp, npq);   /* move string from wnodes into ayNodeBuf */

      sqlWritten+= nodeSize(npq)+80;  /* estimate the lenght of the fixed
                                         SQL part to be 80 bytes */ 
      if (sqlWritten > sqlMax) {   /* submitt buffer, if too large already */
         sqlExec(); sqlWritten=0;
      };

    }; /* loop over the Update nodes */
    if (sqlWritten) sqlExec();

    if (ayTrace>4) printNodeArray("ayWrite/after update:ayNodeBuf",ayNodeBuf);
    if (ayTrace>4) printNodeArray("ayWrite/after update:wnodes",wnodes);

  }; /* if anything to update */
  
  /* deletes ------------------- */
  sqlWritten=0;
  if (arrayMax(ayDeleteSet)) {
    for (j=0; j<arrayMax(ayDeleteSet); j++) {
      npp=arrp(ayNodeBuf, arr(ayDeleteSet,j,int), ayNodeType);
      i=npp->dataType;
      if (i>=typeNameCnt)
	fprintf(stderr,"ayWrite: datatype %d not supported\n",i);
      dbfcmd(dbproc,
      "delete from nodes_%s where key=%d and thisNode=%d ",
      typeNames[i], npp->key, npp->thisNode);

      sqlWritten+= nodeSize(npp)+40;
      if (sqlWritten>sqlMax) {   /* submitt buffer, if too large already */
        sqlExec(); sqlWritten=0;
      };
    }; /* loop over the delete nodes */
      
    if (sqlWritten) sqlExec();
    deleteSet();  /* from ayNodeBuf */

  } /* if anything to delete */

  /* inserts ------------------- */
  if (arrayMax(ayInsertSet)) {
    sqlWritten=0;
    for (j=0; j<arrayMax(ayInsertSet); j++) {
      npp=arrp(wnodes, arr(ayInsertSet,j,int), ayNodeType);
      i=npp->dataType;
      dbfcmd(dbproc,
      "insert into nodes_%s values (%d, %d, %d, %d, %d, %d, %d,",
      typeNames[i], npp->key, npp->thisNode, npp->rightNode, npp->downNode,
      npp->classNo, npp->tagNo, npp->valPos);
      dbAddData(dbproc, npp);
      dbfcmd(dbproc, ") ");

      sqlWritten+= nodeSize(npp)+80;
      if (sqlWritten>sqlMax) {   /* submitt buffer, if too large already */
         sqlExec(); sqlWritten=0;
      };
    }; /* loop over the insert nodes */
      
    if (sqlWritten) sqlExec();

    /* insert these nodes into ayNodeBuf */
    /* currently this is implemented very dumb,
       shifts large amounts of memory around */
    if (ayTrace>4) printNodeArray("write/before insert:ayNodeBuf",ayNodeBuf);
    if (ayTrace>4) printNodeArray("write/before insert:wnodes",wnodes);
    for (j=0; j<arrayMax(ayInsertSet); j++) {
      npp=arrp(wnodes, arr(ayInsertSet,j,int), ayNodeType);
      if (!arrayInsert(ayNodeBuf, npp, ayOrder))
        messcrash("writeNode: j=%d alread there",j);
      /* find the place where insert occured */
      if (!arrayFind(ayNodeBuf, npp, &i, ayOrder))
        messcrash("ayWrite/insert: arrayFind fails");
      aySize += nodeSize(npp); 
      StringSteal(npp, arrp(ayNodeBuf,i,ayNodeType));  /* also the string,
           this destroys npp */
    };
    if (ayTrace>4) printNodeArray("write/after insert:ayNodeBuf",ayNodeBuf);
    if (ayTrace>4) printNodeArray("write/after insert:wnodes",wnodes);
    /* this is an inefficent algorithm;
       better alternatives: 
         run the same merge algorithm again as for determining the set
      or fold the delete and the insert together
    */

  }; /* if ayInsertSet not emtpy */

  reaper();    /* truncate the cache */
  if (ayTrace>2) printNodeArray("ayWrite at return: wnodes", wnodes);
  if (ayTrace>2) printNodeArray("ayWrite at return: ayNodeBuf", ayNodeBuf);

  return;
};

/*-------------------------------------------------------------------*/
int ayDelete(Array ayNodes) { 
  /* delete objects whose keys are in ayNodes
  input: ayNodes, where at least node.key is filled
  postcondition: denoted objects are deleted from buffer and database
  */
  int i, p, q, rowcount;
  ayNodeType *npp;                            /* node pointer in ayNodes */
  KEY curKey;
  arrayMax(ayDeleteSet)=0;                 /* clear index into ayNodeBuf */
  if (!ayNodes || arrayMax(ayNodes)<=0) 
    /* nothing to delete */ return 1;

  if (ayTrace>2) printNodeArray("ayDelete input:", ayNodes);

  /* determine the delete set and delete in database */
  curKey=0;                                            /* no key seen yet */
  for (p=0; p<arrayMax(ayNodes); p++) {      /* loop over the input nodes */
    npp=arrp(ayNodes,p,ayNodeType);
    if (curKey != npp->key) {                               /* a new key */
      curKey = npp->key;

      /* locate all nodes of key in ayNodeBuf and mark them for deletion */
      npp->thisNode =0;                              /* just to be sure */
      if (arrayFind(ayNodeBuf, npp, &q, ayOrder)) {  /* buffered object */
        while ( q<arrayMax(ayNodeBuf) 
           &&   arrp(ayNodeBuf,q,ayNodeType)->key==curKey) {
          /* mark for deletion */
          array(ayDeleteSet, arrayMax(ayDeleteSet), int)=q++; 
        };
      };  /* buffered object */

      /* delete in database */
      rowcount=0;  
      dbfcmd(dbproc,"exec delete_key %d",curKey);  
	/* call stored proc that deletes */
      dbsqlexec(dbproc);
      while ((result_code = dbresults(dbproc)) != NO_MORE_RESULTS) {
	if (DBCURCMD(dbproc)!=1)               /* can only be one command */
	   messcrash("ayRead: DBCURCMD(dbproc)=%d",DBCURCMD(dbproc));
	if (result_code == SUCCEED) {         /* results must be consumed */
	  while (dbnextrow(dbproc) != NO_MORE_ROWS) {};
	  rowcount=DBCOUNT(dbproc);
	};
      };
      if (rowcount==0) {  
	fprintf(stderr,"ayDelete: key=%d not found in database\n", (int) curKey);
	return 1;
      };
      /* end of delete in database */

    };  /* a new key */
  }; /* loop over the nodes in ayNodes */

  if (ayTrace>2) printDeleteSet();

  deleteSet();   /* delete in ayNodeBuf */
  return 0;
};

/*-------------------------------------------------------------------*/
int ayCount(int classNo) {
  /* count object on class classNo */
  return 0;
}; 

/*-------------------------------------------------------------------*/



/* To compile this file, type
 *    gcc -o ts -Wall -DALPHA -I../wh ts.c -L../bin.ALPHA_4_GCC -lace -lfree -lm -lpthreads
 */

/*
 * This is a test program to test phreading in acedb
 * There is a multitude of printout to allow checking 
 * at every stage
 */


#include "ac.h"
#include "mytime.h"
#include "pthread.h"
#include "query.h"

#define MAXK 255
#define MAXN 10
#define MAXDIRTY 10

static int disk[MAXK] ;
static int lock[MAXK] ;
static int dirty[MAXK] ;
static Associator cache = 0 ;
static Array mutexArray = 0 ;
static BOOL  cleanUpRunning = FALSE ;
static int nDirty, nCached ;

static void* cleanCache (void *v) ;

static AC_DB DB = 0 ;
AC_OBJ ac_key2obj (AC_DB db, KEY key, BOOL fillIt, AC_HANDLE handle) ;

static BOOL getMutex (int nn)
{
  pthread_mutex_t *mm = 0 ;

  if (!mutexArray)
    mutexArray = arrayCreate (64, pthread_mutex_t) ;
  mm = arrayp (mutexArray, nn, pthread_mutex_t) ;
  if (!mm)
    pthread_mutex_init (mm, 0) ;
  pthread_mutex_lock (mm) ;
  
  return TRUE ;
}

static BOOL releaseMutex (int nn)
{
  pthread_mutex_t *mm = 0 ;

  if (!mutexArray)
    mutexArray = arrayCreate (64, pthread_mutex_t) ;
  mm = arrayp (mutexArray, nn, pthread_mutex_t) ;
  if (mm)
    pthread_mutex_unlock (mm) ;

  return TRUE ;
}

static int getFromDisk (KEY key)
{
  sleep ( ( 6 + (key % 5)) % 5) ;
  return disk[key] ;
}

static BOOL writeToDisk (KEY key, int data)
{
  sleep ( ( 6 + (key % 5)) % 5) ;
  disk[key]  = data ;
  if (key < 10)
    printf ("Writing key = %d data = %d \n", key, data) ; 
  return TRUE ;
}

static BOOL getFromCache (KEY key, int *datap)
{
  BOOL ok = FALSE ;
  void *v ;

  getMutex (1) ;
  ok = assFind (cache, assVoid(key), &v) ;
  *datap = assInt(v) ;
  releaseMutex(1) ;
  return ok ;
}


static BOOL writeToCache (KEY key, int data)
{
  BOOL ok = FALSE, old ;
  int ii ;
  void *v ;

  getMutex (1) ;
  old = assFind (cache, assVoid(key), &v) ;
  if (!old)
    nCached++ ;
  if ( !old || assInt(v) != data)
    {
      dirty[key] = 1 ; nDirty++ ;
      ok = assInsert (cache, assVoid(key), assVoid(data)) ;
    }
  if (!cleanUpRunning && nDirty > MAXDIRTY)
    { 
      pthread_t cleanUpThread ;
      ii = pthread_create(&cleanUpThread, 0, cleanCache, 0) ;
      if (ii)
	printf("Error %d  creating cleanUpThread\n", ii) ;
      else
	{
	  cleanUpRunning = TRUE ;
	  pthread_detach (cleanUpThread) ;
	}
    }
  releaseMutex(1) ;
    
  return ok ;
}

static void* cleanCache (void *v)
{
  KEY key ;
  int i = MAXK, data ;

  getMutex (1) ;
  while (i--)
    {
      if (dirty[i] && !lock[i])
	{
	  key = i + 1 ;
	  releaseMutex(1) ;
	  getMutex (2 + key) ; /* avoid 1 == controlMutex */
	  if (assFind (cache, assVoid(key), &v))
	    {
	      data = assInt(v) ;
	      writeToDisk (key, data) ;
	    }
	  dirty[i] = 0 ;
	  releaseMutex (2 + key) ;
	  getMutex (1) ;
	  lock[i] = 0 ;
	  nDirty-- ;
	}
    }
  cleanUpRunning = FALSE ;
  releaseMutex(1) ; 
  return 0 ;
}

static BOOL lockKey (KEY key)
{
  BOOL ok = FALSE ;
  int n ;

  getMutex (1) ;
  n = lock[key] ;

  if (!n) { ok = TRUE ; lock[key] = 1 ; }
  releaseMutex (1) ;
  return ok ;
}

static BOOL unLockKey (KEY key)
{
  BOOL ok = FALSE ;
  int n ;

  getMutex (1) ;
  n = lock[key] ;

  if (!n) messcrash ("bad call to unlock key") ;
  ok = TRUE ; lock[key] = 0 ; 
  releaseMutex (1) ;
  return ok ;
}

static int getObj (KEY key)
{
  int data = 0 ;

  getMutex (2 + key) ; /* avoid 1 == controlMutex */
  if (!getFromCache (key, &data))
    {
      data = getFromDisk (key) ;
      writeToCache (key, data) ;
    }
  releaseMutex (2 + key) ;
  lockKey (key) ;
  return data ;  
}

static void destroyObj (KEY key) 
{
  unLockKey (key) ;
}

static void saveObj (KEY key, int data) 
{
  getMutex (2 + key) ; /* avoid 1 == controlMutex */
  if (lock[key])
    {
      writeToCache (key, data) ;
      unLockKey (key) ;
    }
  else
    messcrash ("bad call to saveObj") ; 
  releaseMutex (2 + key) ;
}

static void* doWriteToDisk (void *vp)
{
  KEY key = assInt (vp)  ;
  int data = key + 9 ;

  writeToDisk (key, data) ;
  return assVoid(1) ;
}

static void collectWrite (Array aa)
{
  int i ;
  pthread_t *pp ;
  void *vp ;
  
  for (i = 1 ; i < arrayMax (aa) ; i++)
    {
      pp = arrayp (aa, i, pthread_t) ;
      if (!pthread_join (*pp, &vp)) 
	{ if (assInt(vp) != 1)
	  printf("bad value returned by joined thread %d \n", i) ;
	}
      else
	printf("cannot join phread %d\n", i) ;
      pthread_detach (*pp) ;

    }

  return ;
}

static void* doCheckData (void *vp)
{
  KEY  key = assInt (vp)  ;
  int data ;

  data = getObj (key) ;
  destroyObj (key) ;
  if (data != key + 9)
    printf ("error key = %d data = %d != %d\n",key, data, key+9) ;
  return assVoid(data) ;
}

static void collectCheck (Array aa)
{
  int i ;
  pthread_t *pp ;
  void *vp ;
  
  for (i = 1 ; i < arrayMax (aa) ; i++)
    {
      pp = arrayp (aa, i, pthread_t) ;
      if (!pthread_join (*pp, &vp)) 
	{ if (assInt(vp) != i + 9)
	  printf("bad value returned by joined thread %d \n", i) ;
	}
      else
	printf("cannot join phread %d\n", i) ;
      pthread_detach (*pp) ;

    }

  return ;
}

int main (int argc, char ** argv)
{
  int mode = 0, n, data, loop = 0, ii ;
  KEY key = 0 ;
  Array aa = 0 ; 
  pthread_t *pp ;


  printf("ok\n") ;

  if (argc>1)
    mode = atoi(argv[1]) ;
  else
    { 
      printf ("Usage ts 1 (no threading) (threading)\n") ;
      exit (1) ;
    }
  cache = assCreate () ;

  printf("Start %s\n", timeShowNow()) ; 

  switch (mode)
    {
    case 1:
      n = MAXN ;
      while (n--)
	{ 
	  key = 1 + n ;
	  data = 10 + n ;
	  writeToDisk (key, data) ;
	}
      break ;
    case 2:
      aa = arrayCreate (10000, pthread_t) ;
      pp = arrayp (aa, arrayMax (aa), pthread_t) ; /* for convenience, avoid index zero */

      n = MAXN ;
      while (n--)
	{ 
	  pp = arrayp (aa, arrayMax (aa), pthread_t) ;
	  key = 1 + n ;
	  ii = pthread_create(pp, 0, doWriteToDisk, assVoid (key)) ;
	  if (ii)
	    printf("Error %d  creating thread n = %d \n", ii, n) ;
	}
      collectWrite (aa) ;
      break ;
    }

  printf("Check 1 %s\n", timeShowNow()) ;

ici:

  switch (mode)
    {
    case 1:
      n = MAXN ;
      while (n--)
	{ 
	  key = 1 + n ;
	  data = getObj (key) ;
	  if (data != key + 9)
	    printf ("error key = %d data = %d != %d\n",key, data, key+9) ;
	  destroyObj (key) ;
	}
      break ;
    case 2:
     aa = arrayReCreate (aa, 10000, pthread_t) ;
      pp = arrayp (aa, arrayMax (aa), pthread_t) ; /* for convenience, avoid index zero */

      n = MAXN ;
      while (n--)
	{ 
	  key = 1 + n ;
	  pp = arrayp (aa, key, pthread_t) ;
	  ii = pthread_create(pp, 0, doCheckData, assVoid (key)) ;
	  if (ii)
	    printf("Error %d  creating thread n = %d \n", ii, n) ;
	}
      collectCheck (aa) ;
      break ;
        break ;
    }
  if (!loop++)
    {
      printf("Check 2 %s\n", timeShowNow()) ;
      goto ici ;
    }
      
  printf("Stop %s\n", timeShowNow()) ;
  return 0 ;
}


static void* doGetLength (void *vp)
{
  AC_OBJ obj = 0 ;
  AC_TABLE tbl = 0 ;
  int  ii = 0 , nSeq = 0 , total = 0 ;
  KEY key = assInt (vp)  ;
  BOOL debug = TRUE ;

  if (debug) printf("doGetLength %s starts \n", ac_key_name(key) ) ;
  if (TRUE)
    { 
      obj = ac_key2obj (DB, key, TRUE, 0) ;
      tbl = obj ? ac_tag_table (obj, "DNA", 0) : 0 ;
      if (tbl && tbl->rows)
	{
	  total += ac_table_int (tbl, 0, 1, 0) ;
	  nSeq++ ;
	}
      ac_free (tbl) ;
      ac_free (obj) ;
    }

  if (debug) printf("doGetLength %s done, n = %d \n", ac_key_name(key), ii ) ;

  return assVoid (ii) ;
}


static int collectLength (void)
{
  int i, n = 0 ;
  pthread_t *pp ;
  void *vp ;
  Array aa = arrayCreate (60, pthread_t) ;
  
  for (i = 0 ; i < arrayMax (aa) ; i++)
    {
      pp = arrayp (aa, i, pthread_t) ;
      if (!pthread_join (*pp, &vp))
	n += assInt (vp) ;  
      else
	printf("cannot join phread %d\n", i) ;
      pthread_detach (*pp) ;

    }
  arrayDestroy (aa) ;
  return n ;
}

static int getlength (KEY key, BOOL multi)
{
  /*  pthread_attr_t tha; */
  pthread_t *pp ;
  int n = 0 ;
  Array aa = arrayCreate (10000, pthread_t) ;

  if (!aa) 
    { aa = arrayCreate (10000, pthread_t) ;
    array (aa, arrayMax (aa), pthread_t) ; /* for convenience, avoid index zero */
    }
  pp = arrayp (aa, arrayMax (aa), pthread_t) ;
  printf("call phread %d on %s\n", arrayMax(aa) - 1, ac_key_name(key) ) ;
  if (multi)
    pthread_create(pp, 0, doGetLength, assVoid (key)) ;
  else
    n = assInt (doGetLength (assVoid (key))) ;
  return n ;
}


#ifdef JUNK
  there is alreadya main function above
int main2 (int argc, char **argv)
{ 
  KEY ks ;
  AC_OBJ obj = 0 ;
  AC_TABLE tbl = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  int n, ii = 0 , nSeq = 0 , total = 0, mode = 1 ;
  KEY key ;
  char *error = 0 ;
  BOOL debug = FALSE ;
      
  DB = ac_open_db (argc>1 ? argv[1]:".", &error) ;
  if (!DB)
    {
      fprintf (stderr, "Cannot open database %s : %s", argc>1 ? argv[1]:".", error) ;
      exit (1) ;
    }
  if (argc>2)
    mode = atoi(argv[2]) ;

  switch (mode)
    {
      case 1:
	ks = query(0, "FIND Sequence") ;
	if (!keySetMax(ks))
	  { 
	    fprintf (stderr, "No Sequence in this database, sorry\n%s\n") ;
	    exit (1) ;
	  }
	
	for (n = 0 ; n < keySetMax(ks) ; n++)
	  { 
	    key = keySet(ks,n) ;
	    obj = ac_key2obj (DB, key, h) ;
	    tbl = obj ? ac_tag_table (obj, "DNA", h) : 0 ;
	    if (tbl && tbl->rows)
	      {
		total += ac_table_int (tbl, 0, 1, 0) ;
		nSeq++ ;
	      }
	    ac_free (obj) ;
	    if (debug) printf("%d\t%d\t%s\n",ii,total,ac_key_name(key)) ;
	  }

	break ;
 
    case 2:  /* multi phreaded test */
    case 3:
	ks = query(0, "FIND Sequence") ;
	if (!keySetMax(ks))
	  { 
	    printf("No Sequence in this database, sorry\n%s\n",aceErrorMessage(0)) ;
	    return 1 ;
	  }
	
	for (n = 0 ; n < keySetMax(ks) ; n++)
	  { 
	    key = keySet(ks,n) ;
	    if (keyFindTag (key, str2tag("DNA")))
	      {
		printf ("ok %s\n", ac_key_name(key)) ;
		nSeq++ ;
		total += getlength (key, mode == 3 ? TRUE : FALSE ) ;
	      }
	    else
	      printf ("no tag dna in %s\n", aceName(key)) ;
	  }
	if (mode == 3)
	  total += collectLength () ;
	
	break ;
    }
  printf("%d Sequences, %d bases, on average %g base per sequence\n",
	 nSeq,total, (float)(total)/(nSeq ? nSeq : 1)) ;

  printf ("closing\n") ;
  fflush (stdout) ;
  
  ac_close (DB) ;
  return 0 ;
}
#endif

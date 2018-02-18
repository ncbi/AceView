#include    "graph.h"
#include    <vstd.h>

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  smlEDIT Functionality : Simple Memory Line Editor
_/
_/  This functionality does nothing with graphics,
_/  it only manages memory buffer for lines of the editor
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    Flags for marking the editor or lines.

_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

#define smlEDITDIRTY    0x00000001 /* editor content was modified */
#define smlLINEUPDATE   0x10000000 /* line content was modified */
#define smlLINESHIFT    0x20000000 /* line was shifted in the list */

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    smlLINENUM      get the pointer to the line with asolute number of num
    smlEditLinePtr  get the pointer to the line with ordinal number of num
    smlLINATRIB     get the pointer to the attribute buffer 

_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
#define smlLINENUM(_v_sml,_v_num)  ((smlLINE *) (((char * ) ((_v_sml) ->lineList)) +(sizeof(smlLINE) +((_v_sml) ->lineBufSize) *2) *(_v_num)  )) 
#define smlEditLinePtr(_v_sml,_v_num)  smlLINENUM((_v_sml) ,(smlLINENUM((_v_sml) ,(_v_num)) ->orderIndex)) 
#define smlLINATRIB(_v_sml,_v_sl)  (((smlLINE*) (_v_sl)) ->buf+((_v_sml) ->lineBufSize)  ) 




/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    smlLINE     line information structure 
    smlEDIT     sml editor infromation structure 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

typedef struct _smlLINE
{
  int     orderIndex ; /* the absolute number of the line having the current ordinal number */
  int     freeNum ;    /* the free lines stack cell */
  int     nextIndex,prevIndex ; /* the next and the previous lines indices */
  
  int     len ;        /* the length of the buffer */
  int     srvFlags ;   /* line service flags */
  
  void *  usrInfo ;    /* user specific pointer */
  int     usrFlags ;   /* user specific flags */
  
  char    buf[4] ;     /* it is 4 to be aligned */
} smlLINE ;

typedef struct _smlEDIT
{
  smlLINE  *  lineList ;     /* line array pointer */  
  
  int         lineMaxSize ;  /* the maximum size of the line */
  int         lineBufSize ;  /* thr max size aligned to paragraph */
  int         lineBufQuant ; /* every time we allocate memory for this many new lines */
  
  int         lineCnt,linePlace ; /* how many lines we have and how much is the maximum lpace */
  int         lineFreeCount ;  /* how many free lines do we have */
  int         lineFirstNum ;   /* the first line in the list */
  
  int         bufLen ;         /* the total length of the buffer */
  
  int         editFlags ;      /* editor's general flags */
  char        separ[1] ;       /* line wrapping separators */
} smlEDIT ;


/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  
  smlEditGetFreeLine(smlEDIT * sml,int markAsBusy) 
  
  gets a free line from the free line stack. If there 
  is no more lines it allocates more. If markAsBusy is not 
  set the line stays free.
  
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


static int smlEditGetFreeLine (smlEDIT * sml,int markAsBusy) 
{
  smlLINE *   curLine ;
  int         curNum,sizeLine,iL ;
  void *      newBufr ;
  
  /* ************** reallocate the space ************ */                 
  if (!sml->lineFreeCount) 
    {    /* if there are no free cells allocate them */
      /* compute the new quantized size for smlLINE structures */
      sml->linePlace = (( (sml->lineCnt) /sml->lineBufQuant) +1) *sml->lineBufQuant ;
      /* reallocate the space */
      sizeLine = sizeof (smlLINE) +sml->lineBufSize*2 ; /* for attribute buffer also */
      newBufr = messalloc (sizeLine* (sml->linePlace+1))  ;
      if (!newBufr) return 0 ;
      memset (newBufr,0,sizeLine* (sml->linePlace+1))  ;
      
      /* copy the old elements and replace the buffer */
      if (sml->lineList) 
	{
	  memmove (newBufr, (void *) sml->lineList,sml->lineCnt*sizeLine)  ;
	  messfree (sml->lineList)  ;
        }
      sml->lineList = (smlLINE *) newBufr ;
      
      for (iL = sml->linePlace-1 ;iL >= sml->lineCnt ;iL--) 
	{
	  /* keep the free free lines info */
	  curLine = smlLINENUM (sml,sml->lineFreeCount)  ; 
	  curLine->freeNum = iL ;
	  sml->lineFreeCount++ ;
        }
    }
  
    /* take the top free line */
  curLine = smlLINENUM (sml, (sml->lineFreeCount-1))  ;
  curNum = curLine->freeNum ;
  
  /* remove this line from free lines list */
  if (markAsBusy) 
    {
      curLine->freeNum = -1 ;
      sml->lineFreeCount-- ;
    }
  
  /* mark it as dead end */
  curLine = smlLINENUM (sml,curNum)  ;
  /*  curLine->orderIndex = sml->lineCnt ; */ /* do not do this - even free lines contain order information */
  curLine->nextIndex = -1 ;
  curLine->prevIndex = -1 ;
  curLine->len = 0 ;
  /*    curLine->buf[0] = 0 ;*/ /* do not do this bacuase we might 
	free lines while hoping to copy 
	  their content */
  
  curLine->usrInfo = 0 ;
  curLine->srvFlags = 0 ;
  sml->lineCnt++ ;
  
  return curNum ;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    smlEditFreeLine (smlEDIT * sml,int num) 

    frees the line and puts it into free lines stack. 
    The line is marked as busy and can be requested 
    by a smlEditGetFreeLine function.

_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static void smlEditFreeLine (smlEDIT * sml,int num) 
{
  smlLINE * curLine ;
  
  /* mark the line as free */
  curLine = smlLINENUM (sml,sml->lineFreeCount)  ;
  curLine->freeNum = num ;
  sml->lineFreeCount++ ;
  sml->lineCnt-- ;
  
  curLine = smlLINENUM (sml,num)  ;
  curLine->prevIndex = -1 ;
  curLine->nextIndex = -1 ;    
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    smlEditChangeLine (smlEDIT * sml,char * buf,int bufSize,int startRow,int startCol,int endRow,int endCol) 

    replaces the content of the buffer starting from 
    (startRow,startCol)  to (endRow,endCol)  by bufSize 
    characters of a given buffer (buf) . 
    
    If bufSize is zero bufSize is taken as strlen (buf) . 
    If buf is zero the specified range is deleted. 
    If any of startRow, or endRow is equal to -1 they are 
    taken as the ordinal number of the last row. 
    If startCol or endCol is -1 they are taken as the last 
    column in corresponding lines.


    can be used to replace, delete, insert any number of 
    lines into any position in the buffer.
     
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

int smlEditChangeLine (smlEDIT * sml,char * buf,int bufSize,int startRow,int startCol,int endRow,int endCol) 
{
  int         iB,iL,insertCnt,originalBufr,initCntLines ;
  char        ch,* curBuf ;
  int         curNum,nxtNum,prvNum,freeNum,ordFixNum ;
  smlLINE *   curLine,*tmpLine ;
  char        delBuf[2] ;
  
  /* fix the size of buffer if not defined */ 
  if (!buf) 
    {buf = delBuf ;delBuf[0] = 0xFF ;bufSize = 1 ;}
  if (!bufSize) bufSize = strlen (buf)  ;
  
  initCntLines = sml->lineCnt ; /* remember the number of lines at the beginning */
  
  if (!sml->lineCnt) 
    {
      startRow = smlEditGetFreeLine (sml,1)  ; /* ask for a single line to allocate space */
      endRow = startRow ;endCol = 0 ;startCol = 0 ;
      ordFixNum = 0 ;
      curLine = smlLINENUM (sml,0)  ;
    }
  else 
    {

      /* adjust the size and position to start */
      if (startRow == -1 || startRow >= sml->lineCnt) startRow = sml->lineCnt-1 ;
      if (endRow == -1 || endRow >= sml->lineCnt) endRow = sml->lineCnt-1 ;
      
      /* correct the order of rows */
      if (startRow>endRow) 
	{
	  curNum = endRow ;endRow = startRow ;startRow = curNum ;
	  curNum = endCol ;endCol = startCol ;startCol = curNum ;
        }

      /* get the real order of lines containing this positions in our lineList */
      ordFixNum = startRow ; /* we start fixing line sort orders from here */
      startRow = smlLINENUM (sml,startRow) ->orderIndex ;
      endRow = smlLINENUM (sml,endRow) ->orderIndex ;

      /* correct misspositions */    
      curLine = smlLINENUM (sml,startRow)  ;
      if (startCol == -1 || startCol >= curLine->len) 
	{ /* put on the end of the line */
	  for (startCol = 0 ;startCol<curLine->len ;startCol++) if (curLine->buf[startCol] == '\n') break ; /* look for '\n' */
        }
      
      curLine = smlLINENUM (sml,endRow)  ;
      if (endCol == -1 || endCol >= curLine->len) 
	{
	  for (endCol = 0 ;endCol<curLine->len ;endCol++) if (curLine->buf[endCol] == '\n') break ; /* look for '\n' */
	  if (curLine->buf[endCol] == '\n') endCol++ ;
        }

      /* correct the column positions */
      if (startRow == endRow && startCol>endCol) 
	{curNum = endCol ;endCol = startCol ;startCol = curNum ;}
    }
        
  
  /* get the next of endrow line to glue after our insertion */    
  nxtNum = curLine->nextIndex ; /* by the way curLine was pointing to enrow */
  
  /*  if there is something left on the endRow line 
        transfer the leftover to a free cell an make it 
        the nextNum to be glued after insertion */
  if ((iL = curLine->len-endCol) >0) 
    { /* the length of the endRow line */
      nxtNum = smlEditGetFreeLine (sml,1)  ; /* we need to glue this line after our insertion */
      tmpLine = smlLINENUM (sml,nxtNum)  ;
      
      memmove (tmpLine->buf,curLine->buf+endCol,iL)  ;tmpLine->buf[iL] = 0 ; /* transfer the leftover */
      tmpLine->len = iL ; /* fix the length */
      
      /* copy properties of the endRow line */
      tmpLine->nextIndex = curLine->nextIndex ;
      tmpLine->prevIndex = curLine->prevIndex ;
      tmpLine->usrInfo = curLine->usrInfo ;
      tmpLine->srvFlags = curLine->srvFlags ;
    }
  
  /* free lines inside of the selection in the opposite order */
  for (freeNum = endRow ;freeNum != startRow ;) 
    { 
      curLine = smlLINENUM (sml,freeNum)  ;
      prvNum = freeNum ;
      freeNum = curLine->prevIndex ; /* change the following line currently ought to be freed  */
      smlEditFreeLine (sml,prvNum)  ;   /* mark the line as free */        
      sml->bufLen -= curLine->len ;  /* adjust the total buffer size */
    }
  
  /* set the beginning of the selection replacement */
  curNum = startRow ;
  curLine = smlLINENUM (sml,curNum)  ;
  if (startCol >= sml->lineMaxSize) curBuf = 0 ; /* this forces line change */
  else curBuf = curLine->buf+startCol ;  /* continue in the current line */
  iL = startCol ;
  insertCnt = 0 ;
  originalBufr = 1 ;
  
  if (curBuf) 
    {
      insertCnt++ ;
      curLine->srvFlags |= smlLINEUPDATE ; /* mark the previous line to be updated */
      curLine->len = iL ;
      curBuf[iL] = 0 ;
    }
  
    for (iB = 0 ;iB<bufSize ;iB++) 
      {
        ch = buf[iB] ;
        
        if (!curBuf) 
	  {    /* occupy new line buffer */
            /* occupy another free line  */
            prvNum = curNum ;
            curNum = smlEditGetFreeLine (sml,1)  ;
            curLine = smlLINENUM (sml,curNum)  ;
            

            /* link the newly added line to the previous line list */            
            curLine->prevIndex = prvNum ;
            smlLINENUM (sml,prvNum) ->nextIndex = curNum ;
            
            /* change to the next buffer */
            curBuf = curLine->buf ;iL = 0 ;   
            insertCnt++ ;    /* one more line added */
            
            ordFixNum++ ; /* ordinal number of the current line */
	    
            curLine->srvFlags |= smlLINEUPDATE ; /* mark the previous line to be updated */
            /* set the ordinal number and updated flag */
            tmpLine = smlLINENUM (sml,ordFixNum)  ;
            if (tmpLine->orderIndex != curNum) 
	      {
                tmpLine->orderIndex = curNum ;
                curLine->srvFlags |= smlLINESHIFT ; /* mark the line as shifted */
	      }            
	  }
	
        /* copy the content to our line buffer */
        if ((unsigned char) ch< (unsigned char) 0xFF) 
	  {
            *curBuf = ch ; 
            iL++ ;curBuf++ ;
	  }
        
        /* replace the source buffer */
        if (iB == bufSize-1) 
	  { /* check if this was the last byte of the SOURCE buffer */
            /* if (nxtNum<0) ch = '\n' ; */
            if ((originalBufr || ch != '\n') && nxtNum >= 0) 
	      { /* paragraph is not yet over ? */
                
                tmpLine = smlLINENUM (sml,nxtNum)  ; /*  the next line is the source buffer */
                iB = -1 ;buf = tmpLine->buf ;bufSize = tmpLine->len ;
		
                prvNum = nxtNum ;  /* remember this line */
                nxtNum = tmpLine->nextIndex ; /* choose the next to be glued candidate */
                smlEditFreeLine (sml,prvNum)  ;    /* free this line */
		
                originalBufr = 0 ;/* now we are switching to the end of the paragraph */
	      }
        }
        
        if (ch == '\n' || iL >= sml->lineMaxSize) 
	  { /* terminate the Line */
            *curBuf = 0 ;      
            sml->bufLen += (iL-curLine->len)  ; /* adjust the total buffer size */
            curLine->len = iL ;                /* adjust the length of current string */
            curBuf = 0 ;                       /* this line buffer is finished */
	    
	  }
      }
    
    /* link the next line after the position where we have broken the list */
    curLine->nextIndex = nxtNum ; /* pointing the newcomer to the continuation */    
    if (nxtNum != -1) 
      { /* there was something after this block ? */
        smlLINENUM (sml,nxtNum) ->prevIndex = curNum ;
      }
    
    /* lines were shifted we need to reoder the rest of lines */
    if (initCntLines != sml->lineCnt) 
      {
        for (ordFixNum++ ;nxtNum>0 ;ordFixNum++) 
	  {
            tmpLine = smlLINENUM (sml,ordFixNum)  ;
            curLine = smlLINENUM (sml,nxtNum)  ;
            
            if (tmpLine->orderIndex != nxtNum) 
	      {
                tmpLine->orderIndex = nxtNum ; 
                curLine->srvFlags |= smlLINESHIFT ; /* mark the line as shifted */
	      }
            
            nxtNum = curLine->nextIndex ;
	  }
	
      }
    
    sml->editFlags |= smlEDITDIRTY ;
    return insertCnt ;
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    void smlEditMiniMax (smlEDIT * sml, int * sy,int * sx, int * fy,int * fx) 

    fixes the minal and maximal positions according to the 
    content of our strings buffer. It arranges so *sy comes 
    before *fy and *sx comes before *fx.

    fx and sx parameters can be zero, that means do not 
    adjust columns. 

    If lines are out of linelength or out of lineCnt they 
    will be properly adjusted. 
    
    -1 for any of this *fx, *sx, *fy, *sy parameters means the 
    last possible row or column.
     
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

void smlEditMiniMax (smlEDIT * sml, int * sy,int * sx, int * fy,int * fx) 
{
  smlLINE * sl ;
  int     curNum ;
  
  if (!sml->lineCnt) 
    {
      (*sx)  = (*sy)  ;
      if (fx && fy)  (*fx)  = (*fy)  ;
      return ;
    }
  
  /* fix out of buffer rows */
  if ((*sy)  == -1 || (*sy)  >= sml->lineCnt)  (*sy)  = sml->lineCnt-1 ;
  if ((*fy)  == -1 || (*fy)  >= sml->lineCnt)  (*fy)  = sml->lineCnt-1 ;
  
  /* correct the order of rows */
  if ((*sy) > (*fy)) 
    {
      curNum = (*fy)  ; (*fy)  = (*sy)  ; (*sy)  = curNum ;
        if (fx && fy) curNum = (*fx)  ; (*fx)  = (*sx)  ; (*sx)  = curNum ;
    }
  
  if (!sx || !fx) return ;
  
  /* fix out of line columns */
  sl = smlEditLinePtr (sml, (*sy))  ;
  if ((*sx)  == -1 || (*sx) >sl->len)  (*sx)  = sl->len ;
  sl = smlEditLinePtr (sml, (*fy))  ;
  if ((*fx)  == -1 || (*fx) >sl->len)  (*fx)  = sl->len ;
  
  /* fix reversed columns */
    if ((*sy)  == (*fy)  && (*sx) > (*fx)) 
      {curNum = (*fx)  ; (*fx)  = (*sx)  ; (*sx)  = curNum ;}        
    
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    int smlEditEnumLines (smlEDIT * sml,smlLINEENUMFunction enFunc,void * param,int startRow,int startCol,int endRow,int endCol) 

    Calls a given enFunc callback function with a give parameter 
    for all the lines in the range of  startRow, startCol,endRow,endCol. 

    enumFunc can be zero. (to count the size of the selection buffer) 
    callback has a prototype of int smlLINEENUMFunction (smlEDIT * sml,smlLINE * sl,int absPos,int startCol, int endCol, int row,void *param) 
    sl is the current line slLINE *
    startCol and endCol show the selection in the given line.
    absPos is its absolute position of the startCol in the total selection buffer
    param is the parameter passed to smlEditEnumLines
    returns the size of the selection buffer.

    
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
typedef int (* smlLINEENUMFunction)  (smlEDIT * sml,smlLINE * sl,int absPos,int startCol, int endCol, int row,void *param)  ;
int smlEditEnumLines (smlEDIT * sml,smlLINEENUMFunction enFunc,void * param,int startRow,int startCol,int endRow,int endCol) 
{
  smlLINE * sl ;
  int     iL,st,ed,cntAll ;
  
  if (!sml->lineCnt) return 0 ;
  /* adjust the size and position to start */
  smlEditMiniMax (sml,&startRow,&startCol,&endRow,&endCol)  ;
  
  for (cntAll = 0,iL = startRow ;iL <= endRow ;iL++) 
    {
      sl = smlLINENUM (sml,iL)  ;
      sl = smlLINENUM (sml,sl->orderIndex)  ;
      
      /* set up the beginning and the end of scanning in the current line */
      if (iL>startRow) st = 0 ;else st = startCol ;
      if (iL<endRow) ed = sl->len ;else ed = endCol ;
      if (enFunc) 
	{
	  if (!enFunc (sml,sl,cntAll,st,ed,iL,param)) break ;
        }
      cntAll += ed-st ;
    }
  return cntAll ;
}

int _smlEditGetBuffer (smlEDIT * sml,smlLINE * sl,int absPos,int startCol, int endCol, int row,void *param) 
{
  memmove (( (char *) param+absPos) ,sl->buf+startCol,endCol-startCol)  ;
  return 1 ;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    char * smlEditGetBuffer (smlEDIT * sml,int * bufSize,int startRow,int startCol,int endRow,int endCol) 

    allocates and copies the content of the buffer in the given range 
    startRow, startCol,endRow,endCol. 

    returns the resulting buffer.
    If bufSize is not zero the size of the buffer is copied to *bufSize

    
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
    
char * smlEditGetBuffer (smlEDIT * sml,int * bufSize,int startRow,int startCol,int endRow,int endCol) 
{
  char *  resultBuf ;
  int     size, size2 ;
  
  /* get the size */
  size = smlEditEnumLines (sml,0,0,startRow,startCol,endRow,endCol)  ;if (!size) return 0 ;
  /* allocate a buffer */
  resultBuf = messalloc (size+1)  ;if (!resultBuf) return 0 ;
  /* copy the buffer */
  size2 = smlEditEnumLines (sml,_smlEditGetBuffer,resultBuf,startRow,startCol,endRow,endCol)  ;
  if (size2 != size) messcrash ("inconsistency in smlEditGetBuffer") ;
  resultBuf[size] = 0 ;
  
  if (bufSize) *bufSize = size ;
  return resultBuf ;
}    


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    void smlEditPrintBuffer (smlEDIT * sml,int isOrdered,int cnt)  
    
    prints cnt lines of original (isOrdered = 0)  or ordered ( = 1)  
    presentation of the lines list.
    

_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/ 

void smlEditPrintBuffer (smlEDIT * sml,int isOrdered,int cnt) 
{
  int i ;
  smlLINE * l ;
  
  printf (
	  "lineCnt = %d linePlace = %d lineFirstNum = %d bufLen = %d\n"
	  ,sml->lineCnt,sml->linePlace,sml->lineFirstNum,sml->bufLen)  ;
  
  if (!cnt) cnt = sml->lineCnt ;
  if (isOrdered) i = sml->lineFirstNum ;else i = 0 ;
  for (i = 0 ;i<cnt ;i++) 
    {
      if (isOrdered) 
	{
	  l = smlLINENUM (sml,i)  ;
	  l = smlLINENUM (sml,l->orderIndex)  ;
        }
      else 
	{
	  l = smlLINENUM (sml,i)  ;
	  if (l->prevIndex == -1 && l->nextIndex == -1) cnt++ ;
        }
      
      printf ("i = %d flg = %x ord = %d prv = %d nxt = %d len = %d | buf = %s%s"
	      ,i ,l->srvFlags, l->orderIndex,l->prevIndex,l->nextIndex,l->len,l->buf,strchr (l->buf,'\n')  ? "" : "\n")  ;
    }
    printf ("Free Items\n")  ;
    
    i = sml->lineFreeCount-10 ;if (i<0) i = 0 ;
    for ( ;i<sml->lineFreeCount ;i++) 
      {
        l = smlLINENUM (sml,i)  ;
        printf ("%d ",l->freeNum)  ;
      }
    printf ("\n")  ;
    
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
    smlEDIT * smlEditInit (char * initText,int maxLineSize,char * separ) 
    void smlEditFree (smlEDIT * sml) 

    construction and destruction of the smlEDITOR with 
    initial iniText and wrapping size maxLineSize. 
    If separator is given lines are wrapped at one fo the symbols in *separ
    otherwise they are wrapped at the end of the line character.
     
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

smlEDIT * smlEditInit (char * initText,int maxLineSize,char * separ) 
{
  smlEDIT *   sml ;
  int         len ;
  
  /* set initial parameters */
  if (!maxLineSize) maxLineSize = 256 ;
  len = sizeof (smlEDIT) +maxLineSize ;
  if (separ) len += strlen (separ)  ;
  
  sml = (smlEDIT *) messalloc (len)  ;if (!sml) return 0 ;
  memset (sml,0,len)  ;
  sml->lineMaxSize = maxLineSize ;
  sml->lineBufSize = ((maxLineSize/0x10) +1) *0x10 ; /* to be aligned */
  sml->lineBufQuant = 256 ;
  if (separ) strcpy (sml->separ,separ)  ;
  
  if (initText && !smlEditChangeLine (sml,initText,0,0,0,0,0)) 
    return 0 ;
  
  return sml ;
}

void smlEditFree (smlEDIT * sml) 
{
  if (!sml) return ;
  
  if (sml->lineList) messfree (sml->lineList)  ;
  messfree (sml)  ;
}


/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  
  char smlEditFixPosition (smlEDIT * sml, smlLINE * * psl, int * px, int * py,int isJmpNext) 

  adjusts the given position to our buffers content.
  out-of-line columns or out-of-buffer lines are brought 
  inside. -1 for position means the last row or columnt 
  correspondingly .
  
  if isJmpNext columns more than the line lenght are set 
  to the beginning of the next line.
  
  returns the character at the adjusted position
    
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

char smlEditFixPosition (smlEDIT * sml, smlLINE * * psl, int * px, int * py,int isJmpNext) 
{
  smlLINE * sl ;
  int     curx = *px ; 
  int     cury = *py ;
  int     len,yWasZero ;
  
  if (!sml->lineCnt) 
    {*px = *py = 0 ;return 0 ;}
  if (cury <= 0) 
    {cury = 0 ;yWasZero = 1 ;}else yWasZero = 0 ;
  
  /* fix the Y position */
  if (cury >= sml->lineCnt) cury = sml->lineCnt-1 ;
  if (curx<0 && cury>0) cury-- ;
  
  sl = smlEditLinePtr (sml,cury)  ;
  len = sl->len ; /* the length without '\n' */
  if ((!isJmpNext)  && len && sl->buf[len-1] == '\n') len-- ; 
  
  if (curx >= len) 
    {
      if (isJmpNext && cury<sml->lineCnt-1) 
	{
	  curx = 0 ;
	  cury++ ;sl = smlEditLinePtr (sml,cury)  ;
        }
      else curx = len ;
      
      if (curx >= sml->lineMaxSize) curx = len-1 ;
    }
  if (curx<0 && !yWasZero) curx = len-1 ;if (curx<0) curx = 0 ;
  *px = curx ; 
  *py = cury ;
  if (psl) *psl = sl ;
  return sl->buf[curx] ;
}



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    smlEditMarkLine (smlEDIT* sml,int sRow, int fRow,int flag) 
    smlEditUnmarkLine (smlEDIT* sml,int sRow, int fRow,int flag) 

    Marks or unmarks all lines from within a given range sRow,fRow
    by a flag.

_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static int _smlEditMarkLine (smlEDIT * sml,smlLINE * sl,int absPos,int startCol, int endCol, int row,void *param) 
{ (sl->srvFlags)  |= (* ((int *) param))  ;return 1 ;}

static int _smlEditUnmarkLine (smlEDIT * sml,smlLINE * sl,int absPos,int startCol, int endCol, int row,void *param) 
{ (sl->srvFlags)  &= (0xFFFFFFFF^ (* ((int *) param)) )  ;return 1 ;}

void smlEditMarkLine (smlEDIT* sml,int sRow, int fRow,int flag) 
{
    smlEditEnumLines (sml,_smlEditMarkLine, (void*)  (&flag) ,sRow,0,fRow,-1)  ;
}

void smlEditUnmarkLine (smlEDIT* sml,int sRow, int fRow,int flag) 
{
    smlEditEnumLines (sml,_smlEditUnmarkLine, (void*)  (&flag) ,sRow,0,fRow,-1)  ;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    end of smEDITor section

_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


















#include    <gmleditor.h>


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Xtra-Graphics Functionality : 
_/
_/  draws fancy staff or eazes drawing of other primitives
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* 
    draws dashed line of given coordinartes where dash length  =  lng
*/
static void graphPatternedXorLine (float lng,float sx,float sy,float fx,float fy) 
{
  float leen,l,prvl,cosx,cosy ;
  int     i ;
  
  if (lng == 0) 
    {
      graphXorLine (sx,sy,sx,sy)  ;
      return ;
    }
  cosx = fx-sx ;cosy = fy-sy ;
  leen = sqrt ( cosx*cosx+cosy*cosy )  ;
  cosx /= leen ;cosy /= leen ;
  
  for (i = 0,l = 0,prvl = 0 ;l<leen ;i++) 
    {
      l += lng ;
      if (l> (leen)) l = leen ;
      if (! (i%2)) graphXorLine (sx+prvl*cosx,sy+prvl*cosy,sx+l*cosx,sy+l*cosy)  ;
      prvl = l ;
    }
}

/* 
    draws colored (foreGround = fg, bkGround = bg)  polygon (empty or filled)  
    coordintes are taken from arr, scaled (sclx,sclx)  and shifted (ofx,ofy) 
*/

static void _gmlDrawPolygon (float * arr,int fg, int bg,int filled, int dim,float ofx,float ofy, float sclx, float scly) 
{
  Array   temp ;
  int     i ;
  float   x,y ;
  
  temp  =  arrayCreate ((2*dim) , float)  ;
  for (i = 0 ;i<2*dim ;i++) 
    {
      x = arr[i]*sclx+ofx ;
      array (temp, i, float)   =  x ;
      i++ ;
      y = arr[i]*scly+ofy ;
      array (temp, i, float)   =  y ;
    }
  
  if (filled) 
    {
      graphColor (bg)  ;
      graphPolygon (temp)  ;
    }
  graphColor (fg)  ;
  graphLineSegs (temp)  ;
  
  arrayDestroy (temp)  ;
}

























/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  gmlEDIT STRUCTURES AND MACROS 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


#define gmlEDITsetCallbackFlag(_v_gml,_v_flag)  ((_v_gml) ->callBackFlags)  |= ((0x00000001) << (_v_flag)) 

/* scroller attributes */
#define gmlVERTRSCRL     0x00000001 /* vertical right side */
#define gmlVERTLSCRL     0x00000002 
#define gmlVERTSCRL      (gmlVERTRSCRL|gmlVERTLSCRL) 
#define gmlHORZBSCRL     0x00000004 /* horizontal bottom side */
#define gmlHORZTSCRL     0x00000008
#define gmlHORZSCRL      (gmlHORZBSCRL|gmlHORZTSCRL) 

#define gmlCHNGSCRL     0x00000001 /* changed limits or repositioned */
#define gmlREPOSSCRL    0x00000002 /* repositioned */
#define gmlONSCRL       0x00000004 /* currently grabbed  */
#define gmlDRGBOX       0x00000008 /* XOr box is drawn */


#define gmlLINETOREDRAW     0x00000001 /* lines marked to redraw */

/*
    scroller information structure 
*/

typedef struct _gmlSCRL
{
  int     style,status ;   /* attrributes and the state */
  double  minVal,maxVal,curVal,scale ; /* ranges and the position */
  float   sx,sy,ext ;  /* original geomtery */
  
  /* coordinaes for primitives */
  float   elevSx,elevSy,elevFx,elevFy,oldW ;
  float   lineSx,lineSy,lineFx,lineFy ;
  float   lExt,drawScale,drawShift,diffVal ;
  float   minFigSx,minFigSy,maxFigSx,maxFigSy ;
  float   minFigFx,minFigFy,maxFigFx,maxFigFy ;
  
  float   dragX,dragY ; /* current drag coodinates */
  
} gmlSCRL ;



typedef struct _graphMultiLineEdit
{
  int     style ; /* the style of our editor */
  gmlEDITCONFIG cfg ; /* configuration */
  smlEDIT  *  sml ; /*smlEDITOR buffer holding our content */
  
  /* coordinates */
  float   sx,sy,cx,cy ;
  float   textCx,textCy ;
  
  /* curent cursor and selection */
  int     cursorX,cursorY,selectX,selectY,cursorWishX ;
  
  
  int     updatedLastLine ;    /* different dynamic information */
  
  /* callbacks */
  int     callBackFlags ;
  gmlEDITCallback callbackFunc[gmlEDITlastCALLBACK] ;
  void *          callbackParam[gmlEDITlastCALLBACK] ;
  
  /*scrollers info */
  gmlSCRL scrlVert,scrlHoriz ;
  int     boxScrlV,boxScrlH ;
  
  /* a little magic from Jean to Vahan - boxes graphs and so on */
  int *   magic ; 
  Graph   subGraph ;
  int     boxesCounter,maxBoxesCount ;
  struct  
  {int num ;}boxInfo[1] ;
  
} gmlEDIT ;















/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  gmlSCROLL Functionality : Graphics Memory Line Scroller
_/
_/  This is to be implemented from gmlEDIT 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/* 
    initilizes the scroller ot the given type styl with 
    given coordinates and given extent (ext) 

    return's the width it occupies.
*/

static float gmlScrollInit (gmlEDIT * gml,gmlSCRL * gsc,int styl,int sx,int sy,int ext) 
{
  float   wid ;
  /* scrollers settings */
  gsc->style = styl ;
  gsc->ext = ext ;
  
  wid = gml->cfg.scrlWidth ;
  if (styl&gmlHORZSCRL) wid *= gml->cfg.propVert_to_Horiz ;
  
  if (styl&gmlVERTRSCRL) 
    {
      gsc->sx = sx-wid ;gsc->sy = sy ;
    }
  else if (styl&gmlVERTLSCRL) 
    {
      gsc->sx = sx ;gsc->sy = sy ;
    }
  else if (styl&gmlHORZBSCRL) 
    {
      gsc->sx = sx ;gsc->sy = sy-wid ;
    }
  else if (styl&gmlHORZTSCRL) 
    {
      gsc->sx = sx ;gsc->sy = sy ;
    }
  return wid ;
}


/*    draws xored elevator for drag mode */
static void gmlDrawXorScroll (gmlEDIT * gml,gmlSCRL * gsc) 
{
  if (!gsc->style) return ;
  graphColor (gml->cfg.xorBoxColor)  ;
  
  graphPatternedXorLine (gml->cfg.selPattern, gsc->elevSx,gsc->elevSy,gsc->elevFx,gsc->elevSy)  ;
  graphPatternedXorLine (gml->cfg.selPattern, gsc->elevFx,gsc->elevSy,gsc->elevFx,gsc->elevFy)  ;
  graphPatternedXorLine (gml->cfg.selPattern, gsc->elevFx,gsc->elevFy,gsc->elevSx,gsc->elevFy)  ;
  graphPatternedXorLine (gml->cfg.selPattern, gsc->elevSx,gsc->elevFy,gsc->elevSx,gsc->elevSy)  ;
}

/* 
    set's the range (minV, maxV) ,  the scale and 
    the coordinate of the scroller curV
*/
        
static void gmlScrollSet (gmlEDIT * gml,gmlSCRL * gsc,int minV, int maxV,int curV,int scale) 
{
  float   diffVal ;
  
  if (!gsc->style) return ;
  
  /* adjust to limits */
  if (maxV<minV) maxV = minV ;if (curV<minV) curV = minV ;if (curV>maxV) curV = maxV ;
  /* did anything change */
  if (gsc->drawScale && gsc->minVal  ==  minV && gsc->maxVal  ==  maxV && gsc->curVal  ==  curV ) return ;
  if (gsc->curVal  !=  curV ) gsc->status |= gmlREPOSSCRL ;
  
  /* accept the new ones */
  gsc->minVal = minV ;
  gsc->curVal = curV ;
  gsc->maxVal = maxV ;
  gsc->scale = scale ;
  
  gsc->status |= gmlCHNGSCRL ; /* set the flag saying we changed something */
  
  /* comute the cordinates for primitves */
  diffVal = gsc->maxVal-gsc->minVal ;
  /* the length of elevator's possible shift positions */
  gsc->lExt = gsc->ext-2*gml->cfg.scrlArrowLen ; 
  
  /* the length of the elevator */
  if (diffVal) 
    {
      /*gsc->drawScale = (gsc->scale) *gsc->lExt/diffVal ;*/
      gsc->drawScale = gsc->scale*gsc->lExt/ (diffVal+gsc->scale)  ;
      if (gsc->drawScale>gsc->lExt) gsc->drawScale = gsc->lExt ;if (gsc->drawScale<0) gsc->drawScale = 0 ;
    }else gsc->drawScale = gsc->lExt ;
  
  /* offset where start drawing the elevator */
  if (diffVal) 
    {
      gsc->drawShift = (gsc->curVal-gsc->minVal) * (gsc->lExt-gsc->drawScale) /diffVal ;
      if (gsc->drawShift> (gsc->lExt-gsc->drawScale)) gsc->drawShift = (gsc->lExt-gsc->drawScale)  ;if (gsc->drawShift<0) gsc->drawShift = 0 ;
    }else gsc->drawShift = 0 ;
  
  
    if (gsc->style&gmlVERTSCRL)    /* vertical */
      
      {  
        gsc->lineSx = gsc->lineFx = gsc->sx+gml->cfg.scrlWidth/2 ;
        gsc->lineSy = gsc->sy+gml->cfg.scrlArrowLen ;gsc->lineFy = gsc->lineSy+gsc->lExt ;
        
        gsc->elevSx = gsc->lineSx-gml->cfg.scrlWidth/2 ;gsc->elevFx = gsc->elevSx+gml->cfg.scrlWidth ;
        gsc->elevSy = gsc->lineSy+gsc->drawShift ;gsc->elevFy = gsc->elevSy+gsc->drawScale ;
        gsc->minFigSx = gsc->sx ;gsc->minFigSy = gsc->sy ;
        gsc->maxFigSx = gsc->sx ;gsc->maxFigSy = gsc->sy+gsc->ext-gml->cfg.scrlArrowLen ;
	
        gsc->minFigFx = gsc->minFigSx+gml->cfg.scrlWidth ;gsc->minFigFy = gsc->minFigSy+gml->cfg.scrlArrowLen ;
        gsc->maxFigFx = gsc->maxFigSx+gml->cfg.scrlWidth ;gsc->maxFigFy = gsc->maxFigSy+gml->cfg.scrlArrowLen ;
	
        if (gsc->status&gmlREPOSSCRL) gmlEDITsetCallbackFlag (gml,gmlEDITchangeVSCROLL)  ;
      }
    else if (gsc->style&gmlHORZSCRL)  /* horizontal */
      
      {  
        gsc->lineSy = gsc->lineFy = gsc->sy+gml->cfg.scrlWidth*gml->cfg.propVert_to_Horiz/2 ;
        gsc->lineSx = gsc->sx+gml->cfg.scrlArrowLen/gml->cfg.propVert_to_Horiz ;gsc->lineFx = gsc->lineSx+gsc->lExt ;
        
        gsc->elevSy = gsc->lineSy-gml->cfg.scrlWidth*gml->cfg.propVert_to_Horiz/2 ;gsc->elevFy = gsc->elevSy+gml->cfg.scrlWidth*gml->cfg.propVert_to_Horiz ;
        gsc->elevSx = gsc->lineSx+gsc->drawShift ;gsc->elevFx = gsc->elevSx+gsc->drawScale ;
        gsc->minFigSx = gsc->sx ;gsc->minFigSy = gsc->sy ;
        gsc->maxFigSx = gsc->sx+gsc->ext-gml->cfg.scrlArrowLen/gml->cfg.propVert_to_Horiz ;gsc->maxFigSy = gsc->sy ;
	
        gsc->minFigFx = gsc->minFigSx+gml->cfg.scrlArrowLen/gml->cfg.propVert_to_Horiz ;gsc->minFigFy = gsc->minFigSy+gml->cfg.scrlWidth*gml->cfg.propVert_to_Horiz ;
        gsc->maxFigFx = gsc->maxFigSx+gml->cfg.scrlArrowLen/gml->cfg.propVert_to_Horiz ;gsc->maxFigFy = gsc->maxFigSy+gml->cfg.scrlWidth*gml->cfg.propVert_to_Horiz ;
	
        if (gsc->status&gmlREPOSSCRL) gmlEDITsetCallbackFlag (gml,gmlEDITchangeHSCROLL)  ;
      }
}


/* 
    reacts on mouse event at given coordinates 
    returns 0 if the event was not treated 
    otherwise 1.
*/

static int gmlScrollMouseAction (gmlEDIT * gml,gmlSCRL * gsc,double x, double y,int button) 
{
  double    chg = 0 ;
  
  if (!gsc->style) return 0 ;
  if ((! (gsc->status&gmlONSCRL))  && ( x<gsc->minFigSx || x>gsc->maxFigFx || y<gsc->minFigSy || y>gsc->maxFigFy )) 
    return 0 ;
  
  switch (button) 
    {
      
    case MIDDLE_DOWN:
      gsc->status |= gmlONSCRL ;  /* set drag mode */
      gsc->status &= (0xFFFFFF^gmlDRGBOX)  ;
      gsc->dragX = gsc->drawScale/2 ; /*remember where we have dragged it */ 
      gsc->dragY = gsc->drawScale/2 ;
      
    case MIDDLE_DRAG:
    case LEFT_DRAG:
      {
	
        if (! (gsc->status&gmlONSCRL)) return 0 ;
        if (gsc->dragX == -1  && gsc->dragY == -1) return 1 ;
        
        /*compute the new scroller position */
        if (gsc->style&gmlHORZSCRL) 
	  chg = gsc->minVal+ ((x-gsc->dragX) -gsc->lineSx) * (gsc->maxVal-gsc->minVal) / (gsc->lExt-gsc->drawScale) -gsc->curVal ;
        else if (gsc->style&gmlVERTSCRL)             
	  chg = gsc->minVal+ ((y-gsc->dragY) -gsc->lineSy) * (gsc->maxVal-gsc->minVal) / (gsc->lExt-gsc->drawScale) -gsc->curVal ;            
        /* hide the scroller's XORed elevator */
        if (gsc->status&gmlDRGBOX) gmlDrawXorScroll (gml,gsc)  ;
        
        /* change the scroller position */
        gmlScrollSet (gml,gsc,gsc->minVal,gsc->maxVal,gsc->curVal+chg,gsc->scale)  ;
        /* show the scroller's XORed elevator */
        if (button != MIDDLE_DOWN) 
	  {gmlDrawXorScroll (gml,gsc)  ;gsc->status |= gmlDRGBOX ;}
      }
      
    return 1 ;
    case MIDDLE_UP:
    case LEFT_UP:
      { /* releasing the scroller */
        if (! (gsc->status&gmlONSCRL)) return 0 ;
        if (gsc->status&gmlDRGBOX) gmlDrawXorScroll (gml,gsc)  ;
        gsc->status &= (0xFFFFFFFF^gmlONSCRL)  ;
        gsc->status &= (0xFFFFFFFF^gmlDRGBOX)  ;
      }
      return 1 ; 
    case LEFT_DOWN:
      {
        gsc->status |= gmlONSCRL ;  /* set drag mode */
        gsc->status &= (0xFFFFFF^gmlDRGBOX)  ;
        gsc->dragX = -1 ;gsc->dragY = -1 ; /*remember where we have dragged it */ 
        if (x >= gsc->elevSx && x <= gsc->elevFx && y >= gsc->elevSy && y <= gsc->elevFy ) 
	  {
            gsc->dragX = x-gsc->elevSx ; /*remember where we have dragged it */ 
            gsc->dragY = y-gsc->elevSy ;
            return 1 ;
	  }
        /* minFig clicked ? - step down */
        else if (x >= gsc->minFigSx && x <= gsc->minFigFx && y >= gsc->minFigSy && y <= gsc->minFigFy ) 
            chg = -gml->cfg.scrlStepSize ;
        /* maxFig clicked ? - step down */
        else if (x >= gsc->maxFigSx && x <= gsc->maxFigFx && y >= gsc->maxFigSy && y <= gsc->maxFigFy ) 
	  chg = +gml->cfg.scrlStepSize ;
	
        else 
	  { /* page up or down ? */
            if (gsc->style&gmlHORZSCRL && y >= gsc->elevSy && y <= gsc->elevFy) 
	      {
                if (x >= gsc->lineSx && x <= gsc->elevSx )  /* page down */
		  chg = -gml->cfg.scrlPageSize*gsc->scale ;
                else if (x >= gsc->elevFx && x <= gsc->lineFx )  /* page up */
		  chg = +gml->cfg.scrlPageSize*gsc->scale ;
	      }
            else if (gsc->style&gmlVERTSCRL && x >= gsc->elevSx && x <= gsc->elevFx) 
	      {
                if (y >= gsc->lineSy && y <= gsc->elevSy )  /* page down */
		  chg = -gml->cfg.scrlPageSize*gsc->scale ;
                else if (y >= gsc->elevFy && y <= gsc->lineFy )  /* page up */
		  chg = +gml->cfg.scrlPageSize*gsc->scale ;
	      }
	  }
        /* scroller position changed ? update it */
        if (chg) gmlScrollSet (gml,gsc,gsc->minVal,gsc->maxVal,gsc->curVal+chg,gsc->scale)  ;
      }
      break ;
    default :break ;
    }    
  
  return 1 ;
}

/* 
    draws the scroller with (dragMode = 0)  or 
    without (dragMode = 1)  the elevator 
*/

static void gmlDrawScroll (gmlEDIT * gml,gmlSCRL * gsc,int dragMode) 
{
  int     box,minFigDim,maxFigDim,oldW ;
  static  float downTriangle[]  = 
  { 0.0,0.0,   1.0,0.0,   0.5,1.0,   0.0,0.0 } ;
  static  float upTriangle[]    = 
  { 0.5,0.0,   1.0,1.0,   0.0,1.0,   0.5,0.0 } ;
  static  float leftTriangle[]  = 
  { 0.0,0.5,   1.0,0.0,   1.0,1.0,   0.0,0.5 } ;
  static  float rightTriangle[] = 
  { 0.0,0.0,   1.0,0.5,   0.0,1.0,   0.0,0.0 } ;
  float   * minFig,*maxFig ;
  
  if (!gsc->style) return ;    
  
  /* determine the figures for arrows to draw */
  if (gsc->style&gmlVERTSCRL)    /* vertical */
    
    {  
      minFig = upTriangle ;minFigDim = vDim (upTriangle) /2 ;
      maxFig = downTriangle ;maxFigDim = vDim (downTriangle) /2 ;
    }
  else if (gsc->style&gmlHORZSCRL)  /* horizontal */
    
    {  
      minFig = leftTriangle ;minFigDim = vDim (leftTriangle) /2 ;
      maxFig = rightTriangle ;maxFigDim = vDim (rightTriangle) /2 ;
    }else return ;
  
  /* draw */    
  box = graphBoxStart ()  ;graphBoxColor (box,gml->cfg.scrlLineColor,TRANSPARENT)  ;
  
  graphColor (gml->cfg.scrlLineColor)  ;
  oldW = graphLinewidth (gml->cfg.scrlLineWidth)  ;
  graphLine (gsc->lineSx,gsc->lineSy,gsc->lineFx,gsc->lineFy)  ;
        
  /* elevator cable - line */
  graphLinewidth (oldW)  ;
  
  /* elevator */
  /*        if (!dragMode) 
	    {*/
  graphColor (gml->cfg.scrlElevatorColor)  ;
  graphFillRectangle (gsc->elevSx,gsc->elevSy,gsc->elevFx,gsc->elevFy)  ;
  graphColor (gml->cfg.scrlLineColor)  ;
  graphRectangle (gsc->elevSx,gsc->elevSy,gsc->elevFx,gsc->elevFy)  ;
  /*        }*/
  
  /* arrows */ 
  _gmlDrawPolygon (minFig,gml->cfg.scrlLineColor,gml->cfg.scrlArrColor,1,minFigDim,gsc->minFigSx,gsc->minFigSy
		   ,gsc->minFigFx-gsc->minFigSx,gsc->minFigFy-gsc->minFigSy)  ;
  _gmlDrawPolygon (maxFig,gml->cfg.scrlLineColor,gml->cfg.scrlArrColor,1,maxFigDim,gsc->maxFigSx,gsc->maxFigSy
		   ,gsc->minFigFx-gsc->minFigSx,gsc->minFigFy-gsc->minFigSy)  ;
  
  graphBoxEnd ()  ;        
  
}





















/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  gmlEDIT Functionality : Graphics Memory Line Editor
_/
_/  This functionality implements smlEDIT functionality 
_/  and graphical interface
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/



static int gmlMagic  =  12345 ;

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Connection to xace graphic system 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
static void gmlLeftDown (double x, double y)  ;
static void gmlLeftDrag (double x, double y)  ;
static void gmlLeftUp (double x, double y)  ;
static void gmlMiddleUp (double x, double y)  ;
static void gmlMiddleDrag (double x, double y)  ;
static void gmlMiddleDown (double x, double y)  ;
static void gmlKeyboard (unsigned char key)  ;

/* creates  an Ace DB interface to gmlEDIT */
static Graph gmlCreateSubGraph (gmlEDIT *gml) 
{
  Graph old  =  graphActive ()  ;
  int graphBusyCursorAll (int)  ;
  
  gml->subGraph  =  graphSubCreate (TEXT_FIT, gml->sx, gml->sy, gml->cx, gml->cy)  ;
  
  /* register events */  
  /* graphRegister (IDLE, gmlIdle)  ;*/
  graphRegister (LEFT_DOWN, gmlLeftDown)  ;
  graphRegister (LEFT_DRAG, gmlLeftDrag)  ;
  graphRegister (LEFT_UP, gmlLeftUp)  ;
  graphRegister (MIDDLE_UP, gmlMiddleUp)  ;
  graphRegister (MIDDLE_DRAG, gmlMiddleDrag)  ;
  graphRegister (MIDDLE_DOWN, gmlMiddleDown)  ;    
  graphRegister (KEYBOARD, gmlKeyboard)  ;
  
  graphAssociate (&gmlMagic, gml)  ; 
  graphActivate (old)  ;
  graphBusyCursorAll (FALSE)  ;
  
  return gml->subGraph ;
}

/* get's gml structure for active graph */
static gmlEDIT * gmlGetGML (void ) 
{
  gmlEDIT * gml ;
  if (graphAssFind (&gmlMagic, &gml)  && gml->magic  ==  &gmlMagic) return gml ;
  return 0 ;
}











/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Constructor and desctructor
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/*  creates an editor with initial content 
    <initTxt> occupying a rectangle sx,sy,sx+cx,sy+cy
    Returns the handle for this editor.
    
    if wrap  == -1 wrap is taken such to fit the given geometry.
    if (wrap == 0)  wrap is taken to 1024.
*/
void * gmlEditorInit (char * initTxt,float sx, float sy, int cx, int cy,int wrap,int style,gmlEDITCONFIG * gConf) 
{
  gmlEDIT *   gml ;
  int         len ;
  
  /* allocate  */
  len = sizeof (gmlEDIT) +sizeof (int) *cy ; 
  gml = (gmlEDIT *) messalloc (len)  ;if (!gml) return 0 ;
  memset ((void*) gml,0,len)  ;
  
  /* init config */
  gmlEditorSetConfig (gml,gConf)  ;
  
  /* set geometry */    
  gml->style = style ;
  gml->sx = sx ;gml->sy = sy ;gml->cx = cx ;gml->cy = cy ;
  gml->textCx = cx ;
  gml->textCy = cy ;
  
  /* aceDB graphics preparations */
  gml->maxBoxesCount = 2000 ; /* we clean up graph and boxes after these many steps */
  gml->magic  =  &gmlMagic ;
  gmlCreateSubGraph (gml)  ;
  
  /* set up the scrollers */
  if ((gml->style&gmlEDIT_VSCROLL)) gml->textCx -= gmlScrollInit (gml,&gml->scrlVert,gmlVERTRSCRL,gml->textCx,0,gml->textCy)  ;
  if ((gml->style&gmlEDIT_HSCROLL)) gml->textCy -= gmlScrollInit (gml,&gml->scrlHoriz,gmlHORZBSCRL,0,gml->textCy,gml->textCx)  ;
  
  /* initialize smlEDIT */
  if (wrap == -1) wrap = gml->textCx-1 ;if (!wrap) wrap = 1024 ;
  gml->sml = smlEditInit (initTxt,wrap," \t")  ;if (!gml->sml) 
    {messfree (gml)  ;return 0 ;}
  
    /* reposition  the scroll */
  gmlEditorSetScrollPos (gml,0,0)  ;
  
  gml->updatedLastLine = gml->textCy ;
  
  return (void *) gml ;
}

/*  
    destroys previously created editor <gmlHandle> object 
*/
void gmlEditorDestroy (void * gmlHandle) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  smlEditFree (gml->sml)  ;
  messfree (gml)  ;
  
}
    












/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Cursor functions 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static int gmlFixSelection (gmlEDIT * gml,int cury,int curx,double * sely,double * selx) 
{   
  int dir ;
    
  /* determine the direction of what comes after what */
  if ( cury < (int)  (*sely)  || 
       ((int) cury  ==  (int)  (*sely)   && curx < (int)  (*selx)) ) 
    dir = 1 ;
  else dir = -1 ;
  
  if (selx) 
    {
      if (dir>0 )  (*selx)  = (double)  (( (int)  (*selx)) +1)  ;
      /*        else if ( cury  ==  (int)  (*sely)  )  (*selx)  = (double)  (( (int)  (*selx)) -1)  ;*/
    }
  return dir ;
}
/* 
   draws XORed  the selection are 
*/
static void gmlDrawSelection (gmlEDIT * gml) 
{
  float sX,sY,fX,fY,tmp ;
  int dir,oneline,maxx ;
  
  /* we skip drawing if no selection or no need to show it */
  if (gml->style& (gmlEDIT_SLINECURSOR|gmlEDIT_MLINECURSOR)) return ;
  if (gml->cursorX == gml->selectX && gml->cursorY == gml->selectY ) return ;
  
  /* the maximum column number */
  maxx = gml->textCx ;
  if (maxx>gml->sml->lineMaxSize) maxx = gml->sml->lineMaxSize ;
  
  /* get and adjust selection coordinates */
  sX = gml->cursorX-gml->scrlHoriz.curVal ;sY = gml->cursorY-gml->scrlVert.curVal ;
  fX = gml->selectX-gml->scrlHoriz.curVal ;fY = gml->selectY-gml->scrlVert.curVal ;
  if (sY>fY) 
    { /* correct the Y order */
      tmp = sX ;sX = fX ;fX = tmp ;
      tmp = sY ;sY = fY ;fY = tmp ;
    }
  if (sY == fY && sX>fX) 
    { /* correct the X order */
      tmp = sX ;sX = fX ;fX = tmp ;
    }
  
  /* determine the direction of what comes after what */
  if ( gml->cursorY < gml->selectY || 
       (gml->cursorY  ==  gml->selectY  && gml->cursorX < gml->selectX)) 
    dir = 1 ;
  else dir = -1 ;

  if (sY == fY) oneline = 1 ;else oneline = 0 ;
  if (dir<1) 
    {sY += 1 ;}
  else 
    {fY += 1 ;} 

  graphColor (gml->cfg.xorBoxColor)  ;
  
  /* draw to vertical lines on sY and fy*/
  graphPatternedXorLine (gml->cfg.selPattern,sX,sY,sX,sY+dir)  ;
  graphPatternedXorLine (gml->cfg.selPattern,fX,fY,fX,fY-dir)  ;
  
  if (oneline) 
    { /* one line selection ? */
      graphPatternedXorLine (gml->cfg.selPattern,sX,sY,fX,sY)  ;
      graphPatternedXorLine (gml->cfg.selPattern,sX,fY,fX,fY)  ;
    }else 
      {
        /* horizontal lines */
        if (dir<1) 
	  {
            graphPatternedXorLine (gml->cfg.selPattern,0,sY,sX,sY)  ;
            graphPatternedXorLine (gml->cfg.selPattern,fX,fY,maxx,fY)  ;
            graphPatternedXorLine (gml->cfg.selPattern,0,fY+1,fX,fY+1)  ;
            graphPatternedXorLine (gml->cfg.selPattern,sX,sY-1,maxx,sY-1)  ;
	  }else 
	    {
	      graphPatternedXorLine (gml->cfg.selPattern,sX,sY,maxx,sY)  ;
	      graphPatternedXorLine (gml->cfg.selPattern,0,fY,fX,fY)  ; 
            graphPatternedXorLine (gml->cfg.selPattern,0,sY+1,sX,sY+1)  ;
            graphPatternedXorLine (gml->cfg.selPattern,fX,fY-1,maxx,fY-1)  ;
	    }
        /* vertical borders */
        if (dir>0) 
	  {
	    graphPatternedXorLine (gml->cfg.selPattern,maxx,sY,maxx,fY-1)  ;
	    graphPatternedXorLine (gml->cfg.selPattern,0,sY+1,0,fY)  ;
	  }else 
	    {
	      graphPatternedXorLine (gml->cfg.selPattern,maxx,sY-1,maxx,fY)  ;
	      graphPatternedXorLine (gml->cfg.selPattern,0,sY,0,fY+1)  ;
	    }
      }
}

/* 
   sets the cursor and selection to a given position 
   after adjsuting it to our buffer content 
*/
static void gmlEditorSetCursor (void * gmlHandle, int cury,int curx, int sely,int selx) 
{
  gmlEDIT *   gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  
  /* adjust the positions */
  smlEditFixPosition (gml->sml,0,&curx,&cury,0)  ;
  smlEditFixPosition (gml->sml,0,&selx,&sely,0)  ;
  
  /* if this is LINESELECT editor */
  if (gml->style& (gmlEDIT_SLINECURSOR|gmlEDIT_MLINECURSOR)) 
    {
      
      if (cury == gml->cursorY)  /* set double click if selection did not change */
	gmlEDITsetCallbackFlag (gml,gmlEDITactionDBLCLK)  ;
      
      /* select the whole lines */
      smlEditMiniMax (gml->sml,&cury,&curx,&sely,&selx)  ;
      curx = 0 ;selx = smlEditLinePtr (gml->sml,sely) ->len-1 ; 
      /* or only a single line */        
      if (gml->style&gmlEDIT_SLINECURSOR) sely = cury ; 
    }
  
  /* any changes in the position ? */
  if ( curx != gml->cursorX || cury != gml->cursorY || 
       selx != gml->selectX || sely != gml->selectY ) 
    {
      
      gmlDrawSelection (gml)  ; /* hide the selection */
      
      /* mark the old lines for redraw */
      smlEditMarkLine (gml->sml,gml->cursorY,gml->selectY,gmlLINETOREDRAW)  ;

      /* set the new cursor */    
      gml->cursorX = curx ;
      gml->cursorY = cury ;        
      
      /* set the new selector */
      smlEditFixPosition (gml->sml,0,&selx,&sely,0)  ;
      gml->selectX = selx ;
      gml->selectY = sely ;
      
      /* mark the new lines for redraw */
      smlEditMarkLine (gml->sml,cury,sely,gmlLINETOREDRAW)  ;
      
      /* draw the selection XORed */
      gmlDrawSelection (gml)  ; 
      
      /* set up a callback for cursor or selection change */
      gmlEDITsetCallbackFlag (gml,gmlEDITchangeCURSEL)  ;
    }        
  
}
















/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Drawing functionality 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
/* 
    draws a text line at given coordinates.
    Clips to the visible window.
*/
static int _gmlDrawTextLine (gmlEDIT * gml,char * src,int len,float x, float y) 
{
  char        ch,cr = 0 ;
  int         crlen ;
  
  if (!len) return 0 ;
  
  /* remove invisible lines */
  if (x+len <= 0)  return 0 ; 
  
  /* line is cut from the left ? */
  if (x<0) 
    {src -= (int)  (x)  ;len += x ;x = 0 ;}
  /* line is cut from the right of visible area ? */
  if (gml->textCx<x+len+1) 
    {len = gml->textCx-x ;}
  
  /* avoid printing more than the length of the line and '\n' */
  ch = src[len] ;src[len] = 0 ; 
  if (src[len-1] == '\n') 
    {crlen = len-1 ;cr = src[crlen] ;src[crlen] = ' ' ;}
  
  /* draw the line */
  graphText (src,x,y)  ;
  
  /* restore cuts */
  if (cr) 
    {src[crlen] = cr ;}
  src[len] = ch ;
  
  return 1 ;
}

/* 
    draw the given (num)  line with at the given 
    coordinates taking care about the selection
    and the cursor
*/
    
static int gmlDrawTextLine (gmlEDIT * gml,int num,int x, int y) 
{
  smlLINE * sl ;
  int     sX,sY,fX,fY,tmp,box ;
  int     pos1 = 0,pos2 = 0,pos3 = 0,pos4 = 0 ;
  char    cursorSym[2] ;
  
  /* not a nempty line */
  if (num >= 0) 
    {
      /* get thew line info pointer */
      sl = smlEditLinePtr (gml->sml,num)  ;if (!sl) return 0 ;
      
      /* get and adjust selection coordinates */
      sX = gml->cursorX ;sY = gml->cursorY ;fX = gml->selectX ;fY = gml->selectY ;
      smlEditMiniMax (gml->sml,&sY,&sX,&fY,&fX)  ;
      
      /* determine selection coordinates */
      pos4 = sl->len ;    
      if (num<sY || num>fY) 
	{}
      if (num>sY && num<fY) 
	{pos2 = 0 ;pos3 = sl->len ;} 
      if (num == sY) 
	{pos2 = sX ;pos3 = sl->len ;}
        if (num == fY) 
	  {pos3 = fX ;}
	
    }
  
  /* start drawing  */
  if (gml->boxInfo[y].num) graphBoxClear (gml->boxInfo[y].num)  ;
  gml->boxInfo[y].num = graphBoxStart ()  ;
  gml->boxesCounter++ ;/*printf ("lin ")  ;*/
    
    
  if (pos2>pos1) 
    { /* before the selection */
      box = graphBoxStart ()  ;graphBoxColor (box,gml->cfg.textColor,gml->cfg.textBgColor)  ;
      gml->boxesCounter++ ;/*printf ("beg ")  ;*/
      
      _gmlDrawTextLine (gml,sl->buf+pos1,pos2-pos1,x+pos1-gml->scrlHoriz.curVal,y)  ;
      graphBoxEnd ()  ;
    }
  if (pos3>pos2) 
    { /* the selection */
      box = graphBoxStart ()  ;graphBoxColor (box,gml->cfg.selectColor,gml->cfg.selectBgColor)  ;
      gml->boxesCounter++ ;/*printf ("sel ")  ;*/
      
      _gmlDrawTextLine (gml,sl->buf+pos2,pos3-pos2,x+pos2-gml->scrlHoriz.curVal,y)  ;
      
      graphBoxEnd ()  ;
    }
  if (pos4>pos3) 
    { /* after the selection */
      box = graphBoxStart ()  ;graphBoxColor (box,gml->cfg.textColor,gml->cfg.textBgColor)  ;
      gml->boxesCounter++ ;/*printf ("end ")  ;*/
      
      _gmlDrawTextLine (gml,sl->buf+pos3,pos4-pos3,x+pos3-gml->scrlHoriz.curVal,y)  ;
      
      graphBoxEnd ()  ;
    }
  
  
  /* clean after the end of the line */
    if (pos4-gml->scrlHoriz.curVal<gml->textCx) 
      {tmp = pos4-gml->scrlHoriz.curVal ;if (tmp<0) tmp = 0 ;}
    else tmp = floor (gml->textCx)  ;
    if (tmp != gml->textCx) 
      {
	graphColor (gml->cfg.bgColor)  ;
        graphFillRectangle (tmp,y,gml->textCx,y+1)  ;
      }
    
    /* draw the cursor */
    if ( pos3 <= pos2 && /* no selection */
	 sl && num == gml->cursorY && (! (gml->style& (gmlEDIT_MLINECURSOR|gmlEDIT_SLINECURSOR)) )  ) 
      {
        box = graphBoxStart ()  ;graphBoxColor (box,gml->cfg.cursorColor,gml->cfg.cursorBgColor)  ;
        gml->boxesCounter++ ;/*printf ("cur ")  ;*/
        
        /* extract character at cursor position  */
        cursorSym[0] = sl->buf[gml->cursorX] ;cursorSym[1] = 0 ;
        if (cursorSym[0] == '\n') cursorSym[0] = '>' ;
        /* draw it */
        _gmlDrawTextLine (gml,cursorSym,1,x+gml->cursorX-gml->scrlHoriz.curVal,y)  ;
	
        graphBoxEnd ()  ;
      }
    
    /* terminate this line drawing */
    graphBoxEnd ()  ;
    graphBoxDraw (gml->boxInfo[y].num,TRANSPARENT,TRANSPARENT)  ;

    
    return 1 ;
}


/* 
   redraws all the flagged lines.
   or if line flag == 0 redraws all.
   if scrollers were repositioned - redraws all.
*/
static void gmlEditorDrawUpdate (void * gmlHandle,int flag) 
{
  int         iL,y,endY ;
  gmlEDIT *   gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  
  
  /* what do we do - if do at all ? */
  if (gml->style&gmlEDIT_DONOTDRAW) return ;
  if (gml->scrlHoriz.status&gmlREPOSSCRL || gml->scrlVert.status&gmlREPOSSCRL) 
    flag = 0 ; /* redraw all lines if scrolls changed the position */
  
  /* if too many boxes - clean up them */
  if (gml->boxesCounter >= gml->maxBoxesCount) 
    {
      gml->boxesCounter = 0 ;
      gmlEditorDraw (gmlHandle)  ;
      /*        printf ("redraw issued \n")  ; */
      return ;
    }
  
  
  /* scan through all the visible lines */
  endY = (gml->scrlVert.curVal+gml->textCy)  ;
  if (endY>gml->sml->lineCnt) endY = gml->sml->lineCnt ;
  /* here we go */ 
  for (iL = gml->scrlVert.curVal,y = 0 ;iL<endY ;iL++,y++) 
    { 
      
      if ( flag &&  /* if flag is asked but this line is not flagged - skip it */
	   (! ((smlEditLinePtr (gml->sml,iL) ->srvFlags) &flag))   ) continue ;
      
      /* draw the line */
      gmlDrawTextLine (gml,iL,0,y)  ;
      
      if (flag)  /* unmark the flag as done */
	smlEditUnmarkLine (gml->sml,iL,iL,flag)  ;
    }
  
  /* erase (draw empty)  rectangle to cover possibly disappeared lines since previous drawing */
  if (iL<gml->updatedLastLine) 
    {
      int box = graphBoxStart ()  ;
      gml->boxesCounter++ ;/*printf ("empt ")  ;*/
      
      graphColor (gml->cfg.bgColor)  ;
      graphFillRectangle (0,y,gml->textCx,gml->updatedLastLine+1)  ;
        
      graphBoxEnd ()  ;
      graphBoxDraw (box,TRANSPARENT,TRANSPARENT)  ;
    }
  gml->updatedLastLine = y ;
  
  
    /* draw the vertical scroller */
  if ((gml->style&gmlEDIT_VSCROLL)  && (gml->scrlVert.status&gmlCHNGSCRL || gml->boxScrlV == 0)) 
    {
      if (gml->boxScrlV) graphBoxClear (gml->boxScrlV)  ;
      gml->boxScrlV = graphBoxStart ()  ;
      gml->boxesCounter += 2 ;/*printf ("scrl V ")  ;*/
      
      gmlDrawScroll (gml,& (gml->scrlVert) ,gml->scrlHoriz.status&gmlONSCRL)  ;
      gml->scrlVert.status &= (0xFFFFFFFF^ (gmlCHNGSCRL|gmlREPOSSCRL))  ;
      
      graphBoxEnd ()  ;
      graphBoxDraw (gml->boxScrlV,gml->cfg.bgColor,gml->cfg.bgColor)  ;
    }
  /* and the horizontal scroller */    
  if ((gml->style&gmlEDIT_HSCROLL)  && (gml->scrlHoriz.status&gmlCHNGSCRL || gml->boxScrlH == 0)) 
    {
      if (gml->boxScrlH) graphBoxClear (gml->boxScrlH)  ;
      gml->boxScrlH = graphBoxStart ()  ;
      gml->boxesCounter += 2 ;/*printf ("scrl H ")  ;*/
      
      gmlDrawScroll (gml,& (gml->scrlHoriz) ,gml->scrlHoriz.status&gmlONSCRL)  ;
      gml->scrlHoriz.status &= (0xFFFFFFFF^ (gmlCHNGSCRL|gmlREPOSSCRL))  ;
      
      graphBoxEnd ()  ;
      graphBoxDraw (gml->boxScrlH,gml->cfg.bgColor,gml->cfg.bgColor)  ;
    }   
}


/*  
    draws editor on the display 
*/
void  gmlEditorDraw (void * gmlHandle) 
{
  int         box ;
  gmlEDIT *   gml ;
  Graph       old ;
  
  gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;

  /* activate the current graph */
  old  =  graphActive ()  ;graphActivate (gml->subGraph)  ;
  graphClear ()  ;
  
  /*clean all the boxes drawn before */
  for (box = 0 ;box<gml->cy ;box++) gml->boxInfo[box].num = 0 ; 
  gml->boxScrlV = 0 ;gml->boxScrlH = 0 ; 
  gml->boxesCounter = 0 ;
  box  =  graphBoxStart ()  ;
  
  /* update all */
  gml->updatedLastLine = gml->textCy ;
  gmlEditorDrawUpdate (gml,0)  ;
  
  /* finalyze */
  /*    graphTextBounds (gml->cx,gml->cy)  ; */
  graphBoxEnd ()  ;
  graphRedraw ()  ; 
  graphActivate (old)  ;
}














/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  EVENT functions 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
/* 
    calls all the marked callbacks 
*/

static int gmlEditorCallbackProcs (gmlEDIT * gml) 
{
  int    iType,mask ;
  
  /* call all neccessary callbacks */
  for (iType = 0 ;iType<gmlEDITlastCALLBACK ;iType++) 
    {
      
      /* the event callback flag is set and there is a registered callback for it ? */
      mask = ((0x00000001) << (iType))  ;
      if ( gml->callBackFlags&mask && gml->callbackFunc[iType]) 
	{
	  
	  /* callback and clear the flag */
	  gml->callbackFunc[iType] ((void *) gml,iType,gml->callbackParam[iType])  ;
	  gml->callBackFlags &= (0xFFFFFFFF^mask)  ;
        }
    }
  return 0 ;
}



#include    "key.h"
static void gmlKeyboard (unsigned char key) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  smlLINE * sl ;
  char    cont[12] ;
  int     curx,cury,selx,sely,nochange,i ;
  int     edSx,edSy,edFx,edFy ;
  int     isSelection = 0,cntLines,typeKey = 0 ;
  char    ch ;
  
  if (!gml) return ;
  
  /* the status of the editor  */
  nochange = 1 ;
  cont[0] = 0 ;cont[1] = 0 ;
  cntLines = gml->sml->lineCnt ;
  
  /* cursor, selection and edition coordinates */
  curx = gml->cursorX ;cury = gml->cursorY ;
  selx = gml->selectX ;sely = gml->selectY ;
  sl = smlEditLinePtr (gml->sml,cury)  ;
  edSx = gml->cursorX ;edSy = gml->cursorY ;edFx = gml->selectX ;edFy = gml->selectY ;
  smlEditMiniMax (gml->sml,&edSy,&edSx,&edFy,&edFx)  ;
  if (selx != curx || sely != cury) isSelection = 1 ;
  else isSelection = 0 ;
  /*printf ("key = %d\n",key)  ;*/
  /* what do we do now ?*/
  switch (key)       
    {
      
      /* cursor movement keys */
      
      /*case LEFT_KEY:*/
    case 3: 
      curx-- ;
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,1)  ;
      gml->cursorWishX = -1 ;
      break ;
      /*case RIGHT_KEY:*/
    case 4: 
      curx++ ;
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,1)  ;
      gml->cursorWishX = -1 ;
      break ;
      /*case UP_KEY:*/
    case 1: 
      cury-- ;
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,0)  ;
      curx = gml->cursorWishX ;
      break ;
      /*case DOWN_KEY:*/
    case 2: 
      cury++ ;
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,0)  ;
      curx = gml->cursorWishX ;
      break ;
      /*case HOME_KEY:*/
    case 5:
      for (ch = -1 ;ch ;) 
	{ /* look for the beginning of the paragraph */
	  curx = -1 ; 
	  ch = smlEditFixPosition (gml->sml,&sl,&curx,&cury,1)  ;
	  if (curx == 0 && cury == 0) 
	    {curx-- ;break ;}
	  if (ch == '\n') 
	    {curx++ ;break ;}
	}
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,1)  ;
      gml->cursorWishX = -1 ;
      break ;
      /*case END_KEY:*/
    case 6:
      for (ch = -1 ;ch ;) 
	{ /* look for the end of the paragraph */
	  curx = sl->len ;
	  ch = smlEditFixPosition (gml->sml,&sl,&curx,&cury,0)  ;
	  if (ch == '\n') 
	    {break ;}
	  cury++ ;
	}
      gml->cursorWishX = -1 ;
      break ;
      
      /*case INSERT_KEY:*/
		/* case 9:
		   if (gml->style&gmlEDIT_OVERWRITE) gml->style &= (0xFFFFFFFF^gmlEDIT_OVERWRITE)  ;
		   else gml->style |= gmlEDIT_OVERWRITE ;
		   
		   break ;
		*/            
    case 23: /* ctrl W */
      gmlEditorSetCurSel (gml,0,0,gml->sml->lineCnt,-1)  ;
      gmlEditorDrawUpdate (gml,-1)  ;
      return ;
      break ;
      
      /* editing keys */
      
    case '\r':case '\n':
      cont[0] = '\n' ;
      cury++ ;curx = 0 ;
      nochange = 0 ;
      break ;
    case '\t':
      for (i = 0 ;i<gml->cfg.tabSize && i<vDim (cont)  ;i++) cont[i] = ' ' ;cont[i] = 0 ;
      curx += i ;
      nochange = 0 ;
      break ;
    case DELETE_KEY:
      cont[0] = 0xFF ;
      if (!isSelection) edFx = edSx+1 ;
      curx = edSx ;cury = edSy ;
      nochange = 0 ;
      break ;
    case BACKSPACE_KEY:
      cont[0] = 0xFF ;
      curx-- ;
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,1)  ;
      if (!isSelection) 
	{
	  edSx = curx ;edSy = cury ;
	  edFx = curx+1 ;edFy = cury ;
	}
      else 
	{
	  curx = edSx ;cury = edSy ;
	}
      
      nochange = 0 ;
      break ; 
      
    default:
      
      /* ascii symbols */
      if (key >= 32 && key<127) 
	{
	  cont[0] = key ;
	  nochange = 0 ;
	  typeKey = 1 ;
	}
      if ((gml->style&gmlEDIT_OVERWRITE)  && edFx == edSx) edFx++ ;  /* overwrite mode */
      gml->cursorWishX = -1 ;
      break ;
    }
  
  
  
  /* content has been changed */
  if (! (gml->style&gmlEDIT_READONLY)  && !nochange) 
    {
      gmlEditorSetBufferReplace (gml,cont,0,edSy,edSx,edFy,edFx)  ;
    }

  if (typeKey) 
    { /* this letter move the cursor by one ? */
      curx++ ;
      smlEditFixPosition (gml->sml,&sl,&curx,&cury,1)  ;
    }
  
  /* change cursor position */
  gmlEditorSetCursor (gml, cury, curx,cury,curx)  ;
  if (gml->cursorWishX == -1) gml->cursorWishX = gml->cursorX ;
  
  /* adjust the scrolls centering the position of insertion */
  if (gml->cursorX<gml->scrlHoriz.curVal || gml->cursorX >= gml->scrlHoriz.curVal+gml->textCx) 
    gmlScrollSet (gml,&gml->scrlHoriz,0,gml->sml->lineMaxSize-gml->textCx,gml->cursorX-gml->textCx/2,gml->textCx)  ;
  if (gml->cursorY<gml->scrlVert.curVal || gml->cursorY >= gml->scrlVert.curVal+gml->textCy) 
    gmlScrollSet (gml,&gml->scrlVert,0,gml->sml->lineCnt-gml->textCy,gml->cursorY-gml->textCy/2,gml->textCy)  ;
  else if (cntLines != gml->sml->lineCnt) 
    {
      gmlScrollSet (gml,&gml->scrlVert,0,gml->sml->lineCnt-gml->textCy,gml->scrlVert.curVal,gml->textCy)  ;
    }
  
    /* callback ? redraw ? */
  gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */
  gmlEditorDrawUpdate (gml,-1)  ;
  
  /*printf ("%x %d '%c'\n",key,key,key)  ;*/
  
  return ;
}


static void gmlLeftDown (double x, double y) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  int     trt ;
  
  /* first try to interpret by scrollers */
  trt = gmlScrollMouseAction (gml,&gml->scrlVert,x,y,LEFT_DOWN)  ;
  if (!trt) trt = gmlScrollMouseAction (gml,&gml->scrlHoriz,x,y,LEFT_DOWN)  ;
  if (!trt) 
    { /* and then by cursor/selector */
      gmlEditorSetCursor (gml,
			  (int)  (y+gml->scrlVert.curVal) , (int)  (x+gml->scrlHoriz.curVal) ,
			  (int)  (y+gml->scrlVert.curVal) , (int)  (x+gml->scrlHoriz.curVal))  ;
      gml->cursorWishX = gml->cursorX ;
    }
  
  gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */
  gmlEditorDrawUpdate (gml,-1)  ;
}



static void gmlLeftDrag (double x, double y) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  int     trt = 0 ;
  
  /* first try to interpret by scrollers */
  if ((gml->style&gmlEDIT_VSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlVert,x,y,LEFT_DRAG)  ;
  if (!trt && (gml->style&gmlEDIT_HSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlHoriz,x,y,LEFT_DRAG)  ;
  if (!trt) 
    { /* and then by cursor/selector */
      
      y += +gml->scrlVert.curVal ;x += gml->scrlHoriz.curVal ;
      gmlFixSelection (gml,gml->cursorY,gml->cursorX,&y,&x)  ;
      
      gmlEditorSetCursor (gml,gml->cursorY,gml->cursorX, (int) y, (int) x)  ;
      gml->cursorWishX = gml->cursorX ;
    }
  
    gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */    
}


static void gmlLeftUp (double x, double y) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  char    * bufr ;
  int     trt = 0 ;
  
  /* first try to interpret by scrollers */
  if ((gml->style&gmlEDIT_VSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlVert,x,y,LEFT_UP)  ;
  if (!trt && (gml->style&gmlEDIT_HSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlHoriz,x,y,LEFT_UP)  ;
  if (!trt) 
    { /* and then by cursor/selector */
      gmlDrawSelection (gml)  ; /* hide the selection */
    }
  
  /* get the  selection buffer and copy it into Xselection clipboard buffer  */
  bufr = smlEditGetBuffer (gml->sml,0,gml->cursorY,gml->cursorX,gml->selectY,gml->selectX)  ;
  if (bufr) 
    {
      graphPostBuffer (bufr)  ; /* write to screen cut/paste buffer */
      messfree (bufr)  ;
    }
  
  gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */    
  gmlEditorDrawUpdate (gml,-1)  ;
}

static void gmlMiddleDown (double x, double y) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  int     trt = 0 ;
  
  /* first try to interpret by scrollers */
  if ((gml->style&gmlEDIT_VSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlVert,x,y,MIDDLE_DOWN)  ;
  if (!trt && (gml->style&gmlEDIT_HSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlHoriz,x,y,MIDDLE_DOWN)  ;
  if (trt) 
    {/* and then by cursor/selector */
      gmlEditorDrawUpdate (gml,-1)  ;
    }
  gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */
  
}
static void gmlMiddleDrag (double x, double y) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  int     trt = 0 ;
  
  /* first try to interpret by scrollers */
  if ((gml->style&gmlEDIT_VSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlVert,x,y,MIDDLE_DRAG)  ;
  if (!trt && (gml->style&gmlEDIT_HSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlHoriz,x,y,MIDDLE_DRAG)  ;
  
  if ((gml->style&gmlEDIT_FASTSCROLL)) 
    {
      gmlEditorDrawUpdate (gml,gmlLINETOREDRAW)  ;
      gml->scrlVert.status &= (0xFFFFFFFF^gmlDRGBOX)  ;
    }
  gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */    
}

static void gmlMiddleUp (double x, double y) 
{
  gmlEDIT * gml  =  gmlGetGML ()  ;
  char    * bufr ;
  int     trt = 0 ;
  
  /* first try to interpret by scrollers */
  if ((gml->style&gmlEDIT_VSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlVert,x,y,MIDDLE_UP)  ;
  if (!trt && (gml->style&gmlEDIT_HSCROLL)) trt = gmlScrollMouseAction (gml,&gml->scrlHoriz,x,y,MIDDLE_UP)  ;
  
  if (!trt && ! (gml->style&gmlEDIT_READONLY)) 
    { /* and then by cursor/selector */
      gmlDrawSelection (gml)  ; /* hide the selector XOR box  */
      bufr = graphPasteBuffer ()  ;
      gmlEditorSetCursor (gml, (int)  (y+gml->scrlVert.curVal) , (int)  (gml->scrlHoriz.curVal+x) , (int)  (gml->scrlVert.curVal+y) , (int)  (gml->scrlHoriz.curVal+x))  ;
      if (bufr) 
	{
	  gmlEditorSetBufferReplace (gml,bufr,0,gml->cursorY,gml->cursorX,gml->selectY,gml->selectX)  ;
	  
	  gmlEditorSetScrollPos (gml,gml->scrlHoriz.curVal,gml->scrlVert.curVal)  ;
        }
    }
  
  gmlEditorCallbackProcs (gml)  ;    /* call all marked callbacks */
  gmlEditorDrawUpdate (gml,-1)  ;
}







/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  public access functions 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/*
_/_/_/_/_/_/_/_/
_/
_/  cursor and selection functions
_/
*/
void gmlEditorGetCurSel (void * gmlHandle,int * cury,int * curx,int * sely, int * selx) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  if (curx) *curx = gml->cursorX ;
  if (cury) *cury = gml->cursorY ;
  if (selx) *selx = gml->selectX ;
  if (sely) *sely = gml->selectY ;    
}

void gmlEditorSetCurSel (void * gmlHandle,int cury,int curx,int sely, int selx) 
{
  char * bufr ;
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  gmlEditorSetCursor (gmlHandle,cury,curx,sely,selx)  ;
  
  /* get the  selection buffer and copy it into Xselection clipboard buffer  */
  bufr = smlEditGetBuffer (gml->sml,0,gml->cursorY,gml->cursorX,gml->selectY,gml->selectX)  ;
  if (bufr) 
    {
      graphPostBuffer (bufr)  ; /* write to screen cut/paste buffer */
      messfree (bufr)  ;
    }
  
}

void gmlEditorHideCurSel (void * gmlHandle) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  
  smlEditMarkLine (gml->sml,gml->cursorY,gml->selectY,gmlLINETOREDRAW)  ;
  gml->cursorX = gml->cursorY = gml->selectX = gml->selectY = -1 ;
  
}


/*
  _/_/_/_/_/_/_/_/
  _/
  _/  scroller Functions
  _/
*/
void gmlEditorGetScroll (void * gmlHandle,int * scrx,int * scry) 
{
  gmlEDIT *gml = (gmlEDIT *) gmlHandle ;
  if (gml)
    {
      if (scry) *scry = gml->scrlVert.curVal ;
      if (scrx) *scrx = gml->scrlHoriz.curVal ;
    }
}

void gmlEditorSetScrollPos (void * gmlHandle,int x,int y) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;
  
  if (gml)
    {
      if ((gml->style & gmlEDIT_VSCROLL)) 
	gmlScrollSet (gml,&gml->scrlVert,0,gml->sml->lineCnt-gml->textCy,y,gml->textCy)  ;
      if ((gml->style & gmlEDIT_HSCROLL)) 
	gmlScrollSet (gml,&gml->scrlHoriz,0,gml->sml->lineMaxSize-gml->textCx,x,gml->textCx)  ;
    }
}

void gmlEditorSetScrollCenter (void * gmlHandle,int x,int y) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;
  
  if (gml)
    {
      gmlScrollSet (gml,&gml->scrlHoriz,0,gml->sml->lineMaxSize-gml->textCx,x-gml->textCx/2,gml->textCx)  ;
      gmlScrollSet (gml,&gml->scrlVert,0,gml->sml->lineCnt-gml->textCy,y-gml->textCy/2,gml->textCy)  ;
    }
}


/*
_/_/_/_/_/_/_/_/
_/
_/  buffer content reading access 
_/
*/
int gmlEditorGetLineCnt (void * gmlHandle,int * pLineCnt) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return 0 ;
  if (pLineCnt) *pLineCnt = gml->sml->lineCnt ;
  return gml->sml->lineCnt ;
}

char * gmlEditorGetBufferSelection (void * gmlHandle,int * bufSize,int startRow,int startCol,int endRow,int endCol) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return 0 ;
  return smlEditGetBuffer (gml->sml,bufSize,startRow,startCol,endRow,endCol)  ;
}

char * gmlEditorGetLinePtr (void * gmlHandle,int startRow) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return 0 ;
  /* adjust the position */
  smlEditMiniMax (gml->sml,&startRow,0,&startRow,0)  ;
  return smlEditLinePtr (gml->sml,startRow) ->buf ;
}

/*
  _/_/_/_/_/_/_/_/
  _/
  _/  buffer content writing access 
_/
*/
void    gmlEditorSetBufferReplace (void * gmlHandle,char * buf,int bufSize,int startRow,int startCol,int endRow,int endCol) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  
  smlEditChangeLine (gml->sml,buf,bufSize,startRow,startCol,endRow,endCol)  ;
  if (gml->sml->editFlags&smlEDITDIRTY) gmlEDITsetCallbackFlag (gml,gmlEDITchangeCONTENT)  ;
  
  /* reset these since they might be invalid */
  gmlEditorSetScrollPos (gmlHandle,gml->scrlHoriz.curVal,gml->scrlVert.curVal)  ;
  gmlEditorSetCursor (gmlHandle,gml->cursorY,gml->cursorX,gml->selectY,gml->selectX)  ;
  
  if (gml->style & gmlEDIT_AUTOREDRAW) 
    gmlEditorUpdate (gmlHandle)  ;
  
}


/*
  _/_/_/_/_/_/_/_/
  _/
  _/  buffer update functions
  _/
*/
void gmlEditorUpdate (void * gmlHandle) 
{
  Graph       old ;
  gmlEDIT *   gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  
  old  =  graphActive ()  ;graphActivate (gml->subGraph)  ;
  
  gml->updatedLastLine = gml->textCy ;    
  gmlEditorDrawUpdate (gmlHandle,-1)  ;
  
  graphActivate (old)  ;
}

int gmlEditorDirty (void * gmlHandle,int isDirty) 
{
  int     oldDirty ;
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return 0 ;
  
  if (gml->sml->editFlags&smlEDITDIRTY) oldDirty = 1 ;else oldDirty = 0 ;

  if (isDirty != -1) 
    {
      if (isDirty) gml->sml->editFlags |= smlEDITDIRTY ;
      else gml->sml->editFlags &= (0xFFFFFFFF^smlEDITDIRTY)  ;
    }
  return oldDirty ;
}

/*
_/_/_/_/_/_/_/_/
_/
_/  callback functions
_/
*/
gmlEDITCallback gmlRegisterCallback (void * gmlHandle,gmlEDITCallbackTYPE functionType,gmlEDITCallback func,void * param) 
{
  gmlEDITCallback oldFun ;
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return 0 ;
  
  oldFun = gml->callbackFunc[functionType] ;
  gml->callbackFunc[functionType] = func ;
  gml->callbackParam[functionType] = param ;
  
  return oldFun ;
}


/*
  _/_/_/_/_/_/_/_/
  _/
  _/  configuration functions
  _/
*/

gmlEDITCONFIG gmlEditDeafultConfig = 
{
  gmlEDIT_VSCROLL ,/* style */
  
  WHITE           ,/* bgColor */
  LIGHTGRAY,BLUE  ,/* selectBgColor,selectColor */
  RED,YELLOW      ,/* cursorBgColor,cursorColor */
  WHITE,BLUE      ,/* textBgColor,textColor */
  YELLOW,BLACK,LIGHTBLUE, /* scrlArrColor,scrlLineColor,scrlElevatorColor */
  BLUE,                /* xorBoxColor */    
  
  4,0.3,          /* tabSize , selPattern */
  2.,1.5,3./8.,   /* scrlArrowLen,scrlWidth ,scrlLineWidth */
  1.,0.5,         /* scrlStepSize,scrlPageSize */
  8./13.          /* propVert_to_Horiz */
} ;

void gmlEditorSetConfig (void * gmlHandle,gmlEDITCONFIG * gConf) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  
  if (!gConf) 
    { gml->cfg = gmlEditDeafultConfig ;return ; }
  
  gml->cfg = *gConf ;
  
  if (gml->cfg.style == -1) gml->cfg.style = gmlEditDeafultConfig.style ;
  if (gml->cfg.bgColor == -1) gml->cfg.bgColor = gmlEditDeafultConfig.bgColor ;
  if (gml->cfg.selectBgColor == -1) gml->cfg.selectBgColor = gmlEditDeafultConfig.selectBgColor ;
  if (gml->cfg.selectColor == -1) gml->cfg.selectColor = gmlEditDeafultConfig.selectColor ;
  if (gml->cfg.cursorBgColor == -1) gml->cfg.cursorBgColor = gmlEditDeafultConfig.cursorBgColor ;
  if (gml->cfg.cursorColor == -1) gml->cfg.cursorColor = gmlEditDeafultConfig.cursorColor ;
  if (gml->cfg.textBgColor == -1) gml->cfg.textBgColor = gmlEditDeafultConfig.textBgColor ;
  if (gml->cfg.textColor == -1) gml->cfg.textColor = gmlEditDeafultConfig.textColor ;
  if (gml->cfg.scrlArrColor == -1) gml->cfg.scrlArrColor = gmlEditDeafultConfig.scrlArrColor ;
  if (gml->cfg.scrlLineColor == -1) gml->cfg.scrlLineColor = gmlEditDeafultConfig.scrlLineColor ;
  if (gml->cfg.scrlElevatorColor == -1) gml->cfg.scrlElevatorColor = gmlEditDeafultConfig.scrlElevatorColor ;
  
  if (gml->cfg.tabSize <= 0) gml->cfg.tabSize = gmlEditDeafultConfig.tabSize ;
  if (gml->cfg.selPattern <= 0) gml->cfg.selPattern = gmlEditDeafultConfig.selPattern ;
  if (gml->cfg.scrlArrowLen <= 0) gml->cfg.scrlArrowLen = gmlEditDeafultConfig.scrlArrowLen ;
  if (gml->cfg.scrlWidth <= 0) gml->cfg.scrlWidth = gmlEditDeafultConfig.scrlWidth ;
  if (gml->cfg.scrlStepSize <= 0) gml->cfg.scrlStepSize = gmlEditDeafultConfig.scrlStepSize ;
  if (gml->cfg.scrlPageSize <= 0) gml->cfg.scrlPageSize = gmlEditDeafultConfig.scrlPageSize ;
  if (gml->cfg.scrlLineWidth <= 0) gml->cfg.scrlLineWidth = gmlEditDeafultConfig.scrlLineWidth ;
  if (gml->cfg.propVert_to_Horiz <= 0) gml->cfg.propVert_to_Horiz = gmlEditDeafultConfig.propVert_to_Horiz ;
  
}


void gmlEditorSetStyle (void * gmlHandle,int style,int doSet) 
{
  gmlEDIT * gml = (gmlEDIT *) gmlHandle ;if (!gml) return ;
  if (doSet) gml->cfg.style |= style ;
  else gml->cfg.style &= (0xFFFFFFFF^style)  ;
  gml->style = gml->cfg.style ;
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  gmlEDIT Functionality : graphical Multi Line Editor
_/
_/  This provides graphical interface to smlEDIT 
_/  
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* styles */
#define gmlEDIT_HSCROLL     0x00000001
#define gmlEDIT_VSCROLL     0x00000002
#define gmlEDIT_FASTSCROLL  0x00000004
#define gmlEDIT_AUTOREDRAW  0x00000008
#define gmlEDIT_MLINECURSOR 0x00000010
#define gmlEDIT_SLINECURSOR 0x00000020
#define gmlEDIT_READONLY    0x00000040
#define gmlEDIT_OVERWRITE   0x00000080
#define gmlEDIT_DONOTDRAW   0x00000100

/* config */
typedef struct _gmlEDITCONFIG{
    int     style;
    
    /* colors */
    int     bgColor;
    int     selectBgColor,selectColor;
    int     cursorBgColor,cursorColor;
    int     textBgColor,textColor;
    int     scrlArrColor,scrlLineColor,scrlElevatorColor;
    int     xorBoxColor;

    /* miscelaneous parameters */
    int     tabSize;
    float   selPattern; /* the lengt of dashes in drag-selection line */

    /* scroller */
    float   scrlArrowLen,scrlWidth,scrlLineWidth;
    float   scrlStepSize,scrlPageSize;
    float   propVert_to_Horiz;
    
    }gmlEDITCONFIG;
extern gmlEDITCONFIG gmlEditDeafultConfig;
void gmlEditorSetConfig(void * gmlHandle,gmlEDITCONFIG * gConf);
void gmlEditorSetStyle(void * gmlHandle,int style,int doSet);



/*  contruction and destruction */
void * gmlEditorInit(char * initTxt,float sx, float sy, int cx, int cy,int wrap,int style,gmlEDITCONFIG * gConf);
void    gmlEditorDestroy(void * gmlHandle); /*  destroys previously created editor <gmlHandle> object */

/* draw and update functions */
int     gmlEditorDirty(void * gmlHandle,int isDirty);
void    gmlEditorDraw(void * gmlHandle); /*  draws editor on the display */
void    gmlEditorUpdate(void * gmlHandle);


/* cursor and selection */
#define gmlEDIT_LASTLINE    (-1)
#define gmlEDIT_LASTCOLUMN  (-1)
void gmlEditorSetCurSel(void * gmlHandle, int cury, int curx,int sely,int selx);
void gmlEditorGetCurSel(void * gmlHandle,int * cury,int * curx,int * selx, int * sely);
void gmlEditorHideCurSel(void * gmlHandle);

/* scroller */
void gmlEditorSetScrollPos(void * gmlHandle,int x,int y);
void gmlEditorSetScrollCenter(void * gmlHandle,int x,int y);
void gmlEditorGetScroll(void * gmlHandle,int * scrx,int * scry);

/* buffer manipulation */
int     gmlEditorGetLineCnt(void * gmlHandle,int * pLineCnt);
char *  gmlEditorGetLinePtr(void * gmlHandle,int startRow);

/* buffer reading */
char *  gmlEditorGetBufferSelection(void * gmlHandle,int * bufSize,int startRow,int startCol,int endRow,int endCol);
#define gmlEditorGetBuffer(_v_gmlHandle,_v_bufSize) gmlEditorGetBufferSelection((_v_gmlHandle),(_v_bufSize),0,0,-1,-1)
#define gmlEditorGetLineBuffer(_v_gmlHandle,_v_cury) gmlEditorGetBufferSelection((_v_gmlHandle),0,(_v_cury),0,(_v_cury),-1)

/* buffer writing */
void    gmlEditorSetBufferReplace(void * gmlHandle,char * buf,int bufSize,int startRow,int startCol,int endRow,int endCol);
#define gmlEditorSetBuffer(_v_gmlHandle,_v_buf,_v_bufSize)  gmlEditorSetBufferReplace((_v_gmlHandle),(_v_buf),(_v_bufSize),0,0,-1,-1)
#define gmlEditorAppendBuffer(_v_gmlHandle,_v_buf,_v_bufSize)   gmlEditorSetBufferReplace((_v_gmlHandle),(_v_buf),(_v_bufSize),-1,-1,-1,-1)
#define gmlEditorCut(_v_gmlHandle,_v_startRow,_v_startCol,_v_endRow,_v_endCol)  gmlEditorSetBufferReplace((_v_gmlHandle),0,0,(_v_startRow),(_v_startCol),(_v_endRow),(_v_endCol))
#define gmlEditorClean(_v_gmlHandle)         gmlEditorSetBufferReplace((_v_gmlHandle),0,0,0,0,-1,-1)

/* callbacks */
typedef enum    {
    gmlEDITchangeCURSEL=0,
    gmlEDITchangeCONTENT,
    gmlEDITchangeHSCROLL,
    gmlEDITchangeVSCROLL,    
    gmlEDITactionDBLCLK,
    gmlEDITlastCALLBACK
} gmlEDITCallbackTYPE;
    
typedef int (* gmlEDITCallback)(void * gmlHandle,gmlEDITCallbackTYPE message,void * param);
gmlEDITCallback gmlRegisterCallback(void * gmlHandle,gmlEDITCallbackTYPE functionType,gmlEDITCallback func,void * param);




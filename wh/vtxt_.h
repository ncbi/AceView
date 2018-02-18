typedef int vFILE ;
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  FILE OS DEPENDENTS 
	_/                        
	*/
#ifdef S_IRGRP  /* windows systems */
  #define vFILEREAD (S_IRUSR|S_IXUSR|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)
  #define vFILEWRITE S_IWUSR
  #define vFILEDIRSEPARATOR "\\"
  #define vdirCreate(vdirnm) mkdir(vdirnm,S_IRUSR|S_IWUSR|S_IXUSR|S_IRGRP|S_IXGRP)
#else /* Unix systems */
  char *getcwd (char *buffer, int maxlen );
  #define vFILEREAD _S_IREAD
  #define vFILEWRITE _S_IWRITE
  #define vFILEDIRSEPARATOR "/"
  #define vdirCreate(vdirnm) mkdir(vdirnm)
#endif
    
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  FILE RELATED CONSTANTS 
	_/                        
	*/
	#define vFILEDEFAULT (vFILEWRITE|vFILEREAD)
	#define vFILEBLOCKSIZE	(16*1024)


vFILE vtxtFileCreate (const char *fileName, int mode)   ;
vFILE vtxtFileOpen (const char *fileName, int flag, int mode)   ;
int vtxtFileRead (vFILE fileHandle, char *Buffer, int BuffSize)   ;
char *vtxtFileReadUntil (vFILE hFile,int * pPos,char * lookFor, int minOfs, AC_HANDLE h) ;
int vtxtFileWrite (vFILE fileHandle, const char *Buffer, int BuffSize) ;
void vtxtFileClose (vFILE fileHandle) ;
int vtxtFileSetPos (vFILE fileHandle, myoff_t offset) ;
int vtxtFileGetPos (vFILE fileHandle) ;
BOOL vtxtFileRemove (char *fileName) ;
int vtxtFileGetLength (vFILE fileHandle) ;
BOOL vtxtFileIsEndReached (vFILE fileHandle) ;
char *vtxtFileMakePath (const char *dnam,  const char *flnm, const char *ext, AC_HANDLE h) ;
int vtxtFileIsReadeable (const char *fileName) ;
char *vtxtFileGetContent (const char *fileName, int *lengthp, AC_HANDLE h) ;
int vtxtFileSetContent (const char *fileName, const char *contentBuf) ;
char *vtxtFileAppend (const char *fileName, char *contentBuf) ;
char *vtxtFileExecution (const char *tmplt, const char *commandLine, const char *inputBuf,const char *stdoutFile, AC_HANDLE h) ;
char *vtxtDirGetCurrent (AC_HANDLE h) ;
char *vtxtDirList (char *dnam, int * pcnt) ;

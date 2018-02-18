#ifndef VSTD_H_DEFINED

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdarg.h>
#include <ctype.h>
#include <regex.h>
#include <sys/ioctl.h>
#include <time.h>


#ifdef WIN32
    #include <process.h>
    #include <Winsock.h>
	#include <io.h>
#else
	#include <limits.h>
	#include <unistd.h>
    #include <pthread.h>
    #include <sys/mman.h> 
    #include <sys/socket.h> 
    #include <netdb.h>
	#if ! defined(LINUX) && !defined(OPTERON) && !defined(ALPHA) && !defined (CYGWIN)
	    #include  <sys/filio.h> 
	#endif
	#include <dirent.h>
#endif


#ifdef __cplusplus
    extern "C" {
#endif

#include "../wvsim/vdef.h"





/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ Memory functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  Basics
	_/                        
	*/


#include "regular.h"


	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/
	_/  miscelaneous gunctionality
	_/
	*/
		int vmemHashMem(void * mem, int len, int bits, int isDiff);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/
	_/  extensible memory
	_/
	*/

	typedef struct tagMEM {
		int curPos,cnt;
		int	piece,curSize;
		int	flags;
		char * buf;
	}vMEX;
	#define	vMEXSETZERO	0x01
	#define	vMEXALIGN	0x02

	int vmexReplace(vMEX * mem,int pozReplace,void * add,int sizeAdd,int sizeDel);
	#define vmexInsert(v_mem,v_pos,v_add,v_sizeAdd) vmexReplace((v_mem),(v_pos),(v_add),(v_sizeAdd),0)
	#define vmexAppend(v_mem,v_add,v_sizeAdd) vmexReplace((v_mem),-1,(v_add),(v_sizeAdd),0)
	#define vmexDelete(v_mem,v_pos,v_sizeDel) vmexReplace((v_mem),(v_pos),0,0,(v_sizeDel))
	#define vmexExpand(v_mem,v_sizeRequired) vmexReplace((v_mem),-1,0,((v_sizeRequired)-((v_mem)->curPos)),0)
	#define	vmexEmpty(v_mem) ( (v_mem)->curPos=(int)(messfree((v_mem)->buf)) )
	#define vmexAppendString(v_mem,v_add) vmexAppend((v_mem),(v_add),strlen(v_add))
	#define	vmexArr(v_mem,v_itm,v_type) (((v_type *)((v_mem)->buf))[(v_itm)])


	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  MEMORY MAPPING FUNCTIONS
	_/                        
	*/
		char * vmemMapFile(char * fileName,unsigned long * sizeMem, int * fileHandle);
		void vmemUnmapFile(char * mapAddress,int size,int fileHandle);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  MEMORY DICTIONARY FUNCTIONALITY
	_/                        
	*/

		typedef struct tagMEG {
			int	dim,cnt,bits;
			int	collisions;
			vMEX tbl;
			vMEX dat;
		}vMEG;
	
		typedef struct tagMEGITM{
			int	hashTo; /* reserved - do not touch */
			int idOfs,idLen;
			int dataOfs,dataLen;
			void * usrPtr;
		}vMEGITM;

		int vmegFindBin(vMEG * meg,int * pIndex,void * id,int lenId);
		int vmegAddBin(vMEG * meg,int * pIndex,void * id,int lenId, void * data,int lenData);

		void vmegEmpty(vMEG * meg);

		#define vmegAdd(v_meg,v_fnd,v_name,v_data,v_lenData) vmegAddBin((v_meg),(v_fnd),(v_name),(strlen(v_name)+1),(v_data),(v_lenData))
		#define vmegFind(v_meg,v_fnd,v_name) vmegFindBin((v_meg),(v_fnd),(v_name),(strlen(v_name)+1))

		#define	vmegItm(v_meg,v_iFnd) (&(((vMEGITM *)(v_meg)->tbl.buf)[(v_iFnd)]))

		#define	vmegId(v_meg,v_iFnd) ((void *)((v_meg)->dat.buf + vmegItm((v_meg),(v_iFnd))->idOfs))
		#define	vmegIdLen(v_meg,v_iFnd) (vmegItm((v_meg),(v_iFnd))->idLen)
		#define	vmegData(v_meg,v_iFnd) ((void *)((v_meg)->dat.buf + vmegItm((v_meg),(v_iFnd))->dataOfs))
		#define	vmegDataLen(v_meg,v_iFnd) (vmegItm((v_meg),(v_iFnd))->dataLen)

		#define	vmegUsrPtr(v_meg,v_iFnd) ((void *)(vmegItm((v_meg),(v_iFnd))->usrPtr))
		


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ String Functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  STRING AND BUFFER RELATED CONSTANTS 
	_/                        
	*/

		#define	_00				"\0\0" 
		#define	_0				"\0" 
		#define vSTRENDLINE     "\r\n"_00
		#define vSTRBLANK       " \t\r\n"_00
		#define vSTRWILDCARD    "*?"_00
		#define vSTRCOMMENT		";!"_00
		#define	vSTRMAXNAME      (1024)
		#define	vSTRCASELO		 (1)
		#define	vSTRCASEUP		 (2)

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  BASIC FUNCTIONS
	_/                        
	*/

		char * vstrGlueList(int cnt, char * firstPtr, ... );   /* CHECK : mem */
	        #define vstrExpand(vrn_str,vrn_add) vmexExpand(&((vrn_str)->mem),((vrn_str)->mem.curPos)+(vrn_add)) 
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  DOUBLE ZERO TERMINATED STRING FUNCTIONS
	_/                        
	*/
		char * vstrNextList(char * cont,int cnt);
		int vstrCntList(char * cont);
		int vstrIsOnList(char * src,char * list,int * numfnd,int isCaseInSensitive);
		int vstr00(char * src);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  STRING SEARCH FUNCTIONS
	_/                        
	*/
		char * vstrSearchSubstring(char * src,char * find,int occurence,char * stopFind,int isCaseInSensitive);
		int vstrCompareUntilSymbol(char * str1,char * str2,char * symblist,int isCaseInSensitive);
		char * vstrStructuredSearch(char * src,char * begin,char * end,int * pst,int * pfn);
		char * vstrSkipWords(char * src, int num,char * separators);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  STRING EDITING FUNCTIONS
	_/                        
	*/
		char * vstrFindReplaceSymbols(char * src,char * dst,int sizeDest,char * find, char * replacement,int maxTags,int isMatch,int isSkipMultiple);
		char * vstrFindReplaceStrings(char * src,char * dst,int sizeDest,char * find, char * replacement,int maxTags,int isCaseInSensitive);
		char * vstrClnMarkup(char * src,char * dst,int sizeDest,char * startTag, char * endTag,char * replacement,int maxTags,int inside,int isCaseInSensitive);
		int vstrCleanEnds(char * src,char * dst,int sizeDest,char * find,int isMatch);
		char * vstrSeparateSubstringFromString(char * src,char * dst,int sizeDest,int nextStp,char * nextSepar, int isCaseInSensitive);
		int vstrCopyUntilSymbol(char * strsrc,char * strdst,int sizeDest,char * symblist);
		int vstrHungarianText(char * TextStr,char * TextDest,int sizeDest,int CaseType,int IsName,int IsRemIntBlanks);
		int vstrCase(char * TextStr,char * TextDest,int sizeDest,int CaseType);
		#define vstrUpperCase(vrn_str) vstrCase((vrn_str),(vrn_str),0,vSTRCASEUP)
		#define vstrLowerCase(vrn_str) vstrCase((vrn_str),(vrn_str),0,vSTRCASELO)
										
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  EXTENDED PRINTF AND SSCANF FUNCTIONS
	_/                        
	*/
		/*
			int 		"i[%i][[=defVal][>numLo][<numUp]][;]"
			double 		"r[%i][[=defVal][>numLo][<numUp]][;]"
			int enum 	"e[=defVal][^val1]^[val2][;]"
			int flag	"f[=defVal][|val1]|[val2][;]"
			int logic	"l[%=0][^lval1][^lval2][;]";
			char *		"t[%s][=defVal][;]"

				def 		- default value
				numLo 		- the lowest acceptible value 
				numHi 		- the highest acceptible value 
				%format 	- format specifications as in printf funcitons
				; 			- additional space line after;
				enu1,enu2	- textual names of values 0,1,2,3,4...
				flg1,flg2	- textual names of values 0x00,0x01,0x02,0x04,0x08
				lval1,lval2	- can be any of ^FALSE^TRUE^OFF^ON^NO^YES^CANCEL^OK^0^1

		*/
		int vstrExtendedSScanf(char * textScan,char * KeyFormat,void * pVal);
		char *  vstrExtendedSPrintf(char * RetBuf,char * KeyFormat,void * pVal);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  DYNAMIC SPRINTF FUNCTIONS
	_/                        
	*/
		int vstrPrintfSizeOfArgList(char * formatDescription,va_list marker);

		char * vstrDynaPrintfArgList(char * formatDescription,va_list marker);
		char * vstrDynaPrintf(char * formatDescription , ...);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  Extensible string buffer Functions
	_/                        
	*/
				   
		typedef struct tagSTR {
			int mode;
			char * outFile;
			void (* funcCallback)(void *,void *,char *); 
			void * paramCallback;
			int	hFile;
			vMEX mem;
		}vSTR;
		
		/* modes for  */
		#define vSTRECHO 0x01
		#define vSTRDISCARD 0x02

		int vstrPrintfArgList(vSTR * str,char * formatDescription,va_list marker);
		int vstrPrintf(vSTR * str,char * formatDescription,...);
		int vstrPrintfWrapped(vSTR * str,char * src,char * separ,int charrayLen,int caseChar, int maxnum);
		#define	vstrEmpty(v_str) (vmexEmpty(&(v_str)->mem))
		#define vstrPtr(v_str) ((v_str)->mem.buf)
		#define vstrLen(v_str) ((v_str)->mem.curPos)


	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  Regular expression functions 
	_/                        
	*/
		/* int vstrSearchRegexp(char * src, char * rexpLine, int cntmatch, int * pmatch,int isCaseInSentitive);*/



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ File I/O functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  FILE OS DEPENDENTS 
	_/                        
	*/
	#ifdef S_IRGRP          /* Unix systems */
	    #define vFILEREAD (S_IRUSR|S_IXUSR|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)
	    #define vFILEWRITE S_IWUSR
    	#define vFILEDIRSEPARATOR "\\"

	#else /* windows systems */
	    char *getcwd( char *buffer, int maxlen );
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

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  DIRECTORY FUNCTIONS
	_/                        
	*/
		char * vdirGetCurrent(char * WorkDir,int maxBuf);
		char * vdirList(char * dnam,int * pcnt);



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ Online Functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

	#define vONLINEBLOCKSIZE	(16*1024)
	#ifndef MAXHOSTNAMELEN
		#define MAXHOSTNAMELEN	1024
	#endif


	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  HTTP/CGI BUFFER MANIPULATION
	_/                        
	*/
		#define vONLINEMAXURLSIZE (16*1024)

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ Sorting Functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
	typedef int (* vSORTFunction)(void * pVal,int i, int j);
	#define	vSORTNUMBER	0x01
	#define	vSORTORDER	0x02

	void vsortQuickCallbackIndex(int toDO,void  * pVal,int n,int * istack,int * ind,vSORTFunction isCondFunc);



/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ Text analysis Functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
		#define	vTEXTMAXLINE		 (16*1024)
		#define	vTEXTMAXFILESIZE	 (1024*256) 

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  Text tree functions 
	_/                        
	*/

		enum	{vTEXTTREELNK=0,vTEXTTREELNKNXT,vTEXTTREELNKPRV};
		typedef struct _vTEXTTREEWRD{
			int		score,cnt;
			vMEX	ref[3]; 
			}vTEXTTREEWRD;

		#define	vTEXTTREEBEGINWRD		0
		#define	vTEXTTREEENDWRD			1
		
		#define vtextTreeWordLnk(v_dic,v_itm)		(((vTEXTTREEWRD * )vmegData((v_dic),(v_itm))))
		#define vtextTreeWordScore(v_dic,v_itm)		(((vTEXTTREEWRD * )vmegData((v_dic),(v_itm)))->score)
		#define vtextTreeWordCnt(v_dic,v_itm)		(((vTEXTTREEWRD * )vmegData((v_dic),(v_itm)))->cnt)
		#define vtextTreeWordRef(v_dic,v_itm,v_rf)	(&(((vTEXTTREEWRD * )vmegData((v_dic),(v_itm)))->ref[v_rf]))

		void vtextTreeLinkAttach(vMEX * ref,int  ix);
		void vtextTreeCreate(char * in,vMEG * dicStc, vMEG * dicWrd,int sortType,vCallbackLogical filterFunc);
		void vtextTreeEmpty(vMEG * dic);
		void vtextTreeSort(vMEG * dic,int lnk,vMEG * dicSort,int sortType);
		void vtextTreePrint(vSTR * out, char * title,vMEG * dic1, vMEG * dic2,int what, int * indSort);

	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  Text Dictionary functions 
	_/                        
	*/

		char * vtextDicSentenceCompose(vMEG * dicWrd,int * pWrdList, int maxwrd, char * pStc, int * pScore);
		int vtextDicSentenceConsume(char * stc,vMEG * dicWrd,int * pCnt,int sizeDat);
                void vtextDicSortIndex(vMEG * dic,int sortType, int * pInd, int maxInd) ;
	/*
	_/_/_/_/_/_/_/_/_/_/_/_/_/
	_/                        
	_/  Text Comparison functions 
	_/                        
	*/

		int vtxtCompareSentences(char * seq1,char * seq2,vMEG * dicRes);

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ Keyword-File analysis Functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

	

#ifdef __cplusplus
	}
#endif
#define VSTD_H_DEFINED
#endif




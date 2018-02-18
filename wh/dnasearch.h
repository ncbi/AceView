#ifndef _dnasearchlib
#define _dnasearchlib


#define PARSEERRORMSGLEN 256

typedef struct {
  int first, last;
  int priority;
  char msg[PARSEERRORMSGLEN];
} DnaLexer_ParseErrorInfo;

#ifdef __cplusplus
extern "C" {
#endif
  extern DnaLexer_ParseErrorInfo dnasearch_parseerrorinfo;
  enum ObjTypes { SOBJ_SHORTAPPROXMATCH, SOBJ_LONGAPPROXMATCH, SOBJ_SPACES, 
		  SOBJ_EPSILON, SOBJ_AND, QSOBJ_COMPLETEQUERY, QSOBJ_SCOREQUERY, QSOBJ_ROOT, SOBJ_POSITION };

  typedef void *DnaSearch_Query;
  typedef int (*DnaSearch_InterruptProc)(int n_bytes_read);
  typedef char *(*DnaSearch_InfoToStringProc)(void *info);
  
  typedef struct {
    int type;   /* if < 0 : global result, else detail of a result. */
    int x1, x2; /* coordinates in the data */
    union {
      struct {
	int p1,p2;  /* coordinates of the part matching in the query string */
	int info;   /* pointer to object-specific data about the match. use xDnaSearch_Result.infos[info]*/
	DnaSearch_InfoToStringProc InfoToString;
      } details;
      float score; /* valid if type < 0 */
    } params;
  } DnaSearch_ResultItem;
  
  typedef struct {
    DnaSearch_ResultItem *items;
    int itemcount;
    int occurences;
    char *infos;
  } DnaSearch_Result;
  
  DnaSearch_Query DnaSearch_CompileQuery(char *query,int optimize);
  void DnaSearch_InitSearch(DnaSearch_Query compiledquery, int savebest);
  void DnaSearch_SearchBlock(DnaSearch_Query compiledquery, char *data, int size);
  void DnaSearch_SearchString(DnaSearch_Query compiledquery, char *data);
  int DnaSearch_SearchBlockInterrupt(DnaSearch_Query compiledquery, char *data, int size, int atomicsize, DnaSearch_InterruptProc testproc);
  DnaSearch_Result *DnaSearch_TerminateSearch(DnaSearch_Query compiledquery);
  void DnaSearch_PrintResult(FILE *output, DnaSearch_Result *result, int  verboselevel);
  void DnaSearch_DestroyResultProc(DnaSearch_Result *result);
  void DnaSearch_DestroyQueryProc(DnaSearch_Query compiledquery);
#define DnaSearch_DestroyResult(__result) {DnaSearch_DestroyResultProc(__result);__result=NULL;}
#define DnaSearch_DestroyQuery(__query) {DnaSearch_DestroyQueryProc(__query);__query=NULL;}
  
#ifdef __cplusplus
}
#endif

#endif

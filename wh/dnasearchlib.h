#ifndef _dnasearchlib
#define _dnasearchlib

#ifdef __cplusplus
extern "C" {
#endif

typedef void *CompiledQuery;
typedef int (*Int_Proc_Int)(int bytesread);

typedef struct {
  int type;  /* if < 0 : global result, else detail of a result. */
  int x1, x2; /* coordinates in the data */
  union {
    struct {
      int p1,p2; /* coordinates of the part matching in the query string */
    } details;
    float score; /* valid if type < 0 */
  } params;
} DnaSearch_ResultItem;

typedef struct {
  DnaSearch_ResultItem *items;
  int itemcount;
  int occurences;
} DnaSearch_Result;

CompiledQuery DnaSearch_CompileQuery(char *query);
void DnaSearch_InitSearch(CompiledQuery compiledquery);
void DnaSearch_SearchBlock(CompiledQuery compiledquery, char *data, int size);
void DnaSearch_SearchString(CompiledQuery compiledquery, char *data);
int DnaSearch_SearchBlockInterrupt(CompiledQuery compiledquery, char *data, int size, int atomicsize, Int_Proc_Int testproc);
DnaSearch_Result *DnaSearch_TerminateSearch(CompiledQuery compiledquery);
void DnaSearch_PrintResult(FILE *output, DnaSearch_Result *result, int  verboselevel);
void DnaSearch_DestroyResult(DnaSearch_Result *result);
void DnaSearch_DestroyQuery(CompiledQuery compiledquery);

#ifdef __cplusplus
}
#endif


#endif

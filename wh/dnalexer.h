#ifndef _dnalexer
#define _dnalexer

#define PARSEERRORMSGLEN 256

typedef struct {
  int first, last;
  int priority;
  char msg[PARSEERRORMSGLEN];
} DnaLexer_ParseErrorInfo;

#ifdef __cplusplus
extern "C" DnaLexer_ParseErrorInfo parseerrorinfo;
#else
extern DnaLexer_ParseErrorInfo parseerrorinfo;
#endif

#ifdef _regexptree
//for internal use
Regexp_Node *DnaLexer_ParseString(char *str);

#endif

#endif

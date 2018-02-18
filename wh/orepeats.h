 
#ifndef OREPEAT_DEF
#define OREPEAT_DEF

typedef struct
{ 
  int phase ;
  Array oLengths, profile ;
} OREPEAT ;

OREPEAT *oligoRepeatRegister (OREPEAT *oligoRepeats, int x1, int x2, KEY key, AC_HANDLE h ) ;
BOOL oligoRepeatDump (FILE* f, Stack s, KEY k) ;
BOOL oligoRepeatParse (int level, KEY key) ;
void oligoRepeatRC (OREPEAT *or, int max) ;
void oligoRepeatDoDestroy (OREPEAT *oligoRepeats) ;
#define oligoRepeatDestroy(_or) (oligoRepeatDoDestroy(_or), (_or)=0)

#endif

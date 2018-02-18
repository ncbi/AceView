#ifndef ACINSIDE_DEF
#define ACINSIDE_DEF
#define ACE_4_7

/*
* Here we complete the incomplete types from ac.h.
*
* Because they are handled via handleAlloc/Free, it is mandatory
* to start each struct on a basic type, not on another allocated
* memory like an Array or and AC_*
*/

struct ac_db 
{
  int magic ;
  AC_HANDLE handle ;
  Stack command_stack ;
  AceCommand look ;
};

struct ac_keyset
{
  int magic ;
  AC_DB db ;
  KEYSET ks ;
};

struct ac_object
{
  int magic ;
  AC_DB db ;
  OBJ obj ;
  KEY key ;
  BOOL isEmpty ;
  char buf25[25] ; /* to store ac_tag_printable result */
};

struct ac_iter
{
  int magic ;
  AC_DB db ;
  int max, curr ;
  KEYSET ks ;
  AC_OBJ ac_obj ;
};


#endif

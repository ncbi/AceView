/* private include file supporting the ac library */

#ifndef AC__DEF
#define AC__DEF

/*
* all of these objects begin with the field
*	int ac_magic
*
* This field is used internally by the library to recognize different
* kinds of objects.
*/

#define MAGIC_BASE 0x41635f30	/* "Ac_0" */

#define	MAGIC_AC_DB		0
#define MAGIC_AC_OBJECT		1
#define MAGIC_AC_ITERATOR	2
#define MAGIC_AC_TABLE		3
#define MAGIC_AC_KEYSET		4

#define MAGIC_AC_MAX		4

/* these functions are private to the library and should not be
 * used in the applications
 */

AC_TABLE ac_subtable (AC_TABLE table, int row, int col, AC_HANDLE handle) ;
AC_TABLE ac_db_empty_table (AC_DB db, int rows, int cols, AC_HANDLE handle) ;
void ac_table_copydown(AC_TABLE table, int row, int col) ;
AC_OBJ ac_key2obj (AC_DB db, KEY key, BOOL fillIt, AC_HANDLE handle) ;
const char *ac_key_name (KEY key) ;



#endif /* AC__DEF */

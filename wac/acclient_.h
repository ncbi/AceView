#include <wh/dict.h>
#include <wh/bitset.h>


/*
* Define each of these according to which types of database
* servers you want support for.  These macros are only used in
* aceclient_XXX.c where XXX is one of the transport types.  They
* are NEVER used anywhere else, so we can keep the code from turning
* into a tangle of ifdefs.
*
*
* Because they are handled via handleAlloc/Free, it is mandatory
* to start each struct on a basic type, not on another allocated
* memory like an Array or and AC_*
*/


/* #define AC_HAVE_RPC		 conventional RPC ace server */
#define AC_HAVE_SOCKET		/* Sanger "socket server" */
#define AC_HAVE_SOCKET_JEAN	/* Sanger "socket server" in Jean's code base*/
#define AC_HAVE_ACETCP		/* Mark S acetcp protocol */

struct ac_db
	{
	  int		magic;
		/*
		* magic number to identify it as AC_DB for ac_free
		*/

          int date, refresh_date ;
	        /*
		 * transaction clock, to know if objects are obsolete
		 */
	  void (*ac_partial_command) (AC_DB db, const char *command, 
			unsigned char **response, int *response_length, 
			unsigned char **response_free, int *encore);
		/*
		* which function to use to issue a command to the database
		*	response is filled in with the answer
		*	response_length is the length
		*	response_free is a value to pass to free() after you are
		*		done with response
		*	*encore is true if there is more data
		*/

	  void (*lazy_command)(AC_DB db, const char *command );
		/*
		* sends a command, but does not return a response.
		* The transport is free to delay the command if
		* it wants, but commands will always be executed
		* in the order they are sent.
		*/

	  void (*close_transport)( AC_DB t );
		/*
		* function to use to close the connection to the database
		*/

 	  void * db_access;
		/*
		* this block of memory is private to a particular type
		* of transport.  It contains file descriptors, context,
		* whatever, that is needed to talk to the database server.
		*/

	  int 	swap_needed;
		/*
		* indicates the relation between the endianness of the
		* server and client.
		*/
	  int keysetId ; 
		/* 
	   	* a way to get unique identidiers for all keysets created
	   	* during the existence of the db handle
	   	*/
	  int 	maxTag;
		/*
		* tag_names[x] is the name of tag number x.  All the names
		* for x >= 0 and x < num_tags exist, even if some may be ""
		*/
	  
	  int nActiveClients ;
	        /*
		 * to see if the server is too buzy
		 */
	    
	  DICT *classeDict ;
	 	 /* the class names */

	  DICT *tag_dict;
		/*
		* tag_dict is a dictionary of all tag names known by this
		* database.
		* 	dictName(db->tag_dict, i ) is the name of tag # i
	   	*/

	  Array flags[256] ;
	  	/*
	   	* per class Array of flags
               	* 0x1: empty obj
	  	*/

		/*
		* a simplistic caching system:
		*
		* cached_names is a dict of object names for each class.
		* cached_objects is an associator that associates the char *name
		*	in the dict with the actual AC_OBJ
		*/
	  DICT *	cached_names[256];
	  Associator 	cached_objects[256];
	};


/*
* about ac_object and ac_real_object:
*
* These two data structures are split because that is the only way to
* make it work with the AC_HANDLEs used by halloc().  Ordinarily, I
* would have all the data in ac_object and just decrement the reference count
* when freeing it.  BUT, the handle library does not understand that -- when
* it frees something, it MUST have some bit of memory to deallocate.  To
* deal with that, we allocate and free a small data structure for each 
* _object_reference_.  Now we name that data structure "ac_object" and
* rename the original to "ac_real_object".
*
* The only significance of this from the user's view will be that the same 
* object can have different addresses:
*	AC_OBJ a,b;
* 	a = ac_get_obj (db,"x","y", 0);
* 	b = ac_get_obj (db,"x","y", 0);
*	assert(a != b);
*/

struct ac_object
	{
	int magic;
	struct ac_real_object *x;
	char buf25[25] ; /* to store ac_tag_printable results */
	};

struct ac_real_object
	{
	  int refcount;
	  int get_date ;
	  AC_DB db ;
	  int filled;
		/*
		* true if the data for this object has been loaded.  In
		* that case, one of table or a_data will have a value.
		*/
	  AC_TABLE table;
		/*
		* table contains all of the data of a type B object.  If
		* table is NULL, then the data has not been loaded yet.  An
		* object that contains no tags has a 0 by 0 table.  Note the
		* significant implication that a type A object has an empty
		* table.
		*/
	  unsigned char *a_data;
	  int a_data_len;
		/*
		* a_data is a pointer to type A data for this object.
		*
		* Note that every object has both type A and type B data 
		* fields but nothing that says which is valid.  Note that the
		* interface definition for this library does not require it.
		* ( Example:  Does a type A object contain tag "foo"?  No, 
		* therefore the library is correct in saying that the tag is 
		* not present.  )
		* 
		*/
	  BitSet hasTag ;
		/*
		* bit(hasTag,N) is true if tag number N is present somewhere
		* in the object data.  Of course, you can only belive this if
		* the table field is set.
		*/
	  Array tagLine ;   
		/* 
		* array of struct tagLine describing where the tag occurs
		* in this object.
		*/
	  KEY key ;
                /*
		 * a unique identifier valid during the life of the AC_DB connection
		 * useful to compare objects constructed in different ways
		 * if to AC_OBJ have the same key, they correspond to
		 * the same object on the AC_DB server
		 */
	  int classe_number;
		/*
		* class of this object as a small integer
		*/
	  const char *classe;	    
		/*
		* class of this object as a name
		*/
	  const char *name;	    
		/* 
		* name of this object 
		*/
	};

/*
* tagLine is used to find a particular tag in the table of an object.
* It only locates the tag that should be found by ac_tag_* functions.
* If the tag appears at multiple locations in the tree, the more
* deeply nested instances do not have their locations saved.
*/
struct tagLine
	{
	int row;
	int col;
	int tag_number;
	const char *tag_name;
	};


/*
* remember the game we play with ac_object and ac_real_object?  We do
* the same thing with ac_keyset.
*/

struct ac_keyset
	{
	int magic;
	struct ac_real_keyset *x;
	};

struct ac_real_keyset
	{
	  int id ;          /* the identifier on the server side */
	  int max ;         /* number of elements in keyset */
	  int refcount;
	  AC_DB db ;
	};

struct ac_iter
	{
	  int magic ;
	  struct ac_real_keyset *real_keyset;
		/*
		* the keyset we are iterating over.
		*/
	  int next ;        
		/* 
		* index of the next object to be returned by ac_next.
		*/
	  int ra_next;
		/*
		* index of the next object in the read-ahead buffer
		*/
	  int ra_max;
		/*
		* number of objects in the read-ahead buffer
		*/
#define AC_ITER_OBJECT_READAHEAD 10
	  AC_OBJ obj[AC_ITER_OBJECT_READAHEAD];
		/*
		* objects that we have read during our read-ahead
		*/
	  int fillhint ;
		/*
		* true means the user suggests that we load object 
		* data during iteration.  false means that the user
		* suggests that we do NOT load object data during
		* iteration.  Obeying this hint can help efficiency.
		*/
	  AC_DB db ;
		/*
		* db descriptor - copied from keyset 
		*/
	  int id ;
		/*
		* keyset id - copied from keyset
		*/
	  int max ;
		/*
		* number of elements in keyset - copied from keyset 
		*/
	  AC_HANDLE user_handle;
		/*
		* handle that the user used to create this iterator.
		* objects returned from ac_next() are on this handle.
		*/
	} ;


/*
* These functions are defined here so we can use them in the test program.
*/
unsigned char *ac_partial_command(AC_DB db, char *command, int *response_length, int *encore);

void ac_dump(unsigned char *cp, int n);


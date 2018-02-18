char *pgetenv (char *a);
	/*
	* returns environment variable a, or "(null)" if
	* that variable does not exist.
	*/

char *logPrefix (char *desc);
	/*
	* returns a string prefix to use when printing to stderr.
	* In any CGI program, use
	*	fprintf(stderr,"%s ...", logPrefix("what_error"), ... 
	*
	* desc must NOT contain any spaces, or it will confuse the
	* log analyzer
	*
	* ALWAYS use this prefix so that the log analyzer can find
	* our messages.
	*
	* The returned string came from messprintf(), so you must copy
	* it if you pass it to any Ace functions.  You can free it
	* with messfree(), but most likely you will exit the cgi
	* shortly after you print it so you don't need to free it.
	*/


/*
* av_db_desc is a description of one line from the db_config file,
* broken out into fields.
*
* If you have a pointer to one of these, you can just free() it.
*/
struct av_db_desc
	{
	char db_version[100];
		/* 
		* name of the magic number in AceView/CHANGE
		* for example, "v9"
		*/
	char db_name[100];
		/*
		* this is AceView's name for the database: worm, human, 31, etc
		*/
	char db_access[100];
		/*
		* this is the name to pass to ac_open to open it
		*/
	char db_species_char[100]; 
		/*
		* this is a single character describing how we treat this
		* species:
		*	h	human
		*	w 	worm
		*/
	char db_species_title[100];
		/*
		* this is the name of the species to use for display
		*/
	char db_title[100];
		/*
		* this is the title of the database to use for display
		*/
	};


struct av_db_desc *
get_db_desc(char *dbname, char **msg);
	/*
	* look up a database and return a description of it.
	*/

/*
* do we need a thing to list the databases?  not yet, but maybe later.
*/

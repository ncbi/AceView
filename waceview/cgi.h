extern void url_decode_inplace(char *s);	
	/* 
	* Decode a URLencoded string in place.  The decoded string
	* length is always <= the encoded length, so it does not
	* matter that it overwrites the original string.
	*
	* Do not try to use this function on string constants.
	*/

extern char *html_encode(const char *s);
	/*
	* encode a string so you can print it in an html document
	* without worrying about whether it my have html tags in it.
	*/

char *url_encode(const char *s);
	/*
	* encode a string so you can use it in a URL
	*/

extern char *cgi_arg(char *name);
	/*
	* retrieve the named parameter. return NULL if not present.
	* obviously, this does not allow you to use the same parameter
	* multiple times (like "x.cgi?a=1&a=2&a=3") so don't do that.
	*/

extern char *cgi_clean_arg(char *name);
	/*
	* same, but sets to blank < > scrip SCRIPT Script
	*/

extern int cgi_arg_exists(char *name);
	/*
	* retrieve the named parameter. return 0 if not present.
	* return 1 if present at least once
	*/

extern char *cgi_cookie(char *name);
	/*
	* retrieve the named cookie.  return NULL if not present.
	*/

extern int cgi_init(char *debugger_file_name);
	/*
	* fetch the cgi parameters - debugger_file_name is the name
	* of a file that stores the cgi input state.  If you are running
	* in a web server, cgi_init() saves the debugging state in
	* that file.  If not, it loads the cgi parameters from that
	* file.
	*
	* returns 1 for success, -1 for malloc error, or 0 for incorrect
	* passing of cgi parameters.
	*/

void
cgi_call(char *name, void *arg, void (*callback)(char *value, void *arg));
	/*
	* calls the callback with the value for each parameter with the
	* given name.  this is so you can use x.cgi?a=1&a=2&a=3
	*/

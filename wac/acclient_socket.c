#if defined(AC_HAVE_SOCKET_JEAN) || defined(AC_HAVE_SOCKET)


#ifdef AC_HAVE_SOCKET_JEAN
#include "wac/acclient_socket_lib.c"
#endif

struct socket_descriptor
{
  AceConnection aceconn_connection ;
};

static void ac_partial_command_socket(AC_DB db, const char *command, 
				      unsigned char **response, 
	int *response_length, 	
	unsigned char **response_free, 
	int *encore)
{
  struct socket_descriptor *d;
  
  d = (struct socket_descriptor *) db->db_access;
  
  AceConnRequest(
		 d->aceconn_connection,
		 command,
		 (void **) response,
		 response_length );
  
  *response_free = *response;
  
  *encore = 0;
}


/*
* the socket server transport does not have lazy commands, so we fake it.
*/
static void ac_lazy_command_socket( AC_DB db, const char *command )
{
  int len, encore;
  unsigned char *response, *r_free;
  ac_partial_command_socket(db, command, &response, &len, &r_free, &encore);
  free(r_free);
}


void ac_close_socket(AC_DB db)
{
  struct socket_descriptor *d;

  d = (struct socket_descriptor *) db->db_access;

  AceConnDisconnect(d->aceconn_connection);
  AceConnDestroy(d->aceconn_connection);
  messfree(d);
}


static char *ac_open_socket(AC_DB db, const char *protocol, const char *host, int port, const char *other )
{
  struct socket_descriptor *d;
  char *user, *password;
  AceConnStatus status ;
  AceConnection connection = NULL ;
  int timeout = 0;

  while (! isdigit((int)*protocol) )
        protocol++;
  if (*protocol == '/')
	protocol++;

  timeout = atoi(protocol);
  if ( timeout != 0 )
    timeout = 300;

  d = (struct socket_descriptor*) halloc(sizeof(struct socket_descriptor), 0);

  /*
  * the database specifier is "acetcp:host:port" or "a:host:port"
  * but thanks to HTTP, people want to write "acetcp://host:port".
  * I will let them.
  */
  while (*host == '/') host++;

  /*
  * a blank host name means the local host
  */
  if (*host == '\0')
    host = "localhost";

  /*
  * pluck the user/password out of the extra info from the database name
  */
  user = NULL;
  if (! other || !other[0] )
	{
	/*
	* the user field was missing or blank; try to get a username/password
	* from a config file.
	*
	* The environment variable ACEDB_USERNAME can either be a username/password
	* directly, or if it starts with a '/' it is a file that contains a list
	* of usernames & passwords for various servers.
	*/
	char *fname;
	fname = getenv("ACEDB_USERNAME");
	if (fname)
		{
		if ( *fname != '/' )
			user = strnew(fname, 0);
		else
			{
			FILE *f;
			f = fopen(fname,"r");
			if (f)
				{
				char b[100], f_host[100], f_user[100], f_port[100];
				while (fgets(b,sizeof(b),f))
					{
					int n;
					char *s;
					if (b[0] == '/')
						continue;
					s = strchr(b,'\n');
					if (s) *s = 0;
					n = sscanf(b,"%s%s%s",f_host, f_port, f_user);
					if (n != 3)
						continue;
					if ((f_host[0] != '*') && (strcmp(host, f_host) != 0))
						continue;
					if ((f_port[0] != '*') && (port != atoi(f_port)))
						continue;
					user = strnew(f_user, 0);
					break;
					}
				fclose(f);
				}
			}
		}
	}
  else
  	user = strnew(other,0);

  if (!user || !user[0])
	return "must give user/password";

  password = strchr(user, '/');
  if (password)
	*password++ = 0;
  else
	password = "";

  /*
   * it takes two calls to create the connection?
   */
  status = AceConnCreate(&connection, host, port, user, password, timeout);
  if (status != ACECONN_OK)
	goto bad;

  status = AceConnConnect(connection) ;
  if (status != ACECONN_OK)
	goto bad;

  d->aceconn_connection = connection;

  db->ac_partial_command = ac_partial_command_socket;
  db->close_transport = ac_close_socket;
  db->lazy_command = ac_lazy_command_socket;

  db->db_access = d;

  messfree(user);

  return NULL;

bad:
  if (user)
	messfree(user);
  if (connection)
	{
	char *error;
	error = connection->last_errmsg;
	 /* free connection */ 
	return error;
	}

  return "cannot contact server";
}

#else

static char *ac_open_socket(AC_DB db, char *protocol, char *host, int port, char *extras )
{
  return "000 Socket Server not supported in this library";
}

#endif

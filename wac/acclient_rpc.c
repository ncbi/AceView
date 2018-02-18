#ifdef AC_HAVE_RPC

#include <wh/aceclient.h>

struct rpc_db_connection
{
  ace_handle *server;
};


static void ac_partial_command_rpc(AC_DB db, const char *command, 
				   unsigned char **response, int *response_length, unsigned char **response_free, int *encore)
{
  int error;
	
  *encore = 0;
	
  /* bug: somebody should look at askServerBinary and try to 
   * see what it means when when you call it with *encore == 3
   */

  error = askServerBinary(
			  ((struct rpc_db_connection *)db->db_access)->server, command,
			  response, response_length,
			  encore, 0);

  *response_free = *response;

  if (error > 0)
    *response = 0;

  return;
}


/*
* the RPC transport does not have lazy commands, so we fake it.
*/
static void ac_lazy_command_rpc( AC_DB db, const char *command )
{
  int len, encore;
  unsigned char *response, *r_free;
  ac_partial_command_rpc(db, command, &response, &len, &r_free, &encore);
  free(r_free);
}

static void ac_close_rpc(AC_DB db)
{
  struct rpc_db_connection *r;
  r = db->db_access;
  closeServer(r->server);
  messfree(r);
}

static char *ac_open_rpc(AC_DB db, const char *protocol, const char *host, int port )
{
  struct rpc_db_connection *r;
  int timeout = 0;

  while (! isdigit((int)*protocol) )
        protocol++;
  if (*protocol == '/')
	protocol++;

  timeout = atoi(protocol);
  if ( timeout != 0 )
    timeout = 300;

  r = halloc(sizeof(struct rpc_db_connection),0);

  while  (*host == '/') host++;

  if (*host == '\0')
    host = "localhost";

  db->ac_partial_command = ac_partial_command_rpc;
  db->close_transport = ac_close_rpc;
  db->lazy_command = ac_lazy_command_rpc;

  db->db_access = r;
  r->server = openServer(host, port, timeout ) ;

  if (!r->server)
    return "000 database open failed";
  return 0;
}

#else

static char *ac_open_rpc(AC_DB db, const char *protocol, const char *host, int port )
{
  return "000 RPC not supported in this library";
}

#endif

#ifdef AC_HAVE_ACETCP

#include <netinet/tcp.h>
#include <wtcp/tcp_connect.h>

struct acetcp_descriptor
{
  int	fd;
  ITR in ;
  OTR out ;
  char *server_challenge;
};

static unsigned char *acetcp_response = 0, *acetcp_rp;
static int acetcp_response_length;


/*
* This callback processes a single transaction that was just
* received.  In this case, we save it in some global variables
* so that ac_partial_command_acetcp() can look at the value.
*
* incoming_transaction() cannot just return the data that it
* received because in some contexts it may receive more than
* one command.  In this case, we limit it to just one, so we
* know this trick of storing the result in a global is safe.
*/
void acetcp_process_incoming( char *s, int n, void *callback_cookie )
{
  acetcp_response_length = n-11;
  acetcp_response = (unsigned char *) malloc(n+1);
  if (! acetcp_response)
    messcrash("malloc failed");
  acetcp_rp = acetcp_response + 11;	/* skip length, point at type */
  memcpy(acetcp_response, s, n);
}


/*
* send a command that expects no response
*/
static void ac_lazy_command_acetcp(AC_DB db, const char *command)
{
  struct acetcp_descriptor *d;

  d = (struct acetcp_descriptor *) db->db_access;

  outgoing_transaction_send( d->out, 'n', command, strlen(command)+1);
}


/*
* acetcp variant of ac_partial_command - this is compatible with
* all the other transport's partial_command functions
*/

static void ac_partial_command_acetcp(AC_DB db, const char *command, 
				      unsigned char **response, 
					int *response_length, 
			unsigned char **response_free, 
	int *encore)
{
  struct acetcp_descriptor *d;

  d = (struct acetcp_descriptor *) db->db_access;

  * encore = 0;

  /*
  * recognize the special command "encore" and make it ask for an
  * 'encore' transaction instead of a command transaction.  This is a
  * yucky way to recognize it, but it is compatible with all the other 
  * transports.
  */
  if (strcmp(command,"encore") == 0)
  	outgoing_transaction_send( d->out, 'e', command, strlen(command)+1);
  else if (strcasecmp(command,"comment") == 0)
  	outgoing_transaction_send( d->out, 'n', command, strlen(command)+1);
  else
  	outgoing_transaction_send( d->out, 't', command, strlen(command)+1);

  acetcp_response = 0;

  /*
  * eat up some data - keep trying to read responses until we get either
  * EOF or some actual data
  */
  while (( incoming_transaction(d->in, d->fd, 1, acetcp_process_incoming, db->db_access ) == 0)
	 && (acetcp_response == 0))
    ;

  if (command && !strcasecmp (command, "quit"))   /* 214_08_04 trying to remove a non reproducible parasite messerror on quit */
    return ;
  if (command && !strncasecmp (command, "shutdown", 8))   /* 214_08_04 trying to remove a non reproducible parasite messerror on quit */
    return ;
  if (! acetcp_response)
    messcrash("no response from server - dropped connection?");

  /* 
   * if you implement authentication, you would check it here and
   * move acetcp_rp past the authentication block.
   */

  /*
   * switch on transaction type
   */
  switch (* acetcp_rp)
    {
    case 'C':
      /*
       * We do not expect the server to issue multiple challenges.
       */
      messcrash("server challenged us again");

    case 'a':
      /*
       * The server sends us an auth type message when it changes the
       * type of authentication it uses.  If you implement authentication,
       * you would switch the auth type here and go back to process more
       * messages.
       */
      messcrash("server selected auth type - not implememented in preliminary library ");

    case 'E':
      /*
       * This is a major error about the transport, not something like
       * "found 0 objects" or "unrecognized command".  This is no
       * continuing after this.
       */
      messcrash("server sent back error - %s",acetcp_rp+1);

    case 'R':
      /*
       * a real response.  check it for validity and 
       */

      acetcp_rp++;		/* skip transaction type */
      *encore = *acetcp_rp++;	/* fetch encore value (only R has encore) */
      acetcp_response_length -= 2;

      *response = acetcp_rp;
      *response_length = acetcp_response_length;
      *response_free = acetcp_response;
      acetcp_response = 0 ;

      return;

    default:
      messcrash("unrecognized transaction type %02x %c",*acetcp_rp, *acetcp_rp);
    }

  return;	/* not reached, but compiler doesn't know about messcrash */
}


static void ac_close_acetcp(AC_DB db)
{
  struct acetcp_descriptor *d;
  d = (struct acetcp_descriptor *) db->db_access;
  messfree(d->server_challenge);
  incoming_transaction_free(d->in);
  outgoing_transaction_free(d->out);
  close(d->fd);
  messfree (db->db_access) ;
  if (acetcp_response)
    free (acetcp_response) ;
}

/*
* connect to the acetcp server
*/

static char *ac_open_acetcp(AC_DB db, const char *protocol, const char *host, int port )
{
  struct acetcp_descriptor *d;
  static char dyn_error[200];
  int timeout = 0;

  while (! isdigit((int)*protocol) )
	protocol++;
  if (*protocol == '/')
	protocol++;

  timeout = atoi(protocol);
  if ( timeout == 0 )
    timeout = 300;
  /* unfortunately, we don't use the timeout in acetcp */

  d = (struct acetcp_descriptor*) halloc(sizeof(struct acetcp_descriptor), 0);

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

  d->fd = connect_socket(host, port, TRUE);

  if (d->fd < 0)
    return "000 connect to database server failed";

  d->in = incoming_transaction_create();
  d->out = outgoing_transaction_create(d->fd);

  db->db_access = d;

  db->close_transport = ac_close_acetcp;
  db->ac_partial_command = ac_partial_command_acetcp;
  db->lazy_command = ac_lazy_command_acetcp;

  /*
   * consume the challenge that the server will send us on connect.  if
   * you implement authentication, store the challenge to use in later
   * packets.
   */
  acetcp_response = 0;

  while (( incoming_transaction(d->in, d->fd, 1, acetcp_process_incoming, db->db_access ) == 0)
	 && (acetcp_response == 0))
    ;

  if (! acetcp_rp)
    return "server disconnect";

  switch (* acetcp_rp)
    {
    case 'C':
      {
	int k, nActiveClients = 0 ;
	
	d->server_challenge = strnew((char *)acetcp_rp,0);
	k = sscanf (d->server_challenge,"C%d-", &nActiveClients) ;
	if (k >= 1)
	  db->nActiveClients = 1 + nActiveClients ;
	}
      break;
    case 'E':
      strncpy(dyn_error,(char *)acetcp_rp+1, 199);
      dyn_error[199] = 0;
      return dyn_error;
    default:
      return "server sent back unrecognizable response";
    }

  return 0;
}

#else

static char *ac_open_acetcp(AC_DB db, const char *protocol, const char *host, int port )
{
  return "000 acetcp not supported in this library";
}

#endif

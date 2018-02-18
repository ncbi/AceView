#include <unistd.h>
#include <sys/time.h>
#include <sys/uio.h>

#include <wh/regular.h>
#include <wh/array.h>

/*
 * a record for gathering incoming transactions
 */
struct incoming_transaction_rec
{
  int size;
  int awaiting;
  int arraysize;
  int eofcount;
  Array data;
  struct authentication_type *authentication;
};

/*
 * a record for constructing outgoing transactions
 *
 */
struct outgoing_transaction_rec
{
  int	fd;
  struct authentication_type *authentication;
};

#include <wtcp/tcp_connect.h>
BOOL TCPWEBSERVER = FALSE ;

/*
 * an authentication record - How are individual transactions checked
 * to see that they come from who we think they are from.  This is not
 * related to whether the user gave a valid username/password or whether
 * they connected from an authorized host.  
 */

enum authentication_list
  { 
    auth_none, 
    /*
     * no authentication beyond that it comes on the
     * right socket.
     */
    auth_user_md5
    /*
     * MD5 signature authentication.
     */
  };

struct authentication_type
{
  enum authentication_list type;
  union
  {
    struct {
      char *challenge;
      char *password;
    }
    md5_info;
  } x;
};

struct authentication_type * new_auth()
{
  struct authentication_type *a;
  a = (struct authentication_type *) halloc(sizeof(struct authentication_type), 0);
  a->type = auth_none;
  return a;
}

void set_md5_authentication_info( struct authentication_type *a, char *challenge, char *password)
{
  a->x.md5_info.challenge = strnew(challenge, 0);
  a->x.md5_info.password = strnew(password, 0);
  a->type = auth_user_md5;
}


/*
*****
* incoming data
*****
*/


ITR incoming_transaction_create(void)
{
  struct incoming_transaction_rec *t 
    = (struct incoming_transaction_rec *) halloc(sizeof(struct incoming_transaction_rec), 0);
  t->size = 0;
  t->awaiting = 0;
  t->arraysize = 0;
  t->eofcount = 0;
  t->data = arrayCreate(1000, char);
  t->authentication = new_auth();
  return (ITR) t;
}

void incoming_transaction_free (ITR t)
{
  arrayDestroy(t->data);
  messfree(t->authentication);
  messfree(t);
}





/* weird, was is this function not used ?
static int xttime = 0;
static int p_read(int x, void *y, int z)
{
  int n;
  struct timeval tvs, tve;
  gettimeofday(&tvs,0);
  n = read(x,y,z);
  gettimeofday(&tve,0);
  xttime += (tve.tv_sec -tvs.tv_sec ) * 1000000 + tve.tv_usec - tvs.tv_usec;
  return n;
}
*/

int incoming_transaction(ITR t, int fd, int limit_one,
			 void (*process_incoming)(char *s, int n, void *callback_cookie), void *callback_cookie )
{

#define BUFFER_SIZE 49152
  int n;

  /*
   * extend the buffer if necessary
   */
  if (t->size + BUFFER_SIZE > t->arraysize)
    {
      t->arraysize += BUFFER_SIZE;
      (void) array(t->data, t->arraysize, char);
      /* printf("grow\n"); */
	
    }

  /*
   * read some data directly into the buffer
   */
  n = read(fd, & array( t->data, t->size, char), BUFFER_SIZE);
  if (n <= 0)
    {
      if (t->eofcount++ > 3)
	return 1;
      else
	return 0;
    }
  /*
   * printf("socket read: (%d) %d\n",fd,n);
   * ac_dump( & array( t->data, t->size, char), n);
  */
  t->size += n;

  t->eofcount = 0;

  /*
   * Now that we have some data, see if we have a complete
   * transaction.  As long as we do, process them.
   */
  for (;;)
    {
      unsigned char *start;

      /*
       * a transaction has at least 3 bytes of length at the
       * begining
       */
      if (t->size < 3)
	return 0;

      /*
       * get a pointer so we don't need to keep calling the
       * array functions.
       */
      start = arrp( t->data, 0, unsigned char );
	
      /*
       * If t->awaiting is 0, we don't know how many bytes we are
       * expecting.  The length is at the front of the buffer, so
       * we can just pick it up.
       */
      if (t->awaiting == 0)
	{
	  t->awaiting = start[0] | ( start[1] << 8 ) | ( start[2] << 16 );
	}

      /*
       * if the size of the data we have is smaller than the length
       * of the transaction, we can't process it.  We return, and
       * when more data is available, we will be called again to 
       * again read from the buffer.
       */
      if (arrayExists (t->data) && ! strncmp(arrp(t->data, 0,char), "GET ", 4))
	t->awaiting = t->size = t->size + 1 ;
      if (t->size < t->awaiting)
	return 0;

      /*
       * Ok, there is a whole transaction ready.  Process it and then
       * copy rest of the buffer down to the beginning.
       *
       * You could try to process multiple transactions in place before 
       * copying the remainder of the buffer.  That might save some
       * time, but I wouldn't expect much payoff for the extra work and
       * risk of bugs.
       */

      (process_incoming)(arrp(t->data, 0, char), t->awaiting, callback_cookie);

      memmove( arrp(t->data, 0, char), 
	       arrp(t->data, t->awaiting, char ),
	       t->size - t->awaiting );

      t->size -= t->awaiting;
      t->awaiting = 0;

      /*
       * a way to limit the number of transactions processed per call
       * to this function:
       */
      if (limit_one)
	return 0;
    }
}


/*
*****
* outgoing tcp transactions
*****
*/


void outgoing_transaction_free(struct outgoing_transaction_rec *t)
{
	messfree(t->authentication);
 	messfree(t);
}

struct outgoing_transaction_rec *
outgoing_transaction_create(int fd)
{
  OTR t ;
  t = (OTR) halloc(sizeof(struct outgoing_transaction_rec), 0);
  t->fd = fd;
  t->authentication = new_auth();
  return t;
}


int outgoing_transaction_send (OTR t,int type, const char *message, int message_len )
{
  int total_len;
  char b[4];
  total_len = message_len + 3 + 1;
  b[0] = total_len;
  b[1] = total_len >> 8;
  b[2] = total_len >> 16;
  /* no signatures yet. */
  b[3] = type;
  /*
   * for now, I'm just going to write the message and block on
   * the write if it is too big.  That means that if the write
   * fails, we have a broken pipe.  (Note that we ignore sigpipe
   * in wrpc/tcp_server_hooks.c)
   */
#define HAVE_WRITEV
#if defined(HAVE_WRITEV)
  {
    /*
     * writev allows us to write all the different segments
     * without copying it into a single buffer
     */
    struct iovec iov[2];
    int n;
    iov[0].iov_base = b;
    iov[0].iov_len = 4;
    iov[1].iov_base = (char *)message;
    iov[1].iov_len = message_len;
    if (! TCPWEBSERVER && type != 'G')
      {
	n = writev(t->fd, iov, 2);
	if (n != message_len + 4)
	  return 1;
      }
    else if (TCPWEBSERVER && type == 'G')
      {
	int n1 ;
	char crlf[3] ;
	crlf[0] = 0x0d ; 
	crlf[1] = 0x0a ; 
	crlf[2] = 0 ;
	char *cp = messprintf ("HTTP/1.1 200 OK%s"
			       "Server: acedb%s"
			       "Accept-Ranges: bytes%s"
			       "Content-Length: %d%s"
			       "Connection : close%s"
			       "Content-type: text/html; charste=UTF-8%s"
			       "%s"
			       , crlf, crlf, crlf, message_len, crlf, crlf, crlf, crlf
			       ) ;
	n1 = strlen(cp) ;
	if (1 && write(t->fd, cp, n1) != n1)
	  return 1 ;
	if (write(t->fd,message, message_len) != message_len)
	  return 1;
      }
  }
#else
  /*
   * this is what you do if your system doesn't have writev
   */
  if (! TCPWEBSERVER && type != 'G')
    {
      if (write(t->fd,b,4) != 4)
	return 1;
    }
  else
    {
      char crlf[3] ;
      crlf[0] = 0x0d ; 
      crlf[1] = 0x0a ; 
      crlf[2] = 0 ;
      char *cp = messprintf ("HTTP/1.1 200 OK%s"
			     "Server: acedb%s"
			     "Accept-Ranges: bytes%s"
			     "Content-Length: %s"
			     "Connection : close%s"
			     "Content-type: text/html; charste=UTF-8%s"
			     "%s"
			     , crlf, crlf, crlf, message_len, crlf, crlf, crlf, crfl
			     ) ;
	n1 = strlen(cp) ;
	if (1 && write(t->fd, cp, n1) != n1)
	  return 1 ;
    }

  if (write(t->fd,message, message_len) != message_len)
    return 1;
#endif

  return 0;
}


#if 0
<!--
self.moveTo(0,0)
     self.resizeTo(screen.availWidth,screen.availHeight)
     //-->

#endif

     

#include "regular.h"
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <netinet/tcp.h>
#include <sys/ioctl.h>
#if !defined(LINUX) && !defined(OPTERON) && !defined (CYGWIN) && !defined (ALPHA)
 #include <sys/filio.h>
#endif

int connect_socket (const char *host, int port, BOOL no_delay)
{
  int fd, n, one = 1 ;
  struct sockaddr_in there;
  struct hostent *h;
  
  h = gethostbyname(host);
  if (!h)
    {
      fprintf(stderr,"no such host\n");
      return -1;
    }
  
  fd=socket(AF_INET, SOCK_STREAM, 6);
  /*
   * Strictly speaking, you should use getprotobyname to find
   * the protocol number for TCP, but only a moron would design
   * a system that doesn't use 6 for TCP.
   */
  
  if (fd < 0)
    {
      perror("socket");
      return -1;
    }
  
  n=49152;
  setsockopt(fd, SOL_SOCKET, SO_SNDBUF, &n, sizeof(n));
  setsockopt(fd, SOL_SOCKET, SO_RCVBUF, &n, sizeof(n));
  
  memset(&there,0,sizeof(there));
  there.sin_family = AF_INET;
  there.sin_port = htons(port);
  memcpy(&there.sin_addr.s_addr, h->h_addr_list[0],sizeof(there.sin_addr.s_addr));
  
  if (connect(fd, (struct sockaddr *)&there, sizeof(there)) < 0)
    {
      perror("connect");
      return -1;
    }
  
  /*
   * it is important to set TCP_NODELAY on the socket.  if you
   * do not, transactions will be VERY slow.  Without TCP_NODELAY,
   * the system hopes you will write more data to the socket
   * very soon, so it does NOT transmit the data immediately -- it
   * wants to put the results of your next write in the same
   * packet.
   *
   * 6 is what getprotobyname returns for "tcp" - technically, I'm
   * cheating but TCP is _always_ protocol number 6 because that
   * is the IP protocol number for TCP.
   */
  if (no_delay && setsockopt (fd, 6, TCP_NODELAY , &one, sizeof(int)) < 0)
    {
      /*
       * this setsockopt can never fail, so it is ok to crash if it does.
       */
      if (fd)
	{ close (fd) ; return -1 ; }
    }

  return fd;
}

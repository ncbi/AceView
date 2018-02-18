#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "../wtcp/tcp_connect.h"

extern int listen_socket(int port) ;
/* extern void exit (unsigned int s) ; */
#if defined(ALPHA)
typedef unsigned long   socklen_t;      /* 64-bits */
#endif

void daemon_main(int daemon_port, int daemon_fork, int (*daemon_service)(int) )
{
int fd;
int s;
socklen_t n;
int ecount;

fd = listen_socket(daemon_port);
ecount = 0;

for (n = 1; n < daemon_fork; n++)
	{
	s = fork();
	if (s < 0)
		perror("fork");
	else if (s == 0)
		break;
	}

/* start */
for (;;)
	{
	struct sockaddr_in there;
	n = (socklen_t) sizeof(there);
	s = accept(fd, (struct sockaddr *) &there, &n);
	if (s < 0)
		{
		perror("accept");
		if (ecount++ > 5)
			exit(1);
		}
	else
		{
		ecount = 0;
		(*daemon_service)(s);
		close(s);
		}
	}
}

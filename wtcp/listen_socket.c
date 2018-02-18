#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stdio.h>
#include <string.h>

int listen_socket(int port)
{
int fd,n;
struct sockaddr_in here;
	
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

n=1;
setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &n, sizeof(n));

n=49152;
setsockopt(fd, SOL_SOCKET, SO_SNDBUF, &n, sizeof(n));
setsockopt(fd, SOL_SOCKET, SO_RCVBUF, &n, sizeof(n));

memset(&here,0,sizeof(here));
here.sin_family = AF_INET;
here.sin_port = htons(port);
here.sin_addr.s_addr = INADDR_ANY;

if (bind(fd, (struct sockaddr *)&here, sizeof(here)) < 0)
	{
	perror("bind");
	return -1;
	}

if (listen(fd, 5) < 0)
	{
	perror("listen");
	return -1;
	}

return fd;
}


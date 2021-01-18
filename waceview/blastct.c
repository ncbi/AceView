#include "regular.h"
#include <wtcp/tcp_connect.h>
#include <wsocket7/acesocket_.h>

int main (int argc, char **argv)
{
  int fd;
  char buf[1000];
  char blastd_args[100];
  char sequence[1000] ;
  FILE *fp;
  char *blast_host = "localhost";
  
  if (argv[1])
    blast_host = argv[1];
  
  setbuf(stdout,NULL);
  
  printf("enter args: ");
  fgets(blastd_args, 99, stdin);
  
  printf("enter sequence: ");
  fgets(sequence, 999, stdin);
  
  fd = connect_socket(blast_host,3001, FALSE);
  
  if (fd < 0)
    {
      printf("<h1>Blast search server is down</h1>\n");
      printf("This is unexpected.  Try again later.\n");
      exit(1);
    }
  
  /*
   * this interface with the server is documented in README.  If you change
   * it, update that file.  The server is in blastd.c.
   */
  write(fd,blastd_args,sizeof(blastd_args));
  
  write(fd,sequence,strlen(sequence));
  write(fd,"\n",1);
  shutdown(fd, 1);
  
  fp = fdopen(fd,"r");
  
  while (fgets(buf,sizeof(buf),fp))
    printf("%s",buf);
  
  exit(0);
}

/*  $Id: test.client.c,v 1.1.1.1 2002/07/19 20:23:32 sienkiew Exp $  */
/*  Last edited: Nov 25 19:59 1992 (mieg) */
/*
** this program send message to server
*/

#include <stdio.h>
#include <string.h>
#include "nace_com.h"


void main( argc, argv )
     int argc;
     char *argv[];
{
  net_data data_s;   /* data for server */
  net_data *data_c;  /* result from server*/
  int   erc;
  int   time_out;
  CLIENT *cl;
  
  
  
  if ( argc < 3 ) 
    {
      fprintf(stderr, "usage : %s host message  [time_out]\n",argv[0]);
      exit(1);
    }
  

  bzero( (char *)&data_s, sizeof(data_s) );

  data_s.cmd.cmd_len = (u_int)strlen(argv[2])+1;
  data_s.cmd.cmd_val = (u_char *)strdup(argv[2]); 
  if ( data_s.cmd.cmd_val == NULL ) 
    {
      perror ("strdup");
      exit(1);
    }
  if ( argc > 3 ) 
    {
      sscanf(argv[3],"%i",&time_out);
    }
  else
    {
      time_out = 0;
    }
  
  cl = open_ace_server (argv[1], NET_ACE , NET_ACE_VERS );
  if (cl == NULL)
    exit (1);
  erc = call_ace_server (cl,&data_s,&data_c,time_out);
  if ( erc < 0) 
    {
      fprintf(stderr, "call_ace_server return with : %i\n",erc);
      exit(1);
    }
  printf("server_message : len = %i\n"
         "               : val = %s\n",data_c->cmd.cmd_len, data_c->cmd.cmd_val);

  close_ace_server(cl);
  exit(0);
}


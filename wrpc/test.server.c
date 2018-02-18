/*  $Id: test.server.c,v 1.1.1.1 2002/07/19 20:23:32 sienkiew Exp $  */
/*  Last edited: Mar 24 12:51 1993 (mieg) */

/*
** simple server which make  toupper() and return
** message back.
*/

#include <errno.h>
#include <string.h>
#include "nace_com.h"

int r_toupper ( data_in , data_out )
     net_data *data_in;
     net_data *data_out;
{
  int i;
  data_out->cmd.cmd_len = data_in->cmd.cmd_len;
  data_out->cmd.cmd_val = (u_char *)strdup(data_in->cmd.cmd_val);
  for (i = 0; i < data_out->cmd.cmd_len ; i++)
    {
      data_out->cmd.cmd_val[i] = toupper( (char) data_out->cmd.cmd_val[i]);
    }
  
  return(0);
}


int main ()
{
  int  ret;
  ret = wait_for_client ( &r_toupper );
  printf("return with  %i and errno %i \n",ret, errno);
  
}

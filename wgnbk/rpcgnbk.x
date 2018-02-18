/* $Id: rpcgnbk.x,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */

/*
** This file gets processed by rpcgen to make machine specific 
** rpc library hooks
**
** gnbk_data is for transfer from client to server
** gnbk_reponse is the canonical rpcgen union.
*/

/* 
   question:
      set by client: a buffer containing the request

   reponse: 
      set by server: a buffer containing the answer

*/

#define RPC_PORT rpc_port

struct gnbk_data {
  string   question <>;
  opaque   reponse <> ;
  };
 
union gnbk_reponse switch ( int errnumber) {
case 0:
  gnbk_data    res_data;
default:
  void;
};

/*
** Please don't change this !!!
*/

program RPC_GNBK {
  version RPC_GNBK_VERS   {
    gnbk_reponse  GNBK_SERVER(gnbk_data) = 1;
  } = 1;
} = RPC_PORT ;

/* const RPC_ANYSOCK = rpc_socket; should be defined already */

/********* end of file **********/



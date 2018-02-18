#ifndef TCP_CONNECT_DEF
#define TCP_CONNECT_DEF

#if !defined(LINUX) && !defined(OPTERON) && !defined (CYGWIN) && !defined (ALPHA)
 #include <sys/filio.h>
#else
 #include <sys/ioctl.h>
#endif
#if defined (CYGWIN) &&  !defined (FIONREAD)
#include <sys/termios.h>
#if  ! defined (FIONREAD) && defined (TIOCINQ)
#define FIONREAD TIOCINQ
#endif
#endif

typedef struct outgoing_transaction_rec *OTR ;
typedef struct incoming_transaction_rec *ITR ;

int connect_socket (const char *host, int port, int no_delay) ;
ITR incoming_transaction_create(void) ;
void incoming_transaction_free(ITR t) ;
void outgoing_transaction_free(OTR t) ;
int incoming_transaction(ITR t, int fd, int limit_one,
			 void (*process_incoming)(char *s, int n, void *callback_cookie), void *callback_cookie ) ;


OTR outgoing_transaction_create (int fd) ;
int outgoing_transaction_send (OTR t, int type, const char *message, int message_len ) ;
void daemon_main(int daemon_port, int daemon_fork, int (*daemon_service)(int) ) ;

#endif

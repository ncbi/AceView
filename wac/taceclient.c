#include <signal.h>

#include <wac/ac.h>
#include <array.h>

void usage()
{
fprintf(stderr,"usage:\n");
fprintf(stderr,"accmd database [ -prompt xx ] [ -timeout seconds ] [ -ace_in | -ace_out ]\n");
fprintf(stderr,"         database is one of\n") ;
fprintf(stderr,"              r:host:port example r:vesta:2000100, for oldstyle rpc server\n") ;
fprintf(stderr,"              a:host:port example a:vesta:12345, Mark's simple tcp tgifacemblyserver\n") ;
fprintf(stderr,"              s:hots:port:user:passwd example s:kaa:12345:anonymous, Sanger tcp taceserver\n") ;
fprintf(stderr,"	-prompt xx	set the prompt\n");
fprintf(stderr,"	-ace_in		read standard input as a .ace file to parse into \n");
fprintf(stderr,"			the database \n");
fprintf(stderr,"	-ace_out	suppress the prompt and command output.  the special \n");
fprintf(stderr,"			command 'write' displays the result of the last \n");
fprintf(stderr,"			command issued.  For example, you might: \n");
fprintf(stderr,"				accmd db -ace_out << EOF \n");
fprintf(stderr,"				find object \n");
fprintf(stderr,"				query tag \n");
fprintf(stderr,"				show -a \n");
fprintf(stderr,"				write \n");
fprintf(stderr,"				EOF \n");
fprintf(stderr,"			to extract a .ace file from the database \n");
fprintf(stderr,"	-timeout sec	exit if open takes longer than <sec> seconds. \n");
fprintf(stderr,"exit codes are:\n");
fprintf(stderr,"0	ok\n");
fprintf(stderr,"1	database open failed immediately \n");
fprintf(stderr,"2	database open failed after timeout \n");
fprintf(stderr,"3	server disconnect during command \n");
fprintf(stderr,"4	unrecognized parameters \n");
fprintf(stderr,"5	shutdown failed\n");
exit(4);
}

static void timed_out(int x)
{
write(2,"timeout\n",8);
exit(2);
}

int main(int argc,char **argv)
{
	AC_DB db;
	const char *err ;
	char *prompt, *dbname, *s, *s1, *s2, *s3 ;
	unsigned char *response = 0 ;
	int ace_out = 0, ace_in = 0, ace_form = 0 ;
	char b[2000];
	int len ;
	int timeout = 0;

	argv++;

	if (!argv[0])
		usage();

	dbname = argv[0];

	argv++;
	
	prompt=halloc(strlen(dbname)+20, 0);
	sprintf(prompt,"acedb@%s>",dbname);
	while (argv[0])
		{
		  if (strcasecmp(argv[0],"-help") == 0 || strcasecmp(argv[0],"--help") == 0 || strcasecmp(argv[0],"-h") == 0)
		    usage () ;
		if (strcasecmp(argv[0],"-prompt") == 0 || strcasecmp(argv[0],"-p") == 0)
			{
			argv++;
			if (!argv[0])
				usage();
			prompt=argv[0];
			}
		else if (strcasecmp(argv[0],"-ace_out") == 0 || strcasecmp(argv[0],"-aceout") == 0)
			{
			ace_out = 1;
			prompt="";
			}
		else if (strcasecmp(argv[0],"-ace_in") == 0 || strcasecmp(argv[0],"-acein") == 0)
			ace_in = 1;
		else if (strcasecmp(argv[0],"-form") == 0)
			ace_form = 1;
		else if (strcasecmp(argv[0],"-timeout") == 0 || strcasecmp(argv[0],"-time_out") == 0)
			{
			argv++;
			if (!argv[0])
				usage();
			timeout = atoi(argv[0]);
			}
		else
			usage();
		argv++;
		}

	if (timeout)
		{
		signal(SIGALRM, timed_out);
		alarm(timeout);
		}
	db = ac_open_db ( dbname, &err );
	if (timeout)
		alarm(0);
	if (!db)
		{
		fprintf(stderr,"database open error: %s\n",err);
		exit(1);
		}

	if (ace_in)
		{
		/*
		* ace_in is a special case - all of stdin will be an ace file
		*/
		Stack st;
		AC_KEYSET ks;
		const char *err;
		ks = 0;

		st = stackCreate(10000);

		b[sizeof(b)-1]=0;
		while (fgets(b,sizeof(b)-1,stdin))
			catText(st, b);

		ac_parse( db, stackText(st,0), &err, &ks, NULL_HANDLE);
		/*
		  if (!ok)
		  fprintf(stderr,"parse error %s \n", err ? err : "");
		  else
		*/
		fprintf(stderr,"%d objects parsed \n", ac_keyset_count(ks));
		if (err && *err)
		  fprintf(stderr,"%s",err);
		exit(0);
		}

	if (ace_form)
	  {
	    while (fgets(b,sizeof(b),stdin))
	      {
		s = b ;
		s1 = strchr(b,'\n');
		if (s1) 
		  *s1 = 0;
		while (s && *s)
		  {
		    s1 = strstr(s,"</") ;
		    if (s1) 
		      *s1 = 0;
		    if (s[0])
		    printf ("%s", s) ; 
		    if (s1)
		      {
			s = s1 + 2 ;
			s1 = strstr(s,"/>") ;
			if (! s1)
			  { /* we did not recognize a query, we print as is */
			    printf ("%s", s-2) ;
			    s = 0 ;
			  }
			else
			  {
			    *s1 = 0 ;
			    /* perform all the ; separated queries at once */
			    while (s)
			      {
				s3 = strstr(s,";;") ;
				if (s3)
				  *s3 = 0 ;
				if (*s)
				  {
				    if (response)
				      messfree(response);
				    response = ac_command(db, s, &len, 0);
				    if (!response)
				      {
					printf("server disconnect\n");
					exit(3);
				      }
				  }
				s = s3 ? s3 + 2 : 0 ;
			      }
			    s2 = strstr ((char *)response, "// Found ") ;
			    s3 = s2 ? strstr (s2, " objects") : 0 ;
			    if (s3)
			      {
				*s3 = 0 ;
				printf("%s", s2+9) ;
			      }
			    s = s1 + 2 ;
			  }
		      }
		    else
		      s = 0 ;
		  }
		printf ("\n") ;
	      }
	    exit (0) ;
	  }

	response = 0;

	for (;;)
		{
		char b[2000];
		int len = 0 ;

		printf("%s",prompt);
		fflush(stdout);
		if (! fgets(b,sizeof(b),stdin))
			break;
		s = strchr(b,'\n');
		if (s) 
			*s = 0;
		if (strcmp(b,"quit") == 0)
			{
			/* not really necessary, but people want to type 'quit' */
			break;
			}
		else if (strcmp(b,"write") == 0)
			{
			  if (len) write(1,response,len);
			}
		else if (strncmp(b,"shutdown",8) == 0)
			{
			response = ac_command(db, b, &len, 0);
			if (response && *response && strstr((char*)response,"Sorry") == 0)
				{
				/*
				* the server would not let us shut it down.
				*/
				write(1,response,len);
				exit(5);
				}
			else
			  {
			    write(1,"shutdown\n",10) ;
			    exit(0);
			  }
			}
		else
			{
			if (response)
				messfree(response);
			response = ac_command(db, b, &len, 0);
			if (!response)
				{
				printf("server disconnect\n");
				exit(3);
				}
			if (! strncasecmp ((char *)response, "shutdown", 8))
			  {
			    write(1,response,len);
			    exit (0) ;
			  }
			if (! ace_out)
				write(1,response,len);
			}
		}
	exit(0);
}

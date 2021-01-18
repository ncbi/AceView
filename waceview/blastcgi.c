#include <regular.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/socket.h>
#include <ctype.h>
#include <waceview/cgi.h>
#include <wtcp/tcp_connect.h> 

/*
* notes:
* - This is a cgi.  It leaks memory.
* - There are no nice error messages for users who don't use the form to submit.
*/

static char *search_type="p";
static char *evalue="";
static char dbchars[32]="+";
static char *fromdb;


/*
* from the sequence entered, try to guess what type of a
* search it is.  The choices are:
*/

static char *guess_search_type(char *s)
{
  static char *dna_chars = "tcagn";
  static char *pep_chars = "acdefghiklmnpqrstvwyx*zbu";
  int flags;
  char ch;
  
  flags = 0;
  
 newline:
  /*
   * if it is a fasta name line, skip over the whole line
   */
  if (*s == '>')
    {
      while (*s && *s != '\n')
	s++;
      goto newline;
    }
  
 next_char:
  /*
   * examine a character
   */
  ch = *s;
  if (isupper((int)ch)) 
    ch = tolower(ch);
  if (ch == '\0')
    goto end_str;
  else if (strchr(dna_chars, ch))
    flags |= 1;
  else if (strchr(pep_chars, ch))
    flags |= 2;
  else if (ch == '\n')
    {
      s++;
      goto newline;
    }
  else if (isspace((int)ch))
    ;
  else
    return "?";
  s++;
  goto next_char;
  
 end_str:
  /*
   * If it contains bogus characters, it is an error.
   * If it contains peptide characters (but no bogus characters), 
   *	it is a peptide even if it contains DNA.
   * If it contains DNA characters (but not bogus and not peptide),
   *	it is DNA
   * else we didn't see any valid characters and it is an error
   */
  if (flags & 4)
    return "?";	/* bogus characters */
  if (flags & 2)
    return "p";	/* contains peptide characters */
  if (flags & 1)
    return "n";	/* contains dna characters */
  
  return "?";		/* no valid characters */
}


/*
 * output the link back to the main AceView page
 */
static void link_aceview()
{
  printf("<a href='index.html?%s",fromdb);
  printf("' target='_top'><IMG src='images/aceview.gif' width=140 ");
  printf(" border=2></a><br>");
}

/*
 * output the link to a new blast query
 */
static void link_blast()
{
  printf("<a href='blast.html?");
  printf("%s",dbchars);
  printf(",%s",evalue);
  printf(",%s",fromdb);
  printf("'>New BLAST search</a><p>\n");
}

/*
 * Create the link back into av.cgi, for a particular mrna
 */
static void cook_link(char *s)
{
char *cp, *dbl, *url_e ;
  switch (*s++)
    {
    case 'h':
      dbl = "human";
      break;
    case 'o':
      dbl = "human";
      break;
    case 'a':
      dbl = "ara";
      break;
    case 'w':
      dbl="worm";
      break;
    case 'm':
      dbl="mouse";
      break;
    case 'r':
      dbl="rat";
      break;
    default:
      dbl="xxx";
      break;
    }
  /* the name looks like 'hGENE|TAG|mRNA|stuff' we wish to call back on mRNA */
  cp = strchr (s, '|') ;
  url_e = url_encode(s);
  if (cp) 
    {
      *cp = 0 ;
      url_e = url_encode(s);
      s = cp+1 ; cp = strchr (s, '|') ;
      if (cp)
	{
	  s = cp+1 ; cp = strchr (s, '|') ;
	  if (cp) *cp = 0 ;
	}
    }
  url_e = url_encode(s);
  printf("%-5s ", dbl);
  printf("<a href='av.cgi?db=%s&c=mrna&q=%s'>",dbl, url_e);
  printf("%s</a>", url_e);
}

void errorfunction(int error)
{
  printf("<HTML><TITLE>BLAST Search Results</TITLE><BODY "
	 "BGCOLOR='#FFFFFF' LINK='#0000FF' VLINK='#660099' ALINK='#660099'>"
	 "<A name=top href='https://www.ncbi.nlm.nih.gov' target='_top'><IMG src='images/ncbi_logo.gif' width=130 "
	 "border=2></A><br>\n");
  link_aceview();
  printf(
	 "<h2>Error</h2><ul>\n");
  if (error & 2)
    printf("<li>You must select at least one species to search\n");
  if (error & 1)
    printf("<li>The Evalue must be a number\n");	
  if (error & 4)
	printf("<li>You must enter a sequence to search for\n");
  if (error & 8)
    printf("<li>Unexpected database problem - try again later\n");
  if (error & 16)
    printf("<li>You entered invalid characters in your sequence\n");
  printf("</ul><p>\n");
  link_blast();
  exit(0);
}

int main()
{
  int fd;
  char buf[1000];
  char blastd_args[100];
  char *s, *sequence, *t;
  double d;
  int anydb;
  FILE *fp;
  int error;
  char blast_host[100];
  int blast_port=0;
  
  setbuf(stdout,NULL);
  
  blast_host[0] = 0;
  fp = fopen("blast_host","r");
  if (fp)
    {
      if (fgets(blast_host,sizeof(blast_host),fp))
	{
	  s = strchr(blast_host,'\n');
	  *s=0;
	}
      fclose(fp);
    }
  if (blast_host[0])
    {
      char *s;
      s = strchr(blast_host,' ');
      if (s)
	{
	  *s++ = 0;
	  blast_port = atoi(s);
	  if (! blast_port)
	    blast_port = 3001;
	}
    }
  else
    {
      strcpy(blast_host,"ace01");
      blast_port=3001;
    }
  
  printf("content-type: text/html\n\n");
  
  putenv("REQUEST_METHOD=POST") ;
  cgi_init("blast.log");

  error = 0;
  
  /*
   * get the evalue
   */
  evalue = cgi_clean_arg("evalue");
  if (evalue)
    {
      if (sscanf(evalue,"%lf",&d) != 1)
	error |= 1;
      sprintf(blastd_args,"%g",d);
      evalue = strdup(blastd_args);
    }
  else if (0)
    exit(0) ;

  strcat(blastd_args," ");
  
  /*
   * find which type of search - p or n
   */
  fromdb = cgi_clean_arg("fromdb");
  if (! fromdb)
    fromdb = "human";
  
  
  /*
   * the sequence
   */
  sequence = cgi_clean_arg("sequence");
  if (0)
    {
      printf("evalue=%s\n", evalue ? evalue : "NULL") ;
      printf("sequence==%s\n", sequence ? sequence : "NULL") ;
    }
  if (! sequence)
    exit(0);
  
  if (! *sequence)
    error |= 4;
  
  /*
   * guess search type from sequence
   */
  search_type = guess_search_type(sequence);
  
  if (*search_type == '?')
    error |= 16;
  
  strcat(blastd_args,search_type);
  strcat(blastd_args," ");
  
  /*
   * find which databases they want to search
   */
  anydb=0;
  if (cgi_clean_arg("db_worm"))
    {
      anydb=1;
      strcat(blastd_args,"n.worm_");
      strcat(dbchars,"w");
    }
  if (cgi_clean_arg("db_human"))
    {
      anydb=1;
      strcat(blastd_args,"n.human_");
      strcat(dbchars,"h");
    }
  if (cgi_clean_arg("db_mouse"))
    {
      anydb=1;
      strcat(blastd_args,"n.mouse_");
      strcat(dbchars,"m");
    }
  if (cgi_clean_arg("db_rat"))
    {
      anydb=1;
      strcat(blastd_args,"n.rat_");
      strcat(dbchars,"r");
    }
  if (cgi_clean_arg("db_omim"))
    {
      anydb=1;
      strcat(blastd_args,"n.omim_");
      strcat(dbchars,"h");
    }
  if (cgi_clean_arg("db_ara"))
    {
      anydb=1;
      strcat(blastd_args,"n.ara_");
      strcat(dbchars,"a");
    }
  
  if (! anydb)
    error |= 2;
  
  if (error)
    errorfunction(error);
  if (0)  printf("<br>hello error=%d<br>\n", error) ;  
  /*
   *
   */
  fd = connect_socket(blast_host,blast_port,0) ;
  
  if (fd < 0)
    {
      printf("<h1>Blast search server (%s:%d) is down</h1>\n", blast_host, blast_port);
      printf("This is unexpected. Please complain to mieg@ncbi.nlm.nih.gov, merci.\n");
      exit(0);
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
  
#define PREFIX(b,s) ( strncmp(b,s,sizeof(s)-1) == 0 )
  
  /*
   * up to and including "^<BODY ", output lines as-is
   */
  while (fgets(buf,sizeof(buf),fp) && ! PREFIX(buf,"<BODY "))
    {
      if (strcmp(buf,"dberror\n") == 0)
	errorfunction(8);
      printf("%s",buf);
    }
  
  printf("%s",buf);
  
  /*
   * insert our page header
   */
  printf( "<A name=top href='https://www.ncbi.nlm.nih.gov' target='_top'><IMG src='images/ncbi_logo.gif' width=130 border=2></A><br>"
	  );
  link_aceview();
  link_blast();
  
 another:
  /*
   * output up to "Sequences producing significant alignments"
   */
  while ( fgets(buf,sizeof(buf),fp) && ! PREFIX(buf,"Sequences producing significant alignments") )
    printf("%s",buf);
  
  printf("%s",buf);
  
  /*
   * skip blank line
   */
  if (!fgets(buf,sizeof(buf),fp))
    printf("%s",buf);
  
  /*
   * Until we get to the </PRE>, tag the first word of each line
   */
  while (fgets(buf,sizeof(buf),fp) && ! PREFIX(buf,"</PRE>"))
    {
      s = strchr(buf,' ');
      if (s)
	{
	  *s++ = 0;
	  cook_link(buf);
	  printf(" %s",s);
	}
      else
	cook_link(buf);
    }
  
  printf("\n");
  
  /*
   * at each line "^><a name = "
   *        tag from last '>' to end of line
   */
  
  while (fgets(buf,sizeof(buf),fp))
    {
      if (PREFIX(buf,"<PRE><b>BLASTN 2.2.5 [Nov-16-2002]</b>"))
	goto another;
      if (PREFIX(buf,"<b>TBLASTN"))
	goto another;
      if (PREFIX(buf,"<b>BLASTN"))
	goto another;
      if (PREFIX(buf,"><a name = "))
	{
	  s = strrchr(buf,'>');
	  if (s)
	    {
	      *s++=0;
	      printf("%s>",buf);
	      t=strchr(s,' ');
	      if (t)
		{
		  *t++ = 0;
		  cook_link(s);
		  printf(" %s",t);
		}
	      else
		{
		  t = strchr(s,'\n');
		  if (t) *t=0;
		  cook_link(s);
		}
	    }
	  else
	    printf("->s null\n");
	}
      else
	printf("%s",buf);
    }
  return 0;
}


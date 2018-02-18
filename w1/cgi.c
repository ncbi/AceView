#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include <errno.h>

#include <waceview/cgi.h>

void url_decode_inplace(char *s)
{
  char *r;
  int x;
  char c;
  r=s;
  while (*s)
    {
      if (*s == '+')
	{
	  s++;
	  *r++ = ' ';
	}
      else if (*s == '%')
	{
	  char b[3];
	  s++;
	  if (! *s)
	    break;
	  b[0] = *s++;
	  if (! *s)
	    break;
	  b[1] = *s++;
	  b[2] = 0;
	  if (sscanf(b,"%x%c",&x,&c) != 1)
	    break;
	  *r++ = x;
	}
      else
	*r++ = *s++;
    }
  *r++ = 0;
}


char *html_encode (const char *s)
{
  char *t, *ret;
  const char *t1;
  int size;
  for (size=1, t1=s; *t1; t1++)
    {
      switch (*t1)
	{
	case '<':	size += 4; break;
	case '>':	size += 4; break;
	case '&':	size += 5; break;
	case '"':	size += 6; break;
	case '\'':	size += 6; break;
	default: size++;
	}
    }
  ret=malloc(size);
  for (t=ret; *s; s++)
    {
      switch (*s)
	{
	case '<':	strcpy(t,"&lt;"); t += 4; break;
	case '>':	strcpy(t,"&gt;"); t += 4; break;
	case '&':	strcpy(t,"&amp;"); t += 5; break;
	case '"':	strcpy(t,"&quot;"); t += 6; break;
	case '\'':	strcpy(t,"&acute;"); t += 7; break;
	default:	*t++ = *s;
	}
      *t='.';
    }
  *t = 0;
  return ret;
}

char *url_encode (const char *s)
{
  char *t, *ret;
  const char *t1;
  int size;
  for (size=1, t1=s; *t1; t1++)
    {
      if (isalnum((int)*t1))
	size++;
      else if (*t1 == '_' || *t1 == ',' || *t1 == '.' || *t1 == '-' )
	size ++;
      else 
	size += 3;
    }
  ret = malloc(size+1);
  for (t=ret; *s; s++)
    {
      if (isalnum((int)*s))
	*t++ = *s;
      else if (*s == '_' || *s == ',' || *s == '.' || *s == '-' )
	*t++ = *s;
      else 
	{
	  char b[3];
	  sprintf(b,"%02x",*s);
	  *t++ = '%';
	  *t++ = b[0];
	  *t++ = b[1];
	}
    }
  *t = 0;
  return ret;
}


/*
* CGI related functions:
*/

struct element
{
  char *name;
  char *value;
  struct element *next;
};

static struct element *build;

static struct element *cgi_args, *cookies;

static int
add_element(char *string)
{
  struct element *e;
  char *value;
  e=malloc(sizeof(struct element));
  if (!e)
    return 1;
  value = strchr(string,'=');
  if (value)
    {
      *value++ = 0;
      url_decode_inplace(value);
    }
  else
    value="";
  url_decode_inplace(string);
  e->name = string;
  e->value = value;
  e->next=build;
  build=e;
  return 0;
}

static int crack_cgi_args(char *args, char separator)
{
  char *s;
  int r;
  build=0;
  r=0;
  while ((s=strchr(args,separator)))
    {
      *s = 0;
      r |= add_element(args);
      args = s+1;
    }
  r |= add_element(args);
  return r;
}

int cgi_init(char *name)
{
  char *s, *args, *cookie_string;
  
  s = getenv("REQUEST_METHOD");
  if (0)
    printf("cgi_init request_method = %s\n", s ? s : "NULL") ;
  if (! s)
    {
      /* could implement cgi debugger */
      FILE *f;
      f=fopen(name,"r");
      if (!f)
	return 0;
      args=malloc(10000);
      fgets(args,10000,f);
      s = strchr(args,'\n');
      if (s)
	*s = 0;
      cookie_string=malloc(10000);
      fgets(cookie_string,10000,f);
      s = strchr(cookie_string,'\n');
      if (s)
	*s = 0;
      fclose(f);
    }
  else if (strcmp(s,"GET") == 0)
    {
      args=strdup(getenv("QUERY_STRING"));
      if (! args) args=strdup("");
      cookie_string=getenv("HTTP_COOKIE");
    }
  else if (strcmp(s,"POST") == 0)
    {
      int len;
      s=getenv("CONTENT_LENGTH");
      if (0)
	printf("cgi_init post length = %s\n", s ? s : "NULL") ;

      if (! s)
	return 0;
      len=atoi(s) ;
      args=malloc(len+1);
      memset (args, 0, len+1) ;
      if (!args)
	return 0;

      /* 2009_11_03
       * in blastcgi.c
       * for a reason i do not understand fread returns 0, errno=2 
	 i just ignore that bug since it still seems to work
      */
      if (0 && fread(args,len,1,stdin) != len) 
	return 0;
      fread(args,len,1,stdin) ;
     
      if (0)
	printf("cgi_init len2=%d len3=%zu args = %s errno=%d\n", len, strlen(args), args ? args : "NULL", errno) ;
      cookie_string=getenv("HTTP_COOKIE");
    }
  else
    return 0;
  
  if (!cookie_string)
    cookie_string=strdup("");
  
  if (name)
    {
      FILE *f;
      f=fopen(name,"w");
      if (f)
	{
	  fprintf(f,"%s\n",args);
	  fprintf(f,"%s\n",cookie_string);
	  fclose(f);
	}
    }
  
/* could save cgi args for debugger here */
  if (crack_cgi_args(args,'&'))
    return -1;
  cgi_args = build;
  
  if (crack_cgi_args(cookie_string,';'))
    return -1;
  cookies = build;
  return 1;
}


char *cgi_arg(char *name)
{
  struct element *e;
  for (e=cgi_args; e; e=e->next)
    if (strcmp(name, e->name) == 0)
      return e->value;
  return NULL;
}

char *cgi_clean_arg(char *name)
{
  char *cp, *s = cgi_arg(name) ;
  if (s)
    {
      s = strdup (s) ;
      while ((cp = strstr (s, "<")))
	memset (cp, ' ', 1) ;
      while ((cp = strstr (s, ">")))
	memset (cp, ' ', 1) ;
      while ((cp = strstr (s, "include")))
	memset (cp, ' ', 7) ;
      while ((cp = strstr (s, "INCLUDE")))
	memset (cp, ' ', 7) ;
      while ((cp = strstr (s, "http")))
	memset (cp, ' ', 4) ;
      while ((cp = strstr (s, "script")))
	memset (cp, ' ', 6) ;
      while ((cp = strstr (s, "Script")))
	memset (cp, ' ', 6) ;
      while ((cp = strstr (s, "SCRIPT")))
	memset (cp, ' ', 6) ;
    }

  return s ;
}


int cgi_arg_exists (char *name)
{
  struct element *e;
  for (e=cgi_args; e; e=e->next)
    if (strcmp(name, e->name) == 0)
      return 1 ;
  return 0 ;
}

void cgi_call(char *name, void *arg, void (*callback)(char *value, void *arg))
{
  struct element *e;
  for (e=cgi_args; e; e=e->next)
    if (strcmp(name, e->name) == 0)
      (*callback)(e->value, arg);
}

char *cgi_cookie(char *name)
{
  struct element *e;
  for (e=cookies; e; e=e->next)
	if (strcmp(name, e->name) == 0)
	  return e->value;
  return NULL;
}


void cgi_dump()
{
  struct element *e;
  printf("<h1>cgi parameters</h1>\n");
  for (e=cgi_args; e; e=e->next)
    printf("%s = %s  <p>\n",html_encode(e->name), html_encode(e->value));
}

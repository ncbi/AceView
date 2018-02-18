/*
 * User notes:
 * 
 * see AceView/help/pfam_help.html
 * 
 * 
 * Theory of operation:
 * 
 * If this CGI is run without parameters or with the "clear" parameter, it
 * produces the first page form, which has a blank query.
 * 
 * When you submit a query, it creates a keyset that matches your query and
 * stores it in a "context".  (The context is really a keyset file on the
 * server.)  When you refine the query, it makes various further queries
 * against that keyset and creates a new context.
 * 
 * You MUST create a new context for each query result, because you cannot know
 * if the old context is still in use -- the user may have type ^N at IE or
 * used "open in new window".
 * 
 * When you refine a query, it loads the keyset from the context, performs some
 * number of queries on it, then saves the result in a new context.
 * 
 * Results are displayed by using acCommand to create a table using tablemaker,
 * and parsing the resulting output lines.
 * 
 */

/*
 * Note: This program leaks some of the memory it allocates.  A CGI will only
 * run for a short time, so it does not much matter if it leaks a few KB.
 */

#include <stdio.h>
#include <unistd.h>
#include <ctype.h>

#include <waceview/cgi.h>
#include <wh/regular.h>
#include <waceview/avlib.h>
#include <wac/ac.h>

struct av_db_desc * my_db_desc;

/*
 * query context is currently stored as a temp file on the database server.
 * This is the directory where we keep the context.  You need to set
 * something up to delete files from this directory pretty regularly
 */

#define CTXDIR "/home/mieg/ee/SERVER/QUERY_TMP"

static char           *db;		/* the database we are querying */

static AC_DB          db_descriptor;	/* descriptor for database, using new acec */

static char           *term;

static AC_KEYSET db_context_keyset = 0;

static char *context_name;

static void create_context()
{
  char            host[100], day[50];
  time_t          t;
  struct tm      *tp;
  host[0] = '?';
  host[1] = 0;
  gethostname(host, sizeof(host));
  t = time(0);
  tp = gmtime(&t);
  strftime(day, 50, "%j%H", tp);
  context_name = strdup(messprintf("ctx-%s-%s-%d", day, host, getpid()));
}


static void load_context()
{
  char           *q;
  for (q = context_name; *q; q++)
    if (!isalnum((int) *q) && (*q != '.') && (*q != '_') && (*q != '-'))
      exit(1);
  if (strchr(context_name, '/'))
    exit(1);
  q = strdup(messprintf( CTXDIR "/%s", context_name));
  db_context_keyset = ac_read_keyset (db_descriptor, q, NULL);
  if (! db_context_keyset) {
    printf("<h1>Query results expired </h1>\n");
    printf("<p>Query results are only stored for about an hour after you make ");
    printf("the query.  Click <a href='pfam.cgi'>here</a> to start a new query\n");
    exit(0);
  }
}

static void save_context()
{
  char *q;
  create_context();
  q = strdup(messprintf( CTXDIR "/%s", context_name));
  ac_keyset_write ( db_context_keyset, q );
}

static void db_select(char *db)
{
  /*
   * hack for now
   */
  char           *hs = "", *ms = "", *ws = "", *as = "";
  printf("<select name='db'>\n");
  if (strcmp(db, "human") == 0)
    hs = "selected";
  if (strcmp(db, "mouse") == 0)
    ms = "selected";
  if (strcmp(db, "worm") == 0)
    ws = "selected";
  if (strcmp(db, "ara") == 0)
    as = "selected";
  printf("<option value='human' %s>Human</option>\n", hs); 
  printf("<option value='mouse' %s>Mouse</option>\n", ms); 
  printf("<option value='worm' %s>C elegans</option>\n", ws);
  printf("<option value='ara' %s>A thaliana</option>\n", as);
  printf("</select>\n");
}

static char *db_name()
{
  /*
   * hack for now
   */
  if (strcmp(db, "human") == 0)
    return "Human";
  if (strcmp(db, "mouse") == 0)
    return "Mouse";
  if (strcmp(db, "worm") == 0)
    return "C elegans";
  if (strcmp(db, "ara") == 0)
    return "A thaliana";
  return "?";
}

static void page_top(char *query, char *db)
{
  printf(
	 "<html>\n"
	 "<head><title>PFAM search</title></head>\n"
	 "<body bgcolor='white'>\n"
	 "<font face='verdana'>\n"
	 "<A name=top href='/' target='_top'>"
	 "<IMG src='/corehtml/left.GIF' width=130 border=0>"
	 "</A>"
	 "<br>\n"
	 );
  printf(
	 "<a href='index.html?%s' target='_top'><IMG src='images/aceview.gif' width=140 border=0></a>"
	 "<br>\n",
	 db
	 );
}

static void page_first_time()
{
  printf(
	 "\n"
	 "<FORM method=get action='pfam.cgi'>\n"
	 );
  
  printf(
	 "\n"
	 "<h2>Pfam Protein Family search\n"
	 "</h2>\n"
	 );
  printf("<a href='HelpQuery.html#pfam_search'>Help</a>\n") ; /* ,my_db_desc->db_version; */
  printf(
	 "<p>Search for word: "
	 );
  
  printf(
	 "<INPUT type=text name=term size=40>\n");
  
  printf(
	 "in the name, definition, or description of all Pfam motifs and related Gene Ontology terms. "
	 "<p>For example, search PF00988, PH, dehydrogenase, sugar metabolism, or mitotic spindle."
	 "<p>in:\n"
	 );
  db_select(db);
  printf(
	 "<br>\n"
	 "<INPUT type=submit name=search value='Search'>\n"
	 );
  
  printf("</form>\n");
  printf("</body></html>\n");
}

static int             f_row = 1;	/* first row to display */
static int             n_row = 50;	/* how many rows to display */

static void table_top(char *name, char *num, char *def, char *prod, char *genes, int n_rows, int tot_prod, int tot_gene)
{
  
  printf("<a href='pfam.cgi?db=%s'>New search</a>\n", db);
  printf("<br><a href='HelpQuery.html#pfam_search'>Help</a>") ; /* ,my_db_desc->db_version); */
  
  printf(
	 "<h2>Search Results Found %d Pfam motifs in %s</h2>\n", n_rows, db_name());
  
  printf(
	 "<form method=get action='pfam.cgi'>\n"
	 );
  
  printf("<p>");
  printf("<input type=hidden name=context value='%s'>\n", context_name);
  printf("<input type=hidden name=db value='%s'>\n", db);
  
  printf("<input type=submit name=search value='Refine results'><p>"
	 "Query: <INPUT type=text name=term size=40><br>\n"
	 );
  
  /*
   * Danielle wants to suppress this input box if there are fewer than
   * 10 rows of data.  You also need it if the number of rows is
   * greater than the number we are displaying, even if the number we
   * display is less than 10.  You can't use just n_rows > n_row
   * because you want to allow them to reduce the number of rows if
   * they set it to 300 or something.
   */
  if ((n_rows > 10) || (n_rows > n_row)) {
    printf("<p>Show <input type=text size=5 name=n_row value=%d>\n", n_row);
    printf(" rows starting at row <input type=text size=6 name=f_row value=%d ><br>\n", f_row);
    printf("<input type=submit name=prev value='Previous Page'> \n");
    printf("<input type=submit name=next value='Next Page'> \n");
  }
  printf("<p>");
  
  printf(
	 "<p>\n"
	 "<table width=100%% border>\n"
	 "\n"
	 );
  
  
  printf(
	 "<tr>\n"
	 );
  printf("<th valign='bottom'>remove<br>checked</th>\n");
  
  /*
   * the nbsp prevents some versions of internet explorer from putting
   * a line break between the check box and the text.
   */
  printf("<th align='left'>");
  printf("<INPUT type='radio' name='n_name' value='y' checked>&nbsp;include<BR>");
  printf("<INPUT type='radio' name='n_name' value='n'>&nbsp;exclude<BR>");
  printf("<INPUT type=text name=a_name size=10><br>");
  printf("</th>\n");
  
  printf("<th align='left'>");
  printf("<INPUT type='radio' name='n_num' value='y' checked>&nbsp;include<BR>");
  printf("<INPUT type='radio' name='n_num' value='n'>&nbsp;exclude<BR>");
  printf("<INPUT type=text name=a_num  size=8 ><br>");
  printf("</th>\n");
  
  printf("<th align='left'>");
  printf("<INPUT type='radio' name='n_def' value='y' checked>&nbsp;include<BR>");
  printf("<INPUT type='radio' name='n_def' value='n'>&nbsp;exclude<BR>");
  printf("<INPUT type=text name=a_def  size=50 ><br>");
  printf("</th>\n");
  
  printf("<th valign='bottom'>");
  printf("<INPUT type='radio' name='n_prod' value='g' checked>&nbsp;&gt;<BR>");
  printf("<INPUT type='radio' name='n_prod' value='e' >&nbsp;=<BR>");
  printf("<INPUT type='radio' name='n_prod' value='l' >&nbsp;&lt;<BR>");
  printf("<INPUT type=text name=a_prod size=5 >");
  printf("</th>\n");
  
  printf("<th valign='bottom'>");
  printf("<INPUT type='radio' name='n_gene' value='g' checked>&nbsp;&gt;<BR>");
  printf("<INPUT type='radio' name='n_gene' value='e' >&nbsp;=<BR>");
  printf("<INPUT type='radio' name='n_gene' value='l' >&nbsp;&lt;<BR>");
  printf("<INPUT type=text name=a_gene size=5 >");
  printf("</th>\n");
  
  printf(
	 "</tr>\n"
	 "\n"
	 );
  printf(
	 "<tr>\n"
	 "<th> &nbsp; </th>\n"
	 "<th>%d<br>Motifs</th>\n"
	 "<th>Accession</th>\n"
	 "<th>Definition</th>\n"
	 "<th>%d<br>Products</th>\n"
	 "<th>%d<br>Genes</th>\n"
	 "</tr>\n",
	 n_rows,
	 tot_prod, tot_gene
	 );
}

static void announce_download()
{
  
  printf("<p>You can download this table and related data.<br>");
  printf(
	 "<br>Download the table of results "
	 "<a href=\"pfam.cgi?db=%s&download=tt&context=%s\">as tab delimited file</a>\n", db, context_name);
  
  printf("<br>Download list of product names "
	 "<a href=\"pfam.cgi?db=%s&download=ap&context=%s\">as AceDB file</a>\n", db, context_name);
  
  printf("<br>Download list of gene names "
	 "<a href=\"pfam.cgi?db=%s&download=ag&context=%s\">as AceDB file</a>\n", db, context_name);
  
  printf("<br><a href='pfam_downloads.%s.html'>About Pfam query downloads</a>\n",my_db_desc->db_version);
}

static void clean_query(char *s)
{
  for (; *s; s++) {
    switch (*s) {
    case '"':
    case '\'':
    case '%':
      *s = '?';
      break;
    }
  }
}


static unsigned char *get_table (void)
{
  unsigned char           *r;
  
  load_context();
  r = ac_command(db_descriptor, "table -active -n webpfamg", NULL, NULL);
  return r;
}

static void general_query (char *user_query, int find)
{
  char           *q;
  char           *next_query;
  
  /*
   * this requirement came in after all the code was written:  If the
   * user types two or more words, do NOT look for that string.
   * Instead, look for the first word, then look for the second word
   * independently etc, refinining the search each way.  That is "ice
   * cream" matches "Bailey's Irish Cream tastes good on ice".
   */
 hack_loop:
  
  next_query = strchr(user_query, ' ');
  if (next_query)
    *next_query++ = 0;
  
  /*
   * */
  if (find) {
    if (*user_query)
      {
	char           *q1, *q2, *q3, *q4, *q5;
	q1 = strdup(messprintf(
			       "{ find pfam ; ( IS \"*%s*\" ) || ( accession = \"*%s*\" ) || "
			       "( definition = \"*%s*\" ) || ( pfam_comment = \"*%s*\") || ( interpro_comment = \"*%s*\") }",
			       user_query, user_query, user_query, user_query, user_query));
	q2 = strdup(messprintf(
			       "{ find go_b ; IS \"*%s*\" ; follow go ; follow pfam }",
			       user_query));
	q3 = strdup(messprintf(
			       "{ find go_m ; IS \"*%s*\" ; follow go ; follow pfam }",
			       user_query));
	q4 = strdup(messprintf(
			       "{ find go_c ; IS \"*%s*\" ; follow go ; follow pfam }",
			       user_query));
	
	q5 = strdup(messprintf(
			       "{ find paper ; pfam ; ( title = \"*%s*\" || author = \"*%s*\" ) ; follow pfam } ",
			       user_query, user_query));
	
	q = strdup(messprintf("%s SETOR %s SETOR %s SETOR %s SETOR %s", q1, q2, q3, q4, q5));
      }
    else
      {
	q = strdup("find pfam") ;
      }
    
  } 
  else 
    {
      char           *q1, *q2 ;
      /*
       * in this case, we have a keyset full of pfams
       */
      q1 = strdup(messprintf(
			     " ( IS \"*%s*\" ) || ( accession = \"*%s*\" ) || "
			     " ( definition = \"*%s*\" ) || ( pfam_comment = \"*%s*\") || "
			     " ( interpro_comment = \"*%s*\")  || "
			     " ( COUNT {>GO ; Type:2 == \"*%s*\" } > 0 ) ",
			     user_query, user_query, user_query, user_query, user_query, user_query));
      
      q2 = strdup(messprintf(
			     "follow paper ; pfam ; ( title = \"*%s*\" || author = \"*%s*\" ) ; follow pfam ",
			     user_query, user_query));
      
      q = strdup(messprintf("{ %s } SETOR { %s }", q1, q2));
    }
  
  if (db_context_keyset)
    db_context_keyset = ac_ksquery_keyset( db_context_keyset, q, NULL);
  else
    db_context_keyset = ac_dbquery_keyset( db_descriptor, q, NULL);
  
  if (next_query && *next_query) {
    find = 0;
    user_query = next_query;
    goto hack_loop;
  }
}


char *numeric_filter(char *q)
{
  char *s;
  s = cgi_clean_arg(q);
  if (!s)
    return ">";
  if (*s == 'g')
    return ">";
  if (*s == 'e')
    return "=";
  if (*s == 'l')
    return "<";
  return ">";
}


static char not_char(char *q)
{
  char *s;
  s = cgi_clean_arg(q);
  if (!s)
    return ' ';
  if (*s == 'n')
    return '!';
  return ' ';
}

static void filter_query()
{
  char           *q, *s, *t;
  char            not;
  
  q = "";
  
  s = cgi_clean_arg("a_name");
  if (s && *s) {
    s = strdup(s);
    not = not_char("n_name");
    clean_query(s);
    q = strdup(messprintf("%s ( %c ( IS \"*%s*\" ) ) ; ", q, not, s));
  }
  s = cgi_clean_arg("a_num");
  if (s && *s) {
    s = strdup(s);
    clean_query(s);
    not = not_char("n_num");
    q = strdup(messprintf("%s ( %c ( accession = \"*%s*\" ) ) ; ", q, not, s));
  }
  s = cgi_clean_arg("a_def");
  if (s && *s) {
    s = strdup(s);
    clean_query(s);
    not = not_char("n_def");
    q = strdup(messprintf("%s ( %c ( definition = \"*%s*\" ) ) ; ", q, not, s));
  }
  s = cgi_clean_arg("a_prod");
  if (s && *s) {  
    s = strdup(s);
    clean_query(s);
 
    t = numeric_filter("n_prod");
    if (t)
      q = strdup(messprintf("%s ( COUNT product %s %s ) ; ", q, t, s));
  }
  s = cgi_clean_arg("a_gene");
  if (s && *s) {
    s = strdup(s);
    clean_query(s);
    t = numeric_filter("n_gene");
    if (t)
      q = strdup(messprintf("%s ( COUNT gene %s %s ) ; ", q, t, s));
  }
  if (*q) {
    if (db_context_keyset)
      db_context_keyset = ac_ksquery_keyset(db_context_keyset, q, NULL);
    else
      db_context_keyset = ac_dbquery_keyset(db_descriptor, q, NULL);
  }
}

static void download (char *download)
{
  unsigned char *r ;
  const char *err ;
  char *q ;
  
  if (!context_name)
    return;
  
  db_descriptor = ac_open_db (my_db_desc->db_access , &err) ;
  if (! db_descriptor) 
    {
      printf("<h1>error opening database: %s\n", err);
      return;
    }

  load_context();
  switch (download[0])
    {
    case 'a':
      switch (download[1]) {
      case 'p':
	q = "show -a product";
	break;
      case 'g':
	q = "show -a gene";
	break;
      default:
	return;
      }
      printf("// %s\n", q);
      load_context();
      r = ac_command(db_descriptor, q, NULL, NULL);
      printf("%s\n", r);
      return;
    case 't':
      switch (download[1]) {
      case 't':
	r = get_table();
	printf("%s\n", r);
	return;
      }
    }
}


static char *html_safe(char *s)
{
	static char *save;
	if (save)
		free(save);
	if (! s || ! *s )
		{
		save = 0;
		return "&nbsp;";
		}
	save = html_encode(s);
	return save;
}

static char *url_safe(char *s)
{
	static char *save;
	if (save)
		free(save);
	if (! s || ! *s)
		{
		save = 0;
		return "";
		}
	save = url_encode(s);
	return save;
}

static void display_current()
{
  AC_HANDLE h = ac_new_handle () ;
  int             n_rows;
  int row;
  int tot_prod, tot_gene;
  int prod, gene, end_row;
  AC_TABLE	tb;
  AC_KEYSET	ks_p, ks_g;
  
  tb = ac_tablemaker_table ( db_descriptor, "webpfamg",
			     db_context_keyset, ac_tablemaker_db, NULL, NULL, NULL, 0);
  
  n_rows = tb->rows;
  
  ks_p = ac_ksquery_keyset(db_context_keyset, "follow product",  NULL);
  
  tot_prod = ac_keyset_count(ks_p);
  
  ks_g = ac_ksquery_keyset (ks_p, "follow genebox",  NULL);
  
  tot_gene = ac_keyset_count(ks_g);
  
  table_top("", "", "", "", "", n_rows, tot_prod, tot_gene);
  
  end_row = f_row + n_row;
  if (tb->rows+1 < end_row)
    end_row = tb->rows+1;
  
  for (row = f_row; row < end_row; row++)
    {
      char *target, *cp;
      char *html, *url;
      printf("<tr>\n");
      
      /*
       * row number, checkbox
       * pfam name
       */
      
      url = html = strnew (ac_table_printable(tb, row-1, 0, NULL), h) ;
      html = html_safe(html);
      url = url_safe(url);
      target = strdup(url);
      
      printf("<td>%d <input type=checkbox name=del value=%s></td>",
	     row, url);
      printf("<td>%s</td>", html);
      
      /*
       * accession number
       */
      html = url = strnew (ac_table_printable(tb, row-1, 1, NULL), h) ;
      cp = html + strlen(html) - 1 ;
      while (cp > html && *cp != '.') cp-- ;
      if (cp > html) *cp = 0 ;  /* cut out version number */
      html = html_safe(html);
      url = url_safe(url);
      
      printf("<td><a href=http://www.sanger.ac.uk/cgi-bin/Pfam/getacc?%s>%s</a></td>", url, html);
      
      /*
       * definition
       */
      html = strnew (ac_table_printable(tb, row-1, 2, NULL), h) ;
      html = html_safe(html);
      printf("<td>%s</td>",html);
      
      /*
       * products
       */
      prod = ac_table_int(tb, row-1, 3, 0);
      printf("<td>%d</td>",prod);
      
      /*
       * genes
       */
      gene = ac_table_int(tb, row-1, 4, 0);
      printf("<td><a href=\"av.cgi?db=%s&a=list&c=pfam&q=%s\">%d</a></db>", db, target, gene);
      
      printf("</tr>\n");
    }
  
  printf("</table>\n");
  
  printf("<input type=submit name=search value='Refine results'><p>");
  
  announce_download();
  ac_free (h) ;
}


/*
 * Delete a row from the results.
 * 
 * It might be nice to avoid doing a database transaction for each row we
 * delete, but I think this is fast enough and uncommon enough that it does
 * not really matter.
 */
static void keydel(char *value, void *arg)
{
  char           *s;
  
  s = messprintf("! IS %s", value);
  db_context_keyset = ac_ksquery_keyset(db_context_keyset, s, NULL) ;
}

static int func_search(char *user_query)
{
  const char *err = 0 ;
  int             find;
  
  printf("Content-type: text/html\n\n");
  
  db_descriptor = ac_open_db (my_db_desc->db_access , &err );
  if ( ! db_descriptor )
    {
      printf("<h1>error opening database: %s\n", err);
      return 0;
    }
  if (context_name)
    {
      load_context();
      find = 0;
    } 
  else 
    {
      create_context();
      db_context_keyset = ac_new_keyset( db_descriptor, NULL );
      find = 1;
    }
  
  clean_query(user_query);
  
  general_query(user_query, find);
  
  filter_query();
  
  cgi_call("del", 0, keydel);
  
  save_context();
  
  page_top("", db);
  
  display_current();
  
  return 0;
}


int main()
{
  char           *s;
  int             n;
  setbuf(stdout, NULL);
  n = cgi_init("pfamcgi.log");
  if (n < 0)
    {
      printf("Content-type: text/html\n\n");
      printf("<h1>cgi_init -> %d\n</h1>", n);
      exit(1);
    }
  db = cgi_clean_arg("db");
  if (!db)
    db = "worm";
  
  /*
   * does not really open the database - just find it in the
   * list so we can get a "v8" or "v9" to put in URL for the
   * help file.
   */
  my_db_desc = get_db_desc(db, &s);
  if (!my_db_desc) {
    printf("content-type: text/html\n\n");
    printf("database error - %s\n",s);
    exit(1);
  }
#if 0
  printf("content-type: text/plain\n\n");
  printf("version %s\n\n", my_db_desc->db_version);
  printf("name    %s\n\n", my_db_desc->db_name);
  printf("access  %s\n\n", my_db_desc->db_access);
  printf("spec_c  %s\n\n", my_db_desc->db_species_char);
  printf("s_titl  %s\n\n", my_db_desc->db_species_title);
  printf("db_t    %s\n\n", my_db_desc->db_title);
#endif
  
  
  s = cgi_clean_arg("n_row");
  if (s)
    
    n_row = atoi(s);
  if (n_row < 1)
    n_row = 1;
  
  s = cgi_clean_arg("f_row");
  if (s)
    f_row = atoi(s);
  if (f_row < 1)
    f_row = 1;
  
  term = cgi_clean_arg("term");
  if (term)
    {
      term = strdup(term);
      clean_query(term);
    }
  
  context_name = cgi_clean_arg("context");
  if (context_name)
    {
      context_name = strdup(context_name);
      clean_query(context_name);
    }
  
  if (cgi_clean_arg("clear")) {
    printf("Content-type: text/html\n\n");
    page_top("", db);
    page_first_time();
    return 0;
  }
  s = cgi_clean_arg("search");
  if (s)
    {
      return func_search (term) ;
    }
  
  s = cgi_clean_arg("prev");
  if (s) 
    {
      f_row = f_row - n_row;
      if (f_row < 1)
	f_row = 1;
      func_search (term) ;
      return 0;
    }
  s = cgi_clean_arg("next");
  if (s) 
    {
      f_row = f_row + n_row;
      func_search (term) ;
      return 0;
    }
  s = cgi_clean_arg("download");
  if (s)
    {
      printf("Content-type: text/plain\n\n");
      download (term);
      return 0;
    }
  if (term)
    return func_search(term);
  
  printf("Content-type: text/html\n\n");
  page_top("", db);
  page_first_time();
  return 0;
}

/*
 * User notes:
 * 
 * see AceView/help/zone_help.html
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

static struct av_db_desc * my_db_desc ;


/*
 * query context is currently stored as a temp file on the database server.
 * This is the directory where we keep the context.  You need to set
 * something up to delete files from this directory pretty regularly
 */

#define CTXDIR "/net/vesta/a/mieg/SERVER/QUERY_TMP"


static char *db ;		/* the database we are querying */

static AC_DB          db_descriptor ;	/* descriptor for database, using new acec */

static AC_KEYSET db_context_keyset = 0 ;

static char *context_name ;

/*******************************************************************/
/*******************************************************************/

typedef struct tblStruct { 
  char nam [32] ; /* the name of the table in the aceserver */
  int cols ;      /* the number of columns of the table */
  char *title[30] ; /* the title of the columns */  /*** should be dynamic */
  AC_TABLE tbl ;
} TBL ;

/**************/

static TBL *tblCreate (char *nam)
{
  TBL *tbl = messalloc (sizeof (struct tblStruct)) ;

  if (1) 
    {
      strcpy (tbl->nam, nam) ;
      tbl->cols = 4 ;

      tbl->title[0] = "Gene" ;
      tbl->title[1] = "" ;
      tbl->title[2] = "" ;
      tbl->title[3] = "" ;
      tbl->title[4] = "" ;
      tbl->title[5] = "" ;
      tbl->title[6] = "" ;
      tbl->title[7] = "" ;
      tbl->title[8] = "" ;
      tbl->title[9] = "" ;
      tbl->title[10] = "" ;
    }
  return tbl ;
} /*  */

/*******************************************************************/
/*******************************************************************/

static void create_context (void)
{
  char            host[100], day[50] ;
  time_t          t ;
  struct tm      *tp ;
  host[0] = '?' ;
  host[1] = 0 ;
  gethostname (host, sizeof (host)) ;
  t = time (0) ;
  tp = gmtime (&t) ;
  strftime (day, 50, "%j%H", tp) ;
  context_name = strdup (messprintf ("ctz-%s-%s-%d", day, host, getpid ())) ;
} /*  */

/*******************************************************************/

static void load_context (void)
{
  char *q ;
  for (q = context_name ; *q ; q++)
    if (!isalnum ( (int) *q) && (*q != '.') && (*q != '_') && (*q != '-'))
      exit (1) ;
  if (strchr (context_name, '/'))
    exit (1) ;
  q = strdup (messprintf ( CTXDIR "/%s", context_name)) ;
  db_context_keyset = ac_read_keyset (db_descriptor, q, NULL) ;
  if (! db_context_keyset) 
    {
      printf ("<h1>Query results expired </h1>\n") ;
      printf ("<p>Query results are only stored for about an hour after you make ") ;
      printf ("the query.  Click <a href='zone.cgi'>here</a> to start a new query\n") ;
      exit (0) ;
    }
} /*  */

/*******************************************************************/

static void save_context (void)
{
  char *q ;
  create_context () ;
  q = strdup (messprintf ( CTXDIR "/%s", context_name)) ;
  ac_keyset_write ( db_context_keyset, q ) ;
} /*  */

/*******************************************************************/
/*******************************************************************/

static void db_select (char *db)
{
  /*
   * hack for now
   */
  char *hs = "", *ws = "", *as = "" ;
  printf ("<select name='db'>\n") ;
  if (strcmp (db, "human") == 0)
    hs = "selected" ;
  if (strcmp (db, "worm") == 0)
    ws = "selected" ;
  if (strcmp (db, "ara") == 0)
    as = "selected" ;
  printf ("<option value='human' %s>Human</option>\n", hs) ; 
  printf ("<option value='worm' %s>C elegans</option>\n", ws) ;
  printf ("<option value='ara' %s>A thaliana</option>\n", as) ;
  printf ("</select>\n") ;
} /*  */

/*******************************************************************/

static char *db_name (void)
{
  /*
   * hack for now
   */
  if (strcmp (db, "human") == 0)
    return "Human" ;
  if (strcmp (db, "worm") == 0)
    return "C elegans" ;
  if (strcmp (db, "ara") == 0)
    return "A thaliana" ;
  return "?" ;
} /*  */

/*******************************************************************/
/*******************************************************************/

static void page_top (char *query, char *db)
{
  printf (
	  "<html>\n"
	  "<head><title>ZONE search</title></head>\n"
	  "<body bgcolor='white'>\n"
	  "<font face='verdana'>\n"
	  "<A name=top href='/' target='_top'>"
	  "<IMG src='/corehtml/left.GIF' width=130 border=0>"
	  "</A>"
	  "<br>\n"
	  ) ;
  printf (
	  "<a href='index.html?%s' target='_top'><IMG src='images/aceview.gif' width=140 border=0></a>"
	  "<br>\n",
	  db
	  ) ;
} /*  */

/*******************************************************************/

static void page_first_time (void)
{
  printf (
	  "\n"
	  "<FORM method=get action='zone.cgi'>\n"
	  ) ;
  
  printf (
	  "\n"
	  "<h2>Zone search\n"
	  "</h2>\n"
	  ) ;
  printf ("<a href='help.%s.html#zone_search'>Help</a>\n",my_db_desc->db_version) ;
  printf (
	  "<p>Search chromosome: "
	  ) ;
  
  printf (
	  "<INPUT type=text name=chrom size=3>\n") ;
  
  printf (" from bp ") ;
  printf (
	  "<INPUT type=text name=mini size=9>\n") ;
  
   printf (" to bp ") ;
   printf (
	  "<INPUT type=text name=maxi size=9>\n") ;
  
  printf (
	  "The chromosome and the mini and maxi coordinates in bp."
	  "<p>For example, search 3, 1415926 5353239"
	  "<p>in:\n"
	  ) ;
  db_select (db) ;
  printf (
	  "<br>\n"
	  "<INPUT type=submit name=search value='Search'>\n"
	  ) ;
  
  printf ("</form>\n") ;
  printf ("</body></html>\n") ;
} /*  */

/*******************************************************************/

static int f_row = 1 ;	/* first row to display */
static int n_row = 50 ;	/* how many rows to display */

static void table_top (TBL *tbl, /* char *name, char *num, char *def, char *prod, char *genes, */
		       int tot_prod, int tot_gene)
{
  int i ;
  
  printf ("<a href='zone.cgi?db=%s'>New search</a>\n", db) ;
  printf ("<br><a href='help.%s.html#zone_search'>Help</a>",my_db_desc->db_version) ;
  
  printf (
	  "<h2>Search Results Found %d genes in %s</h2>\n", tbl->tbl->rows, db_name ()) ;
  
  printf (
	  "<form method=get action='zone.cgi'>\n"
	  ) ;
  
  printf ("<p>") ;
  printf ("<input type=hidden name=context value='%s'>\n", context_name) ;
  printf ("<input type=hidden name=db value='%s'>\n", db) ;
  
  printf ("<input type=submit name=search value='Refine results'><p>"
	  "Query: <INPUT type=text name=chrom size=3><br>\n"
	  ) ;
  
  /*
   * Danielle wants to suppress this input box if there are fewer than
   * 10 rows of data.  You also need it if the number of rows is
   * greater than the number we are displaying, even if the number we
   * display is less than 10.  You can't use just n_rows > n_row
   * because you want to allow them to reduce the number of rows if
   * they set it to 300 or something.
   */
  if ( (tbl->tbl->rows > 10) || (tbl->tbl->rows > n_row)) 
    {
      printf ("<p>Show <input type=text size=5 name=n_row value=%d>\n", n_row) ;
      printf (" rows starting at row <input type=text size=6 name=f_row value=%d ><br>\n", f_row) ;
      printf ("<input type=submit name=prev value='Previous Page'> \n") ;
      printf ("<input type=submit name=next value='Next Page'> \n") ;
    }
  printf ("<p>") ;
  
  printf (
	  "<p>\n"
	  "<table width=100%% border>\n"
	  "\n"
	  ) ;
  
  
  printf (
	  "<tr>\n"
	  ) ;
  printf ("<th valign='bottom'>remove<br>checked</th>\n") ;
  
  /*
   * the nbsp prevents some versions of internet explorer from putting
   * a line break between the check box and the text.
   */
  printf ("<th align='left'>") ;
  printf ("<INPUT type='radio' name='n_name' value='y' checked>&nbsp ;include<BR>") ;
  printf ("<INPUT type='radio' name='n_name' value='n'>&nbsp ;exclude<BR>") ;
  printf ("<INPUT type=text name=a_name size=10><br>") ;
  printf ("</th>\n") ;
  
  printf ("<th align='left'>") ;
  printf ("<INPUT type='radio' name='n_num' value='y' checked>&nbsp ;include<BR>") ;
  printf ("<INPUT type='radio' name='n_num' value='n'>&nbsp ;exclude<BR>") ;
  printf ("<INPUT type=text name=a_num  size=8 ><br>") ;
  printf ("</th>\n") ;
  
  printf ("<th align='left'>") ;
  printf ("<INPUT type='radio' name='n_def' value='y' checked>&nbsp ;include<BR>") ;
  printf ("<INPUT type='radio' name='n_def' value='n'>&nbsp ;exclude<BR>") ;
  printf ("<INPUT type=text name=a_def  size=50 ><br>") ;
  printf ("</th>\n") ;
  
  printf ("<th valign='bottom'>") ;
  printf ("<INPUT type='radio' name='n_prod' value='g' checked>&nbsp ;&gt ;<BR>") ;
  printf ("<INPUT type='radio' name='n_prod' value='e' >&nbsp ;=<BR>") ;
  printf ("<INPUT type='radio' name='n_prod' value='l' >&nbsp ;&lt ;<BR>") ;
  printf ("<INPUT type=text name=a_prod size=5 >") ;
  printf ("</th>\n") ;
  
  printf ("<th valign='bottom'>") ;
  printf ("<INPUT type='radio' name='n_gene' value='g' checked>&nbsp ;&gt ;<BR>") ;
  printf ("<INPUT type='radio' name='n_gene' value='e' >&nbsp ;=<BR>") ;
  printf ("<INPUT type='radio' name='n_gene' value='l' >&nbsp ;&lt ;<BR>") ;
  printf ("<INPUT type=text name=a_gene size=5 >") ;
  printf ("</th>\n") ;
  
  printf (
	  "</tr>\n"
	  "\n"
	  ) ;
  printf (
	  "<tr>\n"
	  "<th> &nbsp ; </th>\n"
	  "<th>%d<br>%s</th>\n"
	  , tbl->tbl->rows, tbl->title[0]) ;

  for (i = 1 ; i < tbl->cols ; i++)
    printf ("<th>%s</th>\n", tbl->title[i]) ;

  printf ("</tr>\n") ;
} /*  */

/*******************************************************************/

static void announce_download (void)
{
  
  printf ("<p>You can download this table and related data.<br>") ;
  printf (
	  "<br>Download the table of results <a href=\"zone.cgi?db=%s&download=tt&context=%s\">as tab delimited file</a>\n"
	  , db, context_name) ;
  
  /*
    printf ("<br>Download list of product names "
    "<a href=\"zone.cgi?db=%s&download=ap&context=%s\">as AceDB file</a>\n", db, context_name) ;
    
    printf ("<br>Download list of gene names "
    "<a href=\"zone.cgi?db=%s&download=ag&context=%s\">as AceDB file</a>\n", db, context_name) ;
  */

  printf ("<br><a href='zone_downloads.%s.html'>About Zone query downloads</a>\n",my_db_desc->db_version) ;
} /*  */

/*******************************************************************/

static void clean_query (char *s)
{
  for ( ; *s ; s++) 
    {
      switch (*s) 
	{
	case '"':
	case '\'':
	case '%':
	  *s = '?' ;
	  break ;
	}
    }
} /* clean_query */

/*******************************************************************/

static void general_query (char *chrom, int mini, int maxi)
{
  char *q ;

  if (db_context_keyset)
    {
      q = messprintf ("IS \"%s*\" && IS <= \"%s%d\" && IS >= \"%s%d\" "
		      , chrom, mini, maxi) ;
      db_context_keyset = ac_ksquery_keyset ( db_context_keyset, q, NULL) ;
    }
  else
    {
      q = messprintf ("Find NewName IS \"%s*\" && IS <= \"%s%d\" && IS >= \"%s%d\" "
		      , chrom, mini, maxi) ;
      db_context_keyset = ac_dbquery_keyset ( db_descriptor, q, NULL) ;
    }

} /* general_query */

/*******************************************************************/

static char *numeric_filter (char *q)
{
  char *s ;
  s = cgi_clean_arg (q) ;
  if (!s)
    return ">" ;
  if (*s == 'g')
    return ">" ;
  if (*s == 'e')
    return "=" ;
  if (*s == 'l')
    return "<" ;
  return ">" ;
} /* numeric_filter */

/*******************************************************************/

static char not_char (char *q)
{
  char *s = cgi_clean_arg (q) ;
  if (s && *s == 'n')
    return '!' ;
  return ' ' ;
} /* not_char */

/*******************************************************************/

static void filter_query (void)
{
  char *q, *s, *t ;
  char            not ;
  
  q = "" ;
  
  s = cgi_clean_arg ("a_name") ;
  if (s && *s) 
    {
      s = strdup (s) ;
      not = not_char ("n_name") ;
      clean_query (s) ;
      q = strdup (messprintf ("%s ( %c ( IS \"*%s*\" ) ) ; ", q, not, s)) ;
    }
  s = cgi_clean_arg ("a_num") ;
  if (s && *s) 
    {
      s = strdup (s) ;
      clean_query (s) ;
      not = not_char ("n_num") ;
      q = strdup (messprintf ("%s ( %c ( accession = \"*%s*\" ) ) ; ", q, not, s)) ;
    }
  s = cgi_clean_arg ("a_def") ;
  if (s && *s) 
    {
      s = strdup (s) ;
      clean_query (s) ;
      not = not_char ("n_def") ;
      q = strdup (messprintf ("%s ( %c ( definition = \"*%s*\" ) ) ; ", q, not, s)) ;
    }
  s = cgi_clean_arg ("a_prod") ;
  if (s && *s) 
    {
      t = numeric_filter ("n_prod") ;
      if (t)
	q = strdup (messprintf ("%s ( COUNT product %s %s ) ; ", q, t, s)) ;
    }
  s = cgi_clean_arg ("a_gene") ;
  if (s && *s) 
    {
      t = numeric_filter ("n_gene") ;
      if (t)
	q = strdup (messprintf ("%s ( COUNT gene %s %s ) ; ", q, t, s)) ;
    }
  if (*q) 
    {
      if (db_context_keyset)
	db_context_keyset = ac_ksquery_keyset (db_context_keyset, q, NULL) ;
      else
	db_context_keyset = ac_dbquery_keyset (db_descriptor, q, NULL) ;
    }
} /*  */

/*******************************************************************/

static void download (TBL *tbl, char *download)
{
  unsigned char *r ;
  const char *err ;
  char *q ;
  
  if (!context_name)
    return ;
  
  db_descriptor = ac_open_db (my_db_desc->db_access , &err) ;
  if (! db_descriptor) 
    {
      printf ("<h1>error opening database: %s\n", err) ;
      return ;
    }
  
  load_context () ;
  switch (download[0])
    
    {
    case 'a':
      switch (download[1]) 
	{
	case 'p':
	  q = "show -a product" ;
	  break ;
	case 'g':
	  q = "show -a gene" ;
	  break ;
	default:
	  return ;
	}
      printf ("// %s\n", q) ;
      r = ac_command (db_descriptor, q, NULL, NULL) ;
      printf ("%s\n", r) ;
      return ;
    case 't':
      switch (download[1]) 
	{
	case 't':
	  /*
	    ac_tablemaker_table (db_descriptor, tbl->nam,
		       initial_keyset, ac_tablemaker_db, 0, 0) ;
	  */  
	  q = messprintf ("table -active -n %s", tbl->nam) ;
	  r = ac_command (db_descriptor, q, NULL, NULL) ;

	  printf ("%s\n", r) ;
	  return ;
	}
    }
} /* download */

/*******************************************************************/

static char *html_safe (const char *s)
{
  static char *save ;
  if (save)
    free (save) ;
  if (! s || ! *s )
    
    {
      save = 0 ;
      return "&nbsp ;" ;
    }
  save = html_encode (s) ;
  return save ;
} /*  */

/*******************************************************************/

static char *url_safe (const char *s)
{
  static char *save ;
  if (save)
    free (save) ;
  if (! s || ! *s)
    
    {
      save = 0 ;
      return "" ;
    }
  save = url_encode (s) ;
  return save ;
} /* url_safe */

/*******************************************************************/

static void display_current (TBL *tbl)
{
  int row ;
  int tot_prod, tot_gene ;
  int prod, gene, end_row ;
  AC_KEYSET	ks_p, ks_g ;
  
  tbl->tbl = ac_tablemaker_table ( db_descriptor, tbl->nam,
			     db_context_keyset, ac_tablemaker_db, NULL, NULL, NULL, 0) ;
  
  ks_p = ac_new_keyset (db_descriptor, 0) ;
  /* ac_ksquery_keyset (db_context_keyset, "follow product",  NULL) ; */
  
  tot_prod = ac_keyset_count (ks_p) ;
  
  ks_g = ac_copy_keyset (db_context_keyset, 0) ;
  
  tot_gene = ac_keyset_count (ks_g) ;
  
  table_top (tbl, tot_prod, tot_gene) ;
  
  end_row = f_row + n_row ;
  if (tbl->tbl->rows+1 < end_row)
    end_row = tbl->tbl->rows+1 ;
  
  for (row = f_row ; row < end_row ; row++)
    
    {
      char *target ;
      const char *html, *url ;
      printf ("<tr>\n") ;
      
      /*
       * row number, checkbox
       * zone name
       */
      
      url = html = ac_table_printable (tbl->tbl, row-1, 0, NULL) ;
      html = html_safe (html) ;
      url = url_safe (url) ;
      target = strdup (url) ;
      
      printf ("<td>%d <input type=checkbox name=del value=%s></td>",
	      row, url) ;
      printf ("<td>%s</td>", html) ;
      
      /*
       * accession number
       */
      url = html = ac_table_printable (tbl->tbl, row-1, 1, NULL) ;
      html = html_safe (html) ;
      url = url_safe (url) ;
      
      printf ("<td><a href=http://www.sanger.ac.uk/cgi-bin/Zone/getacc?%s>%s</a></td>", url, html) ;
      
      /*
       * definition
       */
      html = ac_table_printable (tbl->tbl, row-1, 2, NULL) ;
      html = html_safe (html) ;
      printf ("<td>%s</td>",html) ;
      
      /*
       * products
       */
      prod = ac_table_int (tbl->tbl, row-1, 3, 0) ;
      printf ("<td>%d</td>",prod) ;
      
      /*
       * genes
       */
      gene = ac_table_int (tbl->tbl, row-1, 4, 0) ;
      printf ("<td><a href=\"av.cgi?db=%s&a=list&c=zone&q=%s\">%d</a></db>", db, target, gene) ;
      
      printf ("</tr>\n") ;
    }
  
  printf ("</table>\n") ;
  
  printf ("<input type=submit name=search value='Refine results'><p>") ;
  
  announce_download () ;
} /* display_current */

/*******************************************************************/
/*
 * Delete a row from the results.
 * 
 * It might be nice to avoid doing a database transaction for each row we
 * delete, but I think this is fast enough and uncommon enough that it does
 * not really matter.
 */
static void keydel (char *value, void *arg)
{
  char *s ;
  
  s = messprintf ("! IS %s", value) ;
  db_context_keyset = ac_ksquery_keyset (db_context_keyset, s, NULL) ;
} /* keydel */

/*******************************************************************/
/*******************************************************************/

static int func_search (TBL *tbl, char *chrom, int mini, int maxi)
{
  const char *err ;
  
  printf ("Content-type: text/html\n\n") ;
  
  db_descriptor = ac_open_db ( my_db_desc->db_access , &err ) ;
  if ( ! db_descriptor ) 
    {
      printf ("<h1>error opening database: %s\n", err) ;
      return 0 ;
    }
  if (context_name) 
    {
      load_context () ;
    } 
  else 
    {
      create_context () ;
      db_context_keyset = ac_new_keyset (db_descriptor, NULL) ;
    }
  
  /* clean_query (user_query) ; */
  
  general_query (chrom, mini, maxi) ;
  
  filter_query () ;
  
  cgi_call ("del", 0, keydel) ;
  
  save_context () ;
  
  page_top ("", db) ;
  
  display_current (tbl) ;
  
  return 0 ;
} /* func_search */

/*******************************************************************/

int main ()
{
  char *s ;
  int n ;
  TBL *tbl = 0 ;
  char *chrom ;
  int mini=0, maxi=0 ;

  setbuf (stdout, NULL) ;
  n = cgi_init ("zonecgi.log") ;
  if (n < 0) 
    {
      printf ("Content-type: text/html\n\n") ;
      printf ("<h1>cgi_init -> %d\n</h1>", n) ;
      exit (1) ;
    }
 
  db = cgi_clean_arg ("db") ;
  if (!db)
    db = "worm" ;
  
  tbl = tblCreate ("webzone") ;
  /*
   * does not really open the database - just find it in the
   * list so we can get a "v8" or "v9" to put in URL for the
   * help file.
   */
  my_db_desc = get_db_desc (db, &s) ;
  if (!my_db_desc) 
    {
      printf ("content-type: text/html\n\n") ;
      printf ("database error - %s\n",s) ;
      exit (1) ;
    }
#if 0
  cgi_dump () ;
 printf ("content-type: text/plain\n\n") ;
  printf ("version %s\n\n", my_db_desc->db_version) ;
  printf ("name    %s\n\n", my_db_desc->db_name) ;
  printf ("access  %s\n\n", my_db_desc->db_access) ;
  printf ("spec_c  %s\n\n", my_db_desc->db_species_char) ;
  printf ("s_titl  %s\n\n", my_db_desc->db_species_title) ;
  printf ("db_t    %s\n\n", my_db_desc->db_title) ;
#endif
  
  
  s = cgi_clean_arg ("n_row") ;
  if (s)
    n_row = atoi (s) ;
  if (n_row < 1)
    n_row = 1 ;
  
  s = cgi_clean_arg ("f_row") ;
  if (s)
    f_row = atoi (s) ;
  if (f_row < 1)
    f_row = 1 ;
  
  chrom = cgi_clean_arg ("chrom") ;

  s = cgi_clean_arg ("mini") ;
  if (s)
    mini = atoi (s) ;
  if (mini < 1)
    mini = 1 ;

  s = cgi_clean_arg ("maxi") ;
  if (s)
    maxi = atoi (s) ;
  if (maxi < 1)
    maxi = 1 ;
  
  context_name = cgi_clean_arg ("context") ;
  
  if (cgi_clean_arg ("clear")) 
    {
      printf ("Content-type: text/html\n\n") ;
      page_top ("", db) ;
      page_first_time () ;
      return 0 ;
    }
  s = cgi_clean_arg ("search") ;
  if (0 && s)
    return func_search (tbl, chrom, mini, maxi) ;
  
  
  s = cgi_clean_arg ("prev") ;
  if (s) 
    {
      f_row = f_row - n_row ;
      if (f_row < 1)
	f_row = 1 ;
      func_search (tbl, chrom, mini, maxi) ;
      return 0 ;
    }
  s = cgi_clean_arg ("next") ;
  if (s) 
    {
      f_row = f_row + n_row ;
      func_search (tbl, chrom, mini, maxi) ;
      return 0 ;
    }
  s = cgi_clean_arg ("download") ;
  if (s) 
    {
      printf ("Content-type: text/plain\n\n") ;
      download (tbl, s) ;
      return 0 ;
    }
  if (chrom)
    return func_search (tbl, chrom, mini, maxi) ;
  
  printf ("Content-type: text/html\n\n") ;
  page_top ("", db) ;
  page_first_time () ;
  return 0 ;
} /* main */

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

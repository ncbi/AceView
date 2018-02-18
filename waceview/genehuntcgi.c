/*
 * File genehuntcgi.c
 * must be compiled as genehunt.cgi and placed in the
 * AceView directory to be found by the ncbi_banner.html
 * 
 * Users guide: see AceView/help/gene_hunt_help.html
 * 
 * Programmers guide:
 * 
 * This CGI export an .html document made of 2 parts
 *  -top part: a selector with text entry fields to let
 *   the user ask for a set of genes
 *  -bottom part: a possible long table with the genes
 * 
 * The table has to be precomputed on the acedb server and
 * exported in the AceView directory under the name
 *   -GeneProperties.txt
 * The table definition: GeneProperties.def
 * must be in synch with the menu in this code
 *
 * The code looks for GeneProperties.dat
 * If this file is absent it is created 
 * GeneProperties.dat is a compiled view of GeneProperties.txt
 *
 * Then the code parses the cgi parameters and reads GeneProperties.dat
 * keeping only the lines fitting the user choice
 * 
 * Finally a table slice is exported
 */
/*

#include <unistd.h>
#include <wh/regular.h>
*/

#include <ctype.h>
#include <waceview/cgi.h>
#include <waceview/avlib.h>
#include <vtxt.h>
#include <aceio.h>

static Stack s = 0 ;
static Stack sFlag = 0 ;
static Array tt = 0 ;
static char *chroms[] = { "*", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", 0} ;
int selectedChrom = 0 ;

/**************************************************************/
/**************************************************************/
/**************************************************************/
/*
  "3PAP"  "5_32358576"    "54545" 7       1       4       NULL    "Genetic"
  20      "5p13.3"        NULL    125965  "3PAP.aDec03"     2325    2310    128
*/
/*
  display: 0/1/2, 0:no_show, 1 show, 2:show and hot-link this column in the exported table
  select:  0/1, present in the hand selector with text-entry fields
  type: 0=number, 1=text 2=BOOL 3=codon 4=chrom 5=from 6=to 100=menu
  value: the value selected by the user (number or stackMark for a text)
  title: used at op of table and as URL parameter
  longTitle: use to name the slector field
*/

typedef struct selectStruct { int display, select, type, value ; char *title, *longTitle, *prompt ; } MM ;
static MM menu[] = {
  { 2, 0, 1, 0, "gene", "Gene identifier", "*" },
  { 1, 1, 104, 0, "newName", "Chromosome", "0: any, 1,2,3...X,Y" },
  { 0, 1, 5, 0, "a1", "...from", "0: any, 1,2,... minimal position in Mb (million bp)" },  
  { 0, 1, 6, 0, "a2", "...to", "0: any, 1,2,... maximal position in Mb (million bp)" },  
  { 1, 1, 1, 0, "locusId", "LocusId", "*=optional,**=known,-1=none,*5123*=desired value" },
  { 1, 1, 0, 0, "nMrna", "Number of variants", "0,1,2... : desiredminimal value" },
  { 1, 1, 0, 0, "nNm", "Number of NM", "0:optional; 1,2... : desired minimal value, -1:none"   },
  { 0, 1, 0, 0, "nPaper", "Number of papers" , "0:optional; 1,2... : desired minimal value, -1:none"   },
  { 1, 1, 0, 0, "nPfam", "Number of PFAM motifs" , "0:optional; 1,2... : desired minimal value, -1:none"   },
  { 1, 0, 2, 0, "genetic", "Genetic phenotype" , "0,y,n : 0=optional, y=something known, n=unknown" },
  { 1, 1, 0, 0, "nGt_ag", "Number of gt_ag introns",  "0:optional; 1,2... : desired minimal value, -1:none"},
  { 0, 0, 1, 0, "chromBand", "Chromosome band", "*=optional,**=known, 3Q*=desired value" },
  { 0, 1, 2, 0, "cloud", "Cloud gene", "0,y,n : 0=optional, y=cloud gene, n=normal gene" },
  { 0, 0, 0, 0, "covers", "Extent on genome", "0,1,2, ... : Minimal extent in bp" },
  { 0, 0, 1, 0, "mrnaName" , "Name of variant", "*" },
  { 1, 1, 3, 0, "orf", "Length of ORF (AA)" , "0,1,2, ... : Minimal number of amino acids" },
  { 1, 1, 3, 0, "cds", "Length of CDS (AA)" , "0,1,2, ... : Minimal number of amino acids" },
  { 1, 0, 0, 0, "nClones", "Number of supporting clones", "0,1,2.. : desired minimal value."  },
  { 1, 1, 1, 0, "title" , "Title", "*, *atpase*" },
  { 0, 0, 0, 0, 0} 
} ;

/*****************************************************************/

static void selectSubMenu (vTXT layerContent, MM *mm, int jj)
{
  int nn = 0 ;

  if (jj) vtxtPrintf (layerContent, ",\n ") ; 

  if (!nn++)
    vtxtPrint (layerContent, "\"<table bgcolor='#ebfff5' border=2>") ;
  if (FALSE)
    vtxtPrintf (layerContent, "  <tr><td><a href='#%s'>%s</a></td></tr>", mm->title, mm->longTitle) ;
  else
    vtxtPrintf (layerContent, "  <tr><td> %s</td></tr>", mm->longTitle) ;

  if (nn) 
    vtxtPrint (layerContent, "</table>\"\n") ;
  else
    vtxtPrint (layerContent, "\"\"") ;
} /* gmpJumpPointShowSubMenu */

/**************************************************************/

static void selectDoShow (vTXT blk, vTXT layerContent)
{
  /*
   * hack for now
   */
  int i, jj = 0 ;
  MM *mm = menu - 1 ;

  /*
   * the nbsp prevents some versions of internet explorer from putting
   * a line break between the check box and the text.
   */
  vtxtPrint (blk, "<script language=\'JavaScript\'>") ;
  vtxtPrint (blk, "document.write(layerCreate ('gmpPopUpLayer','',''));") ;
  vtxtPrint (blk, "</script>\n") ;
  vtxtPrint (blk, "  <font size=12>") ;
  vtxtPrint (blk, "    <table  bgcolor='#dffff2' border=2 width=180>\n") ;
  
  vtxtPrint (layerContent, "<script language = 'JavaScript'>\n<!--\n") ;
  vtxtPrint (layerContent, "var gmpPopUpMenu = new Array (\n") ; /* open gmpMenu Array */
  
  vtxtPrint (blk, "<tr bgcolor=#afafff><td>Feature</td><td>Value</td><td>Possible values</td></tr>") ;
  for (mm = menu ; mm->title ; mm++)
    if (mm->select)
      {
	vtxtPrint (blk, "<tr>") ;
	vtxtPrintf (blk, "      <tr><td  onMouseOver='gmpPopUp (event,%d);'>", jj) ;
	vtxtPrint (blk, mm->longTitle) ;
	vtxtPrint (blk, "</td>\n") ;
	selectSubMenu (layerContent, mm, jj++) ; /* fill gmpMenu Array */


	switch (mm->type)
	  {
	  case 0: /* NUMBER */
	  case 3: /* CODON */
	  case 5: /* a1 */
	  case 6: /* a2 */
	    vtxtPrintf(blk, "<td><INPUT type=text name='%s' size=10 value=%d></td>", mm->title, mm->value);
	    break ;
	  case 2:
	    vtxtPrintf(blk, "<td><INPUT type=text name='%s' size=10 value=%s></td>"
		   , mm->title
		   , mm->value ? (mm->value > 0  ?"y" : "n") : "0"
		   ) ;
	    break ;
	  case 1: /* TEXT */
	  case 4: /* chrom */
	    vtxtPrintf(blk, "<td><INPUT type=text name='%s' size=10 value=%s></td>"
		   , mm->title, mm->value ? (mm->value != -1 ? stackText (sFlag,mm->value) : "-1") : "*");
	    break ;
	  case 104: /* chrom */	
	    vtxtPrint (blk, "<td><SELECT name=newName>") ;
	    for (i = 0 ; chroms[i] ; i++)
	      vtxtPrintf(blk, "<OPTION value=%s %s>%s</OPTION>"
			 , chroms[i]
			 , i == selectedChrom ? "selected" : ""
			 , chroms[i]) ;
	    vtxtPrintf(blk, "</SELECT></td>") ;
	    break ;
	  }
	vtxtPrintf (blk, "<td>%s</td>\n", mm->prompt) ;
	vtxtPrint (blk, "</tr>\n") ;
      }    
  vtxtPrint (layerContent, ");\n") ; /* close gmpMenu Array */
  vtxtPrint (layerContent, "//-->\n</script>\n") ;

  vtxtPrint (blk, "<table>\n") ;   
  vtxtPrint (blk, "  </font>\n") ;


} /* selectDoShow */

/**************************************************************/

static void selectShow (void)
{
  vTXT blk = vtxtCreate () ;
  vTXT layerContent = vtxtCreate () ;

  selectDoShow (blk, layerContent) ;

  /* create top anchor */
  printf ("<a href='#' name='top'></a>") ;

  printf ("%s\n", vtxtPtr (layerContent)) ; /* first declare gmpMenu */

  printf ("<table border=0>\n") ;
  printf ("  <tr>\n") ;
  printf ("    <td>\n      ") ;
  printf ("%s", vtxtPtr (blk)) ;          /* then use it, hence we need the two vTXT */
  printf ("\n    </td>\n") ;
  printf ("    <td>\n      ") ;
  printf ("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;") ;
  printf ("\n    </td>\n") ;

  printf ("  </tr>\n") ;
  printf ("</table>\n") ;

  vtxtDestroy (layerContent) ;
  vtxtDestroy (blk) ;
} /* selectShow */

/**************************************************************/
/* parse flags parameter assumed to look like
   0.0.1.0.23.0.542.1
   with one integer for each tunable flag
*/

static void flagSet (void)
{
  MM *mm = menu ;
  char *cp ;
  int n, position = 0 ;

  sFlag = stackCreate (1000) ;
  stackTextOnly (sFlag) ;
  pushText (sFlag, "0") ;

  for (mm = menu ; mm->title ; mm++)
    if (mm->select)
      {
	mm->value = 0 ; 
	if ((cp = cgi_clean_arg (mm->title)))
	  {
	    switch (mm->type)
	      {
	      case 0: /* NUMBER */
	      case 3: /* CODON */
		if (sscanf (cp, "%d", &n) == 1)
		  mm->value = n ;   
		break ;
	      case 2:/* BOOL */
		if (!strcasecmp (cp, "y")) /* BOOL */
		  mm->value = 1 ;
		else if (!strcasecmp (cp, "n"))
		  mm->value = -1 ;
		break ;
	      case 1: /* TEXT */
	      case 4: /* CHROM */
		if (!strcmp (cp, "") || !strcmp (cp, "*") || !strcmp (cp, "0")) ;
		else if (!strcmp (cp, "-1"))
		  mm->value = -1 ;
		else
		  {
		    position = mm->value = stackMark (sFlag) ;
		    pushText (sFlag, cp) ;
		  }
		break ;
	      case 104: /* CHROM */
		if (!strcmp (cp, "") || !strcmp (cp, "*") || !strcmp (cp, "0") ||
		    !strcmp (cp, "-1"))
		  selectedChrom = mm->value = -1 ;
		else
		  {
		    int i ;

		    position = mm->value = stackMark (sFlag) ;
		    pushText (sFlag, cp) ;
		    for (i = 0 ; chroms[i] ; i++)
		      if (!strcasecmp (cp, chroms[i]))
			selectedChrom = i ;
		  }
		break ;
	      case 5: /* a1 */
	      case 6: /* a2 */
		if (position && sscanf (cp, "%d", &n) == 1)
		  mm->value = n ;   
		break ;
	      }
	  }
      }
} /* flagSet */

/**************************************************************/

static int tableLoad (void)
{
  int nn = 0, line = 0, lastNn, ok, k, pos ;
  char buffer [1024] ;
  ACEIN fi ;
  MM *mm ;

  fi = aceInCreateFromFile ("GeneProperties.txt", "r", 0, 0) ;
  if (!fi)
    {
      printf ("cannot open GeneProperties.txt, sorry<br>") ;
      return 0 ;
    }
  s = stackCreate (100000000) ;
  stackTextOnly (s) ;

  pushText (s, "0") ; /* avoid zero */
  tt = arrayCreate (200000, KEY) ;
  while (aceInCard (fi))
    {
      line++ ; 
      lastNn = arrayMax (tt) ;
      pos = -1 ;
      for (mm = menu, ok = TRUE ; ok && mm->title ; mm++) /* grab NCOL entries */
	{
	  if (mm->type == 5) /* a1 */
	    {
	      if (mm->value > 0 && (pos == -1 || mm->value * 1000000 > pos))
		ok = FALSE ; 
	      continue ; /* do not scan the table, do not increment nn */
	    }
	  else if (mm->type == 6) /* a1 */
	    {
	      if (mm->value > 0 && (pos == -1 || mm->value * 1000000  < pos))
		ok = FALSE ; 
	      continue ; /* do not scan the table, do not increment nn */
	    }
	 if (mm->type == 0 ||  /* NUMBER */
	      mm->type == 3) /* CODON */
	    {
	      if (! aceInInt (fi, &k))
		{ nn = lastNn ; goto abort ; }
	      if (k == UT_NON_INT) 
		k = 0 ; /* parsed a "NULL" */
	      if (mm->select && ((mm->value > (mm->type == 3 ? k/3 : k)) || (k && mm->value == -1)))
		ok = FALSE ;  /* do not break, finish parsing the line */
	    }
	  else  /* TEXT BOOL CHROM */
	    {
	      if (! aceInWordScan (fi, buffer, 1024))
		{ nn = lastNn ; goto abort ; }
	      if (! strcmp (buffer, "NULL"))
		k = 0 ;
	      else
		{
		  k = stackMark (s) ;
		  pushText (s, buffer) ;
		}
	      
	      if (ok && mm->select && mm->value)
		{
		  switch (mm->type)
		    {
		    case 1: /* TEXT */
		      if (
			  (mm->value == -1 && k > 0) ||
			  (mm->value != -1 &&
			   (!k || !pickMatch(buffer, stackText (sFlag, mm->value))))
			  )
			ok = FALSE ;
		      break ;
		    case 2: /* BOOL */
		      if ((!k && mm->value == 1) || (k && mm->value == -1))
			ok = FALSE ;
		      break ; 
		    case 4: /* chrom */
		      if (k)
			{
			  char cc, *cp = buffer ;
			  while (*cp && *cp != '|' && *cp != '_') cp++ ;
			  cc = *cp ; *cp = 0 ;
			  if (strcasecmp (buffer, stackText (sFlag, mm->value)))
			    ok = FALSE ;
			  if (cc && sscanf (cp+1, "%d", &pos) == 1) ;
			  else pos = -1 ;
			  *cp = cc ;
			}
		      break ;
		    case 104: /* chrom menu */
		      if (k && selectedChrom)
			{
			  char cc, *cp = buffer ;
			  while (*cp && *cp != '|' && *cp != '_') cp++ ;
			  cc = *cp ; *cp = 0 ;
			  if (strcasecmp (buffer, chroms[selectedChrom]))
			    ok = FALSE ;
			  *cp = cc ;
			}
		      break ;
		    }
		}
	    }
	  array (tt, nn++, int)  = k ; /* always increment the cell number nn */
	}
      if (!ok)
	nn = arrayMax (tt) = lastNn ;
    }
 abort:
  aceInDestroy (fi) ;
  for (mm = menu, k = 0 ; mm->title ; mm++) k++ ;
  return nn/k ;
} /* tableLoad */

/**************************************************************/

static void tableShow (char *db)
{
  MM *mm ;
  int nn = 0, bb ;
  int line = 0 ;
  char *cp ;

  printf ("<p><table>") ; 

  printf ("<tr bgcolor=#afafff>") ;
  for (mm = menu ; mm->title ; mm++) 
    if (mm->display)
      printf ("<td>%s</td>", mm->title) ;
  printf ("</tr>") ;

  while (nn < arrayMax (tt))
    {
      if ((line++) % 2)
	printf ("<tr bgcolor=#efefff>") ;
      else
	printf ("<tr bgcolor=#ffffff>") ;

      for (mm = menu ; mm->title ; mm++)
	{
	  if (mm->type == 5 || mm->type == 6)
	    continue ; /* do not increment nn */
	  switch (mm->display)
	    {
	    case 0:
	      break ;
	    case 1:
	      switch (mm->type)
		{
		case 0: /* NUMBER */
		  printf ("<td>%d</td>", arr (tt, nn, int)) ;
		  break ;
		case 1: /* TEXT */
		case 4: /* CHROM */
		  printf ("<td>%s</td>", stackText (s, arr (tt, nn, int))) ;
		  break ;
		case 2: /* BOOL */
		  bb = arr (tt, nn, int) ;
		  printf ("<td>%s</td>", bb  ? (bb > 0 ? "yes" : "no") : "0") ;
		  break ;
		case 3: /* CODON */
		  printf ("<td>%d</td>", arr (tt, nn, int)/3) ;
		  break ;
		}
	      break ;
	    case 2:
	      cp = stackText (s, arr (tt, nn, int)) ;
	      printf ("<td><a href=\"//www.humangenes.org/av.cgi?db=%s&q=%s\">%s</a></td>"
		      , db, cp, cp) ;
	      break ;
	    }
	  nn++ ;
	}
      printf ("</tr>") ;
      if (line > 30)
	break ;
    }
  printf ("</table><p>") ;
} /* tableShow */

/**************************************************************/
/**************************************************************/

static void page_top (char *db)
{
  printf(
	 "Content-type: text/html\n\n"
	 "<html>\n"
	 "<head>\n<title>Gene Hunt</title>\n"
	 "<script type=\'text/javascript\' src='jscript/basic.js'></script>"
	 "<script type=\'text/javascript\' src='jscript/main2.js'></script>"
	 "</head>\n"
	 "<body bgcolor='white'>\n"
	 "<font face='verdana'>\n"
	 "<A name=top href='/' target='_top'>"
	 "<IMG src='/corehtml/left.GIF' width=130 border=0>"
	 "</A>"
	 "<br>\n"
	 );
  printf ("<script language=\'JavaScript\'>\n"
          "<!--\n"
          "gInitAll() ;\n"
          "//-->\n</script>\n") ;

  printf(
	 "<a href='index.html?%s' target='_top'><IMG src='images/aceview.gif' width=140 border=0></a>"
	 "<br>\n",
	 db
	 );
  printf(
	 "<h2>Query by Features: \n"
	 "<a href='HelpQuery.html#genehunt_search'>Help</a>\n"
	 "</h2>\n"
	 );

  printf(
	 "\n"
	 "<FORM method=get action='genehunt.cgi'>\n"
	 );
  
  selectShow() ;
  printf(
	 "<br>\n"
	 "<INPUT type=submit name=search value='Search'>\n"
	 );
  

  printf("</FORM>\n") ;
} /* page_top */

/**************************************************************/

static void page_end (void)
{
  printf ("</body>\n</html>\n") ;
} /* page_end */

/**************************************************************/
/**************************************************************/

static void table_top(char *db, char *flags, char *name, char *num, char *def, char *prod, char *genes, int n_rows, int tot_prod, int tot_gene)
{
  
  int n_row = 12, f_row = 10 ;

  printf("<a href='genehunt.cgi?db=%s'>New search</a>\n", db);
  printf("<br><a href='HelpQuery.html#genehunt_search'>Help</a>") ; /* ,my_db_desc->db_version); */
  
  printf(
	 "<h2>Search Results Found %d Genehunt motifs in %s</h2>\n", n_rows, db) ;
  
  printf(
	 "<form method=get action='genehunt.cgi'>\n"
	 );
  
  printf("<p>");
  printf("<input type=hidden name=flags value='%s'>\n", flags);
  printf("<input type=hidden name=db value='%s'>\n", db);
  
  printf("<input type=submit name=search value='refine'><p>"
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

/**************************************************************/

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

/**************************************************************/

static char *numeric_filter(char *q)
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

/**************************************************************/

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

/**************************************************************/

static char *html_safe (const char *s)
{
  static char *save;
  if (save)
    free(save);
  if (! s || ! *s )
    {
      save = 0;
      return "&nbsp;";
    }
  save = html_encode (s);
  return save;
}

/**************************************************************/

static char *url_safe(const char *s)
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

/**************************************************************/
/**************************************************************/

int main (int argc, char **argv)
{
  char *db ;
  int n, nn ;
  Array tt = 0 ;


  /*
  * cgi debugging hack - you can pass in the query string as argv[1]
  * and the CGI system will believe it is in a real web server
  */
  if ( argv[1] && ( getenv("REQUEST_METHOD") == NULL ) )
    {
    char *s;
    s = argv[1];
    if (*s == '?')
      s++;
    s=messprintf("QUERY_STRING=%s",s);
    s=strdup(s);		/* putenv needs the string forever */
    putenv(strdup(s));
    putenv("REQUEST_METHOD=GET");
    }

   n = cgi_init ("genehunt_cgi.log");
  if (n < 0)
    {
      printf("Content-type: text/html\n\n");
      printf("<h1>cgi_init -> %d\n</h1>", n);
      exit(1);
    }

  /* get or initialise the flags */
  db = cgi_clean_arg ("db") ;
  if (!db) db = "human" ;

  if (0) printf("Content-type: text/html\n\n");

  flagSet () ;
  page_top (db); /* calls selectShow() */

  nn = tableLoad () ;

  if (0)
    {
      char *cp = getenv("QUERY_STRING") ;
      char *cq =  cgi_clean_arg ("nMrna") ;
      
      printf ("param:  %s<br>\n", cp ? cp : "empty") ;
      printf ("nMrna ->  %s<br>\n", cq ? cq : "empty") ;
    }
  printf ("<h3>Loaded %d genes</h3><br>\n", nn) ;
  if (nn) tableShow (db) ;
  page_end () ;
  
  arrayDestroy (tt) ;
  return 0;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/

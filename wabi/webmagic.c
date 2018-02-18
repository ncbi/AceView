/*  File: webmagic.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              linked in with graphical acembly programs
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 10:39 1998 (fw)
 * Created: Mon Apr 18 17:41:24 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: webmagic.c,v 1.8 2017/03/18 15:31:02 mieg Exp $ */

/*
#define ARRAY_CHECK  
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "acembly.h"
#include "freeout.h"
#include "parse.h"
#include "a.h"
#include "bs.h"

#include "display.h"
#include "../wh/menu.h"
#include "plot.h"
#include "key.h"
#include "sysclass.h"
#include "session.h"
#include "pick.h"
#include "cdna.h"
#include "map.h"		/* opaque decl of LOOK */
#include "basecall.h"		/* completes LookStruct -> TRACELOOK */
#include "fmap.h"		/* public fmap function headers */



typedef struct baStuct { 
  AC_HANDLE h ;
  AC_DB db ; 
  AC_TABLE tbl ;
  const char *outFileName ;
  DICT *dict ;
} BA ;

/*************************************************************************************/
/*************************************************************************************/

static void magicHtmlHeader (vTXT txt)
{
  if (1)
    {  /* header */
      vtxtPrintf (txt, 
	       "<html>\n"
	       "  <head>\n"
	       "    <!--Load the AJAX API-->\n"
	       "    <script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n"
		  /*	       "    <script type=\"text/javascript\" src=\"jquery-1.6.2.min.js\"></script>\n" */
	       "    <script type=\"text/javascript\">\n"
	       "    \n"
	       "    // Load the Visualization API and the piechart package.\n"
	       "    google.load('visualization', '1.0', {'packages':['controls']});\n"
	       "    google.load('visualization', '1', {'packages':['corechart']});\n"
	       "    google.load('visualization', '1', {'packages':['table']});\n"
	       "      \n"
	       "    // Set a callback to run when the Google Visualization API is loaded.\n"
	       "    google.setOnLoadCallback(drawDashBoard);\n"
	       "      \n"
	       "    function drawDashBoard() {\n"
	       "              \n"
	       ) ;
    }
} /* magicHtmlHeader */
  
/*************************************************************************************/
/*************************************************************************************/

static void magicTable2jason (vTXT txt, AC_TABLE tbl, const char *types, const char **titres, char *nam)
{
  int ic, ir ;
  const char **cpp ;
  const char *tt ;
  char *typ[256] ;

  typ['f'] = "number" ;
  typ['i'] = "number" ;
  typ['s'] = "string" ;

  if (1) 
    {  /* table */
       vtxtPrintf (txt, 
		"var %s = {\n"
		"  \"cols\": [\n"
		   , nam
		) ;

       for (cpp = titres, tt = types, ic = 0 ; *tt && *cpp ; cpp++, tt++, ic++)
	 vtxtPrintf (txt, "        %s{\"id\":\"%c\",\"label\":\"%s\",\"type\":\"%s\"}\n"
		     , ic ? "," : ""
		     , 'A' + ic
		     , *cpp
		     , typ[(int)*tt]
		     ) ;
	   

       vtxtPrintf (txt, 
		"],\n"
		"\"rows\": ["
		) ;

       for (ir = 0 ; ir < tbl->rows ; tt++, ir++)
	 { 
	   vtxtPrint (txt, "   {\"c\":[") ;

	   for (tt = types, ic = 0 ; ic < tbl->cols ; ic++, tt++)
	     {
	       if (ic)  vtxtPrint (txt, ",") ;
	       switch (*tt)
		 {
		 case 'i':
		   vtxtPrintf (txt, "{\"v\":%d}", ac_table_int (tbl, ir, ic, 0)) ;
		   break ;
		 case 'f':
		   vtxtPrintf (txt, "{\"v\":%g}", ac_table_float (tbl, ir, ic, 0)) ;
		   break ;
		 case 's':
		   vtxtPrintf (txt, "{\"v\":\"%s\"}"
			       , ac_table_printable (tbl, ir, ic, "") 
			       ) ;
		   break ;
		 default:
		   messcrash ("wrong type %c in  magicTable2jason", *tt) ;
		 }
	     }

	   vtxtPrintf (txt, "]}") ;
	   if (ir < tbl->rows - 1)
	     vtxtPrintf (txt, ",") ;
	   vtxtPrintf (txt, "\n") ;
	 }
       
      vtxtPrintf (txt, 
		"]\n"
		"} ;\n"
		) ;
    }
} /* magicTable2jason */

/*


      vtxtPrintf (txt, 
		  "      // Create our data table out of JSON data loaded from server.\n"
		  "      var %s = new google.visualization.DataTable(jsonData);\n"
		  "\n"
		  , nam
		  ) ;

      if (0)
	vtxtPrintf (txt, 
		    "      // Instantiate and draw our chart, passing in some options.\n"
		    "      var chart = new google.visualization.PieChart(document.getElementById('chart_div'));\n"
		    , nam
		    ) ;
 
      if (dashboard)
	vtxtPrintf (txt, 
		    "      // Instantiate and draw our dashboard, passing in some options.\n"
		    "      var chart = new google.visualization.Dashboard(document.getElementById('%s'));\n"
		    , dashboard
		    ) ;
      if (line_graph)
	vtxtPrintf (txt, 
		    "      var chart3 = new google.visualization.ChartWrapper ({\n"
		    "         'chartType' : 'LineChart',\n"
		    "         'containerId' : '%s',\n"
		    "         'options' : {\n"
		    "            'width': 300,\n"
		    "            'height': 300\n"
		    "            }\n"
		    "         }) ;\n"
		    , line_graph
		    , nam
		    ) ;
      if (table_graph)
	vtxtPrintf (txt, 
		    "      var chart3 = new google.visualization.ChartWrapper ({\n"
		    "         'chartType' : 'Table',\n"
		    "         'containerId' : '%s',\n"
		    "         'options' : {\n"
		    "            'width': 300,\n"
		    "            'height': 300\n"
		    "            }\n"
		    "         }) ;\n"
		    , table_graph
		    , nam
		    ) ;
      if (dashboard && line_graph && table_graph)
	vtxtPrintf (txt, 
		    "      %s.bind (%s,%s);\n"
		    , dashboard, line_graph, table_graph
		    ) ; 
      if (dashboard)
	vtxtPrintf (txt, 
		    "      %s.draw(%s);\n"
		    , dashboard, nam
		    ) ;

  
	vtxtPrintf (txt, 
		    "      var chart2 = new google.visualization.ControlWrapperTable({\n"
		    "      chart2.draw(%s, null) ;\n"
		    , dashboard
		    , nam
		    ) ;
*/

 /*************************************************************************************/

static void magicTableDisplay (vTXT txt, char *dashboard, char *line_graph, char *table_graph)
{
  AC_HANDLE h = ac_new_handle () ;

  vtxtPrintf (txt, 
	      "    <div id=\"%s\">\n"
	      , dashboard
	      ) ;
  if (1) 
    {  /* table */
      if (line_graph)
	vtxtPrintf (txt, 
		    "    <div id=\"%s\"></div>\n"
		    , line_graph
		    ) ;
      if (table_graph)
	vtxtPrintf (txt, 
		    "    <div id=\"%s\"></div>\n"
		    , table_graph
		    ) ;
    }

  vtxtPrintf (txt, "</div>\n") ;
  ac_free (h) ;
  return ;
} /* test1 */

/*************************************************************************************/
/*************************************************************************************/
/* display firs last base aligned for one run */

static BOOL magicAli_run (AC_DB db, vTXT txt, const char *run) 
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  const char *errors ;
  const char *titres[] = {"Cycle", "First base aligned F1",  0 };


  
  magicHtmlHeader (txt) ;
  
  qq = hprintf (h, "select x,f1, f2, r1, r2 from r in class \"Ali\" where r==\"%s\", f in r->Letter_profile, xx in f[1], x in r->First_base_aligned where x == xx, y in r->Last_base_aligned where y == x, f1 in x[1], f2 in f1[1], r1 in y[1], r2 in r1[1]", run) ;

  fprintf(stderr, "%s\n",qq) ;
  tbl = ac_bql_table (db, qq, 0,0,&errors,h) ;
  magicTable2jason (txt, tbl, "sffff", titres, "profile1") ;
  
  vtxtPrintf (txt, "\n}\n</script>\n</head>\n\n<body>\n") ;
  
  vtxtPrintf (txt, 
	      "    <h3>First and last base aligned per cycle in run %s</h3>\n"
	      , run
	      ) ;
  
  magicTableDisplay (txt, "div_dashboard", "div_line1", "div_table1") ;

  vtxtPrint (txt, "  </body>\n") ;
  
  ac_free (h) ;
  return TRUE ;
} /* magicATCG_run */

/*
select x,f1, f2, r1, r2 from r in class "Ali" where r="Ghs434", f in r->Letter_profile, xx in f[1], x in r->First_base_aligned where x = xx, y in r->Last_base_aligned where y = x, f1 in x[1], f2 in f1[1], r1 in y[1], r2 in r1[1]

select r,xx, x from r in class "Ali" where r="Ghs434", f in r->Letter_profile, xx in f[1], x in r->First_base_aligned where x = xx

select r,xx,x from r in class "Ali" where r="Ghs434", f in r->Letter_profile, xx in f[1], x in r->First_base_aligned 
*/


/*************************************************************************************/
/*************************************************************************************/
/* display the profile for one run */

static BOOL magicATCG_run (AC_DB db, vTXT txt, const char *run) 
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  const char *errors ;
  const char *titres[] = {"x","Cycle", "A", "T", "G", "C", "N", "quality", 0 };    
  
  magicHtmlHeader (txt) ;
  
  qq = hprintf (h, "select x,x,a,t,g,c,n,q from r in class \"Ali\" where r==\"%s\", f in r->Letter_profile, x in f[1], a in x[1], t in a[1], g in t[1], c in g[1], n in c[1], f1 in r->Quality_profile where f1 == f, x1 in f1[1] where x1==x, q in x1[1]", run) ;
  fprintf(stderr, "%s\n",qq) ;
  tbl = ac_bql_table (db, qq, 0,0,&errors,h) ;
  magicTable2jason (txt, tbl, "siffffff", titres, "ATGC") ;

  
  vtxtPrintf (txt, 
	     "        // Create a dashboard.\n"
"        var dashboard = new google.visualization.Dashboard(\n"
"          document.getElementById('dashboard_div'));\n"
"\n"
"        // Create a range slider, passing some options\n"
"        var cycleSlider = new google.visualization.ControlWrapper({\n"
"          'controlType': 'NumberRangeFilter',\n"
"          'containerId': 'filter_div',\n"
"          'options': {\n"
"            'filterColumnLabel': 'Cycle',\n"
"            'ui': {'labelStacking' : 'vertical'}\n"
"          }\n"
"        });\n"
"\n"
"\n"
"        var chart3 = new google.visualization.ChartWrapper ({\n"
"         'chartType' : 'LineChart',\n"
"         'containerId' : 'div_line1',\n"
"         'options' : {\n"
"            'filterColumnLabel': 'Cycle',\n"
"            'width': 600,\n"
"            'height': 300\n"
"            }\n"
"    , 'view': {'columns': [0,2, 3,4,5,6,7]}\n"
	      /*, 'colors': ['black','green','red','black',blue','gray','cerise']\n" */
"         }) ;\n"
"         var chart4 = new google.visualization.ChartWrapper ({\n"
"         'chartType' : 'Table',\n"
"         'containerId' : 'div_tbl1',\n"
"\n"
"         'options' : {\n"
"            'filterColumnLabel': 'Cycle',\n"
"            'width': '300px',\n"
"            'height': '200px'\n"
"            }\n"
"    , 'view': {'columns': [1,2, 3,4,5,6,7]}\n"
"         }) ;\n"
" \n"
" var barChart = new google.visualization.ChartWrapper({\n"
"    'chartType': 'BarChart',\n"
"    'containerId': 'div_bar1',\n"
"    'options': {\n"
"      'filterColumnLabel': 'Cycle',\n"
"      'width': 400,\n"
"      'height': 300,\n"
"      'chartArea': {top: 0, right: 0, bottom: 0}\n"
"    },\n"
"    // Configure the barchart to use columns 2 (City) and 3 (Population)\n"
"    'view': {'columns': [2, 3]}\n"
"  });\n"
"\n"
"        // Establish dependencies, declaring that 'filter' drives 'pieChart',\n"
"        // so that the pie chart will only display entries that are let through\n"
"        // given the chosen slider range.\n"
"        dashboard.bind([cycleSlider], [chart3,chart4]);\n"
"\n"
"        // Draw the dashboard.\n"
"        dashboard.draw(ATGC);\n"
	      /* , {'colors': ['black','green','red','black',blue','gray','cerise']} */
"      }\n"
"    </script>\n"
"  </head>\n"
"\n"
"  <body>\n"
"<h3>Run %s</h3>\n"
"    <!--Div that will hold the dashboard-->\n"
"    <div id='dashboard_div'>\n"
"      <!--Divs that will hold each control and chart-->\n"
"      <div id='filter_div'></div>\n"
"      <div id='chart_div'></div>\n"

"<table>\n"
"  <tr><td>\n"
"        <div id='div_line1'></div>\n"
"      </td>\n"
"      <td>\n"
"      <div id='div_tbl1'></div>\n"
"</td>\n"
"</tr>\n"
"</table>\n"
"      <div id='div_bar1'></div>\n"
"    </div>\n"
	  , run  ) ;


  vtxtPrint (txt, "  </body>\n") ;
  
  ac_free (h) ;
  return TRUE ;
} /* magicATCG_run */

/*************************************************************************************/
/*************************************************************************************/

static BOOL magicATCG_list (AC_DB db, vTXT txt)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  const char *errors ;
  int ir ;
  const char *titres[] = {"Run", "Number of reads", "kb", "Read length", "Show profile", "Show error types", 0 } ;
  
  magicHtmlHeader (txt) ;
  
  qq = hprintf (h, "select a,t,kb,ln from r in class \"Run\" where exists_tag r->Is_run, a in r->Ali, t in a->Accepted[3], kb in t[2], ln in kb[2], L in a->Letter_profile, E in a->Error_profile") ;
  tbl = ac_bql_table (db, qq, 0,0,&errors,h) ;
  
  if (0)
    for (ir = 0 ; ir < tbl->rows ; ir++)
      vtxtPrintf (txt,"<a href='http:%s'>run %s</a></br></body>"
		  , ac_table_printable (tbl, ir, 0, "")
		  , ac_table_printable (tbl, ir, 0, "")
		  ) ;
  
  magicTable2jason (txt, tbl, "sfff", titres, "runinfo") ; 
  
  vtxtPrintf (txt, "\n}\n</script>\n</head>\n\n<body>\n") ;
  
  vtxtPrintf (txt,"<html><body>\n") ;
  
  vtxtPrintf (txt, 
	      "    <h3>Run information</h3>\n"
	      ) ;
  magicTableDisplay (txt,  "div_dashboard", "div_line0", "div_table0") ;
  vtxtPrint (txt, "  </body>\n") ;
  
  vtxtPrintf (txt,"</body>") ;

  ac_free (h) ;
  return TRUE ;
} /* magicATCG_list */

/*************************************************************************************/
/*************************************************************************************/

static BOOL magicProcessParseRequest (char *command, char *qq[])
{
  char *cp, *cq, *q1 = 0, *q2 = 0 ;

  cp = strstr (command, "magic/") ;
  if (cp)
    {
      cp += 6 ;
      cq = strstr (cp, " ") ;
      if (cq) *cq = 0 ;
      cq = strstr (cp, "/") ;
      if (cq) *cq = 0 ;
      q1 = cp ;
      if (cq)
	{
	  cp = cq+1 ;
	  cq = strstr (cp, " ") ;
	  if (cq) *cq = 0 ;
	  cq = strstr (cp, "/") ;
	  if (cq) *cq = 0 ;
	  q2 = cp ;
	} 
      else
	{
	  cq = strstr (cp, " ") ;
	  if (cq) *cq = 0 ;
	  cq = strstr (cp, "/") ;
	  if (cq) *cq = 0 ;
	}
    }
  qq[1] = q1 ;
  qq[2] = q2 ;

  return TRUE ;
} /* magicProcessParseRequest */

/*************************************************************************************/

static BOOL  magicProcess (vTXT txt, char **qq, AC_HANDLE h)
{
  static AC_DB db = 0 ;
  const char *errors ;

  if (! db)
    db = ac_open_db ("local", &errors) ;
  if (!db) 
    messcrash ("cannot open the local database, no idea why") ;

  /* actual work */
  else if (! strcmp (qq[1], "ATGC"))
    {
      if (qq[2])
	magicATCG_run (db, txt, qq[2]) ;
      else
	magicATCG_list (db, txt) ;
    }
  else if (! strcmp (qq[1], "Ali"))
    {
      if (qq[2])
	magicAli_run (db, txt, qq[2]) ;
      /*      else
	      magicAli_list (db, txt) ;
      */
    } 
    
  vtxtPrintf (txt,"</html>\n") ;
  return TRUE ;
}

/*************************************************************************************/

BOOL magicProcessWebRequest (vTXT txt, char *command)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = TRUE ;
  char *qq[3] ;

  memset (qq, 0, sizeof (qq)) ;
  magicProcessParseRequest (command, qq) ;
  if (! qq[1])
    vtxtPrintf (txt,"<html><body>Bad request : %s <br></body>", command) ;
  else
    magicProcess (txt, qq, h) ;
  /*
bql   select  x,a,t,g,c,n,e from r in class "Ali" where r=="Ghs434" , f in r->Letter_profile, x in f[1], a in x[1], t in a[1], g in t[1], c in g[1], n in c[1], y in r->Error_position where y = x, e in y[1]
  */

  ac_free (h) ;
  return ok ;
}

/*************************************************************************************/
/* [-o file] type run :: export an html page using google API */
KEYSET magicProcessAcedbRequest (KEYSET ksOld)
{
  AC_HANDLE h = ac_new_handle () ;
  int n = 0 ;
  const char *ccp ;
  char *qq[3] ;
  KEYSET ks = 0 ;
  ACEOUT ao = 0 ;
  vTXT txt = 0 ;

  memset (qq, 0, sizeof (qq)) ;

  while (n < 2 && (ccp = freeword ()))
    {
      if (!strcmp (ccp, "-o"))
	{
	  ccp = freeword () ;
	  if (ccp)
	    {
	      ao = aceOutCreate (ccp, 0, 0, h) ;
	      if (! ao)
		{
		  messout ("Sorry, cannot open output file %s",  ccp) ;
		  goto done ;
		}
	    }
	  else
	    {
	      messout ("missing file name after -o, usage : magic [-o file_name] type run") ;
	      goto done ;
	    }
	  continue ;
	}
      qq[++n] = strnew (ccp, h) ;
    }
  if (! ao)
    ao = aceOutCreate (0, 0, 0, h) ;
  if (n)
    {
      txt = vtxtHandleCreate (h) ;
      magicProcess (txt, qq, h) ;
      if (vtxtPtr (txt))
	aceOut (ao, vtxtPtr (txt)) ; 
    }
  else
    messout ("Empty request, usage : magic [-o file_name] type run") ;

 done:
  ac_free (h) ;
  return ks ? ks : keySetCopy (ksOld) ;

} /* magicProcessAcedbRequest */
/* Usage
magic ATGC Rhs531
magic -o /home/mieg/MAC/Desktop/toto.html Ali  Rhs531
magic -o /home/mieg/MAC/Desktop/toto.html ATGC Rhs531
bql select x,a,t,g,c,n from r in class "Ali" where r="Ghs531", f in r->Letter_profile, x in f[1], a in x[1], t in a[1], g in t[1], c in g[1], n in c[1], y in r->Error_position where y == x, e in y[1]
*/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/





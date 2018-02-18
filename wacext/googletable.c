/*  File: googletable.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
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
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.aceview.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */

#include "ac.h"
#include "acedb.h"
/* #define ARRAY_CHECK */

typedef struct baStuct { 
  AC_HANDLE h ;
  AC_DB db ; 
  AC_TABLE tbl ;
  const char *outFileName ;
  DICT *dict ;
} BA ;

/*************************************************************************************/
/*************************************************************************************/

static void test1 (BA *ba, const char **titres)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao  = aceOutCreate (ba->outFileName, ".html", FALSE, h) ;
  int ir ;
  const char **cpp ;

  if (1)
    {  /* header */
      aceOutf (ao, 
	       "<html>\n"
	       "  <head>\n"
	       "    <!--Load the AJAX API-->\n"
	       "    <script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n"
	       "    <script type=\"text/javascript\" src=\"jquery-1.6.2.min.js\"></script>\n"
	       "    <script type=\"text/javascript\">\n"
	       "    \n"
	       "    // Load the Visualization API and the piechart package.\n"
	       "    google.load('visualization', '1', {'packages':['corechart']});\n"
	       "    google.load('visualization', '1', {'packages':['table']});\n"
	       "      \n"
	       "    // Set a callback to run when the Google Visualization API is loaded.\n"
	       "    google.setOnLoadCallback(drawChart);\n"
	       "      \n"
	       "    function drawChart() {\n"
	       "              \n"
	       ) ;
    }
  
  if (1) 
    {  /* table */
       aceOutf (ao, 
		"var jsonData2 = {\n"
		"  \"cols\": [\n"
		"        {\"id\":\"\",\"label\":\"Topping\",\"pattern\":\"\",\"type\":\"string\"},\n"
		"        {\"id\":\"\",\"label\":\"Slices\",\"pattern\":\"\",\"type\":\"number\"},\n"
		"        {\"id\":\"\",\"label\":\"Prices\",\"type\":\"number\"},\n"
		"        {\"id\":\"\",\"label\":\"Weight\",\"type\":\"number\"}\n"
		"],\n"
		"\"rows\": ["
		"   {\"c\":[{\"v\":\"Mushrooms\",\"f\":null},{\"v\":3,\"f\":null},{\"v\":12},{\"v\":43}]},"
		"    {\"c\":[{\"v\":\"Onions\",\"f\":null},{\"v\":1,\"f\":null},{\"v\":15},{\"v\":32}]},"
		"   {\"c\":[{\"v\":\"Olives\",\"f\":null},{\"v\":1,\"f\":null},{\"v\":13},{\"v\":15}]},"
		"   {\"c\":[{\"v\":\"Zucchini\",\"f\":null},{\"v\":1,\"f\":null},{\"v\":40},{\"v\":17}]},"
		"   {\"c\":[{\"v\":\"Pepperoni\",\"f\":null},{\"v\":2,\"f\":null},{\"v\":20},{\"v\":40}]}\n"
		"]\n"
		"} ;\n"
		) ;

       aceOutf (ao, 
		"var jsonData = {\n"
		"  \"cols\": [\n"
		) ;


       cpp = titres ;
	 aceOutf (ao, "        {\"id\":\"\",\"label\":\"%s\",\"pattern\":\"\",\"type\":\"string\"},\n", *cpp) ;
       for (++cpp ; *cpp ; cpp++)
	 aceOutf (ao, "        {\"id\":\"\",\"label\":\"%s\",\"type\":\"number\"},\n", *cpp) ;


       aceOutf (ao, 
		"],\n"
		"\"rows\": ["
		) ;

       for (ir = 0 ; ir < ba->tbl->rows ; ir++)
	 {
	   aceOutf (ao, 
		    "   {\"c\":[{\"v\":\"%s\"},{\"v\":%g},{\"v\":%g}]}"
		    , ac_table_printable (ba->tbl, ir, 0, "") 
		    , ac_table_float (ba->tbl, ir, 1, 0) 
		    , ac_table_float (ba->tbl, ir, 2, 0) 
		    ) ;
	   if (ir < ba->tbl->rows - 1)
	     aceOutf (ao, ",") ;
	   aceOutf (ao, "\n") ;
	 }
       
      aceOutf (ao, 
		"]\n"
		"} ;\n"
		) ;


      
       aceOutf (ao, 
		"      // Create our data table out of JSON data loaded from server.\n"
		"      var data = new google.visualization.DataTable(jsonData);\n"
		"\n"
		"      // Instantiate and draw our chart, passing in some options.\n"
		"      var chart = new google.visualization.PieChart(document.getElementById('chart_div'));\n"
		"      chart.draw(data, {width: 400, height: 240});\n"
		"      var chart2 = new google.visualization.Table(document.getElementById('chart_div2'));\n"
		"      chart2.draw(data, null) ;\n"
		"      var chart3 = new google.visualization.LineChart(document.getElementById('chart_div3'));\n"
		"      chart3.draw(data, null) ;\n"
		"\n"
		"      var atgc = new google.visualization.DataTable() ;\n"
		"      atgc.addColumn ('number','position') ;\n"
		"      atgc.addColumn ('number','A') ;\n"
		"      atgc.addColumn ('number','T') ;\n"
		"      atgc.addColumn ('number','G') ;\n"
		"      atgc.addColumn ('number','C') ;\n"
		"\n"
		"      atgc.addRows ([\n"
		"        [1, 25,25,25,23]\n"
		"        , [2, 25,25,25,23]\n"
		"        , [3, 25,25,25,23]\n"
		"        , [4, 25,25,25,23]\n"
		"        , [5, 20,25,35,23]\n"
		"        , [6, 25,25,25,23]\n"
		"        , [7, 25,12,25,23]\n"
		"        , [8, 25,25,25,23]\n"
		"        , [9, 25,25,25,23]\n"
		"        , [10, 25,25,25,28]\n"
		"      ], false) ;\n"
		"      var chart4 = new google.visualization.LineChart(document.getElementById('chart_div4'));\n"
		"      chart4.draw(atgc, null) ;\n"
		"\n"
		"      }\n"
		"\n"
		"    </script>\n"
		"  </head>\n"
		"\n"
		"  <body>\n"
		"    <!--Div that will hold the pie chart-->\n"
		"    <div id=\"chart_div\"></div>\n"
		"    <div id=\"chart_div2\"></div>\n"
		"    <div id=\"chart_div3\"></div>\n"
		"    <div id=\"chart_div4\"></div>\n"
		"  </body>\n"
		"</html>\n"
		) ;
    }
  ac_free (h) ;
  return ;
} /* test1 */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// googletable.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//   to create tables and other google web pages\n"
	    "// Help\n"
	    "//   -help : export this document\n"
	    "// Example\n"
	    "//    googletable -db MetaDB -project droso -o HTML/droso -export t\n"
	    "// Connect to a server\n"
	    "//   -db <MetaDB> : the name of a MAGIC MetaDB directory\n"
	    "//   -project <name> : the name of a project, typically $MAGIC\n"
	    "// Output\n"
	    "//   -o <file_prefix> : all results will be called file_prefix.*\n"
	    "//   -export [t...] : list of desired tables\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  BA ba ;
  AC_HANDLE h = ac_new_handle () ;
  const char *dbName = 0 ;
  char *errors = 0 ;

  memset (&ba, 0, sizeof (BA)) ;
  ba.h = h ;

  if (argc == 1)
    usage (0) ;

  ba.dict = dictHandleCreate (100, h) ;

  /* optional arguments */

  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  getCmdLineOption (&argc, argv, "-o", &ba.outFileName) ;

  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  if (argc > 1)
    {
      fprintf (stderr, "Unknown argument %s, try -help", argv[argc-1]) ;
      exit (1) ;
    }

  /* actions */

  fprintf (stderr, "Start %s\n", timeShowNow ()) ;

  if (dbName)
    {
      const char *errors ;
      ba.db = ac_open_db (dbName, &errors);
      if (! ba.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }

  /* actual work */
  if (ba.db)
    {
      if (0)
	{
	  const char *titres[] = {"Run", "Reads", "kb", 0 };     
	  ba.tbl = ac_bql_table (ba.db, "select r,t,kb from r in class \"Ali\", s in r->accepted, t in s[2], kb in t[2]", 0,0,&errors,h) ;
	  test1 (&ba, titres) ;
	}
      if (1)
	{
	  const char *titres[] = {"Cycle", "A", "T", "G", "C", "N", 0 };     
	  ba.tbl = ac_bql_table (ba.db, "select x,a,t,g,c,n from r in class \"Ali\" where r==\"exact_ReverseStrand\", f in r->Letter_profile, x in f[1], a in x[7], t in a[1], g in t[1], c in g[1], n in c[1]", 0,0,&errors,h) ;
	  test1 (&ba, titres) ;
	}
    }


  /* clean up */
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }
  
  /* report the global counters */

  ac_free (ba.db) ;
  fprintf (stderr, "Done %s\n", timeShowNow ()) ;
  ac_free (h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


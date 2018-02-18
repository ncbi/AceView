#include "vcgi.h"

/****************************************************************************/
/****************************************************************************/
/* export an H1 title */
void vcgiH1 (vTXT blkp, char *text)
{
  if (text && *text)
    vtxtPrintf (blkp, "<h1>%s</h1>\n", text) ;
} /* vcgiH1 */

/****************************************************************************/
/* start a new html document it with vcgiHtmlEnd () */
void vcgiHtmlStart (vTXT blkp)
{ 
  vtxtPrint (blkp,  "Content-type: text/html\n\n") ; /* the 2 \n are crucial */
} /* vcgiHtmlStart */

void vcgiHtmlEnd (vTXT blkp)
{
  vtxtPrint (blkp, "</HTML>\n") ;
} /* vcgiHtmlEnd */

/****************************************************************************/
/* start a new form, close it with vcgiFormEnd () */
void vcgiFormStart (vTXT blkp, char *action)
{

  vtxtPrint (blkp, "<!-- Form generated using  AceView-cgi library -->");
  vtxtPrint (blkp, "<form method=\"POST\" enctype=\"multipart/form-data\" ");

  vtxtPrintf (blkp
	      , " action=\"%s\">\n"
	      , action
	      );
} /* vcgiFormStart */

/****************************************************************************/

void vcgiFormEnd (vTXT blkp)
{
  vtxtPrint (blkp, "</form>\n" );
} /* vcgiFormEnd */

/****************************************************************************/
/* create a text entry box of length len and initial content dflt */
void vcgiTextEntry (vTXT blkp, char *nam, int len, char *dflt)
{
  vtxtPrintf (blkp
	      , "<input type=\"text\" name=\"%s\">%s\n"
	      , nam
	      , dflt ? dflt : ""
	      ) ;
} /* vcgiTextEntry */

/****************************************************************************/
/* create a multi line text area entry box of initial content dflt */
void vcgiTextAreaEntry (vTXT blkp, char *nam, int width, int height, char *dflt)
{
  if (!nam || !*nam || width < 1 || height <1 ||
      width > 200 || height > 200)
    messcrash ("Bad call of function vcgiTextAreaEntry") ;
  vtxtPrintf (blkp
	      ,  "<textarea NAME=\"%s\" ROWS=%d COLS=%d>\n%s\n</textarea>\n"
	      , nam
	      , height
	      , width
	      , dflt ? dflt : ""
	      ) ;
} /* vcgiTextAreaEntry */

/****************************************************************************/
/* create a integer entry box of initial content dflt */
void vcgiIntEntry (vTXT blkp, char *nam, int dflt)
{
  vtxtPrintf (blkp
	      ,  "<input type=\"text\" name=\"%s\" value=\"%d\">\n"
	      , nam
	      , dflt
	      ) ;
} /* vcgiIntEntry */

/****************************************************************************/
/* create a integer entry box of initial content dflt */
void vcgiFloatEntry (vTXT blkp, char *nam, float dflt)
{
  vtxtPrintf (blkp
	      ,  "<input type=\"text\" name=\"%s\" value=\"%g\">\n"
	      , nam
	      , dflt
	      ) ;
}  /* vcgiFloatEntry */

/****************************************************************************/
/* create a single check box */
void vcgiCheckBox (vTXT blkp, char *nam, BOOL dflt)
{
  vtxtPrintf (blkp
	      , "<input type=\"checkbox\" name=\"%s\" %s>%s\n"
	      , nam
	      , dflt ? "checked" : "unchecked"
	      , nam
	      ) ;
} /* vcgiCheckBox */

/****************************************************************************/
/* create a menu box allowing a single selection */
/* we limit all list to 100, to avoid infinite loops */
void vcgiSelect (vTXT blkp, char *nam, char **options)
{
  char **opt ;
  
  vtxtPrintf (blkp
	      , "<select name=\"%s\">\n"
	      , nam
	      ) ;
  for (opt = options ; *opt ; opt++)
    vtxtPrintf (blkp, "<option value=\"%s\">%s\n"
		, *opt, *opt) ;
  vtxtPrint (blkp, "</select>\n") ;
} /* vcgiSelect */

/****************************************************************************/
/* create a menu box allowing a multiple selection */
void vcgiMultiSelect (vTXT blkp, char *nam, char **options)
{
  char **opt ;
  int i ;
   
  vtxtPrintf (blkp
	      , "<select name=\"%s\" multiple>\n"
	      , nam
	      ) ;
  for (i = 0, opt = options ; i < 100 && *opt ; i++, opt++)
    vtxtPrintf (blkp, "<option value=\"%s\">%s\n"
		, *opt, *opt) ;
  vtxtPrint (blkp, "</select>\n") ;
} /* vcgiMultiSelect */

/****************************************************************************/
/* create a set of mutually exclusive radio buttons */
void vcgiRadioButtons (vTXT blkp, char *nam, char **options)
{
  char **opt ;
  int i ;

  for (i = 0, opt = options ; i < 100 && *opt ; i++, opt++)
    vtxtPrintf (blkp
		, "<input type=\"radio\" name=\"%s\" value=\"%s\">%s\n"
		, nam, *opt, *opt) ;
} /* vcgiRadioButtons */

/****************************************************************************/
/* create a set of non exclusive buttons */
void vcgiMultiButtons (vTXT blkp, char *nam, char **options)
{
  char **opt ;
  int i ;
  
  for (i = 0, opt = options ; i < 100 && *opt ; i++, opt++)
    vtxtPrintf (blkp
		, "<input type=\"checkbox\" name=\"%s\" value=\"%s\">%s\n"
		, nam, *opt, *opt) ;
} /* vcgiMultiButtons */

/****************************************************************************/
/****************************************************************************/


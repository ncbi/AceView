#include "wac/ac.h"
#include "vcgi.h"

static void vcgiFormTest (vTXT blkp, char *action)
{ 
  vtxtMarkup (blkp) ; /* from now on, use html markup in following commands */
  vcgiFormStart (blkp, action);
  vtxtBreak (blkp) ;

  vtxtPrint (blkp, "Text Field containing Plaintext") ;
  vtxtBreak (blkp) ;
  vcgiTextEntry (blkp, "name", 12, "Your name") ;
  vtxtEmptyLine (blkp, 1) ;

  vtxtPrint (blkp, "Multiple-Line Text Field") ;
  vtxtBreak (blkp) ;
  vcgiTextAreaEntry (blkp, "address", 40, 4, "Default contents go here.") ;
  vtxtEmptyLine (blkp, 1) ;
  
  vtxtPrint (blkp, "Checkbox") ; 
  vtxtBreak (blkp) ;
  vcgiCheckBox (blkp, "Hungry", TRUE) ;
  vtxtEmptyLine (blkp, 1) ;
  
  vtxtPrint (blkp, "Float entry");
  vtxtBreak (blkp) ;
  vcgiFloatEntry (blkp, "temperature", 98.6) ;
  vtxtPrint (blkp, "Blood Temperature (80.0-120.0)");
  vtxtEmptyLine (blkp, 1) ;
  
  vtxtPrint (blkp, "Int entry");
  vtxtBreak (blkp) ;
  vcgiFloatEntry (blkp, "frogs", 1) ;
  vtxtPrint (blkp, "Frogs Eaten");
  vtxtEmptyLine (blkp, 1) ;
  
  vtxtPrint (blkp, "Selector");
  vtxtBreak (blkp) ;
  {
    char *colors [] = { "Bleu", "Blanc", "Rouge", 0 } ;
    vcgiSelect (blkp, "colors", colors) ;
  }
  vtxtEmptyLine (blkp, 1) ;

  vtxtPrint (blkp, "Multi-selector");
  vtxtBreak (blkp) ;
  {
    char *flavors [] = { "up", "down", "strange", "charm", "truth", "beauty", 0 } ;
    vcgiMultiSelect (blkp, "colors", flavors) ;
  }
  vtxtEmptyLine (blkp, 1) ;

  vtxtPrint (blkp, "Radio buttons");
  vtxtBreak (blkp) ;
  vtxtPrint (blkp, "How old is your baby");
  {
    char *ages [] = { "1", "2", "3", 0 } ;
    vcgiRadioButtons (blkp, "age", ages) ;
  }
  vtxtEmptyLine (blkp, 1) ;

  vtxtPrint (blkp, "Multi buttons");
  vtxtBreak (blkp) ;
  vtxtPrint (blkp, "Do you have boys and girls");
  {
    char *genders [] = { "Boy", "Girl", 0 } ;
    vcgiMultiButtons (blkp, "gender", genders) ;
  }
  vtxtEmptyLine (blkp, 1) ;

  vcgiFormEnd (blkp);
}

/****************************************************************************/

int main (int argc, char **argv)
{
  char *action = "toto.cgi" ;
  AC_HANDLE h = ac_new_handle () ;
  vTXT blkp = vtxtHandleCreate (h) ; 

  vcgiHtmlStart (blkp) ;
  vtxtPrint (blkp
	     , "<HEAD>\n"
	     "    <TITLE>cgic test</TITLE>\n  </HEAD>\n"
	     "  <BODY>\n"
	     ) ;

  vcgiH1 (blkp, "this is a H1 title") ;
  vcgiFormTest (blkp, action) ;
  vtxtPrint (blkp, "  </BODY>\n</HTML>\n") ;
  vcgiHtmlEnd (blkp) ;
  printf ("%s", vtxtPtr (blkp)) ;
  ac_free (h) ;

  return 0 ;
}

/****************************************************************************/
/****************************************************************************/

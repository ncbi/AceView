#ifndef vCGI_H_DEF
#define vCGI_H_DEF

#include <vtxt.h>
/* export an H1 title */
void vcgiH1 (vTXT blkp, char *text) ;
/* start a new html document it with vcgiHtmlEnd () */
void vcgiHtmlStart (vTXT blkp) ;
void vcgiHtmlEnd (vTXT blkp) ;
/* start a new form, close it with vcgiFormEnd () */
void vcgiFormStart (vTXT blkp, char *actionScrpt) ;
void vcgiFormEnd (vTXT blkp) ;
/* create a text entry box of length len and initial content dflt */
void vcgiTextEntry (vTXT blkp, char *nam, int len, char *dflt) ;
/* create a multi line text area entry box of initial content dflt */
void vcgiTextAreaEntry (vTXT blkp, char *nam, int width, int height, char *dflt) ;
/* create a integer entry box of initial content dflt */
void vcgiIntEntry (vTXT blkp, char *nam, int dflt) ;
/* create a integer entry box of initial content dflt */
void vcgiFloatEntry (vTXT blkp, char *nam, float dflt) ;
/* create a single check box */
void vcgiCheckBox (vTXT blkp, char *nam, BOOL dflt) ;
/* create a menu box allowing a single selection */
/* options must be NULL terminated */
/* we limit all list to 100, to avoid infinite loops in case the NULL is missing */
void vcgiSelect (vTXT blkp, char *nam, char **options) ;
/* create a menu box allowing a multiple selection */
void vcgiMultiSelect (vTXT blkp, char *nam, char **options) ;
/* create a set of mutually exclusive radio buttons */
void vcgiRadioButtons (vTXT blkp, char *nam, char **options) ;
/* create a set of non exclusive buttons */
void vcgiMultiButtons (vTXT blkp, char *nam, char **options) ;

#endif

/*  File: html.h
 *  Authors: Martin Senger (m.senger@dkfz-heidelberg.de)
 *           Petr Kocab    (p.kocab@dkfz-heidelberg.de)
 *
 *-------------------------------------------------------------------
 * Description:
 *    Header file for generating HTML documents.
 *
 * HISTORY:
 *    Created: Sun Apr 16 1995 (senger)
 *-------------------------------------------------------------------
 */

/* $Id: html.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

/******************************************************************************
*
*   --- constants, literals --- 
*
******************************************************************************/
/* --- types of HTML lists --- */
#define HTML_LIST_BULLET       1
#define HTML_LIST_NUMBER       2
#define HTML_LIST_DICT         3

/******************************************************************************
*
*   --- function prototypes --- 
*
******************************************************************************/
FILE* set_html_output (FILE *output);

void html_title (int underline, char* header, char* title);
void html_end (void);
void html (char* format, ...);
void html_button (char* value, char* url);
void html_start_list (int type);
void html_list (char* format, ...);
void html_end_list (void);
void html_line (void);

char* html_bold (char* text);
char* html_link (char* text, char* url);
char* html_name (char* name);
char* html_esc  (char* str);

/* graphwin32_.cpp : implementation file
 *  Author: Richard Bruskiewich
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 *
 * (A superset of) ACEDB GraphTypes are encapsulated within WIN32 Framework C++ objects
 * in such a manner that each possible graph type has a CMultiDocTemplate defined for
 * GraphType-specific IDR_xxxx, CDocument and CView (or similar) derived classes.
 * The specifics for each of these framework classes are defined in separate class definition and
 * implementation files (suitably named) and in the WinAce.rc resource file.
 *
 * A listing of these GraphTypes and Classes is as follows,in enum order of declaration
 * with ???? standing for the source filenames where '*' == 'Doc' and 'View', and the
 * CMultiDocTemplate designated CDocument file extension is .ac? where ? is a unique character 
 * specific to the given GraphType;  the file extension is only used internally to the
 * program, to communicate the GraphType to the CMultiDocTemplate:
 *
 *	GraphType				IDR_xxxx						Document/
 *	(????*.cpp)				(.ac?)							ViewClass(parent view class)
 *  ______________________________________________________________________________________
 *
 *  PLAIN					IDR_PLAIN_GRAPH_TYPE			CPlainGraphDoc/
 *	(Plain*.cpp)			(.acn)							CPlainGraphView(CView)
 *
 *	TEXT_SCROLL				IDR_TEXT_SCROLL_GRAPH_TYPE		CTextScrollGraphDoc/
 *	(TextScroll*.cpp) 		(.acs)							CTextScrollGraphView(CScrollView)
 *
 *	TEXT_FIT				IDR_TEXT_FIT_GRAPH_TYPE			CTextFitGraphDoc/
 *	(TextFit*.cpp)			(.acf)							CTextFitGraphView(CView)
 *
 *	MAP_SCROLL				IDR_MAP_SCROLL_GRAPH_TYPE		CMapScrollGraphDoc/
 *	(MapScroll*.cpp)		(.acm)							CMapScrollGraphView(CScrollView)
 *
 *	PIXEL_SCROLL			IDR_PIXEL_SCROLL_GRAPH_TYPE		CPixelScrollGraphDoc/
 *	(PixelScroll*.cpp)		(.acx)							CPixelScrollGraphView(CScrollView)
 *
 *	TEXT_FULL_SCROLL		IDR_TEXT_FULL_SCROLL_GRAPH_TYPE	CTextFullScrollGraphDoc/
 *	(TextFullScroll*.cpp)	(.acu)							CTextFullScrollGraphView(CScrollView)
 *
 *	PIXEL_FIT				IDR_PIXEL_FIT_GRAPH_TYPE		CPixelFitGraphDoc/
 *	(PixelFit*.cpp)			(.acp)							CPixelFitGraphView(CView)
 *
 *	TEXT_FULL_EDIT			IDR_TEXT_FULL_EDIT_GRAPH_TYPE	CTextFullEditGraphDoc/
 *	(TextFullEdit*.cpp)	(.act)								CTextFullEditGraphView(CEditView)
 *
 *  An additional pseudo-GraphType is the .ace file display, which (contrary to the other
 *  GraphTypes) is externally accessible by the user in the "Open/New" file menus, has the following
 *  pseudo-GraphType specifics:
 *
 *	 	  n/a				IDR_ACE_FILE_TYPE 				CAceFileDoc/
 *	(AceFile*.cpp)			(.ace)							CAceFileView(CEditView)
 *
 * HISTORY:
 * Last edited: Jun 11 03:16 1996 (rbrusk):
 *	-	modified and extended for new WIN32 graph implementation
 * Created: Jul 10 22:50 1995 (rbrusk)
 *-------------------------------------------------------------------*/
/*  Last edited: Jul 10 22:50 1995 (rbrusk) */

// $Id: graphwin32_.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $

/* graphwin32_.cpp
     #include "graph_.h before this file

*/

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "graphwin32_.h" 


GRAPH_DOC_DATA AceGraphDesc[NUMGRAPHTYPES] =
 {
	{IDR_PLAIN_GRAPH_TYPE, ".acn"},
 
 	{IDR_TEXT_SCROLL_GRAPH_TYPE,".acs"},
 
 	{IDR_TEXT_FIT_GRAPH_TYPE, ".acf"},

 	{IDR_MAP_SCROLL_GRAPH_TYPE, ".acm"},
 
 	{IDR_PIXEL_SCROLL_GRAPH_TYPE, ".acx"},
 
 	{IDR_TEXT_FULL_SCROLL_GRAPH_TYPE, ".acu"},
 
 	{IDR_PIXEL_FIT_GRAPH_TYPE, ".acp"},

 	{IDR_TEXT_FULL_EDIT_GRAPH_TYPE, ".act"}
 
} ;
 
  
/****** end of file ******/
 
 

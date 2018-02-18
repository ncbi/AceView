/* $Id: graphloop.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
/* graphloop.h : header file
 *  Author: Richard Bruskiewich
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:   Component of WIN32 emulation of graphLoop()
 * Exported:  class CGraphLoop
 * HISTORY:
 * Last edited:
 * Created: Aug 25 11:56:08 1995 (rmb)
 *-------------------------------------------------------------------
 */

/////////////////////////////////////////////////////////////////////////////
// CGraphLoop thread

class CGraphLoop : public CObject
{
	DECLARE_DYNAMIC(CGraphLoop)

// Attributes
public:
	CGraphLoop();           // protected constructor used by dynamic creation

	// semaphore to flag RECursive graphLoopReturn or graphDestroy invocations
	BOOL recGLR ;

	BOOL isBlocking ;
	int retval ;
	Graph_ theGraph ;			// gActive at time of graphLoop() call
	HWND theGraphWnd ;			// the associated MS Windows HWND handle
	VoidRoutine destroyFunc ; 	// gActive associated DESTROY callback

// Operations
public:
	virtual ~CGraphLoop();

// Overrides

// Implementation
};

/////////////////////////////////////////////////////////////////////////////
 

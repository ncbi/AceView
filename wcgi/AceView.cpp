/*  $Id: AceView.cpp,v 1.4 2009/09/07 02:10:11 mieg Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Denis Vakatov, 
 *          Modifed by Jean Thierry-Mieg: wrapper for an external C cgi 
 *
 * File Description:
 *   CGI/FastCGI application, serving as a wrapper for external C code constructing an html page
 *
 *   USAGE:  sample.cgi?message=Some+Message
 *
 *   NOTE:
 *     1) needs the externally compiled aceviewlib module and libraries
 *
 */

#include <ncbi_pch.hpp>
#include <cgi/cgiapp.hpp>
#include <cgi/cgictx.hpp>

// C_html_page_generator, a C function, exports a complete html page to stdout
extern "C" int C_html_page_generator (void) ;
extern "C" char *C_html_extra_message_generator (int) ;

using namespace ncbi;

/////////////////////////////////////////////////////////////////////////////
//  CCgiSampleApplication::
//

class CCgiSampleApplication : public CCgiApplication
{
public:
  virtual void Init(void);
  virtual int  ProcessRequest(CCgiContext& ctx);
};


void CCgiSampleApplication::Init()
{
  // Standard CGI framework initialization
  CCgiApplication::Init();
}


int CCgiSampleApplication::ProcessRequest(CCgiContext& ctx)
{
  // Given "CGI context", get access to its "HTTP response" sub-object
  CCgiResponse&      response = ctx.GetResponse();
  
  // Call the legacy C code to create the HTML page and flush it
  response.out() <<  C_html_page_generator() ;
  response.Flush();
  
  // Open a diagnostic stream for the web stats
  // 2009_08_07 i remove this because, since introduced on april 17, we have extra errors on web log
  // CDiagContext_Extra extra = GetDiagContext().Extra() ;
  // extra.Print("action", C_html_extra_message_generator(0)) ;
  // extra.Flush();

  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//  MAIN
//

int main(int argc, const char* argv[])
{
    GetDiagContext().SetOldPostFormat(false); // Switch to the new log format
    // 2009_08_11 
    // we used to have
    //    int result = CCgiSampleApplication().AppMain(argc, argv, 0, eDS_Default, "", "AceView");
    //  this seems to create a 'no config file' error, i try to set the last arg to NULL
    //  but if i do that, the applog no longer registers the stats
    int result = CCgiSampleApplication().AppMain(argc, argv, 0, eDS_Default, "", "AceView");
    _TRACE("back to normal diags");
    return result;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


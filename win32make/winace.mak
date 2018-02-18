# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# $Id: winace.mak,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

!IF "$(CFG)" == ""
CFG=WinAce - Win32 Release
!MESSAGE No configuration specified.  Defaulting to WinAce - Win32 Release.
!ENDIF 

!IF "$(CFG)" != "WinAce - Win32 Release" && "$(CFG)" != "WinAce - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "winace.mak" CFG="WinAce - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "WinAce - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "WinAce - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "WinAce - Win32 Debug"
RSC=rc.exe
MTL=mktyplib.exe
CPP=cl.exe

!IF  "$(CFG)" == "WinAce - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "winrel"
# PROP BASE Intermediate_Dir "winrel"
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "c:\acedb\bin\winrel"
# PROP Intermediate_Dir "c:\acedb\bin\winrel"
OUTDIR=c:\acedb\bin\winrel
INTDIR=c:\acedb\bin\winrel

ALL : "$(OUTDIR)\winace.exe" "$(OUTDIR)\winace.bsc"

CLEAN : 
	-@erase "c:\acedb\bin\winrel\winace.bsc"
	-@erase "c:\acedb\bin\winrel\fontpreference.sbr"
	-@erase "c:\acedb\bin\winrel\fmapgene.sbr"
	-@erase "c:\acedb\bin\winrel\freesubs.sbr"
	-@erase "c:\acedb\bin\winrel\chrono.sbr"
	-@erase "c:\acedb\bin\winrel\gmapdata.sbr"
	-@erase "c:\acedb\bin\winrel\fmapfeatures.sbr"
	-@erase "c:\acedb\bin\winrel\gmapmarkercol.sbr"
	-@erase "c:\acedb\bin\winrel\windialogs.sbr"
	-@erase "c:\acedb\bin\winrel\graphwin.sbr"
	-@erase "c:\acedb\bin\winrel\plainview.sbr"
	-@erase "c:\acedb\bin\winrel\sysclass.sbr"
	-@erase "c:\acedb\bin\winrel\parse.sbr"
	-@erase "c:\acedb\bin\winrel\align.sbr"
	-@erase "c:\acedb\bin\winrel\vmapdrag.sbr"
	-@erase "c:\acedb\bin\winrel\sprdmap.sbr"
	-@erase "c:\acedb\bin\winrel\timesubs.sbr"
	-@erase "c:\acedb\bin\winrel\dnacpt.sbr"
	-@erase "c:\acedb\bin\winrel\colourbox.sbr"
	-@erase "c:\acedb\bin\winrel\check.sbr"
	-@erase "c:\acedb\bin\winrel\fpdisp.sbr"
	-@erase "c:\acedb\bin\winrel\preferences.sbr"
	-@erase "c:\acedb\bin\winrel\fmapsequence.sbr"
	-@erase "c:\acedb\bin\winrel\pepdisp.sbr"
	-@erase "c:\acedb\bin\winrel\hexcode.sbr"
	-@erase "c:\acedb\bin\winrel\cmapdisp.sbr"
	-@erase "c:\acedb\bin\winrel\bump.sbr"
	-@erase "c:\acedb\bin\winrel\viewedit.sbr"
	-@erase "c:\acedb\bin\winrel\sprdop.sbr"
	-@erase "c:\acedb\bin\winrel\session.sbr"
	-@erase "c:\acedb\bin\winrel\alignment.sbr"
	-@erase "c:\acedb\bin\winrel\gfcode.sbr"
	-@erase "c:\acedb\bin\winrel\tags.sbr"
	-@erase "c:\acedb\bin\winrel\banner.sbr"
	-@erase "c:\acedb\bin\winrel\cmultiass.sbr"
	-@erase "c:\acedb\bin\winrel\pepseqcol.sbr"
	-@erase "c:\acedb\bin\winrel\win32.sbr"
	-@erase "c:\acedb\bin\winrel\status.sbr"
	-@erase "c:\acedb\bin\winrel\vmapdata2.sbr"
	-@erase "c:\acedb\bin\winrel\mapscrollview.sbr"
	-@erase "c:\acedb\bin\winrel\freeout.sbr"
	-@erase "c:\acedb\bin\winrel\fitview.sbr"
	-@erase "c:\acedb\bin\winrel\dict.sbr"
	-@erase "c:\acedb\bin\winrel\sprddata.sbr"
	-@erase "c:\acedb\bin\winrel\mapcontrol.sbr"
	-@erase "c:\acedb\bin\winrel\textscrollview.sbr"
	-@erase "c:\acedb\bin\winrel\fmapcontrol.sbr"
	-@erase "c:\acedb\bin\winrel\embl.sbr"
	-@erase "c:\acedb\bin\winrel\pepgraphcol.sbr"
	-@erase "c:\acedb\bin\winrel\win32process.sbr"
	-@erase "c:\acedb\bin\winrel\lexalpha.sbr"
	-@erase "c:\acedb\bin\winrel\cgraph.sbr"
	-@erase "c:\acedb\bin\winrel\picksubs.sbr"
	-@erase "c:\acedb\bin\winrel\linkdate.sbr"
	-@erase "c:\acedb\bin\winrel\bsubs.sbr"
	-@erase "c:\acedb\bin\winrel\biblio.sbr"
	-@erase "c:\acedb\bin\winrel\qbedisp.sbr"
	-@erase "c:\acedb\bin\winrel\objcache.sbr"
	-@erase "c:\acedb\bin\winrel\winace.sbr"
	-@erase "c:\acedb\bin\winrel\acedbprofile.sbr"
	-@erase "c:\acedb\bin\winrel\keyset.sbr"
	-@erase "c:\acedb\bin\winrel\pephomolcol.sbr"
	-@erase "c:\acedb\bin\winrel\blxview.sbr"
	-@erase "c:\acedb\bin\winrel\randsubs.sbr"
	-@erase "c:\acedb\bin\winrel\nicedump.sbr"
	-@erase "c:\acedb\bin\winrel\display.sbr"
	-@erase "c:\acedb\bin\winrel\blocksub.sbr"
	-@erase "c:\acedb\bin\winrel\longtext.sbr"
	-@erase "c:\acedb\bin\winrel\sprdctrl.sbr"
	-@erase "c:\acedb\bin\winrel\msgwindow.sbr"
	-@erase "c:\acedb\bin\winrel\afxtrace.sbr"
	-@erase "c:\acedb\bin\winrel\filsubs.sbr"
	-@erase "c:\acedb\bin\winrel\dotterKarlin.sbr"
	-@erase "c:\acedb\bin\winrel\userinterface.sbr"
	-@erase "c:\acedb\bin\winrel\winfilquery.sbr"
	-@erase "c:\acedb\bin\winrel\dotter.sbr"
	-@erase "c:\acedb\bin\winrel\graphcon.sbr"
	-@erase "c:\acedb\bin\winrel\acedbprofileview.sbr"
	-@erase "c:\acedb\bin\winrel\quovadis.sbr"
	-@erase "c:\acedb\bin\winrel\caceprintpage.sbr"
	-@erase "c:\acedb\bin\winrel\splashbox.sbr"
	-@erase "c:\acedb\bin\winrel\gmapconvert.sbr"
	-@erase "c:\acedb\bin\winrel\messubs.sbr"
	-@erase "c:\acedb\bin\winrel\peptide.sbr"
	-@erase "c:\acedb\bin\winrel\bstools.sbr"
	-@erase "c:\acedb\bin\winrel\graphbox.sbr"
	-@erase "c:\acedb\bin\winrel\colcontrol.sbr"
	-@erase "c:\acedb\bin\winrel\mainpick.sbr"
	-@erase "c:\acedb\bin\winrel\action.sbr"
	-@erase "c:\acedb\bin\winrel\pmapdisp.sbr"
	-@erase "c:\acedb\bin\winrel\graphview.sbr"
	-@erase "c:\acedb\bin\winrel\geldisp.sbr"
	-@erase "c:\acedb\bin\winrel\gmapdatacol.sbr"
	-@erase "c:\acedb\bin\winrel\pixelfitview.sbr"
	-@erase "c:\acedb\bin\winrel\pixelscrollview.sbr"
	-@erase "c:\acedb\bin\winrel\class.sbr"
	-@erase "c:\acedb\bin\winrel\ksetdisp.sbr"
	-@erase "c:\acedb\bin\winrel\acefileview.sbr"
	-@erase "c:\acedb\bin\winrel\lexsubs4.sbr"
	-@erase "c:\acedb\bin\winrel\win32menusubs.sbr"
	-@erase "c:\acedb\bin\winrel\dbprofile.sbr"
	-@erase "c:\acedb\bin\winrel\pmapconvert.sbr"
	-@erase "c:\acedb\bin\winrel\method.sbr"
	-@erase "c:\acedb\bin\winrel\model.sbr"
	-@erase "c:\acedb\bin\winrel\update.sbr"
	-@erase "c:\acedb\bin\winrel\treedisp.sbr"
	-@erase "c:\acedb\bin\winrel\keysetdump.sbr"
	-@erase "c:\acedb\bin\winrel\win32lib.sbr"
	-@erase "c:\acedb\bin\winrel\asubs.sbr"
	-@erase "c:\acedb\bin\winrel\opp.sbr"
	-@erase "c:\acedb\bin\winrel\dump.sbr"
	-@erase "c:\acedb\bin\winrel\gmapremarkcol.sbr"
	-@erase "c:\acedb\bin\winrel\plot.sbr"
	-@erase "c:\acedb\bin\winrel\bssubs.sbr"
	-@erase "c:\acedb\bin\winrel\graphramp.sbr"
	-@erase "c:\acedb\bin\winrel\winmain.sbr"
	-@erase "c:\acedb\bin\winrel\fullscrollview.sbr"
	-@erase "c:\acedb\bin\winrel\gd.sbr"
	-@erase "c:\acedb\bin\winrel\bsdumps.sbr"
	-@erase "c:\acedb\bin\winrel\help.sbr"
	-@erase "c:\acedb\bin\winrel\bstree.sbr"
	-@erase "c:\acedb\bin\winrel\graphwin32_.sbr"
	-@erase "c:\acedb\bin\winrel\lexsubs.sbr"
	-@erase "c:\acedb\bin\winrel\queryexe.sbr"
	-@erase "c:\acedb\bin\winrel\vmapdisp.sbr"
	-@erase "c:\acedb\bin\winrel\acefiledoc.sbr"
	-@erase "c:\acedb\bin\winrel\metab.sbr"
	-@erase "c:\acedb\bin\winrel\stdafx.sbr"
	-@erase "c:\acedb\bin\winrel\griddisp.sbr"
	-@erase "c:\acedb\bin\winrel\mainfrm.sbr"
	-@erase "c:\acedb\bin\winrel\drawdisp.sbr"
	-@erase "c:\acedb\bin\winrel\querydisp.sbr"
	-@erase "c:\acedb\bin\winrel\vmapphys.sbr"
	-@erase "c:\acedb\bin\winrel\oldhelp.sbr"
	-@erase "c:\acedb\bin\winrel\graphsub.sbr"
	-@erase "c:\acedb\bin\winrel\heap.sbr"
	-@erase "c:\acedb\bin\winrel\sprddisplay.sbr"
	-@erase "c:\acedb\bin\winrel\gmapposnegcol.sbr"
	-@erase "c:\acedb\bin\winrel\newkey.sbr"
	-@erase "c:\acedb\bin\winrel\menu.sbr"
	-@erase "c:\acedb\bin\winrel\gmapdisp.sbr"
	-@erase "c:\acedb\bin\winrel\gmapintervalcol.sbr"
	-@erase "c:\acedb\bin\winrel\win32print.sbr"
	-@erase "c:\acedb\bin\winrel\translate.sbr"
	-@erase "c:\acedb\bin\winrel\acedialogs.sbr"
	-@erase "c:\acedb\bin\winrel\disknew.sbr"
	-@erase "c:\acedb\bin\winrel\gmapphys.sbr"
	-@erase "c:\acedb\bin\winrel\dnasubs.sbr"
	-@erase "c:\acedb\bin\winrel\textfitview.sbr"
	-@erase "c:\acedb\bin\winrel\gmaplocuscol.sbr"
	-@erase "c:\acedb\bin\winrel\call.sbr"
	-@erase "c:\acedb\bin\winrel\graphselect.sbr"
	-@erase "c:\acedb\bin\winrel\querybuild.sbr"
	-@erase "c:\acedb\bin\winrel\arraysub.sbr"
	-@erase "c:\acedb\bin\winrel\winace.exe"
	-@erase "c:\acedb\bin\winrel\graphwin32_.obj"
	-@erase "c:\acedb\bin\winrel\lexsubs.obj"
	-@erase "c:\acedb\bin\winrel\queryexe.obj"
	-@erase "c:\acedb\bin\winrel\vmapdisp.obj"
	-@erase "c:\acedb\bin\winrel\acefiledoc.obj"
	-@erase "c:\acedb\bin\winrel\metab.obj"
	-@erase "c:\acedb\bin\winrel\stdafx.obj"
	-@erase "c:\acedb\bin\winrel\griddisp.obj"
	-@erase "c:\acedb\bin\winrel\mainfrm.obj"
	-@erase "c:\acedb\bin\winrel\drawdisp.obj"
	-@erase "c:\acedb\bin\winrel\querydisp.obj"
	-@erase "c:\acedb\bin\winrel\vmapphys.obj"
	-@erase "c:\acedb\bin\winrel\oldhelp.obj"
	-@erase "c:\acedb\bin\winrel\graphsub.obj"
	-@erase "c:\acedb\bin\winrel\heap.obj"
	-@erase "c:\acedb\bin\winrel\sprddisplay.obj"
	-@erase "c:\acedb\bin\winrel\gmapposnegcol.obj"
	-@erase "c:\acedb\bin\winrel\newkey.obj"
	-@erase "c:\acedb\bin\winrel\menu.obj"
	-@erase "c:\acedb\bin\winrel\gmapdisp.obj"
	-@erase "c:\acedb\bin\winrel\gmapintervalcol.obj"
	-@erase "c:\acedb\bin\winrel\win32print.obj"
	-@erase "c:\acedb\bin\winrel\translate.obj"
	-@erase "c:\acedb\bin\winrel\acedialogs.obj"
	-@erase "c:\acedb\bin\winrel\disknew.obj"
	-@erase "c:\acedb\bin\winrel\gmapphys.obj"
	-@erase "c:\acedb\bin\winrel\dnasubs.obj"
	-@erase "c:\acedb\bin\winrel\textfitview.obj"
	-@erase "c:\acedb\bin\winrel\gmaplocuscol.obj"
	-@erase "c:\acedb\bin\winrel\call.obj"
	-@erase "c:\acedb\bin\winrel\graphselect.obj"
	-@erase "c:\acedb\bin\winrel\querybuild.obj"
	-@erase "c:\acedb\bin\winrel\arraysub.obj"
	-@erase "c:\acedb\bin\winrel\fontpreference.obj"
	-@erase "c:\acedb\bin\winrel\fmapgene.obj"
	-@erase "c:\acedb\bin\winrel\freesubs.obj"
	-@erase "c:\acedb\bin\winrel\chrono.obj"
	-@erase "c:\acedb\bin\winrel\gmapdata.obj"
	-@erase "c:\acedb\bin\winrel\fmapfeatures.obj"
	-@erase "c:\acedb\bin\winrel\gmapmarkercol.obj"
	-@erase "c:\acedb\bin\winrel\windialogs.obj"
	-@erase "c:\acedb\bin\winrel\graphwin.obj"
	-@erase "c:\acedb\bin\winrel\plainview.obj"
	-@erase "c:\acedb\bin\winrel\sysclass.obj"
	-@erase "c:\acedb\bin\winrel\parse.obj"
	-@erase "c:\acedb\bin\winrel\align.obj"
	-@erase "c:\acedb\bin\winrel\vmapdrag.obj"
	-@erase "c:\acedb\bin\winrel\sprdmap.obj"
	-@erase "c:\acedb\bin\winrel\timesubs.obj"
	-@erase "c:\acedb\bin\winrel\dnacpt.obj"
	-@erase "c:\acedb\bin\winrel\colourbox.obj"
	-@erase "c:\acedb\bin\winrel\check.obj"
	-@erase "c:\acedb\bin\winrel\fpdisp.obj"
	-@erase "c:\acedb\bin\winrel\preferences.obj"
	-@erase "c:\acedb\bin\winrel\fmapsequence.obj"
	-@erase "c:\acedb\bin\winrel\pepdisp.obj"
	-@erase "c:\acedb\bin\winrel\hexcode.obj"
	-@erase "c:\acedb\bin\winrel\cmapdisp.obj"
	-@erase "c:\acedb\bin\winrel\bump.obj"
	-@erase "c:\acedb\bin\winrel\viewedit.obj"
	-@erase "c:\acedb\bin\winrel\sprdop.obj"
	-@erase "c:\acedb\bin\winrel\session.obj"
	-@erase "c:\acedb\bin\winrel\alignment.obj"
	-@erase "c:\acedb\bin\winrel\gfcode.obj"
	-@erase "c:\acedb\bin\winrel\tags.obj"
	-@erase "c:\acedb\bin\winrel\banner.obj"
	-@erase "c:\acedb\bin\winrel\cmultiass.obj"
	-@erase "c:\acedb\bin\winrel\pepseqcol.obj"
	-@erase "c:\acedb\bin\winrel\win32.obj"
	-@erase "c:\acedb\bin\winrel\status.obj"
	-@erase "c:\acedb\bin\winrel\vmapdata2.obj"
	-@erase "c:\acedb\bin\winrel\mapscrollview.obj"
	-@erase "c:\acedb\bin\winrel\freeout.obj"
	-@erase "c:\acedb\bin\winrel\fitview.obj"
	-@erase "c:\acedb\bin\winrel\dict.obj"
	-@erase "c:\acedb\bin\winrel\sprddata.obj"
	-@erase "c:\acedb\bin\winrel\mapcontrol.obj"
	-@erase "c:\acedb\bin\winrel\textscrollview.obj"
	-@erase "c:\acedb\bin\winrel\fmapcontrol.obj"
	-@erase "c:\acedb\bin\winrel\embl.obj"
	-@erase "c:\acedb\bin\winrel\pepgraphcol.obj"
	-@erase "c:\acedb\bin\winrel\win32process.obj"
	-@erase "c:\acedb\bin\winrel\lexalpha.obj"
	-@erase "c:\acedb\bin\winrel\cgraph.obj"
	-@erase "c:\acedb\bin\winrel\picksubs.obj"
	-@erase "c:\acedb\bin\winrel\linkdate.obj"
	-@erase "c:\acedb\bin\winrel\bsubs.obj"
	-@erase "c:\acedb\bin\winrel\biblio.obj"
	-@erase "c:\acedb\bin\winrel\qbedisp.obj"
	-@erase "c:\acedb\bin\winrel\objcache.obj"
	-@erase "c:\acedb\bin\winrel\winace.obj"
	-@erase "c:\acedb\bin\winrel\acedbprofile.obj"
	-@erase "c:\acedb\bin\winrel\keyset.obj"
	-@erase "c:\acedb\bin\winrel\pephomolcol.obj"
	-@erase "c:\acedb\bin\winrel\blxview.obj"
	-@erase "c:\acedb\bin\winrel\randsubs.obj"
	-@erase "c:\acedb\bin\winrel\nicedump.obj"
	-@erase "c:\acedb\bin\winrel\display.obj"
	-@erase "c:\acedb\bin\winrel\blocksub.obj"
	-@erase "c:\acedb\bin\winrel\longtext.obj"
	-@erase "c:\acedb\bin\winrel\sprdctrl.obj"
	-@erase "c:\acedb\bin\winrel\msgwindow.obj"
	-@erase "c:\acedb\bin\winrel\afxtrace.obj"
	-@erase "c:\acedb\bin\winrel\filsubs.obj"
	-@erase "c:\acedb\bin\winrel\dotterKarlin.obj"
	-@erase "c:\acedb\bin\winrel\userinterface.obj"
	-@erase "c:\acedb\bin\winrel\winfilquery.obj"
	-@erase "c:\acedb\bin\winrel\dotter.obj"
	-@erase "c:\acedb\bin\winrel\graphcon.obj"
	-@erase "c:\acedb\bin\winrel\acedbprofileview.obj"
	-@erase "c:\acedb\bin\winrel\quovadis.obj"
	-@erase "c:\acedb\bin\winrel\caceprintpage.obj"
	-@erase "c:\acedb\bin\winrel\splashbox.obj"
	-@erase "c:\acedb\bin\winrel\gmapconvert.obj"
	-@erase "c:\acedb\bin\winrel\messubs.obj"
	-@erase "c:\acedb\bin\winrel\peptide.obj"
	-@erase "c:\acedb\bin\winrel\bstools.obj"
	-@erase "c:\acedb\bin\winrel\graphbox.obj"
	-@erase "c:\acedb\bin\winrel\colcontrol.obj"
	-@erase "c:\acedb\bin\winrel\mainpick.obj"
	-@erase "c:\acedb\bin\winrel\action.obj"
	-@erase "c:\acedb\bin\winrel\pmapdisp.obj"
	-@erase "c:\acedb\bin\winrel\graphview.obj"
	-@erase "c:\acedb\bin\winrel\geldisp.obj"
	-@erase "c:\acedb\bin\winrel\gmapdatacol.obj"
	-@erase "c:\acedb\bin\winrel\pixelfitview.obj"
	-@erase "c:\acedb\bin\winrel\pixelscrollview.obj"
	-@erase "c:\acedb\bin\winrel\class.obj"
	-@erase "c:\acedb\bin\winrel\ksetdisp.obj"
	-@erase "c:\acedb\bin\winrel\acefileview.obj"
	-@erase "c:\acedb\bin\winrel\lexsubs4.obj"
	-@erase "c:\acedb\bin\winrel\win32menusubs.obj"
	-@erase "c:\acedb\bin\winrel\dbprofile.obj"
	-@erase "c:\acedb\bin\winrel\pmapconvert.obj"
	-@erase "c:\acedb\bin\winrel\method.obj"
	-@erase "c:\acedb\bin\winrel\model.obj"
	-@erase "c:\acedb\bin\winrel\update.obj"
	-@erase "c:\acedb\bin\winrel\treedisp.obj"
	-@erase "c:\acedb\bin\winrel\keysetdump.obj"
	-@erase "c:\acedb\bin\winrel\win32lib.obj"
	-@erase "c:\acedb\bin\winrel\asubs.obj"
	-@erase "c:\acedb\bin\winrel\opp.obj"
	-@erase "c:\acedb\bin\winrel\dump.obj"
	-@erase "c:\acedb\bin\winrel\gmapremarkcol.obj"
	-@erase "c:\acedb\bin\winrel\plot.obj"
	-@erase "c:\acedb\bin\winrel\bssubs.obj"
	-@erase "c:\acedb\bin\winrel\graphramp.obj"
	-@erase "c:\acedb\bin\winrel\winmain.obj"
	-@erase "c:\acedb\bin\winrel\fullscrollview.obj"
	-@erase "c:\acedb\bin\winrel\gd.obj"
	-@erase "c:\acedb\bin\winrel\bsdumps.obj"
	-@erase "c:\acedb\bin\winrel\help.obj"
	-@erase "c:\acedb\bin\winrel\bstree.obj"
	-@erase "c:\acedb\bin\winrel\winace.res"
	-@erase "c:\acedb\bin\winrel\winace.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D "ACEDB4" /D "_MBCS" /Fr /YX /c
CPP_PROJ=/nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fp"$(INTDIR)/winace.pch" /YX\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=c:\acedb\bin\winrel/
CPP_SBRS=c:\acedb\bin\winrel/
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/winace.bsc" 
BSC32_SBRS= \
	"$(INTDIR)/fontpreference.sbr" \
	"$(INTDIR)/fmapgene.sbr" \
	"$(INTDIR)/freesubs.sbr" \
	"$(INTDIR)/chrono.sbr" \
	"$(INTDIR)/gmapdata.sbr" \
	"$(INTDIR)/fmapfeatures.sbr" \
	"$(INTDIR)/gmapmarkercol.sbr" \
	"$(INTDIR)/windialogs.sbr" \
	"$(INTDIR)/graphwin.sbr" \
	"$(INTDIR)/plainview.sbr" \
	"$(INTDIR)/sysclass.sbr" \
	"$(INTDIR)/parse.sbr" \
	"$(INTDIR)/align.sbr" \
	"$(INTDIR)/vmapdrag.sbr" \
	"$(INTDIR)/sprdmap.sbr" \
	"$(INTDIR)/timesubs.sbr" \
	"$(INTDIR)/dnacpt.sbr" \
	"$(INTDIR)/colourbox.sbr" \
	"$(INTDIR)/check.sbr" \
	"$(INTDIR)/fpdisp.sbr" \
	"$(INTDIR)/preferences.sbr" \
	"$(INTDIR)/fmapsequence.sbr" \
	"$(INTDIR)/pepdisp.sbr" \
	"$(INTDIR)/hexcode.sbr" \
	"$(INTDIR)/cmapdisp.sbr" \
	"$(INTDIR)/bump.sbr" \
	"$(INTDIR)/viewedit.sbr" \
	"$(INTDIR)/sprdop.sbr" \
	"$(INTDIR)/session.sbr" \
	"$(INTDIR)/alignment.sbr" \
	"$(INTDIR)/gfcode.sbr" \
	"$(INTDIR)/tags.sbr" \
	"$(INTDIR)/banner.sbr" \
	"$(INTDIR)/cmultiass.sbr" \
	"$(INTDIR)/pepseqcol.sbr" \
	"$(INTDIR)/win32.sbr" \
	"$(INTDIR)/status.sbr" \
	"$(INTDIR)/vmapdata2.sbr" \
	"$(INTDIR)/mapscrollview.sbr" \
	"$(INTDIR)/freeout.sbr" \
	"$(INTDIR)/fitview.sbr" \
	"$(INTDIR)/dict.sbr" \
	"$(INTDIR)/sprddata.sbr" \
	"$(INTDIR)/mapcontrol.sbr" \
	"$(INTDIR)/textscrollview.sbr" \
	"$(INTDIR)/fmapcontrol.sbr" \
	"$(INTDIR)/embl.sbr" \
	"$(INTDIR)/pepgraphcol.sbr" \
	"$(INTDIR)/win32process.sbr" \
	"$(INTDIR)/lexalpha.sbr" \
	"$(INTDIR)/cgraph.sbr" \
	"$(INTDIR)/picksubs.sbr" \
	"$(INTDIR)/linkdate.sbr" \
	"$(INTDIR)/bsubs.sbr" \
	"$(INTDIR)/biblio.sbr" \
	"$(INTDIR)/qbedisp.sbr" \
	"$(INTDIR)/objcache.sbr" \
	"$(INTDIR)/winace.sbr" \
	"$(INTDIR)/acedbprofile.sbr" \
	"$(INTDIR)/keyset.sbr" \
	"$(INTDIR)/pephomolcol.sbr" \
	"$(INTDIR)/blxview.sbr" \
	"$(INTDIR)/randsubs.sbr" \
	"$(INTDIR)/nicedump.sbr" \
	"$(INTDIR)/display.sbr" \
	"$(INTDIR)/blocksub.sbr" \
	"$(INTDIR)/longtext.sbr" \
	"$(INTDIR)/sprdctrl.sbr" \
	"$(INTDIR)/msgwindow.sbr" \
	"$(INTDIR)/afxtrace.sbr" \
	"$(INTDIR)/filsubs.sbr" \
	"$(INTDIR)/dotterKarlin.sbr" \
	"$(INTDIR)/userinterface.sbr" \
	"$(INTDIR)/winfilquery.sbr" \
	"$(INTDIR)/dotter.sbr" \
	"$(INTDIR)/graphcon.sbr" \
	"$(INTDIR)/acedbprofileview.sbr" \
	"$(INTDIR)/quovadis.sbr" \
	"$(INTDIR)/caceprintpage.sbr" \
	"$(INTDIR)/splashbox.sbr" \
	"$(INTDIR)/gmapconvert.sbr" \
	"$(INTDIR)/messubs.sbr" \
	"$(INTDIR)/peptide.sbr" \
	"$(INTDIR)/bstools.sbr" \
	"$(INTDIR)/graphbox.sbr" \
	"$(INTDIR)/colcontrol.sbr" \
	"$(INTDIR)/mainpick.sbr" \
	"$(INTDIR)/action.sbr" \
	"$(INTDIR)/pmapdisp.sbr" \
	"$(INTDIR)/graphview.sbr" \
	"$(INTDIR)/geldisp.sbr" \
	"$(INTDIR)/gmapdatacol.sbr" \
	"$(INTDIR)/pixelfitview.sbr" \
	"$(INTDIR)/pixelscrollview.sbr" \
	"$(INTDIR)/class.sbr" \
	"$(INTDIR)/ksetdisp.sbr" \
	"$(INTDIR)/acefileview.sbr" \
	"$(INTDIR)/lexsubs4.sbr" \
	"$(INTDIR)/win32menusubs.sbr" \
	"$(INTDIR)/dbprofile.sbr" \
	"$(INTDIR)/pmapconvert.sbr" \
	"$(INTDIR)/method.sbr" \
	"$(INTDIR)/model.sbr" \
	"$(INTDIR)/update.sbr" \
	"$(INTDIR)/treedisp.sbr" \
	"$(INTDIR)/keysetdump.sbr" \
	"$(INTDIR)/win32lib.sbr" \
	"$(INTDIR)/asubs.sbr" \
	"$(INTDIR)/opp.sbr" \
	"$(INTDIR)/dump.sbr" \
	"$(INTDIR)/gmapremarkcol.sbr" \
	"$(INTDIR)/plot.sbr" \
	"$(INTDIR)/bssubs.sbr" \
	"$(INTDIR)/graphramp.sbr" \
	"$(INTDIR)/winmain.sbr" \
	"$(INTDIR)/fullscrollview.sbr" \
	"$(INTDIR)/gd.sbr" \
	"$(INTDIR)/bsdumps.sbr" \
	"$(INTDIR)/help.sbr" \
	"$(INTDIR)/bstree.sbr" \
	"$(INTDIR)/graphwin32_.sbr" \
	"$(INTDIR)/lexsubs.sbr" \
	"$(INTDIR)/queryexe.sbr" \
	"$(INTDIR)/vmapdisp.sbr" \
	"$(INTDIR)/acefiledoc.sbr" \
	"$(INTDIR)/metab.sbr" \
	"$(INTDIR)/stdafx.sbr" \
	"$(INTDIR)/griddisp.sbr" \
	"$(INTDIR)/mainfrm.sbr" \
	"$(INTDIR)/drawdisp.sbr" \
	"$(INTDIR)/querydisp.sbr" \
	"$(INTDIR)/vmapphys.sbr" \
	"$(INTDIR)/oldhelp.sbr" \
	"$(INTDIR)/graphsub.sbr" \
	"$(INTDIR)/heap.sbr" \
	"$(INTDIR)/sprddisplay.sbr" \
	"$(INTDIR)/gmapposnegcol.sbr" \
	"$(INTDIR)/newkey.sbr" \
	"$(INTDIR)/menu.sbr" \
	"$(INTDIR)/gmapdisp.sbr" \
	"$(INTDIR)/gmapintervalcol.sbr" \
	"$(INTDIR)/win32print.sbr" \
	"$(INTDIR)/translate.sbr" \
	"$(INTDIR)/acedialogs.sbr" \
	"$(INTDIR)/disknew.sbr" \
	"$(INTDIR)/gmapphys.sbr" \
	"$(INTDIR)/dnasubs.sbr" \
	"$(INTDIR)/textfitview.sbr" \
	"$(INTDIR)/gmaplocuscol.sbr" \
	"$(INTDIR)/call.sbr" \
	"$(INTDIR)/graphselect.sbr" \
	"$(INTDIR)/querybuild.sbr" \
	"$(INTDIR)/arraysub.sbr"

"$(OUTDIR)\winace.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /map /machine:I386 /nodefaultlib:"libc.lib" /nodefaultlib:" libcmt.lib" /nodefaultlib:" msvcrtd.lib" /nodefaultlib:" libcd.lib" /nodefaultlib:" libcmtd.lib"
# SUBTRACT LINK32 /verbose /nodefaultlib
LINK32_FLAGS=/nologo /version:4.31 /entry:"WinMainCRTStartup"\
 /subsystem:windows /incremental:no /pdb:"$(OUTDIR)/winace.pdb"\
 /map:"$(INTDIR)/winace.map" /machine:I386 /nodefaultlib:"libc.lib"\
 /nodefaultlib:" libcmt.lib" /nodefaultlib:" msvcrtd.lib"\
 /nodefaultlib:" libcd.lib" /nodefaultlib:" libcmtd.lib"\
 /out:"$(OUTDIR)/winace.exe" 
LINK32_OBJS= \
	"$(INTDIR)/graphwin32_.obj" \
	"$(INTDIR)/lexsubs.obj" \
	"$(INTDIR)/queryexe.obj" \
	"$(INTDIR)/vmapdisp.obj" \
	"$(INTDIR)/acefiledoc.obj" \
	"$(INTDIR)/metab.obj" \
	"$(INTDIR)/stdafx.obj" \
	"$(INTDIR)/griddisp.obj" \
	"$(INTDIR)/mainfrm.obj" \
	"$(INTDIR)/drawdisp.obj" \
	"$(INTDIR)/querydisp.obj" \
	"$(INTDIR)/vmapphys.obj" \
	"$(INTDIR)/oldhelp.obj" \
	"$(INTDIR)/graphsub.obj" \
	"$(INTDIR)/heap.obj" \
	"$(INTDIR)/sprddisplay.obj" \
	"$(INTDIR)/gmapposnegcol.obj" \
	"$(INTDIR)/newkey.obj" \
	"$(INTDIR)/menu.obj" \
	"$(INTDIR)/gmapdisp.obj" \
	"$(INTDIR)/gmapintervalcol.obj" \
	"$(INTDIR)/win32print.obj" \
	"$(INTDIR)/translate.obj" \
	"$(INTDIR)/acedialogs.obj" \
	"$(INTDIR)/disknew.obj" \
	"$(INTDIR)/gmapphys.obj" \
	"$(INTDIR)/dnasubs.obj" \
	"$(INTDIR)/textfitview.obj" \
	"$(INTDIR)/gmaplocuscol.obj" \
	"$(INTDIR)/call.obj" \
	"$(INTDIR)/graphselect.obj" \
	"$(INTDIR)/querybuild.obj" \
	"$(INTDIR)/arraysub.obj" \
	"$(INTDIR)/fontpreference.obj" \
	"$(INTDIR)/fmapgene.obj" \
	"$(INTDIR)/freesubs.obj" \
	"$(INTDIR)/chrono.obj" \
	"$(INTDIR)/gmapdata.obj" \
	"$(INTDIR)/fmapfeatures.obj" \
	"$(INTDIR)/gmapmarkercol.obj" \
	"$(INTDIR)/windialogs.obj" \
	"$(INTDIR)/graphwin.obj" \
	"$(INTDIR)/plainview.obj" \
	"$(INTDIR)/sysclass.obj" \
	"$(INTDIR)/parse.obj" \
	"$(INTDIR)/align.obj" \
	"$(INTDIR)/vmapdrag.obj" \
	"$(INTDIR)/sprdmap.obj" \
	"$(INTDIR)/timesubs.obj" \
	"$(INTDIR)/dnacpt.obj" \
	"$(INTDIR)/colourbox.obj" \
	"$(INTDIR)/check.obj" \
	"$(INTDIR)/fpdisp.obj" \
	"$(INTDIR)/preferences.obj" \
	"$(INTDIR)/fmapsequence.obj" \
	"$(INTDIR)/pepdisp.obj" \
	"$(INTDIR)/hexcode.obj" \
	"$(INTDIR)/cmapdisp.obj" \
	"$(INTDIR)/bump.obj" \
	"$(INTDIR)/viewedit.obj" \
	"$(INTDIR)/sprdop.obj" \
	"$(INTDIR)/session.obj" \
	"$(INTDIR)/alignment.obj" \
	"$(INTDIR)/gfcode.obj" \
	"$(INTDIR)/tags.obj" \
	"$(INTDIR)/banner.obj" \
	"$(INTDIR)/cmultiass.obj" \
	"$(INTDIR)/pepseqcol.obj" \
	"$(INTDIR)/win32.obj" \
	"$(INTDIR)/status.obj" \
	"$(INTDIR)/vmapdata2.obj" \
	"$(INTDIR)/mapscrollview.obj" \
	"$(INTDIR)/freeout.obj" \
	"$(INTDIR)/fitview.obj" \
	"$(INTDIR)/dict.obj" \
	"$(INTDIR)/sprddata.obj" \
	"$(INTDIR)/mapcontrol.obj" \
	"$(INTDIR)/textscrollview.obj" \
	"$(INTDIR)/fmapcontrol.obj" \
	"$(INTDIR)/embl.obj" \
	"$(INTDIR)/pepgraphcol.obj" \
	"$(INTDIR)/win32process.obj" \
	"$(INTDIR)/lexalpha.obj" \
	"$(INTDIR)/cgraph.obj" \
	"$(INTDIR)/picksubs.obj" \
	"$(INTDIR)/linkdate.obj" \
	"$(INTDIR)/bsubs.obj" \
	"$(INTDIR)/biblio.obj" \
	"$(INTDIR)/qbedisp.obj" \
	"$(INTDIR)/objcache.obj" \
	"$(INTDIR)/winace.obj" \
	"$(INTDIR)/acedbprofile.obj" \
	"$(INTDIR)/keyset.obj" \
	"$(INTDIR)/pephomolcol.obj" \
	"$(INTDIR)/blxview.obj" \
	"$(INTDIR)/randsubs.obj" \
	"$(INTDIR)/nicedump.obj" \
	"$(INTDIR)/display.obj" \
	"$(INTDIR)/blocksub.obj" \
	"$(INTDIR)/longtext.obj" \
	"$(INTDIR)/sprdctrl.obj" \
	"$(INTDIR)/msgwindow.obj" \
	"$(INTDIR)/afxtrace.obj" \
	"$(INTDIR)/filsubs.obj" \
	"$(INTDIR)/dotterKarlin.obj" \
	"$(INTDIR)/userinterface.obj" \
	"$(INTDIR)/winfilquery.obj" \
	"$(INTDIR)/dotter.obj" \
	"$(INTDIR)/graphcon.obj" \
	"$(INTDIR)/acedbprofileview.obj" \
	"$(INTDIR)/quovadis.obj" \
	"$(INTDIR)/caceprintpage.obj" \
	"$(INTDIR)/splashbox.obj" \
	"$(INTDIR)/gmapconvert.obj" \
	"$(INTDIR)/messubs.obj" \
	"$(INTDIR)/peptide.obj" \
	"$(INTDIR)/bstools.obj" \
	"$(INTDIR)/graphbox.obj" \
	"$(INTDIR)/colcontrol.obj" \
	"$(INTDIR)/mainpick.obj" \
	"$(INTDIR)/action.obj" \
	"$(INTDIR)/pmapdisp.obj" \
	"$(INTDIR)/graphview.obj" \
	"$(INTDIR)/geldisp.obj" \
	"$(INTDIR)/gmapdatacol.obj" \
	"$(INTDIR)/pixelfitview.obj" \
	"$(INTDIR)/pixelscrollview.obj" \
	"$(INTDIR)/class.obj" \
	"$(INTDIR)/ksetdisp.obj" \
	"$(INTDIR)/acefileview.obj" \
	"$(INTDIR)/lexsubs4.obj" \
	"$(INTDIR)/win32menusubs.obj" \
	"$(INTDIR)/dbprofile.obj" \
	"$(INTDIR)/pmapconvert.obj" \
	"$(INTDIR)/method.obj" \
	"$(INTDIR)/model.obj" \
	"$(INTDIR)/update.obj" \
	"$(INTDIR)/treedisp.obj" \
	"$(INTDIR)/keysetdump.obj" \
	"$(INTDIR)/win32lib.obj" \
	"$(INTDIR)/asubs.obj" \
	"$(INTDIR)/opp.obj" \
	"$(INTDIR)/dump.obj" \
	"$(INTDIR)/gmapremarkcol.obj" \
	"$(INTDIR)/plot.obj" \
	"$(INTDIR)/bssubs.obj" \
	"$(INTDIR)/graphramp.obj" \
	"$(INTDIR)/winmain.obj" \
	"$(INTDIR)/fullscrollview.obj" \
	"$(INTDIR)/gd.obj" \
	"$(INTDIR)/bsdumps.obj" \
	"$(INTDIR)/help.obj" \
	"$(INTDIR)/bstree.obj" \
	"$(INTDIR)/winace.res"

"$(OUTDIR)\winace.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "windebug"
# PROP BASE Intermediate_Dir "windebug"
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "c:\acedb\bin\windebug"
# PROP Intermediate_Dir "c:\acedb\bin\windebug"
OUTDIR=c:\acedb\bin\windebug
INTDIR=c:\acedb\bin\windebug

ALL : "$(OUTDIR)\winace.exe" "..\bin\winace\windebug\winace.bsc"

CLEAN : 
	-@erase "c:\acedb\bin\windebug\vc40.pdb"
	-@erase "c:\acedb\bin\windebug\vc40.idb"
	-@erase "..\bin\winace\windebug\winace.bsc"
	-@erase "c:\acedb\bin\windebug\graphwin32_.sbr"
	-@erase "c:\acedb\bin\windebug\cgraph.sbr"
	-@erase "c:\acedb\bin\windebug\newkey.sbr"
	-@erase "c:\acedb\bin\windebug\alignment.sbr"
	-@erase "c:\acedb\bin\windebug\sprddata.sbr"
	-@erase "c:\acedb\bin\windebug\cmultiass.sbr"
	-@erase "c:\acedb\bin\windebug\caceprintpage.sbr"
	-@erase "c:\acedb\bin\windebug\pixelfitview.sbr"
	-@erase "c:\acedb\bin\windebug\pepseqcol.sbr"
	-@erase "c:\acedb\bin\windebug\vmapdata2.sbr"
	-@erase "c:\acedb\bin\windebug\qbedisp.sbr"
	-@erase "c:\acedb\bin\windebug\sprddisplay.sbr"
	-@erase "c:\acedb\bin\windebug\dump.sbr"
	-@erase "c:\acedb\bin\windebug\metab.sbr"
	-@erase "c:\acedb\bin\windebug\plot.sbr"
	-@erase "c:\acedb\bin\windebug\disknew.sbr"
	-@erase "c:\acedb\bin\windebug\graphramp.sbr"
	-@erase "c:\acedb\bin\windebug\linkdate.sbr"
	-@erase "c:\acedb\bin\windebug\chrono.sbr"
	-@erase "c:\acedb\bin\windebug\objcache.sbr"
	-@erase "c:\acedb\bin\windebug\vmapdisp.sbr"
	-@erase "c:\acedb\bin\windebug\win32menusubs.sbr"
	-@erase "c:\acedb\bin\windebug\textfitview.sbr"
	-@erase "c:\acedb\bin\windebug\help.sbr"
	-@erase "c:\acedb\bin\windebug\gd.sbr"
	-@erase "c:\acedb\bin\windebug\graphselect.sbr"
	-@erase "c:\acedb\bin\windebug\blocksub.sbr"
	-@erase "c:\acedb\bin\windebug\vmapphys.sbr"
	-@erase "c:\acedb\bin\windebug\gmapremarkcol.sbr"
	-@erase "c:\acedb\bin\windebug\dnacpt.sbr"
	-@erase "c:\acedb\bin\windebug\sprdmap.sbr"
	-@erase "c:\acedb\bin\windebug\fpdisp.sbr"
	-@erase "c:\acedb\bin\windebug\fullscrollview.sbr"
	-@erase "c:\acedb\bin\windebug\peptide.sbr"
	-@erase "c:\acedb\bin\windebug\bstools.sbr"
	-@erase "c:\acedb\bin\windebug\gmapdisp.sbr"
	-@erase "c:\acedb\bin\windebug\pepdisp.sbr"
	-@erase "c:\acedb\bin\windebug\hexcode.sbr"
	-@erase "c:\acedb\bin\windebug\heap.sbr"
	-@erase "c:\acedb\bin\windebug\acefileview.sbr"
	-@erase "c:\acedb\bin\windebug\sprdop.sbr"
	-@erase "c:\acedb\bin\windebug\geldisp.sbr"
	-@erase "c:\acedb\bin\windebug\menu.sbr"
	-@erase "c:\acedb\bin\windebug\gmapphys.sbr"
	-@erase "c:\acedb\bin\windebug\gfcode.sbr"
	-@erase "c:\acedb\bin\windebug\banner.sbr"
	-@erase "c:\acedb\bin\windebug\session.sbr"
	-@erase "c:\acedb\bin\windebug\preferences.sbr"
	-@erase "c:\acedb\bin\windebug\class.sbr"
	-@erase "c:\acedb\bin\windebug\status.sbr"
	-@erase "c:\acedb\bin\windebug\parse.sbr"
	-@erase "c:\acedb\bin\windebug\graphbox.sbr"
	-@erase "c:\acedb\bin\windebug\align.sbr"
	-@erase "c:\acedb\bin\windebug\freesubs.sbr"
	-@erase "c:\acedb\bin\windebug\win32print.sbr"
	-@erase "c:\acedb\bin\windebug\gmapdata.sbr"
	-@erase "c:\acedb\bin\windebug\acedialogs.sbr"
	-@erase "c:\acedb\bin\windebug\freeout.sbr"
	-@erase "c:\acedb\bin\windebug\fitview.sbr"
	-@erase "c:\acedb\bin\windebug\model.sbr"
	-@erase "c:\acedb\bin\windebug\check.sbr"
	-@erase "c:\acedb\bin\windebug\querybuild.sbr"
	-@erase "c:\acedb\bin\windebug\colcontrol.sbr"
	-@erase "c:\acedb\bin\windebug\lexsubs4.sbr"
	-@erase "c:\acedb\bin\windebug\gmapposnegcol.sbr"
	-@erase "c:\acedb\bin\windebug\sysclass.sbr"
	-@erase "c:\acedb\bin\windebug\gmaplocuscol.sbr"
	-@erase "c:\acedb\bin\windebug\timesubs.sbr"
	-@erase "c:\acedb\bin\windebug\winmain.sbr"
	-@erase "c:\acedb\bin\windebug\opp.sbr"
	-@erase "c:\acedb\bin\windebug\treedisp.sbr"
	-@erase "c:\acedb\bin\windebug\win32lib.sbr"
	-@erase "c:\acedb\bin\windebug\dotterKarlin.sbr"
	-@erase "c:\acedb\bin\windebug\win32.sbr"
	-@erase "c:\acedb\bin\windebug\fmapfeatures.sbr"
	-@erase "c:\acedb\bin\windebug\biblio.sbr"
	-@erase "c:\acedb\bin\windebug\fmapcontrol.sbr"
	-@erase "c:\acedb\bin\windebug\fontpreference.sbr"
	-@erase "c:\acedb\bin\windebug\winace.sbr"
	-@erase "c:\acedb\bin\windebug\dbprofile.sbr"
	-@erase "c:\acedb\bin\windebug\keysetdump.sbr"
	-@erase "c:\acedb\bin\windebug\keyset.sbr"
	-@erase "c:\acedb\bin\windebug\pepgraphcol.sbr"
	-@erase "c:\acedb\bin\windebug\bump.sbr"
	-@erase "c:\acedb\bin\windebug\blxview.sbr"
	-@erase "c:\acedb\bin\windebug\gmapmarkercol.sbr"
	-@erase "c:\acedb\bin\windebug\mainfrm.sbr"
	-@erase "c:\acedb\bin\windebug\fmapsequence.sbr"
	-@erase "c:\acedb\bin\windebug\display.sbr"
	-@erase "c:\acedb\bin\windebug\oldhelp.sbr"
	-@erase "c:\acedb\bin\windebug\acedbprofileview.sbr"
	-@erase "c:\acedb\bin\windebug\tags.sbr"
	-@erase "c:\acedb\bin\windebug\queryexe.sbr"
	-@erase "c:\acedb\bin\windebug\pephomolcol.sbr"
	-@erase "c:\acedb\bin\windebug\filsubs.sbr"
	-@erase "c:\acedb\bin\windebug\dotter.sbr"
	-@erase "c:\acedb\bin\windebug\gmapintervalcol.sbr"
	-@erase "c:\acedb\bin\windebug\dict.sbr"
	-@erase "c:\acedb\bin\windebug\bsubs.sbr"
	-@erase "c:\acedb\bin\windebug\griddisp.sbr"
	-@erase "c:\acedb\bin\windebug\drawdisp.sbr"
	-@erase "c:\acedb\bin\windebug\lexalpha.sbr"
	-@erase "c:\acedb\bin\windebug\graphsub.sbr"
	-@erase "c:\acedb\bin\windebug\dnasubs.sbr"
	-@erase "c:\acedb\bin\windebug\mapcontrol.sbr"
	-@erase "c:\acedb\bin\windebug\acefiledoc.sbr"
	-@erase "c:\acedb\bin\windebug\winfilquery.sbr"
	-@erase "c:\acedb\bin\windebug\picksubs.sbr"
	-@erase "c:\acedb\bin\windebug\messubs.sbr"
	-@erase "c:\acedb\bin\windebug\embl.sbr"
	-@erase "c:\acedb\bin\windebug\action.sbr"
	-@erase "c:\acedb\bin\windebug\querydisp.sbr"
	-@erase "c:\acedb\bin\windebug\pixelscrollview.sbr"
	-@erase "c:\acedb\bin\windebug\randsubs.sbr"
	-@erase "c:\acedb\bin\windebug\gmapconvert.sbr"
	-@erase "c:\acedb\bin\windebug\mapscrollview.sbr"
	-@erase "c:\acedb\bin\windebug\nicedump.sbr"
	-@erase "c:\acedb\bin\windebug\win32process.sbr"
	-@erase "c:\acedb\bin\windebug\textscrollview.sbr"
	-@erase "c:\acedb\bin\windebug\longtext.sbr"
	-@erase "c:\acedb\bin\windebug\sprdctrl.sbr"
	-@erase "c:\acedb\bin\windebug\afxtrace.sbr"
	-@erase "c:\acedb\bin\windebug\method.sbr"
	-@erase "c:\acedb\bin\windebug\arraysub.sbr"
	-@erase "c:\acedb\bin\windebug\update.sbr"
	-@erase "c:\acedb\bin\windebug\fmapgene.sbr"
	-@erase "c:\acedb\bin\windebug\gmapdatacol.sbr"
	-@erase "c:\acedb\bin\windebug\translate.sbr"
	-@erase "c:\acedb\bin\windebug\acedbprofile.sbr"
	-@erase "c:\acedb\bin\windebug\call.sbr"
	-@erase "c:\acedb\bin\windebug\graphcon.sbr"
	-@erase "c:\acedb\bin\windebug\msgwindow.sbr"
	-@erase "c:\acedb\bin\windebug\graphwin.sbr"
	-@erase "c:\acedb\bin\windebug\bssubs.sbr"
	-@erase "c:\acedb\bin\windebug\quovadis.sbr"
	-@erase "c:\acedb\bin\windebug\pmapconvert.sbr"
	-@erase "c:\acedb\bin\windebug\vmapdrag.sbr"
	-@erase "c:\acedb\bin\windebug\windialogs.sbr"
	-@erase "c:\acedb\bin\windebug\mainpick.sbr"
	-@erase "c:\acedb\bin\windebug\cmapdisp.sbr"
	-@erase "c:\acedb\bin\windebug\pmapdisp.sbr"
	-@erase "c:\acedb\bin\windebug\bstree.sbr"
	-@erase "c:\acedb\bin\windebug\bsdumps.sbr"
	-@erase "c:\acedb\bin\windebug\viewedit.sbr"
	-@erase "c:\acedb\bin\windebug\plainview.sbr"
	-@erase "c:\acedb\bin\windebug\asubs.sbr"
	-@erase "c:\acedb\bin\windebug\splashbox.sbr"
	-@erase "c:\acedb\bin\windebug\lexsubs.sbr"
	-@erase "c:\acedb\bin\windebug\ksetdisp.sbr"
	-@erase "c:\acedb\bin\windebug\colourbox.sbr"
	-@erase "c:\acedb\bin\windebug\stdafx.sbr"
	-@erase "c:\acedb\bin\windebug\graphview.sbr"
	-@erase "c:\acedb\bin\windebug\userinterface.sbr"
	-@erase "c:\acedb\bin\windebug\winace.exe"
	-@erase "c:\acedb\bin\windebug\arraysub.obj"
	-@erase "c:\acedb\bin\windebug\update.obj"
	-@erase "c:\acedb\bin\windebug\fmapgene.obj"
	-@erase "c:\acedb\bin\windebug\gmapdatacol.obj"
	-@erase "c:\acedb\bin\windebug\translate.obj"
	-@erase "c:\acedb\bin\windebug\acedbprofile.obj"
	-@erase "c:\acedb\bin\windebug\call.obj"
	-@erase "c:\acedb\bin\windebug\graphcon.obj"
	-@erase "c:\acedb\bin\windebug\msgwindow.obj"
	-@erase "c:\acedb\bin\windebug\graphwin.obj"
	-@erase "c:\acedb\bin\windebug\bssubs.obj"
	-@erase "c:\acedb\bin\windebug\quovadis.obj"
	-@erase "c:\acedb\bin\windebug\pmapconvert.obj"
	-@erase "c:\acedb\bin\windebug\vmapdrag.obj"
	-@erase "c:\acedb\bin\windebug\windialogs.obj"
	-@erase "c:\acedb\bin\windebug\mainpick.obj"
	-@erase "c:\acedb\bin\windebug\cmapdisp.obj"
	-@erase "c:\acedb\bin\windebug\pmapdisp.obj"
	-@erase "c:\acedb\bin\windebug\bstree.obj"
	-@erase "c:\acedb\bin\windebug\bsdumps.obj"
	-@erase "c:\acedb\bin\windebug\viewedit.obj"
	-@erase "c:\acedb\bin\windebug\plainview.obj"
	-@erase "c:\acedb\bin\windebug\asubs.obj"
	-@erase "c:\acedb\bin\windebug\splashbox.obj"
	-@erase "c:\acedb\bin\windebug\lexsubs.obj"
	-@erase "c:\acedb\bin\windebug\ksetdisp.obj"
	-@erase "c:\acedb\bin\windebug\colourbox.obj"
	-@erase "c:\acedb\bin\windebug\stdafx.obj"
	-@erase "c:\acedb\bin\windebug\graphview.obj"
	-@erase "c:\acedb\bin\windebug\userinterface.obj"
	-@erase "c:\acedb\bin\windebug\graphwin32_.obj"
	-@erase "c:\acedb\bin\windebug\cgraph.obj"
	-@erase "c:\acedb\bin\windebug\newkey.obj"
	-@erase "c:\acedb\bin\windebug\alignment.obj"
	-@erase "c:\acedb\bin\windebug\sprddata.obj"
	-@erase "c:\acedb\bin\windebug\cmultiass.obj"
	-@erase "c:\acedb\bin\windebug\caceprintpage.obj"
	-@erase "c:\acedb\bin\windebug\pixelfitview.obj"
	-@erase "c:\acedb\bin\windebug\pepseqcol.obj"
	-@erase "c:\acedb\bin\windebug\vmapdata2.obj"
	-@erase "c:\acedb\bin\windebug\qbedisp.obj"
	-@erase "c:\acedb\bin\windebug\sprddisplay.obj"
	-@erase "c:\acedb\bin\windebug\dump.obj"
	-@erase "c:\acedb\bin\windebug\metab.obj"
	-@erase "c:\acedb\bin\windebug\plot.obj"
	-@erase "c:\acedb\bin\windebug\disknew.obj"
	-@erase "c:\acedb\bin\windebug\graphramp.obj"
	-@erase "c:\acedb\bin\windebug\linkdate.obj"
	-@erase "c:\acedb\bin\windebug\chrono.obj"
	-@erase "c:\acedb\bin\windebug\objcache.obj"
	-@erase "c:\acedb\bin\windebug\vmapdisp.obj"
	-@erase "c:\acedb\bin\windebug\win32menusubs.obj"
	-@erase "c:\acedb\bin\windebug\textfitview.obj"
	-@erase "c:\acedb\bin\windebug\help.obj"
	-@erase "c:\acedb\bin\windebug\gd.obj"
	-@erase "c:\acedb\bin\windebug\graphselect.obj"
	-@erase "c:\acedb\bin\windebug\blocksub.obj"
	-@erase "c:\acedb\bin\windebug\vmapphys.obj"
	-@erase "c:\acedb\bin\windebug\gmapremarkcol.obj"
	-@erase "c:\acedb\bin\windebug\dnacpt.obj"
	-@erase "c:\acedb\bin\windebug\sprdmap.obj"
	-@erase "c:\acedb\bin\windebug\fpdisp.obj"
	-@erase "c:\acedb\bin\windebug\fullscrollview.obj"
	-@erase "c:\acedb\bin\windebug\peptide.obj"
	-@erase "c:\acedb\bin\windebug\bstools.obj"
	-@erase "c:\acedb\bin\windebug\gmapdisp.obj"
	-@erase "c:\acedb\bin\windebug\pepdisp.obj"
	-@erase "c:\acedb\bin\windebug\hexcode.obj"
	-@erase "c:\acedb\bin\windebug\heap.obj"
	-@erase "c:\acedb\bin\windebug\acefileview.obj"
	-@erase "c:\acedb\bin\windebug\sprdop.obj"
	-@erase "c:\acedb\bin\windebug\geldisp.obj"
	-@erase "c:\acedb\bin\windebug\menu.obj"
	-@erase "c:\acedb\bin\windebug\gmapphys.obj"
	-@erase "c:\acedb\bin\windebug\gfcode.obj"
	-@erase "c:\acedb\bin\windebug\banner.obj"
	-@erase "c:\acedb\bin\windebug\session.obj"
	-@erase "c:\acedb\bin\windebug\preferences.obj"
	-@erase "c:\acedb\bin\windebug\class.obj"
	-@erase "c:\acedb\bin\windebug\status.obj"
	-@erase "c:\acedb\bin\windebug\parse.obj"
	-@erase "c:\acedb\bin\windebug\graphbox.obj"
	-@erase "c:\acedb\bin\windebug\align.obj"
	-@erase "c:\acedb\bin\windebug\freesubs.obj"
	-@erase "c:\acedb\bin\windebug\win32print.obj"
	-@erase "c:\acedb\bin\windebug\gmapdata.obj"
	-@erase "c:\acedb\bin\windebug\acedialogs.obj"
	-@erase "c:\acedb\bin\windebug\freeout.obj"
	-@erase "c:\acedb\bin\windebug\fitview.obj"
	-@erase "c:\acedb\bin\windebug\model.obj"
	-@erase "c:\acedb\bin\windebug\check.obj"
	-@erase "c:\acedb\bin\windebug\querybuild.obj"
	-@erase "c:\acedb\bin\windebug\colcontrol.obj"
	-@erase "c:\acedb\bin\windebug\lexsubs4.obj"
	-@erase "c:\acedb\bin\windebug\gmapposnegcol.obj"
	-@erase "c:\acedb\bin\windebug\sysclass.obj"
	-@erase "c:\acedb\bin\windebug\gmaplocuscol.obj"
	-@erase "c:\acedb\bin\windebug\timesubs.obj"
	-@erase "c:\acedb\bin\windebug\winmain.obj"
	-@erase "c:\acedb\bin\windebug\opp.obj"
	-@erase "c:\acedb\bin\windebug\treedisp.obj"
	-@erase "c:\acedb\bin\windebug\win32lib.obj"
	-@erase "c:\acedb\bin\windebug\dotterKarlin.obj"
	-@erase "c:\acedb\bin\windebug\win32.obj"
	-@erase "c:\acedb\bin\windebug\fmapfeatures.obj"
	-@erase "c:\acedb\bin\windebug\biblio.obj"
	-@erase "c:\acedb\bin\windebug\fmapcontrol.obj"
	-@erase "c:\acedb\bin\windebug\fontpreference.obj"
	-@erase "c:\acedb\bin\windebug\winace.obj"
	-@erase "c:\acedb\bin\windebug\dbprofile.obj"
	-@erase "c:\acedb\bin\windebug\keysetdump.obj"
	-@erase "c:\acedb\bin\windebug\keyset.obj"
	-@erase "c:\acedb\bin\windebug\pepgraphcol.obj"
	-@erase "c:\acedb\bin\windebug\bump.obj"
	-@erase "c:\acedb\bin\windebug\blxview.obj"
	-@erase "c:\acedb\bin\windebug\gmapmarkercol.obj"
	-@erase "c:\acedb\bin\windebug\mainfrm.obj"
	-@erase "c:\acedb\bin\windebug\fmapsequence.obj"
	-@erase "c:\acedb\bin\windebug\display.obj"
	-@erase "c:\acedb\bin\windebug\oldhelp.obj"
	-@erase "c:\acedb\bin\windebug\acedbprofileview.obj"
	-@erase "c:\acedb\bin\windebug\tags.obj"
	-@erase "c:\acedb\bin\windebug\queryexe.obj"
	-@erase "c:\acedb\bin\windebug\pephomolcol.obj"
	-@erase "c:\acedb\bin\windebug\filsubs.obj"
	-@erase "c:\acedb\bin\windebug\dotter.obj"
	-@erase "c:\acedb\bin\windebug\gmapintervalcol.obj"
	-@erase "c:\acedb\bin\windebug\dict.obj"
	-@erase "c:\acedb\bin\windebug\bsubs.obj"
	-@erase "c:\acedb\bin\windebug\griddisp.obj"
	-@erase "c:\acedb\bin\windebug\drawdisp.obj"
	-@erase "c:\acedb\bin\windebug\lexalpha.obj"
	-@erase "c:\acedb\bin\windebug\graphsub.obj"
	-@erase "c:\acedb\bin\windebug\dnasubs.obj"
	-@erase "c:\acedb\bin\windebug\mapcontrol.obj"
	-@erase "c:\acedb\bin\windebug\acefiledoc.obj"
	-@erase "c:\acedb\bin\windebug\winfilquery.obj"
	-@erase "c:\acedb\bin\windebug\picksubs.obj"
	-@erase "c:\acedb\bin\windebug\messubs.obj"
	-@erase "c:\acedb\bin\windebug\embl.obj"
	-@erase "c:\acedb\bin\windebug\action.obj"
	-@erase "c:\acedb\bin\windebug\querydisp.obj"
	-@erase "c:\acedb\bin\windebug\pixelscrollview.obj"
	-@erase "c:\acedb\bin\windebug\randsubs.obj"
	-@erase "c:\acedb\bin\windebug\gmapconvert.obj"
	-@erase "c:\acedb\bin\windebug\mapscrollview.obj"
	-@erase "c:\acedb\bin\windebug\nicedump.obj"
	-@erase "c:\acedb\bin\windebug\win32process.obj"
	-@erase "c:\acedb\bin\windebug\textscrollview.obj"
	-@erase "c:\acedb\bin\windebug\longtext.obj"
	-@erase "c:\acedb\bin\windebug\sprdctrl.obj"
	-@erase "c:\acedb\bin\windebug\afxtrace.obj"
	-@erase "c:\acedb\bin\windebug\method.obj"
	-@erase "c:\acedb\bin\windebug\winace.res"
	-@erase "c:\acedb\bin\windebug\winace.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR /YX /c
CPP_PROJ=/nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/winace.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=c:\acedb\bin\windebug/
CPP_SBRS=c:\acedb\bin\windebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"\acedb\bin\winace\windebug\winace.bsc"
BSC32_FLAGS=/nologo /o"\acedb\bin\winace\windebug\winace.bsc" 
BSC32_SBRS= \
	"$(INTDIR)/graphwin32_.sbr" \
	"$(INTDIR)/cgraph.sbr" \
	"$(INTDIR)/newkey.sbr" \
	"$(INTDIR)/alignment.sbr" \
	"$(INTDIR)/sprddata.sbr" \
	"$(INTDIR)/cmultiass.sbr" \
	"$(INTDIR)/caceprintpage.sbr" \
	"$(INTDIR)/pixelfitview.sbr" \
	"$(INTDIR)/pepseqcol.sbr" \
	"$(INTDIR)/vmapdata2.sbr" \
	"$(INTDIR)/qbedisp.sbr" \
	"$(INTDIR)/sprddisplay.sbr" \
	"$(INTDIR)/dump.sbr" \
	"$(INTDIR)/metab.sbr" \
	"$(INTDIR)/plot.sbr" \
	"$(INTDIR)/disknew.sbr" \
	"$(INTDIR)/graphramp.sbr" \
	"$(INTDIR)/linkdate.sbr" \
	"$(INTDIR)/chrono.sbr" \
	"$(INTDIR)/objcache.sbr" \
	"$(INTDIR)/vmapdisp.sbr" \
	"$(INTDIR)/win32menusubs.sbr" \
	"$(INTDIR)/textfitview.sbr" \
	"$(INTDIR)/help.sbr" \
	"$(INTDIR)/gd.sbr" \
	"$(INTDIR)/graphselect.sbr" \
	"$(INTDIR)/blocksub.sbr" \
	"$(INTDIR)/vmapphys.sbr" \
	"$(INTDIR)/gmapremarkcol.sbr" \
	"$(INTDIR)/dnacpt.sbr" \
	"$(INTDIR)/sprdmap.sbr" \
	"$(INTDIR)/fpdisp.sbr" \
	"$(INTDIR)/fullscrollview.sbr" \
	"$(INTDIR)/peptide.sbr" \
	"$(INTDIR)/bstools.sbr" \
	"$(INTDIR)/gmapdisp.sbr" \
	"$(INTDIR)/pepdisp.sbr" \
	"$(INTDIR)/hexcode.sbr" \
	"$(INTDIR)/heap.sbr" \
	"$(INTDIR)/acefileview.sbr" \
	"$(INTDIR)/sprdop.sbr" \
	"$(INTDIR)/geldisp.sbr" \
	"$(INTDIR)/menu.sbr" \
	"$(INTDIR)/gmapphys.sbr" \
	"$(INTDIR)/gfcode.sbr" \
	"$(INTDIR)/banner.sbr" \
	"$(INTDIR)/session.sbr" \
	"$(INTDIR)/preferences.sbr" \
	"$(INTDIR)/class.sbr" \
	"$(INTDIR)/status.sbr" \
	"$(INTDIR)/parse.sbr" \
	"$(INTDIR)/graphbox.sbr" \
	"$(INTDIR)/align.sbr" \
	"$(INTDIR)/freesubs.sbr" \
	"$(INTDIR)/win32print.sbr" \
	"$(INTDIR)/gmapdata.sbr" \
	"$(INTDIR)/acedialogs.sbr" \
	"$(INTDIR)/freeout.sbr" \
	"$(INTDIR)/fitview.sbr" \
	"$(INTDIR)/model.sbr" \
	"$(INTDIR)/check.sbr" \
	"$(INTDIR)/querybuild.sbr" \
	"$(INTDIR)/colcontrol.sbr" \
	"$(INTDIR)/lexsubs4.sbr" \
	"$(INTDIR)/gmapposnegcol.sbr" \
	"$(INTDIR)/sysclass.sbr" \
	"$(INTDIR)/gmaplocuscol.sbr" \
	"$(INTDIR)/timesubs.sbr" \
	"$(INTDIR)/winmain.sbr" \
	"$(INTDIR)/opp.sbr" \
	"$(INTDIR)/treedisp.sbr" \
	"$(INTDIR)/win32lib.sbr" \
	"$(INTDIR)/dotterKarlin.sbr" \
	"$(INTDIR)/win32.sbr" \
	"$(INTDIR)/fmapfeatures.sbr" \
	"$(INTDIR)/biblio.sbr" \
	"$(INTDIR)/fmapcontrol.sbr" \
	"$(INTDIR)/fontpreference.sbr" \
	"$(INTDIR)/winace.sbr" \
	"$(INTDIR)/dbprofile.sbr" \
	"$(INTDIR)/keysetdump.sbr" \
	"$(INTDIR)/keyset.sbr" \
	"$(INTDIR)/pepgraphcol.sbr" \
	"$(INTDIR)/bump.sbr" \
	"$(INTDIR)/blxview.sbr" \
	"$(INTDIR)/gmapmarkercol.sbr" \
	"$(INTDIR)/mainfrm.sbr" \
	"$(INTDIR)/fmapsequence.sbr" \
	"$(INTDIR)/display.sbr" \
	"$(INTDIR)/oldhelp.sbr" \
	"$(INTDIR)/acedbprofileview.sbr" \
	"$(INTDIR)/tags.sbr" \
	"$(INTDIR)/queryexe.sbr" \
	"$(INTDIR)/pephomolcol.sbr" \
	"$(INTDIR)/filsubs.sbr" \
	"$(INTDIR)/dotter.sbr" \
	"$(INTDIR)/gmapintervalcol.sbr" \
	"$(INTDIR)/dict.sbr" \
	"$(INTDIR)/bsubs.sbr" \
	"$(INTDIR)/griddisp.sbr" \
	"$(INTDIR)/drawdisp.sbr" \
	"$(INTDIR)/lexalpha.sbr" \
	"$(INTDIR)/graphsub.sbr" \
	"$(INTDIR)/dnasubs.sbr" \
	"$(INTDIR)/mapcontrol.sbr" \
	"$(INTDIR)/acefiledoc.sbr" \
	"$(INTDIR)/winfilquery.sbr" \
	"$(INTDIR)/picksubs.sbr" \
	"$(INTDIR)/messubs.sbr" \
	"$(INTDIR)/embl.sbr" \
	"$(INTDIR)/action.sbr" \
	"$(INTDIR)/querydisp.sbr" \
	"$(INTDIR)/pixelscrollview.sbr" \
	"$(INTDIR)/randsubs.sbr" \
	"$(INTDIR)/gmapconvert.sbr" \
	"$(INTDIR)/mapscrollview.sbr" \
	"$(INTDIR)/nicedump.sbr" \
	"$(INTDIR)/win32process.sbr" \
	"$(INTDIR)/textscrollview.sbr" \
	"$(INTDIR)/longtext.sbr" \
	"$(INTDIR)/sprdctrl.sbr" \
	"$(INTDIR)/afxtrace.sbr" \
	"$(INTDIR)/method.sbr" \
	"$(INTDIR)/arraysub.sbr" \
	"$(INTDIR)/update.sbr" \
	"$(INTDIR)/fmapgene.sbr" \
	"$(INTDIR)/gmapdatacol.sbr" \
	"$(INTDIR)/translate.sbr" \
	"$(INTDIR)/acedbprofile.sbr" \
	"$(INTDIR)/call.sbr" \
	"$(INTDIR)/graphcon.sbr" \
	"$(INTDIR)/msgwindow.sbr" \
	"$(INTDIR)/graphwin.sbr" \
	"$(INTDIR)/bssubs.sbr" \
	"$(INTDIR)/quovadis.sbr" \
	"$(INTDIR)/pmapconvert.sbr" \
	"$(INTDIR)/vmapdrag.sbr" \
	"$(INTDIR)/windialogs.sbr" \
	"$(INTDIR)/mainpick.sbr" \
	"$(INTDIR)/cmapdisp.sbr" \
	"$(INTDIR)/pmapdisp.sbr" \
	"$(INTDIR)/bstree.sbr" \
	"$(INTDIR)/bsdumps.sbr" \
	"$(INTDIR)/viewedit.sbr" \
	"$(INTDIR)/plainview.sbr" \
	"$(INTDIR)/asubs.sbr" \
	"$(INTDIR)/splashbox.sbr" \
	"$(INTDIR)/lexsubs.sbr" \
	"$(INTDIR)/ksetdisp.sbr" \
	"$(INTDIR)/colourbox.sbr" \
	"$(INTDIR)/stdafx.sbr" \
	"$(INTDIR)/graphview.sbr" \
	"$(INTDIR)/userinterface.sbr"

"..\bin\winace\windebug\winace.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386
# ADD LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /profile /map /debug /machine:I386 /nodefaultlib:"libc.lib" /nodefaultlib:" libcmt.lib" /nodefaultlib:" msvcrtd.lib" /nodefaultlib:" libcd.lib" /nodefaultlib:" libcmtd.lib"
# SUBTRACT LINK32 /verbose /nodefaultlib
LINK32_FLAGS=/nologo /version:4.31 /entry:"WinMainCRTStartup"\
 /subsystem:windows /profile /map:"$(INTDIR)/winace.map" /debug /machine:I386\
 /nodefaultlib:"libc.lib" /nodefaultlib:" libcmt.lib"\
 /nodefaultlib:" msvcrtd.lib" /nodefaultlib:" libcd.lib"\
 /nodefaultlib:" libcmtd.lib" /out:"$(OUTDIR)/winace.exe" 
LINK32_OBJS= \
	"$(INTDIR)/arraysub.obj" \
	"$(INTDIR)/update.obj" \
	"$(INTDIR)/fmapgene.obj" \
	"$(INTDIR)/gmapdatacol.obj" \
	"$(INTDIR)/translate.obj" \
	"$(INTDIR)/acedbprofile.obj" \
	"$(INTDIR)/call.obj" \
	"$(INTDIR)/graphcon.obj" \
	"$(INTDIR)/msgwindow.obj" \
	"$(INTDIR)/graphwin.obj" \
	"$(INTDIR)/bssubs.obj" \
	"$(INTDIR)/quovadis.obj" \
	"$(INTDIR)/pmapconvert.obj" \
	"$(INTDIR)/vmapdrag.obj" \
	"$(INTDIR)/windialogs.obj" \
	"$(INTDIR)/mainpick.obj" \
	"$(INTDIR)/cmapdisp.obj" \
	"$(INTDIR)/pmapdisp.obj" \
	"$(INTDIR)/bstree.obj" \
	"$(INTDIR)/bsdumps.obj" \
	"$(INTDIR)/viewedit.obj" \
	"$(INTDIR)/plainview.obj" \
	"$(INTDIR)/asubs.obj" \
	"$(INTDIR)/splashbox.obj" \
	"$(INTDIR)/lexsubs.obj" \
	"$(INTDIR)/ksetdisp.obj" \
	"$(INTDIR)/colourbox.obj" \
	"$(INTDIR)/stdafx.obj" \
	"$(INTDIR)/graphview.obj" \
	"$(INTDIR)/userinterface.obj" \
	"$(INTDIR)/graphwin32_.obj" \
	"$(INTDIR)/cgraph.obj" \
	"$(INTDIR)/newkey.obj" \
	"$(INTDIR)/alignment.obj" \
	"$(INTDIR)/sprddata.obj" \
	"$(INTDIR)/cmultiass.obj" \
	"$(INTDIR)/caceprintpage.obj" \
	"$(INTDIR)/pixelfitview.obj" \
	"$(INTDIR)/pepseqcol.obj" \
	"$(INTDIR)/vmapdata2.obj" \
	"$(INTDIR)/qbedisp.obj" \
	"$(INTDIR)/sprddisplay.obj" \
	"$(INTDIR)/dump.obj" \
	"$(INTDIR)/metab.obj" \
	"$(INTDIR)/plot.obj" \
	"$(INTDIR)/disknew.obj" \
	"$(INTDIR)/graphramp.obj" \
	"$(INTDIR)/linkdate.obj" \
	"$(INTDIR)/chrono.obj" \
	"$(INTDIR)/objcache.obj" \
	"$(INTDIR)/vmapdisp.obj" \
	"$(INTDIR)/win32menusubs.obj" \
	"$(INTDIR)/textfitview.obj" \
	"$(INTDIR)/help.obj" \
	"$(INTDIR)/gd.obj" \
	"$(INTDIR)/graphselect.obj" \
	"$(INTDIR)/blocksub.obj" \
	"$(INTDIR)/vmapphys.obj" \
	"$(INTDIR)/gmapremarkcol.obj" \
	"$(INTDIR)/dnacpt.obj" \
	"$(INTDIR)/sprdmap.obj" \
	"$(INTDIR)/fpdisp.obj" \
	"$(INTDIR)/fullscrollview.obj" \
	"$(INTDIR)/peptide.obj" \
	"$(INTDIR)/bstools.obj" \
	"$(INTDIR)/gmapdisp.obj" \
	"$(INTDIR)/pepdisp.obj" \
	"$(INTDIR)/hexcode.obj" \
	"$(INTDIR)/heap.obj" \
	"$(INTDIR)/acefileview.obj" \
	"$(INTDIR)/sprdop.obj" \
	"$(INTDIR)/geldisp.obj" \
	"$(INTDIR)/menu.obj" \
	"$(INTDIR)/gmapphys.obj" \
	"$(INTDIR)/gfcode.obj" \
	"$(INTDIR)/banner.obj" \
	"$(INTDIR)/session.obj" \
	"$(INTDIR)/preferences.obj" \
	"$(INTDIR)/class.obj" \
	"$(INTDIR)/status.obj" \
	"$(INTDIR)/parse.obj" \
	"$(INTDIR)/graphbox.obj" \
	"$(INTDIR)/align.obj" \
	"$(INTDIR)/freesubs.obj" \
	"$(INTDIR)/win32print.obj" \
	"$(INTDIR)/gmapdata.obj" \
	"$(INTDIR)/acedialogs.obj" \
	"$(INTDIR)/freeout.obj" \
	"$(INTDIR)/fitview.obj" \
	"$(INTDIR)/model.obj" \
	"$(INTDIR)/check.obj" \
	"$(INTDIR)/querybuild.obj" \
	"$(INTDIR)/colcontrol.obj" \
	"$(INTDIR)/lexsubs4.obj" \
	"$(INTDIR)/gmapposnegcol.obj" \
	"$(INTDIR)/sysclass.obj" \
	"$(INTDIR)/gmaplocuscol.obj" \
	"$(INTDIR)/timesubs.obj" \
	"$(INTDIR)/winmain.obj" \
	"$(INTDIR)/opp.obj" \
	"$(INTDIR)/treedisp.obj" \
	"$(INTDIR)/win32lib.obj" \
	"$(INTDIR)/dotterKarlin.obj" \
	"$(INTDIR)/win32.obj" \
	"$(INTDIR)/fmapfeatures.obj" \
	"$(INTDIR)/biblio.obj" \
	"$(INTDIR)/fmapcontrol.obj" \
	"$(INTDIR)/fontpreference.obj" \
	"$(INTDIR)/winace.obj" \
	"$(INTDIR)/dbprofile.obj" \
	"$(INTDIR)/keysetdump.obj" \
	"$(INTDIR)/keyset.obj" \
	"$(INTDIR)/pepgraphcol.obj" \
	"$(INTDIR)/bump.obj" \
	"$(INTDIR)/blxview.obj" \
	"$(INTDIR)/gmapmarkercol.obj" \
	"$(INTDIR)/mainfrm.obj" \
	"$(INTDIR)/fmapsequence.obj" \
	"$(INTDIR)/display.obj" \
	"$(INTDIR)/oldhelp.obj" \
	"$(INTDIR)/acedbprofileview.obj" \
	"$(INTDIR)/tags.obj" \
	"$(INTDIR)/queryexe.obj" \
	"$(INTDIR)/pephomolcol.obj" \
	"$(INTDIR)/filsubs.obj" \
	"$(INTDIR)/dotter.obj" \
	"$(INTDIR)/gmapintervalcol.obj" \
	"$(INTDIR)/dict.obj" \
	"$(INTDIR)/bsubs.obj" \
	"$(INTDIR)/griddisp.obj" \
	"$(INTDIR)/drawdisp.obj" \
	"$(INTDIR)/lexalpha.obj" \
	"$(INTDIR)/graphsub.obj" \
	"$(INTDIR)/dnasubs.obj" \
	"$(INTDIR)/mapcontrol.obj" \
	"$(INTDIR)/acefiledoc.obj" \
	"$(INTDIR)/winfilquery.obj" \
	"$(INTDIR)/picksubs.obj" \
	"$(INTDIR)/messubs.obj" \
	"$(INTDIR)/embl.obj" \
	"$(INTDIR)/action.obj" \
	"$(INTDIR)/querydisp.obj" \
	"$(INTDIR)/pixelscrollview.obj" \
	"$(INTDIR)/randsubs.obj" \
	"$(INTDIR)/gmapconvert.obj" \
	"$(INTDIR)/mapscrollview.obj" \
	"$(INTDIR)/nicedump.obj" \
	"$(INTDIR)/win32process.obj" \
	"$(INTDIR)/textscrollview.obj" \
	"$(INTDIR)/longtext.obj" \
	"$(INTDIR)/sprdctrl.obj" \
	"$(INTDIR)/afxtrace.obj" \
	"$(INTDIR)/method.obj" \
	"$(INTDIR)/winace.res"

"$(OUTDIR)\winace.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

################################################################################
# Begin Target

# Name "WinAce - Win32 Release"
# Name "WinAce - Win32 Debug"

!IF  "$(CFG)" == "WinAce - Win32 Release"

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\acedb\w1\timesubs.c
DEP_CPP_TIMES=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\menu.c
DEP_CPP_MENU_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\menu_.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\menu.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\menu.obj" : $(SOURCE) $(DEP_CPP_MENU_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\menu.sbr" : $(SOURCE) $(DEP_CPP_MENU_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\menu.obj" : $(SOURCE) $(DEP_CPP_MENU_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\menu.sbr" : $(SOURCE) $(DEP_CPP_MENU_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\bump.c
DEP_CPP_BUMP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\dict.c
DEP_CPP_DICT_=\
	"\acedb\wh\dict.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\filsubs.c
DEP_CPP_FILSU=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\call.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\mydirent.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	"\acedb\wh\mystdlib.h"\
	".\..\wh\ibmdir.h"\
	
NODEP_CPP_FILSU=\
	".\..\wh\aixfs\dir.h"\
	".\..\wh\jfs\dir.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\arraysub.c
DEP_CPP_ARRAY=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\help.c
DEP_CPP_HELP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	".\..\wgd\gd.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\help.sbr" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\help.sbr" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\messubs.c
DEP_CPP_MESSU=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\messubs.obj" : $(SOURCE) $(DEP_CPP_MESSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messubs.sbr" : $(SOURCE) $(DEP_CPP_MESSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\messubs.obj" : $(SOURCE) $(DEP_CPP_MESSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messubs.sbr" : $(SOURCE) $(DEP_CPP_MESSU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\call.c
DEP_CPP_CALL_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\heap.c
DEP_CPP_HEAP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\heap.obj" : $(SOURCE) $(DEP_CPP_HEAP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\heap.sbr" : $(SOURCE) $(DEP_CPP_HEAP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\heap.obj" : $(SOURCE) $(DEP_CPP_HEAP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\heap.sbr" : $(SOURCE) $(DEP_CPP_HEAP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\freesubs.c
DEP_CPP_FREES=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\randsubs.c

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\display.c
DEP_CPP_DISPL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	
NODEP_CPP_DISPL=\
	".\..\w1\MPWIncludes.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\display.obj" : $(SOURCE) $(DEP_CPP_DISPL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\display.sbr" : $(SOURCE) $(DEP_CPP_DISPL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\display.obj" : $(SOURCE) $(DEP_CPP_DISPL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\display.sbr" : $(SOURCE) $(DEP_CPP_DISPL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\parse.c
DEP_CPP_PARSE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bstree.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\parse.sbr" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\parse.sbr" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\newkey.c
DEP_CPP_NEWKE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\bstree.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\newkey.obj" : $(SOURCE) $(DEP_CPP_NEWKE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\newkey.sbr" : $(SOURCE) $(DEP_CPP_NEWKE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\newkey.obj" : $(SOURCE) $(DEP_CPP_NEWKE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\newkey.sbr" : $(SOURCE) $(DEP_CPP_NEWKE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\keyset.c
DEP_CPP_KEYSE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keyset.sbr" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keyset.sbr" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\chrono.c
DEP_CPP_CHRON=\
	"\acedb\wh\acedb.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\chrono.sbr" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\chrono.sbr" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\dotter.c
DEP_CPP_DOTTE=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\dotter_.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\blxview.h"\
	".\..\WH\dotter.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dotter.obj" : $(SOURCE) $(DEP_CPP_DOTTE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dotter.sbr" : $(SOURCE) $(DEP_CPP_DOTTE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dotter.obj" : $(SOURCE) $(DEP_CPP_DOTTE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dotter.sbr" : $(SOURCE) $(DEP_CPP_DOTTE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\freeout.c
DEP_CPP_FREEO=\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w2\graphcon.c
DEP_CPP_GRAPH=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\graph.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphcon.obj" : $(SOURCE) $(DEP_CPP_GRAPH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphcon.sbr" : $(SOURCE) $(DEP_CPP_GRAPH) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphcon.obj" : $(SOURCE) $(DEP_CPP_GRAPH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphcon.sbr" : $(SOURCE) $(DEP_CPP_GRAPH) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w2\graphsub.c
DEP_CPP_GRAPHS=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphsub.obj" : $(SOURCE) $(DEP_CPP_GRAPHS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphsub.sbr" : $(SOURCE) $(DEP_CPP_GRAPHS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphsub.obj" : $(SOURCE) $(DEP_CPP_GRAPHS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphsub.sbr" : $(SOURCE) $(DEP_CPP_GRAPHS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\mainpick.c
DEP_CPP_MAINP=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\mainpick.obj" : $(SOURCE) $(DEP_CPP_MAINP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mainpick.sbr" : $(SOURCE) $(DEP_CPP_MAINP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\mainpick.obj" : $(SOURCE) $(DEP_CPP_MAINP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mainpick.sbr" : $(SOURCE) $(DEP_CPP_MAINP) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\update.c
DEP_CPP_UPDAT=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\update.obj" : $(SOURCE) $(DEP_CPP_UPDAT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\update.sbr" : $(SOURCE) $(DEP_CPP_UPDAT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\update.obj" : $(SOURCE) $(DEP_CPP_UPDAT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\update.sbr" : $(SOURCE) $(DEP_CPP_UPDAT) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\picksubs.c
DEP_CPP_PICKS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\picksubs.sbr" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\picksubs.sbr" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\longtext.c
DEP_CPP_LONGT=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\longtext.sbr" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\longtext.sbr" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\dump.c
DEP_CPP_DUMP_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	
NODEP_CPP_DUMP_=\
	".\..\w4\MPWIncludes.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dump.sbr" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dump.sbr" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\status.c
DEP_CPP_STATU=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\status.obj" : $(SOURCE) $(DEP_CPP_STATU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\status.sbr" : $(SOURCE) $(DEP_CPP_STATU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\status.obj" : $(SOURCE) $(DEP_CPP_STATU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\status.sbr" : $(SOURCE) $(DEP_CPP_STATU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\linkdate.c

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\linkdate.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\linkdate.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\model.c
DEP_CPP_MODEL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk.h"\
	
NODEP_CPP_MODEL=\
	".\..\w4\newmodel.c"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\model.sbr" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\model.sbr" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\queryexe.c
DEP_CPP_QUERY=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\queryexe.sbr" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\queryexe.sbr" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapcontrol.c
DEP_CPP_FMAPC=\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\key.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fmapcontrol.obj" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapcontrol.sbr" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapcontrol.obj" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapcontrol.sbr" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapgene.c
DEP_CPP_FMAPG=\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\dna.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fmapgene.obj" : $(SOURCE) $(DEP_CPP_FMAPG) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapgene.sbr" : $(SOURCE) $(DEP_CPP_FMAPG) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapgene.obj" : $(SOURCE) $(DEP_CPP_FMAPG) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapgene.sbr" : $(SOURCE) $(DEP_CPP_FMAPG) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapfeatures.c
DEP_CPP_FMAPF=\
	"\acedb\wh\fmap.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\bump.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	".\..\WH\dotter.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fmapfeatures.obj" : $(SOURCE) $(DEP_CPP_FMAPF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapfeatures.sbr" : $(SOURCE) $(DEP_CPP_FMAPF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapfeatures.obj" : $(SOURCE) $(DEP_CPP_FMAPF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapfeatures.sbr" : $(SOURCE) $(DEP_CPP_FMAPF) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapsequence.c
DEP_CPP_FMAPS=\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fmapsequence.obj" : $(SOURCE) $(DEP_CPP_FMAPS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapsequence.sbr" : $(SOURCE) $(DEP_CPP_FMAPS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapsequence.obj" : $(SOURCE) $(DEP_CPP_FMAPS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapsequence.sbr" : $(SOURCE) $(DEP_CPP_FMAPS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmaplocuscol.c
DEP_CPP_GMAPL=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmaplocuscol.obj" : $(SOURCE) $(DEP_CPP_GMAPL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmaplocuscol.sbr" : $(SOURCE) $(DEP_CPP_GMAPL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmaplocuscol.obj" : $(SOURCE) $(DEP_CPP_GMAPL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmaplocuscol.sbr" : $(SOURCE) $(DEP_CPP_GMAPL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\mapcontrol.c
DEP_CPP_MAPCO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\mapcontrol.obj" : $(SOURCE) $(DEP_CPP_MAPCO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mapcontrol.sbr" : $(SOURCE) $(DEP_CPP_MAPCO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\mapcontrol.obj" : $(SOURCE) $(DEP_CPP_MAPCO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mapcontrol.sbr" : $(SOURCE) $(DEP_CPP_MAPCO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\metab.c
DEP_CPP_METAB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\metab.obj" : $(SOURCE) $(DEP_CPP_METAB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\metab.sbr" : $(SOURCE) $(DEP_CPP_METAB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\metab.obj" : $(SOURCE) $(DEP_CPP_METAB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\metab.sbr" : $(SOURCE) $(DEP_CPP_METAB) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\drawdisp.c
DEP_CPP_DRAWD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\menu.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\drawdisp.obj" : $(SOURCE) $(DEP_CPP_DRAWD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\drawdisp.sbr" : $(SOURCE) $(DEP_CPP_DRAWD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\drawdisp.obj" : $(SOURCE) $(DEP_CPP_DRAWD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\drawdisp.sbr" : $(SOURCE) $(DEP_CPP_DRAWD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapdisp.c
DEP_CPP_GMAPD=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapdisp.obj" : $(SOURCE) $(DEP_CPP_GMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapdisp.sbr" : $(SOURCE) $(DEP_CPP_GMAPD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapdisp.obj" : $(SOURCE) $(DEP_CPP_GMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapdisp.sbr" : $(SOURCE) $(DEP_CPP_GMAPD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\cmapdisp.c
DEP_CPP_CMAPD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\cmapdisp.obj" : $(SOURCE) $(DEP_CPP_CMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\cmapdisp.sbr" : $(SOURCE) $(DEP_CPP_CMAPD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\cmapdisp.obj" : $(SOURCE) $(DEP_CPP_CMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\cmapdisp.sbr" : $(SOURCE) $(DEP_CPP_CMAPD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\vmapdrag.c
DEP_CPP_VMAPD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\vmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\vmapdrag.obj" : $(SOURCE) $(DEP_CPP_VMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapdrag.sbr" : $(SOURCE) $(DEP_CPP_VMAPD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\vmapdrag.obj" : $(SOURCE) $(DEP_CPP_VMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapdrag.sbr" : $(SOURCE) $(DEP_CPP_VMAPD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapmarkercol.c
DEP_CPP_GMAPM=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapmarkercol.obj" : $(SOURCE) $(DEP_CPP_GMAPM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapmarkercol.sbr" : $(SOURCE) $(DEP_CPP_GMAPM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapmarkercol.obj" : $(SOURCE) $(DEP_CPP_GMAPM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapmarkercol.sbr" : $(SOURCE) $(DEP_CPP_GMAPM) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\method.c
DEP_CPP_METHO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\method.obj" : $(SOURCE) $(DEP_CPP_METHO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\method.sbr" : $(SOURCE) $(DEP_CPP_METHO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\method.obj" : $(SOURCE) $(DEP_CPP_METHO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\method.sbr" : $(SOURCE) $(DEP_CPP_METHO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapremarkcol.c
DEP_CPP_GMAPR=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapremarkcol.obj" : $(SOURCE) $(DEP_CPP_GMAPR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapremarkcol.sbr" : $(SOURCE) $(DEP_CPP_GMAPR) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapremarkcol.obj" : $(SOURCE) $(DEP_CPP_GMAPR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapremarkcol.sbr" : $(SOURCE) $(DEP_CPP_GMAPR) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\griddisp.c
DEP_CPP_GRIDD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\griddisp.obj" : $(SOURCE) $(DEP_CPP_GRIDD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\griddisp.sbr" : $(SOURCE) $(DEP_CPP_GRIDD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\griddisp.obj" : $(SOURCE) $(DEP_CPP_GRIDD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\griddisp.sbr" : $(SOURCE) $(DEP_CPP_GRIDD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapintervalcol.c
DEP_CPP_GMAPI=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapintervalcol.obj" : $(SOURCE) $(DEP_CPP_GMAPI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapintervalcol.sbr" : $(SOURCE) $(DEP_CPP_GMAPI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapintervalcol.obj" : $(SOURCE) $(DEP_CPP_GMAPI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapintervalcol.sbr" : $(SOURCE) $(DEP_CPP_GMAPI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\pmapconvert.c
DEP_CPP_PMAPC=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\pmap_.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	
NODEP_CPP_PMAPC=\
	".\..\WH\pmap_trace_.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pmapconvert.obj" : $(SOURCE) $(DEP_CPP_PMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pmapconvert.sbr" : $(SOURCE) $(DEP_CPP_PMAPC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pmapconvert.obj" : $(SOURCE) $(DEP_CPP_PMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pmapconvert.sbr" : $(SOURCE) $(DEP_CPP_PMAPC) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\dnacpt.c
DEP_CPP_DNACP=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dnacpt.obj" : $(SOURCE) $(DEP_CPP_DNACP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnacpt.sbr" : $(SOURCE) $(DEP_CPP_DNACP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dnacpt.obj" : $(SOURCE) $(DEP_CPP_DNACP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnacpt.sbr" : $(SOURCE) $(DEP_CPP_DNACP) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapconvert.c
DEP_CPP_GMAPC=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapconvert.obj" : $(SOURCE) $(DEP_CPP_GMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapconvert.sbr" : $(SOURCE) $(DEP_CPP_GMAPC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapconvert.obj" : $(SOURCE) $(DEP_CPP_GMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapconvert.sbr" : $(SOURCE) $(DEP_CPP_GMAPC) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fpdisp.c
DEP_CPP_FPDIS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\fingerp.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fpdisp.obj" : $(SOURCE) $(DEP_CPP_FPDIS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fpdisp.sbr" : $(SOURCE) $(DEP_CPP_FPDIS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fpdisp.obj" : $(SOURCE) $(DEP_CPP_FPDIS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fpdisp.sbr" : $(SOURCE) $(DEP_CPP_FPDIS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\vmapdisp.c
DEP_CPP_VMAPDI=\
	"\acedb\wh\vmap.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\vmapdisp.obj" : $(SOURCE) $(DEP_CPP_VMAPDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapdisp.sbr" : $(SOURCE) $(DEP_CPP_VMAPDI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\vmapdisp.obj" : $(SOURCE) $(DEP_CPP_VMAPDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapdisp.sbr" : $(SOURCE) $(DEP_CPP_VMAPDI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\pmapdisp.c
DEP_CPP_PMAPD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\pmap_.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	
NODEP_CPP_PMAPD=\
	".\..\WH\pmap_trace_.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pmapdisp.obj" : $(SOURCE) $(DEP_CPP_PMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pmapdisp.sbr" : $(SOURCE) $(DEP_CPP_PMAPD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pmapdisp.obj" : $(SOURCE) $(DEP_CPP_PMAPD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pmapdisp.sbr" : $(SOURCE) $(DEP_CPP_PMAPD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\geldisp.c
DEP_CPP_GELDI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\geldisp.obj" : $(SOURCE) $(DEP_CPP_GELDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\geldisp.sbr" : $(SOURCE) $(DEP_CPP_GELDI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\geldisp.obj" : $(SOURCE) $(DEP_CPP_GELDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\geldisp.sbr" : $(SOURCE) $(DEP_CPP_GELDI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\biblio.c
DEP_CPP_BIBLI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\biblio.sbr" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\biblio.sbr" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapposnegcol.c
DEP_CPP_GMAPP=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapposnegcol.obj" : $(SOURCE) $(DEP_CPP_GMAPP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapposnegcol.sbr" : $(SOURCE) $(DEP_CPP_GMAPP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapposnegcol.obj" : $(SOURCE) $(DEP_CPP_GMAPP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapposnegcol.sbr" : $(SOURCE) $(DEP_CPP_GMAPP) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\tags.c
DEP_CPP_TAGS_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tags.sbr" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tags.sbr" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\sysclass.c
DEP_CPP_SYSCL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sysclass.sbr" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sysclass.sbr" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\quovadis.c
DEP_CPP_QUOVA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\menu.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\syn.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\quovadis.sbr" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\quovadis.sbr" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\class.c
DEP_CPP_CLASS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\menu.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\class.sbr" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\class.sbr" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\disptype.h

!IF  "$(CFG)" == "WinAce - Win32 Release"

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\classes.h

!IF  "$(CFG)" == "WinAce - Win32 Release"

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\tags.h

!IF  "$(CFG)" == "WinAce - Win32 Release"

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\sysclass.h

!IF  "$(CFG)" == "WinAce - Win32 Release"

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\systags.h

!IF  "$(CFG)" == "WinAce - Win32 Release"

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprdop.c
DEP_CPP_SPRDO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdop.sbr" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdop.sbr" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprdctrl.c
DEP_CPP_SPRDC=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\sprdctrl.obj" : $(SOURCE) $(DEP_CPP_SPRDC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdctrl.sbr" : $(SOURCE) $(DEP_CPP_SPRDC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprdctrl.obj" : $(SOURCE) $(DEP_CPP_SPRDC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdctrl.sbr" : $(SOURCE) $(DEP_CPP_SPRDC) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\ksetdisp.c
DEP_CPP_KSETD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dump.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\ksetdisp.obj" : $(SOURCE) $(DEP_CPP_KSETD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\ksetdisp.sbr" : $(SOURCE) $(DEP_CPP_KSETD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\ksetdisp.obj" : $(SOURCE) $(DEP_CPP_KSETD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\ksetdisp.sbr" : $(SOURCE) $(DEP_CPP_KSETD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprdmap.c
DEP_CPP_SPRDM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\sprdmap.obj" : $(SOURCE) $(DEP_CPP_SPRDM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdmap.sbr" : $(SOURCE) $(DEP_CPP_SPRDM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprdmap.obj" : $(SOURCE) $(DEP_CPP_SPRDM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdmap.sbr" : $(SOURCE) $(DEP_CPP_SPRDM) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\qbedisp.c
DEP_CPP_QBEDI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\query_.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\key.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\qbedisp.obj" : $(SOURCE) $(DEP_CPP_QBEDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\qbedisp.sbr" : $(SOURCE) $(DEP_CPP_QBEDI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\qbedisp.obj" : $(SOURCE) $(DEP_CPP_QBEDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\qbedisp.sbr" : $(SOURCE) $(DEP_CPP_QBEDI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprddisplay.c
DEP_CPP_SPRDD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\sprddisplay.obj" : $(SOURCE) $(DEP_CPP_SPRDD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddisplay.sbr" : $(SOURCE) $(DEP_CPP_SPRDD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprddisplay.obj" : $(SOURCE) $(DEP_CPP_SPRDD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddisplay.sbr" : $(SOURCE) $(DEP_CPP_SPRDD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\querybuild.c
DEP_CPP_QUERYB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\query_.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\menu.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\querybuild.obj" : $(SOURCE) $(DEP_CPP_QUERYB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\querybuild.sbr" : $(SOURCE) $(DEP_CPP_QUERYB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\querybuild.obj" : $(SOURCE) $(DEP_CPP_QUERYB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\querybuild.sbr" : $(SOURCE) $(DEP_CPP_QUERYB) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprddata.c
DEP_CPP_SPRDDA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\sprddata.obj" : $(SOURCE) $(DEP_CPP_SPRDDA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddata.sbr" : $(SOURCE) $(DEP_CPP_SPRDDA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprddata.obj" : $(SOURCE) $(DEP_CPP_SPRDDA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddata.sbr" : $(SOURCE) $(DEP_CPP_SPRDDA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\treedisp.c
DEP_CPP_TREED=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\treedisp.obj" : $(SOURCE) $(DEP_CPP_TREED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\treedisp.sbr" : $(SOURCE) $(DEP_CPP_TREED) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\treedisp.obj" : $(SOURCE) $(DEP_CPP_TREED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\treedisp.sbr" : $(SOURCE) $(DEP_CPP_TREED) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bsubs.c
DEP_CPP_BSUBS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\b_.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsubs.sbr" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsubs.sbr" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bsdumps.c
DEP_CPP_BSDUM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsdumps.sbr" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsdumps.sbr" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\check.c
DEP_CPP_CHECK=\
	"\acedb\wh\acedb.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\check.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\check.sbr" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\check.sbr" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\peptide.c
DEP_CPP_PEPTI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\peptide.sbr" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\peptide.sbr" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bstools.c
DEP_CPP_BSTOO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bstree.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstools.sbr" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstools.sbr" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bstree.c
DEP_CPP_BSTRE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\b_.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\key.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstree.sbr" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstree.sbr" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\nicedump.c
DEP_CPP_NICED=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nicedump.sbr" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nicedump.sbr" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\dnasubs.c
DEP_CPP_DNASU=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnasubs.sbr" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnasubs.sbr" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\asubs.c
DEP_CPP_ASUBS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\a_.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\asubs.sbr" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\asubs.sbr" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bssubs.c
DEP_CPP_BSSUB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\dump.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\check.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bssubs.sbr" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bssubs.sbr" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\action.c
DEP_CPP_ACTIO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\action.obj" : $(SOURCE) $(DEP_CPP_ACTIO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\action.sbr" : $(SOURCE) $(DEP_CPP_ACTIO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\action.obj" : $(SOURCE) $(DEP_CPP_ACTIO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\action.sbr" : $(SOURCE) $(DEP_CPP_ACTIO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\querydisp.c
DEP_CPP_QUERYD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\query_.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\igdevent.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\pick.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\querydisp.obj" : $(SOURCE) $(DEP_CPP_QUERYD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\querydisp.sbr" : $(SOURCE) $(DEP_CPP_QUERYD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\querydisp.obj" : $(SOURCE) $(DEP_CPP_QUERYD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\querydisp.sbr" : $(SOURCE) $(DEP_CPP_QUERYD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\blocksub.c
DEP_CPP_BLOCK=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk__.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blocksub.sbr" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blocksub.sbr" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\objcache.c
DEP_CPP_OBJCA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\cache_.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\objcache.sbr" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\objcache.sbr" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\disknew.c
DEP_CPP_DISKN=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\byteswap.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk_.h"\
	
NODEP_CPP_DISKN=\
	".\..\w5\macfilemgrsubs.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\disknew.sbr" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\disknew.sbr" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\lexsubs4.c
DEP_CPP_LEXSU=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs4.sbr" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs4.sbr" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\keysetdump.c
DEP_CPP_KEYSET=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keysetdump.sbr" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keysetdump.sbr" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\lexalpha.c
DEP_CPP_LEXAL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexalpha.sbr" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexalpha.sbr" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\lexsubs.c
DEP_CPP_LEXSUB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs.sbr" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs.sbr" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gmapdata.c
DEP_CPP_GMAPDA=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapdata.obj" : $(SOURCE) $(DEP_CPP_GMAPDA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapdata.sbr" : $(SOURCE) $(DEP_CPP_GMAPDA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapdata.obj" : $(SOURCE) $(DEP_CPP_GMAPDA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapdata.sbr" : $(SOURCE) $(DEP_CPP_GMAPDA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\embl.c
DEP_CPP_EMBL_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\embl.obj" : $(SOURCE) $(DEP_CPP_EMBL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\embl.sbr" : $(SOURCE) $(DEP_CPP_EMBL_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\embl.obj" : $(SOURCE) $(DEP_CPP_EMBL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\embl.sbr" : $(SOURCE) $(DEP_CPP_EMBL_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gfcode.c
DEP_CPP_GFCOD=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gfcode.obj" : $(SOURCE) $(DEP_CPP_GFCOD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gfcode.sbr" : $(SOURCE) $(DEP_CPP_GFCOD) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gfcode.obj" : $(SOURCE) $(DEP_CPP_GFCOD) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gfcode.sbr" : $(SOURCE) $(DEP_CPP_GFCOD) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\vmapphys.c
DEP_CPP_VMAPP=\
	"\acedb\wh\vmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\vmapphys.obj" : $(SOURCE) $(DEP_CPP_VMAPP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapphys.sbr" : $(SOURCE) $(DEP_CPP_VMAPP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\vmapphys.obj" : $(SOURCE) $(DEP_CPP_VMAPP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapphys.sbr" : $(SOURCE) $(DEP_CPP_VMAPP) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\vmapdata2.c
DEP_CPP_VMAPDA=\
	"\acedb\wh\vmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\vmapdata2.obj" : $(SOURCE) $(DEP_CPP_VMAPDA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapdata2.sbr" : $(SOURCE) $(DEP_CPP_VMAPDA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\vmapdata2.obj" : $(SOURCE) $(DEP_CPP_VMAPDA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vmapdata2.sbr" : $(SOURCE) $(DEP_CPP_VMAPDA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\align.c
DEP_CPP_ALIGN=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\align.obj" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\align.sbr" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\align.obj" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\align.sbr" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\translate.c
DEP_CPP_TRANS=\
	"\acedb\wh\iupac.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\translate.obj" : $(SOURCE) $(DEP_CPP_TRANS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\translate.sbr" : $(SOURCE) $(DEP_CPP_TRANS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\translate.obj" : $(SOURCE) $(DEP_CPP_TRANS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\translate.sbr" : $(SOURCE) $(DEP_CPP_TRANS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gmapphys.c
DEP_CPP_GMAPPH=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapphys.obj" : $(SOURCE) $(DEP_CPP_GMAPPH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapphys.sbr" : $(SOURCE) $(DEP_CPP_GMAPPH) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapphys.obj" : $(SOURCE) $(DEP_CPP_GMAPPH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapphys.sbr" : $(SOURCE) $(DEP_CPP_GMAPPH) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gmapdatacol.c
DEP_CPP_GMAPDAT=\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gmapdatacol.obj" : $(SOURCE) $(DEP_CPP_GMAPDAT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapdatacol.sbr" : $(SOURCE) $(DEP_CPP_GMAPDAT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gmapdatacol.obj" : $(SOURCE) $(DEP_CPP_GMAPDAT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gmapdatacol.sbr" : $(SOURCE) $(DEP_CPP_GMAPDAT) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\blxview.c
DEP_CPP_BLXVI=\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WH\dotter.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\acedb.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\blxview.obj" : $(SOURCE) $(DEP_CPP_BLXVI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blxview.sbr" : $(SOURCE) $(DEP_CPP_BLXVI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\blxview.obj" : $(SOURCE) $(DEP_CPP_BLXVI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blxview.sbr" : $(SOURCE) $(DEP_CPP_BLXVI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winace.cpp
DEP_CPP_WINAC=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\array.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\graph.h"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\fontpreference.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\winace.obj" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winace.sbr" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winace.rc

!IF  "$(CFG)" == "WinAce - Win32 Release"


"$(INTDIR)\winace.res" : $(SOURCE) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" $(SOURCE)


!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"


"$(INTDIR)\winace.res" : $(SOURCE) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\stdafx.cpp
DEP_CPP_STDAF=\
	"\acedb\win32h\stdafx.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\mainfrm.cpp
DEP_CPP_MAINF=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\mainfrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mainfrm.sbr" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphwin32_.cpp
DEP_CPP_GRAPHW=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\graphwin32_.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\graphwin32_.obj" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphwin32_.sbr" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\afxtrace.cpp
DEP_CPP_AFXTR=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fp"$(INTDIR)/winace.pch" /YX\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\afxtrace.obj" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\afxtrace.sbr" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\afxtrace.obj" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\afxtrace.sbr" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32print.cpp
DEP_CPP_WIN32=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\win32print.obj" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32print.sbr" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\msgwindow.cpp
DEP_CPP_MSGWI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\msgwindow.obj" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\msgwindow.sbr" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32.c
DEP_CPP_WIN32_=\
	"\acedb\win32h\win32.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32.obj" : $(SOURCE) $(DEP_CPP_WIN32_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32.sbr" : $(SOURCE) $(DEP_CPP_WIN32_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\win32.obj" : $(SOURCE) $(DEP_CPP_WIN32_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32.sbr" : $(SOURCE) $(DEP_CPP_WIN32_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\cmultiass.cpp
DEP_CPP_CMULT=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\cmultiass.obj" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\cmultiass.sbr" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32process.cpp
DEP_CPP_WIN32P=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\win32process.obj" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32process.sbr" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32lib.cpp
DEP_CPP_WIN32L=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\fontpreference.h"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\win32lib.obj" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32lib.sbr" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wgd\gd.c
DEP_CPP_GD_Cd4=\
	".\..\wgd\gd.h"\
	".\..\wgd\mtables.c"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\gd.obj" : $(SOURCE) $(DEP_CPP_GD_Cd4) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gd.sbr" : $(SOURCE) $(DEP_CPP_GD_Cd4) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gd.obj" : $(SOURCE) $(DEP_CPP_GD_Cd4) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gd.sbr" : $(SOURCE) $(DEP_CPP_GD_Cd4) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wstaden\opp.c
DEP_CPP_OPP_C=\
	"\acedb\wh\opp.h"\
	"\acedb\wh\seq.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WH\mach-io.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\opp.obj" : $(SOURCE) $(DEP_CPP_OPP_C) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\opp.sbr" : $(SOURCE) $(DEP_CPP_OPP_C) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\opp.obj" : $(SOURCE) $(DEP_CPP_OPP_C) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\opp.sbr" : $(SOURCE) $(DEP_CPP_OPP_C) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\mapscrollview.cpp
DEP_CPP_MAPSC=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\mapscrollview.obj" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mapscrollview.sbr" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\textfitview.cpp
DEP_CPP_TEXTF=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	"\acedb\win32h\FitView.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\textfitview.obj" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\textfitview.sbr" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\textscrollview.cpp
DEP_CPP_TEXTS=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\textscrollview.obj" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\textscrollview.sbr" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\pixelfitview.cpp
DEP_CPP_PIXEL=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	"\acedb\win32h\FitView.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\pixelfitview.obj" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pixelfitview.sbr" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\pixelscrollview.cpp
DEP_CPP_PIXELS=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\pixelscrollview.obj" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pixelscrollview.sbr" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\acefiledoc.cpp
DEP_CPP_ACEFI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\acefiledoc.obj" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acefiledoc.sbr" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\acefileview.cpp
DEP_CPP_ACEFIL=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\acefileview.obj" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acefileview.sbr" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\fullscrollview.cpp
DEP_CPP_FULLS=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\fullscrollview.obj" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fullscrollview.sbr" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\cgraph.cpp
DEP_CPP_CGRAP=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\cgraph.obj" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\cgraph.sbr" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphview.cpp
DEP_CPP_GRAPHV=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\key.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\graphview.obj" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphview.sbr" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphwin.cpp
DEP_CPP_GRAPHWI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\graphwin32_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\graphwin.obj" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphwin.sbr" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32menusubs.cpp
DEP_CPP_WIN32M=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\menu_.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\win32menusubs.obj" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32menusubs.sbr" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphselect.cpp
DEP_CPP_GRAPHSE=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\graphselect.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\graphselect.obj" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphselect.sbr" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\windialogs.cpp
DEP_CPP_WINDI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\windialogs.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\windialogs.obj" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\windialogs.sbr" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winfilquery.cpp
DEP_CPP_WINFI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\winfilquery.obj" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winfilquery.sbr" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W6\alignment.c
DEP_CPP_ALIGNM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\alignment.sbr" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\alignment.sbr" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W9\dotterKarlin.c
DEP_CPP_DOTTER=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dotterKarlin.obj" : $(SOURCE) $(DEP_CPP_DOTTER) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dotterKarlin.sbr" : $(SOURCE) $(DEP_CPP_DOTTER) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dotterKarlin.obj" : $(SOURCE) $(DEP_CPP_DOTTER) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dotterKarlin.sbr" : $(SOURCE) $(DEP_CPP_DOTTER) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W9\hexcode.c
DEP_CPP_HEXCO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\method.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\hexcode.obj" : $(SOURCE) $(DEP_CPP_HEXCO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\hexcode.sbr" : $(SOURCE) $(DEP_CPP_HEXCO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\hexcode.obj" : $(SOURCE) $(DEP_CPP_HEXCO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\hexcode.sbr" : $(SOURCE) $(DEP_CPP_HEXCO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\caceprintpage.cpp
DEP_CPP_CACEP=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\caceprintpage.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\caceprintpage.obj" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\caceprintpage.sbr" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\preferences.cpp
DEP_CPP_PREFE=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\fontpreference.h"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\preferences.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\preferences.sbr" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\dbprofile.cpp
DEP_CPP_DBPRO=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\dbprofile.obj" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dbprofile.sbr" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\userinterface.cpp
DEP_CPP_USERI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\userinterface.obj" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\userinterface.sbr" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\acedbprofileview.cpp
DEP_CPP_ACEDB=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\acedbprofileview.obj" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acedbprofileview.sbr" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\acedbprofile.cpp
DEP_CPP_ACEDBP=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\acedbprofile.obj" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acedbprofile.sbr" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\fontpreference.cpp
DEP_CPP_FONTP=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\fontpreference.obj" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fontpreference.sbr" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\WIN32\splashbox.cpp
DEP_CPP_SPLAS=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\splashbox.obj" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\splashbox.sbr" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pepseqcol.c
DEP_CPP_PEPSE=\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pepseqcol.obj" : $(SOURCE) $(DEP_CPP_PEPSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepseqcol.sbr" : $(SOURCE) $(DEP_CPP_PEPSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepseqcol.obj" : $(SOURCE) $(DEP_CPP_PEPSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepseqcol.sbr" : $(SOURCE) $(DEP_CPP_PEPSE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pephomolcol.c
DEP_CPP_PEPHO=\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pephomolcol.obj" : $(SOURCE) $(DEP_CPP_PEPHO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pephomolcol.sbr" : $(SOURCE) $(DEP_CPP_PEPHO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pephomolcol.obj" : $(SOURCE) $(DEP_CPP_PEPHO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pephomolcol.sbr" : $(SOURCE) $(DEP_CPP_PEPHO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pepgraphcol.c
DEP_CPP_PEPGR=\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pepgraphcol.obj" : $(SOURCE) $(DEP_CPP_PEPGR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepgraphcol.sbr" : $(SOURCE) $(DEP_CPP_PEPGR) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepgraphcol.obj" : $(SOURCE) $(DEP_CPP_PEPGR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepgraphcol.sbr" : $(SOURCE) $(DEP_CPP_PEPGR) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pepdisp.c
DEP_CPP_PEPDI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pepdisp.obj" : $(SOURCE) $(DEP_CPP_PEPDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepdisp.sbr" : $(SOURCE) $(DEP_CPP_PEPDI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepdisp.obj" : $(SOURCE) $(DEP_CPP_PEPDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepdisp.sbr" : $(SOURCE) $(DEP_CPP_PEPDI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W2\colcontrol.c
DEP_CPP_COLCO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\colcontrol_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\colcontrol.obj" : $(SOURCE) $(DEP_CPP_COLCO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\colcontrol.sbr" : $(SOURCE) $(DEP_CPP_COLCO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\colcontrol.obj" : $(SOURCE) $(DEP_CPP_COLCO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\colcontrol.sbr" : $(SOURCE) $(DEP_CPP_COLCO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W2\viewedit.c
DEP_CPP_VIEWE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\bs.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\colcontrol_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\viewedit.obj" : $(SOURCE) $(DEP_CPP_VIEWE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\viewedit.sbr" : $(SOURCE) $(DEP_CPP_VIEWE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\viewedit.obj" : $(SOURCE) $(DEP_CPP_VIEWE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\viewedit.sbr" : $(SOURCE) $(DEP_CPP_VIEWE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\acedialogs.c
DEP_CPP_ACEDI=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\acedialogs.obj" : $(SOURCE) $(DEP_CPP_ACEDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acedialogs.sbr" : $(SOURCE) $(DEP_CPP_ACEDI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\acedialogs.obj" : $(SOURCE) $(DEP_CPP_ACEDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acedialogs.sbr" : $(SOURCE) $(DEP_CPP_ACEDI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winmain.c

!IF  "$(CFG)" == "WinAce - Win32 Release"

DEP_CPP_WINMA=\
	"\acedb\win32h\win32.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\menu.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\winmain.obj" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winmain.sbr" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

DEP_CPP_WINMA=\
	"\acedb\win32h\win32.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\menu.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\winmain.obj" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winmain.sbr" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W4\banner.c
DEP_CPP_BANNE=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\array.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\banner.sbr" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\banner.sbr" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphbox.cpp
DEP_CPP_GRAPHB=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\graphbox.obj" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphbox.sbr" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\plainview.cpp
DEP_CPP_PLAIN=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\plainview.obj" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\plainview.sbr" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W1\oldhelp.c
DEP_CPP_OLDHE=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oldhelp.sbr" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oldhelp.sbr" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W4\session.c
DEP_CPP_SESSI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\whooks\disptype.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\session.sbr" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\session.sbr" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\colourbox.c
DEP_CPP_COLOU=\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\colourbox.obj" : $(SOURCE) $(DEP_CPP_COLOU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\colourbox.sbr" : $(SOURCE) $(DEP_CPP_COLOU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\colourbox.obj" : $(SOURCE) $(DEP_CPP_COLOU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\colourbox.sbr" : $(SOURCE) $(DEP_CPP_COLOU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W2\graphramp.c

!IF  "$(CFG)" == "WinAce - Win32 Release"

DEP_CPP_GRAPHR=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\graph.h"\
	
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphramp.obj" : $(SOURCE) $(DEP_CPP_GRAPHR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphramp.sbr" : $(SOURCE) $(DEP_CPP_GRAPHR) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

DEP_CPP_GRAPHR=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\graph.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphramp.obj" : $(SOURCE) $(DEP_CPP_GRAPHR) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphramp.sbr" : $(SOURCE) $(DEP_CPP_GRAPHR) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W6\plot.c
DEP_CPP_PLOT_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "WinAce - Win32 Release"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "ACEDB" /D\
 "ACEDB4" /D "_MBCS" /Fr"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\plot.obj" : $(SOURCE) $(DEP_CPP_PLOT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\plot.sbr" : $(SOURCE) $(DEP_CPP_PLOT_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "WinAce - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\plot.obj" : $(SOURCE) $(DEP_CPP_PLOT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\plot.sbr" : $(SOURCE) $(DEP_CPP_PLOT_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\WIN32\fitview.cpp
DEP_CPP_FITVI=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\wh\graph_.h"\
	".\..\WIN32H\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\WIN32H\CMultiAss.h"\
	".\..\WIN32H\win32menus.h"\
	".\..\WIN32H\graphbox.h"\
	".\..\WIN32H\msgwindow.h"\
	".\..\WIN32H\graphview.h"\
	".\..\WIN32H\graphwin.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\fitview.obj" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fitview.sbr" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"
   $(BuildCmds)

# End Source File
# End Target
# End Project
################################################################################
 
 

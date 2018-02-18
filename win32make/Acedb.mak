# $Id: Acedb.mak,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $
# Microsoft Developer Studio Generated NMAKE File, Format Version 4.20
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

!IF "$(CFG)" == ""
CFG=acedb - Win32 WinTace debug
!MESSAGE No configuration specified.  Defaulting to acedb - Win32 WinTace\
 debug.
!ENDIF 

!IF "$(CFG)" != "acedb - Win32 Release" && "$(CFG)" != "acedb - Win32 Debug" &&\
 "$(CFG)" != "acedb - Win32 Acembly Debug" && "$(CFG)" !=\
 "acedb - Win32 Acembly Release" && "$(CFG)" != "acedb - Win32 AceClient Debug"\
 && "$(CFG)" != "acedb - Win32 AceServer Debug" && "$(CFG)" !=\
 "acedb - Win32 AceServer Release" && "$(CFG)" !=\
 "acedb - Win32 AceClient Release" && "$(CFG)" !=\
 "acedb - Win32 NetClient Debug" && "$(CFG)" !=\
 "acedb - Win32 NetClient Release" && "$(CFG)" != "acedb - Win32 WinTace debug"\
 && "$(CFG)" != "acedb - Win32 WinTace Release"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "Acedb.mak" CFG="acedb - Win32 WinTace debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "acedb - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 Acembly Debug" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 Acembly Release" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 AceClient Debug" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 AceServer Debug" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 AceServer Release" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 AceClient Release" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 NetClient Debug" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 NetClient Release" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 WinTace debug" (based on "Win32 (x86) Application")
!MESSAGE "acedb - Win32 WinTace Release" (based on "Win32 (x86) Application")
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
# PROP Target_Last_Scanned "acedb - Win32 Debug"
MTL=mktyplib.exe
RSC=rc.exe
CPP=cl.exe

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "winrel"
# PROP BASE Intermediate_Dir "winrel"
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "\acedb\bin\winrel"
# PROP Intermediate_Dir "\acedb\bin\winrel"
OUTDIR=\acedb\bin\winrel
INTDIR=\acedb\bin\winrel

ALL : "$(OUTDIR)\winace.exe"

CLEAN : 
	-@erase "$(INTDIR)\acedbprofile.obj"
	-@erase "$(INTDIR)\acedbprofileview.obj"
	-@erase "$(INTDIR)\acedialogs.obj"
	-@erase "$(INTDIR)\acefiledoc.obj"
	-@erase "$(INTDIR)\acefileview.obj"
	-@erase "$(INTDIR)\action.obj"
	-@erase "$(INTDIR)\afxtrace.obj"
	-@erase "$(INTDIR)\align.obj"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\blxview.obj"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\caceprintpage.obj"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\cgraph.obj"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\cmapdisp.obj"
	-@erase "$(INTDIR)\cmultiass.obj"
	-@erase "$(INTDIR)\colcontrol.obj"
	-@erase "$(INTDIR)\colourbox.obj"
	-@erase "$(INTDIR)\dbprofile.obj"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\display.obj"
	-@erase "$(INTDIR)\dnacpt.obj"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dotter.obj"
	-@erase "$(INTDIR)\dotterKarlin.obj"
	-@erase "$(INTDIR)\drawdisp.obj"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\embl.obj"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\fitview.obj"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\fmapblast.obj"
	-@erase "$(INTDIR)\fmapcontrol.obj"
	-@erase "$(INTDIR)\fmapfeatures.obj"
	-@erase "$(INTDIR)\fmapgene.obj"
	-@erase "$(INTDIR)\fmapmenes.obj"
	-@erase "$(INTDIR)\fmaposp.obj"
	-@erase "$(INTDIR)\fmapsequence.obj"
	-@erase "$(INTDIR)\fontpreference.obj"
	-@erase "$(INTDIR)\forest.obj"
	-@erase "$(INTDIR)\forestdisplay.obj"
	-@erase "$(INTDIR)\fpdisp.obj"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\fullscrollview.obj"
	-@erase "$(INTDIR)\gd.obj"
	-@erase "$(INTDIR)\geldisp.obj"
	-@erase "$(INTDIR)\gfcode.obj"
	-@erase "$(INTDIR)\gmapconvert.obj"
	-@erase "$(INTDIR)\gmapdata.obj"
	-@erase "$(INTDIR)\gmapdatacol.obj"
	-@erase "$(INTDIR)\gmapdisp.obj"
	-@erase "$(INTDIR)\gmapintervalcol.obj"
	-@erase "$(INTDIR)\gmaplocuscol.obj"
	-@erase "$(INTDIR)\gmapmarkercol.obj"
	-@erase "$(INTDIR)\gmapphys.obj"
	-@erase "$(INTDIR)\gmapposnegcol.obj"
	-@erase "$(INTDIR)\gmapremarkcol.obj"
	-@erase "$(INTDIR)\graphbox.obj"
	-@erase "$(INTDIR)\graphcon.obj"
	-@erase "$(INTDIR)\graphramp.obj"
	-@erase "$(INTDIR)\graphselect.obj"
	-@erase "$(INTDIR)\graphsub.obj"
	-@erase "$(INTDIR)\graphview.obj"
	-@erase "$(INTDIR)\graphwin.obj"
	-@erase "$(INTDIR)\graphwin32_.obj"
	-@erase "$(INTDIR)\griddisp.obj"
	-@erase "$(INTDIR)\heap.obj"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\hexcode.obj"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\ksetdisp.obj"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\mainfrm.obj"
	-@erase "$(INTDIR)\mainpick.obj"
	-@erase "$(INTDIR)\mapcontrol.obj"
	-@erase "$(INTDIR)\mapscrollview.obj"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\menu.obj"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messubs.obj"
	-@erase "$(INTDIR)\metab.obj"
	-@erase "$(INTDIR)\method.obj"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\msgwindow.obj"
	-@erase "$(INTDIR)\newkey.obj"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\nqcpass1.obj"
	-@erase "$(INTDIR)\nqcpass2.obj"
	-@erase "$(INTDIR)\nqctools.obj"
	-@erase "$(INTDIR)\o2m.obj"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\oxgriddisp.obj"
	-@erase "$(INTDIR)\oxhomlist.obj"
	-@erase "$(INTDIR)\pairmapdisp.obj"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\pepactivezonecol.obj"
	-@erase "$(INTDIR)\pepdisp.obj"
	-@erase "$(INTDIR)\pepfeaturecol.obj"
	-@erase "$(INTDIR)\pepgraphcol.obj"
	-@erase "$(INTDIR)\pephomolcol.obj"
	-@erase "$(INTDIR)\pepseqcol.obj"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\pixelfitview.obj"
	-@erase "$(INTDIR)\pixelscrollview.obj"
	-@erase "$(INTDIR)\plainview.obj"
	-@erase "$(INTDIR)\plot.obj"
	-@erase "$(INTDIR)\pmapconvert.obj"
	-@erase "$(INTDIR)\pmapdisp.obj"
	-@erase "$(INTDIR)\preferences.obj"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\qbedisp.obj"
	-@erase "$(INTDIR)\querybuild.obj"
	-@erase "$(INTDIR)\querydisp.obj"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\specg.obj"
	-@erase "$(INTDIR)\splashbox.obj"
	-@erase "$(INTDIR)\sprdctrl.obj"
	-@erase "$(INTDIR)\sprddata.obj"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprddisplay.obj"
	-@erase "$(INTDIR)\sprdmap.obj"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\status.obj"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\textfitview.obj"
	-@erase "$(INTDIR)\textscrollview.obj"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\tqdebug.obj"
	-@erase "$(INTDIR)\translate.obj"
	-@erase "$(INTDIR)\treedisp.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\userinterface.obj"
	-@erase "$(INTDIR)\viewedit.obj"
	-@erase "$(INTDIR)\vmapdata2.obj"
	-@erase "$(INTDIR)\vmapdisp.obj"
	-@erase "$(INTDIR)\vmapdrag.obj"
	-@erase "$(INTDIR)\vmapphys.obj"
	-@erase "$(INTDIR)\win32.obj"
	-@erase "$(INTDIR)\win32lib.obj"
	-@erase "$(INTDIR)\win32menusubs.obj"
	-@erase "$(INTDIR)\win32print.obj"
	-@erase "$(INTDIR)\win32process.obj"
	-@erase "$(INTDIR)\win32thread.obj"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\winace.obj"
	-@erase "$(INTDIR)\winace.res"
	-@erase "$(INTDIR)\winacemail.obj"
	-@erase "$(INTDIR)\windialogs.obj"
	-@erase "$(INTDIR)\winfilquery.obj"
	-@erase "$(INTDIR)\winmain.obj"
	-@erase "$(OUTDIR)\winace.exe"
	-@erase "$(OUTDIR)\winace.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /YX /c
# SUBTRACT CPP /Fr
CPP_PROJ=/nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\winrel/
CPP_SBRS=.\.
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"\acedb\bin\winrel/winace.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/winace.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /map /machine:I386 /nodefaultlib:"libcmt.lib" /out:"\acedb\bin\winrel/winace.exe"
# SUBTRACT LINK32 /verbose /nodefaultlib
LINK32_FLAGS=/nologo /version:4.31 /entry:"WinMainCRTStartup"\
 /subsystem:windows /incremental:no /pdb:"$(OUTDIR)/winace.pdb"\
 /map:"$(INTDIR)/winace.map" /machine:I386 /nodefaultlib:"libcmt.lib"\
 /out:"$(OUTDIR)/winace.exe" 
LINK32_OBJS= \
	"$(INTDIR)\acedbprofile.obj" \
	"$(INTDIR)\acedbprofileview.obj" \
	"$(INTDIR)\acedialogs.obj" \
	"$(INTDIR)\acefiledoc.obj" \
	"$(INTDIR)\acefileview.obj" \
	"$(INTDIR)\action.obj" \
	"$(INTDIR)\afxtrace.obj" \
	"$(INTDIR)\align.obj" \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\blxview.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\caceprintpage.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\cgraph.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\cmapdisp.obj" \
	"$(INTDIR)\cmultiass.obj" \
	"$(INTDIR)\colcontrol.obj" \
	"$(INTDIR)\colourbox.obj" \
	"$(INTDIR)\dbprofile.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\display.obj" \
	"$(INTDIR)\dnacpt.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dotter.obj" \
	"$(INTDIR)\dotterKarlin.obj" \
	"$(INTDIR)\drawdisp.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\embl.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\fitview.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\fmapblast.obj" \
	"$(INTDIR)\fmapcontrol.obj" \
	"$(INTDIR)\fmapfeatures.obj" \
	"$(INTDIR)\fmapgene.obj" \
	"$(INTDIR)\fmapmenes.obj" \
	"$(INTDIR)\fmaposp.obj" \
	"$(INTDIR)\fmapsequence.obj" \
	"$(INTDIR)\fontpreference.obj" \
	"$(INTDIR)\forest.obj" \
	"$(INTDIR)\forestdisplay.obj" \
	"$(INTDIR)\fpdisp.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\fullscrollview.obj" \
	"$(INTDIR)\gd.obj" \
	"$(INTDIR)\geldisp.obj" \
	"$(INTDIR)\gfcode.obj" \
	"$(INTDIR)\gmapconvert.obj" \
	"$(INTDIR)\gmapdata.obj" \
	"$(INTDIR)\gmapdatacol.obj" \
	"$(INTDIR)\gmapdisp.obj" \
	"$(INTDIR)\gmapintervalcol.obj" \
	"$(INTDIR)\gmaplocuscol.obj" \
	"$(INTDIR)\gmapmarkercol.obj" \
	"$(INTDIR)\gmapphys.obj" \
	"$(INTDIR)\gmapposnegcol.obj" \
	"$(INTDIR)\gmapremarkcol.obj" \
	"$(INTDIR)\graphbox.obj" \
	"$(INTDIR)\graphcon.obj" \
	"$(INTDIR)\graphramp.obj" \
	"$(INTDIR)\graphselect.obj" \
	"$(INTDIR)\graphsub.obj" \
	"$(INTDIR)\graphview.obj" \
	"$(INTDIR)\graphwin.obj" \
	"$(INTDIR)\graphwin32_.obj" \
	"$(INTDIR)\griddisp.obj" \
	"$(INTDIR)\heap.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\hexcode.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\ksetdisp.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\mainfrm.obj" \
	"$(INTDIR)\mainpick.obj" \
	"$(INTDIR)\mapcontrol.obj" \
	"$(INTDIR)\mapscrollview.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\menu.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\messubs.obj" \
	"$(INTDIR)\metab.obj" \
	"$(INTDIR)\method.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\msgwindow.obj" \
	"$(INTDIR)\newkey.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\nqcpass1.obj" \
	"$(INTDIR)\nqcpass2.obj" \
	"$(INTDIR)\nqctools.obj" \
	"$(INTDIR)\o2m.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\oxgriddisp.obj" \
	"$(INTDIR)\oxhomlist.obj" \
	"$(INTDIR)\pairmapdisp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\pepactivezonecol.obj" \
	"$(INTDIR)\pepdisp.obj" \
	"$(INTDIR)\pepfeaturecol.obj" \
	"$(INTDIR)\pepgraphcol.obj" \
	"$(INTDIR)\pephomolcol.obj" \
	"$(INTDIR)\pepseqcol.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\pixelfitview.obj" \
	"$(INTDIR)\pixelscrollview.obj" \
	"$(INTDIR)\plainview.obj" \
	"$(INTDIR)\plot.obj" \
	"$(INTDIR)\pmapconvert.obj" \
	"$(INTDIR)\pmapdisp.obj" \
	"$(INTDIR)\preferences.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\qbedisp.obj" \
	"$(INTDIR)\querybuild.obj" \
	"$(INTDIR)\querydisp.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\specg.obj" \
	"$(INTDIR)\splashbox.obj" \
	"$(INTDIR)\sprdctrl.obj" \
	"$(INTDIR)\sprddata.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprddisplay.obj" \
	"$(INTDIR)\sprdmap.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\status.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\textfitview.obj" \
	"$(INTDIR)\textscrollview.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\tqdebug.obj" \
	"$(INTDIR)\translate.obj" \
	"$(INTDIR)\treedisp.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\userinterface.obj" \
	"$(INTDIR)\viewedit.obj" \
	"$(INTDIR)\vmapdata2.obj" \
	"$(INTDIR)\vmapdisp.obj" \
	"$(INTDIR)\vmapdrag.obj" \
	"$(INTDIR)\vmapphys.obj" \
	"$(INTDIR)\win32.obj" \
	"$(INTDIR)\win32lib.obj" \
	"$(INTDIR)\win32menusubs.obj" \
	"$(INTDIR)\win32print.obj" \
	"$(INTDIR)\win32process.obj" \
	"$(INTDIR)\win32thread.obj" \
	"$(INTDIR)\win32util.obj" \
	"$(INTDIR)\winace.obj" \
	"$(INTDIR)\winace.res" \
	"$(INTDIR)\winacemail.obj" \
	"$(INTDIR)\windialogs.obj" \
	"$(INTDIR)\winfilquery.obj" \
	"$(INTDIR)\winmain.obj"

"$(OUTDIR)\winace.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "windebug"
# PROP BASE Intermediate_Dir "windebug"
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "\acedb\bin\windebug"
# PROP Intermediate_Dir "\acedb\bin\windebug"
OUTDIR=\acedb\bin\windebug
INTDIR=\acedb\bin\windebug

ALL : "$(OUTDIR)\winaced.exe" "$(OUTDIR)\winaced.bsc"

CLEAN : 
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\acedbprofile.obj"
	-@erase "$(INTDIR)\acedbprofile.sbr"
	-@erase "$(INTDIR)\acedbprofileview.obj"
	-@erase "$(INTDIR)\acedbprofileview.sbr"
	-@erase "$(INTDIR)\acedialogs.obj"
	-@erase "$(INTDIR)\acedialogs.sbr"
	-@erase "$(INTDIR)\acefiledoc.obj"
	-@erase "$(INTDIR)\acefiledoc.sbr"
	-@erase "$(INTDIR)\acefileview.obj"
	-@erase "$(INTDIR)\acefileview.sbr"
	-@erase "$(INTDIR)\action.obj"
	-@erase "$(INTDIR)\action.sbr"
	-@erase "$(INTDIR)\afxtrace.obj"
	-@erase "$(INTDIR)\afxtrace.sbr"
	-@erase "$(INTDIR)\align.obj"
	-@erase "$(INTDIR)\align.sbr"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\alignment.sbr"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\arraysub.sbr"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\asubs.sbr"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\banner.sbr"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\biblio.sbr"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\blocksub.sbr"
	-@erase "$(INTDIR)\blxview.obj"
	-@erase "$(INTDIR)\blxview.sbr"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bsdumps.sbr"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bssubs.sbr"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstools.sbr"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bstree.sbr"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bsubs.sbr"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\bump.sbr"
	-@erase "$(INTDIR)\caceprintpage.obj"
	-@erase "$(INTDIR)\caceprintpage.sbr"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\call.sbr"
	-@erase "$(INTDIR)\cgraph.obj"
	-@erase "$(INTDIR)\cgraph.sbr"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\check.sbr"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\chrono.sbr"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\class.sbr"
	-@erase "$(INTDIR)\cmapdisp.obj"
	-@erase "$(INTDIR)\cmapdisp.sbr"
	-@erase "$(INTDIR)\cmultiass.obj"
	-@erase "$(INTDIR)\cmultiass.sbr"
	-@erase "$(INTDIR)\colcontrol.obj"
	-@erase "$(INTDIR)\colcontrol.sbr"
	-@erase "$(INTDIR)\colourbox.obj"
	-@erase "$(INTDIR)\colourbox.sbr"
	-@erase "$(INTDIR)\dbprofile.obj"
	-@erase "$(INTDIR)\dbprofile.sbr"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\dict.sbr"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\disknew.sbr"
	-@erase "$(INTDIR)\display.obj"
	-@erase "$(INTDIR)\display.sbr"
	-@erase "$(INTDIR)\dnacpt.obj"
	-@erase "$(INTDIR)\dnacpt.sbr"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dnasubs.sbr"
	-@erase "$(INTDIR)\dotter.obj"
	-@erase "$(INTDIR)\dotter.sbr"
	-@erase "$(INTDIR)\dotterKarlin.obj"
	-@erase "$(INTDIR)\dotterKarlin.sbr"
	-@erase "$(INTDIR)\drawdisp.obj"
	-@erase "$(INTDIR)\drawdisp.sbr"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\dump.sbr"
	-@erase "$(INTDIR)\embl.obj"
	-@erase "$(INTDIR)\embl.sbr"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\filsubs.sbr"
	-@erase "$(INTDIR)\fitview.obj"
	-@erase "$(INTDIR)\fitview.sbr"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\flag.sbr"
	-@erase "$(INTDIR)\fmapblast.obj"
	-@erase "$(INTDIR)\fmapblast.sbr"
	-@erase "$(INTDIR)\fmapcontrol.obj"
	-@erase "$(INTDIR)\fmapcontrol.sbr"
	-@erase "$(INTDIR)\fmapfeatures.obj"
	-@erase "$(INTDIR)\fmapfeatures.sbr"
	-@erase "$(INTDIR)\fmapgene.obj"
	-@erase "$(INTDIR)\fmapgene.sbr"
	-@erase "$(INTDIR)\fmapmenes.obj"
	-@erase "$(INTDIR)\fmapmenes.sbr"
	-@erase "$(INTDIR)\fmaposp.obj"
	-@erase "$(INTDIR)\fmaposp.sbr"
	-@erase "$(INTDIR)\fmapsequence.obj"
	-@erase "$(INTDIR)\fmapsequence.sbr"
	-@erase "$(INTDIR)\fontpreference.obj"
	-@erase "$(INTDIR)\fontpreference.sbr"
	-@erase "$(INTDIR)\forest.obj"
	-@erase "$(INTDIR)\forest.sbr"
	-@erase "$(INTDIR)\forestdisplay.obj"
	-@erase "$(INTDIR)\forestdisplay.sbr"
	-@erase "$(INTDIR)\fpdisp.obj"
	-@erase "$(INTDIR)\fpdisp.sbr"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freeout.sbr"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\freesubs.sbr"
	-@erase "$(INTDIR)\fullscrollview.obj"
	-@erase "$(INTDIR)\fullscrollview.sbr"
	-@erase "$(INTDIR)\gd.obj"
	-@erase "$(INTDIR)\gd.sbr"
	-@erase "$(INTDIR)\geldisp.obj"
	-@erase "$(INTDIR)\geldisp.sbr"
	-@erase "$(INTDIR)\gfcode.obj"
	-@erase "$(INTDIR)\gfcode.sbr"
	-@erase "$(INTDIR)\gmapconvert.obj"
	-@erase "$(INTDIR)\gmapconvert.sbr"
	-@erase "$(INTDIR)\gmapdata.obj"
	-@erase "$(INTDIR)\gmapdata.sbr"
	-@erase "$(INTDIR)\gmapdatacol.obj"
	-@erase "$(INTDIR)\gmapdatacol.sbr"
	-@erase "$(INTDIR)\gmapdisp.obj"
	-@erase "$(INTDIR)\gmapdisp.sbr"
	-@erase "$(INTDIR)\gmapintervalcol.obj"
	-@erase "$(INTDIR)\gmapintervalcol.sbr"
	-@erase "$(INTDIR)\gmaplocuscol.obj"
	-@erase "$(INTDIR)\gmaplocuscol.sbr"
	-@erase "$(INTDIR)\gmapmarkercol.obj"
	-@erase "$(INTDIR)\gmapmarkercol.sbr"
	-@erase "$(INTDIR)\gmapphys.obj"
	-@erase "$(INTDIR)\gmapphys.sbr"
	-@erase "$(INTDIR)\gmapposnegcol.obj"
	-@erase "$(INTDIR)\gmapposnegcol.sbr"
	-@erase "$(INTDIR)\gmapremarkcol.obj"
	-@erase "$(INTDIR)\gmapremarkcol.sbr"
	-@erase "$(INTDIR)\graphbox.obj"
	-@erase "$(INTDIR)\graphbox.sbr"
	-@erase "$(INTDIR)\graphcon.obj"
	-@erase "$(INTDIR)\graphcon.sbr"
	-@erase "$(INTDIR)\graphramp.obj"
	-@erase "$(INTDIR)\graphramp.sbr"
	-@erase "$(INTDIR)\graphselect.obj"
	-@erase "$(INTDIR)\graphselect.sbr"
	-@erase "$(INTDIR)\graphsub.obj"
	-@erase "$(INTDIR)\graphsub.sbr"
	-@erase "$(INTDIR)\graphview.obj"
	-@erase "$(INTDIR)\graphview.sbr"
	-@erase "$(INTDIR)\graphwin.obj"
	-@erase "$(INTDIR)\graphwin.sbr"
	-@erase "$(INTDIR)\graphwin32_.obj"
	-@erase "$(INTDIR)\graphwin32_.sbr"
	-@erase "$(INTDIR)\griddisp.obj"
	-@erase "$(INTDIR)\griddisp.sbr"
	-@erase "$(INTDIR)\heap.obj"
	-@erase "$(INTDIR)\heap.sbr"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\help.sbr"
	-@erase "$(INTDIR)\hexcode.obj"
	-@erase "$(INTDIR)\hexcode.sbr"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keyset.sbr"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\keysetdump.sbr"
	-@erase "$(INTDIR)\ksetdisp.obj"
	-@erase "$(INTDIR)\ksetdisp.sbr"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexalpha.sbr"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs.sbr"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\lexsubs4.sbr"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\linkdate.sbr"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\longtext.sbr"
	-@erase "$(INTDIR)\mainfrm.obj"
	-@erase "$(INTDIR)\mainfrm.sbr"
	-@erase "$(INTDIR)\mainpick.obj"
	-@erase "$(INTDIR)\mainpick.sbr"
	-@erase "$(INTDIR)\mapcontrol.obj"
	-@erase "$(INTDIR)\mapcontrol.sbr"
	-@erase "$(INTDIR)\mapscrollview.obj"
	-@erase "$(INTDIR)\mapscrollview.sbr"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\matchtable.sbr"
	-@erase "$(INTDIR)\menu.obj"
	-@erase "$(INTDIR)\menu.sbr"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messclean.sbr"
	-@erase "$(INTDIR)\messubs.obj"
	-@erase "$(INTDIR)\messubs.sbr"
	-@erase "$(INTDIR)\metab.obj"
	-@erase "$(INTDIR)\metab.sbr"
	-@erase "$(INTDIR)\method.obj"
	-@erase "$(INTDIR)\method.sbr"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\model.sbr"
	-@erase "$(INTDIR)\msgwindow.obj"
	-@erase "$(INTDIR)\msgwindow.sbr"
	-@erase "$(INTDIR)\newkey.obj"
	-@erase "$(INTDIR)\newkey.sbr"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\nicedump.sbr"
	-@erase "$(INTDIR)\nqcpass1.obj"
	-@erase "$(INTDIR)\nqcpass1.sbr"
	-@erase "$(INTDIR)\nqcpass2.obj"
	-@erase "$(INTDIR)\nqcpass2.sbr"
	-@erase "$(INTDIR)\nqctools.obj"
	-@erase "$(INTDIR)\nqctools.sbr"
	-@erase "$(INTDIR)\o2m.obj"
	-@erase "$(INTDIR)\o2m.sbr"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\objcache.sbr"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\oldhelp.sbr"
	-@erase "$(INTDIR)\oxgriddisp.obj"
	-@erase "$(INTDIR)\oxgriddisp.sbr"
	-@erase "$(INTDIR)\oxhomlist.obj"
	-@erase "$(INTDIR)\oxhomlist.sbr"
	-@erase "$(INTDIR)\pairmapdisp.obj"
	-@erase "$(INTDIR)\pairmapdisp.sbr"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\parse.sbr"
	-@erase "$(INTDIR)\pepactivezonecol.obj"
	-@erase "$(INTDIR)\pepactivezonecol.sbr"
	-@erase "$(INTDIR)\pepdisp.obj"
	-@erase "$(INTDIR)\pepdisp.sbr"
	-@erase "$(INTDIR)\pepfeaturecol.obj"
	-@erase "$(INTDIR)\pepfeaturecol.sbr"
	-@erase "$(INTDIR)\pepgraphcol.obj"
	-@erase "$(INTDIR)\pepgraphcol.sbr"
	-@erase "$(INTDIR)\pephomolcol.obj"
	-@erase "$(INTDIR)\pephomolcol.sbr"
	-@erase "$(INTDIR)\pepseqcol.obj"
	-@erase "$(INTDIR)\pepseqcol.sbr"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\peptide.sbr"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\picksubs.sbr"
	-@erase "$(INTDIR)\pixelfitview.obj"
	-@erase "$(INTDIR)\pixelfitview.sbr"
	-@erase "$(INTDIR)\pixelscrollview.obj"
	-@erase "$(INTDIR)\pixelscrollview.sbr"
	-@erase "$(INTDIR)\plainview.obj"
	-@erase "$(INTDIR)\plainview.sbr"
	-@erase "$(INTDIR)\plot.obj"
	-@erase "$(INTDIR)\plot.sbr"
	-@erase "$(INTDIR)\pmapconvert.obj"
	-@erase "$(INTDIR)\pmapconvert.sbr"
	-@erase "$(INTDIR)\pmapdisp.obj"
	-@erase "$(INTDIR)\pmapdisp.sbr"
	-@erase "$(INTDIR)\preferences.obj"
	-@erase "$(INTDIR)\preferences.sbr"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\prefsubs.sbr"
	-@erase "$(INTDIR)\qbedisp.obj"
	-@erase "$(INTDIR)\qbedisp.sbr"
	-@erase "$(INTDIR)\querybuild.obj"
	-@erase "$(INTDIR)\querybuild.sbr"
	-@erase "$(INTDIR)\querydisp.obj"
	-@erase "$(INTDIR)\querydisp.sbr"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\queryexe.sbr"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\quovadis.sbr"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\randsubs.sbr"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\session.sbr"
	-@erase "$(INTDIR)\specg.obj"
	-@erase "$(INTDIR)\specg.sbr"
	-@erase "$(INTDIR)\splashbox.obj"
	-@erase "$(INTDIR)\splashbox.sbr"
	-@erase "$(INTDIR)\sprdctrl.obj"
	-@erase "$(INTDIR)\sprdctrl.sbr"
	-@erase "$(INTDIR)\sprddata.obj"
	-@erase "$(INTDIR)\sprddata.sbr"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprddef.sbr"
	-@erase "$(INTDIR)\sprddisplay.obj"
	-@erase "$(INTDIR)\sprddisplay.sbr"
	-@erase "$(INTDIR)\sprdmap.obj"
	-@erase "$(INTDIR)\sprdmap.sbr"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\sprdop.sbr"
	-@erase "$(INTDIR)\status.obj"
	-@erase "$(INTDIR)\status.sbr"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\stdafx.sbr"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\sysclass.sbr"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\table.sbr"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\tags.sbr"
	-@erase "$(INTDIR)\textfitview.obj"
	-@erase "$(INTDIR)\textfitview.sbr"
	-@erase "$(INTDIR)\textscrollview.obj"
	-@erase "$(INTDIR)\textscrollview.sbr"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\timesubs.sbr"
	-@erase "$(INTDIR)\tqdebug.obj"
	-@erase "$(INTDIR)\tqdebug.sbr"
	-@erase "$(INTDIR)\translate.obj"
	-@erase "$(INTDIR)\translate.sbr"
	-@erase "$(INTDIR)\treedisp.obj"
	-@erase "$(INTDIR)\treedisp.sbr"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\update.sbr"
	-@erase "$(INTDIR)\userinterface.obj"
	-@erase "$(INTDIR)\userinterface.sbr"
	-@erase "$(INTDIR)\vc40.idb"
	-@erase "$(INTDIR)\vc40.pdb"
	-@erase "$(INTDIR)\viewedit.obj"
	-@erase "$(INTDIR)\viewedit.sbr"
	-@erase "$(INTDIR)\vmapdata2.obj"
	-@erase "$(INTDIR)\vmapdata2.sbr"
	-@erase "$(INTDIR)\vmapdisp.obj"
	-@erase "$(INTDIR)\vmapdisp.sbr"
	-@erase "$(INTDIR)\vmapdrag.obj"
	-@erase "$(INTDIR)\vmapdrag.sbr"
	-@erase "$(INTDIR)\vmapphys.obj"
	-@erase "$(INTDIR)\vmapphys.sbr"
	-@erase "$(INTDIR)\win32.obj"
	-@erase "$(INTDIR)\win32.sbr"
	-@erase "$(INTDIR)\win32lib.obj"
	-@erase "$(INTDIR)\win32lib.sbr"
	-@erase "$(INTDIR)\win32menusubs.obj"
	-@erase "$(INTDIR)\win32menusubs.sbr"
	-@erase "$(INTDIR)\win32print.obj"
	-@erase "$(INTDIR)\win32print.sbr"
	-@erase "$(INTDIR)\win32process.obj"
	-@erase "$(INTDIR)\win32process.sbr"
	-@erase "$(INTDIR)\win32thread.obj"
	-@erase "$(INTDIR)\win32thread.sbr"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\win32util.sbr"
	-@erase "$(INTDIR)\winace.obj"
	-@erase "$(INTDIR)\winace.res"
	-@erase "$(INTDIR)\winace.sbr"
	-@erase "$(INTDIR)\winacemail.obj"
	-@erase "$(INTDIR)\winacemail.sbr"
	-@erase "$(INTDIR)\windialogs.obj"
	-@erase "$(INTDIR)\windialogs.sbr"
	-@erase "$(INTDIR)\winfilquery.obj"
	-@erase "$(INTDIR)\winfilquery.sbr"
	-@erase "$(INTDIR)\winmain.obj"
	-@erase "$(INTDIR)\winmain.sbr"
	-@erase "$(OUTDIR)\winaced.bsc"
	-@erase "$(OUTDIR)\winaced.exe"
	-@erase "$(OUTDIR)\winaced.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR /YX /c
CPP_PROJ=/nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\windebug/
CPP_SBRS=\acedb\bin\windebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"\acedb\bin\windebug/winaced.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/winaced.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\acedbprofile.sbr" \
	"$(INTDIR)\acedbprofileview.sbr" \
	"$(INTDIR)\acedialogs.sbr" \
	"$(INTDIR)\acefiledoc.sbr" \
	"$(INTDIR)\acefileview.sbr" \
	"$(INTDIR)\action.sbr" \
	"$(INTDIR)\afxtrace.sbr" \
	"$(INTDIR)\align.sbr" \
	"$(INTDIR)\alignment.sbr" \
	"$(INTDIR)\arraysub.sbr" \
	"$(INTDIR)\asubs.sbr" \
	"$(INTDIR)\banner.sbr" \
	"$(INTDIR)\biblio.sbr" \
	"$(INTDIR)\blocksub.sbr" \
	"$(INTDIR)\blxview.sbr" \
	"$(INTDIR)\bsdumps.sbr" \
	"$(INTDIR)\bssubs.sbr" \
	"$(INTDIR)\bstools.sbr" \
	"$(INTDIR)\bstree.sbr" \
	"$(INTDIR)\bsubs.sbr" \
	"$(INTDIR)\bump.sbr" \
	"$(INTDIR)\caceprintpage.sbr" \
	"$(INTDIR)\call.sbr" \
	"$(INTDIR)\cgraph.sbr" \
	"$(INTDIR)\check.sbr" \
	"$(INTDIR)\chrono.sbr" \
	"$(INTDIR)\class.sbr" \
	"$(INTDIR)\cmapdisp.sbr" \
	"$(INTDIR)\cmultiass.sbr" \
	"$(INTDIR)\colcontrol.sbr" \
	"$(INTDIR)\colourbox.sbr" \
	"$(INTDIR)\dbprofile.sbr" \
	"$(INTDIR)\dict.sbr" \
	"$(INTDIR)\disknew.sbr" \
	"$(INTDIR)\display.sbr" \
	"$(INTDIR)\dnacpt.sbr" \
	"$(INTDIR)\dnasubs.sbr" \
	"$(INTDIR)\dotter.sbr" \
	"$(INTDIR)\dotterKarlin.sbr" \
	"$(INTDIR)\drawdisp.sbr" \
	"$(INTDIR)\dump.sbr" \
	"$(INTDIR)\embl.sbr" \
	"$(INTDIR)\filsubs.sbr" \
	"$(INTDIR)\fitview.sbr" \
	"$(INTDIR)\flag.sbr" \
	"$(INTDIR)\fmapblast.sbr" \
	"$(INTDIR)\fmapcontrol.sbr" \
	"$(INTDIR)\fmapfeatures.sbr" \
	"$(INTDIR)\fmapgene.sbr" \
	"$(INTDIR)\fmapmenes.sbr" \
	"$(INTDIR)\fmaposp.sbr" \
	"$(INTDIR)\fmapsequence.sbr" \
	"$(INTDIR)\fontpreference.sbr" \
	"$(INTDIR)\forest.sbr" \
	"$(INTDIR)\forestdisplay.sbr" \
	"$(INTDIR)\fpdisp.sbr" \
	"$(INTDIR)\freeout.sbr" \
	"$(INTDIR)\freesubs.sbr" \
	"$(INTDIR)\fullscrollview.sbr" \
	"$(INTDIR)\gd.sbr" \
	"$(INTDIR)\geldisp.sbr" \
	"$(INTDIR)\gfcode.sbr" \
	"$(INTDIR)\gmapconvert.sbr" \
	"$(INTDIR)\gmapdata.sbr" \
	"$(INTDIR)\gmapdatacol.sbr" \
	"$(INTDIR)\gmapdisp.sbr" \
	"$(INTDIR)\gmapintervalcol.sbr" \
	"$(INTDIR)\gmaplocuscol.sbr" \
	"$(INTDIR)\gmapmarkercol.sbr" \
	"$(INTDIR)\gmapphys.sbr" \
	"$(INTDIR)\gmapposnegcol.sbr" \
	"$(INTDIR)\gmapremarkcol.sbr" \
	"$(INTDIR)\graphbox.sbr" \
	"$(INTDIR)\graphcon.sbr" \
	"$(INTDIR)\graphramp.sbr" \
	"$(INTDIR)\graphselect.sbr" \
	"$(INTDIR)\graphsub.sbr" \
	"$(INTDIR)\graphview.sbr" \
	"$(INTDIR)\graphwin.sbr" \
	"$(INTDIR)\graphwin32_.sbr" \
	"$(INTDIR)\griddisp.sbr" \
	"$(INTDIR)\heap.sbr" \
	"$(INTDIR)\help.sbr" \
	"$(INTDIR)\hexcode.sbr" \
	"$(INTDIR)\keyset.sbr" \
	"$(INTDIR)\keysetdump.sbr" \
	"$(INTDIR)\ksetdisp.sbr" \
	"$(INTDIR)\lexalpha.sbr" \
	"$(INTDIR)\lexsubs.sbr" \
	"$(INTDIR)\lexsubs4.sbr" \
	"$(INTDIR)\linkdate.sbr" \
	"$(INTDIR)\longtext.sbr" \
	"$(INTDIR)\mainfrm.sbr" \
	"$(INTDIR)\mainpick.sbr" \
	"$(INTDIR)\mapcontrol.sbr" \
	"$(INTDIR)\mapscrollview.sbr" \
	"$(INTDIR)\matchtable.sbr" \
	"$(INTDIR)\menu.sbr" \
	"$(INTDIR)\messclean.sbr" \
	"$(INTDIR)\messubs.sbr" \
	"$(INTDIR)\metab.sbr" \
	"$(INTDIR)\method.sbr" \
	"$(INTDIR)\model.sbr" \
	"$(INTDIR)\msgwindow.sbr" \
	"$(INTDIR)\newkey.sbr" \
	"$(INTDIR)\nicedump.sbr" \
	"$(INTDIR)\nqcpass1.sbr" \
	"$(INTDIR)\nqcpass2.sbr" \
	"$(INTDIR)\nqctools.sbr" \
	"$(INTDIR)\o2m.sbr" \
	"$(INTDIR)\objcache.sbr" \
	"$(INTDIR)\oldhelp.sbr" \
	"$(INTDIR)\oxgriddisp.sbr" \
	"$(INTDIR)\oxhomlist.sbr" \
	"$(INTDIR)\pairmapdisp.sbr" \
	"$(INTDIR)\parse.sbr" \
	"$(INTDIR)\pepactivezonecol.sbr" \
	"$(INTDIR)\pepdisp.sbr" \
	"$(INTDIR)\pepfeaturecol.sbr" \
	"$(INTDIR)\pepgraphcol.sbr" \
	"$(INTDIR)\pephomolcol.sbr" \
	"$(INTDIR)\pepseqcol.sbr" \
	"$(INTDIR)\peptide.sbr" \
	"$(INTDIR)\picksubs.sbr" \
	"$(INTDIR)\pixelfitview.sbr" \
	"$(INTDIR)\pixelscrollview.sbr" \
	"$(INTDIR)\plainview.sbr" \
	"$(INTDIR)\plot.sbr" \
	"$(INTDIR)\pmapconvert.sbr" \
	"$(INTDIR)\pmapdisp.sbr" \
	"$(INTDIR)\preferences.sbr" \
	"$(INTDIR)\prefsubs.sbr" \
	"$(INTDIR)\qbedisp.sbr" \
	"$(INTDIR)\querybuild.sbr" \
	"$(INTDIR)\querydisp.sbr" \
	"$(INTDIR)\queryexe.sbr" \
	"$(INTDIR)\quovadis.sbr" \
	"$(INTDIR)\randsubs.sbr" \
	"$(INTDIR)\session.sbr" \
	"$(INTDIR)\specg.sbr" \
	"$(INTDIR)\splashbox.sbr" \
	"$(INTDIR)\sprdctrl.sbr" \
	"$(INTDIR)\sprddata.sbr" \
	"$(INTDIR)\sprddef.sbr" \
	"$(INTDIR)\sprddisplay.sbr" \
	"$(INTDIR)\sprdmap.sbr" \
	"$(INTDIR)\sprdop.sbr" \
	"$(INTDIR)\status.sbr" \
	"$(INTDIR)\stdafx.sbr" \
	"$(INTDIR)\sysclass.sbr" \
	"$(INTDIR)\table.sbr" \
	"$(INTDIR)\tags.sbr" \
	"$(INTDIR)\textfitview.sbr" \
	"$(INTDIR)\textscrollview.sbr" \
	"$(INTDIR)\timesubs.sbr" \
	"$(INTDIR)\tqdebug.sbr" \
	"$(INTDIR)\translate.sbr" \
	"$(INTDIR)\treedisp.sbr" \
	"$(INTDIR)\update.sbr" \
	"$(INTDIR)\userinterface.sbr" \
	"$(INTDIR)\viewedit.sbr" \
	"$(INTDIR)\vmapdata2.sbr" \
	"$(INTDIR)\vmapdisp.sbr" \
	"$(INTDIR)\vmapdrag.sbr" \
	"$(INTDIR)\vmapphys.sbr" \
	"$(INTDIR)\win32.sbr" \
	"$(INTDIR)\win32lib.sbr" \
	"$(INTDIR)\win32menusubs.sbr" \
	"$(INTDIR)\win32print.sbr" \
	"$(INTDIR)\win32process.sbr" \
	"$(INTDIR)\win32thread.sbr" \
	"$(INTDIR)\win32util.sbr" \
	"$(INTDIR)\winace.sbr" \
	"$(INTDIR)\winacemail.sbr" \
	"$(INTDIR)\windialogs.sbr" \
	"$(INTDIR)\winfilquery.sbr" \
	"$(INTDIR)\winmain.sbr"

"$(OUTDIR)\winaced.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386
# ADD LINK32 /nologo /version:4.31 /subsystem:windows /profile /map /debug /machine:I386 /out:"\acedb\bin\windebug/winaced.exe"
# SUBTRACT LINK32 /verbose /nodefaultlib
LINK32_FLAGS=/nologo /version:4.31 /subsystem:windows /profile\
 /map:"$(INTDIR)/winaced.map" /debug /machine:I386 /out:"$(OUTDIR)/winaced.exe" 
LINK32_OBJS= \
	"$(INTDIR)\acedbprofile.obj" \
	"$(INTDIR)\acedbprofileview.obj" \
	"$(INTDIR)\acedialogs.obj" \
	"$(INTDIR)\acefiledoc.obj" \
	"$(INTDIR)\acefileview.obj" \
	"$(INTDIR)\action.obj" \
	"$(INTDIR)\afxtrace.obj" \
	"$(INTDIR)\align.obj" \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\blxview.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\caceprintpage.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\cgraph.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\cmapdisp.obj" \
	"$(INTDIR)\cmultiass.obj" \
	"$(INTDIR)\colcontrol.obj" \
	"$(INTDIR)\colourbox.obj" \
	"$(INTDIR)\dbprofile.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\display.obj" \
	"$(INTDIR)\dnacpt.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dotter.obj" \
	"$(INTDIR)\dotterKarlin.obj" \
	"$(INTDIR)\drawdisp.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\embl.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\fitview.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\fmapblast.obj" \
	"$(INTDIR)\fmapcontrol.obj" \
	"$(INTDIR)\fmapfeatures.obj" \
	"$(INTDIR)\fmapgene.obj" \
	"$(INTDIR)\fmapmenes.obj" \
	"$(INTDIR)\fmaposp.obj" \
	"$(INTDIR)\fmapsequence.obj" \
	"$(INTDIR)\fontpreference.obj" \
	"$(INTDIR)\forest.obj" \
	"$(INTDIR)\forestdisplay.obj" \
	"$(INTDIR)\fpdisp.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\fullscrollview.obj" \
	"$(INTDIR)\gd.obj" \
	"$(INTDIR)\geldisp.obj" \
	"$(INTDIR)\gfcode.obj" \
	"$(INTDIR)\gmapconvert.obj" \
	"$(INTDIR)\gmapdata.obj" \
	"$(INTDIR)\gmapdatacol.obj" \
	"$(INTDIR)\gmapdisp.obj" \
	"$(INTDIR)\gmapintervalcol.obj" \
	"$(INTDIR)\gmaplocuscol.obj" \
	"$(INTDIR)\gmapmarkercol.obj" \
	"$(INTDIR)\gmapphys.obj" \
	"$(INTDIR)\gmapposnegcol.obj" \
	"$(INTDIR)\gmapremarkcol.obj" \
	"$(INTDIR)\graphbox.obj" \
	"$(INTDIR)\graphcon.obj" \
	"$(INTDIR)\graphramp.obj" \
	"$(INTDIR)\graphselect.obj" \
	"$(INTDIR)\graphsub.obj" \
	"$(INTDIR)\graphview.obj" \
	"$(INTDIR)\graphwin.obj" \
	"$(INTDIR)\graphwin32_.obj" \
	"$(INTDIR)\griddisp.obj" \
	"$(INTDIR)\heap.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\hexcode.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\ksetdisp.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\mainfrm.obj" \
	"$(INTDIR)\mainpick.obj" \
	"$(INTDIR)\mapcontrol.obj" \
	"$(INTDIR)\mapscrollview.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\menu.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\messubs.obj" \
	"$(INTDIR)\metab.obj" \
	"$(INTDIR)\method.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\msgwindow.obj" \
	"$(INTDIR)\newkey.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\nqcpass1.obj" \
	"$(INTDIR)\nqcpass2.obj" \
	"$(INTDIR)\nqctools.obj" \
	"$(INTDIR)\o2m.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\oxgriddisp.obj" \
	"$(INTDIR)\oxhomlist.obj" \
	"$(INTDIR)\pairmapdisp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\pepactivezonecol.obj" \
	"$(INTDIR)\pepdisp.obj" \
	"$(INTDIR)\pepfeaturecol.obj" \
	"$(INTDIR)\pepgraphcol.obj" \
	"$(INTDIR)\pephomolcol.obj" \
	"$(INTDIR)\pepseqcol.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\pixelfitview.obj" \
	"$(INTDIR)\pixelscrollview.obj" \
	"$(INTDIR)\plainview.obj" \
	"$(INTDIR)\plot.obj" \
	"$(INTDIR)\pmapconvert.obj" \
	"$(INTDIR)\pmapdisp.obj" \
	"$(INTDIR)\preferences.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\qbedisp.obj" \
	"$(INTDIR)\querybuild.obj" \
	"$(INTDIR)\querydisp.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\specg.obj" \
	"$(INTDIR)\splashbox.obj" \
	"$(INTDIR)\sprdctrl.obj" \
	"$(INTDIR)\sprddata.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprddisplay.obj" \
	"$(INTDIR)\sprdmap.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\status.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\textfitview.obj" \
	"$(INTDIR)\textscrollview.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\tqdebug.obj" \
	"$(INTDIR)\translate.obj" \
	"$(INTDIR)\treedisp.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\userinterface.obj" \
	"$(INTDIR)\viewedit.obj" \
	"$(INTDIR)\vmapdata2.obj" \
	"$(INTDIR)\vmapdisp.obj" \
	"$(INTDIR)\vmapdrag.obj" \
	"$(INTDIR)\vmapphys.obj" \
	"$(INTDIR)\win32.obj" \
	"$(INTDIR)\win32lib.obj" \
	"$(INTDIR)\win32menusubs.obj" \
	"$(INTDIR)\win32print.obj" \
	"$(INTDIR)\win32process.obj" \
	"$(INTDIR)\win32thread.obj" \
	"$(INTDIR)\win32util.obj" \
	"$(INTDIR)\winace.obj" \
	"$(INTDIR)\winace.res" \
	"$(INTDIR)\winacemail.obj" \
	"$(INTDIR)\windialogs.obj" \
	"$(INTDIR)\winfilquery.obj" \
	"$(INTDIR)\winmain.obj"

"$(OUTDIR)\winaced.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "acedb__"
# PROP BASE Intermediate_Dir "acedb__"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "\acedb\bin\acemblydebug"
# PROP Intermediate_Dir "\acedb\bin\acemblydebug"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\acemblydebug
INTDIR=\acedb\bin\acemblydebug

ALL : "$(OUTDIR)\winacembly.exe" "$(OUTDIR)\winacembly.bsc"

CLEAN : 
	-@erase "$(INTDIR)\abifix.obj"
	-@erase "$(INTDIR)\abifix.sbr"
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\acedbprofile.obj"
	-@erase "$(INTDIR)\acedbprofile.sbr"
	-@erase "$(INTDIR)\acedbprofileview.obj"
	-@erase "$(INTDIR)\acedbprofileview.sbr"
	-@erase "$(INTDIR)\acedialogs.obj"
	-@erase "$(INTDIR)\acedialogs.sbr"
	-@erase "$(INTDIR)\acefiledoc.obj"
	-@erase "$(INTDIR)\acefiledoc.sbr"
	-@erase "$(INTDIR)\acefileview.obj"
	-@erase "$(INTDIR)\acefileview.sbr"
	-@erase "$(INTDIR)\acemblyhook.obj"
	-@erase "$(INTDIR)\acemblyhook.sbr"
	-@erase "$(INTDIR)\action.obj"
	-@erase "$(INTDIR)\action.sbr"
	-@erase "$(INTDIR)\afxtrace.obj"
	-@erase "$(INTDIR)\afxtrace.sbr"
	-@erase "$(INTDIR)\align.obj"
	-@erase "$(INTDIR)\align.sbr"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\alignment.sbr"
	-@erase "$(INTDIR)\aligntools.obj"
	-@erase "$(INTDIR)\aligntools.sbr"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\arraysub.sbr"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\asubs.sbr"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\banner.sbr"
	-@erase "$(INTDIR)\basecall.obj"
	-@erase "$(INTDIR)\basecall.sbr"
	-@erase "$(INTDIR)\basecallstat.obj"
	-@erase "$(INTDIR)\basecallstat.sbr"
	-@erase "$(INTDIR)\basepad.obj"
	-@erase "$(INTDIR)\basepad.sbr"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\biblio.sbr"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\blocksub.sbr"
	-@erase "$(INTDIR)\blxview.obj"
	-@erase "$(INTDIR)\blxview.sbr"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bsdumps.sbr"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bssubs.sbr"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstools.sbr"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bstree.sbr"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bsubs.sbr"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\bump.sbr"
	-@erase "$(INTDIR)\caceprintpage.obj"
	-@erase "$(INTDIR)\caceprintpage.sbr"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\call.sbr"
	-@erase "$(INTDIR)\cgraph.obj"
	-@erase "$(INTDIR)\cgraph.sbr"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\check.sbr"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\chrono.sbr"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\class.sbr"
	-@erase "$(INTDIR)\cmapdisp.obj"
	-@erase "$(INTDIR)\cmapdisp.sbr"
	-@erase "$(INTDIR)\cmultiass.obj"
	-@erase "$(INTDIR)\cmultiass.sbr"
	-@erase "$(INTDIR)\colcontrol.obj"
	-@erase "$(INTDIR)\colcontrol.sbr"
	-@erase "$(INTDIR)\colourbox.obj"
	-@erase "$(INTDIR)\colourbox.sbr"
	-@erase "$(INTDIR)\command.obj"
	-@erase "$(INTDIR)\command.sbr"
	-@erase "$(INTDIR)\dbprofile.obj"
	-@erase "$(INTDIR)\dbprofile.sbr"
	-@erase "$(INTDIR)\defcpt.obj"
	-@erase "$(INTDIR)\defcpt.sbr"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\dict.sbr"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\disknew.sbr"
	-@erase "$(INTDIR)\display.obj"
	-@erase "$(INTDIR)\display.sbr"
	-@erase "$(INTDIR)\dnaalign.obj"
	-@erase "$(INTDIR)\dnaalign.sbr"
	-@erase "$(INTDIR)\dnacpt.obj"
	-@erase "$(INTDIR)\dnacpt.sbr"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dnasubs.sbr"
	-@erase "$(INTDIR)\dotter.obj"
	-@erase "$(INTDIR)\dotter.sbr"
	-@erase "$(INTDIR)\dotterKarlin.obj"
	-@erase "$(INTDIR)\dotterKarlin.sbr"
	-@erase "$(INTDIR)\drawdisp.obj"
	-@erase "$(INTDIR)\drawdisp.sbr"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\dump.sbr"
	-@erase "$(INTDIR)\embl.obj"
	-@erase "$(INTDIR)\embl.sbr"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filqng.sbr"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\filsubs.sbr"
	-@erase "$(INTDIR)\fitview.obj"
	-@erase "$(INTDIR)\fitview.sbr"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\flag.sbr"
	-@erase "$(INTDIR)\fmapblast.obj"
	-@erase "$(INTDIR)\fmapblast.sbr"
	-@erase "$(INTDIR)\fmapcontrol.obj"
	-@erase "$(INTDIR)\fmapcontrol.sbr"
	-@erase "$(INTDIR)\fmapfeatures.obj"
	-@erase "$(INTDIR)\fmapfeatures.sbr"
	-@erase "$(INTDIR)\fmapgene.obj"
	-@erase "$(INTDIR)\fmapgene.sbr"
	-@erase "$(INTDIR)\fmaposp.obj"
	-@erase "$(INTDIR)\fmaposp.sbr"
	-@erase "$(INTDIR)\fmapsequence.obj"
	-@erase "$(INTDIR)\fmapsequence.sbr"
	-@erase "$(INTDIR)\fmaptrace.obj"
	-@erase "$(INTDIR)\fmaptrace.sbr"
	-@erase "$(INTDIR)\fontpreference.obj"
	-@erase "$(INTDIR)\fontpreference.sbr"
	-@erase "$(INTDIR)\forest.obj"
	-@erase "$(INTDIR)\forest.sbr"
	-@erase "$(INTDIR)\forestdisplay.obj"
	-@erase "$(INTDIR)\forestdisplay.sbr"
	-@erase "$(INTDIR)\fpdisp.obj"
	-@erase "$(INTDIR)\fpdisp.sbr"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freeout.sbr"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\freesubs.sbr"
	-@erase "$(INTDIR)\fullscrollview.obj"
	-@erase "$(INTDIR)\fullscrollview.sbr"
	-@erase "$(INTDIR)\gd.obj"
	-@erase "$(INTDIR)\gd.sbr"
	-@erase "$(INTDIR)\geldisp.obj"
	-@erase "$(INTDIR)\geldisp.sbr"
	-@erase "$(INTDIR)\gfcode.obj"
	-@erase "$(INTDIR)\gfcode.sbr"
	-@erase "$(INTDIR)\gmapconvert.obj"
	-@erase "$(INTDIR)\gmapconvert.sbr"
	-@erase "$(INTDIR)\gmapdata.obj"
	-@erase "$(INTDIR)\gmapdata.sbr"
	-@erase "$(INTDIR)\gmapdatacol.obj"
	-@erase "$(INTDIR)\gmapdatacol.sbr"
	-@erase "$(INTDIR)\gmapdisp.obj"
	-@erase "$(INTDIR)\gmapdisp.sbr"
	-@erase "$(INTDIR)\gmapintervalcol.obj"
	-@erase "$(INTDIR)\gmapintervalcol.sbr"
	-@erase "$(INTDIR)\gmaplocuscol.obj"
	-@erase "$(INTDIR)\gmaplocuscol.sbr"
	-@erase "$(INTDIR)\gmapmarkercol.obj"
	-@erase "$(INTDIR)\gmapmarkercol.sbr"
	-@erase "$(INTDIR)\gmapphys.obj"
	-@erase "$(INTDIR)\gmapphys.sbr"
	-@erase "$(INTDIR)\gmapposnegcol.obj"
	-@erase "$(INTDIR)\gmapposnegcol.sbr"
	-@erase "$(INTDIR)\gmapremarkcol.obj"
	-@erase "$(INTDIR)\gmapremarkcol.sbr"
	-@erase "$(INTDIR)\graphbox.obj"
	-@erase "$(INTDIR)\graphbox.sbr"
	-@erase "$(INTDIR)\graphcon.obj"
	-@erase "$(INTDIR)\graphcon.sbr"
	-@erase "$(INTDIR)\graphramp.obj"
	-@erase "$(INTDIR)\graphramp.sbr"
	-@erase "$(INTDIR)\graphselect.obj"
	-@erase "$(INTDIR)\graphselect.sbr"
	-@erase "$(INTDIR)\graphsub.obj"
	-@erase "$(INTDIR)\graphsub.sbr"
	-@erase "$(INTDIR)\graphview.obj"
	-@erase "$(INTDIR)\graphview.sbr"
	-@erase "$(INTDIR)\graphwin.obj"
	-@erase "$(INTDIR)\graphwin.sbr"
	-@erase "$(INTDIR)\graphwin32_.obj"
	-@erase "$(INTDIR)\graphwin32_.sbr"
	-@erase "$(INTDIR)\griddisp.obj"
	-@erase "$(INTDIR)\griddisp.sbr"
	-@erase "$(INTDIR)\heap.obj"
	-@erase "$(INTDIR)\heap.sbr"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\help.sbr"
	-@erase "$(INTDIR)\hexcode.obj"
	-@erase "$(INTDIR)\hexcode.sbr"
	-@erase "$(INTDIR)\intrinsictree.obj"
	-@erase "$(INTDIR)\intrinsictree.sbr"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keyset.sbr"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\keysetdump.sbr"
	-@erase "$(INTDIR)\ksetdisp.obj"
	-@erase "$(INTDIR)\ksetdisp.sbr"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexalpha.sbr"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs.sbr"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\lexsubs4.sbr"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\linkdate.sbr"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\longtext.sbr"
	-@erase "$(INTDIR)\mach-io.obj"
	-@erase "$(INTDIR)\mach-io.sbr"
	-@erase "$(INTDIR)\mainfrm.obj"
	-@erase "$(INTDIR)\mainfrm.sbr"
	-@erase "$(INTDIR)\mainpick.obj"
	-@erase "$(INTDIR)\mainpick.sbr"
	-@erase "$(INTDIR)\mapcontrol.obj"
	-@erase "$(INTDIR)\mapcontrol.sbr"
	-@erase "$(INTDIR)\mapscrollview.obj"
	-@erase "$(INTDIR)\mapscrollview.sbr"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\matchtable.sbr"
	-@erase "$(INTDIR)\menu.obj"
	-@erase "$(INTDIR)\menu.sbr"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messclean.sbr"
	-@erase "$(INTDIR)\messubs.obj"
	-@erase "$(INTDIR)\messubs.sbr"
	-@erase "$(INTDIR)\metab.obj"
	-@erase "$(INTDIR)\metab.sbr"
	-@erase "$(INTDIR)\method.obj"
	-@erase "$(INTDIR)\method.sbr"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\model.sbr"
	-@erase "$(INTDIR)\msgwindow.obj"
	-@erase "$(INTDIR)\msgwindow.sbr"
	-@erase "$(INTDIR)\myNetwork.obj"
	-@erase "$(INTDIR)\myNetwork.sbr"
	-@erase "$(INTDIR)\newkey.obj"
	-@erase "$(INTDIR)\newkey.sbr"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\nicedump.sbr"
	-@erase "$(INTDIR)\nqcpass1.obj"
	-@erase "$(INTDIR)\nqcpass1.sbr"
	-@erase "$(INTDIR)\nqcpass2.obj"
	-@erase "$(INTDIR)\nqcpass2.sbr"
	-@erase "$(INTDIR)\nqctools.obj"
	-@erase "$(INTDIR)\nqctools.sbr"
	-@erase "$(INTDIR)\o2m.obj"
	-@erase "$(INTDIR)\o2m.sbr"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\objcache.sbr"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\oldhelp.sbr"
	-@erase "$(INTDIR)\opp.obj"
	-@erase "$(INTDIR)\opp.sbr"
	-@erase "$(INTDIR)\oxgriddisp.obj"
	-@erase "$(INTDIR)\oxgriddisp.sbr"
	-@erase "$(INTDIR)\oxhomlist.obj"
	-@erase "$(INTDIR)\oxhomlist.sbr"
	-@erase "$(INTDIR)\pairmapdisp.obj"
	-@erase "$(INTDIR)\pairmapdisp.sbr"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\parse.sbr"
	-@erase "$(INTDIR)\pepactivezonecol.obj"
	-@erase "$(INTDIR)\pepactivezonecol.sbr"
	-@erase "$(INTDIR)\pepdisp.obj"
	-@erase "$(INTDIR)\pepdisp.sbr"
	-@erase "$(INTDIR)\pepfeaturecol.obj"
	-@erase "$(INTDIR)\pepfeaturecol.sbr"
	-@erase "$(INTDIR)\pepgraphcol.obj"
	-@erase "$(INTDIR)\pepgraphcol.sbr"
	-@erase "$(INTDIR)\pephomolcol.obj"
	-@erase "$(INTDIR)\pephomolcol.sbr"
	-@erase "$(INTDIR)\pepseqcol.obj"
	-@erase "$(INTDIR)\pepseqcol.sbr"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\peptide.sbr"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\picksubs.sbr"
	-@erase "$(INTDIR)\pixelfitview.obj"
	-@erase "$(INTDIR)\pixelfitview.sbr"
	-@erase "$(INTDIR)\pixelscrollview.obj"
	-@erase "$(INTDIR)\pixelscrollview.sbr"
	-@erase "$(INTDIR)\plainview.obj"
	-@erase "$(INTDIR)\plainview.sbr"
	-@erase "$(INTDIR)\plot.obj"
	-@erase "$(INTDIR)\plot.sbr"
	-@erase "$(INTDIR)\pmapconvert.obj"
	-@erase "$(INTDIR)\pmapconvert.sbr"
	-@erase "$(INTDIR)\pmapdisp.obj"
	-@erase "$(INTDIR)\pmapdisp.sbr"
	-@erase "$(INTDIR)\preferences.obj"
	-@erase "$(INTDIR)\preferences.sbr"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\prefsubs.sbr"
	-@erase "$(INTDIR)\qbedisp.obj"
	-@erase "$(INTDIR)\qbedisp.sbr"
	-@erase "$(INTDIR)\querybuild.obj"
	-@erase "$(INTDIR)\querybuild.sbr"
	-@erase "$(INTDIR)\querydisp.obj"
	-@erase "$(INTDIR)\querydisp.sbr"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\queryexe.sbr"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\quovadis.sbr"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\randsubs.sbr"
	-@erase "$(INTDIR)\regression.obj"
	-@erase "$(INTDIR)\regression.sbr"
	-@erase "$(INTDIR)\seq.obj"
	-@erase "$(INTDIR)\seq.sbr"
	-@erase "$(INTDIR)\seqIOSCF.obj"
	-@erase "$(INTDIR)\seqIOSCF.sbr"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\session.sbr"
	-@erase "$(INTDIR)\specg.obj"
	-@erase "$(INTDIR)\specg.sbr"
	-@erase "$(INTDIR)\splashbox.obj"
	-@erase "$(INTDIR)\splashbox.sbr"
	-@erase "$(INTDIR)\sprdctrl.obj"
	-@erase "$(INTDIR)\sprdctrl.sbr"
	-@erase "$(INTDIR)\sprddata.obj"
	-@erase "$(INTDIR)\sprddata.sbr"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprddef.sbr"
	-@erase "$(INTDIR)\sprddisplay.obj"
	-@erase "$(INTDIR)\sprddisplay.sbr"
	-@erase "$(INTDIR)\sprdmap.obj"
	-@erase "$(INTDIR)\sprdmap.sbr"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\sprdop.sbr"
	-@erase "$(INTDIR)\status.obj"
	-@erase "$(INTDIR)\status.sbr"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\stdafx.sbr"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\sysclass.sbr"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\table.sbr"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\tags.sbr"
	-@erase "$(INTDIR)\textfitview.obj"
	-@erase "$(INTDIR)\textfitview.sbr"
	-@erase "$(INTDIR)\textscrollview.obj"
	-@erase "$(INTDIR)\textscrollview.sbr"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\timesubs.sbr"
	-@erase "$(INTDIR)\topology.obj"
	-@erase "$(INTDIR)\topology.sbr"
	-@erase "$(INTDIR)\tqdebug.obj"
	-@erase "$(INTDIR)\tqdebug.sbr"
	-@erase "$(INTDIR)\trace.obj"
	-@erase "$(INTDIR)\trace.sbr"
	-@erase "$(INTDIR)\translate.obj"
	-@erase "$(INTDIR)\translate.sbr"
	-@erase "$(INTDIR)\treedisp.obj"
	-@erase "$(INTDIR)\treedisp.sbr"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\update.sbr"
	-@erase "$(INTDIR)\userinterface.obj"
	-@erase "$(INTDIR)\userinterface.sbr"
	-@erase "$(INTDIR)\vc40.idb"
	-@erase "$(INTDIR)\vc40.pdb"
	-@erase "$(INTDIR)\viewedit.obj"
	-@erase "$(INTDIR)\viewedit.sbr"
	-@erase "$(INTDIR)\vmapdata2.obj"
	-@erase "$(INTDIR)\vmapdata2.sbr"
	-@erase "$(INTDIR)\vmapdisp.obj"
	-@erase "$(INTDIR)\vmapdisp.sbr"
	-@erase "$(INTDIR)\vmapdrag.obj"
	-@erase "$(INTDIR)\vmapdrag.sbr"
	-@erase "$(INTDIR)\vmapphys.obj"
	-@erase "$(INTDIR)\vmapphys.sbr"
	-@erase "$(INTDIR)\win32.obj"
	-@erase "$(INTDIR)\win32.sbr"
	-@erase "$(INTDIR)\win32lib.obj"
	-@erase "$(INTDIR)\win32lib.sbr"
	-@erase "$(INTDIR)\win32menusubs.obj"
	-@erase "$(INTDIR)\win32menusubs.sbr"
	-@erase "$(INTDIR)\win32print.obj"
	-@erase "$(INTDIR)\win32print.sbr"
	-@erase "$(INTDIR)\win32process.obj"
	-@erase "$(INTDIR)\win32process.sbr"
	-@erase "$(INTDIR)\win32thread.obj"
	-@erase "$(INTDIR)\win32thread.sbr"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\win32util.sbr"
	-@erase "$(INTDIR)\winace.obj"
	-@erase "$(INTDIR)\winace.res"
	-@erase "$(INTDIR)\winace.sbr"
	-@erase "$(INTDIR)\winacemail.obj"
	-@erase "$(INTDIR)\winacemail.sbr"
	-@erase "$(INTDIR)\windialogs.obj"
	-@erase "$(INTDIR)\windialogs.sbr"
	-@erase "$(INTDIR)\winfilquery.obj"
	-@erase "$(INTDIR)\winfilquery.sbr"
	-@erase "$(INTDIR)\winmain.obj"
	-@erase "$(INTDIR)\winmain.sbr"
	-@erase "$(OUTDIR)\winacembly.bsc"
	-@erase "$(OUTDIR)\winacembly.exe"
	-@erase "$(OUTDIR)\winacembly.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "acedb4" /D "_AFXDLL" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR /YX /c
CPP_PROJ=/nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\acemblydebug/
CPP_SBRS=\acedb\bin\acemblydebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"\acedb\bin\acemblydebug/winacembly.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/winacembly.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\abifix.sbr" \
	"$(INTDIR)\acedbprofile.sbr" \
	"$(INTDIR)\acedbprofileview.sbr" \
	"$(INTDIR)\acedialogs.sbr" \
	"$(INTDIR)\acefiledoc.sbr" \
	"$(INTDIR)\acefileview.sbr" \
	"$(INTDIR)\acemblyhook.sbr" \
	"$(INTDIR)\action.sbr" \
	"$(INTDIR)\afxtrace.sbr" \
	"$(INTDIR)\align.sbr" \
	"$(INTDIR)\alignment.sbr" \
	"$(INTDIR)\aligntools.sbr" \
	"$(INTDIR)\arraysub.sbr" \
	"$(INTDIR)\asubs.sbr" \
	"$(INTDIR)\banner.sbr" \
	"$(INTDIR)\basecall.sbr" \
	"$(INTDIR)\basecallstat.sbr" \
	"$(INTDIR)\basepad.sbr" \
	"$(INTDIR)\biblio.sbr" \
	"$(INTDIR)\blocksub.sbr" \
	"$(INTDIR)\blxview.sbr" \
	"$(INTDIR)\bsdumps.sbr" \
	"$(INTDIR)\bssubs.sbr" \
	"$(INTDIR)\bstools.sbr" \
	"$(INTDIR)\bstree.sbr" \
	"$(INTDIR)\bsubs.sbr" \
	"$(INTDIR)\bump.sbr" \
	"$(INTDIR)\caceprintpage.sbr" \
	"$(INTDIR)\call.sbr" \
	"$(INTDIR)\cgraph.sbr" \
	"$(INTDIR)\check.sbr" \
	"$(INTDIR)\chrono.sbr" \
	"$(INTDIR)\class.sbr" \
	"$(INTDIR)\cmapdisp.sbr" \
	"$(INTDIR)\cmultiass.sbr" \
	"$(INTDIR)\colcontrol.sbr" \
	"$(INTDIR)\colourbox.sbr" \
	"$(INTDIR)\command.sbr" \
	"$(INTDIR)\dbprofile.sbr" \
	"$(INTDIR)\defcpt.sbr" \
	"$(INTDIR)\dict.sbr" \
	"$(INTDIR)\disknew.sbr" \
	"$(INTDIR)\display.sbr" \
	"$(INTDIR)\dnaalign.sbr" \
	"$(INTDIR)\dnacpt.sbr" \
	"$(INTDIR)\dnasubs.sbr" \
	"$(INTDIR)\dotter.sbr" \
	"$(INTDIR)\dotterKarlin.sbr" \
	"$(INTDIR)\drawdisp.sbr" \
	"$(INTDIR)\dump.sbr" \
	"$(INTDIR)\embl.sbr" \
	"$(INTDIR)\filqng.sbr" \
	"$(INTDIR)\filsubs.sbr" \
	"$(INTDIR)\fitview.sbr" \
	"$(INTDIR)\flag.sbr" \
	"$(INTDIR)\fmapblast.sbr" \
	"$(INTDIR)\fmapcontrol.sbr" \
	"$(INTDIR)\fmapfeatures.sbr" \
	"$(INTDIR)\fmapgene.sbr" \
	"$(INTDIR)\fmaposp.sbr" \
	"$(INTDIR)\fmapsequence.sbr" \
	"$(INTDIR)\fmaptrace.sbr" \
	"$(INTDIR)\fontpreference.sbr" \
	"$(INTDIR)\forest.sbr" \
	"$(INTDIR)\forestdisplay.sbr" \
	"$(INTDIR)\fpdisp.sbr" \
	"$(INTDIR)\freeout.sbr" \
	"$(INTDIR)\freesubs.sbr" \
	"$(INTDIR)\fullscrollview.sbr" \
	"$(INTDIR)\gd.sbr" \
	"$(INTDIR)\geldisp.sbr" \
	"$(INTDIR)\gfcode.sbr" \
	"$(INTDIR)\gmapconvert.sbr" \
	"$(INTDIR)\gmapdata.sbr" \
	"$(INTDIR)\gmapdatacol.sbr" \
	"$(INTDIR)\gmapdisp.sbr" \
	"$(INTDIR)\gmapintervalcol.sbr" \
	"$(INTDIR)\gmaplocuscol.sbr" \
	"$(INTDIR)\gmapmarkercol.sbr" \
	"$(INTDIR)\gmapphys.sbr" \
	"$(INTDIR)\gmapposnegcol.sbr" \
	"$(INTDIR)\gmapremarkcol.sbr" \
	"$(INTDIR)\graphbox.sbr" \
	"$(INTDIR)\graphcon.sbr" \
	"$(INTDIR)\graphramp.sbr" \
	"$(INTDIR)\graphselect.sbr" \
	"$(INTDIR)\graphsub.sbr" \
	"$(INTDIR)\graphview.sbr" \
	"$(INTDIR)\graphwin.sbr" \
	"$(INTDIR)\graphwin32_.sbr" \
	"$(INTDIR)\griddisp.sbr" \
	"$(INTDIR)\heap.sbr" \
	"$(INTDIR)\help.sbr" \
	"$(INTDIR)\hexcode.sbr" \
	"$(INTDIR)\intrinsictree.sbr" \
	"$(INTDIR)\keyset.sbr" \
	"$(INTDIR)\keysetdump.sbr" \
	"$(INTDIR)\ksetdisp.sbr" \
	"$(INTDIR)\lexalpha.sbr" \
	"$(INTDIR)\lexsubs.sbr" \
	"$(INTDIR)\lexsubs4.sbr" \
	"$(INTDIR)\linkdate.sbr" \
	"$(INTDIR)\longtext.sbr" \
	"$(INTDIR)\mach-io.sbr" \
	"$(INTDIR)\mainfrm.sbr" \
	"$(INTDIR)\mainpick.sbr" \
	"$(INTDIR)\mapcontrol.sbr" \
	"$(INTDIR)\mapscrollview.sbr" \
	"$(INTDIR)\matchtable.sbr" \
	"$(INTDIR)\menu.sbr" \
	"$(INTDIR)\messclean.sbr" \
	"$(INTDIR)\messubs.sbr" \
	"$(INTDIR)\metab.sbr" \
	"$(INTDIR)\method.sbr" \
	"$(INTDIR)\model.sbr" \
	"$(INTDIR)\msgwindow.sbr" \
	"$(INTDIR)\myNetwork.sbr" \
	"$(INTDIR)\newkey.sbr" \
	"$(INTDIR)\nicedump.sbr" \
	"$(INTDIR)\nqcpass1.sbr" \
	"$(INTDIR)\nqcpass2.sbr" \
	"$(INTDIR)\nqctools.sbr" \
	"$(INTDIR)\o2m.sbr" \
	"$(INTDIR)\objcache.sbr" \
	"$(INTDIR)\oldhelp.sbr" \
	"$(INTDIR)\opp.sbr" \
	"$(INTDIR)\oxgriddisp.sbr" \
	"$(INTDIR)\oxhomlist.sbr" \
	"$(INTDIR)\pairmapdisp.sbr" \
	"$(INTDIR)\parse.sbr" \
	"$(INTDIR)\pepactivezonecol.sbr" \
	"$(INTDIR)\pepdisp.sbr" \
	"$(INTDIR)\pepfeaturecol.sbr" \
	"$(INTDIR)\pepgraphcol.sbr" \
	"$(INTDIR)\pephomolcol.sbr" \
	"$(INTDIR)\pepseqcol.sbr" \
	"$(INTDIR)\peptide.sbr" \
	"$(INTDIR)\picksubs.sbr" \
	"$(INTDIR)\pixelfitview.sbr" \
	"$(INTDIR)\pixelscrollview.sbr" \
	"$(INTDIR)\plainview.sbr" \
	"$(INTDIR)\plot.sbr" \
	"$(INTDIR)\pmapconvert.sbr" \
	"$(INTDIR)\pmapdisp.sbr" \
	"$(INTDIR)\preferences.sbr" \
	"$(INTDIR)\prefsubs.sbr" \
	"$(INTDIR)\qbedisp.sbr" \
	"$(INTDIR)\querybuild.sbr" \
	"$(INTDIR)\querydisp.sbr" \
	"$(INTDIR)\queryexe.sbr" \
	"$(INTDIR)\quovadis.sbr" \
	"$(INTDIR)\randsubs.sbr" \
	"$(INTDIR)\regression.sbr" \
	"$(INTDIR)\seq.sbr" \
	"$(INTDIR)\seqIOSCF.sbr" \
	"$(INTDIR)\session.sbr" \
	"$(INTDIR)\specg.sbr" \
	"$(INTDIR)\splashbox.sbr" \
	"$(INTDIR)\sprdctrl.sbr" \
	"$(INTDIR)\sprddata.sbr" \
	"$(INTDIR)\sprddef.sbr" \
	"$(INTDIR)\sprddisplay.sbr" \
	"$(INTDIR)\sprdmap.sbr" \
	"$(INTDIR)\sprdop.sbr" \
	"$(INTDIR)\status.sbr" \
	"$(INTDIR)\stdafx.sbr" \
	"$(INTDIR)\sysclass.sbr" \
	"$(INTDIR)\table.sbr" \
	"$(INTDIR)\tags.sbr" \
	"$(INTDIR)\textfitview.sbr" \
	"$(INTDIR)\textscrollview.sbr" \
	"$(INTDIR)\timesubs.sbr" \
	"$(INTDIR)\topology.sbr" \
	"$(INTDIR)\tqdebug.sbr" \
	"$(INTDIR)\trace.sbr" \
	"$(INTDIR)\translate.sbr" \
	"$(INTDIR)\treedisp.sbr" \
	"$(INTDIR)\update.sbr" \
	"$(INTDIR)\userinterface.sbr" \
	"$(INTDIR)\viewedit.sbr" \
	"$(INTDIR)\vmapdata2.sbr" \
	"$(INTDIR)\vmapdisp.sbr" \
	"$(INTDIR)\vmapdrag.sbr" \
	"$(INTDIR)\vmapphys.sbr" \
	"$(INTDIR)\win32.sbr" \
	"$(INTDIR)\win32lib.sbr" \
	"$(INTDIR)\win32menusubs.sbr" \
	"$(INTDIR)\win32print.sbr" \
	"$(INTDIR)\win32process.sbr" \
	"$(INTDIR)\win32thread.sbr" \
	"$(INTDIR)\win32util.sbr" \
	"$(INTDIR)\winace.sbr" \
	"$(INTDIR)\winacemail.sbr" \
	"$(INTDIR)\windialogs.sbr" \
	"$(INTDIR)\winfilquery.sbr" \
	"$(INTDIR)\winmain.sbr"

"$(OUTDIR)\winacembly.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /profile /map /debug /machine:I386 /nodefaultlib:"libc.lib" /nodefaultlib:"libcmt.lib" /nodefaultlib:"msvcrtd.lib" /nodefaultlib:"libcd.lib" /nodefaultlib:"libcmtd.lib"
# SUBTRACT BASE LINK32 /verbose /nodefaultlib
# ADD LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /profile /map /debug /machine:I386 /nodefaultlib:"libcmtd.lib" /out:"\acedb\bin\acemblydebug/winacembly.exe"
# SUBTRACT LINK32 /verbose /nodefaultlib
LINK32_FLAGS=/nologo /version:4.31 /entry:"WinMainCRTStartup"\
 /subsystem:windows /profile /map:"$(INTDIR)/winacembly.map" /debug\
 /machine:I386 /nodefaultlib:"libcmtd.lib" /out:"$(OUTDIR)/winacembly.exe" 
LINK32_OBJS= \
	"$(INTDIR)\abifix.obj" \
	"$(INTDIR)\acedbprofile.obj" \
	"$(INTDIR)\acedbprofileview.obj" \
	"$(INTDIR)\acedialogs.obj" \
	"$(INTDIR)\acefiledoc.obj" \
	"$(INTDIR)\acefileview.obj" \
	"$(INTDIR)\acemblyhook.obj" \
	"$(INTDIR)\action.obj" \
	"$(INTDIR)\afxtrace.obj" \
	"$(INTDIR)\align.obj" \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\aligntools.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\basecall.obj" \
	"$(INTDIR)\basecallstat.obj" \
	"$(INTDIR)\basepad.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\blxview.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\caceprintpage.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\cgraph.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\cmapdisp.obj" \
	"$(INTDIR)\cmultiass.obj" \
	"$(INTDIR)\colcontrol.obj" \
	"$(INTDIR)\colourbox.obj" \
	"$(INTDIR)\command.obj" \
	"$(INTDIR)\dbprofile.obj" \
	"$(INTDIR)\defcpt.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\display.obj" \
	"$(INTDIR)\dnaalign.obj" \
	"$(INTDIR)\dnacpt.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dotter.obj" \
	"$(INTDIR)\dotterKarlin.obj" \
	"$(INTDIR)\drawdisp.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\embl.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\fitview.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\fmapblast.obj" \
	"$(INTDIR)\fmapcontrol.obj" \
	"$(INTDIR)\fmapfeatures.obj" \
	"$(INTDIR)\fmapgene.obj" \
	"$(INTDIR)\fmaposp.obj" \
	"$(INTDIR)\fmapsequence.obj" \
	"$(INTDIR)\fmaptrace.obj" \
	"$(INTDIR)\fontpreference.obj" \
	"$(INTDIR)\forest.obj" \
	"$(INTDIR)\forestdisplay.obj" \
	"$(INTDIR)\fpdisp.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\fullscrollview.obj" \
	"$(INTDIR)\gd.obj" \
	"$(INTDIR)\geldisp.obj" \
	"$(INTDIR)\gfcode.obj" \
	"$(INTDIR)\gmapconvert.obj" \
	"$(INTDIR)\gmapdata.obj" \
	"$(INTDIR)\gmapdatacol.obj" \
	"$(INTDIR)\gmapdisp.obj" \
	"$(INTDIR)\gmapintervalcol.obj" \
	"$(INTDIR)\gmaplocuscol.obj" \
	"$(INTDIR)\gmapmarkercol.obj" \
	"$(INTDIR)\gmapphys.obj" \
	"$(INTDIR)\gmapposnegcol.obj" \
	"$(INTDIR)\gmapremarkcol.obj" \
	"$(INTDIR)\graphbox.obj" \
	"$(INTDIR)\graphcon.obj" \
	"$(INTDIR)\graphramp.obj" \
	"$(INTDIR)\graphselect.obj" \
	"$(INTDIR)\graphsub.obj" \
	"$(INTDIR)\graphview.obj" \
	"$(INTDIR)\graphwin.obj" \
	"$(INTDIR)\graphwin32_.obj" \
	"$(INTDIR)\griddisp.obj" \
	"$(INTDIR)\heap.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\hexcode.obj" \
	"$(INTDIR)\intrinsictree.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\ksetdisp.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\mach-io.obj" \
	"$(INTDIR)\mainfrm.obj" \
	"$(INTDIR)\mainpick.obj" \
	"$(INTDIR)\mapcontrol.obj" \
	"$(INTDIR)\mapscrollview.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\menu.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\messubs.obj" \
	"$(INTDIR)\metab.obj" \
	"$(INTDIR)\method.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\msgwindow.obj" \
	"$(INTDIR)\myNetwork.obj" \
	"$(INTDIR)\newkey.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\nqcpass1.obj" \
	"$(INTDIR)\nqcpass2.obj" \
	"$(INTDIR)\nqctools.obj" \
	"$(INTDIR)\o2m.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\opp.obj" \
	"$(INTDIR)\oxgriddisp.obj" \
	"$(INTDIR)\oxhomlist.obj" \
	"$(INTDIR)\pairmapdisp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\pepactivezonecol.obj" \
	"$(INTDIR)\pepdisp.obj" \
	"$(INTDIR)\pepfeaturecol.obj" \
	"$(INTDIR)\pepgraphcol.obj" \
	"$(INTDIR)\pephomolcol.obj" \
	"$(INTDIR)\pepseqcol.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\pixelfitview.obj" \
	"$(INTDIR)\pixelscrollview.obj" \
	"$(INTDIR)\plainview.obj" \
	"$(INTDIR)\plot.obj" \
	"$(INTDIR)\pmapconvert.obj" \
	"$(INTDIR)\pmapdisp.obj" \
	"$(INTDIR)\preferences.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\qbedisp.obj" \
	"$(INTDIR)\querybuild.obj" \
	"$(INTDIR)\querydisp.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\regression.obj" \
	"$(INTDIR)\seq.obj" \
	"$(INTDIR)\seqIOSCF.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\specg.obj" \
	"$(INTDIR)\splashbox.obj" \
	"$(INTDIR)\sprdctrl.obj" \
	"$(INTDIR)\sprddata.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprddisplay.obj" \
	"$(INTDIR)\sprdmap.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\status.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\textfitview.obj" \
	"$(INTDIR)\textscrollview.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\topology.obj" \
	"$(INTDIR)\tqdebug.obj" \
	"$(INTDIR)\trace.obj" \
	"$(INTDIR)\translate.obj" \
	"$(INTDIR)\treedisp.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\userinterface.obj" \
	"$(INTDIR)\viewedit.obj" \
	"$(INTDIR)\vmapdata2.obj" \
	"$(INTDIR)\vmapdisp.obj" \
	"$(INTDIR)\vmapdrag.obj" \
	"$(INTDIR)\vmapphys.obj" \
	"$(INTDIR)\win32.obj" \
	"$(INTDIR)\win32lib.obj" \
	"$(INTDIR)\win32menusubs.obj" \
	"$(INTDIR)\win32print.obj" \
	"$(INTDIR)\win32process.obj" \
	"$(INTDIR)\win32thread.obj" \
	"$(INTDIR)\win32util.obj" \
	"$(INTDIR)\winace.obj" \
	"$(INTDIR)\winace.res" \
	"$(INTDIR)\winacemail.obj" \
	"$(INTDIR)\windialogs.obj" \
	"$(INTDIR)\winfilquery.obj" \
	"$(INTDIR)\winmain.obj"

"$(OUTDIR)\winacembly.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "acedb_0"
# PROP BASE Intermediate_Dir "acedb_0"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "\acedb\bin\acemblyrel"
# PROP Intermediate_Dir "\acedb\bin\acemblyrel"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\acemblyrel
INTDIR=\acedb\bin\acemblyrel

ALL : "$(OUTDIR)\winacembly.exe"

CLEAN : 
	-@erase "$(INTDIR)\abifix.obj"
	-@erase "$(INTDIR)\acedbprofile.obj"
	-@erase "$(INTDIR)\acedbprofileview.obj"
	-@erase "$(INTDIR)\acedialogs.obj"
	-@erase "$(INTDIR)\acefiledoc.obj"
	-@erase "$(INTDIR)\acefileview.obj"
	-@erase "$(INTDIR)\acemblyhook.obj"
	-@erase "$(INTDIR)\action.obj"
	-@erase "$(INTDIR)\afxtrace.obj"
	-@erase "$(INTDIR)\align.obj"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\aligntools.obj"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\basecall.obj"
	-@erase "$(INTDIR)\basecallstat.obj"
	-@erase "$(INTDIR)\basepad.obj"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\blxview.obj"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\caceprintpage.obj"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\cgraph.obj"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\cmapdisp.obj"
	-@erase "$(INTDIR)\cmultiass.obj"
	-@erase "$(INTDIR)\colcontrol.obj"
	-@erase "$(INTDIR)\colourbox.obj"
	-@erase "$(INTDIR)\command.obj"
	-@erase "$(INTDIR)\dbprofile.obj"
	-@erase "$(INTDIR)\defcpt.obj"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\display.obj"
	-@erase "$(INTDIR)\dnaalign.obj"
	-@erase "$(INTDIR)\dnacpt.obj"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dotter.obj"
	-@erase "$(INTDIR)\dotterKarlin.obj"
	-@erase "$(INTDIR)\drawdisp.obj"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\embl.obj"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\fitview.obj"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\fmapblast.obj"
	-@erase "$(INTDIR)\fmapcontrol.obj"
	-@erase "$(INTDIR)\fmapfeatures.obj"
	-@erase "$(INTDIR)\fmapgene.obj"
	-@erase "$(INTDIR)\fmaposp.obj"
	-@erase "$(INTDIR)\fmapsequence.obj"
	-@erase "$(INTDIR)\fmaptrace.obj"
	-@erase "$(INTDIR)\fontpreference.obj"
	-@erase "$(INTDIR)\forest.obj"
	-@erase "$(INTDIR)\forestdisplay.obj"
	-@erase "$(INTDIR)\fpdisp.obj"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\fullscrollview.obj"
	-@erase "$(INTDIR)\gd.obj"
	-@erase "$(INTDIR)\geldisp.obj"
	-@erase "$(INTDIR)\gfcode.obj"
	-@erase "$(INTDIR)\gmapconvert.obj"
	-@erase "$(INTDIR)\gmapdata.obj"
	-@erase "$(INTDIR)\gmapdatacol.obj"
	-@erase "$(INTDIR)\gmapdisp.obj"
	-@erase "$(INTDIR)\gmapintervalcol.obj"
	-@erase "$(INTDIR)\gmaplocuscol.obj"
	-@erase "$(INTDIR)\gmapmarkercol.obj"
	-@erase "$(INTDIR)\gmapphys.obj"
	-@erase "$(INTDIR)\gmapposnegcol.obj"
	-@erase "$(INTDIR)\gmapremarkcol.obj"
	-@erase "$(INTDIR)\graphbox.obj"
	-@erase "$(INTDIR)\graphcon.obj"
	-@erase "$(INTDIR)\graphramp.obj"
	-@erase "$(INTDIR)\graphselect.obj"
	-@erase "$(INTDIR)\graphsub.obj"
	-@erase "$(INTDIR)\graphview.obj"
	-@erase "$(INTDIR)\graphwin.obj"
	-@erase "$(INTDIR)\graphwin32_.obj"
	-@erase "$(INTDIR)\griddisp.obj"
	-@erase "$(INTDIR)\heap.obj"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\hexcode.obj"
	-@erase "$(INTDIR)\intrinsictree.obj"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\ksetdisp.obj"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\mach-io.obj"
	-@erase "$(INTDIR)\mainfrm.obj"
	-@erase "$(INTDIR)\mainpick.obj"
	-@erase "$(INTDIR)\mapcontrol.obj"
	-@erase "$(INTDIR)\mapscrollview.obj"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\menu.obj"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messubs.obj"
	-@erase "$(INTDIR)\metab.obj"
	-@erase "$(INTDIR)\method.obj"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\msgwindow.obj"
	-@erase "$(INTDIR)\myNetwork.obj"
	-@erase "$(INTDIR)\newkey.obj"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\nqcpass1.obj"
	-@erase "$(INTDIR)\nqcpass2.obj"
	-@erase "$(INTDIR)\nqctools.obj"
	-@erase "$(INTDIR)\o2m.obj"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\opp.obj"
	-@erase "$(INTDIR)\oxgriddisp.obj"
	-@erase "$(INTDIR)\oxhomlist.obj"
	-@erase "$(INTDIR)\pairmapdisp.obj"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\pepactivezonecol.obj"
	-@erase "$(INTDIR)\pepdisp.obj"
	-@erase "$(INTDIR)\pepfeaturecol.obj"
	-@erase "$(INTDIR)\pepgraphcol.obj"
	-@erase "$(INTDIR)\pephomolcol.obj"
	-@erase "$(INTDIR)\pepseqcol.obj"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\pixelfitview.obj"
	-@erase "$(INTDIR)\pixelscrollview.obj"
	-@erase "$(INTDIR)\plainview.obj"
	-@erase "$(INTDIR)\plot.obj"
	-@erase "$(INTDIR)\pmapconvert.obj"
	-@erase "$(INTDIR)\pmapdisp.obj"
	-@erase "$(INTDIR)\preferences.obj"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\qbedisp.obj"
	-@erase "$(INTDIR)\querybuild.obj"
	-@erase "$(INTDIR)\querydisp.obj"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\regression.obj"
	-@erase "$(INTDIR)\seq.obj"
	-@erase "$(INTDIR)\seqIOSCF.obj"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\specg.obj"
	-@erase "$(INTDIR)\splashbox.obj"
	-@erase "$(INTDIR)\sprdctrl.obj"
	-@erase "$(INTDIR)\sprddata.obj"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprddisplay.obj"
	-@erase "$(INTDIR)\sprdmap.obj"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\status.obj"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\textfitview.obj"
	-@erase "$(INTDIR)\textscrollview.obj"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\topology.obj"
	-@erase "$(INTDIR)\tqdebug.obj"
	-@erase "$(INTDIR)\trace.obj"
	-@erase "$(INTDIR)\translate.obj"
	-@erase "$(INTDIR)\treedisp.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\userinterface.obj"
	-@erase "$(INTDIR)\viewedit.obj"
	-@erase "$(INTDIR)\vmapdata2.obj"
	-@erase "$(INTDIR)\vmapdisp.obj"
	-@erase "$(INTDIR)\vmapdrag.obj"
	-@erase "$(INTDIR)\vmapphys.obj"
	-@erase "$(INTDIR)\win32.obj"
	-@erase "$(INTDIR)\win32lib.obj"
	-@erase "$(INTDIR)\win32menusubs.obj"
	-@erase "$(INTDIR)\win32print.obj"
	-@erase "$(INTDIR)\win32process.obj"
	-@erase "$(INTDIR)\win32thread.obj"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\winace.obj"
	-@erase "$(INTDIR)\winace.res"
	-@erase "$(INTDIR)\winacemail.obj"
	-@erase "$(INTDIR)\windialogs.obj"
	-@erase "$(INTDIR)\winfilquery.obj"
	-@erase "$(INTDIR)\winmain.obj"
	-@erase "$(OUTDIR)\winacembly.exe"
	-@erase "$(OUTDIR)\winacembly.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "acedb" /D "WIN32" /D "_WINDOWS" /D "acedb4" /D "_AFXDLL" /D "_MBCS" /YX /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /YX /c
# SUBTRACT CPP /Fr
CPP_PROJ=/nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\acemblyrel/
CPP_SBRS=.\.
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"\acedb\bin\acemblyrel/winacembly.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/winacembly.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /map /machine:I386 /nodefaultlib:"libc.lib" /nodefaultlib:"libcmt.lib" /nodefaultlib:"msvcrtd.lib" /nodefaultlib:"libcd.lib" /nodefaultlib:"libcmtd.lib"
# SUBTRACT BASE LINK32 /verbose /nodefaultlib
# ADD LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /map /machine:I386 /nodefaultlib:"libcmt.lib" /out:"\acedb\bin\acemblyrel/winacembly.exe"
# SUBTRACT LINK32 /verbose /nodefaultlib
LINK32_FLAGS=/nologo /version:4.31 /entry:"WinMainCRTStartup"\
 /subsystem:windows /incremental:no /pdb:"$(OUTDIR)/winacembly.pdb"\
 /map:"$(INTDIR)/winacembly.map" /machine:I386 /nodefaultlib:"libcmt.lib"\
 /out:"$(OUTDIR)/winacembly.exe" 
LINK32_OBJS= \
	"$(INTDIR)\abifix.obj" \
	"$(INTDIR)\acedbprofile.obj" \
	"$(INTDIR)\acedbprofileview.obj" \
	"$(INTDIR)\acedialogs.obj" \
	"$(INTDIR)\acefiledoc.obj" \
	"$(INTDIR)\acefileview.obj" \
	"$(INTDIR)\acemblyhook.obj" \
	"$(INTDIR)\action.obj" \
	"$(INTDIR)\afxtrace.obj" \
	"$(INTDIR)\align.obj" \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\aligntools.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\basecall.obj" \
	"$(INTDIR)\basecallstat.obj" \
	"$(INTDIR)\basepad.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\blxview.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\caceprintpage.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\cgraph.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\cmapdisp.obj" \
	"$(INTDIR)\cmultiass.obj" \
	"$(INTDIR)\colcontrol.obj" \
	"$(INTDIR)\colourbox.obj" \
	"$(INTDIR)\command.obj" \
	"$(INTDIR)\dbprofile.obj" \
	"$(INTDIR)\defcpt.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\display.obj" \
	"$(INTDIR)\dnaalign.obj" \
	"$(INTDIR)\dnacpt.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dotter.obj" \
	"$(INTDIR)\dotterKarlin.obj" \
	"$(INTDIR)\drawdisp.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\embl.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\fitview.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\fmapblast.obj" \
	"$(INTDIR)\fmapcontrol.obj" \
	"$(INTDIR)\fmapfeatures.obj" \
	"$(INTDIR)\fmapgene.obj" \
	"$(INTDIR)\fmaposp.obj" \
	"$(INTDIR)\fmapsequence.obj" \
	"$(INTDIR)\fmaptrace.obj" \
	"$(INTDIR)\fontpreference.obj" \
	"$(INTDIR)\forest.obj" \
	"$(INTDIR)\forestdisplay.obj" \
	"$(INTDIR)\fpdisp.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\fullscrollview.obj" \
	"$(INTDIR)\gd.obj" \
	"$(INTDIR)\geldisp.obj" \
	"$(INTDIR)\gfcode.obj" \
	"$(INTDIR)\gmapconvert.obj" \
	"$(INTDIR)\gmapdata.obj" \
	"$(INTDIR)\gmapdatacol.obj" \
	"$(INTDIR)\gmapdisp.obj" \
	"$(INTDIR)\gmapintervalcol.obj" \
	"$(INTDIR)\gmaplocuscol.obj" \
	"$(INTDIR)\gmapmarkercol.obj" \
	"$(INTDIR)\gmapphys.obj" \
	"$(INTDIR)\gmapposnegcol.obj" \
	"$(INTDIR)\gmapremarkcol.obj" \
	"$(INTDIR)\graphbox.obj" \
	"$(INTDIR)\graphcon.obj" \
	"$(INTDIR)\graphramp.obj" \
	"$(INTDIR)\graphselect.obj" \
	"$(INTDIR)\graphsub.obj" \
	"$(INTDIR)\graphview.obj" \
	"$(INTDIR)\graphwin.obj" \
	"$(INTDIR)\graphwin32_.obj" \
	"$(INTDIR)\griddisp.obj" \
	"$(INTDIR)\heap.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\hexcode.obj" \
	"$(INTDIR)\intrinsictree.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\ksetdisp.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\mach-io.obj" \
	"$(INTDIR)\mainfrm.obj" \
	"$(INTDIR)\mainpick.obj" \
	"$(INTDIR)\mapcontrol.obj" \
	"$(INTDIR)\mapscrollview.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\menu.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\messubs.obj" \
	"$(INTDIR)\metab.obj" \
	"$(INTDIR)\method.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\msgwindow.obj" \
	"$(INTDIR)\myNetwork.obj" \
	"$(INTDIR)\newkey.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\nqcpass1.obj" \
	"$(INTDIR)\nqcpass2.obj" \
	"$(INTDIR)\nqctools.obj" \
	"$(INTDIR)\o2m.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\opp.obj" \
	"$(INTDIR)\oxgriddisp.obj" \
	"$(INTDIR)\oxhomlist.obj" \
	"$(INTDIR)\pairmapdisp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\pepactivezonecol.obj" \
	"$(INTDIR)\pepdisp.obj" \
	"$(INTDIR)\pepfeaturecol.obj" \
	"$(INTDIR)\pepgraphcol.obj" \
	"$(INTDIR)\pephomolcol.obj" \
	"$(INTDIR)\pepseqcol.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\pixelfitview.obj" \
	"$(INTDIR)\pixelscrollview.obj" \
	"$(INTDIR)\plainview.obj" \
	"$(INTDIR)\plot.obj" \
	"$(INTDIR)\pmapconvert.obj" \
	"$(INTDIR)\pmapdisp.obj" \
	"$(INTDIR)\preferences.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\qbedisp.obj" \
	"$(INTDIR)\querybuild.obj" \
	"$(INTDIR)\querydisp.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\regression.obj" \
	"$(INTDIR)\seq.obj" \
	"$(INTDIR)\seqIOSCF.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\specg.obj" \
	"$(INTDIR)\splashbox.obj" \
	"$(INTDIR)\sprdctrl.obj" \
	"$(INTDIR)\sprddata.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprddisplay.obj" \
	"$(INTDIR)\sprdmap.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\status.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\textfitview.obj" \
	"$(INTDIR)\textscrollview.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\topology.obj" \
	"$(INTDIR)\tqdebug.obj" \
	"$(INTDIR)\trace.obj" \
	"$(INTDIR)\translate.obj" \
	"$(INTDIR)\treedisp.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\userinterface.obj" \
	"$(INTDIR)\viewedit.obj" \
	"$(INTDIR)\vmapdata2.obj" \
	"$(INTDIR)\vmapdisp.obj" \
	"$(INTDIR)\vmapdrag.obj" \
	"$(INTDIR)\vmapphys.obj" \
	"$(INTDIR)\win32.obj" \
	"$(INTDIR)\win32lib.obj" \
	"$(INTDIR)\win32menusubs.obj" \
	"$(INTDIR)\win32print.obj" \
	"$(INTDIR)\win32process.obj" \
	"$(INTDIR)\win32thread.obj" \
	"$(INTDIR)\win32util.obj" \
	"$(INTDIR)\winace.obj" \
	"$(INTDIR)\winace.res" \
	"$(INTDIR)\winacemail.obj" \
	"$(INTDIR)\windialogs.obj" \
	"$(INTDIR)\winfilquery.obj" \
	"$(INTDIR)\winmain.obj"

"$(OUTDIR)\winacembly.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "acedb___"
# PROP BASE Intermediate_Dir "acedb___"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "\acedb\bin\aceclientdebug"
# PROP Intermediate_Dir "\acedb\bin\aceclientdebug"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\aceclientdebug
INTDIR=\acedb\bin\aceclientdebug

ALL : "$(OUTDIR)\aceclientd.exe" "$(OUTDIR)\aceclient.bsc"

CLEAN : 
	-@erase "$(INTDIR)\aceclient.obj"
	-@erase "$(INTDIR)\aceclient.sbr"
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\arraysub.sbr"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\bump.sbr"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\call.sbr"
	-@erase "$(INTDIR)\dceclientlib.obj"
	-@erase "$(INTDIR)\dceclientlib.sbr"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\dict.sbr"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filqng.sbr"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\filsubs.sbr"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freeout.sbr"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\freesubs.sbr"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messclean.sbr"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\prefsubs.sbr"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\randsubs.sbr"
	-@erase "$(INTDIR)\rpcace_c.obj"
	-@erase "$(INTDIR)\rpcace_c.sbr"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\stdafx.sbr"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\timesubs.sbr"
	-@erase "$(INTDIR)\vc40.idb"
	-@erase "$(INTDIR)\vc40.pdb"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\win32util.sbr"
	-@erase "$(OUTDIR)\aceclient.bsc"
	-@erase "$(OUTDIR)\aceclientd.exe"
	-@erase "$(OUTDIR)\aceclientd.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR /c
# SUBTRACT CPP /YX
CPP_PROJ=/nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\aceclientdebug/
CPP_SBRS=\acedb\bin\aceclientdebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceclientdebug/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\windebug/winace.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\aceclientdebug/aceclient.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/aceclient.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\aceclient.sbr" \
	"$(INTDIR)\arraysub.sbr" \
	"$(INTDIR)\bump.sbr" \
	"$(INTDIR)\call.sbr" \
	"$(INTDIR)\dceclientlib.sbr" \
	"$(INTDIR)\dict.sbr" \
	"$(INTDIR)\filqng.sbr" \
	"$(INTDIR)\filsubs.sbr" \
	"$(INTDIR)\freeout.sbr" \
	"$(INTDIR)\freesubs.sbr" \
	"$(INTDIR)\messclean.sbr" \
	"$(INTDIR)\prefsubs.sbr" \
	"$(INTDIR)\randsubs.sbr" \
	"$(INTDIR)\rpcace_c.sbr" \
	"$(INTDIR)\stdafx.sbr" \
	"$(INTDIR)\timesubs.sbr" \
	"$(INTDIR)\win32util.sbr"

"$(OUTDIR)\aceclient.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 /nologo /version:4.31 /subsystem:windows /profile /map /debug /machine:I386 /out:"\acedb\bin\windebug/winace.exe"
# SUBTRACT BASE LINK32 /verbose /nodefaultlib
# ADD LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386 /out:"\acedb\bin\aceclientdebug/aceclientd.exe"
LINK32_FLAGS=rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib\
 winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib\
 uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console\
 /profile /map:"$(INTDIR)/aceclientd.map" /debug /machine:I386\
 /out:"$(OUTDIR)/aceclientd.exe" 
LINK32_OBJS= \
	"$(INTDIR)\aceclient.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\dceclientlib.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\rpcace_c.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\aceclientd.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "acedb__1"
# PROP BASE Intermediate_Dir "acedb__1"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "\acedb\bin\aceserverdebug"
# PROP Intermediate_Dir "\acedb\bin\aceserverdebug"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\aceserverdebug
INTDIR=\acedb\bin\aceserverdebug

ALL : "$(OUTDIR)\aceserverd.exe" "$(OUTDIR)\aceserver.bsc"

CLEAN : 
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\aceserver.obj"
	-@erase "$(INTDIR)\aceserver.sbr"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\alignment.sbr"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\arraysub.sbr"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\asubs.sbr"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\banner.sbr"
	-@erase "$(INTDIR)\basepad.obj"
	-@erase "$(INTDIR)\basepad.sbr"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\biblio.sbr"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\blocksub.sbr"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bsdumps.sbr"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bssubs.sbr"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstools.sbr"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bstree.sbr"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bsubs.sbr"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\bump.sbr"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\call.sbr"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\check.sbr"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\chrono.sbr"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\class.sbr"
	-@erase "$(INTDIR)\command.obj"
	-@erase "$(INTDIR)\command.sbr"
	-@erase "$(INTDIR)\dceserverlib.obj"
	-@erase "$(INTDIR)\dceserverlib.sbr"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\dict.sbr"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\disknew.sbr"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dnasubs.sbr"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\dump.sbr"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filqng.sbr"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\filsubs.sbr"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\flag.sbr"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freeout.sbr"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\freesubs.sbr"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\help.sbr"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keyset.sbr"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\keysetdump.sbr"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexalpha.sbr"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs.sbr"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\lexsubs4.sbr"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\linkdate.sbr"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\longtext.sbr"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\matchtable.sbr"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messclean.sbr"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\model.sbr"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\nicedump.sbr"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\objcache.sbr"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\oldhelp.sbr"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\parse.sbr"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\peptide.sbr"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\picksubs.sbr"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\prefsubs.sbr"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\queryexe.sbr"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\quovadis.sbr"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\randsubs.sbr"
	-@erase "$(INTDIR)\rpcace_s.obj"
	-@erase "$(INTDIR)\rpcace_s.sbr"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\session.sbr"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprddef.sbr"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\sprdop.sbr"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\stdafx.sbr"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\sysclass.sbr"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\table.sbr"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\tags.sbr"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\timesubs.sbr"
	-@erase "$(INTDIR)\vc40.idb"
	-@erase "$(INTDIR)\vc40.pdb"
	-@erase "$(INTDIR)\win32dcelib.obj"
	-@erase "$(INTDIR)\win32dcelib.sbr"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\win32util.sbr"
	-@erase "$(OUTDIR)\aceserver.bsc"
	-@erase "$(OUTDIR)\aceserverd.exe"
	-@erase "$(OUTDIR)\aceserverd.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR /c
# SUBTRACT CPP /YX
CPP_PROJ=/nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\aceserverdebug/
CPP_SBRS=\acedb\bin\aceserverdebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceserverdebug/aceserver.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\windebug/winace.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\aceserverdebug/aceserver.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/aceserver.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\aceserver.sbr" \
	"$(INTDIR)\alignment.sbr" \
	"$(INTDIR)\arraysub.sbr" \
	"$(INTDIR)\asubs.sbr" \
	"$(INTDIR)\banner.sbr" \
	"$(INTDIR)\basepad.sbr" \
	"$(INTDIR)\biblio.sbr" \
	"$(INTDIR)\blocksub.sbr" \
	"$(INTDIR)\bsdumps.sbr" \
	"$(INTDIR)\bssubs.sbr" \
	"$(INTDIR)\bstools.sbr" \
	"$(INTDIR)\bstree.sbr" \
	"$(INTDIR)\bsubs.sbr" \
	"$(INTDIR)\bump.sbr" \
	"$(INTDIR)\call.sbr" \
	"$(INTDIR)\check.sbr" \
	"$(INTDIR)\chrono.sbr" \
	"$(INTDIR)\class.sbr" \
	"$(INTDIR)\command.sbr" \
	"$(INTDIR)\dceserverlib.sbr" \
	"$(INTDIR)\dict.sbr" \
	"$(INTDIR)\disknew.sbr" \
	"$(INTDIR)\dnasubs.sbr" \
	"$(INTDIR)\dump.sbr" \
	"$(INTDIR)\filqng.sbr" \
	"$(INTDIR)\filsubs.sbr" \
	"$(INTDIR)\flag.sbr" \
	"$(INTDIR)\freeout.sbr" \
	"$(INTDIR)\freesubs.sbr" \
	"$(INTDIR)\help.sbr" \
	"$(INTDIR)\keyset.sbr" \
	"$(INTDIR)\keysetdump.sbr" \
	"$(INTDIR)\lexalpha.sbr" \
	"$(INTDIR)\lexsubs.sbr" \
	"$(INTDIR)\lexsubs4.sbr" \
	"$(INTDIR)\linkdate.sbr" \
	"$(INTDIR)\longtext.sbr" \
	"$(INTDIR)\matchtable.sbr" \
	"$(INTDIR)\messclean.sbr" \
	"$(INTDIR)\model.sbr" \
	"$(INTDIR)\nicedump.sbr" \
	"$(INTDIR)\objcache.sbr" \
	"$(INTDIR)\oldhelp.sbr" \
	"$(INTDIR)\parse.sbr" \
	"$(INTDIR)\peptide.sbr" \
	"$(INTDIR)\picksubs.sbr" \
	"$(INTDIR)\prefsubs.sbr" \
	"$(INTDIR)\queryexe.sbr" \
	"$(INTDIR)\quovadis.sbr" \
	"$(INTDIR)\randsubs.sbr" \
	"$(INTDIR)\rpcace_s.sbr" \
	"$(INTDIR)\session.sbr" \
	"$(INTDIR)\sprddef.sbr" \
	"$(INTDIR)\sprdop.sbr" \
	"$(INTDIR)\stdafx.sbr" \
	"$(INTDIR)\sysclass.sbr" \
	"$(INTDIR)\table.sbr" \
	"$(INTDIR)\tags.sbr" \
	"$(INTDIR)\timesubs.sbr" \
	"$(INTDIR)\win32dcelib.sbr" \
	"$(INTDIR)\win32util.sbr"

"$(OUTDIR)\aceserver.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 /nologo /version:4.31 /subsystem:windows /profile /map /debug /machine:I386 /out:"\acedb\bin\windebug/winace.exe"
# SUBTRACT BASE LINK32 /verbose /nodefaultlib
# ADD LINK32 rpcrt4.lib rpcns4.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386 /out:"\acedb\bin\aceserverdebug/aceserverd.exe"
LINK32_FLAGS=rpcrt4.lib rpcns4.lib /nologo /version:4.31 /subsystem:console\
 /profile /map:"$(INTDIR)/aceserverd.map" /debug /machine:I386\
 /out:"$(OUTDIR)/aceserverd.exe" 
LINK32_OBJS= \
	"$(INTDIR)\aceserver.obj" \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\basepad.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\command.obj" \
	"$(INTDIR)\dceserverlib.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\rpcace_s.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32dcelib.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\aceserverd.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "acedb__2"
# PROP BASE Intermediate_Dir "acedb__2"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "\acedb\bin\aceserverrel"
# PROP Intermediate_Dir "\acedb\bin\aceserverrel"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\aceserverrel
INTDIR=\acedb\bin\aceserverrel

ALL : "$(OUTDIR)\aceserver.exe"

CLEAN : 
	-@erase "$(INTDIR)\aceserver.obj"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\basepad.obj"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\command.obj"
	-@erase "$(INTDIR)\dceserverlib.obj"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\rpcace_s.obj"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\win32dcelib.obj"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(OUTDIR)\aceserver.exe"
	-@erase "$(OUTDIR)\aceserver.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /YX /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /c
# SUBTRACT CPP /Fr /YX
CPP_PROJ=/nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\aceserverrel/
CPP_SBRS=.\.
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceserverrel/aceserver.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\winrel/winace.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\aceserverrel/aceserver.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/aceserver.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /map /machine:I386 /nodefaultlib:"libcmt.lib" /out:"\acedb\bin\winrel/winace.exe"
# SUBTRACT BASE LINK32 /verbose /nodefaultlib
# ADD LINK32 rpcrt4.lib rpcns4.lib /nologo /version:4.31 /subsystem:console /map /machine:I386 /out:"\acedb\bin\aceserverrel/aceserver.exe"
# SUBTRACT LINK32 /pdb:none
LINK32_FLAGS=rpcrt4.lib rpcns4.lib /nologo /version:4.31 /subsystem:console\
 /incremental:no /pdb:"$(OUTDIR)/aceserver.pdb" /map:"$(INTDIR)/aceserver.map"\
 /machine:I386 /out:"$(OUTDIR)/aceserver.exe" 
LINK32_OBJS= \
	"$(INTDIR)\aceserver.obj" \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\basepad.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\command.obj" \
	"$(INTDIR)\dceserverlib.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\rpcace_s.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32dcelib.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\aceserver.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "acedb__0"
# PROP BASE Intermediate_Dir "acedb__0"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "\acedb\bin\aceclientrel"
# PROP Intermediate_Dir "\acedb\bin\aceclientrel"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\aceclientrel
INTDIR=\acedb\bin\aceclientrel

ALL : "$(OUTDIR)\aceclient.exe" ".\rpcace_c.c" ".\rpcace.h"

CLEAN : 
	-@erase "$(INTDIR)\aceclient.obj"
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\dceclientlib.obj"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\rpcace_c.obj"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(OUTDIR)\aceclient.exe"
	-@erase "$(OUTDIR)\aceclient.map"
	-@erase ".\rpcace.h"
	-@erase ".\rpcace_c.c"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /YX /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /c
# SUBTRACT CPP /Fr /YX
CPP_PROJ=/nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\aceclientrel/
CPP_SBRS=.\.
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceclientrel/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\winrel/winace.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\aceclientrel/aceclient.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/aceclient.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 /nologo /version:4.31 /entry:"WinMainCRTStartup" /subsystem:windows /map /machine:I386 /nodefaultlib:"libcmt.lib" /out:"\acedb\bin\winrel/winace.exe"
# SUBTRACT BASE LINK32 /verbose /nodefaultlib
# ADD LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /map /machine:I386 /out:"\acedb\bin\aceclientrel/aceclient.exe"
# SUBTRACT LINK32 /pdb:none
LINK32_FLAGS=rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib\
 winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib\
 uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console\
 /incremental:no /pdb:"$(OUTDIR)/aceclient.pdb" /map:"$(INTDIR)/aceclient.map"\
 /machine:I386 /out:"$(OUTDIR)/aceclient.exe" 
LINK32_OBJS= \
	"$(INTDIR)\aceclient.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\dceclientlib.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\rpcace_c.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\aceclient.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "acedb___"
# PROP BASE Intermediate_Dir "acedb___"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "\acedb\bin\netclientdebug"
# PROP Intermediate_Dir "\acedb\bin\netclientdebug"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\netclientdebug
INTDIR=\acedb\bin\netclientdebug

ALL : "$(OUTDIR)\netclientd.exe" "$(OUTDIR)\netclient.bsc"

CLEAN : 
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\arraysub.sbr"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\bump.sbr"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\call.sbr"
	-@erase "$(INTDIR)\dceclientlib.obj"
	-@erase "$(INTDIR)\dceclientlib.sbr"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\dict.sbr"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filqng.sbr"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\filsubs.sbr"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freeout.sbr"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\freesubs.sbr"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messclean.sbr"
	-@erase "$(INTDIR)\netclient.obj"
	-@erase "$(INTDIR)\netclient.sbr"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\prefsubs.sbr"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\randsubs.sbr"
	-@erase "$(INTDIR)\rpcace_c.obj"
	-@erase "$(INTDIR)\rpcace_c.sbr"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\stdafx.sbr"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\timesubs.sbr"
	-@erase "$(INTDIR)\vc40.idb"
	-@erase "$(INTDIR)\vc40.pdb"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\win32util.sbr"
	-@erase "$(OUTDIR)\netclient.bsc"
	-@erase "$(OUTDIR)\netclientd.exe"
	-@erase "$(OUTDIR)\netclientd.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR /c
# SUBTRACT CPP /YX
CPP_PROJ=/nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\netclientdebug/
CPP_SBRS=\acedb\bin\netclientdebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /fo"\acedb\bin\aceclientdebug/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceclientdebug/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\aceclientdebug/aceclient.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\netclientdebug/netclient.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/netclient.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\arraysub.sbr" \
	"$(INTDIR)\bump.sbr" \
	"$(INTDIR)\call.sbr" \
	"$(INTDIR)\dceclientlib.sbr" \
	"$(INTDIR)\dict.sbr" \
	"$(INTDIR)\filqng.sbr" \
	"$(INTDIR)\filsubs.sbr" \
	"$(INTDIR)\freeout.sbr" \
	"$(INTDIR)\freesubs.sbr" \
	"$(INTDIR)\messclean.sbr" \
	"$(INTDIR)\netclient.sbr" \
	"$(INTDIR)\prefsubs.sbr" \
	"$(INTDIR)\randsubs.sbr" \
	"$(INTDIR)\rpcace_c.sbr" \
	"$(INTDIR)\stdafx.sbr" \
	"$(INTDIR)\timesubs.sbr" \
	"$(INTDIR)\win32util.sbr"

"$(OUTDIR)\netclient.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386 /out:"\acedb\bin\aceclientdebug/aceclient.exe"
# ADD LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386 /out:"\acedb\bin\netclientdebug/netclientd.exe"
LINK32_FLAGS=rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib\
 winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib\
 uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console\
 /profile /map:"$(INTDIR)/netclientd.map" /debug /machine:I386\
 /out:"$(OUTDIR)/netclientd.exe" 
LINK32_OBJS= \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\dceclientlib.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\netclient.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\rpcace_c.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\netclientd.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "acedb__0"
# PROP BASE Intermediate_Dir "acedb__0"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "\acedb\bin\netclientrel"
# PROP Intermediate_Dir "\acedb\bin\netclientrel"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\netclientrel
INTDIR=\acedb\bin\netclientrel

ALL : "$(OUTDIR)\netclient.exe" ".\rpcace_c.c" ".\rpcace.h"

CLEAN : 
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\dceclientlib.obj"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\netclient.obj"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\rpcace_c.obj"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(OUTDIR)\netclient.exe"
	-@erase "$(OUTDIR)\netclient.map"
	-@erase ".\rpcace.h"
	-@erase ".\rpcace_c.c"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /c
# SUBTRACT BASE CPP /Fr /YX
# ADD CPP /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /c
# SUBTRACT CPP /Fr /YX
CPP_PROJ=/nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\netclientrel/
CPP_SBRS=.\.
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /fo"\acedb\bin\aceclientrel/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceclientrel/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\aceclientrel/aceclient.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\netclientrel/netclient.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/netclient.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /map /machine:I386 /out:"\acedb\bin\aceclientrel/aceclient.exe"
# SUBTRACT BASE LINK32 /pdb:none
# ADD LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /map /machine:I386 /out:"\acedb\bin\netclientrel/netclient.exe"
# SUBTRACT LINK32 /pdb:none
LINK32_FLAGS=rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib\
 winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib\
 uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console\
 /incremental:no /pdb:"$(OUTDIR)/netclient.pdb" /map:"$(INTDIR)/netclient.map"\
 /machine:I386 /out:"$(OUTDIR)/netclient.exe" 
LINK32_OBJS= \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\dceclientlib.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\netclient.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\rpcace_c.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\netclient.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "acedb___"
# PROP BASE Intermediate_Dir "acedb___"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "\acedb\bin\wintacedebug"
# PROP Intermediate_Dir "\acedb\bin\wintacedebug"
# PROP Target_Dir ""
OUTDIR=\acedb\bin\wintacedebug
INTDIR=\acedb\bin\wintacedebug

ALL : "$(OUTDIR)\wintaced.exe" "$(OUTDIR)\wintaced.bsc"

CLEAN : 
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\alignment.sbr"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\arraysub.sbr"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\asubs.sbr"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\banner.sbr"
	-@erase "$(INTDIR)\basepad.obj"
	-@erase "$(INTDIR)\basepad.sbr"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\biblio.sbr"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\blocksub.sbr"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bsdumps.sbr"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bssubs.sbr"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstools.sbr"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bstree.sbr"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bsubs.sbr"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\bump.sbr"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\call.sbr"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\check.sbr"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\chrono.sbr"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\class.sbr"
	-@erase "$(INTDIR)\command.obj"
	-@erase "$(INTDIR)\command.sbr"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\dict.sbr"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\disknew.sbr"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dnasubs.sbr"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\dump.sbr"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filqng.sbr"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\filsubs.sbr"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\flag.sbr"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freeout.sbr"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\freesubs.sbr"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\help.sbr"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keyset.sbr"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\keysetdump.sbr"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexalpha.sbr"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs.sbr"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\lexsubs4.sbr"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\linkdate.sbr"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\longtext.sbr"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\matchtable.sbr"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\messclean.sbr"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\model.sbr"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\nicedump.sbr"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\objcache.sbr"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\oldhelp.sbr"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\parse.sbr"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\peptide.sbr"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\picksubs.sbr"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\prefsubs.sbr"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\queryexe.sbr"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\quovadis.sbr"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\randsubs.sbr"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\session.sbr"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprddef.sbr"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\sprdop.sbr"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\stdafx.sbr"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\sysclass.sbr"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\table.sbr"
	-@erase "$(INTDIR)\tacemain.obj"
	-@erase "$(INTDIR)\tacemain.sbr"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\tags.sbr"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\timesubs.sbr"
	-@erase "$(INTDIR)\vc40.idb"
	-@erase "$(INTDIR)\vc40.pdb"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(INTDIR)\win32util.sbr"
	-@erase "$(OUTDIR)\wintaced.bsc"
	-@erase "$(OUTDIR)\wintaced.exe"
	-@erase "$(OUTDIR)\wintaced.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR /c
# SUBTRACT CPP /YX
CPP_PROJ=/nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=\acedb\bin\wintacedebug/
CPP_SBRS=\acedb\bin\wintacedebug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x409 /fo"\acedb\bin\aceclientdebug/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceclientdebug/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\aceclientdebug/aceclient.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\wintacedebug/wintaced.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/wintaced.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\alignment.sbr" \
	"$(INTDIR)\arraysub.sbr" \
	"$(INTDIR)\asubs.sbr" \
	"$(INTDIR)\banner.sbr" \
	"$(INTDIR)\basepad.sbr" \
	"$(INTDIR)\biblio.sbr" \
	"$(INTDIR)\blocksub.sbr" \
	"$(INTDIR)\bsdumps.sbr" \
	"$(INTDIR)\bssubs.sbr" \
	"$(INTDIR)\bstools.sbr" \
	"$(INTDIR)\bstree.sbr" \
	"$(INTDIR)\bsubs.sbr" \
	"$(INTDIR)\bump.sbr" \
	"$(INTDIR)\call.sbr" \
	"$(INTDIR)\check.sbr" \
	"$(INTDIR)\chrono.sbr" \
	"$(INTDIR)\class.sbr" \
	"$(INTDIR)\command.sbr" \
	"$(INTDIR)\dict.sbr" \
	"$(INTDIR)\disknew.sbr" \
	"$(INTDIR)\dnasubs.sbr" \
	"$(INTDIR)\dump.sbr" \
	"$(INTDIR)\filqng.sbr" \
	"$(INTDIR)\filsubs.sbr" \
	"$(INTDIR)\flag.sbr" \
	"$(INTDIR)\freeout.sbr" \
	"$(INTDIR)\freesubs.sbr" \
	"$(INTDIR)\help.sbr" \
	"$(INTDIR)\keyset.sbr" \
	"$(INTDIR)\keysetdump.sbr" \
	"$(INTDIR)\lexalpha.sbr" \
	"$(INTDIR)\lexsubs.sbr" \
	"$(INTDIR)\lexsubs4.sbr" \
	"$(INTDIR)\linkdate.sbr" \
	"$(INTDIR)\longtext.sbr" \
	"$(INTDIR)\matchtable.sbr" \
	"$(INTDIR)\messclean.sbr" \
	"$(INTDIR)\model.sbr" \
	"$(INTDIR)\nicedump.sbr" \
	"$(INTDIR)\objcache.sbr" \
	"$(INTDIR)\oldhelp.sbr" \
	"$(INTDIR)\parse.sbr" \
	"$(INTDIR)\peptide.sbr" \
	"$(INTDIR)\picksubs.sbr" \
	"$(INTDIR)\prefsubs.sbr" \
	"$(INTDIR)\queryexe.sbr" \
	"$(INTDIR)\quovadis.sbr" \
	"$(INTDIR)\randsubs.sbr" \
	"$(INTDIR)\session.sbr" \
	"$(INTDIR)\sprddef.sbr" \
	"$(INTDIR)\sprdop.sbr" \
	"$(INTDIR)\stdafx.sbr" \
	"$(INTDIR)\sysclass.sbr" \
	"$(INTDIR)\table.sbr" \
	"$(INTDIR)\tacemain.sbr" \
	"$(INTDIR)\tags.sbr" \
	"$(INTDIR)\timesubs.sbr" \
	"$(INTDIR)\win32util.sbr"

"$(OUTDIR)\wintaced.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386 /out:"\acedb\bin\netclientdebug/netclient.exe"
# ADD LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386 /out:"\acedb\bin\wintacedebug/wintaced.exe"
LINK32_FLAGS=rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib\
 winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib\
 uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console\
 /profile /map:"$(INTDIR)/wintaced.map" /debug /machine:I386\
 /out:"$(OUTDIR)/wintaced.exe" 
LINK32_OBJS= \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\basepad.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\command.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tacemain.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\wintaced.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "acedb__0"
# PROP BASE Intermediate_Dir "acedb__0"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "D:\acedb\bin\wintacerel"
# PROP Intermediate_Dir "D:\acedb\bin\wintacerel"
# PROP Target_Dir ""
OUTDIR=D:\acedb\bin\wintacerel
INTDIR=D:\acedb\bin\wintacerel

ALL : "$(OUTDIR)\wintace.exe"

CLEAN : 
	-@erase "$(INTDIR)\Acedb.pch"
	-@erase "$(INTDIR)\alignment.obj"
	-@erase "$(INTDIR)\arraysub.obj"
	-@erase "$(INTDIR)\asubs.obj"
	-@erase "$(INTDIR)\banner.obj"
	-@erase "$(INTDIR)\basepad.obj"
	-@erase "$(INTDIR)\biblio.obj"
	-@erase "$(INTDIR)\blocksub.obj"
	-@erase "$(INTDIR)\bsdumps.obj"
	-@erase "$(INTDIR)\bssubs.obj"
	-@erase "$(INTDIR)\bstools.obj"
	-@erase "$(INTDIR)\bstree.obj"
	-@erase "$(INTDIR)\bsubs.obj"
	-@erase "$(INTDIR)\bump.obj"
	-@erase "$(INTDIR)\call.obj"
	-@erase "$(INTDIR)\check.obj"
	-@erase "$(INTDIR)\chrono.obj"
	-@erase "$(INTDIR)\class.obj"
	-@erase "$(INTDIR)\command.obj"
	-@erase "$(INTDIR)\dict.obj"
	-@erase "$(INTDIR)\disknew.obj"
	-@erase "$(INTDIR)\dnasubs.obj"
	-@erase "$(INTDIR)\dump.obj"
	-@erase "$(INTDIR)\filqng.obj"
	-@erase "$(INTDIR)\filsubs.obj"
	-@erase "$(INTDIR)\flag.obj"
	-@erase "$(INTDIR)\freeout.obj"
	-@erase "$(INTDIR)\freesubs.obj"
	-@erase "$(INTDIR)\help.obj"
	-@erase "$(INTDIR)\keyset.obj"
	-@erase "$(INTDIR)\keysetdump.obj"
	-@erase "$(INTDIR)\lexalpha.obj"
	-@erase "$(INTDIR)\lexsubs.obj"
	-@erase "$(INTDIR)\lexsubs4.obj"
	-@erase "$(INTDIR)\linkdate.obj"
	-@erase "$(INTDIR)\longtext.obj"
	-@erase "$(INTDIR)\matchtable.obj"
	-@erase "$(INTDIR)\messclean.obj"
	-@erase "$(INTDIR)\model.obj"
	-@erase "$(INTDIR)\nicedump.obj"
	-@erase "$(INTDIR)\objcache.obj"
	-@erase "$(INTDIR)\oldhelp.obj"
	-@erase "$(INTDIR)\parse.obj"
	-@erase "$(INTDIR)\peptide.obj"
	-@erase "$(INTDIR)\picksubs.obj"
	-@erase "$(INTDIR)\prefsubs.obj"
	-@erase "$(INTDIR)\queryexe.obj"
	-@erase "$(INTDIR)\quovadis.obj"
	-@erase "$(INTDIR)\randsubs.obj"
	-@erase "$(INTDIR)\session.obj"
	-@erase "$(INTDIR)\sprddef.obj"
	-@erase "$(INTDIR)\sprdop.obj"
	-@erase "$(INTDIR)\stdafx.obj"
	-@erase "$(INTDIR)\sysclass.obj"
	-@erase "$(INTDIR)\table.obj"
	-@erase "$(INTDIR)\tacemain.obj"
	-@erase "$(INTDIR)\tags.obj"
	-@erase "$(INTDIR)\timesubs.obj"
	-@erase "$(INTDIR)\win32util.obj"
	-@erase "$(OUTDIR)\wintace.exe"
	-@erase "$(OUTDIR)\wintace.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /c
# SUBTRACT BASE CPP /Fr /YX
# ADD CPP /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /c
# SUBTRACT CPP /Fr /YX
CPP_PROJ=/nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=D:\acedb\bin\wintacerel/
CPP_SBRS=.\.
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x409 /fo"\acedb\bin\aceclientrel/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG"
# ADD RSC /l 0x409 /fo"\acedb\bin\aceclientrel/aceclient.res" /i "\acedb\win32h" /i "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo /o"\acedb\bin\aceclientrel/aceclient.bsc"
# ADD BSC32 /nologo /o"\acedb\bin\wintacerel/wintace.bsc"
BSC32_FLAGS=/nologo /o"$(OUTDIR)/wintace.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /map /machine:I386 /out:"\acedb\bin\netclientrel/netclient.exe"
# SUBTRACT BASE LINK32 /pdb:none
# ADD LINK32 rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /map /machine:I386 /out:"\acedb\bin\wintacerel\wintace.exe"
# SUBTRACT LINK32 /pdb:none
LINK32_FLAGS=rpcrt4.lib rpcns4.lib kernel32.lib user32.lib gdi32.lib\
 winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib\
 uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console\
 /incremental:no /pdb:"$(OUTDIR)/wintace.pdb" /map:"$(INTDIR)/wintace.map"\
 /machine:I386 /out:"$(OUTDIR)/wintace.exe" 
LINK32_OBJS= \
	"$(INTDIR)\alignment.obj" \
	"$(INTDIR)\arraysub.obj" \
	"$(INTDIR)\asubs.obj" \
	"$(INTDIR)\banner.obj" \
	"$(INTDIR)\basepad.obj" \
	"$(INTDIR)\biblio.obj" \
	"$(INTDIR)\blocksub.obj" \
	"$(INTDIR)\bsdumps.obj" \
	"$(INTDIR)\bssubs.obj" \
	"$(INTDIR)\bstools.obj" \
	"$(INTDIR)\bstree.obj" \
	"$(INTDIR)\bsubs.obj" \
	"$(INTDIR)\bump.obj" \
	"$(INTDIR)\call.obj" \
	"$(INTDIR)\check.obj" \
	"$(INTDIR)\chrono.obj" \
	"$(INTDIR)\class.obj" \
	"$(INTDIR)\command.obj" \
	"$(INTDIR)\dict.obj" \
	"$(INTDIR)\disknew.obj" \
	"$(INTDIR)\dnasubs.obj" \
	"$(INTDIR)\dump.obj" \
	"$(INTDIR)\filqng.obj" \
	"$(INTDIR)\filsubs.obj" \
	"$(INTDIR)\flag.obj" \
	"$(INTDIR)\freeout.obj" \
	"$(INTDIR)\freesubs.obj" \
	"$(INTDIR)\help.obj" \
	"$(INTDIR)\keyset.obj" \
	"$(INTDIR)\keysetdump.obj" \
	"$(INTDIR)\lexalpha.obj" \
	"$(INTDIR)\lexsubs.obj" \
	"$(INTDIR)\lexsubs4.obj" \
	"$(INTDIR)\linkdate.obj" \
	"$(INTDIR)\longtext.obj" \
	"$(INTDIR)\matchtable.obj" \
	"$(INTDIR)\messclean.obj" \
	"$(INTDIR)\model.obj" \
	"$(INTDIR)\nicedump.obj" \
	"$(INTDIR)\objcache.obj" \
	"$(INTDIR)\oldhelp.obj" \
	"$(INTDIR)\parse.obj" \
	"$(INTDIR)\peptide.obj" \
	"$(INTDIR)\picksubs.obj" \
	"$(INTDIR)\prefsubs.obj" \
	"$(INTDIR)\queryexe.obj" \
	"$(INTDIR)\quovadis.obj" \
	"$(INTDIR)\randsubs.obj" \
	"$(INTDIR)\session.obj" \
	"$(INTDIR)\sprddef.obj" \
	"$(INTDIR)\sprdop.obj" \
	"$(INTDIR)\stdafx.obj" \
	"$(INTDIR)\sysclass.obj" \
	"$(INTDIR)\table.obj" \
	"$(INTDIR)\tacemain.obj" \
	"$(INTDIR)\tags.obj" \
	"$(INTDIR)\timesubs.obj" \
	"$(INTDIR)\win32util.obj"

"$(OUTDIR)\wintace.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

# Name "acedb - Win32 Release"
# Name "acedb - Win32 Debug"
# Name "acedb - Win32 Acembly Debug"
# Name "acedb - Win32 Acembly Release"
# Name "acedb - Win32 AceClient Debug"
# Name "acedb - Win32 AceServer Debug"
# Name "acedb - Win32 AceServer Release"
# Name "acedb - Win32 AceClient Release"
# Name "acedb - Win32 NetClient Debug"
# Name "acedb - Win32 NetClient Release"
# Name "acedb - Win32 WinTace debug"
# Name "acedb - Win32 WinTace Release"

!IF  "$(CFG)" == "acedb - Win32 Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\acedb\w1\timesubs.c
DEP_CPP_TIMES=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\menu.c
DEP_CPP_MENU_=\
	"\acedb\wh\array.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\menu_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\menu.obj" : $(SOURCE) $(DEP_CPP_MENU_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\menu.obj" : $(SOURCE) $(DEP_CPP_MENU_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\bump.c
DEP_CPP_BUMP_=\
	"\acedb\wh\bump.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\dict.c
DEP_CPP_DICT_=\
	"\acedb\wh\array.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\filsubs.c
DEP_CPP_FILSU=\
	"\acedb\wh\array.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\mydirent.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\arraysub.c
DEP_CPP_ARRAY=\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\help.c
DEP_CPP_HELP_=\
	".\..\wgd\gd.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\help.sbr" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\help.sbr" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\messubs.c
DEP_CPP_MESSU=\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\messubs.obj" : $(SOURCE) $(DEP_CPP_MESSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\messubs.obj" : $(SOURCE) $(DEP_CPP_MESSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\call.c
DEP_CPP_CALL_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\heap.c
DEP_CPP_HEAP_=\
	"\acedb\wh\heap.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\heap.obj" : $(SOURCE) $(DEP_CPP_HEAP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\heap.obj" : $(SOURCE) $(DEP_CPP_HEAP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\freesubs.c
DEP_CPP_FREES=\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\randsubs.c

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\display.c
DEP_CPP_DISPL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
NODEP_CPP_DISPL=\
	".\..\w1\MPWIncludes.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\display.obj" : $(SOURCE) $(DEP_CPP_DISPL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\display.obj" : $(SOURCE) $(DEP_CPP_DISPL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\parse.c
DEP_CPP_PARSE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\parse.sbr" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\parse.sbr" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\newkey.c
DEP_CPP_NEWKE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\newkey.obj" : $(SOURCE) $(DEP_CPP_NEWKE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\newkey.obj" : $(SOURCE) $(DEP_CPP_NEWKE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\keyset.c
DEP_CPP_KEYSE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keyset.sbr" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keyset.sbr" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\chrono.c
DEP_CPP_CHRON=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\chrono.sbr" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\chrono.sbr" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\dotter.c
DEP_CPP_DOTTE=\
	"\acedb\wh\array.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\dotter.h"\
	"\acedb\wh\dotter_.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\dotter.obj" : $(SOURCE) $(DEP_CPP_DOTTE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dotter.obj" : $(SOURCE) $(DEP_CPP_DOTTE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\freeout.c
DEP_CPP_FREEO=\
	"\acedb\wh\array.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w2\graphcon.c
DEP_CPP_GRAPH=\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\graphcon.obj" : $(SOURCE) $(DEP_CPP_GRAPH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\graphcon.obj" : $(SOURCE) $(DEP_CPP_GRAPH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w2\graphsub.c
DEP_CPP_GRAPHS=\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\graphsub.obj" : $(SOURCE) $(DEP_CPP_GRAPHS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\graphsub.obj" : $(SOURCE) $(DEP_CPP_GRAPHS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\mainpick.c
DEP_CPP_MAINP=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\mainpick.obj" : $(SOURCE) $(DEP_CPP_MAINP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\mainpick.obj" : $(SOURCE) $(DEP_CPP_MAINP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\update.c
DEP_CPP_UPDAT=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\update.obj" : $(SOURCE) $(DEP_CPP_UPDAT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\update.obj" : $(SOURCE) $(DEP_CPP_UPDAT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\picksubs.c
DEP_CPP_PICKS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\picksubs.sbr" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\picksubs.sbr" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\longtext.c
DEP_CPP_LONGT=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\longtext.sbr" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\longtext.sbr" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\dump.c
DEP_CPP_DUMP_=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
NODEP_CPP_DUMP_=\
	".\..\w4\MPWIncludes.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dump.sbr" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dump.sbr" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\status.c
DEP_CPP_STATU=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\status.obj" : $(SOURCE) $(DEP_CPP_STATU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\status.obj" : $(SOURCE) $(DEP_CPP_STATU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\linkdate.c

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\linkdate.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\linkdate.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\model.c
DEP_CPP_MODEL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
NODEP_CPP_MODEL=\
	".\..\w4\newmodel.c"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\model.sbr" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\model.sbr" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\queryexe.c
DEP_CPP_QUERY=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\queryexe.sbr" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\queryexe.sbr" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapcontrol.c
DEP_CPP_FMAPC=\
	"..\wh\bitset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\matchtable.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\fmapcontrol.obj" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapcontrol.obj" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapcontrol.sbr" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

"$(INTDIR)\fmapcontrol.obj" : $(SOURCE) $(DEP_CPP_FMAPC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapgene.c
DEP_CPP_FMAPG=\
	"..\wh\bitset.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\fmapgene.obj" : $(SOURCE) $(DEP_CPP_FMAPG) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\fmapgene.obj" : $(SOURCE) $(DEP_CPP_FMAPG) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapfeatures.c
DEP_CPP_FMAPF=\
	"..\wh\bitset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dotter.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\fmapfeatures.obj" : $(SOURCE) $(DEP_CPP_FMAPF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\fmapfeatures.obj" : $(SOURCE) $(DEP_CPP_FMAPF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapsequence.c
DEP_CPP_FMAPS=\
	"..\wh\bitset.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\fmapsequence.obj" : $(SOURCE) $(DEP_CPP_FMAPS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\fmapsequence.obj" : $(SOURCE) $(DEP_CPP_FMAPS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmaplocuscol.c
DEP_CPP_GMAPL=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmaplocuscol.obj" : $(SOURCE) $(DEP_CPP_GMAPL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmaplocuscol.obj" : $(SOURCE) $(DEP_CPP_GMAPL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\mapcontrol.c
DEP_CPP_MAPCO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\mapcontrol.obj" : $(SOURCE) $(DEP_CPP_MAPCO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\mapcontrol.obj" : $(SOURCE) $(DEP_CPP_MAPCO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\metab.c
DEP_CPP_METAB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\metab.obj" : $(SOURCE) $(DEP_CPP_METAB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\metab.obj" : $(SOURCE) $(DEP_CPP_METAB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\drawdisp.c
DEP_CPP_DRAWD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\drawdisp.obj" : $(SOURCE) $(DEP_CPP_DRAWD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\drawdisp.obj" : $(SOURCE) $(DEP_CPP_DRAWD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapdisp.c
DEP_CPP_GMAPD=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapdisp.obj" : $(SOURCE) $(DEP_CPP_GMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapdisp.obj" : $(SOURCE) $(DEP_CPP_GMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\cmapdisp.c
DEP_CPP_CMAPD=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\cmapdisp.obj" : $(SOURCE) $(DEP_CPP_CMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\cmapdisp.obj" : $(SOURCE) $(DEP_CPP_CMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\vmapdrag.c
DEP_CPP_VMAPD=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\vmap.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\vmapdrag.obj" : $(SOURCE) $(DEP_CPP_VMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\vmapdrag.obj" : $(SOURCE) $(DEP_CPP_VMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapmarkercol.c
DEP_CPP_GMAPM=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapmarkercol.obj" : $(SOURCE) $(DEP_CPP_GMAPM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapmarkercol.obj" : $(SOURCE) $(DEP_CPP_GMAPM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\method.c
DEP_CPP_METHO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\method.obj" : $(SOURCE) $(DEP_CPP_METHO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\method.obj" : $(SOURCE) $(DEP_CPP_METHO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapremarkcol.c
DEP_CPP_GMAPR=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapremarkcol.obj" : $(SOURCE) $(DEP_CPP_GMAPR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapremarkcol.obj" : $(SOURCE) $(DEP_CPP_GMAPR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\griddisp.c
DEP_CPP_GRIDD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\griddisp.obj" : $(SOURCE) $(DEP_CPP_GRIDD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\griddisp.obj" : $(SOURCE) $(DEP_CPP_GRIDD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapintervalcol.c
DEP_CPP_GMAPI=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapintervalcol.obj" : $(SOURCE) $(DEP_CPP_GMAPI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapintervalcol.obj" : $(SOURCE) $(DEP_CPP_GMAPI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\pmapconvert.c
DEP_CPP_PMAPC=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\pmap_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\pmapconvert.obj" : $(SOURCE) $(DEP_CPP_PMAPC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\pmapconvert.obj" : $(SOURCE) $(DEP_CPP_PMAPC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\dnacpt.c
DEP_CPP_DNACP=\
	"..\wh\bitset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\dnacpt.obj" : $(SOURCE) $(DEP_CPP_DNACP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 0
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dnacpt.obj" : $(SOURCE) $(DEP_CPP_DNACP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapconvert.c
DEP_CPP_GMAPC=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapconvert.obj" : $(SOURCE) $(DEP_CPP_GMAPC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapconvert.obj" : $(SOURCE) $(DEP_CPP_GMAPC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fpdisp.c
DEP_CPP_FPDIS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fingerp.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\fpdisp.obj" : $(SOURCE) $(DEP_CPP_FPDIS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\fpdisp.obj" : $(SOURCE) $(DEP_CPP_FPDIS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\vmapdisp.c
DEP_CPP_VMAPDI=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\vmap.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\vmapdisp.obj" : $(SOURCE) $(DEP_CPP_VMAPDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\vmapdisp.obj" : $(SOURCE) $(DEP_CPP_VMAPDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\pmapdisp.c
DEP_CPP_PMAPD=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\grid.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\pmap_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\pmapdisp.obj" : $(SOURCE) $(DEP_CPP_PMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\pmapdisp.obj" : $(SOURCE) $(DEP_CPP_PMAPD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\geldisp.c
DEP_CPP_GELDI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\restriction.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\geldisp.obj" : $(SOURCE) $(DEP_CPP_GELDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\geldisp.obj" : $(SOURCE) $(DEP_CPP_GELDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\biblio.c
DEP_CPP_BIBLI=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\biblio.sbr" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\biblio.sbr" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\gmapposnegcol.c
DEP_CPP_GMAPP=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapposnegcol.obj" : $(SOURCE) $(DEP_CPP_GMAPP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapposnegcol.obj" : $(SOURCE) $(DEP_CPP_GMAPP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\tags.c
DEP_CPP_TAGS_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tags.sbr" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tags.sbr" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tags.sbr" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\sysclass.c
DEP_CPP_SYSCL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sysclass.sbr" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sysclass.sbr" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sysclass.sbr" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\quovadis.c
DEP_CPP_QUOVA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\syn.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\quovadis.sbr" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\quovadis.sbr" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\quovadis.sbr" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\class.c
DEP_CPP_CLASS=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\class.sbr" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\class.sbr" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\class.sbr" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\disptype.h

!IF  "$(CFG)" == "acedb - Win32 Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\classes.h

!IF  "$(CFG)" == "acedb - Win32 Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\tags.h

!IF  "$(CFG)" == "acedb - Win32 Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\sysclass.h

!IF  "$(CFG)" == "acedb - Win32 Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\whooks\systags.h

!IF  "$(CFG)" == "acedb - Win32 Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprdop.c
DEP_CPP_SPRDO=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdop.sbr" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdop.sbr" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprdctrl.c
DEP_CPP_SPRDC=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\sprdctrl.obj" : $(SOURCE) $(DEP_CPP_SPRDC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprdctrl.obj" : $(SOURCE) $(DEP_CPP_SPRDC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\ksetdisp.c
DEP_CPP_KSETD=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\ksetdisp.obj" : $(SOURCE) $(DEP_CPP_KSETD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\ksetdisp.obj" : $(SOURCE) $(DEP_CPP_KSETD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprdmap.c
DEP_CPP_SPRDM=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\sprdmap.obj" : $(SOURCE) $(DEP_CPP_SPRDM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprdmap.obj" : $(SOURCE) $(DEP_CPP_SPRDM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\qbedisp.c
DEP_CPP_QBEDI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\query_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\qbedisp.obj" : $(SOURCE) $(DEP_CPP_QBEDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\qbedisp.obj" : $(SOURCE) $(DEP_CPP_QBEDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprddisplay.c
DEP_CPP_SPRDD=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\sprddisplay.obj" : $(SOURCE) $(DEP_CPP_SPRDD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprddisplay.obj" : $(SOURCE) $(DEP_CPP_SPRDD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\querybuild.c
DEP_CPP_QUERYB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\query_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\querybuild.obj" : $(SOURCE) $(DEP_CPP_QUERYB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\querybuild.obj" : $(SOURCE) $(DEP_CPP_QUERYB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprddata.c
DEP_CPP_SPRDDA=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\sprddata.obj" : $(SOURCE) $(DEP_CPP_SPRDDA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\sprddata.obj" : $(SOURCE) $(DEP_CPP_SPRDDA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\treedisp.c
DEP_CPP_TREED=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\treedisp.obj" : $(SOURCE) $(DEP_CPP_TREED) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\treedisp.obj" : $(SOURCE) $(DEP_CPP_TREED) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bsubs.c
DEP_CPP_BSUBS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\b_.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsubs.sbr" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsubs.sbr" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bsdumps.c
DEP_CPP_BSDUM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsdumps.sbr" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsdumps.sbr" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\check.c
DEP_CPP_CHECK=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\check.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\check.sbr" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\check.sbr" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\peptide.c
DEP_CPP_PEPTI=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\peptide.sbr" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\peptide.sbr" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bstools.c
DEP_CPP_BSTOO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstools.sbr" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstools.sbr" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bstree.c
DEP_CPP_BSTRE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\b_.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstree.sbr" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstree.sbr" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\nicedump.c
DEP_CPP_NICED=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nicedump.sbr" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nicedump.sbr" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\dnasubs.c
DEP_CPP_DNASU=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnasubs.sbr" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnasubs.sbr" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\asubs.c
DEP_CPP_ASUBS=\
	"\acedb\wh\a.h"\
	"\acedb\wh\a_.h"\
	"\acedb\wh\acache.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\asubs.sbr" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\asubs.sbr" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\bssubs.c
DEP_CPP_BSSUB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\check.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bssubs.sbr" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bssubs.sbr" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\action.c
DEP_CPP_ACTIO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\action.obj" : $(SOURCE) $(DEP_CPP_ACTIO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\action.obj" : $(SOURCE) $(DEP_CPP_ACTIO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\querydisp.c
DEP_CPP_QUERYD=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\igdevent.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\query_.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\querydisp.obj" : $(SOURCE) $(DEP_CPP_QUERYD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\querydisp.obj" : $(SOURCE) $(DEP_CPP_QUERYD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\blocksub.c
DEP_CPP_BLOCK=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blocksub.sbr" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blocksub.sbr" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\objcache.c
DEP_CPP_OBJCA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\cache.h"\
	"\acedb\wh\cache_.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\objcache.sbr" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\objcache.sbr" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\disknew.c
DEP_CPP_DISKN=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\block.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
NODEP_CPP_DISKN=\
	".\..\w5\macfilemgrsubs.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\disknew.sbr" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\disknew.sbr" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\lexsubs4.c
DEP_CPP_LEXSU=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs4.sbr" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs4.sbr" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\keysetdump.c
DEP_CPP_KEYSET=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keysetdump.sbr" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keysetdump.sbr" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\lexalpha.c
DEP_CPP_LEXAL=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexalpha.sbr" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexalpha.sbr" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w5\lexsubs.c
DEP_CPP_LEXSUB=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\lex_bl_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs.sbr" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs.sbr" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gmapdata.c
DEP_CPP_GMAPDA=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapdata.obj" : $(SOURCE) $(DEP_CPP_GMAPDA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapdata.obj" : $(SOURCE) $(DEP_CPP_GMAPDA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\embl.c
DEP_CPP_EMBL_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\embl.obj" : $(SOURCE) $(DEP_CPP_EMBL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\embl.obj" : $(SOURCE) $(DEP_CPP_EMBL_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gfcode.c
DEP_CPP_GFCOD=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gfcode.obj" : $(SOURCE) $(DEP_CPP_GFCOD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gfcode.obj" : $(SOURCE) $(DEP_CPP_GFCOD) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\vmapphys.c
DEP_CPP_VMAPP=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\vmap.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\vmapphys.obj" : $(SOURCE) $(DEP_CPP_VMAPP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\vmapphys.obj" : $(SOURCE) $(DEP_CPP_VMAPP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\vmapdata2.c
DEP_CPP_VMAPDA=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\vmap.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\vmapdata2.obj" : $(SOURCE) $(DEP_CPP_VMAPDA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\vmapdata2.obj" : $(SOURCE) $(DEP_CPP_VMAPDA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\align.c
DEP_CPP_ALIGN=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\align.obj" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\align.obj" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

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
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\translate.obj" : $(SOURCE) $(DEP_CPP_TRANS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\translate.obj" : $(SOURCE) $(DEP_CPP_TRANS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gmapphys.c
DEP_CPP_GMAPPH=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapphys.obj" : $(SOURCE) $(DEP_CPP_GMAPPH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapphys.obj" : $(SOURCE) $(DEP_CPP_GMAPPH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\gmapdatacol.c
DEP_CPP_GMAPDAT=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gmapdatacol.obj" : $(SOURCE) $(DEP_CPP_GMAPDAT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gmapdatacol.obj" : $(SOURCE) $(DEP_CPP_GMAPDAT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w9\blxview.c
DEP_CPP_BLXVI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dotter.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\blxview.obj" : $(SOURCE) $(DEP_CPP_BLXVI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\blxview.obj" : $(SOURCE) $(DEP_CPP_BLXVI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winace.cpp

!IF  "$(CFG)" == "acedb - Win32 Release"

DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /FAcs

"$(INTDIR)\winace.obj" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /FAcs /Fa"$(INTDIR)/" /Fp"$(INTDIR)/Acedb.pch" /YX\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\winace.obj" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\winace.sbr" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\winace.obj" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winace.sbr" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /FAcs
# ADD CPP /FAcs

"$(INTDIR)\winace.obj" : $(SOURCE) $(DEP_CPP_WINAC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /FAcs /Fa"$(INTDIR)/" /Fp"$(INTDIR)/Acedb.pch" /YX\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /FAcs
# ADD CPP /FAcs

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /FAcs
# ADD CPP /FAcs

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	
# ADD BASE CPP /FAcs
# ADD CPP /FAcs

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WINAC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /FAcs
# ADD CPP /FAcs /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winace.rc

!IF  "$(CFG)" == "acedb - Win32 Release"

DEP_RSC_WINACE=\
	"..\win32\res\acelogo.bmp"\
	"..\win32\res\ico00001.ico"\
	"..\win32\res\ico00002.ico"\
	"..\win32\res\ico00003.ico"\
	"..\win32\res\idr_map_.ico"\
	"..\win32\res\idr_pixe.ico"\
	"..\win32\res\idr_text.ico"\
	"..\win32\res\idr_wina.ico"\
	"..\win32\res\KeySetdoc.ico"\
	"..\win32\res\toolbar.bmp"\
	"..\win32\res\winace.ico"\
	"..\win32\res\winace.rc2"\
	

"$(INTDIR)\winace.res" : $(SOURCE) $(DEP_RSC_WINACE) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

DEP_RSC_WINACE=\
	"..\win32\res\acelogo.bmp"\
	"..\win32\res\ico00001.ico"\
	"..\win32\res\ico00002.ico"\
	"..\win32\res\ico00003.ico"\
	"..\win32\res\idr_map_.ico"\
	"..\win32\res\idr_pixe.ico"\
	"..\win32\res\idr_text.ico"\
	"..\win32\res\idr_wina.ico"\
	"..\win32\res\KeySetdoc.ico"\
	"..\win32\res\toolbar.bmp"\
	"..\win32\res\winace.ico"\
	"..\win32\res\winace.rc2"\
	

"$(INTDIR)\winace.res" : $(SOURCE) $(DEP_RSC_WINACE) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

DEP_RSC_WINACE=\
	"..\win32\res\acelogo.bmp"\
	"..\win32\res\ico00001.ico"\
	"..\win32\res\ico00002.ico"\
	"..\win32\res\ico00003.ico"\
	"..\win32\res\idr_map_.ico"\
	"..\win32\res\idr_pixe.ico"\
	"..\win32\res\idr_text.ico"\
	"..\win32\res\idr_wina.ico"\
	"..\win32\res\KeySetdoc.ico"\
	"..\win32\res\toolbar.bmp"\
	"..\win32\res\winace.ico"\
	"..\win32\res\winace.rc2"\
	

"$(INTDIR)\winace.res" : $(SOURCE) $(DEP_RSC_WINACE) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "_DEBUG" /d "_AFXDLL"\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

DEP_RSC_WINACE=\
	"..\win32\res\acelogo.bmp"\
	"..\win32\res\ico00001.ico"\
	"..\win32\res\ico00002.ico"\
	"..\win32\res\ico00003.ico"\
	"..\win32\res\idr_map_.ico"\
	"..\win32\res\idr_pixe.ico"\
	"..\win32\res\idr_text.ico"\
	"..\win32\res\idr_wina.ico"\
	"..\win32\res\KeySetdoc.ico"\
	"..\win32\res\toolbar.bmp"\
	"..\win32\res\winace.ico"\
	"..\win32\res\winace.rc2"\
	

"$(INTDIR)\winace.res" : $(SOURCE) $(DEP_RSC_WINACE) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/winace.res" /i "\acedb\win32h" /i\
 "\acedb\win32" /i "\acedb\wh" /i "\acedb\whooks" /d "NDEBUG" /d "_AFXDLL"\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\stdafx.cpp
DEP_CPP_STDAF=\
	"\acedb\win32h\stdafx.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# ADD BASE CPP /Yc"stdafx.h"
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0
# ADD BASE CPP /Yc"stdafx.h"
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# ADD BASE CPP /Yc"stdafx.h"
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 0

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 0
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# ADD BASE CPP /Yc"stdafx.h"
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# ADD BASE CPP /Yc"stdafx.h"
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\stdafx.sbr" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0
# ADD BASE CPP /Yc"stdafx.h"
# ADD CPP /Yc"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yc"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\stdafx.obj" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Acedb.pch" : $(SOURCE) $(DEP_CPP_STDAF) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\mainfrm.cpp
DEP_CPP_MAINF=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\mainfrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\mainfrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\mainfrm.sbr" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\mainfrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mainfrm.sbr" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\mainfrm.obj" : $(SOURCE) $(DEP_CPP_MAINF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphwin32_.cpp
DEP_CPP_GRAPHW=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\graphwin32_.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\graphwin32_.obj" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphwin32_.obj" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\graphwin32_.sbr" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphwin32_.obj" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphwin32_.sbr" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\graphwin32_.obj" : $(SOURCE) $(DEP_CPP_GRAPHW) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\afxtrace.cpp
DEP_CPP_AFXTR=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\afxtrace.obj" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\afxtrace.obj" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\afxtrace.sbr" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\afxtrace.obj" : $(SOURCE) $(DEP_CPP_AFXTR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32print.cpp
DEP_CPP_WIN32=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\win32print.obj" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32print.obj" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32print.sbr" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\win32print.obj" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32print.sbr" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\win32print.obj" : $(SOURCE) $(DEP_CPP_WIN32) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\msgwindow.cpp
DEP_CPP_MSGWI=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\msgwindow.obj" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\msgwindow.obj" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\msgwindow.sbr" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\msgwindow.obj" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\msgwindow.sbr" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\msgwindow.obj" : $(SOURCE) $(DEP_CPP_MSGWI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32.c
DEP_CPP_WIN32_=\
	"\acedb\win32h\win32.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\win32.obj" : $(SOURCE) $(DEP_CPP_WIN32_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\win32.obj" : $(SOURCE) $(DEP_CPP_WIN32_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\cmultiass.cpp
DEP_CPP_CMULT=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\cmultiass.obj" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\cmultiass.obj" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\cmultiass.sbr" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\cmultiass.obj" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\cmultiass.sbr" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\cmultiass.obj" : $(SOURCE) $(DEP_CPP_CMULT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32process.cpp

!IF  "$(CFG)" == "acedb - Win32 Release"

DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32process.obj" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32process.obj" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32process.sbr" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\win32process.obj" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32process.sbr" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

"$(INTDIR)\win32process.obj" : $(SOURCE) $(DEP_CPP_WIN32P) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_WIN32P=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32lib.cpp
DEP_CPP_WIN32L=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\win32lib.obj" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32lib.obj" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32lib.sbr" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\win32lib.obj" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32lib.sbr" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\win32lib.obj" : $(SOURCE) $(DEP_CPP_WIN32L) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wgd\gd.c
DEP_CPP_GD_Cd6=\
	".\..\wgd\gd.h"\
	".\..\wgd\mtables.c"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\gd.obj" : $(SOURCE) $(DEP_CPP_GD_Cd6) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gd.obj" : $(SOURCE) $(DEP_CPP_GD_Cd6) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gd.sbr" : $(SOURCE) $(DEP_CPP_GD_Cd6) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\gd.obj" : $(SOURCE) $(DEP_CPP_GD_Cd6) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\gd.sbr" : $(SOURCE) $(DEP_CPP_GD_Cd6) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\gd.obj" : $(SOURCE) $(DEP_CPP_GD_Cd6) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wstaden\opp.c
DEP_CPP_OPP_C=\
	"\acedb\wh\mach-io.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\opp.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\seq.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\opp.obj" : $(SOURCE) $(DEP_CPP_OPP_C) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\mapscrollview.cpp
DEP_CPP_MAPSC=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\mapscrollview.obj" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\mapscrollview.obj" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\mapscrollview.sbr" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\mapscrollview.obj" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mapscrollview.sbr" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\mapscrollview.obj" : $(SOURCE) $(DEP_CPP_MAPSC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\textfitview.cpp
DEP_CPP_TEXTF=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\textfitview.obj" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\textfitview.obj" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\textfitview.sbr" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\textfitview.obj" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\textfitview.sbr" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\textfitview.obj" : $(SOURCE) $(DEP_CPP_TEXTF) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\textscrollview.cpp
DEP_CPP_TEXTS=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\textscrollview.obj" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\textscrollview.obj" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\textscrollview.sbr" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\textscrollview.obj" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\textscrollview.sbr" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\textscrollview.obj" : $(SOURCE) $(DEP_CPP_TEXTS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\pixelfitview.cpp
DEP_CPP_PIXEL=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\pixelfitview.obj" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pixelfitview.obj" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\pixelfitview.sbr" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pixelfitview.obj" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pixelfitview.sbr" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\pixelfitview.obj" : $(SOURCE) $(DEP_CPP_PIXEL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\pixelscrollview.cpp
DEP_CPP_PIXELS=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\pixelscrollview.obj" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\pixelscrollview.obj" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\pixelscrollview.sbr" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pixelscrollview.obj" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pixelscrollview.sbr" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\pixelscrollview.obj" : $(SOURCE) $(DEP_CPP_PIXELS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\acefiledoc.cpp
DEP_CPP_ACEFI=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\acefiledoc.obj" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\acefiledoc.obj" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\acefiledoc.sbr" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\acefiledoc.obj" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acefiledoc.sbr" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\acefiledoc.obj" : $(SOURCE) $(DEP_CPP_ACEFI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\acefileview.cpp
DEP_CPP_ACEFIL=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\acefiledoc.h"\
	"\acedb\win32h\acefileview.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\acefileview.obj" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\acefileview.obj" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\acefileview.sbr" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\acefileview.obj" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acefileview.sbr" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\acefileview.obj" : $(SOURCE) $(DEP_CPP_ACEFIL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\fullscrollview.cpp
DEP_CPP_FULLS=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\fullscrollview.obj" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fullscrollview.obj" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\fullscrollview.sbr" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fullscrollview.obj" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fullscrollview.sbr" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\fullscrollview.obj" : $(SOURCE) $(DEP_CPP_FULLS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\cgraph.cpp
DEP_CPP_CGRAP=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\FullScrollView.h"\
	"\acedb\win32h\graphwin32_.h"\
	"\acedb\win32h\MapScrollView.h"\
	"\acedb\win32h\PixelFitView.h"\
	"\acedb\win32h\PixelScrollView.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\TextFitView.h"\
	"\acedb\win32h\TextScrollView.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\cgraph.obj" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\cgraph.obj" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\cgraph.sbr" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\cgraph.obj" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\cgraph.sbr" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\cgraph.obj" : $(SOURCE) $(DEP_CPP_CGRAP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphview.cpp

!IF  "$(CFG)" == "acedb - Win32 Release"

DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

"$(INTDIR)\graphview.obj" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphview.obj" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\graphview.sbr" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphview.obj" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphview.sbr" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

"$(INTDIR)\graphview.obj" : $(SOURCE) $(DEP_CPP_GRAPHV) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_GRAPHV=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32xcolor.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32menusubs.cpp
DEP_CPP_WIN32M=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\menu_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\win32menusubs.obj" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32menusubs.obj" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32menusubs.sbr" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\win32menusubs.obj" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32menusubs.sbr" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\win32menusubs.obj" : $(SOURCE) $(DEP_CPP_WIN32M) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphselect.cpp
DEP_CPP_GRAPHSE=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\graphselect.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\graphselect.obj" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphselect.obj" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\graphselect.sbr" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphselect.obj" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphselect.sbr" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\graphselect.obj" : $(SOURCE) $(DEP_CPP_GRAPHSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\windialogs.cpp
DEP_CPP_WINDI=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	"\acedb\win32h\windialogs.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\windialogs.obj" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\windialogs.obj" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\windialogs.sbr" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\windialogs.obj" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\windialogs.sbr" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\windialogs.obj" : $(SOURCE) $(DEP_CPP_WINDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winfilquery.cpp
DEP_CPP_WINFI=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\winfilquery.obj" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\winfilquery.obj" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\winfilquery.sbr" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\winfilquery.obj" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winfilquery.sbr" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\winfilquery.obj" : $(SOURCE) $(DEP_CPP_WINFI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W6\alignment.c
DEP_CPP_ALIGNM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\alignment.sbr" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\alignment.sbr" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGNM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W9\dotterKarlin.c
DEP_CPP_DOTTER=\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\dotterKarlin.obj" : $(SOURCE) $(DEP_CPP_DOTTER) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\dotterKarlin.obj" : $(SOURCE) $(DEP_CPP_DOTTER) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W9\hexcode.c
DEP_CPP_HEXCO=\
	"..\wh\bitset.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\hexcode.obj" : $(SOURCE) $(DEP_CPP_HEXCO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\hexcode.obj" : $(SOURCE) $(DEP_CPP_HEXCO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\caceprintpage.cpp
DEP_CPP_CACEP=\
	"\acedb\win32h\caceprintpage.h"\
	"\acedb\win32h\stdafx.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\caceprintpage.obj" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\caceprintpage.obj" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\caceprintpage.sbr" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\caceprintpage.obj" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\caceprintpage.sbr" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\caceprintpage.obj" : $(SOURCE) $(DEP_CPP_CACEP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\preferences.cpp
DEP_CPP_PREFE=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\preferences.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\preferences.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\preferences.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\preferences.sbr" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\preferences.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\preferences.sbr" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\preferences.obj" : $(SOURCE) $(DEP_CPP_PREFE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\dbprofile.cpp
DEP_CPP_DBPRO=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\dbprofile.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\dbprofile.obj" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\dbprofile.obj" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\dbprofile.sbr" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dbprofile.obj" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dbprofile.sbr" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\dbprofile.obj" : $(SOURCE) $(DEP_CPP_DBPRO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\userinterface.cpp
DEP_CPP_USERI=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\userinterface.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\userinterface.obj" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\userinterface.obj" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\userinterface.sbr" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\userinterface.obj" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\userinterface.sbr" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\userinterface.obj" : $(SOURCE) $(DEP_CPP_USERI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\acedbprofileview.cpp
DEP_CPP_ACEDB=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\acedbprofileview.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\acedbprofileview.obj" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\acedbprofileview.obj" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\acedbprofileview.sbr" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\acedbprofileview.obj" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acedbprofileview.sbr" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\acedbprofileview.obj" : $(SOURCE) $(DEP_CPP_ACEDB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\acedbprofile.cpp
DEP_CPP_ACEDBP=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\acedbprofile.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\acedbprofile.obj" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\acedbprofile.obj" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\acedbprofile.sbr" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\acedbprofile.obj" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acedbprofile.sbr" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\acedbprofile.obj" : $(SOURCE) $(DEP_CPP_ACEDBP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\Win32\fontpreference.cpp
DEP_CPP_FONTP=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\fontpreference.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\fontpreference.obj" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fontpreference.obj" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\fontpreference.sbr" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fontpreference.obj" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fontpreference.sbr" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\fontpreference.obj" : $(SOURCE) $(DEP_CPP_FONTP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\WIN32\splashbox.cpp
DEP_CPP_SPLAS=\
	"..\win32h\graphloop.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\splashbox.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\splashbox.obj" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\splashbox.obj" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\splashbox.sbr" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\splashbox.obj" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\splashbox.sbr" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\splashbox.obj" : $(SOURCE) $(DEP_CPP_SPLAS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pepseqcol.c
DEP_CPP_PEPSE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\pepseqcol.obj" : $(SOURCE) $(DEP_CPP_PEPSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\pepseqcol.obj" : $(SOURCE) $(DEP_CPP_PEPSE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pephomolcol.c
DEP_CPP_PEPHO=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dotter.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\matchtable.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\pephomolcol.obj" : $(SOURCE) $(DEP_CPP_PEPHO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\pephomolcol.obj" : $(SOURCE) $(DEP_CPP_PEPHO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pepgraphcol.c
DEP_CPP_PEPGR=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\pepgraphcol.obj" : $(SOURCE) $(DEP_CPP_PEPGR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\pepgraphcol.obj" : $(SOURCE) $(DEP_CPP_PEPGR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\pepdisp.c
DEP_CPP_PEPDI=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\pepdisp.obj" : $(SOURCE) $(DEP_CPP_PEPDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\pepdisp.obj" : $(SOURCE) $(DEP_CPP_PEPDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W2\colcontrol.c
DEP_CPP_COLCO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\colcontrol_.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\colcontrol.obj" : $(SOURCE) $(DEP_CPP_COLCO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\colcontrol.obj" : $(SOURCE) $(DEP_CPP_COLCO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W2\viewedit.c
DEP_CPP_VIEWE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\colcontrol_.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\viewedit.obj" : $(SOURCE) $(DEP_CPP_VIEWE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\viewedit.obj" : $(SOURCE) $(DEP_CPP_VIEWE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\acedialogs.c
DEP_CPP_ACEDI=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\acedialogs.obj" : $(SOURCE) $(DEP_CPP_ACEDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\acedialogs.obj" : $(SOURCE) $(DEP_CPP_ACEDI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winmain.c
DEP_CPP_WINMA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\win32h\win32.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\winmain.obj" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\winmain.obj" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winmain.sbr" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# ADD CPP /D "ACEMBLY"
# SUBTRACT CPP /YX

"$(INTDIR)\winmain.obj" : $(SOURCE) $(DEP_CPP_WINMA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W4\banner.c
DEP_CPP_BANNE=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\banner.sbr" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\banner.sbr" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\graphbox.cpp
DEP_CPP_GRAPHB=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\graphbox.obj" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphbox.obj" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\graphbox.sbr" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphbox.obj" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphbox.sbr" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\graphbox.obj" : $(SOURCE) $(DEP_CPP_GRAPHB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\plainview.cpp
DEP_CPP_PLAIN=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\PlainView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\plainview.obj" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\plainview.obj" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\plainview.sbr" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\plainview.obj" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\plainview.sbr" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\plainview.obj" : $(SOURCE) $(DEP_CPP_PLAIN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W1\oldhelp.c
DEP_CPP_OLDHE=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oldhelp.sbr" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oldhelp.sbr" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W4\session.c
DEP_CPP_SESSI=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\byteswap.h"\
	"\acedb\wh\client.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\session.sbr" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\session.sbr" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W7\colourbox.c
DEP_CPP_COLOU=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\colourbox.obj" : $(SOURCE) $(DEP_CPP_COLOU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\colourbox.obj" : $(SOURCE) $(DEP_CPP_COLOU) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W2\graphramp.c
DEP_CPP_GRAPHR=\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\graphramp.obj" : $(SOURCE) $(DEP_CPP_GRAPHR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\graphramp.obj" : $(SOURCE) $(DEP_CPP_GRAPHR) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\W6\plot.c
DEP_CPP_PLOT_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# SUBTRACT CPP /YX

"$(INTDIR)\plot.obj" : $(SOURCE) $(DEP_CPP_PLOT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# SUBTRACT BASE CPP /YX
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

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

"$(INTDIR)\plot.obj" : $(SOURCE) $(DEP_CPP_PLOT_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\WIN32\fitview.cpp
DEP_CPP_FITVI=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\FitView.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\fitview.obj" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\fitview.obj" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\fitview.sbr" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fitview.obj" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fitview.sbr" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\fitview.obj" : $(SOURCE) $(DEP_CPP_FITVI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\Acedb\win32\graphwin.cpp
DEP_CPP_GRAPHWI=\
	"..\win32h\cmultiass.h"\
	"..\win32h\graphbox.h"\
	"..\win32h\graphloop.h"\
	"..\win32h\graphview.h"\
	"..\win32h\graphwin.h"\
	"..\win32h\msgwindow.h"\
	"..\win32h\win32menus.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\graph_.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\cgraph.h"\
	"\acedb\win32h\mainfrm.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	"\acedb\win32h\winace.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\graphwin.obj" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\graphwin.obj" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\graphwin.sbr" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\graphwin.obj" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\graphwin.sbr" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\graphwin.obj" : $(SOURCE) $(DEP_CPP_GRAPHWI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\forest.c
DEP_CPP_FORES=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\forest.obj" : $(SOURCE) $(DEP_CPP_FORES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\forest.obj" : $(SOURCE) $(DEP_CPP_FORES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\forest.sbr" : $(SOURCE) $(DEP_CPP_FORES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\forest.obj" : $(SOURCE) $(DEP_CPP_FORES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\forest.sbr" : $(SOURCE) $(DEP_CPP_FORES) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\forest.obj" : $(SOURCE) $(DEP_CPP_FORES) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\forestdisplay.c

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\forestdisplay.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\forestdisplay.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\forestdisplay.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\forestdisplay.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\forestdisplay.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\forestdisplay.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wabi\basecall.c
DEP_CPP_BASEC=\
	".\..\wh\seqIOSCF.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\basecall.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mach-io.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\seq.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\basecall.obj" : $(SOURCE) $(DEP_CPP_BASEC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\basecall.sbr" : $(SOURCE) $(DEP_CPP_BASEC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\basecall.obj" : $(SOURCE) $(DEP_CPP_BASEC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wabi\fmaptrace.c
DEP_CPP_FMAPT=\
	"..\wh\bitset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\blxview.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\dotter.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmaptrace.obj" : $(SOURCE) $(DEP_CPP_FMAPT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmaptrace.sbr" : $(SOURCE) $(DEP_CPP_FMAPT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\fmaptrace.obj" : $(SOURCE) $(DEP_CPP_FMAPT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wabi\acemblyhook.c
DEP_CPP_ACEMB=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\acemblyhook.obj" : $(SOURCE) $(DEP_CPP_ACEMB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\acemblyhook.sbr" : $(SOURCE) $(DEP_CPP_ACEMB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\acemblyhook.obj" : $(SOURCE) $(DEP_CPP_ACEMB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wabi\abifix.c
DEP_CPP_ABIFI=\
	".\..\wh\seqIOSCF.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\basecall.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mach-io.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\seq.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\abifix.obj" : $(SOURCE) $(DEP_CPP_ABIFI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\abifix.sbr" : $(SOURCE) $(DEP_CPP_ABIFI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\abifix.obj" : $(SOURCE) $(DEP_CPP_ABIFI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wabi\myNetwork.c

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\myNetwork.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\myNetwork.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\myNetwork.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wabi\trace.c
DEP_CPP_TRACE=\
	".\..\wh\seqIOSCF.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\basecall.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mach-io.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\menu.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\seq.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\trace.obj" : $(SOURCE) $(DEP_CPP_TRACE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\trace.sbr" : $(SOURCE) $(DEP_CPP_TRACE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\trace.obj" : $(SOURCE) $(DEP_CPP_TRACE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\aligntools.c
DEP_CPP_ALIGNT=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\aligntools.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\topology.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\aligntools.obj" : $(SOURCE) $(DEP_CPP_ALIGNT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\aligntools.sbr" : $(SOURCE) $(DEP_CPP_ALIGNT) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\aligntools.obj" : $(SOURCE) $(DEP_CPP_ALIGNT) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\basecallstat.c
DEP_CPP_BASECA=\
	"..\wh\bitset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\aligntools.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\topology.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\basecallstat.obj" : $(SOURCE) $(DEP_CPP_BASECA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\basecallstat.sbr" : $(SOURCE) $(DEP_CPP_BASECA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\basecallstat.obj" : $(SOURCE) $(DEP_CPP_BASECA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\defcpt.c
DEP_CPP_DEFCP=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\fingerp.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\topology.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\defcpt.obj" : $(SOURCE) $(DEP_CPP_DEFCP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\defcpt.sbr" : $(SOURCE) $(DEP_CPP_DEFCP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\defcpt.obj" : $(SOURCE) $(DEP_CPP_DEFCP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\dnaalign.c
DEP_CPP_DNAAL=\
	"..\wh\bitset.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\acembly.h"\
	"\acedb\wh\aligntools.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dnaalign.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\plot.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\topology.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dnaalign.obj" : $(SOURCE) $(DEP_CPP_DNAAL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnaalign.sbr" : $(SOURCE) $(DEP_CPP_DNAAL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\dnaalign.obj" : $(SOURCE) $(DEP_CPP_DNAAL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\intrinsictree.c
DEP_CPP_INTRI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\interval.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\topology.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\intrinsictree.obj" : $(SOURCE) $(DEP_CPP_INTRI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\intrinsictree.sbr" : $(SOURCE) $(DEP_CPP_INTRI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\intrinsictree.obj" : $(SOURCE) $(DEP_CPP_INTRI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\regression.c
DEP_CPP_REGRE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regression.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\regression.obj" : $(SOURCE) $(DEP_CPP_REGRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\regression.sbr" : $(SOURCE) $(DEP_CPP_REGRE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\regression.obj" : $(SOURCE) $(DEP_CPP_REGRE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w8\topology.c
DEP_CPP_TOPOL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\topology.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\topology.obj" : $(SOURCE) $(DEP_CPP_TOPOL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\topology.sbr" : $(SOURCE) $(DEP_CPP_TOPOL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\topology.obj" : $(SOURCE) $(DEP_CPP_TOPOL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wstaden\seqIOSCF.c
DEP_CPP_SEQIO=\
	"\acedb\wh\mach-io.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\scf.h"\
	"\acedb\wh\seq.h"\
	"\acedb\wh\seqIOEdit.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\seqIOSCF.obj" : $(SOURCE) $(DEP_CPP_SEQIO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\seqIOSCF.sbr" : $(SOURCE) $(DEP_CPP_SEQIO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\seqIOSCF.obj" : $(SOURCE) $(DEP_CPP_SEQIO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wstaden\seq.c
DEP_CPP_SEQ_C=\
	"\acedb\wh\mach-io.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\seq.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\seq.obj" : $(SOURCE) $(DEP_CPP_SEQ_C) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\seq.sbr" : $(SOURCE) $(DEP_CPP_SEQ_C) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\seq.obj" : $(SOURCE) $(DEP_CPP_SEQ_C) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\acedb\wstaden\mach-io.c"
DEP_CPP_MACH_=\
	"\acedb\wh\mach-io.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\mach-io.obj" : $(SOURCE) $(DEP_CPP_MACH_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\mach-io.sbr" : $(SOURCE) $(DEP_CPP_MACH_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\mach-io.obj" : $(SOURCE) $(DEP_CPP_MACH_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\command.c
DEP_CPP_COMMA=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# ADD CPP /D "ACEMBLY"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\command.sbr" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# ADD CPP /D "ACEMBLY"

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /D "ACEMBLY" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\command.sbr" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\command.sbr" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\basepad.c
DEP_CPP_BASEP=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\basepad.sbr" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\basepad.sbr" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\basepad.sbr" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\pepfeaturecol.c
DEP_CPP_PEPFE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\pepfeaturecol.obj" : $(SOURCE) $(DEP_CPP_PEPFE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepfeaturecol.obj" : $(SOURCE) $(DEP_CPP_PEPFE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepfeaturecol.sbr" : $(SOURCE) $(DEP_CPP_PEPFE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepfeaturecol.obj" : $(SOURCE) $(DEP_CPP_PEPFE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepfeaturecol.sbr" : $(SOURCE) $(DEP_CPP_PEPFE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\pepfeaturecol.obj" : $(SOURCE) $(DEP_CPP_PEPFE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\pepactivezonecol.c
DEP_CPP_PEPAC=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\colcontrol.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\gmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pepdisp.h"\
	"\acedb\wh\peptide.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\pepactivezonecol.obj" : $(SOURCE) $(DEP_CPP_PEPAC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepactivezonecol.obj" : $(SOURCE) $(DEP_CPP_PEPAC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepactivezonecol.sbr" : $(SOURCE) $(DEP_CPP_PEPAC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pepactivezonecol.obj" : $(SOURCE) $(DEP_CPP_PEPAC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pepactivezonecol.sbr" : $(SOURCE) $(DEP_CPP_PEPAC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\pepactivezonecol.obj" : $(SOURCE) $(DEP_CPP_PEPAC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# SUBTRACT BASE CPP /YX
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\winacemail.cpp
DEP_CPP_WINACEM=\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\winacemail.obj" : $(SOURCE) $(DEP_CPP_WINACEM) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\winacemail.obj" : $(SOURCE) $(DEP_CPP_WINACEM) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\winacemail.sbr" : $(SOURCE) $(DEP_CPP_WINACEM) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\winacemail.obj" : $(SOURCE) $(DEP_CPP_WINACEM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\winacemail.sbr" : $(SOURCE) $(DEP_CPP_WINACEM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\winacemail.obj" : $(SOURCE) $(DEP_CPP_WINACEM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD BASE CPP /Yu"stdafx.h"
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\matchtable.c
DEP_CPP_MATCH=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\matchtable.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX /Yc /Yu

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matchtable.sbr" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matchtable.sbr" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matchtable.sbr" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matchtable.sbr" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\matchtable.obj" : $(SOURCE) $(DEP_CPP_MATCH) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w6\sprddef.c
DEP_CPP_SPRDDE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddef.sbr" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddef.sbr" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddef.sbr" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprddef.sbr" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\sprddef.obj" : $(SOURCE) $(DEP_CPP_SPRDDE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapblast.c
DEP_CPP_FMAPB=\
	"..\wh\bitset.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\fmapblast.obj" : $(SOURCE) $(DEP_CPP_FMAPB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapblast.obj" : $(SOURCE) $(DEP_CPP_FMAPB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapblast.sbr" : $(SOURCE) $(DEP_CPP_FMAPB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapblast.obj" : $(SOURCE) $(DEP_CPP_FMAPB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapblast.sbr" : $(SOURCE) $(DEP_CPP_FMAPB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\fmapblast.obj" : $(SOURCE) $(DEP_CPP_FMAPB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmaposp.c
DEP_CPP_FMAPO=\
	"..\wh\bitset.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\fmaposp.obj" : $(SOURCE) $(DEP_CPP_FMAPO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmaposp.obj" : $(SOURCE) $(DEP_CPP_FMAPO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmaposp.sbr" : $(SOURCE) $(DEP_CPP_FMAPO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmaposp.obj" : $(SOURCE) $(DEP_CPP_FMAPO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmaposp.sbr" : $(SOURCE) $(DEP_CPP_FMAPO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\fmaposp.obj" : $(SOURCE) $(DEP_CPP_FMAPO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wjo\oxgriddisp.c
DEP_CPP_OXGRI=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\oxgrid.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\oxgriddisp.obj" : $(SOURCE) $(DEP_CPP_OXGRI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oxgriddisp.obj" : $(SOURCE) $(DEP_CPP_OXGRI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oxgriddisp.sbr" : $(SOURCE) $(DEP_CPP_OXGRI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oxgriddisp.obj" : $(SOURCE) $(DEP_CPP_OXGRI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oxgriddisp.sbr" : $(SOURCE) $(DEP_CPP_OXGRI) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\oxgriddisp.obj" : $(SOURCE) $(DEP_CPP_OXGRI) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wjo\specg.c
DEP_CPP_SPECG=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\oxgrid.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\specg.obj" : $(SOURCE) $(DEP_CPP_SPECG) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX /Yc /Yu

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\specg.obj" : $(SOURCE) $(DEP_CPP_SPECG) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\specg.sbr" : $(SOURCE) $(DEP_CPP_SPECG) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\specg.obj" : $(SOURCE) $(DEP_CPP_SPECG) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\specg.sbr" : $(SOURCE) $(DEP_CPP_SPECG) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\specg.obj" : $(SOURCE) $(DEP_CPP_SPECG) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wjo\pairmapdisp.c
DEP_CPP_PAIRM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\oxgrid.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\pairmapdisp.obj" : $(SOURCE) $(DEP_CPP_PAIRM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pairmapdisp.obj" : $(SOURCE) $(DEP_CPP_PAIRM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pairmapdisp.sbr" : $(SOURCE) $(DEP_CPP_PAIRM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\pairmapdisp.obj" : $(SOURCE) $(DEP_CPP_PAIRM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\pairmapdisp.sbr" : $(SOURCE) $(DEP_CPP_PAIRM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\pairmapdisp.obj" : $(SOURCE) $(DEP_CPP_PAIRM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wjo\oxhomlist.c
DEP_CPP_OXHOM=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\oxgrid.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\oxhomlist.obj" : $(SOURCE) $(DEP_CPP_OXHOM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oxhomlist.obj" : $(SOURCE) $(DEP_CPP_OXHOM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oxhomlist.sbr" : $(SOURCE) $(DEP_CPP_OXHOM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\oxhomlist.obj" : $(SOURCE) $(DEP_CPP_OXHOM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oxhomlist.sbr" : $(SOURCE) $(DEP_CPP_OXHOM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\oxhomlist.obj" : $(SOURCE) $(DEP_CPP_OXHOM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wjo\o2m.c
DEP_CPP_O2M_C=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\oxgrid.h"\
	"\acedb\wh\pmap.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\o2m.obj" : $(SOURCE) $(DEP_CPP_O2M_C) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\o2m.obj" : $(SOURCE) $(DEP_CPP_O2M_C) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\o2m.sbr" : $(SOURCE) $(DEP_CPP_O2M_C) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\o2m.obj" : $(SOURCE) $(DEP_CPP_O2M_C) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\o2m.sbr" : $(SOURCE) $(DEP_CPP_O2M_C) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\o2m.obj" : $(SOURCE) $(DEP_CPP_O2M_C) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wnq\tqdebug.c
DEP_CPP_TQDEB=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\table.h"\
	"\acedb\wh\tq_.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\tqdebug.obj" : $(SOURCE) $(DEP_CPP_TQDEB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tqdebug.obj" : $(SOURCE) $(DEP_CPP_TQDEB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tqdebug.sbr" : $(SOURCE) $(DEP_CPP_TQDEB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tqdebug.obj" : $(SOURCE) $(DEP_CPP_TQDEB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tqdebug.sbr" : $(SOURCE) $(DEP_CPP_TQDEB) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\tqdebug.obj" : $(SOURCE) $(DEP_CPP_TQDEB) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wnq\nqcpass2.c
DEP_CPP_NQCPA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\nqc.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\nqcpass2.obj" : $(SOURCE) $(DEP_CPP_NQCPA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nqcpass2.obj" : $(SOURCE) $(DEP_CPP_NQCPA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nqcpass2.sbr" : $(SOURCE) $(DEP_CPP_NQCPA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nqcpass2.obj" : $(SOURCE) $(DEP_CPP_NQCPA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nqcpass2.sbr" : $(SOURCE) $(DEP_CPP_NQCPA) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\nqcpass2.obj" : $(SOURCE) $(DEP_CPP_NQCPA) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wnq\nqctools.c
DEP_CPP_NQCTO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\nqc.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\nqctools.obj" : $(SOURCE) $(DEP_CPP_NQCTO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nqctools.obj" : $(SOURCE) $(DEP_CPP_NQCTO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nqctools.sbr" : $(SOURCE) $(DEP_CPP_NQCTO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nqctools.obj" : $(SOURCE) $(DEP_CPP_NQCTO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nqctools.sbr" : $(SOURCE) $(DEP_CPP_NQCTO) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\nqctools.obj" : $(SOURCE) $(DEP_CPP_NQCTO) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wnq\table.c
DEP_CPP_TABLE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\java.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\sysclass.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\table.sbr" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\table.sbr" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\table.sbr" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\table.sbr" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 0

"$(INTDIR)\table.obj" : $(SOURCE) $(DEP_CPP_TABLE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wnq\nqcpass1.c
DEP_CPP_NQCPAS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\nqc.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\nqcpass1.obj" : $(SOURCE) $(DEP_CPP_NQCPAS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX /Yc /Yu

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nqcpass1.obj" : $(SOURCE) $(DEP_CPP_NQCPAS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nqcpass1.sbr" : $(SOURCE) $(DEP_CPP_NQCPAS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\nqcpass1.obj" : $(SOURCE) $(DEP_CPP_NQCPAS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nqcpass1.sbr" : $(SOURCE) $(DEP_CPP_NQCPAS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\nqcpass1.obj" : $(SOURCE) $(DEP_CPP_NQCPAS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\rpcace_s.c
DEP_CPP_RPCAC=\
	"..\wdce\rpcace.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\rpcace_s.obj" : $(SOURCE) $(DEP_CPP_RPCAC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\rpcace_s.sbr" : $(SOURCE) $(DEP_CPP_RPCAC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\rpcace_s.obj" : $(SOURCE) $(DEP_CPP_RPCAC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\dceprot.c

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\dceserverlib.c
DEP_CPP_DCESE=\
	"..\wdce\dceprot.c"\
	"..\wdce\rpcace.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dceserverlib.obj" : $(SOURCE) $(DEP_CPP_DCESE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dceserverlib.sbr" : $(SOURCE) $(DEP_CPP_DCESE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\dceserverlib.obj" : $(SOURCE) $(DEP_CPP_DCESE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\rpcace_c.c
DEP_CPP_RPCACE=\
	"..\wdce\rpcace.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\rpcace_c.obj" : $(SOURCE) $(DEP_CPP_RPCACE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\rpcace_c.sbr" : $(SOURCE) $(DEP_CPP_RPCACE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"


"$(INTDIR)\rpcace_c.obj" : $(SOURCE) $(DEP_CPP_RPCACE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\rpcace_c.obj" : $(SOURCE) $(DEP_CPP_RPCACE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\rpcace_c.sbr" : $(SOURCE) $(DEP_CPP_RPCACE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"


"$(INTDIR)\rpcace_c.obj" : $(SOURCE) $(DEP_CPP_RPCACE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\dceclientlib.c
DEP_CPP_DCECL=\
	"..\wdce\dceprot.c"\
	"..\wdce\rpcace.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dceclientlib.obj" : $(SOURCE) $(DEP_CPP_DCECL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dceclientlib.sbr" : $(SOURCE) $(DEP_CPP_DCECL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"


"$(INTDIR)\dceclientlib.obj" : $(SOURCE) $(DEP_CPP_DCECL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\dceclientlib.obj" : $(SOURCE) $(DEP_CPP_DCECL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dceclientlib.sbr" : $(SOURCE) $(DEP_CPP_DCECL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"


"$(INTDIR)\dceclientlib.obj" : $(SOURCE) $(DEP_CPP_DCECL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 0
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\rpcace.idl

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# Begin Custom Build
InputDir=\acedb\wdce
InputPath=\acedb\wdce\rpcace.idl
InputName=rpcace

BuildCmds= \
	$(MSDEVDIR)\bin\midl /align 4 /osf /acf $(InputDir)\$(InputName).acf /out\
                                    $(InputDir) $(InputPath) \
	

"rpcace_c.c" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   $(BuildCmds)

"rpcace.h" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# Begin Custom Build
InputDir=\acedb\wdce
InputPath=\acedb\wdce\rpcace.idl
InputName=rpcace

BuildCmds= \
	$(MSDEVDIR)\bin\midl /align 4 /osf /acf $(InputDir)\$(InputName).acf /out\
                                    $(InputDir) $(InputPath) \
	

"rpcace_c.c" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   $(BuildCmds)

"rpcace.h" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wrpc\aceclient.c
DEP_CPP_ACECL=\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\aceclient.obj" : $(SOURCE) $(DEP_CPP_ACECL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\aceclient.sbr" : $(SOURCE) $(DEP_CPP_ACECL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"


"$(INTDIR)\aceclient.obj" : $(SOURCE) $(DEP_CPP_ACECL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\filqng.c
DEP_CPP_FILQN=\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filqng.sbr" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filqng.sbr" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filqng.sbr" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"


"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filqng.sbr" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"


"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filqng.sbr" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"


"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w1\messclean.c
DEP_CPP_MESSC=\
	"..\w1\messubs.c"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messclean.sbr" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messclean.sbr" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"


"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messclean.sbr" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"


"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messclean.sbr" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"


"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32util.cpp
DEP_CPP_WIN32U=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32util.sbr" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32util.sbr" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32util.sbr" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32util.sbr" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32util.sbr" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32util.sbr" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32util.obj" : $(SOURCE) $(DEP_CPP_WIN32U) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wrpc\aceserver.c

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\aceserver.obj" : $(SOURCE) $(DEP_CPP_ACESE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\aceserver.sbr" : $(SOURCE) $(DEP_CPP_ACESE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

"$(INTDIR)\aceserver.obj" : $(SOURCE) $(DEP_CPP_ACESE) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
DEP_CPP_ACESE=\
	"\acedb\wh\a.h"\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wdce\win32dcelib.cpp
DEP_CPP_WIN32D=\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32.h"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\win32dcelib.obj" : $(SOURCE) $(DEP_CPP_WIN32D) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\win32dcelib.sbr" : $(SOURCE) $(DEP_CPP_WIN32D) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\win32dcelib.obj" : $(SOURCE) $(DEP_CPP_WIN32D) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wrpc\netclient.c
DEP_CPP_NETCL=\
	"\acedb\wh\aceclient.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\netclient.obj" : $(SOURCE) $(DEP_CPP_NETCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\netclient.sbr" : $(SOURCE) $(DEP_CPP_NETCL) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"


"$(INTDIR)\netclient.obj" : $(SOURCE) $(DEP_CPP_NETCL) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\tacemain.c
DEP_CPP_TACEM=\
	"\acedb\wh\a.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\command.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\table.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# PROP Exclude_From_Build 1
# SUBTRACT CPP /YX

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\tacemain.obj" : $(SOURCE) $(DEP_CPP_TACEM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tacemain.sbr" : $(SOURCE) $(DEP_CPP_TACEM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"


"$(INTDIR)\tacemain.obj" : $(SOURCE) $(DEP_CPP_TACEM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\win32\win32thread.cpp

!IF  "$(CFG)" == "acedb - Win32 Release"

DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32thread.obj" : $(SOURCE) $(DEP_CPP_WIN32T) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32thread.obj" : $(SOURCE) $(DEP_CPP_WIN32T) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32thread.sbr" : $(SOURCE) $(DEP_CPP_WIN32T) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE) \
	

"$(INTDIR)\win32thread.obj" : $(SOURCE) $(DEP_CPP_WIN32T) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

"$(INTDIR)\win32thread.sbr" : $(SOURCE) $(DEP_CPP_WIN32T) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

"$(INTDIR)\win32thread.obj" : $(SOURCE) $(DEP_CPP_WIN32T) "$(INTDIR)"\
 "$(INTDIR)\Acedb.pch"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /Yu"stdafx.h" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 1
DEP_CPP_WIN32T=\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\win32h\parsedata.h"\
	"\acedb\win32h\stdafx.h"\
	"\acedb\win32h\win32thread.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
# ADD CPP /Yu"stdafx.h"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w4\prefsubs.c
DEP_CPP_PREFS=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\prefsubs.sbr" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\prefsubs.sbr" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 0

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\prefsubs.sbr" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\prefsubs.sbr" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 0

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 0

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MT /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\prefsubs.sbr" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP Exclude_From_Build 0

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP Exclude_From_Build 0

BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\prefsubs.sbr" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 0

"$(INTDIR)\prefsubs.obj" : $(SOURCE) $(DEP_CPP_PREFS) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\w7\fmapmenes.c
DEP_CPP_FMAPM=\
	"..\wh\bitset.h"\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\fmap.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\map.h"\
	"\acedb\wh\method.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\regular.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\fmapmenes.obj" : $(SOURCE) $(DEP_CPP_FMAPM) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\fmapmenes.obj" : $(SOURCE) $(DEP_CPP_FMAPM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\fmapmenes.sbr" : $(SOURCE) $(DEP_CPP_FMAPM) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\acedb\wnq\flag.c
DEP_CPP_FLAG_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\dict.h"\
	"\acedb\wh\flag.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

!IF  "$(CFG)" == "acedb - Win32 Release"


"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 Debug"

# SUBTRACT CPP /YX

BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/" /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\flag.sbr" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_WINDOWS" /D "ACEDB4" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\flag.sbr" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 Acembly Release"


"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_WINDOWS" /D "ACEDB4" /D\
 "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)/Acedb.pch" /YX /Fo"$(INTDIR)/"\
 /Fd"$(INTDIR)/" /c $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MDd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "ACEDB4" /D "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\flag.sbr" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 AceServer Release"


"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MD /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D "ACEDB4" /D\
 "NON_GRAPHIC" /D "_AFXDLL" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "acedb - Win32 AceClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 NetClient Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace debug"


BuildCmds= \
	$(CPP) /nologo /Zp4 /MTd /Gm /Gi /GR /GX /Zi /Od /Gy /I "\acedb\wh" /I\
 "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "ACEDB" /D "WIN32" /D\
 "_CONSOLE" /D "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /FR"$(INTDIR)/"\
 /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c $(SOURCE) \
	

"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\flag.sbr" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "acedb - Win32 WinTace Release"


"$(INTDIR)\flag.obj" : $(SOURCE) $(DEP_CPP_FLAG_) "$(INTDIR)"
   $(CPP) /nologo /Zp4 /MTd /Gi /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "ACEDB" /D "WIN32" /D "_CONSOLE" /D\
 "NON_GRAPHIC" /D "ACEDB4" /D "_MBCS" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c\
 $(SOURCE)


!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
 

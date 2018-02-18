# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# $Id: wintace.mak,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=WinTace - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to WinTace - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "WinTace - Win32 Release" && "$(CFG)" !=\
 "WinTace - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "wintace.mak" CFG="WinTace - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "WinTace - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "WinTace - Win32 Debug" (based on "Win32 (x86) Console Application")
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
# PROP Target_Last_Scanned "WinTace - Win32 Debug"
RSC=rc.exe
CPP=cl.exe

!IF  "$(CFG)" == "WinTace - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "C:\acedb\bin\WinTace"
# PROP Intermediate_Dir "C:\acedb\bin\WinTace"
# PROP Target_Dir ""
OUTDIR=C:\acedb\bin\WinTace
INTDIR=C:\acedb\bin\WinTace

ALL : "$(OUTDIR)\wintace.exe"

CLEAN : 
	-@erase "C:\acedb\bin\WinTace\wintace.exe"
	-@erase "C:\acedb\bin\WinTace\freeout.obj"
	-@erase "C:\acedb\bin\WinTace\lexalpha.obj"
	-@erase "C:\acedb\bin\WinTace\check.obj"
	-@erase "C:\acedb\bin\WinTace\command.obj"
	-@erase "C:\acedb\bin\WinTace\biblio.obj"
	-@erase "C:\acedb\bin\WinTace\filqng.obj"
	-@erase "C:\acedb\bin\WinTace\picksubs.obj"
	-@erase "C:\acedb\bin\WinTace\randsubs.obj"
	-@erase "C:\acedb\bin\WinTace\nicedump.obj"
	-@erase "C:\acedb\bin\WinTace\bsdumps.obj"
	-@erase "C:\acedb\bin\WinTace\longtext.obj"
	-@erase "C:\acedb\bin\WinTace\arraysub.obj"
	-@erase "C:\acedb\bin\WinTace\quovadis.obj"
	-@erase "C:\acedb\bin\WinTace\oldhelp.obj"
	-@erase "C:\acedb\bin\WinTace\keysetdump.obj"
	-@erase "C:\acedb\bin\WinTace\tags.obj"
	-@erase "C:\acedb\bin\WinTace\dict.obj"
	-@erase "C:\acedb\bin\WinTace\dnasubs.obj"
	-@erase "C:\acedb\bin\WinTace\bssubs.obj"
	-@erase "C:\acedb\bin\WinTace\bstree.obj"
	-@erase "C:\acedb\bin\WinTace\call.obj"
	-@erase "C:\acedb\bin\WinTace\linkdate.obj"
	-@erase "C:\acedb\bin\WinTace\objcache.obj"
	-@erase "C:\acedb\bin\WinTace\class.obj"
	-@erase "C:\acedb\bin\WinTace\alignment.obj"
	-@erase "C:\acedb\bin\WinTace\parse.obj"
	-@erase "C:\acedb\bin\WinTace\blocksub.obj"
	-@erase "C:\acedb\bin\WinTace\model.obj"
	-@erase "C:\acedb\bin\WinTace\asubs.obj"
	-@erase "C:\acedb\bin\WinTace\lexsubs.obj"
	-@erase "C:\acedb\bin\WinTace\keyset.obj"
	-@erase "C:\acedb\bin\WinTace\bump.obj"
	-@erase "C:\acedb\bin\WinTace\chrono.obj"
	-@erase "C:\acedb\bin\WinTace\freesubs.obj"
	-@erase "C:\acedb\bin\WinTace\tacemain.obj"
	-@erase "C:\acedb\bin\WinTace\disknew.obj"
	-@erase "C:\acedb\bin\WinTace\lexsubs4.obj"
	-@erase "C:\acedb\bin\WinTace\dump.obj"
	-@erase "C:\acedb\bin\WinTace\basepad.obj"
	-@erase "C:\acedb\bin\WinTace\sysclass.obj"
	-@erase "C:\acedb\bin\WinTace\timesubs.obj"
	-@erase "C:\acedb\bin\WinTace\filsubs.obj"
	-@erase "C:\acedb\bin\WinTace\bsubs.obj"
	-@erase "C:\acedb\bin\WinTace\sprdop.obj"
	-@erase "C:\acedb\bin\WinTace\banner.obj"
	-@erase "C:\acedb\bin\WinTace\help.obj"
	-@erase "C:\acedb\bin\WinTace\peptide.obj"
	-@erase "C:\acedb\bin\WinTace\bstools.obj"
	-@erase "C:\acedb\bin\WinTace\messclean.obj"
	-@erase "C:\acedb\bin\WinTace\queryexe.obj"
	-@erase "C:\acedb\bin\WinTace\session.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /YX /c
# ADD CPP /nologo /Zp4 /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "ACEDB" /D "NON_GRAPHIC" /c
# SUBTRACT CPP /YX
CPP_PROJ=/nologo /Zp4 /ML /GX /O2 /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "ACEDB" /D\
 "NON_GRAPHIC" /Fo"$(INTDIR)/" /c 
CPP_OBJS=C:\acedb\bin\WinTace/
CPP_SBRS=
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/wintace.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/wintace.pdb" /machine:I386 /out:"$(OUTDIR)/wintace.exe" 
LINK32_OBJS= \
	"$(INTDIR)/freeout.obj" \
	"$(INTDIR)/lexalpha.obj" \
	"$(INTDIR)/check.obj" \
	"$(INTDIR)/command.obj" \
	"$(INTDIR)/biblio.obj" \
	"$(INTDIR)/filqng.obj" \
	"$(INTDIR)/picksubs.obj" \
	"$(INTDIR)/randsubs.obj" \
	"$(INTDIR)/nicedump.obj" \
	"$(INTDIR)/bsdumps.obj" \
	"$(INTDIR)/longtext.obj" \
	"$(INTDIR)/arraysub.obj" \
	"$(INTDIR)/quovadis.obj" \
	"$(INTDIR)/oldhelp.obj" \
	"$(INTDIR)/keysetdump.obj" \
	"$(INTDIR)/tags.obj" \
	"$(INTDIR)/dict.obj" \
	"$(INTDIR)/dnasubs.obj" \
	"$(INTDIR)/bssubs.obj" \
	"$(INTDIR)/bstree.obj" \
	"$(INTDIR)/call.obj" \
	"$(INTDIR)/linkdate.obj" \
	"$(INTDIR)/objcache.obj" \
	"$(INTDIR)/class.obj" \
	"$(INTDIR)/alignment.obj" \
	"$(INTDIR)/parse.obj" \
	"$(INTDIR)/blocksub.obj" \
	"$(INTDIR)/model.obj" \
	"$(INTDIR)/asubs.obj" \
	"$(INTDIR)/lexsubs.obj" \
	"$(INTDIR)/keyset.obj" \
	"$(INTDIR)/bump.obj" \
	"$(INTDIR)/chrono.obj" \
	"$(INTDIR)/freesubs.obj" \
	"$(INTDIR)/tacemain.obj" \
	"$(INTDIR)/disknew.obj" \
	"$(INTDIR)/lexsubs4.obj" \
	"$(INTDIR)/dump.obj" \
	"$(INTDIR)/basepad.obj" \
	"$(INTDIR)/sysclass.obj" \
	"$(INTDIR)/timesubs.obj" \
	"$(INTDIR)/filsubs.obj" \
	"$(INTDIR)/bsubs.obj" \
	"$(INTDIR)/sprdop.obj" \
	"$(INTDIR)/banner.obj" \
	"$(INTDIR)/help.obj" \
	"$(INTDIR)/peptide.obj" \
	"$(INTDIR)/bstools.obj" \
	"$(INTDIR)/messclean.obj" \
	"$(INTDIR)/queryexe.obj" \
	"$(INTDIR)/session.obj"

"$(OUTDIR)\wintace.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "C:\acedb\bin\WinTace"
# PROP Intermediate_Dir "C:\acedb\bin\WinTace"
# PROP Target_Dir ""
OUTDIR=C:\acedb\bin\WinTace
INTDIR=C:\acedb\bin\WinTace

ALL : "$(OUTDIR)\wintace.exe" "$(OUTDIR)\wintace.bsc"

CLEAN : 
	-@erase "C:\acedb\bin\WinTace\vc40.pdb"
	-@erase "C:\acedb\bin\WinTace\vc40.idb"
	-@erase "C:\acedb\bin\WinTace\wintace.bsc"
	-@erase "C:\acedb\bin\WinTace\quovadis.sbr"
	-@erase "C:\acedb\bin\WinTace\oldhelp.sbr"
	-@erase "C:\acedb\bin\WinTace\keysetdump.sbr"
	-@erase "C:\acedb\bin\WinTace\tags.sbr"
	-@erase "C:\acedb\bin\WinTace\dict.sbr"
	-@erase "C:\acedb\bin\WinTace\dnasubs.sbr"
	-@erase "C:\acedb\bin\WinTace\bssubs.sbr"
	-@erase "C:\acedb\bin\WinTace\bstree.sbr"
	-@erase "C:\acedb\bin\WinTace\call.sbr"
	-@erase "C:\acedb\bin\WinTace\linkdate.sbr"
	-@erase "C:\acedb\bin\WinTace\objcache.sbr"
	-@erase "C:\acedb\bin\WinTace\class.sbr"
	-@erase "C:\acedb\bin\WinTace\alignment.sbr"
	-@erase "C:\acedb\bin\WinTace\parse.sbr"
	-@erase "C:\acedb\bin\WinTace\blocksub.sbr"
	-@erase "C:\acedb\bin\WinTace\model.sbr"
	-@erase "C:\acedb\bin\WinTace\asubs.sbr"
	-@erase "C:\acedb\bin\WinTace\lexsubs.sbr"
	-@erase "C:\acedb\bin\WinTace\keyset.sbr"
	-@erase "C:\acedb\bin\WinTace\bump.sbr"
	-@erase "C:\acedb\bin\WinTace\chrono.sbr"
	-@erase "C:\acedb\bin\WinTace\freesubs.sbr"
	-@erase "C:\acedb\bin\WinTace\tacemain.sbr"
	-@erase "C:\acedb\bin\WinTace\disknew.sbr"
	-@erase "C:\acedb\bin\WinTace\lexsubs4.sbr"
	-@erase "C:\acedb\bin\WinTace\dump.sbr"
	-@erase "C:\acedb\bin\WinTace\basepad.sbr"
	-@erase "C:\acedb\bin\WinTace\sysclass.sbr"
	-@erase "C:\acedb\bin\WinTace\timesubs.sbr"
	-@erase "C:\acedb\bin\WinTace\filsubs.sbr"
	-@erase "C:\acedb\bin\WinTace\bsubs.sbr"
	-@erase "C:\acedb\bin\WinTace\sprdop.sbr"
	-@erase "C:\acedb\bin\WinTace\banner.sbr"
	-@erase "C:\acedb\bin\WinTace\help.sbr"
	-@erase "C:\acedb\bin\WinTace\peptide.sbr"
	-@erase "C:\acedb\bin\WinTace\bstools.sbr"
	-@erase "C:\acedb\bin\WinTace\messclean.sbr"
	-@erase "C:\acedb\bin\WinTace\queryexe.sbr"
	-@erase "C:\acedb\bin\WinTace\session.sbr"
	-@erase "C:\acedb\bin\WinTace\freeout.sbr"
	-@erase "C:\acedb\bin\WinTace\lexalpha.sbr"
	-@erase "C:\acedb\bin\WinTace\check.sbr"
	-@erase "C:\acedb\bin\WinTace\command.sbr"
	-@erase "C:\acedb\bin\WinTace\biblio.sbr"
	-@erase "C:\acedb\bin\WinTace\filqng.sbr"
	-@erase "C:\acedb\bin\WinTace\picksubs.sbr"
	-@erase "C:\acedb\bin\WinTace\randsubs.sbr"
	-@erase "C:\acedb\bin\WinTace\nicedump.sbr"
	-@erase "C:\acedb\bin\WinTace\bsdumps.sbr"
	-@erase "C:\acedb\bin\WinTace\longtext.sbr"
	-@erase "C:\acedb\bin\WinTace\arraysub.sbr"
	-@erase "C:\acedb\bin\WinTace\wintace.exe"
	-@erase "C:\acedb\bin\WinTace\freeout.obj"
	-@erase "C:\acedb\bin\WinTace\lexalpha.obj"
	-@erase "C:\acedb\bin\WinTace\check.obj"
	-@erase "C:\acedb\bin\WinTace\command.obj"
	-@erase "C:\acedb\bin\WinTace\biblio.obj"
	-@erase "C:\acedb\bin\WinTace\filqng.obj"
	-@erase "C:\acedb\bin\WinTace\picksubs.obj"
	-@erase "C:\acedb\bin\WinTace\randsubs.obj"
	-@erase "C:\acedb\bin\WinTace\nicedump.obj"
	-@erase "C:\acedb\bin\WinTace\bsdumps.obj"
	-@erase "C:\acedb\bin\WinTace\longtext.obj"
	-@erase "C:\acedb\bin\WinTace\arraysub.obj"
	-@erase "C:\acedb\bin\WinTace\quovadis.obj"
	-@erase "C:\acedb\bin\WinTace\oldhelp.obj"
	-@erase "C:\acedb\bin\WinTace\keysetdump.obj"
	-@erase "C:\acedb\bin\WinTace\tags.obj"
	-@erase "C:\acedb\bin\WinTace\dict.obj"
	-@erase "C:\acedb\bin\WinTace\dnasubs.obj"
	-@erase "C:\acedb\bin\WinTace\bssubs.obj"
	-@erase "C:\acedb\bin\WinTace\bstree.obj"
	-@erase "C:\acedb\bin\WinTace\call.obj"
	-@erase "C:\acedb\bin\WinTace\linkdate.obj"
	-@erase "C:\acedb\bin\WinTace\objcache.obj"
	-@erase "C:\acedb\bin\WinTace\class.obj"
	-@erase "C:\acedb\bin\WinTace\alignment.obj"
	-@erase "C:\acedb\bin\WinTace\parse.obj"
	-@erase "C:\acedb\bin\WinTace\blocksub.obj"
	-@erase "C:\acedb\bin\WinTace\model.obj"
	-@erase "C:\acedb\bin\WinTace\asubs.obj"
	-@erase "C:\acedb\bin\WinTace\lexsubs.obj"
	-@erase "C:\acedb\bin\WinTace\keyset.obj"
	-@erase "C:\acedb\bin\WinTace\bump.obj"
	-@erase "C:\acedb\bin\WinTace\chrono.obj"
	-@erase "C:\acedb\bin\WinTace\freesubs.obj"
	-@erase "C:\acedb\bin\WinTace\tacemain.obj"
	-@erase "C:\acedb\bin\WinTace\disknew.obj"
	-@erase "C:\acedb\bin\WinTace\lexsubs4.obj"
	-@erase "C:\acedb\bin\WinTace\dump.obj"
	-@erase "C:\acedb\bin\WinTace\basepad.obj"
	-@erase "C:\acedb\bin\WinTace\sysclass.obj"
	-@erase "C:\acedb\bin\WinTace\timesubs.obj"
	-@erase "C:\acedb\bin\WinTace\filsubs.obj"
	-@erase "C:\acedb\bin\WinTace\bsubs.obj"
	-@erase "C:\acedb\bin\WinTace\sprdop.obj"
	-@erase "C:\acedb\bin\WinTace\banner.obj"
	-@erase "C:\acedb\bin\WinTace\help.obj"
	-@erase "C:\acedb\bin\WinTace\peptide.obj"
	-@erase "C:\acedb\bin\WinTace\bstools.obj"
	-@erase "C:\acedb\bin\WinTace\messclean.obj"
	-@erase "C:\acedb\bin\WinTace\queryexe.obj"
	-@erase "C:\acedb\bin\WinTace\session.obj"
	-@erase "C:\acedb\bin\WinTace\wintace.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /YX /c
# ADD CPP /nologo /Zp4 /Gm /GX /Zi /Od /I "\acedb\wh" /I "\acedb\whooks" /I "\acedb\win32h" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "ACEDB" /D "NON_GRAPHIC" /FR /c
# SUBTRACT CPP /YX
CPP_PROJ=/nologo /Zp4 /MLd /Gm /GX /Zi /Od /I "\acedb\wh" /I "\acedb\whooks" /I\
 "\acedb\win32h" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "ACEDB" /D\
 "NON_GRAPHIC" /FR"$(INTDIR)/" /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=C:\acedb\bin\WinTace/
CPP_SBRS=C:\acedb\bin\WinTace/
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/wintace.bsc" 
BSC32_SBRS= \
	"$(INTDIR)/quovadis.sbr" \
	"$(INTDIR)/oldhelp.sbr" \
	"$(INTDIR)/keysetdump.sbr" \
	"$(INTDIR)/tags.sbr" \
	"$(INTDIR)/dict.sbr" \
	"$(INTDIR)/dnasubs.sbr" \
	"$(INTDIR)/bssubs.sbr" \
	"$(INTDIR)/bstree.sbr" \
	"$(INTDIR)/call.sbr" \
	"$(INTDIR)/linkdate.sbr" \
	"$(INTDIR)/objcache.sbr" \
	"$(INTDIR)/class.sbr" \
	"$(INTDIR)/alignment.sbr" \
	"$(INTDIR)/parse.sbr" \
	"$(INTDIR)/blocksub.sbr" \
	"$(INTDIR)/model.sbr" \
	"$(INTDIR)/asubs.sbr" \
	"$(INTDIR)/lexsubs.sbr" \
	"$(INTDIR)/keyset.sbr" \
	"$(INTDIR)/bump.sbr" \
	"$(INTDIR)/chrono.sbr" \
	"$(INTDIR)/freesubs.sbr" \
	"$(INTDIR)/tacemain.sbr" \
	"$(INTDIR)/disknew.sbr" \
	"$(INTDIR)/lexsubs4.sbr" \
	"$(INTDIR)/dump.sbr" \
	"$(INTDIR)/basepad.sbr" \
	"$(INTDIR)/sysclass.sbr" \
	"$(INTDIR)/timesubs.sbr" \
	"$(INTDIR)/filsubs.sbr" \
	"$(INTDIR)/bsubs.sbr" \
	"$(INTDIR)/sprdop.sbr" \
	"$(INTDIR)/banner.sbr" \
	"$(INTDIR)/help.sbr" \
	"$(INTDIR)/peptide.sbr" \
	"$(INTDIR)/bstools.sbr" \
	"$(INTDIR)/messclean.sbr" \
	"$(INTDIR)/queryexe.sbr" \
	"$(INTDIR)/session.sbr" \
	"$(INTDIR)/freeout.sbr" \
	"$(INTDIR)/lexalpha.sbr" \
	"$(INTDIR)/check.sbr" \
	"$(INTDIR)/command.sbr" \
	"$(INTDIR)/biblio.sbr" \
	"$(INTDIR)/filqng.sbr" \
	"$(INTDIR)/picksubs.sbr" \
	"$(INTDIR)/randsubs.sbr" \
	"$(INTDIR)/nicedump.sbr" \
	"$(INTDIR)/bsdumps.sbr" \
	"$(INTDIR)/longtext.sbr" \
	"$(INTDIR)/arraysub.sbr"

"$(OUTDIR)\wintace.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:4.31 /subsystem:console /profile /map /debug /machine:I386
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /version:4.31 /subsystem:console /profile\
 /map:"$(INTDIR)/wintace.map" /debug /machine:I386 /out:"$(OUTDIR)/wintace.exe" 
LINK32_OBJS= \
	"$(INTDIR)/freeout.obj" \
	"$(INTDIR)/lexalpha.obj" \
	"$(INTDIR)/check.obj" \
	"$(INTDIR)/command.obj" \
	"$(INTDIR)/biblio.obj" \
	"$(INTDIR)/filqng.obj" \
	"$(INTDIR)/picksubs.obj" \
	"$(INTDIR)/randsubs.obj" \
	"$(INTDIR)/nicedump.obj" \
	"$(INTDIR)/bsdumps.obj" \
	"$(INTDIR)/longtext.obj" \
	"$(INTDIR)/arraysub.obj" \
	"$(INTDIR)/quovadis.obj" \
	"$(INTDIR)/oldhelp.obj" \
	"$(INTDIR)/keysetdump.obj" \
	"$(INTDIR)/tags.obj" \
	"$(INTDIR)/dict.obj" \
	"$(INTDIR)/dnasubs.obj" \
	"$(INTDIR)/bssubs.obj" \
	"$(INTDIR)/bstree.obj" \
	"$(INTDIR)/call.obj" \
	"$(INTDIR)/linkdate.obj" \
	"$(INTDIR)/objcache.obj" \
	"$(INTDIR)/class.obj" \
	"$(INTDIR)/alignment.obj" \
	"$(INTDIR)/parse.obj" \
	"$(INTDIR)/blocksub.obj" \
	"$(INTDIR)/model.obj" \
	"$(INTDIR)/asubs.obj" \
	"$(INTDIR)/lexsubs.obj" \
	"$(INTDIR)/keyset.obj" \
	"$(INTDIR)/bump.obj" \
	"$(INTDIR)/chrono.obj" \
	"$(INTDIR)/freesubs.obj" \
	"$(INTDIR)/tacemain.obj" \
	"$(INTDIR)/disknew.obj" \
	"$(INTDIR)/lexsubs4.obj" \
	"$(INTDIR)/dump.obj" \
	"$(INTDIR)/basepad.obj" \
	"$(INTDIR)/sysclass.obj" \
	"$(INTDIR)/timesubs.obj" \
	"$(INTDIR)/filsubs.obj" \
	"$(INTDIR)/bsubs.obj" \
	"$(INTDIR)/sprdop.obj" \
	"$(INTDIR)/banner.obj" \
	"$(INTDIR)/help.obj" \
	"$(INTDIR)/peptide.obj" \
	"$(INTDIR)/bstools.obj" \
	"$(INTDIR)/messclean.obj" \
	"$(INTDIR)/queryexe.obj" \
	"$(INTDIR)/session.obj"

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

# Name "WinTace - Win32 Release"
# Name "WinTace - Win32 Debug"

!IF  "$(CFG)" == "WinTace - Win32 Release"

!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\tacemain.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_TACEM=\
	"\acedb\wh\acedb.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\spread.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	

"$(INTDIR)\tacemain.obj" : $(SOURCE) $(DEP_CPP_TACEM) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_TACEM=\
	"\acedb\wh\acedb.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\spread.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\tacemain.obj" : $(SOURCE) $(DEP_CPP_TACEM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tacemain.sbr" : $(SOURCE) $(DEP_CPP_TACEM) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\dump.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\disk.h"\
	"\acedb\whooks\disptype.h"\
	
NODEP_CPP_DUMP_=\
	".\..\W4\MPWIncludes.h"\
	

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	".\..\W4\MPWIncludes.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\dump.obj" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dump.sbr" : $(SOURCE) $(DEP_CPP_DUMP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\session.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\session.obj" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\session.sbr" : $(SOURCE) $(DEP_CPP_SESSI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\model.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk.h"\
	
NODEP_CPP_MODEL=\
	".\..\w4\newmodel.c"\
	

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\model.obj" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\model.sbr" : $(SOURCE) $(DEP_CPP_MODEL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\parse.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\parse.obj" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\parse.sbr" : $(SOURCE) $(DEP_CPP_PARSE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\chrono.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_CHRON=\
	"\acedb\wh\acedb.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_CHRON=\
	"\acedb\wh\acedb.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\chrono.obj" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\chrono.sbr" : $(SOURCE) $(DEP_CPP_CHRON) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\longtext.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\longtext.obj" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\longtext.sbr" : $(SOURCE) $(DEP_CPP_LONGT) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\queryexe.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\array.h"\
	

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\queryexe.obj" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\queryexe.sbr" : $(SOURCE) $(DEP_CPP_QUERY) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\alignment.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_ALIGN=\
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
	

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_ALIGN=\
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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\alignment.obj" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\alignment.sbr" : $(SOURCE) $(DEP_CPP_ALIGN) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W7\biblio.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\biblio.obj" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\biblio.sbr" : $(SOURCE) $(DEP_CPP_BIBLI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\WHOOKS\sysclass.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_SYSCL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_SYSCL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\sysclass.obj" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sysclass.sbr" : $(SOURCE) $(DEP_CPP_SYSCL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\WHOOKS\class.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk_.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\class.obj" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\class.sbr" : $(SOURCE) $(DEP_CPP_CLASS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\WHOOKS\tags.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_TAGS_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_TAGS_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\tags.obj" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\tags.sbr" : $(SOURCE) $(DEP_CPP_TAGS_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\WHOOKS\quovadis.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\quovadis.obj" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\quovadis.sbr" : $(SOURCE) $(DEP_CPP_QUOVA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\lexsubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_LEXSU=\
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
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_LEXSU=\
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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\lexsubs.obj" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs.sbr" : $(SOURCE) $(DEP_CPP_LEXSU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\command.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_COMMA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_COMMA=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\dump.h"\
	"\acedb\wh\lex.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\whooks\tags.h"\
	"\acedb\wh\spread.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\biblio.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\command.obj" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\command.sbr" : $(SOURCE) $(DEP_CPP_COMMA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\lexsubs4.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_LEXSUB=\
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
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_LEXSUB=\
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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\lexsubs4.obj" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexsubs4.sbr" : $(SOURCE) $(DEP_CPP_LEXSUB) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\lexalpha.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_LEXAL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_LEXAL=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\a.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\lexalpha.obj" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lexalpha.sbr" : $(SOURCE) $(DEP_CPP_LEXAL) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\disknew.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\array.h"\
	
NODEP_CPP_DISKN=\
	".\..\W5\macfilemgrsubs.h"\
	

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	".\..\W5\macfilemgrsubs.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\disknew.obj" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\disknew.sbr" : $(SOURCE) $(DEP_CPP_DISKN) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\blocksub.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\disk__.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\blocksub.obj" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\blocksub.sbr" : $(SOURCE) $(DEP_CPP_BLOCK) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\objcache.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk__.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\disk__.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk_.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	"\acedb\wh\disk.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\objcache.obj" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\objcache.sbr" : $(SOURCE) $(DEP_CPP_OBJCA) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W5\keysetdump.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_KEYSE=\
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
	"\acedb\wh\array.h"\
	

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_KEYSE=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\pick.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\whooks\sysclass.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\keysetdump.obj" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keysetdump.sbr" : $(SOURCE) $(DEP_CPP_KEYSE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\bsubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\bsubs.obj" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsubs.sbr" : $(SOURCE) $(DEP_CPP_BSUBS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\bstree.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk_.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\bstree.obj" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstree.sbr" : $(SOURCE) $(DEP_CPP_BSTRE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\bstools.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_BSTOO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_BSTOO=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bs_.h"\
	"\acedb\wh\chrono.h"\
	"\acedb\wh\parse.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\bstools.obj" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bstools.sbr" : $(SOURCE) $(DEP_CPP_BSTOO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\bsdumps.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\bsdumps.obj" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bsdumps.sbr" : $(SOURCE) $(DEP_CPP_BSDUM) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\nicedump.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\nicedump.obj" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\nicedump.sbr" : $(SOURCE) $(DEP_CPP_NICED) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\dnasubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\dnasubs.obj" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dnasubs.sbr" : $(SOURCE) $(DEP_CPP_DNASU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\peptide.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\peptide.obj" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\peptide.sbr" : $(SOURCE) $(DEP_CPP_PEPTI) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\check.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_CHECK=\
	"\acedb\wh\acedb.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\check.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_CHECK=\
	"\acedb\wh\acedb.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\check.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\check.obj" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\check.sbr" : $(SOURCE) $(DEP_CPP_CHECK) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\basepad.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_BASEP=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_BASEP=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\lex.h"\
	"\acedb\wh\bs.h"\
	"\acedb\wh\query.h"\
	"\acedb\wh\dna.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\classes.h"\
	"\acedb\whooks\systags.h"\
	"\acedb\wh\a.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\basepad.obj" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\basepad.sbr" : $(SOURCE) $(DEP_CPP_BASEP) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\asubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\asubs.obj" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\asubs.sbr" : $(SOURCE) $(DEP_CPP_ASUBS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\bssubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\wh\bstree.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\bssubs.obj" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bssubs.sbr" : $(SOURCE) $(DEP_CPP_BSSUB) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\picksubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\array.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\picksubs.obj" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\picksubs.sbr" : $(SOURCE) $(DEP_CPP_PICKS) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\linkdate.c

!IF  "$(CFG)" == "WinTace - Win32 Release"


"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"


BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\linkdate.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\linkdate.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\keyset.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_KEYSET=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_KEYSET=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\keyset.obj" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\keyset.sbr" : $(SOURCE) $(DEP_CPP_KEYSET) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W6\sprdop.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

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
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\whooks\disptype.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	"\acedb\wh\bump.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\display.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\whooks\disptype.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\sprdop.obj" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\sprdop.sbr" : $(SOURCE) $(DEP_CPP_SPRDO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W4\banner.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_BANNE=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_BANNE=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\banner.obj" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\banner.sbr" : $(SOURCE) $(DEP_CPP_BANNE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\freesubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_FREES=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_FREES=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\freesubs.obj" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freesubs.sbr" : $(SOURCE) $(DEP_CPP_FREES) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\freeout.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_FREEO=\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_FREEO=\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\freeout.obj" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\freeout.sbr" : $(SOURCE) $(DEP_CPP_FREEO) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\messclean.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_MESSC=\
	".\..\W1\messubs.c"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	

"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_MESSC=\
	".\..\W1\messubs.c"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\messclean.obj" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\messclean.sbr" : $(SOURCE) $(DEP_CPP_MESSC) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\arraysub.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_ARRAY=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_ARRAY=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\arraysub.obj" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\arraysub.sbr" : $(SOURCE) $(DEP_CPP_ARRAY) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\filsubs.c
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
	

!IF  "$(CFG)" == "WinTace - Win32 Release"


"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"


BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\filsubs.obj" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filsubs.sbr" : $(SOURCE) $(DEP_CPP_FILSU) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\timesubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_TIMES=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_TIMES=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\timesubs.obj" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\timesubs.sbr" : $(SOURCE) $(DEP_CPP_TIMES) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\bump.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_BUMP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_BUMP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\bump.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\bump.obj" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\bump.sbr" : $(SOURCE) $(DEP_CPP_BUMP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\randsubs.c

!IF  "$(CFG)" == "WinTace - Win32 Release"


"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"


BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\randsubs.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\randsubs.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\call.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_CALL_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_CALL_=\
	"\acedb\wh\acedb.h"\
	"\acedb\wh\call.h"\
	"\acedb\wh\mytime.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\keyset.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\array.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\call.obj" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\call.sbr" : $(SOURCE) $(DEP_CPP_CALL_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\filqng.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_FILQN=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_FILQN=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\filqng.obj" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\filqng.sbr" : $(SOURCE) $(DEP_CPP_FILQN) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\dict.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_DICT_=\
	"\acedb\wh\dict.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_DICT_=\
	"\acedb\wh\dict.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\regular.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\dict.obj" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\dict.sbr" : $(SOURCE) $(DEP_CPP_DICT_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\help.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_HELP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mydirent.h"\
	"\acedb\wh\dict.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	".\..\wgd\gd.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	".\..\wh\ibmdir.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	
NODEP_CPP_HELP_=\
	".\..\wh\aixfs\dir.h"\
	".\..\wh\jfs\dir.h"\
	

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

DEP_CPP_HELP_=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\mydirent.h"\
	"\acedb\wh\dict.h"\
	{$(INCLUDE)}"\sys\STAT.H"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\graph.h"\
	".\..\wgd\gd.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\mystdlib.h"\
	{$(INCLUDE)}"\sys\TYPES.H"\
	".\..\wh\ibmdir.h"\
	
NODEP_CPP_HELP_=\
	".\..\wh\aixfs\dir.h"\
	".\..\wh\jfs\dir.h"\
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\help.obj" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\help.sbr" : $(SOURCE) $(DEP_CPP_HELP_) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\ACEDB\W1\oldhelp.c

!IF  "$(CFG)" == "WinTace - Win32 Release"

DEP_CPP_OLDHE=\
	"\acedb\wh\regular.h"\
	"\acedb\wh\array.h"\
	"\acedb\wh\graph.h"\
	"\acedb\wh\key.h"\
	"\acedb\wh\freeout.h"\
	"\acedb\wh\session.h"\
	"\acedb\wh\mystdlib.h"\
	"\acedb\wh\disk.h"\
	

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "WinTace - Win32 Debug"

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
	

BuildCmds= \
	$(CPP) $(CPP_PROJ) $(SOURCE) \
	

"$(INTDIR)\oldhelp.obj" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\oldhelp.sbr" : $(SOURCE) $(DEP_CPP_OLDHE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
 

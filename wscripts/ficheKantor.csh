#!/bin/csh -f

if ($3 == "" || ! -d $3) exit 1

set ACEBINDIR=~/ace/bin.$ACEDB_MACHINE
#set ACEBINDIR=~/ace/bin.SOLARIS_4_OPT
set ACESERVER="localhost:200103 "
set COMMANDRUN="$ACEBINDIR/kantor $ACESERVER $4 $5 $6 -out $2.out"
set SERVERRUN="$ACEBINDIR/aceserver $3 200103 3600:3600:0 &"
#echo "$COMMANDRUN"
#echo "$SERVERRUN"

echo "launching: $COMMANDRUN"
set RESULTRUN="`$COMMANDRUN`"
#echo "$RESULTRUN"

set NOCONNECTION="No database connection $ACESERVER"
#echo "$NOCONNECTION"

if ( "$RESULTRUN" == "$NOCONNECTION" ) then
#    echo "starting server"
    echo "$SERVERRUN"
    set RESULTRUN="`$SERVERRUN`"
    sleep 5
    echo "$COMMANDRUN"
    set RESULTRUN = "`$COMMANDRUN`"
endif

if ($1 != "0") kill -USR1 $1

exit 0



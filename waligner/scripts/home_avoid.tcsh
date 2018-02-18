#!/bin/tcsh

\rm bin metaData scripts
ln -s $MAGIC_SRC/bin.$ACEDB_MACHINE bin
ln -s $MAGIC_SRC/waligner/metaData
ln -s $MAGIC_SRC/waligner/scripts

pushd TARGET
set d1=`pwd`
popd
\rm TARGET
ln -s $d1 TARGET


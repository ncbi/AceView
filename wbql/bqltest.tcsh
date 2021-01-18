#!/bin/tcsh -f

echo "The Acedb query language, 2018, test-suite"
echo "Author: Jean Thierry-Mieg"
echo "For questions and information please email mieg@ncbi.nlm.nih.gov"

########################################
# Localize the executable
########################################
set tace=`ls bin*/tace | head | head -1`
set ACEDB_SRC
if (! -x $tace) then
  echo 'FATAL ERROR: the executable tace is not found in bin*/tace'
  echo 'Please compile the code'
  goto done
endif

########################################
# Create a test directory
########################################
# in wbql we ran: tar cf test.tar wspec test.ace test.out.expected
setenv BDIR  ACEDB_QUERY_LANGUAGE_TEST
if (! -d $BDIR) then
  mkdir $BDIR
  pushd $BDIR
  # create the scema in the wspec directory and add the current user
  mkdir database
  tar xf ../wbql/test.wspec.tar
  cp ../wbql/test.ace ../wbql/test.cmd ../wbql/test.out.expected .
  echo `whoami` >> wspec/passwd.wrm
  # initialize the database and parse the test data
  ../$tace . << EOF
y  // initialize the database
    pparse test.ace
    save
    quit
EOF
  popd
endif

########################################
# Run the test and compare to the expected output
########################################

pushd $BDIR
  ../$tace . < test.cmd > test.out
  diff test.out test.out.expected > test.diff
popd

done:
  echo done

if (-e $BDIR/test.diff) then
  set n=`wc -l $BDIR/test.diff | awk '{print $1}'`
  if ($n > 0) then
    echo 'Please check the file the file $BDIR/test.diff, it should have been empty'
    wc -l $BDIR/test.diff
    echo 'If the error is reproducible please report it to mieg@ncbi.nlm.nih.gov'
    echo 'Thank you for testing this code'
  else
    echo "ALL TESTS ARE SUCCESSFUL"
    echo "You may try additional query by typing the command"
    echo 'bin.*/tace ' $BDIR
  endif
endif



#!bin/tcsh -f

# Verify that the code has been installed and is accessible

if (! $?MAGIC_SRC || ! -d $MAGIC_SRC/waligner) then
  echo "please setenv MAGIC_SRC to correspond to the source code directory"
  echo "i.e. the directory containg the subdirectories w1, w2, w3...."
  exit 1
endif

if (! $?ACEDB_MACHINE || ! -d $MAGIC_SRC/bin.$ACEDB_MACHINE) then
   echo "please setenv ACEDB_MACHINE to correspond to your architecture and compile the code (make -k all)"
   exit 1
endif

# verify that the test data are accessible
set species=`ls *.test.ILM.R1.fastq.gz | gawk '{i=index($1,".");n++; if(n==1)z=substr($1,1,i-1);}END{printf("%s",substr($1,1,i-1)); }'`
if (! -e $species.test.ILM.R1.fastq.gz) then
   echo "File $species.rnaSeq.test.fastq.gz is not accessible, please download it."
   exit 1
endif
  
# create a test directory
if (! -d TEST) mkdir TEST
if (! -d TEST)then
  echo "Sorry i cannot create the TEST directory"
  exit 1
endif

echo "Testing for existence of : TARGET.$species.*"
ls -d TARGET.$species.*

set tt=`ls -d TARGET.$species.* | gawk '/gz/{next}{n++;if(n==1)print;}'`
if ($tt == "") set tt=_____
if ($tt == "" || ! -d $tt) then
  echo "Sorry, the TARGET.$species directory is not accessible"
  echo "please download and decompress the TARGET tar.gz file in this directory"
  echo "or create a link pointing to the TARGET dir"
  echo "Then retry, thank you"
  exit 1
endif

###############################
# Configure the test
cd TEST

# run the test suite as MULTICORE
# real data sets should be run as SGE, or the program 
# scripts/submit should be edited and adapted to your harware
setenv MAGIC_SUBMIT MULTICORE

# Intall the TARGET for the correct species by linking to the distribution directory
ln -s ../$tt TARGET

# Configure the directory to point for the compilation directory
$MAGIC_SRC/waligner/scripts/MAGIC init RNA
setenv MAGIC Test
 
# The file names should be declared in a subdirectory possibly using links
# but the name themselves shoudl not start with . or with .. 
mkdir DATA
cd DATA
  ln -s ../../human.test.ILM.R1.fastq.gz
  ln -s ../../human.test.ILM.R2.fastq.gz
cd ..

# declare the Test run using the acedb object oriented syntax (.ace file)
cat << EOF > test.ace
Run Test1
Project Test
RNA
Paired_end
Illumina HiSeq
nonStranded
Title "Test Illumina run with 1 million fragments distributed with the MAGIC pipeline"
Sample "Human_Mixed_tissues_PolyA_ns
Sequencing_laboratory CNL
-D File
File fastq/1 DATA/$species.test.ILM.R1.fastq.gz
File fastq/2 DATA/$species.test.ILM.R2.fastq.gz

EOF

# Parse the declared run in the MAGIC meta-data acedb-database
bin/tacembly MetaDB << EOF
  parse test.ace
  save
  quit
EOF

# Reconfigure the fastq file in fastc format

./MAGIC a0

# Run the alignment
./MAGIC wait ALIGN GENE_EXPRESSION

# Check that the results are correct
# not sure how to do this !!!!


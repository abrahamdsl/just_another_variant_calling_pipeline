#!/bin/bash -x

#
# Written/extended by a.llave@irri.org 16APR2015-0953
#
# A Bash script to align reads to a reference genome using "bwa" Bioinformatics tool.
#
# Derives from:
#   Example BWA alignment script
#   (c) Anna Battenhouse 21 May 2012
#    https://wikis.utexas.edu/display/bioiteam/Example+BWA+alignment+script
#
# (c) 2015 under the C4 Rice Project of the International Rice Research Institute (IRRI).
# Use outside of non-profit research might require legal compliance with IRRI.
# All intellectual property rights are to be determined by/remain with IRRI.

echo "--------------------------------";
echo "`date` : $0 triggered.";
echo " PID: $$";
echo "--------------------------------";
IN_FQ=$1
OUT_PFX=$2
ASSEMBLY=$3
PAIRED=$4
RGROUP=$5
SORT_BAM=$6

POST_BWA_SAM_OUTPUT_FILENAME=${IN_FQ/_1.fastq/.sam};
FIXEDMATE_NAME=${POST_BWA_SAM_OUTPUT_FILENAME/.sam/_fixmate.bam}

# samtools directory
SAMTOOLS_DIR_X="/home/applications/samtools-0.1.19/"

# where is your bwa?
BWA_DIR_X="/usr/bin"

COMPRESS_SAM_PID=-1

PRE_FIXMATE_NAME="";

# show usage if we don't have 4 command-line arguments
if [ "$PAIRED" == "" ]; then
    echo "-----------------------------------------------------------------";
    echo "Align fastq data with bwa, producing a sorted indexed BAM file.";
    echo " ";
    echo "align_bwa.sh in_file out_pfx reference paired(0|1) [0|1]";
    echo " ";
    echo "  in_file   For single-end alignments, path of the input fastq file.";
    echo "            For paired-end alignemtts, path to the the R1 fastq file"
    echo "            which must contain the string '_1.fq' in its name. The ";
    echo "            R2 file is thus expected to be named '_2.fq' automatically. ";
    echo "  out_pfx   Desired prefix of output files.";
    echo "  reference  Supply sequence file in absolute path ";
    echo "  paired    0 = single end alignment; 1 = paired end.";
    echo "  rg        Read group label.";
    echo "  sort_bam  [Optional] Also does fixmate/clean up read pairing information".
    echo "            Default is 0 - do not sort. 1 to sort.";
    echo " ";
    echo "Example:";
    echo "  align_bwa.sh my.fastq mrna_b1_ln1 hg18 0";
    echo "  align_bwa.sh my_L001_1.pe.fq swi6_b2_ln1 sacCer3 1";
    exit 1;
fi
# general function that exits after printing its text argument
#   in a standard format which can be easily grep'd.
err() {
  echo "$1...exiting";
  exit 1; # any non-0 exit code signals an error
}
# function to check return code of programs.
# exits with standard message if code is non-zero;
# otherwise displays completiong message and date.
#   arg 1 is the return code (usually $?)
#   arg2 is text describing what ran
ckRes() {
  if [ "$1" == "0" ]; then
    echo "..Done $2 `date`";
  else
    err "$2 returned non-0 exit code $1";
  fi
}
# function that checks if a file exists
#   arg 1 is the file name
#   arg2 is text describing the file (optional)
ckFile() {
  if [ ! -e "$1" ]; then
    err "$2 File '$1' not found";
  fi
}
# function that checks if a file exists and
#   that it has non-0 length. needed because
#   programs don't always return non-0 return
#   codes, and worse, they also create their
#   output file with 0 length so that just
#   checking for its existence is not enough
#   to ensure the program ran properly
ckFileSz() {
  ckFile $1 $2;
  SZ=`ls -l $1 | awk '{print $5}'`;
  if [ "$SZ" == "0" ]; then
    err "$2 file '$1' is zero length";
  fi
}

# --------------
# Defaulting
# --------------
# if no parameter specified for sorting the sam/bam output, then set to zero.
if [ "$SORT_BAM" == ""  ]; then
        SORT_BAM="0";
else
       if [ \( "$SORT_BAM" == "0" -o "$SORT_BAM" == "1" \) ne 1 ]; then
               echo "Invalid parameter supplied for sort_bam (parameter 6).";
               exit 1;
       fi
       echo "Will be not sorting BAM.";
fi
# If aligning read pairs, find the name of the R2 file
#   based on the R1 filename.
if [ "$PAIRED" == "1" ]; then
    IN_FQ_R1="$IN_FQ";
    IN_FQ_R2=${IN_FQ_R1/_1.fastq/_2.fastq}; # a.llave@irri.org : Modified to suit our needs
fi

# Note: reference indexes created with bwa 0.5.9 are *not*
#   compatible with bwa 0.6 or higher. So we go to some
#   trouble to figure out which BWA version we have
echo "Executing $BWA_DIR_X/bwa";
BWA_VER=`$BWA_DIR_X/bwa 2>&1 | grep Version | awk '{print $2}' | awk -F '.' '{print $1 "." $2}'`;
if [ "$BWA_VER" == "0.5" ]; then BWA_DIR=bwa;
elif [ "$BWA_VER" == "0.6" ]; then BWA_DIR=bwa6;
else BWA_DIR=unknown_bwa_version; fi

REF_PFX="$ASSEMBLY";

# Read group information is part of the SAM/BAM header that desribes
#   what is being aligned. When multiple lanes of data are combined
#   from separate BAM files, read groups provide identification of the
#   source of each read. Variant callers such as GATK depend on
#   having well-defined read group information.
# Here we set the RG variable to be the read group line we want
#   inserted in the header.
READ_GRP=$RGROUP;
RG='@RG\tID:1\tPL:ILLUMINA\tSM:'$READ_GRP'\tDS:ref='$ASSEMBLY',pfx='$REF_PFX

# Display how the program will be run, including
#   defaulted arguments. Do this before running
#   checks so user can see what went wrong.
echo "=================================================================";
echo "$0 - `date`";
if [ "$PAIRED" == "1" ]; then
echo "  fastq read1 file:  $IN_FQ_R1";
echo "  fastq read2 file:  $IN_FQ_R2";
else
echo "  input file:        $IN_FQ";
fi
echo "  output prefix:     $OUT_PFX";
echo "  reference:          $ASSEMBLY";
echo "  bwa version:       $BWA_VER";
echo "  ref prefix:        $REF_PFX";
echo "  read group line:   $RG"
echo "  will sort and fixmate BAM?:    $SORT_BAM";
echo "---------------------------------------------------------";

echo "";
echo "`date`: Sleeping for 10 seconds so that you can check if the above are correct."
echo "";
sleep 10;

# ------------------
# Error Checks
# ------------------
# Make sure the fastq file(s) exist.
# For paired end data, also make sure we have
#   two different files.
if [ "$PAIRED" == "1" ]; then
    ckFile "$IN_FQ_R1" "Fastq read1";
    ckFile "$IN_FQ_R2" "Fastq read2";
    if [ "$IN_FQ_R1" == "$IN_FQ_R2" ]; then
        err "Fastq read1 and read2 files are the same: '$IN_FQ_R1' '$IN_FQ_R2'";
    fi
else
    ckFile "$IN_FQ" "Input fastq";
fi

# Make sure we have found an appropriate reference
#   by checking that one of the standard files exists.
# ckFile "${REF_PFX}.amb" "$ASSEMBLY Reference";

# Make sure version information for our two programs is
#   part of our execution record. This is done by
#   calling the programs with no arguments
echo "---------------------------------------------------------";
echo "Program version information";
echo "---------------------------------------------------------";
$BWA_DIR_X/bwa
$SAMTOOLS_DIR_X/samtools

# echo "Exiting for dry run. Remove when satisfactory.";
# exit 1;

# ------------------
# The actual work!
 # a.llave@irri.org : We utilize multithreading here. However, we limit to 10 CPUs because there are reports that problematic SSAM output has been observed if more than 10 CPUs are used.
# ------------------

SAI_READ1="";
SAI_READ2="";
if [ "$PAIRED" == "0" ]; then
    echo "---------------------------------------------------------";
    echo "Running bwa mem on reads (single-end)";
    echo "Assembly: $ASSEMBLY"
    echo "FASTQ file: $IN_FQ";
    echo "---------------------------------------------------------";
    # at this point, $IN_FQ is just equal to $IN_FQ_R1
    $BWA_DIR_X/bwa mem -R $RG -t 10 $ASSEMBLY $IN_FQ > $POST_BWA_SAM_OUTPUT_FILENAME ; 
else
    echo "---------------------------------------------------------";
    echo "Running bwa mem on read1, read2 ends";
    echo "---------------------------------------------------------";
    $BWA_DIR_X/bwa mem -R $RG -t 10 $ASSEMBLY $IN_FQ_R1 $IN_FQ_R2 > $POST_BWA_SAM_OUTPUT_FILENAME ;
fi

wait;
echo "---------------------------------------------------------";
echo " `date` :  Checking if bwa mem successful...";
echo "---------------------------------------------------------";
ckFileSz $POST_BWA_SAM_OUTPUT_FILENAME;

echo "---------------------------------------------------------";
echo "`date`: Converting to BAM ...";
echo "---------------------------------------------------------";
$SAMTOOLS_DIR_X/samtools view -bS $POST_BWA_SAM_OUTPUT_FILENAME -o $OUT_PFX.bam


echo "---------------------------------------------------------";
echo "`date`: BAM conversion finished. Compressing SAM on a child process but will not wait for it";
echo "---------------------------------------------------------";
tar -czvpf $POST_BWA_SAM_OUTPUT_FILENAME.tar.gz $POST_BWA_SAM_OUTPUT_FILENAME &
COMPRESS_SAM_PID=$!

#################
#
# a.llave@irri.org REMARK 29OCT2012-0944 Sorting and indexing the BAM is only when it was
# 	intended when the script was run.
#
#################
# 	<area id="sort_and_index_bam" >
if [ "$SORT_BAM" == "1"  ]; then 
 echo "---------------------------------------------------------";
 echo "Creating sorted, fixmate-d, indexed bam file";
 echo "---------------------------------------------------------";
 echo " `date` :  Sorting '$OUT_PFX.bam'...";
 PRE_FIXMATE_NAME="$OUT_PFX.pre-fixmate.sorted"
 $SAMTOOLS_DIR_X/samtools sort $OUT_PFX.bam $PRE_FIXMATE_NAME;
 ckRes $? "$SAMTOOLS_DIR_X/samtools sort";
 ckFileSz "$PRE_FIXMATE_NAME.bam";

 # to save space
 echo "Deleting unsorted BAM...";
 rm $OUT_PFX.bam
 if [ ! -f $OUT_PFX.bam  ]; then
        echo "Unsorted BAM deleted.";
 else
        echo "Unsorted BAM still exists???";
 fi
 
 echo "`date` : Indexing '$PRE_FIXMATE_NAME.bam'...";
 $SAMTOOLS_DIR_X/samtools index $PRE_FIXMATE_NAME.bam;
 ckRes $? "$SAMTOOLS_DIR_X/samtools index";
 ckFileSz "$PRE_FIXMATE_NAME.bam.bai";

 echo "---------------------------------------------------------";
 echo "`date`: Cleaning up read pairing information and flags";
 echo "---------------------------------------------------------";
 $SAMTOOLS_DIR_X/samtools fixmate $PRE_FIXMATE_NAME.bam $OUT_PFX.sorted.bam
 wait;

 echo "`date` : Indexing '$OUT_PFX.sorted.bam'...";
 $SAMTOOLS_DIR_X/samtools index $OUT_PFX.sorted.bam;
 ckRes $? "$SAMTOOLS_DIR_X/samtools index";
 ckFileSz "$OUT_PFX.sorted.bam.bai";

 echo "---------------------------------------------------------";
 echo "Collecting alignment statistics";
 echo "---------------------------------------------------------";
 echo "`date` : Running flagstat...";
 $SAMTOOLS_DIR_X/samtools flagstat $OUT_PFX.sorted.bam | tee $OUT_PFX.flagstat.txt
 ckRes $? "$SAMTOOLS_DIR_X/samtools flagstat";
 ckFileSz "$OUT_PFX.flagstat.txt";
# 	</area>
fi

# If we make it here, all went well. Exit with a standard
#   message that can be easily grep'd
echo "---------------------------------------------------------";
echo "All bwa mem tasks completed successfully!";
echo "`date` : PID : $$";
echo ""
echo "If PID $COMPRESS_SAM_PID for compressing SAM file is still running"
echo "please wait for it to finish, otherwise, you can now delete SAM file manually:" 
echo "<<<------IS THERE A PROCESS BELOW HERE?------>>>";
ps aux | grep -v 'grep' | grep $COMPRESS_SAM_PID
echo "<<<------IS THERE A PROCESS ABOVE HERE?-------->>>";
echo "---------------------------------------------------------";
exit 0 ;

#!/bin/sh

echo
echo This is how MIRA extracted 'large contigs' from the total assembly for you:
echo  it made a list of contigs \(see info file ../output_d_info/output_info_largecontigs.txt\)
echo  which reached a certain length \(usually 500, see -MI:lcs=...\) and had at least 1/3 of
echo  the average coverage \(per sequencing technology\).
echo
echo Note that you can redefine what large contigs are for you by simply using any
echo combination of -n, -x, -y and -z parameters of 'miraconvert' instead of only the
echo '-n' parameter as used in this example.
echo
echo You can follow the progress of the conversion in the file "ec.log"
echo

miraconvert  -t caf -t maf -t tcs -t wig -t fasta -n ../output_d_info/output_info_largecontigs.txt -A "--job=genome,denovo,draft,Solexa -NW:cmrnl=no" output_out.maf output_LargeContigs_out >ec.log 2>&1

if [ $? -eq 0 ];then
   rm ec.log
   echo Finished, all done.
else
   tail -50 ec.log
   echo
   echo Ooops, something went wrong. Please consult the file 'ec-log' in this directory.
fi
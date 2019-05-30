#!/usr/bin/env python
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

glimmer_file = sys.argv[1]
fasta_file = sys.argv[2]
output_file = sys.argv[3]

# Read the sequence file
seq_record = SeqIO.parse(open(fasta_file),"fasta").next()

outseqrecords = []
# Read the glimmer file, record by record
for inline in file(glimmer_file):
    if '>' in inline:
        seqname = inline.split()[0][1:]
        outfilename = "%s_g3.tfa" % (seqname)
        continue
    if "orf" not in inline:
        continue
    orfname, sbegin, send, rf, score = inline.strip().split()
    sbegin = int(sbegin)
    send = int(send)
    rf = int(rf)
    # reverse complement
    if sbegin > send:
        sbegin, send = send, sbegin
    sbegin -= 1     # Python indexes start a 0
    score = float(score)
    # split the sequence record
    newseq = seq_record.seq[sbegin:send]
    if rf < 0:
        newseq = newseq.reverse_complement()
    # Add a sequence record to the output
    seqrecord_description = "begin=%d end=%d rf=%d score=%.2f" % (sbegin+1, send, rf, score)
    outseqrecords.append(SeqRecord(newseq,id=seqname+"_"+orfname, description=seqrecord_description))

SeqIO.write(outseqrecords,open(output_file,"w"),"fasta")
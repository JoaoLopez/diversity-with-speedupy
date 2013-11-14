CentC_Analyses
==============
**************************************
*      CentC Project Workflow        *
**************************************
By: Paul Bilinski, UC Davis, 2013

Original, large data files can be found on:
http://figshare.com/account/projects/124


I. Extract all CentC's from the maize reference genome (version 2)
-The original CentC sequences (Genbank + Nagaki et al 2003 Study) are in 

CentC_Seq_Originals.fasta

-Use the reverse complimented sequences for BLASTing.

-Follow instructions in the README_blast_anno.txt

-The DNA distance matrices can be generated using PHYLIP software package, DNADIST
function.  This project was run using PHYLIP 3.69.  To regenerate my dna distance 
matrices, I first performed 7 muscle alignments of the sequences I wanted to have in the
distance matrix (for example, all CentC from chromosome 2 or 5).  This was performed in
Geneious's muscle.  Then the alignment was fed to the dnadist.

-Execute the Spagedi analyses given the distance matrices and text files in the spagedi
manual.  These analyses were run with SPAGeDi 1.3a.

-
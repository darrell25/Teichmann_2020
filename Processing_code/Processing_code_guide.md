Summary of the initial sequence processing code used in this project

Step 1: De-multiplexed fastq files were run through trimmomatic to remove primer sequences. Used the following code:
Need to Input
Output: fastq files with primer sequences removed

Step 2: mothur was run using the "mothur_script.txt" file above in batch mode
Output: shared file, taxonomy file, representative sequence file

Step 3: The function "MakeBlastFasta" from InitialProcessingFunctions.R was run on the representative fasta file ouput by mothur
Output: FastaForBlast.fasta

Step 4: BLAST+ was downloaded and run with the following code:
blastn -db 16SMicrobial -query FastaForBlast.fasta -out BLASTtaxID_OTU.txt -outfmt "6 qseqid pident staxids sscinames" -max_target_seqs 1 -max_hsps 1
Output:BLASTtaxID_OTU.txt

Step 5: The function "AddSpeciesToTaxonomy" from InitialProcessingFunctions.R was run on the taxonomy file output from mothur and the output from BLAST+
Output: FF2_RDP_Blast.taxonomy

Step 6: The function "TransposeTable" from InitialProcessingFunctions.r was run on the shared OTU table output from mothur
Output: FF2_mothur_Transposed.shared

Step 7: a tree file was generated with the program FastTree using the following code:
FastTree -gtr -nt FastaforBlast.fasta > FF2.tree 
Output: FF2.tree

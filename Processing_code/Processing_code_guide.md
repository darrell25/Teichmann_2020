Summary of the initial sequence processing code used in this project

Step 1: De-multiplexed fastq files were run through cutadapt to remove primer sequences. Used the following code:<br/>
for %%G in ("D:\Sequencing\FF2\*.fastq") DO cutadapt -e 0.2 -g ^GTGYCAGCMGCCGCGGTAA -g ^GGACTACNVGGGTWTCTAAT -o D:\Sequencing\FF2\Trimmed\%%~nxG %%~nxG<br/>
Output: fastq files with primer sequences removed

Step 2: mothur was run using the "mothur_script.txt" file above in batch mode<br/>
Output: shared file, taxonomy file, representative sequence file

Step 3: The function "MakeBlastFasta" from InitialProcessingFunctions.R was run on the representative fasta file ouput by mothur<br/>
Output: FastaForBlast.fasta

Step 4: BLAST+ was downloaded and run with the following code:<br/>
blastn -db 16SMicrobial -query FastaForBlast.fasta -out BLASTtaxID_OTU.txt -outfmt "6 qseqid pident staxids sscinames" -max_target_seqs 1 -max_hsps 1<br/>
Output:BLASTtaxID_OTU.txt

Step 5: The function "AddSpeciesToTaxonomy" from InitialProcessingFunctions.R was run on the taxonomy file output from mothur and the output from BLAST+<br/>
Output: FF2_RDP_Blast.taxonomy

Step 6: The function "TransposeTable" from InitialProcessingFunctions.r was run on the shared OTU table output from mothur<br/>
Output: FF2_mothur_Transposed.shared

Step 7: a tree file was generated with the program FastTree using the following code:<br/>
FastTree -gtr -nt FastaforBlast.fasta > FF2.tree<br/> 
Output: FF2.tree

# 5'-end Seq Primary and secondary TSS calling

TSSpredator puts out 3 main factors to calculate transcription start sites; Step Height, Step Factor, Enrichment Factor. 

This code annotates genes for each coordinate given by TSSpredator. 150 nt before the start of every gene are assigned as 5' UTRs. 

Similarly, 150 nt downstream of the end of every gene was assigned as 3'UTR. 

This code also calculates the primary reads according to the maximum step Height by grouping genes and choosing the max step Height. 

To annotate secondary TSS, it calls the second highest maximum Step Height. 

# RA-Eval
DGE Analysis on Ewing sarcoma cells

I was tasked with performing DGE analysis on Ewing sarcoma cells using DESeq2. 
I first explored the data then started to create DEseq2 object with the design formula focusing on condition. 
Then I decied to filter the DESeq2 object by removing the the rows with low gene counts and setting the condition as a factor. 
Then I ran DESeq and created a table from the results of the differently expressed genes. 
I then made a new filtered table from the results by removing the NA values and only showing the genes with a p-value <0.05. 
To visualize the significant results I did a MA plot, PCA plot, and a Dispersion plot. 
The MA plot displays the log fold change vs the average normalized counts with most points residing on the y intercept of 0. 
The PCA plot shows the two separate conditions as two identifiable clusters. 
The Dispersion plot measures variablility in the data and showed that the genes with low read counts show higher dispersion than highly expressed genes. 
To prepare for the creating the heat map I converted the significant results table to a data frame then got the normalized counts and z-score for each row. 
The Heatmap shows the expression changes of the top 10 over and under expressed genes. 
The Volcano Plot shows the most significnat genes towards the top corners. 
Finally I ran KEGG Enrichment Analysis. First I made sure to pull from the KEGG pathway database. 
Then defined the significant genes based on p-value and set the cut off at p-value <= 0.01. 
Then I ran KEGG Enrichment analysis and decided the visualize the results using a barplot. 
Which showed the KEGG pathways in cancer is the most altered in EWSR1-FLI1 knock-down (shEF1) vs control (shCTR).

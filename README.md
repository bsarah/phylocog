# phylocog
Pipeline to create sequence and structure alignments for protein sequences and corresponding phylogenetic trees

Given: protein sequences of orthologous proteins (COGs, https://www.ncbi.nlm.nih.gov/research/cog-project/)

Steps:

- (S) sequence based:

-> (S1) create MSAs for each COG (MAFFT L-INS-I, https://mafft.cbrc.jp/alignment/software/)

-> (S2) create rooted phylogenetic trees (FastTree, http://meta.microbesonline.org/fasttree/)

-> (S3/D6) calculate splits and D-values (modified Fitch:phylocog/treeAnalysis)

-> (S4/D7) plot trees (ete3, phylocog/treeAnalysis)



- (D) structure domain based:

-> (D1) scan protein sequences for ECOD (http://prodata.swmed.edu/ecod/) domains (hmmscan, hmmer.org)

-> (D2) map structural domains to sequences (DomainMapper, https://onlinelibrary.wiley.com/doi/pdf/10.1002/pro.4465)

-> (D3) parse DomainMapper output (phylocog/dmParsing)

-> (D3b) plot structural domain composition (phylocog/plotting)

-> (D4) pairwise alignment (dNWA, see https://github.com/blaueste/bachelor_thesis)

-> (D5) construction of phylogenetic trees (Feng-Doolittle: phylocog/aln2tree)

-> (D6/S3) calculation of splits and D-values (modified Fitch: phylocog/treeAnalysis)

-> (D7/S4) plot trees (ete3, phylocog/treeAnalysis)


- (C) comparison of S-trees and D-trees
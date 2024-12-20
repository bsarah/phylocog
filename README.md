# phylocog
Pipeline to create sequence and structure alignments for protein sequences and corresponding phylogenetic trees

Given: protein sequences of orthologous proteins (COGs, https://www.ncbi.nlm.nih.gov/research/cog-project/)

Steps:

- (S) sequence based:

-> (S1) create MSAs for each COG (MAFFT L-INS-I, https://mafft.cbrc.jp/alignment/software/)

-> (S2) create rooted phylogenetic trees (FastTree, http://meta.microbesonline.org/fasttree/)

-> (S2b) reformat tree to have protIDs and taxIDs on each node or even pattern IDs from structural domains (phylocog/formatting)

-> (S3/D6) calculate splits and D-values (modified Fitch:phylocog/treeAnalysis)

-> (S4/D7 a) plot trees colored based on archaeal and bacterial proteins (ete3, phylocog/treeAnalysis)
-> (S4/D7 b) plot trees colored based on structural patterns (ete3, phylocog/plotting)



- (D) structure domain based:

-> (D1) scan protein sequences for ECOD (http://prodata.swmed.edu/ecod/) domains (hmmscan, hmmer.org)

-> (D2) map structural domains to sequences (DomainMapper, https://onlinelibrary.wiley.com/doi/pdf/10.1002/pro.4465)

-> (D3) parse DomainMapper output (phylocog/dmParsing)
ls

-> (D3b) plot structural domain composition (phylocog/plotting)

-> (D4) pairwise alignment (dNWA, see https://github.com/blaueste/bachelor_thesis)

-> (D5) construction of phylogenetic trees (Feng-Doolittle (e.g. https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle): phylocog/aln2tree)

-> (D6/S3) calculation of splits and D-values (modified Fitch: phylocog/treeAnalysis)

-> (D7/S4 a) plot trees colored based on archaeal and bacterial proteins (ete3, phylocog/treeAnalysis)
-> (D7/S4 b) plot trees colored based on structural patterns (ete3, phylocog/plotting)


- (C) comparison of S-trees and D-trees
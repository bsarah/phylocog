# phylocog
Pipeline to create sequence and structure alignments for protein sequences and corresponding phylogenetic trees

Given:
- protein sequences of orthologous proteins (COGs)

Steps:

- sequence based:
-- S1 create MSAs for each COG (MAFFT L-INS-I)
-- S2 create rooted phylogenetic trees (FastTree)
-- S3 calculate splits and D-values (modified Fitch)
-- S4 plot trees (ete3)

- structure domain based:
-- D1 scan protein sequences for ECOD domains (hmmscan)
-- D2 map structural domains to sequences (DomainMapper)
-- D3 parse DomainMapper output (dommap-parsing)
-- D4 pairwise alignment (dNWA)
-- D5 construction of phylogenetic trees (Feng-Doolittle)
-- D6 calculation of splits and D-values (modified Fitch)
-- D7 plot trees (ete3)


- comparison of S-trees and D-trees
# Amphioxus Conserved Gene Cluster Identification Pipeline
A bioinformatics pipeline to identify conserved gene clusters across multiple amphioxus species (Alref, Bb, Bf, Bj, Bl) using collinearity data from MCScanX results.
## This pipeline consists of three core steps to detect cross-species conserved continuous gene clusters in amphioxus genomes:
step0.Extract conserved collinear gene pairs across target species from MCScanX collinearity files

step1.Identify intra-species continuous gene clusters for each amphioxus species

step2.Integrate intra-species clusters to detect cross-species conserved major clusters

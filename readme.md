This repo contains analyses for processing nanopore sequencing data of anabaena strains sequencing in the Azolla lab at Utrecht University.
Specifically, we're looking for large in/del variants in sequenced strains created by RNA guided transposition.

The _anabaena_/_nostoc_ reference strain: _Nostoc spec_ pcc7120 
Downloaded from: https://www.ncbi.nlm.nih.gov/assembly/GCF_000009705.1/


The analyses documented here include two main approaches. 
First a denovo approach, assembling the anabaena genomes one by one,
and second a reference based approach.

The denovo approach includes:
1. de-novo assembly with flye (dir `denovo`)
2. assembly polishing with medaka (dir `denovo/sample/polished-medaka`)

The reference based approach includes
1. mapping to the ncbi reference with minimap2, then variant calling with medaka
2. mapping to the ncbi reference with ngmlr, then calling structural variants with sniffles

Finally, these last two variant calling steps are also performed against the wild type de-novo sequenced strain.









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
3. assembly annotation with both prokka and bakta
4. mapping of sample and reference strain reads with minimap2
5. locating regions of interest with blat
6. visualisation of all generated data with igv

The reference based approach can take various reference-sample combinations and do:
1. mapping to a reference with minimap2, then variant calling with medaka
2. mapping to a reference with ngmlr, then calling structural variants with sniffles
3. extract fasta files with insertion and deletion SVs
4. locating regions of interest with blat
5. visualisation of all generated data with igv

A mapping table should match up sample names with their appropriate reference like so
```
#ref	samples
WT1	sampleA sampleB sampleC
WT2	sampleX sampleY sampleZ
```



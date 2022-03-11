#!/bin/bash
# This script takes nanopore FastQ data and does several types of variant calling.
# Variant calling is done against a downloaded (ncbi) reference, and a denovo assembled reference.
# All software is supposed to be in dedicated conda environments.
# First, this script takes the raw data and does flye assemblies for all samples.
# When assembled, a snapshot of the assembly graph is generated with Bandage.
# Second, it polishes these flye assemblies with medaka.
# Third, it annotates these pollished assemblies with both prokka and bakta for comparison
# Finally, this script does two forms of variant calling to both the denovo assembled WT and a reference genome.
# The wild type assembly must be defined up front and a ncbi reference must be defined as well.
# Variants are called with medaka again, which includes minimap2 readmapping.
# Structural variants are called with Sniffles.
# To this end, nanopore reads are mapped with a dedicated aligner: ngmlr.

# For visualisation of denovo assemlbes in IGV, I recommend loading a polished genome and the appropriate gff files from bakta and prokka.
# Except for the WT, no read mappings are available.
# For visualisation of variant calling on the denovo WT assembly, I recommend loading the polished WT genome and appropriate GFF files.
# Read mappings are available as minimap2 mappings by medaka and as ngmlr mappings as well. You can load both or just one of the two.
# For visualisation of variant calling on the reference assembly, load the fasta and gff from the reference directory.
# Load bam and vcf files similarly as for the denovo WT assembly.


# define where stuff is:
basedir=/stor/anabaena                                             # this is where we create our output
fqdir=/stor/azolla_sequencing/nanopore/anabaena                    # a dir with nanopore .fastq.gz files to process
ncbi="$basedir"/reference/GCF_000009705.1_ASM970v1_genomic.fna     # a reference genome to call variants on
baktaDB=/stor/scripts/baktaDB                                      # bakta db for annotation
WT_name=WT                                                         # the name of your wild type sample to use as a reference for variant calling
WT="$basedir"/denovo/"$WT_name"/polished-medaka/consensus.fasta    # auto genenerated based on the line above

# A trick to swtich conda environments while using this script, adapt to your particular conda installation.
source /home/laura/miniconda3/etc/profile.d/conda.sh

# get an array of our samples directly from the available sequencing files and check if any samples are found
samples=( $(find "$fqdir" -maxdepth 1 -name '*.fastq.gz' -printf '%P\n') )
if   [ "${#samples[@]}" -lt 1 ]
then echo no samples found
     exit
fi

# for each sample, make a de novo assembly with flye
echo 'Checking if all denovo assemblies are present'
conda activate flye
if     [ ! -d "$basedir"/denovo ]
then   mkdir "$basedir"/denovo
fi
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -d "$basedir"/denovo/"$name" ]
      then   if   [ ! $(command -v flye) ]
             then echo 'flye is not found'
                  exit
             fi
             # assemble with flye expecting a genome of 6.4Mb
             flye --nano-hq "$fqdir/$s"    \
                  --genome-size 6.4M       \
                  --threads $(nproc)       \
                  --scaffold               \
                  --out-dir "$basedir"/denovo/"$name"
             Bandage image "$basedir"/denovo/"$name"/assembly_graph.gfa \
                           "$basedir"/denovo/"$name".png
      fi
done
conda deactivate

# for each sample, polish the assembly with medaka
echo 'Checking if all denovo assemblies are polished'
conda activate medaka
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -d "$basedir"/denovo/"$name/polished-medaka" ]
      then   if   [ ! $(command -v medaka_consensus) ]
             then echo 'medaka is not found'
                  exit
             fi
             medaka_consensus -i "$fqdir/$s"    \
                              -d "$basedir"/denovo/"$name"/assembly.fasta  \
                              -o "$basedir"/denovo/"$name/polished-medaka" \
                              -m r941_min_sup_variant_g507                 \
                              -t 6
      fi
done
conda deactivate

# Annotate all assembled and polished genomes with prokka
echo 'Checking if all polished assemblies are annotated with prokka'
conda activate prokka
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      wd="$basedir"/denovo/"$name"/polished-medaka_prokka-annotation
      if     [ ! -d "$wd" ]
      then   if   [ ! $(command -v prokka) ]
             then echo 'prokka is not found'
                  exit
             fi
             prokka --outdir "$wd"    \
                    --addgenes        \
                    --genus 'nostoc'  \
                    --prefix "$name"  \
                    --kingdom 'Bacteria' \
                    --cpus 0          \
                    --rfam            \
                    "$basedir"/denovo/"$name/polished-medaka/consensus.fasta"
      fi
done
conda deactivate

# Annotate all assembled and polished genomes with bakta
if   [ ! -d "$baktaDB"/amrfinderplus-db ]
then echo 'amrfinderplus-db is not setup correctly, doing that now'
     amrdinamrfinder_update --database "$baktaDB"/amrfinderplus-db
fi

echo 'Checking if all polished assemblies are annotated with bakta'
conda activate bakta
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      wd="$basedir"/denovo/"$name"/polished-medaka_bakta-annotation
      if     [ ! -d "$wd" ]
      then   if   [ ! $(command -v bakta) ]
             then echo 'prokka is not found'
                  exit
             fi
             bakta --output "$wd"    \
                   --db "$baktaDB"   \
                   --genus 'nostoc'  \
                   --prefix "$name"  \
                   --complete        \
                   --threads $(nproc)\
                    "$basedir"/denovo/"$name/polished-medaka/consensus.fasta"
      fi
done
conda deactivate

# call variants with medaka on the wild type de novo assembly
echo 'Checking if all samples are used for variant calling on the denovo assembled WT'
conda activate medaka
wd="$basedir"/haplotypes_WT/medaka
if     [ ! -d   "$wd" ]
then   mkdir -p "$wd"
fi

for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -d "$wd/$name" ]
      then   if   [ ! $(command -v medaka_haploid_variant) ]
             then echo 'cant find medaka'
             fi
             medaka_haploid_variant -i "$fqdir/$s"    \
                                    -r "$WT"          \
                                    -o "$wd/$name"    \
                                    -t $(nproc)       \
                                    -m r941_min_sup_variant_g507
      fi
done
conda deactivate

# map reads with ngmlr on the WT assembled and polished genome for variant calling with sniffles
echo 'Checking if all samples are mapped to the denovo WT with ngmlr'
conda activate nanopore
wd="$basedir"/haplotypes_WT/mapped_ngmlr
if     [ ! -d   "$wd" ]
then   mkdir -p "$wd"
fi

for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -f "$wd/$name".sorted.bam ]
      then   if   [ ! $(command -v ngmlr) ]
             then echo 'cant find ngmlr'
             fi
             ngmlr -q "$fqdir/$s"     \
                   -r "$WT"           \
                   --rg-sm "$name"    \
                   -t $(nproc)        \
                   -x ont             \
                   -o "$wd/$name".sam
             # now process the bam file with samtools for later use and visualisation
             samtools sort -@ 6 -m 9G "$wd/$name".sam \
             | samtools view -b -@ 6 -h               \
             > "$wd/$name".sorted.bam
             samtools index "$wd/$name".sorted.bam
             rm "$wd/$name".sam
      fi
done

# call variants with sniffles on the WT assembled and polished genome
echo 'Checking if all samples are used for structural variant calling on the denovo WT'
wd="$basedir"/haplotypes_WT/sniffles
if     [ ! -d   "$wd" ]
then   mkdir -p "$wd"
fi

for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -f "$wd/$name".vcf ]
      then   if   [ ! $(command -v sniffles) ]
             then echo 'cant find sniffles'
             fi
             sniffles --input "$wd"/../mapped_ngmlr/"$name".sorted.bam  \
                      --reference "$WT"                                 \
                      --snf "$wd/$name".snf                             \
                      --vcf "$wd/$name".vcf
      fi
done
# now combine all sniffles calls in one vcf
if   [ ! -f "$wd"/multi-sample.vcf ]
then sniffles --input "$wd/*.snf" --vcf "$wd"/multi-sample.vcf
fi
conda deactivate

# call variants with medaka on the ncbi reference
echo 'Checking if all samples are used for variant calling on the reference'
conda activate medaka
wd="$basedir"/haplotypes_ref/medaka
if     [ ! -d   "$wd" ]
then   mkdir -p "$wd"
fi

for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastqcc\.gz//g' )
      if     [ ! -d "$wd/$name" ]
      then   if   [ ! $(command -v medaka_haploid_variant) ]
             then echo 'cant find medaka'
             fi
             medaka_haploid_variant -i "$fqdir/$s"    \
                                    -r "$ncbi"        \
                                    -o "$wd/$name"    \
                                    -t $(nproc)       \
                                    -m r941_min_sup_variant_g507
      fi
done
conda deactivate


# map reads with ngmlr on the ncbi reference (for variant calling with sniffles)
echo 'Checking if all samples are mapped to the reference with ngmlr'
conda activate nanopore
wd="$basedir"/haplotypes_ref/mapped_ngmlr
if     [ ! -d   "$wd" ]
then   mkdir -p "$wd"
fi

for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -f "$wd/$name".sorted.bam ]
      then   if   [ ! $(command -v ngmlr) ]
             then echo 'cant find ngmlr'
             fi
             ngmlr -q "$fqdir/$s"    \
                   -r "$ncbi"        \
                   --rg-sm "$name"   \
                   -t $(nproc)       \
                   -x ont            \
                   -o "$wd/$name".sam
             samtools sort -@ 6 -m 9G "$wd/$name".sam \
             | samtools view -b -@ 6 -h               \
             > "$wd/$name".sorted.bam
             samtools index "$wd/$name".sorted.bam
             rm "$wd/$name".sam
      fi
done

# call variants with sniffles on the ncbi reference
echo 'Checking if all samples are used for structural variant calling on the reference'
wd="$basedir"/haplotypes_ref/sniffles
if     [ ! -d   "$wd" ]
then   mkdir -p "$wd"
fi

for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -f "$wd/$name".vcf ]
      then   if   [ ! $(command -v sniffles) ]
             then echo 'cant find sniffles'
             fi
             sniffles --input "$wd"/../mapped_ngmlr/"$name".sorted.bam     \
                      --reference "$ncbi"                           \
                      --snf "$wd/$name".snf                         \
                      --vcf "$wd/$name".vcf
      fi
done
# now combine all sniffles calls in one vcf
if   [ ! -f "$wd"/multi-sample.vcf ]
then sniffles --input "$wd/*.snf" --vcf "$wd"/multi-sample.vcf
fi
conda deactivate
echo 'Script finished'

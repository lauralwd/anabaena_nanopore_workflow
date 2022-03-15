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
condadir=/home/laura/miniconda3                                    # (mini)conda(3) directory
WT="$basedir"/denovo/WT/polished-medaka/consensus.fasta            # the name(s) of your wild type sample to use as a reference for variant calling
CSV15="$basedir"/denovo/CSV15/polished-medaka/consensus.fasta
refs=( "$WT" "$CSV15" "$ncbi" )
ref_names=( 'WT' 'CSV15' 'ncbi' )
maptab="$basedir"/WT_sample.txt

# A trick to swtich conda environments while using this script, adapt to your particular conda installation.
if   [ ! -f "$condadir"/etc/profile.d/conda.sh ]
then echo 'quiting for we need the conda environments in the `envs` directory to proceed'
     exit
else source "$condadir"/etc/profile.d/conda.sh
fi

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

# for each sample, map the sample reads to the polished assembly
echo 'Checking if all denovo assemblies have reads mapped back'
conda activate nanopore
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -f "$basedir"/denovo/"$name/polished-medaka/$name.bam" ]
      then   if   [ ! $(command -v minimap2) ]
             then echo 'minimap2 is not found'
                  exit
             fi
             if   [ ! -f "$basedir"/denovo/"$name"/assembly.fasta.mmi ]
             then minimap2    "$basedir"/denovo/"$name"/assembly.fasta \
                           -d "$basedir"/denovo/"$name"/assembly.fasta.mmi
             fi
             minimap2 -d "$basedir"/denovo/"$name"/assembly.fasta \
                      "$fqdir/$s"           \
                      -x map-ont            \
                      -t 6                  \
                      -Y                    \
                      -a                    \
            | samtools sort -@ 6 -l 9 -m 9G \
            | samtools view -b              \
            > "$basedir"/denovo/"$name/polished-medaka/$name.bam"
      fi
done
conda deactivate
exit

# for each sample, map the WT reads to the polished assembly
echo 'Checking if all denovo assemblies have reads mapped back'
conda activate nanopore
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      refname=$(grep "$name" "$maptab" | grep -v ncbi | cut -f 1 )
      if     [ ! -d "$basedir"/denovo/"$name/polished-medaka/$refname.bam" ]
      then   if   [ ! $(command -v minimap2) ]
             then echo 'minimap2 is not found'
                  exit
             fi
             if   [ ! -f "$basedir"/denovo/"$name"/assembly.fasta.mmi ]
             then minimap2    "$basedir"/denovo/"$name"/assembly.fasta \
                           -d "$basedir"/denovo/"$name"/assembly.fasta.mmi
             fi
             minimap2 -d "$basedir"/denovo/"$name"/assembly.fasta \
                      "$fqdir/$refname".fastq.gz  \
                      -x map-ont            \
                      -t 6                  \
                      -Y                    \
                      -a                    \
            | samtools sort -@ 6 -l 9 -m 9G \
            | samtools view -b              \
            > "$basedir"/denovo/"$name/polished-medaka/$refname.bam"
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
                   --compliant       \
                   --complete        \
                   --keep-contig-headers \
                   --threads $(nproc)\
                    "$basedir"/denovo/"$name/polished-medaka/consensus.fasta"
      fi
done
conda deactivate

# call variants with medaka on all selected reference sequences
conda activate medaka
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo "Checking if all samples are used for variant calling on reference $refname"

     wd="$basedir"/haplotypes_"$refname"/medaka
     if     [ ! -d   "$wd" ]
     then   mkdir -p "$wd"
     fi

     for   s in "${samples[@]}"
     do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
           if     [ ! -d "$wd/$name" ]
           then   if   [ ! $(command -v medaka_haploid_variant) ]
                  then echo 'cant find medaka'
                  elif [ $(grep "$refname" "$maptab" | cut -f 2 | grep "$name" -c ) -eq 1 ]
                  then medaka_haploid_variant -i "$fqdir/$s"    \
                                              -r "${refs[$count]}"    \
                                              -o "$wd/$name"    \
                                              -t $(nproc)       \
                                              -m r941_min_sup_variant_g507
                  fi
           fi
     done
done
conda deactivate
exit

# map reads with ngmlr for variant calling with sniffles
conda activate nanopore
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo "Checking if all samples are mapped to $refname with ngmlr"

     wd="$basedir"/haplotypes_"$refname"/mapped_ngmlr
     if     [ ! -d   "$wd" ]
     then   mkdir -p "$wd"
     fi

     for   s in "${samples[@]}"
     do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
           if     [ ! -f "$wd/$name".sorted.bam ]
           then   if   [ ! $(command -v ngmlr) ]
                  then echo 'cant find ngmlr'
                  elif [ $(grep "$refname" "$maptab" | cut -f 2 | grep "$name" -c ) -eq 1 ]
                  then ngmlr -q "$fqdir/$s"     \
                             -r "${refs[$count]}" \
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
           fi
      done
done

#call variants with sniffles on the reference of of choice
conda acticate sniffles
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo "Checking if all samples are used for structural variant calling on reference $refname"

     wd="$basedir"/haplotypes_"$refname"/sniffles
     if     [ ! -d   "$wd" ]
     then   mkdir -p "$wd"
     fi

     for   s in "${samples[@]}"
     do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
           if     [ ! -f "$wd/$name".vcf ]
           then   if   [ ! $(command -v sniffles) ]
                  then echo 'cant find sniffles'
                  elif [ $(grep "$refname" "$maptab" | cut -f 2 | grep "$name" -c ) -eq 1 ]
                  then sniffles --input "$wd"/../mapped_ngmlr/"$name".sorted.bam  \
                                --reference "${refs[$count]}"                     \
                                --snf "$wd/$name".snf                             \
                                --vcf "$wd/$name".vcf
                  fi

           fi
     done
     # now combine all sniffles calls in one vcf
     if   [ ! -f "$wd"/"$refname"_multi-sample.vcf ]
     then sniffles --input "$wd/"*.snf --vcf "$wd"/"$refname"_multi-sample.vcf
     fi
done
conda deactivate

conda activate igv
# searching for sequences of interest in all reference assemblies
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo "Using blat to annotate particular regions of interest on $refname"

     wd="$basedir"/igv_configs/
     if    [ ! -d   "$wd" ]
     then  mkdir -p "$wd"
     fi
     if    [ ! -f "$wd"/"$refname"_blathits.psl
     then  blat  "${refs[$count]}"              \
                 queries.fasta                  \
                 -t=DNA -q=DNA
                 "$wd"/"$refname"_blathits.psl
     fi
done
conda deactivate

echo 'Script finished'

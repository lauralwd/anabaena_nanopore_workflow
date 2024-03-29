#!/bin/bash
# This script takes nanopore FastQ data and does several types of variant calling.
# Variant calling is done against a downloaded (ncbi) reference, and denovo assembled references.
# A sample-mapping table corresponds the appropriate wild type sample with its treatments

# All software is supposed to be in dedicated conda environments.
# First, this script takes the raw data and does flye assemblies for all samples.
# When assembled, a snapshot of the assembly graph is generated with Bandage.
# Second, it polishes these flye assemblies with medaka.
# The original sample and corresponding WT reads are mapped to the polished assembly for later visualisation.
# Third, it annotates these pollished assemblies with both prokka and bakta for comparison
# Finally, this script does two forms of variant calling to all relevant references.
# The relevant references are defined in a mapping table, corresponding samples to their WT.
# Variants are called with medaka again, which includes minimap2 readmapping.
# Structural variants are called with Sniffles.
# To this end, nanopore reads are mapped with a dedicated aligner: ngmlr.

# For visualisation of denovo assemlbes in IGV, I recommend loading a polished genome and the appropriate gff files from bakta and prokka.
# BAM files for assembly verification are available next to the consensus fasta file.
# For visualisation of variant calling on any reference, I recommend loading the appropriate reference with its GFF files.
# Next, load the VCF files generated by medaka (one file per sample) and sniffles (one file per reference)
# Read mappings are available as minimap2 mappings by medaka and as ngmlr mappings as well. You can load both or just one of the two.

# define where stuff is:
basedir=/stor/anabaena                                             # this is where we create our output
fqdir=/stor/azolla_sequencing/nanopore/anabaena                    # a dir with nanopore .fastq.gz files to process
baktaDB=/stor/scripts/baktaDB                                      # bakta db for annotation
condadir=/home/laura/miniconda3                                    # (mini)conda(3) directory

# define references
ncbi="$basedir"/reference/GCF_000009705.1_ASM970v1_genomic.fna     # a reference genome to call variants on
WT="$basedir"/denovo/WT/polished-medaka/consensus.fasta            # the name(s) of your wild type sample to use as a reference for variant calling
CSV15="$basedir"/denovo/CSV15/polished-medaka/consensus.fasta
CSAM="$basedir"/denovo/CSAM/polished-medaka/consensus.fasta
refs=( "$WT" "$CSV15" "$ncbi" "$CSAM" )
ref_names=( 'WT' 'CSV15' 'ncbi' 'CSAM')
maptab="$basedir"/WT_sample.txt
CPU=$(nproc)

# Coloured text
BLD=$(tput bold)
GRN=$(tput setaf 2)
RED=$(tput setaf 1)
BLU=$(tput setaf 4)
NML=$(tput sgr0)

# check mapping table presence
if   [ ! -f "$maptab" ]
then echo -e "$BLD""$RED""Wildtype to sample mapping table is not found, make sure to set this propperly. $NML"
     echo -e "expecting mapping table at path $maptab"
     exit 1
fi

# functions to use in the script
function checkprog {
if   [ ! $(command -v "$1") ]
then echo -e "$BLD""$RED""Command $1 is not found $NML"
     echo 'Make sure the propper conda environments are installed'
     exit 1
fi
}

function checkwd {
if   [ ! -d "$wd" ]
then mkdir  "$wd"
fi
}

# A trick to swtich conda environments while using this script, adapt to your particular conda installation.
if   [ ! -f "$condadir"/etc/profile.d/conda.sh ]
then echo -e "$BLD""$RED""quiting for we need the conda environments in the `envs` directory to proceed $NML"
     exit
else source "$condadir"/etc/profile.d/conda.sh
fi

# get an array of our samples directly from the available sequencing files and check if any samples are found
samples=( $(find "$fqdir" -maxdepth 1 -name '*.fastq.gz' -printf '%P\n') )
if   [ "${#samples[@]}" -lt 1 ]
then echo -e "$BLD""$RED""no samples found $NML"
     exit
fi

# for each sample, make a de novo assembly with flye
echo -e "$BLD""$GRN""Checking if all denovo assemblies are present $NML"
conda activate flye
wd="$basedir"/denovo
checkwd
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -d "$wd/$name" ]
      then   checkprog flye
             checkprog Bandage
             # assemble with flye expecting a genome of 6.4Mb
             flye --nano-hq "$fqdir/$s"    \
                  --genome-size 6.4M       \
                  --threads "$CPU"         \
                  --scaffold               \
                  --out-dir "$wd/$name"
             Bandage image "$wd/$name"/assembly_graph.gfa \
                           "$wd/$name".png
      fi
done
conda deactivate

# for each sample, polish the assembly with medaka
echo -e "$BLD""$GRN""Checking if all denovo assemblies are polished $NML"
conda activate medaka
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -d "$basedir"/denovo/"$name/polished-medaka" ]
      then   checkprog medaka_consensus
             medaka_consensus -i "$fqdir/$s"    \
                              -d "$basedir"/denovo/"$name"/assembly.fasta  \
                              -o "$basedir"/denovo/"$name/polished-medaka" \
                              -m r941_min_sup_variant_g507                 \
                              -t 6
      fi
done
conda deactivate

# for each sample, map the sample reads to the polished assembly
echo -e "$BLD""$GRN""Checking if all denovo assemblies have sample reads mapped back $NML"
conda activate nanopore
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      if     [ ! -f "$basedir"/denovo/"$name/polished-medaka/$name.bam" ]
      then   checkprog minimap2
             checkprog samtools
             # check if a minimap2 index is already present:
             if   [ ! -f      "$basedir"/denovo/"$name"/polished-medaka/consensus.fasta.mmi ]
             then minimap2    "$basedir"/denovo/"$name"/polished-medaka/consensus.fasta \
                           -d "$basedir"/denovo/"$name"/polished-medaka/consensus.fasta.mmi
             fi
             # run minimap2, then sort the samfile and convert to bam
             minimap2 "$basedir"/denovo/"$name"/polished-medaka/consensus.fasta \
                      "$fqdir/$s"           \
                      -x map-ont            \
                      -t "$CPU"             \
                      -Y                    \
                      -a                    \
             | samtools sort -@ 6 -l 9 -m 9G \
             > "$basedir"/denovo/"$name/polished-medaka/$name.bam"
             # index the bamfile for igv
             samtools index "$basedir"/denovo/"$name/polished-medaka/$name.bam"
      fi
done
conda deactivate

# for each sample, map the WT reads to the polished assembly
echo -e "$BLD""$GRN""Checking if all denovo assemblies have WT reads mapped back $NML"
conda activate nanopore
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      refname=$(grep '\W'"$name" "$maptab" | grep -v ncbi | cut -f 1 )
      if     [ ! -f "$basedir"/denovo/"$name/polished-medaka/$refname.bam" ]
      then   checkprog minimap2
             checkprog samtools
             # check if a minimap2 index is already present:
             if   [ ! -f "$basedir"/denovo/"$name"/assembly.fasta.mmi ]
             then minimap2    "$basedir"/denovo/"$name"/assembly.fasta \
                           -d "$basedir"/denovo/"$name"/assembly.fasta.mmi
             fi
             # run minimap2, then sort the samfile and convert to bam
             minimap2 "$basedir"/denovo/"$name"/polished-medaka/consensus.fasta \
                      "$fqdir/$refname".fastq.gz  \
                      -x map-ont            \
                      -t 6                  \
                      -Y                    \
                      -a                    \
             | samtools sort -@ 6 -l 9 -m 9G \
             > "$basedir"/denovo/"$name/polished-medaka/$refname.bam"
             # index the bamfile for igv
             samtools index "$basedir"/denovo/"$name/polished-medaka/$refname.bam"
      fi
done
conda deactivate

# Annotate all assembled and polished genomes with prokka
echo -e "$BLD""$GRN""Checking if all polished assemblies are annotated with prokka $NML"
conda activate prokka
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      wd="$basedir"/denovo/"$name"/polished-medaka_prokka-annotation
      if     [ ! -d "$wd" ]
      then   checkprog prokka
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
conda activate bakta
if   [ ! -d "$baktaDB"/amrfinderplus-db ]
then echo -e "$BLD""$RED""amrfinderplus-db is not setup correctly, doing that now $NML"
     checkprog amrfinder_update
     amrfinder_update --database "$baktaDB"/amrfinderplus-db
fi

echo -e "$BLD""$GRN""Checking if all polished assemblies are annotated with bakta $nml"
for   s in "${samples[@]}"
do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
      wd="$basedir"/denovo/"$name"/polished-medaka_bakta-annotation
      if     [ ! -d "$wd" ]
      then   checkprog bakta
             bakta --output "$wd"    \
                   --db "$baktaDB"   \
                   --genus 'nostoc'  \
                   --prefix "$name"  \
                   --compliant       \
                   --complete        \
                   --keep-contig-headers \
                   --threads "$CPU"  \
                    "$basedir"/denovo/"$name/polished-medaka/consensus.fasta"
      fi
done
conda deactivate

# call variants with medaka on all selected reference sequences
conda activate medaka
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo -e "$BLD""$GRN""Checking if all samples are used for variant calling on reference $refname $NML"

     wd="$basedir"/haplotypes_"$refname"/medaka
     checkwd

     for   s in "${samples[@]}"
     do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
           if     [ ! -d "$wd/$name" ]
           then   checkprog medaka_haploid_variant
                  # check if this reference-sample combo should be ran,
                  # otherwise continue to the next combo
                  if [ $(grep \^"$refname" "$maptab" | grep "$name" -c ) -eq 1 ]
                  then medaka_haploid_variant -i "$fqdir/$s"    \
                                              -r "${refs[$count]}"    \
                                              -o "$wd/$name"    \
                                              -t "$CPU"         \
                                              -m r941_min_sup_variant_g507
                  fi
           fi
     done
done
conda deactivate

# map reads with ngmlr for variant calling with sniffles
conda activate nanopore
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo -e "$BLD""$GRN""Checking if all samples are mapped to $refname with ngmlr $NML"

     wd="$basedir"/haplotypes_"$refname"/mapped_ngmlr
     checkwd

     for   s in "${samples[@]}"
     do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
           if     [ ! -f "$wd/$name".sorted.bam ]
           then   checkprog ngmlr
                  checkprog samtools
                  # check if this reference-sample combo should be ran,
                  # otherwise continue to the next combo
                  if [ $(grep \^"$refname" "$maptab" | grep "$name" -c ) -eq 1 ]
                  then ngmlr -q "$fqdir/$s"     \
                             -r "${refs[$count]}" \
                             --rg-sm "$name"    \
                             -t "$CPU"          \
                             -x ont             \
                             -o "$wd/$name".sam
                  # now process the bam file with samtools for later use and visualisation
                  samtools sort -@ 6 -m 9G "$wd/$name".sam \
                  > "$wd/$name".sorted.bam
                  samtools index "$wd/$name".sorted.bam
                  rm "$wd/$name".sam
                  fi
           fi
      done
done

#call variants with sniffles on the reference of of choice
for  r in $(seq 1 1 "${#refs[@]}" )
do   count=$(echo "$r -1" | bc)      # correct for 0based counting
     refname="${ref_names[$count]}"  # define refname
     echo -e "$BLD""$GRN""Checking if all samples are used for structural variant calling on reference $refname $NML"

     wd="$basedir"/haplotypes_"$refname"/sniffles
     checkwd

     for   s in "${samples[@]}"
     do    name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
           if     [ ! -f "$wd/$name".vcf ]
           then   checkprog sniffles
                  # check if this reference-sample combo should be ran,
                  # otherwise continue to the next combo
                  if [ $(grep \^"$refname" "$maptab" | grep "$name" -c ) -eq 1 ]
                  then sniffles --input "$wd"/../mapped_ngmlr/"$name".sorted.bam  \
                                --reference "${refs[$count]}"                     \
                                --sample-id "$name"                               \
                                --snf "$wd/$name".snf                             \
                                --vcf "$wd/$name".vcf
                  fi

           fi
     done
     # now combine all sniffles calls in one vcf
     if   [ ! -f "$wd"/"$refname"_multi-sample.vcf ]
     then sniffles --input "$wd/"*.snf --vcf "$wd"/"$refname"_multi-sample.vcf
          # arbitrary fix to remove whole chromosome duplications
          grep -Pv '\tSniffles2\.DUP' "$wd"/"$refname"_multi-sample.vcf > ./temp.txt
          mv temp.txt "$wd"/"$refname"_multi-sample.vcf
          # create a fasta file with the sequences of all insertions
          while read line
          do    name=$(  echo $line | cut -d ' ' -f 1,2,3 | sed 's/\W/_/g')
                seq=$(   echo $line | cut -f 5 -d ' ')
                nonseq=$(echo $seq  | tr -d A | tr -d C | tr -d T | tr -d G | sed 's/\w//g' | wc -c )
                # check if seq is a seq, then print
                if   [ "$nonseq" -eq 1 ]
                then echo \>"$name"
                     echo "$seq"
                fi
                unset name seq nonseq
          done < <(grep -v '#' "$wd"/"$refname"_multi-sample.vcf \
                  | grep PASS \
                  | grep '\.INS\.' \
                  ) \
          > "$wd"/"$refname"_multi-sample_insertions.fasta

          # create a fasta file with the sequences of all deletions
          bedtools getfasta                                             \
                   -fi "${refs[$count]}"                                \
                   -fo "$wd"/"$refname"_multi-sample_deletions.fasta    \
                   -bed <(                                              \
                          grep -v '#' "$wd"/"$refname"_multi-sample.vcf \
                          | grep PASS                                   \
                          | grep '\.DEL\.'                              \
                          | vcf2bed                                     \
                          | cut -f 1,2,9                                \
                          | sed 's/\t[PI].*;END=/\t/g'                  \
                          | cut -d ';' -f 1                             \
                          )
     fi
     echo -e "$BLD""$BLU""In reference $refname the following variants were annotated by sniffles2 $NML"
     echo -e "$BLU"'Insertions:'     $(grep -v '#' "$wd"/"$refname"_multi-sample.vcf | grep PASS | grep '\.INS\.' | wc -l )
     echo -e "$BLU"'Deletions:'      $(grep -v '#' "$wd"/"$refname"_multi-sample.vcf | grep PASS | grep '\.DEL\.' | wc -l )
     echo -e "$BLU"'Recombinations:' $(grep -v '#' "$wd"/"$refname"_multi-sample.vcf | grep PASS | grep '\.BND\.' | wc -l )
done
conda deactivate

conda activate igv
# searching for sequences of interest in all denovo assemblies
echo -e "$BLD""$GRN""hecking blat searches marking sequences of interest in the de-novo assemblies $NML"
wd="$basedir"/igv_configs/
checkwd

for  s in "${samples[@]}"
do   name=$(echo "$s" | sed 's/\.fastq\.gz//g' )
     if    [ ! -f "$wd"/"$name"_blathits.psl ]
     then  checkprog blat
           blat  "$basedir"/denovo/"$name"/polished-medaka/consensus.fasta \
                 reference/queries.fasta                                   \
                 -t=DNA -q=DNA                                             \
                 "$wd"/"$name"_blathits.psl
     fi
done
conda deactivate

# export R markdown page with some statistics on the flye assemblies
if   [ ! -f ./denovo/flyestats.html ]
then Rscript -e "rmarkdown::render('flyestats.Rmd',output_file='./denovo/flyestats.html')"
else echo -e "$GRN""not updating denovo/flyestats html because it already exists. $NML"
fi


echo -e "$BLD""$GRN"'Script finished'"$NML"

# BRAKER3
###### Last updated: May 3, 2024
------------------------------------------------------------------------

## Prerequisites: 

- Soft-masked reference genome
- RNASeq data in .bam format, mapped to the reference genome
- Protein data, preferably downloaded from OrthoDB: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/

We already have the above files prepared. However, when running BRAKER3 on your own data, follow the below steps to prepare your reference genome, RNASeq data and protein data. BRAKER can also be run with only RNASeq data, or with only protein data. For more information on BRAKER, check their github: https://github.com/Gaius-Augustus/BRAKER#f19

1. Softmask the reference genome.
- You can use RepeatMasker and RepeatModeler, or Earl Grey. We'll use results from Earl Grey.

BRAKER doesn't work well with long fasta header names. Check your masked reference genome and confirm that there are no spaces, special symbols, or long names with grep:

```bash
grep ">" <filename>
```

The following code removes everything in fasta headers after a space:

```bash
 cut -d ' ' -f1 <masked_reference> > <outputfile>
```

BRAKER expects you to use hisat2 (https://daehwankimlab.github.io/hisat2/download/) to map RNASeq to your masked reference genome. Activate hisat2:

```bash
module load miniconda3/4.12-pws-472
conda activate hisat2
```

Now index the reference genome with hisat2-build, and then map the index to the .fastq files.

```bash
# Index the reference genome:
hisat2-build <masked-reference> <species-code>
# Map index to .fastq files:
hisat2 -x <species-code> -1 <RNASeq-dataset-1.fastq> -2 <RNASeq-dataset-2.fastq> -S <output.sam>
```

Convert .sam to .bam
```bash
conda deactivate
conda activate samtools # activate environment with samtools
samtools view -bS -o arcto-alignment.bam arcto_alignment.sam #convert .sam to .bam
samtools sort arcto-alignment.bam > arcto-sorted.bam #sort .bam
```

BRAKER has a lot of dependencies, including AUGUSTUS, which also has many dependencies. This makes it a pain to install manually. Because of licensing issues, there isn't an active conda environment kept for BRAKER. Instead, they setup a singularity to run BRAKER. If you have singularity on your HPC (we do), setup BRAKER as follows:

```bash
module load singularity
singularity build braker3.sif docker://teambraker/braker3:latest # this downloads the .sif file which singularity will need to run BRAKER

# download the following config folder from AUGUSTUS - this needs to be in a writable directory to run BRAKER
# We're going to download it from the BRAKER singularity
singularity exec braker3.sif cp -r Augustus/config .

# download three test files to confirm BRAKER is running correctly:
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
```

Edit test files with location of AUGUSTUS config file by adding the below code to all three test files.
```bash
--AUGUSTUS_CONFIG_PATH=<file_path_to_AUGUSTUS_config_directory>
```

Test setup worked by running test files:
```bash
bash test1.sh
bash test2.sh
bash test3.sh
```

Running BRAKER3:

Create a new job script for running BRAKER:

```bash
nano braker.job
```

```bash

#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
#SBATCH --ntasks=12   # number of processor cores (i.e. tasks)
#SBATCH --mail-user=youremail@email.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name= <name>

# Load modules
module load singularity

# Get the SIF file into the path
export BRAKER_SIF=<file_path_to.sif_file>
export AUGUSTUS_CONFIG_PATH=<file_path_to_AUGUSTUS_config_file>

# Run singularity
singularity exec braker3.sif braker.pl \
--genome=data/arcto-renamed.fasta \
--bam=data/arcto-sorted.bam \
--prot_seq=data/Arthropoda.fa \
--workingdir=arcto-grandis \
--gff3 \
--threads=12 \
--species=Arctopsyche-grandis \
--AUGUSTUS_CONFIG_PATH=<file_path_to_AUGUSTUS_config_directory>

```

It's important to understand what each of the above lines is doing. 
--genome references the softmasked reference genome is
--bam references the .bam RNASeq file mapped to your genome
--prot_seq is the protein database downloaded from OrthoDB
--workingdir specifies the working directory all files will be saved into - change this each time you run BRAKER so data isn't overwritten
--gff3 requests output to be in gff3 format
--species is a unique species identifier
--AUGUSTUS_CONFIG_PATH is the file path to the config directory downloaded above

## Examining the output:

```bash
wget https://raw.githubusercontent.com/Gaius-Augustus/GALBA/master/scripts/analyze_exons.py
chmod u+x analyze_exons.py
./analyze_exons.py -f arcto-grandis/braker.gtf
./analyze_exons.py -f arcto-grandis/Augustus/augustus.hints.gtf
./analyze_exons.py -f arcto-grandis/GeneMark-ETP/genemark.gtf
```

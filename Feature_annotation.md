# BRAKER3
###### Last updated: May 3, 2024
------------------------------------------------------------------------

## Prerequisites: 

- Soft-masked reference genome
- RNASeq data in .bam format, mapped to the reference genome
- Protein data, preferably downloaded from OrthoDB: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/

We already have the above files prepared. However, when running BRAKER3 on your own data, follow the below steps to prepare your reference genome, RNASeq data and protein data. BRAKER can also be run with only RNASeq data, or with only protein data. For more information on BRAKER, check their github: https://github.com/Gaius-Augustus/BRAKER#f19

1. Softmask the reference genome.
- Install RepeatModeler and RepeatMasker inside conda environment if not already installed:
```bash
conda install -c bioconda -c conda-forge repeatmodeler
conda install bioconda::repeatmasker
```

- Build a database for Repeat Modeler, this should go fairly quickly.

```bash
BuildDatabase -name <name_of_database> <reference.fasta>
```

Now you can softmask the genome, it's recommended to use multiple threads:
```bash
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
#SBATCH -J "soft-mask"   # job name
#SBATCH --mail-user=youremail@email.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# Activate environment with RepeatModeler and RepeatMasker
source ~/.bashrc
conda activate hfibroin

# Run code
RepeatModeler -database <name_of_database> -threads 24 -LTRStruct
# Name of database should match what was used above
RepeatMasker -threads 24 lib <database>-families.fa -xsmall <reference.fasta>
# -xsmall softmasks the genome
```

BRAKER expects you to use hisat2 to map RNASeq to your masked reference genome. You can find instructions on installing hisat2 here: https://daehwankimlab.github.io/hisat2/download/

```bash
# Index the reference genome:
./hisat2/hisat2-build <masked-reference> <species-code>
# Map index to .fastq files:
./hisat2/hisat2 -x <species-code> -1 <RNASeq-dataset-1.fastq> -2 <RNASeq-dataset-2.fastq> -S <output.sam>
```

Convert .sam to .bam
```bash
conda activate braker3 # activate environment with samtools
samtools view -bS -o arcto-alignment.bam arcto_alignment.sam #convert .sam to .bam
samtools sort arcto-alignment.bam > arcto-sorted.bam #sort .bam
```

BRAKER has a lot of dependencies, including AUGUSTUS, which also has many dependencies. This makes it a pain to install manually. Because of licensing issues, there isn't an active conda environment kept for BRAKER. Instead, they setup a singularity to run BRAKER. If you have singularity on your HPC (we do), setup BRAKER as follows:

```bash
module load singularity
singularity build braker3.sif docker://teambraker/braker3:latest # this downloads the .sif file which singularity will need to run BRAKER

# download the following config folder from AUGUSTUS - this needs to be in a writable directory to run BRAKER
git clone XXX

# download three test files to confirm BRAKER is running correctly:
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
```

Test setup worked by running test files:
```bash
bash test1.sh
bash test2.sh
bash test3.sh
```

Running BRAKER3:

```bash
singularity exec braker3.sif braker.pl \
--genome=data/arcto_4.fasta.masked \
--bam=data/arcto-sorted.bam \
--prot_seq=data/Arthropoda.fa \
--workingdir=arct
--gff3 \
--species=Arcto-test \
--AUGUSTUS_CONFIG_PATH=/nobackup/private/lifesci/fslg_lifesciences/ssamant/braker3/config/
```

It's important to understand what each of the above lines is doing. 
--genome specifies where the softmasked reference genome is
--bam the .bam RNASeq file mapped to your genome
--prot_seq is the protein database downloaded from OrthoDB
--workingdir specifies the working directory all files will be saved into - change this so data isn't overwritten
--gff3 requests output to be in gff3 format
--species is a unique species identifier
--AUGUSTUS_CONFIG_PATH should match what you get when you move into the config folder and run "pwd".


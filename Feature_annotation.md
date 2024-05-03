# BRAKER3
###### Last updated: May 3, 2024
------------------------------------------------------------------------

## Prerequisites: 

- Soft-masked reference genome
- RNASeq data in .bam format, mapped to reference genome
- Protein data, preferably downloaded from OrthoDB: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/

Build a database for Repeat Modeler, this should go fairly quickly.

```bash
BuildDatabase -name arcto arcto_4.fasta
```

Softmask genome, it's recommended to use multiple threads:
```bash
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
#SBATCH -J "soft-mask"   # job name
#SBATCH --mail-user=ssmit038@ucr.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL



# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source ~/.bashrc
conda activate hfibroin # activate environment with RepeatModeler and RepeatMasker
RepeatModeler -database arcto -pa 24 -LTRStruct
RepeatMasker -pa 24 lib arcto-families.fa -xsmall arcto_4.fasta
```

Map RNAseq to masked genome
```bash
./hisat2/hisat2-build arcto_4.fasta.masked arcto-4
./hisat2/hisat2 -x arcto-4 -1 SRR2083574_1.fastq -2 SRR2083574_2.fastq -S arcto_alignment.sam
```

Convert .sam to .bam
```bash
conda activate braker3 # activate environment with samtools
samtools view -bS -o arcto-alignment.bam arcto_alignment.sam
samtools sort arcto-alignment.bam > arcto-sorted.bam
```

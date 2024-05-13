# Whole Genome Assembly

Hifiasm is a tool for assembling genomes, especially for PacBio HiFi reads.

We need to install Hifiasm first to assemble the reads of our organisms' genomes.

We’ll use conda to create an environment and install hifiasm, using the following command:   

```bash
conda create -n hifiasm -c bioconda hifiasm
```
The '-n' option creates a new conda environment named hifiasm, and the "-c" option installs hifiasm from the bioconda channel into the hifiasm environment. 

Press "y" to install hifiasm.  

You can activate the environment with: 

```bash
conda activate hifiasm
```

With the active hifiasm environment, you can check if hifiasm was installed correctly by running the following command: 

```bash
hifiasm --version
```

The output should be the version number of our newly installed hifiasm software, e.g. "0.16.1-r375." 

Then, you can deactivate the active hifiasm environment with: 

```bash
conda deactivate
```

Now, we'll run a batch job (i.e., a task submitted to the supercomputer) with hifiasm.  

We'll first create a job script, a file with parameters and commands to run our task. 

You can create a job script by using BYU's [Job Script Generator](https://rc.byu.edu/documentation/slurm/script-generator). 


Fill out the following parameters and options for your job script:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Limit this job to one node: [select this option]*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Number of processor cores across all nodes: 24* 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Memory per processor: 4 GB* 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Walltime: 24 hours* 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Job name: [add your job’s name]*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Receive email for job events: [click on begin, end, abort]*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Email address: [add your email address]*


Then, scroll down and click on "Copy Script to Clipboard."

Go back to your terminal window.

Use your preferred text editor (vim, vi, nano) to create a text file and name it [name of job].sh 

```bash
nano [name of job].sh 
```

Paste your job script in your newly created text file.  

Scroll to the bottom of your text file, add a few lines of space, and then include the lines below to activate hifiasm and run it along with some running options. 

```bash
source ~/.bashrc 
conda activate hifiasm

hifiasm -o [Prefix of output file] -t 24  [input reads]
```

 The "-t" stands for the number of CPUs or processor cores we’ll use for the genome assembly, which is 24 (previously specified in the job script). 

 Your job script should look something like this: 

```bash
 #!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4096M   # memory per CPU core
#SBATCH -J "genome_assembly "   # job name
#SBATCH --mail-user=mgjijon@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


source ~/.bashrc 
conda activate hifiasm


hifiasm -o caddisfly_genome.asm -t 24 caddisfly_genome.fq.gz
```

 If you'd like to learn more about running hifiasm, here's the [README file](https://github.com/chhylp123/hifiasm).  

Save the changes you made and exit your text file window.  

Now, you can run hifiasm to assemble your organism's genome:  

```bash
sbatch [job name].sh
```
Hit enter and you're done!

To check if your job is running and the state of all your jobs, run the following command:

```bash
squeue -u [username]
```

An "R" should appear next to your hifiasm job number. 

The job should take several hours. 

After the job finish running, several output files will be created, including BED, GFA, and BIN files. 

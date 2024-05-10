# Genome Assembly

Hifiasm is a tool for assembling genomes, especially with PacBio HiFi reads.

We need to install Hifiasm first to assemble the reads of our organisms' genomes.

We’ll use conda to install hifiasm and create a conda environment to host hifiasm.  

```bash
conda install -c bioconda hifiasm
```
Press y to install hifiasm 

You can activate the environment using 

```bash
conda activate hifiasm
```

To run hifiasm, you will first create a job script. 

You can do this by using BYU's [Job Script Generator] (https://rc.byu.edu/documentation/slurm/script-generator). 

Fill out the following parameters and options for your job script:

Limit this job to one node: [select this option]
Number of processor cores across all nodes: 32 
Memory per processor: 4 GB 
Walltime: 48 hours 
Job name: [add your job’s name]
Receive email for job events: [click on begin, end, abort]
Email address: [add your email address]

Then, scroll down and click on "Copy Script to Clipboard."

Go back to your terminal window.

Use your preferred text editor (vim, vi, nano) to create a text file and name it [name of job].sh 

Paste your job script in your newly created text file.  

Scroll to the bottom of your text file and add a few lines of space, and then add the commands to activate hifiasm as well as to run it along with some running options. 

```bash
source ~/.bashrc 
conda activate hifiasm

hifiasm -o [Prefix of output file] -t 32  [input reads]
```

 Th t stands for the number of CPUs or processor cores we’ll be using for the assembly, which is 32 (this was previously specified in the job script, too). 

 If you'd like to learn more about how to run hifiasm, here's the [README file](https://github.com/chhylp123/hifiasm).  

Save the changes you made to the text file.  

Now, you can run hifiasm to assemble our organism's genome:  

```bash
sbatch [job name].sh
```

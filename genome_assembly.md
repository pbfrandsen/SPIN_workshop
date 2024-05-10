Hifiasm is a tool for assembling genomes, especially with PacBio HiFi reads.
We need to install Hifiasm first to assemble the reads of our organisms' genomes.

We’ll use conda to install hifiasm and create a conda environment to host hifiasm.  
```bash
conda install -c bioconda hifiasm
```
Press y to install hifiasm 

you can activate the environment using 

conda activate hifiasm 

#to run hifiasm you will create a job script first. 
#this can create the job script here BYU Job Script Generator 



#for the parameters, you should add:
#Limit this job to one node: [select this option]
#Number of processor cores across all nodes: 32 
#Memory per processor: 4 GB 
#Walltime: 48 hours 
#Job name: [add your job’s name]
#Receive email for job events: [click on begin, end, abort]
#Email address: [add your email address]

#Scroll down and click on 

#use your preferred text editor (vim, vi, nano) to create a text file and name it [name of job].sh 

#paste your job script in the text file 
#at the bottom of your text file, you will include adding the command to run hifiasm as well as the options you will be using while running hifiasm 
#add these lines: 

source ~/.bashrc 
conda activate hifiasm 

hifiasm -o [Prefix of output file] -t 32  [input reads] 

#the t is the number of CPUs or processor cores we’ll be using, which is 24 (previously specified in the job script, too)
#if you want to learn more about how to run this tool, you can read the hifiasm’s README file 

#save the change you made to the text file 

#now you can run hifiasm and assemble our genomes: 

sbatch [job name].sh 

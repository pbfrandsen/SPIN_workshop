# 03: Assembly QC

Now that we have a whole genome assembly, we're going to do some quick checks to see how it looks. The first thing that we will do is estimate the contig statistics with [`QUAST`](https://github.com/ablab/quast). Then we will run `BUSCO` to get a handle on genome completeness.

### QUAST

First, make sure that your assembly file was written appropriately from `hifiasm`. In the last lab, we converted three files from `gfa` format to `fasta` format. The one that we will use from here on is the one that ends in `p_ctg.fasta`. This is the `primary` assembly. Perhaps sometime soon, we'll have more tools that work directly on genome graphs that include the haplotypic diversity in the genome, but most tools currently work with a single fasta file representation of the genome.

Next, we'll go ahead and run `QUAST` on that assembly file. First, load the `conda` module and then activate the `QUAST` environment.

```
module load miniconda3/4.12-pws-472
conda activate quast
```

Then we can simply run quast on the assembly file. It is pretty fast so we can run it interactively and won't need to put it into a job file.

```
quast <assembly_name>.p_ctg.fasta
```

Now it will run on your genome and the results will be added to `quast_results/latest`. Once your run is complete, you can navigate to that folder. There should be a file called `report.pdf`. Go ahead and download that file with `scp` and examine it's contents.

### BUSCO

[`BUSCO`](https://busco.ezlab.org/) was designed to evaluate your genome assembly for the presence of universal single copy orthologs. These gene sets are generated from [`OrthoDB`](https://www.orthodb.org/).

Create a new job script using the [Job Script Generator](https://rc.byu.edu/documentation/slurm/script-generator). Make sure that the job is limited to one node, select 4 processor cores, 4 GB of RAM per processor, and a wall-time of 12 hours. Copy the resulting script into a new job file called `busco.job`.

Copy the `busco_downloads` folder from the following path to your working folder (the one containing your genome assembly file).
```
cp -r ~/fsl_groups/fslg_pws472/compute/lab5/data/busco_downloads ./
```
d. Open the job file and add the following commands to the bottom, separated by line breaks:

```
module load miniconda3/4.12-pws-472
conda activate busco
busco -o plodia_busco -i arctopsyche.asm.p_ctg.fasta -l insecta_odb10 -c 4 -m genome --offline
```

Note, that you should substitute your particular genome assembly file for "arctopsyche.p_ctg.fasta". NOte, we are using the `insecta_odb10` dataset. For some taxa, there are more specific datasets that you could use. For example, there is a `lepidoptera_odb10` dataset and for holometabolous insects, there is a `endoptergyota_odb10` dataset.

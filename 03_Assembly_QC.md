# 03: Assembly QC

Now that we have a whole genome assembly, we're going to do some quick checks to see how it looks. The first thing that we will do is estimate the contig statistics with [`QUAST`](https://github.com/ablab/quast). Then we will run `compleasm` to get a handle on genome completeness.

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

### Compleasm

[`compleasm`](https://github.com/huangnengCSU/compleasm) is a re-implementation of [`BUSCO`](https://busco.ezlab.org/), which was designed to evaluate your genome assembly for the presence of universal single copy orthologs. These gene sets are generated from [`OrthoDB`](https://www.orthodb.org/).

First copy the lineage files from `~/fsl_groups/fslg_nanopore/compute/genomics_workshop_byu_may_2/mb_downloads`. Make sure you copy that folder into the same folder that holds your genome assembly file. For example, if your genome assembly file was in, e.g. `~/compute/genome_workshop`. You could copy it by first navigating to that folder:

```
cd ~/compute/genome_workshop
```

And then copy it over with:

```
cp -r ~/fsl_groups/fslg_nanopore/compute/genomics_workshop_byu_may_2/mb_downloads .
```

Next, you'll want to run `compleasm`. You will want to make a new job file for this from the [BYU job script generator](https://rc.byu.edu/documentation/slurm/script-generator).

Choose 12 threads, 5 GB of RAM per thread and 4 hours.

Copy and paste that into a new file called, e.g. `compleasm.job`. At the bottom of your script add:

```
module load miniconda3/4.12-pws-472
conda activate compleasm
compleasm run -a arctopsyche.p_ctg.fasta -o compleasm_out -l endopterygota -t $SLURM_NPROCS
```

Note, that you should substitute your particular genome assembly file for "arctopsyche.p_ctg.fasta". For caddisflies, the endopterygota gene set is the best set, however, if you have an insect that is not in endopterygota, you might just want to use the "insecta" dataset. To do this, you would substitute "insecta" for "endopterygota" in the job file.
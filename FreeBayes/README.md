**FreeBayes pipeline**
-------------

This FreeBayes SNP calling pipeline follows is similar to what is implemented by the [snippy package](https://github.com/tseemann/snippy) (Seemann, 2018). In brief:  
* *Step 1:* Sequence data is aligned to reference sequences.
* *Step 2:* Clipped alignments and PCR duplicates are removed.
* *Step 3:* Joint SNP calling is carried out.
* *Step 4:* Called SNPs are filtered based on quality thresholds.
* *Step 5:* SNP calls were converted to a user-friendly format.

For additional information on SNP calling parameters, filtering, and conversion steps, see the [snippy package repo](https://github.com/tseemann/snippy) and the pipeline parameters below.

---
**Configuring pipeline**

All of the dependancy programs are included in the snippy conda package. This can be setup following instructions on the [snippy package repo](https://github.com/tseemann/snippy). Alternatively, each dependency can be downloaded individually, and its binary can be added to the PATH environment variable. These include:  
* *bwa mem*
* *samtools*
* *bcftools*
* *bedtools*
* *freebayes / freebayes-parallel*
* *vcflib*
* *vt*
* *snpEff*
* *samclip*
* *seqtk*
* *snp-sites*
* *any2fasta*

---
**Input files**

Three input files are required:
* ***Sequence data****: Raw (or trimmed) fastq or fastq.gz reads.
* ***Reference sequence template file***: This is a fasta file containing all of your SNP regions. You should include ample flanking sequence on either side of your target SNP to ensure raw reads are mapped (e.g. ~75 bp) - this will depend on your input read length. The SNP site can be denoted as an `N`. The fasta sequence ID should be the SNP name. See an example in the `example_inputs` folder.
* ***SNP position file***: This is a tab separated file which includes a column of SNP names and the position of the SNP in your *'Reference sequence template file'*. See an example in the `example_inputs` folder.

---
**Running the pipeline**

Edit the varaibles at the top of the FreeBayes pipeline bash script in the current folder:
`FreeBayes_pipeline.sh`  

Most of the SNP calling and filtering thresholds follow the snippy package, though a few can be altered in this pipeline. Below is a description of the variables that need to be defined: 
* ***FASTQ***: Raw fastq sequences - see above.
* ***REFERENCE***: *'Reference sequence template file'* - see above.
* ***SNP_POSITIONS***: *'SNP position file'* - see above.
* ***SAMPLE***: Sample name to assign to your output files.
* ***MIN_COVERAGE***: Minimum coverage threshold for SNP call.
* ***THREADS***: How many threads the analysis should use.
* ***REMOVE_FILES***: A number of intermediate files are produced - these files are removed if this variable is set to `true`.

In addition, if you downloaded snippy with conda, the relevant conda environment should be activivated and the snippy binaries should be exported:
* ***export PATH***: Exporting the path where the snippy binaries were installed. For example, this could be something like:
`export PATH=/Users/username/anaconda3/envs/snippy_env/bin:$PATH`  
* ***conda activate***: This activates your conda environment. Change the `conda activate` line depending on what you have called your snippy conda environment.

Once all of the variables are set, and the paths/environments are defined, the pipeline can be run with:
```
bash FreeBayes_pipeline.sh
```

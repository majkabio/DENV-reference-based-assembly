# **DENV reference-based assembly**<br/>
This script is a custom pipeline to perform referenced-based assembly of trimmed / filtered FASTQ files <br/>
by **M Galarion**


<br/>

## IMPORTANT DEPENDENCIES!!! <br/>
Make sure the following tools are installed:<br/>
<br/>
&emsp;**minimap2**  https://github.com/lh3/minimap2 <br/>
&emsp;**samtools** https://www.htslib.org/ <br/>
&emsp;**bcftools** https://samtools.github.io/bcftools/ <br/>
&emsp;**medaka**  https://github.com/nanoporetech/medaka <br/>
<br/>
Once these tools are installed, you need to explicitly define the command usage for each tool <br/>
To do this, open the latest version of the bash script in a text editor <br/>
Starting from around line 29, edit the line that corresponds to the tool definition (inside the quotation marks) <br/>
 <br/>
For example:
```

### Define tools and how you call them in your current system

minimap2 = "minimap2"
samtools = "samtools"
bcftools = "bcftools"
medaka = ".~/medaka/venv/bin/activate"

```
The tool definition inside the quotation marks should correspond to how you call them in your current system<br/>
Save the bash file and close<br/>

<br/>

### The script is an executable file written in bash and performs the following steps:<br/>
1. Aligns the reads to a reference file using minimap2 <br/>
2. Generates VCF file using bcftools mpileup <br/>
3. Produces a draft consensus file using bcftools consensus <br/>
4. Extracts mapped reads only using samtools and aligns them back to the draft consensus using medaka_consensus <br/>
5. Medaka generates a polished consensus <br/>
6. Masks low depth positions in the polished consensus using a user-specified minimum threshold <br/>
7. Produces final alignment stats <br/>
<br/>

### Command line usage:
```

bash ref-based-assembly_v2.2.sh -i reads.fastq -r reference.fasta

```
**NOTE:** Always run the script inside the directory where your -i and -r files are located.

<br/>

### To print options and default values:
```

bash ref-based-assembly_v2.2.sh

```
```
Optional parameters:
-t: number of threads (default: 8)
-M: Medaka polishing model (default: r941_min_hac_g507m)
-m: minimum read length (default: 200)
-x: maximum read length (default: none)
-d: minimum read depth (default: 20)
-q: minimum read Q-score (default: 9)
```

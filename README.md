# autoANI

# Installation

No installation is required if dependencies are installed. See these optional steps to save some time running the program if you are going to be the sole person using the software.

For questions about your PATH, see [this link](https://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them) for a decent basic explanation.

1. Edit line 49 to include your e-mail address. Your e-mail address will only be used to query the NCBI databases using the E-utilites. This is to allow NCBI to contact you if there is a problem with your requests or if you are abusing the system. 

If you will be sharing this installation with others, set your EMAIL environment variable. 

```bash
export EMAIL yourname@yourschool.edu
```
OR:
```tcsh
setenv EMAIL yourname@yourschool.edu
```

Otherwise, make the change:

```perl
my $email = '';
```
to
```perl
my $email = 'youremail@yourschool.edu';
```

2. Edit line 44 to the path of your blast executables if they are not in your PATH.

3. If you wish to include autoANI.pl in your path, but don't want to move the ./scripts folder to the same location, put the full path to the scripts folder in line 46. For example:

  * Move autoANI.pl to a location in your PATH (e.g. /usr/bin or /usr/local/bin or ~/bin)
  * Add the full path to the unpacked scripts directory to line 44 (eg. ~/libs/autoANI/scripts)

## Dependencies

Required 

* BLAST+ version 2.2.31+
* Perl version 5.10.1 or higher
* NCBI Entrez Direct scripts - Run edirect-dl.pl from scripts folder. [See here](https://www.ncbi.nlm.nih.gov/books/NBK179288/) for more information
* Perl Modules
  * BioPerl (Bio::Perl) v1.7.000 [(1.007000)](http://search.cpan.org/~cjfields/BioPerl-1.007000/)
  * ~~BioPerl Eutilties (Bio::DB::EUtilities)~~ **NO LONGER REQUIRED**
  * Parallel Fork Manager (Parallel::ForkManager) v1.18 or higher

# Pipeline Description

This software facilitates **A**verage **N**ucleotide **I**dentity calculation for any number of genome comparisons. This is accomplished in several steps:

1. Generate BLAST DBs for all genomes.
2. Divide genomes into 1020nt chunks.
3. Use chunked genomes as queries against full length genomic BLAST DBs.
4. Find shared genomic regions based on cut-offs (70% coverage, 30% identity, by default)
5. Find average identity of all shared genomic regions.
6. Repeat for any number of genomes.

The location of the input files is saved so that you can repeat the analysis with different parameters. Furthermore, you can add new genomes to the analysis and re-use all of the prior calculations. 
The input is a set of complete or draft genome sequences each contained within their own FASTA file.  
The output is a tab-delimited file with all of the pairwise ANI values.

# Usage

This software is designed to be as straightforward to use as possible, with as minimal input as possible. Once the dependencies are installed, running autoANI.pl --help should provide a lengthy help message to remind you of all of the possible options with this software.

The FASTA files have to be formatted properly prior to use with autoANI, excluding any FASTA files that have been downloaded from NCBI (if the FASTA files have a recognizable accession number, they will be detected). This process can be accomplished with the fasta\_format.pl script located in the autoANI/scripts folder. The formatting script expects each genome to be contained in an individual FASTA file, with the format of strainName.fasta. For example, if you were using the Escherichia coli strain K12, the genome file would be named K12.fasta. You would then provide the species name with the -genome flag:

`autoANI/scripts/fasta_format.pl -genome 'Escherichia coli' -overwrite K12.fasta`

If you had several genomes with the appropriate file names in a particular directory, you would do this:

`autoANI/scripts/fasta_format.pl -genome 'Escherichia coli' -overwrite *.fasta`

You must be careful with this rename script, as it will overwrite each file in turn, so preferably you would make a copy of each genome to convert to autoANI format.

The autoANI format is straightforward. For example:

`>gnl|K12|contig1 [Escherichia coli K12]`

The strain identifier is in the second field, and is used to link all of the sequences from this organism together. The full strain name is in square braces in the FASTA description field.

After you have properly formatted FASTA files, you can start a run with this command:

`autoANI.pl -email youremail@yourinstitution.edu -threads 8 *.fasta`

The FASTA files will automatically be split into 1020nt chunks, BLAST DBs for each genome will be generated, and pairwise BLAST searches will be done. Results are stored in a temporary file until all of the searches are complete. Then, the analysis can be finished with `autoANI.pl -finish > ani.out` in the same folder.

You can then sort the tab delimited output file using the autoANI-sort.pl file in the scripts folder. You provide a reference genome name sort start the sorting by, and then the file will be sorted based on ANI values in an iterative fashion (use -help for more information).

# Example Dataset

There is an example dataset in the ./examples folder. Run autoANI.pl from the examples folder like so:

`autoANI.pl -email youremail@yourinstitution.edu -threads 8 ./fna/*.fasta`

This process should only take a few minutes. After that completes, then run:

`autoANI.pl -email youremail@yourinstitution.edu -finish`

To complete the process. Your example output will be printed to STDOUT (your screen).

# Algorithm

To better explain the technical aspects of this script, I am including this section on the algorithm. I hope that even without reading the Perl code, this section will help you understand the specific steps this program uses to calculate ANI values.

## Data Preparation

Before doing any BLAST searches, the software prepares the data in three ways.

1. The nucleotide files are chunked into 1020nt fragments.
2. BLAST databases are generated from the full length genomes.
3. Proper names for the output file are gathered from either NCBI or the input FASTA files (from the headers).

More specifically, these events take place before BLAST searches:

First the script ensures that the proper BLAST version is found, and that blast and makeblastdb are available to the script (generally present in your PATH environment variable). It then generates folders for several of the intermediate files, including chunked FASTA files (./fasta), BLAST output (./blast), and the path to the original input files (./queries). It then saves all of the input files in the queries.txt file in the ./queries folder, and ensures that there are no duplicate file names given in two separate locations.

Next, the files are read and chunked versions are written to ./fasta, with name data written to ani.keys for local files, or accession numbers written to ani.accn.tmp for data from NCBI. BLAST DBs are written to the ./db folder. Using the accession numbers gathered from the headers of the chunked files, the NCBI E-Utilities are queried (using BioPerl modules) to download relevant naming data, including genus, species, and strain names. Some data in NCBI is not saved in the usual format, with Organism='Genus spcies' and Strain='strain name', so several different locations are checked for strain names. The relevant data from NCBI is then also stored in ani.keys.

## BLAST Searches

Once the data is prepared as above, then the BLAST searches can be started. The BLAST settings specified in Goris et al. 2007 are followed, with the 'new' BLAST+ software:

`-X=150 -q=-1 -F=F`

as

`-xdrop_gap 150 -penalty -1 -reward 1 -dust no`

We additionally use the filters `-max_target_seqs 1 and -max_hsps 1` to reduce the amount of output data and only provide the best BLAST hits. Further, the output only contains percent identity, alignment length, number of identical residues, and evalue in the output. Before the algorithm revision in v1.2, we used the straight percent identity from the BLAST search to compute the ANI. Moving forward, we will use the number of identical residues divided by the chunk size, which is what is specified in Goris et al. 2007. To compare to the old algorithm, use -oldani to output the old values. The primary effect of this change is to reduce ANI values, especially those of more distant comparisons. ANI values at or near the species level (~95% ANI) are relatively unaffected by this change.

This script parallelizes the BLAST searches using the Parallel::ForkManager Perl module. This provides an additional speed benefit compared to using -num\_threads in BLAST, as not every step in the BLAST search is threaded.

## Calculating ANI Values

Each line of the BLAST output is searched to identify the query and subject genomes, and then checked to see if the result passes filter: longer than the length (70% of chunk size by default), and greater than the percent ID cutoff (30% by default), then the result is kept and used for ANI calculation. As explained above, the old ANI calculation uses the percent ID given by the BLAST program, while the new calculation divides the number of identical residues by the length. In both cases, the average is taken by adding all of the pairwise identities together and dividing by the number of matches. These values are stored in ani.out.tmp for future reference and so that all of the BLAST searches do not need to be parsed each time the output table is generated.

Output is printed (to STDOUT by default) with the query as each row and the subject as each column. Use the -finish command to generate the table after the BLAST searches are all complete. Genome names are extracted from ani.keys and ANI values are either calculated directly from the BLAST output or taken from the ani.out.tmp file. Output can be sorted with the autoANI-sort.pl script in the ./autoANI/scripts folder.

# Downloading Genomes from NCBI

Use the included download\_genomes.pl script to download genomes from NCBI for inclusion in the pipeline. Use the format `scripts/download_genomes.pl -type ani Pseudomonas-download` to download, for example, the Pseudomonas genomes available in NCBI. Any term from the NCBI taxonomy database is acceptable. The script also provides the count of the number of genomes that will be downloaded, as well as an estimate of the disk usage required after downloading the sequence. Leave the -download flag off to see the number of genomes and size without downloading; this is also useful for checking to see if your taxonomy term is valid.

# History

v2.0.0 - 2016-12-06 - Major update. No longer requires Bio::DB::EUtilities. Not compatible with previous version.

v1.3.0 - 2016-11-15 - Now throw away chunks at end of contigs (those that would be <1kb). Should improve accuracy for highly fragmented genomes.

v1.2.1 - 2016-09-21 - Added better SGE support/notes. SGE support is only available at Oregon State University on the [CGRB infrastructure.](https://shell.cgrb.oregonstate.edu)

v1.2.0 - 2016-09-12 - Fixed ANI calculation to be true to Goris et al 2007. This results in decreased ANI values, with the more distant relationships being more affected than close relationships.

v1.1.0 - 2016-07-15 - Revised ani.out.tmp file generation. First run of -finish will take a while, depending on dataset size. Subsequent runs (on that same dataset) will be fast.

v1.0.0 - 2016-04-29 - First revision released to GitHub.

# Planned Updates

* Remove BioPerl dependency completely. Just need time to do it.
* When NCBI removes GI dependency from E-Utilities, switch GIs for Accessions. Should reduce coding complexity by a small, but significant, amount.

# Description of Included Files

| Filename    | Description  |
| ----------- | ------------ |
| README.md   | This file    |
| LICENSE     | Description of GPL license |
| autoANI.pl  | Main program |
| examples/   | Folder with example files |
| scripts/autoANI-plot\_prep.pl | Helps prepare data for import into R for generating heatmaps |
| scripts/autoANI-sort.pl | Sorts output based on user-provided reference strain name |
| scripts/auto\_edirect.pl | Using accessions as input, downloads strain name and other information |
| scripts/edirect-dl.pl | Downloader for edirect tools from NCBI |
| scripts/fasta\_format.pl | Formats user-generated FASTA files |
| scripts/download\_genomes.pl | Downloads genomes from NCBI |

# Citations

Please cite:

```
Davis II EW, Weisberg AJ, Tabima JF, Grunwald NJ, Chang JH. (2016) Gall-ID: tools for genotyping gall-causing phytopathogenic bacteria. PeerJ 4:e2222. https://doi.org/10.7717/peerj.2222
```

If you use this software. Also, please cite these other publications as integral parts of autoANI:

```
Stajich JE, Block D, Boulez K, Brenner SE, Chervitz SA, Dagdigian C, Fuellen G, Gilbert JG, Korf I, Lapp H, Lehväslaiho H, Matsalla C, Mungall CJ, Osborne BI, Pocock MR, Schattner P, Senger M, Stein LD, Stupka E, Wilkinson MD, Birney E. (2002) The Bioperl toolkit: Perl modules for the life sciences. Genome Res. 12(10):1611-8. https://doi.org/10.1101/gr.361602

Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM. (2007) DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. IJSEM 57(1):81–91. http://doi.org/10.1099/ijs.0.64483-0

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421. https://doi.org/10.1186/1471-2105-10-421
```

# Credits

* [Edward Davis](mailto:davised.dev@gmail.com)
* [Dr. Jeff Chang](mailto:changj@science.oregonstate.edu)

# License

Please see the LICENSE file for more information.

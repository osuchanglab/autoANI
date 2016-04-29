# autoANI


# Installation

No installation is required if dependencies are installed. See these optional steps to save some time running the program if you are going to be the sole person using the software.

For questions about your PATH, see [this link](https://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them) for a decent basic explanation.

1. Edit line 47 to include your e-mail address. Your e-mail address will only be used to query the NCBI databases using the E-utilites. This is to allow NCBI to contact you if there is a problem with your requests or if you are abusing the system.

    Make the change:
    
    ```perl
    my $email = '';
    ```
    to
    ```perl
    my $email = 'youremail@yourschool.edu';
    ```

2. Edit line 42 to the path of your blast executables if they are not in your PATH.

3. If you wish to include autoANI.pl in your path, but don't want to move the ./scripts folder to the same location, put the full path to the scripts folder in line 44. For example:

  * Move autoANI.pl to a location in your PATH (e.g. /usr/bin or /usr/local/bin or ~/bin)
  * Add the full path to the unpacked scripts directory to line 44 (eg. ~/libs/autoANI/scripts)

## Dependencies

Required 

* BLAST+ version 2.2.31+
* Perl version 5.10.1 or higher
* Perl Modules
  * BioPerl (Bio::Perl) v1.6.923, other versions might work but are unsupported
  * BioPerl Eutilities (Bio::DB::EUtilities) v1.73
  * Parallel Fork Manager (Parallel::ForkManager) v1.18

# Pipeline Description

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

# History

v1.0 - 2016-04-29 - First revision released to GitHub.

# Credits

* [Edward Davis](mailto:davised.dev@gmail.com)
* [Dr. Jeff Chang](mailto:changj@science.oregonstate.edu)

# License

Please see the LICENSE file for more information.

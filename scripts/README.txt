
CRISPRdigger Release Notes 
"(version 1.0, 10th Dec. 2015)"

1. INTRODUCTION
===============
This directory contains source code related to CRISPRdigger, which is described in manual.
The purpose of the CRISPRdigger software is to identify CRISPR sequences in the genomes. 
"In fact, the output of this program can be used as input to Argo software(http://www.broadinstitute.org/annotation/argo/) " as a way to make the genomes visual.


2. ABOUT SCRIPTS && softwares
=============================
This directory contains the scripts of about CRISPRdigger software
with the following files:

"scripts/                   contains different script files, as described in the manual."

"softwares_all.rar 			contains all tools which need to be auto installed. Please change its directory into the same folder with scripts after decompressed or unzipped it."
"softwares/                 contains a few tools which need to be auto installed. Other tools can be auto download an installed by AutoconfigCRISPRdigger.pl."
Hint: You just need select one way in two directorys(softwares_all,softwares) to install the needed tools by AutoconfigCRISPRdigger.pl.

AutoconfigCRISPRdigger.pl   is a perl script to auto-install all the softwares which the CRISPRdigger script need. This script needs to be executed in the first step. Please do not change their directory optionally.

CRISPRdigger_manual.pdf    is the complete manual for the CRISPRdigger software.

CRISPRdigger.pl            is a perl script to identify all CRISPRs in one genome or DNA sequence formatted fasta.

batchCRISPR.pl             is a perl script to identify all CRISPRs in a batch of genomes formatted fasta. The script should be put in the same directory with the genomes.
                                
filter-stage-1.prl         is a perl script of RepeatScout software. This script will be substituted for the old one.

GoldenCRISPRs              is a CRISPR golden standard database of all bacteria and archaea.It was obtained from the website: http://crispr.u-psud.fr/crispr/ .                
                    
The parameters of programs are described in the manual.After the CRISPRdigger script was installed successfully, you can execute other scripts. But,you must modify the bioperl library path and some other softwares' path in the right way in the scripts.


3. EXAMPLE FOR RUNNING THE SOFTWARE 
===================================
a. pre-requirement
   Some softwares need to be installed successfully before execute the CRISPRdigger. The envionment variables also were modified by the script automatically.
   The command line like this which you only need to execute in the "scripts" directory:
         perl AutoconfigCRISPRdigger.pl  
   After executing the last script, the tips was given about how to further install the RepeatMasker in the screen or in the file:"softwares/other_config.txt".
   You maybe also need to source your file ".bash_profile" in your home directory, the command like: source .bash_profile

b. Execute the script to discover all CRISPRs
In the "scripts" directory,an genome example(example.fna) is given a demonstration for your testing and debugging purposes.
The command line like this:
         perl CRISPRdigger.pl -i  example.fna 


4. Feedback
===========
 Please contact the development team at: fengfengzhou@gmail.com or ruiquange@gmail.com to submit questions or feedback for us.

 
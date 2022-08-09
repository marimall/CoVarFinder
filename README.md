# CoVarFinder
---This title is under consideration! It's super basic but at the moment I can't think of anything else

Prerequisites

In order to use this code you need to install the following tools:

freebayes: available on https://github.com/freebayes/freebayes
snpEff: available on http://pcingola.github.io/SnpEff/  (also needs java)
Have the snpEff file in the same directory as CoVarFinder, otherwise change it manually in the code, line 42

CoVarFinder can either process aligned .bam files or already processed .vcf files
For already processed files, just use the command

python CoVarFinder --input_file <vcf_file_or_vcf_file_directory> 

While for the .bam files run

python CoVarFinder --input_file <bam_file_or_bam_file_directory> --preprocess True --reference <your_sars_reference_file>


This code is in a very primitive state but it's still very usefull for those who know python basics.Basically the script reads all the genomic mutations, translates them in aminoacid alterations and compares them with all the defining mutations of all Variants of Concern. The final output is a file containing all unique mutations,shared mutations and mutations that were not reported on my all_mutations file. In the future I want to expand to other variants as well and to variant frequencies which will be a bitch but I will figure!


Enjoy!

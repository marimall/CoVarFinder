#####################################################
#####################################################
##												   ##
##	 HERE LIES THE MAGNIFICENT CODE WHICH 	   ##
##		ASSIGNS SARS VARIANTS TO BAM FILES		 ##
##		WRITTEN ENTIRELY BY MARIA MALLIAROU		##
##	   	Genius, not a billionaire,				 ##
##			playgirl, philanthropist.			   ##
##												 ##
#####################################################
#####################################################

#-------------Dependencies--------------------------

##check if installed and if not "freebayes", "java" and "have snpEff in directory" 


#--------------Imports-------------------------------

from subprocess import call
import io
import os
try:
    import pandas as pd
except ModuleNotFoundError:
    install(pandas) ###install 
import json
try:
    from Bio import SeqIO
except ModuleNotFoundError:
    install(biopython)
import re


#--------------Functions-----------------------------

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def f2dict (input_fasta):
	'''
	Takes a fasta file and returns a dictionary with the header as key and the sequence as values
	'''
	input_file = open(input_fasta)
	my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
	return my_dict

def fasta_ref(file):
    with open(file, "r") as f:
        lines = [line.strip("\n") for line in f]
    seq = "".join(lines[1:])
    return seq

def gene_calling(vcf_table):
    mutated_genes = []
    for i in vcf_table.index:
        for gene in pos.index:
            if vcf_table.loc[i]["POS"] < pos(gene).split("-")[1] and vcf_table.loc[i]["POS"] > pos(gene).split("-")[0]:
                #print( vcf_table.loc[i]["INFO"].split("|")[9])
                mutated_genes.append("{}:{}".format(gene, vcf_table.loc[i]["INFO"].split("|")[9][2:]) )  
    return mutated_genes

def read_vcf(path):
	'''
	what do you think it does?It reads a vcf file, dah!
	'''
	with open(path, 'r') as f:
		lines = [l for l in f if not l.startswith('##')]
	return pd.read_csv(
		io.StringIO(''.join(lines)),
		dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
			   'QUAL': str, 'FILTER': str, 'INFO': str},
		sep='\t'
	).rename(columns={'#CHROM': 'CHROM'})

def is_aa(aa, aa_table):
	'''
	returns if given string is aminoacid
	'''
	if aa in aa_table:
		return True
	else:
		return False
		

def translate(seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon] if is_aa(codon , table) else "_"
    elif len(seq)%3 == 1 :
        for i in range(0, len(seq) - 1, 3):
            codon = seq[i:i + 3]
            protein += table[codon] if is_aa(codon , table) else "_"
    elif len(seq)%3 == 2 :
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            protein += table[codon] if is_aa(codon , table) else "_"
    return protein	
	
	
def triplet_position(pos):

    if int(pos)%3 == 0:
        return int(int(pos)-2)
    elif int(pos)%3 == 1:
        return int(pos)
    else:
        return int(int(pos)-1)
        

def all_mutations(my_file):     #####################################################3   this needs better files
    with open(my_file, "r") as f:
        lines = [line.strip("\n") for line in f]

    all_muts = {}
    for l in lines:
        if len(l.split(":")) == 3 :
            all_muts[":".join(l.split(":")[0:2])] =  l.split(":")[2].split(",")
    return all_muts
    
    
    
def gene_calling(vcf_file):
    vcf_table = read_vcf(vcf_file)
    mutated_genes = []
    pos = {
        'ORF1a':'266-13468', 'ORF1b':'13465-21555', 'S':'21563-25384', 'ORF3a':'25393-26220', 
        'E':'26245-26472', 'M':'26523-27191', 'ORF6':'27202-27387', 'ORF7a':'27394-27759', 
        'ORF7b':'27756-27887', 'N':'28274-29533', 'ORF10':'29558-29674', 'ORF8':'27894-28259'
        }
        
    for i in vcf_table.index:
        #print(vcf_table.loc[i])
        for gene in pos.keys():
            if int(vcf_table.loc[i]["POS"]) < int(pos[gene].split("-")[1]) and vcf_table.loc[i]["POS"] > int(pos[gene].split("-")[0]):
                #print( vcf_table.loc[i]["INFO"].split("|")[9])
                mutated_genes.append("{}:{}".format(gene, vcf_table.loc[i]["INFO"].split("|")[9][2:]) )  
    return mutated_genes


    


def substitution(mutation, sequence_dict):
    site = re.findall(r'(\d+)([ATGC])>([ATGC])', mutation[1])
    gene = mutation[0]
    pos = triplet_position(site[0][0])
    assert(site[0][1] == sequence_dict[gene][int(site[0][0]) -1])
    temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + site[0][2] + sequence_dict[gene][int(site[0][0]): ]
  
    before = translate(sequence_dict[mutation[0]][pos -1 : pos +3])
    after = translate(temp_seq[pos -1 : pos +3])
    

    return gene + ":" +before + str(int((pos+2)/3)) + after
                
def one_deletion(mutation, sequence_dict):
    site = re.findall(r'(\w+)del([ATGC])', mutation[1])
    gene = mutation[0]
    pos = triplet_position(site[0][0])
    assert(site[0][1] == sequence_dict[gene][int(site[0][0]) -1])
    temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + sequence_dict[gene][int(site[0][0]): ]
  
    before = translate(sequence_dict[mutation[0]][pos -1 : pos +3])
    after = translate(temp_seq[pos -1 : pos +3])
    
    return gene + ":" +before + str(int((pos+2)/3)) + after
    #return gene + ":" +before + str(int((pos+2)/3)) + "del")
    #return "ORF shift"
    
def duplication(mutation, sequence_dict):
    site = re.findall(r'(\w+)dup([ATGC]+)', mutation[1])
    gene = mutation[0]
    pos = triplet_position(site[0][0])
    assert(site[0][1] == sequence_dict[gene][int(site[0][0]) -1])
    temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + site[0][1]+ site[0][1] + sequence_dict[gene][int(site[0][0]): ]
    before = translate(sequence_dict[mutation[0]][pos -1 : pos +3])
    after = translate(temp_seq[pos -1 : pos +3])
    return gene + ":" +before + str(int((pos+2)/3)) + after
    #return "ORF shift"
    
def insertion(mutation, sequence_dict):
    site = re.findall(r'(\w+)_(\w+)ins([ATGC]+)', mutation[1])
    gene = mutation[0]
    pos = triplet_position(site[0][0])
    temp_seq = sequence_dict[gene][0:int(site[0][0])] + site[0][2] + sequence_dict[gene][int(site[0][0]): ]
    before = translate(sequence_dict[mutation[0]][pos -1 : pos +3])
    after = translate(temp_seq[pos -1 : pos +3])
    return gene + ":" +before + str(int((pos+2)/3)) + after
    #return "ORF shift"

def multiple_deletion(mutation, sequence_dict):
    site = re.findall(r'(\w+)_(\w+)del([ATGC+])', mutation[1])
    gene = mutation[0]
    pos1 = triplet_position(site[0][0])
    pos2 = triplet_position(site[0][1])
    #assert(site[0][2] == sequence_dict[gene][int(site[0][0]) -1 :int(site[0][1])])
    #temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + sequence_dict[gene][int(site[0][0]): ]
  
    before = translate(sequence_dict[mutation[0]][pos1 -1 : pos2 +3])
    #after = translate(temp_seq[pos -1 : pos +3])
    
    #return gene + ":" +before + str(int((pos+2)/3)) + after
    return gene + ":" +before[0] + str(int((pos1+2)/3))+ "_" +  before[-1] + str(int((pos2+2)/3)) + "del"
    #return "ORF shift"

def mut_translator(mutation, sequence_dict):
    mut_list = []
    for i in mutation:
        info = i.split(":")
        gene = info[0]
        if re.search(pattern="[ATGC]>[ATGC]", string=info[1]) : 
            mut_list.append(substitution(info, sequence_dict))
        elif re.search(pattern="[0-9]del[ATGC]$", string=info[1]):
            mut_list.append(one_deletion(info, sequence_dict))
        elif re.search(pattern="[1-9]dup[ATGC]+", string=info[1]):
            mut_list.append(duplication(info, sequence_dict))
        elif re.search(pattern="[0-9]+_[0-9]+ins[ATGC]", string=info[1]):
            mut_list.append(insertion(info, sequence_dict))
        elif re.search(pattern="[0-9]+_[0-9]+del[ATGC]+$", string=info[1]):
            mut_list.append(multiple_deletion(info, sequence_dict))
            print(multiple_deletion(info, sequence_dict))
        
        
        else:
            print(info)
            
       
        
        

######################################################

total_mutations = all_mutations("all_mutations.txt")  ###### get all of the mutations reported manually

genes = { i.split(":")[0] : str(f2dict("ncbi_dataset/data/gene.fna")[i].seq) for i in f2dict("ncbi_dataset/data/gene.fna")}  #####get correct positions with corresponding sequence of all SARS-COV-2 genes

muts = gene_calling("vcf_results/Sample1_freebayes_ann.vcf")       #######   this file in order to be correct needs preproccecing with varCalling from freebayes, change ref name (from 2019-nCoV to NC_045512.2 and genomic mutation indentification from snpEff)


sample_mutations = mut_translator(muts,genes)

unique_mutations = {}
shared_mutations = {}
no_hit = {}
for sm in sample_mutations:
    if sm in all_mutations and len(all_mutations[sm]) == 1 and len(all_mutations[sm][0]) > 1 :
        unique_mutations[sm] = all_mutations[sm]
    elif sm in all_mutations and len(all_mutations[sm]) > 1 :
        shared_mutations[sm] = all_mutations[sm]
    else:
        no_hit[sm] = sm
        
print(unique_mutations)
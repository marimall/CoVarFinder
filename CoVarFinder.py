#####################################################
#####################################################
##						   ##
##	 HERE LIES THE MAGNIFICENT CODE WHICH 	   ##
##        ASSIGNS SARS VARIANTS TO BAM FILES	   ##
##	  WRITTEN ENTIRELY BY MARIA MALLIAROU	   ##
##	  					   ##
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
import argparse


#--------------Functions-----------------------------
def preprocess_file(bam, ref, path):
    filename = bam.split("/")[-1][0:-4] if len(bam.split("/")) > 1 else bam[0:-4]
    
    os.system("freebayes -f {} -F 0.1 --pooled-continuous {} > {}/{}_freebayes.vcf".format(ref, bam , path, filename))
    os.system("sed 's/2019-nCoV/NC_045512.2/g' {}/{}_freebayes.vcf > {}/{}_freebayes_sed.vcf".format(path, filename, path, filename))
    os.system("java -jar snpEff/snpEff.jar ann NC_045512.2 {}/{}_freebayes_sed.vcf > {}/{}_snpEff.vcf".format(path, filename, path, filename))
    return "{}/{}_snpEff.vcf".format(path, filename)

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
        

def all_mutations(my_file):     
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
    #return gene + ":" +before + str(int((pos+2)/3)) + "ORF shift"

def multiple_deletion(mutation, sequence_dict):   ###########################this needs fixing 
    print(mutation)
 
    site = re.findall(r'(\w+)_(\w+)del([ATGC]+)', mutation[1])
    gene = mutation[0]
    pos1 = triplet_position(site[0][0])
    pos2 = triplet_position(site[0][1])
    assert(site[0][2] == sequence_dict[gene][int(site[0][0]) -1 :int(site[0][1])])
    #temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + sequence_dict[gene][int(site[0][0]): ]
    before = translate(sequence_dict[mutation[0]][pos1 -1 : pos2 +3])
    if len(site[0][2])% 3 != 0:
        return [gene + ":" + str(int(pos1/3 +1)) + "ORF shift due to deletion"]
    elif pos1 == int(site[0][0]):
        if len(site[0][2]) == 3 :
            return [gene + ":" +before[0] + str(int((pos1+2)/3))+ + "del"]
        else:
            return [  gene + ":" +before[0] + str(int((pos1+2)/3))+ "_" +  before[-1] + str(int((pos2+2)/3)) + "del"]
    else:
        temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + sequence_dict[gene][int(site[0][1]): ]
        if len(site[0][2]) == 3 :
            
            return [ gene + ":" +before[-1] + str(int((pos2+2)/3)) + translate(temp_seq[pos1-1:pos1+3]), gene + ":" +before[0] + str(int((pos1+2)/3) )+ "del"]
        else:
            return [ gene + ":" +before[-1] + str(int((pos2+2)/3)) + translate(temp_seq[pos1-1:pos1+3]), gene + ":" +before[0] + str(int((pos1+2)/3) )+ "_" +  before[-2] + str(int((pos2+2)/3)-1) + "del"]



def deletion_insertion(mutation, sequence_dict):  
    site = re.findall(r'(\w+)_(\w+)del([ATGC]+)ins([ATGC]+)', mutation[1])
    gene = mutation[0]
    pos1 = triplet_position(site[0][0])
    pos2 = triplet_position(site[0][1])
    assert(site[0][2] == sequence_dict[gene][int(site[0][0]) -1 :int(site[0][1])])
    to_return = []
    if len(site[0][2]) == len(site[0][3]):
        temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + site[0][3]+sequence_dict[gene][int(site[0][1]): ]
        before = translate(sequence_dict[mutation[0]][pos1 -1 : pos2 +3])
        after = translate(temp_seq[pos1 -1 : pos2 +3])
        aa_pos =list(range( int((pos1+2)/3) , int((pos2+2)/3) +1 ))
        for i in range(0,len(aa_pos)):
            to_return.append(gene + ":" + before[i] + str(aa_pos[i]) + after[i])
        return to_return
    elif len(site[0][2]) != len(site[0][3]) and (len(site[0][3]) - len(site[0][2])) % 3 == 0:

        temp_seq = sequence_dict[gene][0:int(site[0][0])-1] + site[0][3] +sequence_dict[gene][int(site[0][1]): ]
        shift = int((triplet_position(int(site[0][0]) + len(site[0][3])) + 2 ) /3)
        before = translate(sequence_dict[mutation[0]][pos1 -1 : pos2 +3])
        after = translate(temp_seq[pos1 -1 : pos2 +3])
        aa_pos =list(range( int((pos1+2)/3) , int((pos2+2)/3) +1 ))
        
        aa_subs = [ q1 for q1 in aa_pos if q1<=shift ]
        aa_del = [ q2 for q2 in aa_pos if q2>shift ]
        
        
        for i in range(0,len(aa_subs)):
            to_return.append(gene + ":" + before[i] + str(aa_pos[i]) + after[i])
        to_return.append( gene + ":" + before[len(aa_subs)] + str(aa_del[0]) + "_" + before[-1] + str(aa_del[-1]) + "del"  )
        return to_return
    else:
        return [gene + ":" + str(int(pos1/3 +1)) + "ORF shift"]
    

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
            for d in multiple_deletion(info, sequence_dict):
                print(d)
                mut_list.append(d)
            #print(multiple_deletion(info, sequence_dict))
        elif re.search(pattern="[0-9]+_[0-9]+del[ATGC]+ins[ATGC]+$", string=info[1]):
            for m in deletion_insertion(info, sequence_dict):
                mut_list.append(m)
         
        else:
            pass
    return mut_list
       
def unique_mutations1(sample_muts, all_muts):
    unique_mutations = {}
    for sm in sample_mutations:
        if sm in all_muts and len(all_muts[sm]) == 1 and len(all_muts[sm][0]) > 1 :
            unique_mutations[sm] = all_mutas[sm]
    return 


def get_report(path ,filename, uni, shared, no):
    with open("{}/{}_report.txt".format(path, filename), "w") as f:
        f.write("Final report of sample {}".format(filename) + "\n")
        f.write("Unique mutations found in this sample:" + "\n")
        for i in uni:
            f.write(i+ ":" + " ".join(uni[i])  + "\n")
        f.write("Shared mutations found in this sample:" + "\n")
        for k in shared:
            f.write(k + ":" + " ".join(shared[k]) + "\n")
        f.write("The folllowing mutations are not reported in our database" + "\n")
        for l in no:
            f.write(l + "\n")
        

def core_function(vcf, genes, total_mutations, work) :

    muts = gene_calling(vcf)       #######   this file in order to be correct needs preproccecing with varCalling from freebayes, change ref name (from 2019-nCoV to NC_045512.2 and genomic mutation indentification from snpEff)     
    sample_mutations = mut_translator(muts,genes)
    
    unique_mutations = {}
    shared_mutations = {}
    no_hit = []
    for sm in sample_mutations:
        if sm in total_mutations and len(total_mutations[sm]) == 1 and len(total_mutations[sm][0]) > 1 :
            unique_mutations[sm] = total_mutations[sm]
        elif sm in total_mutations and len(total_mutations[sm]) > 1 :
            shared_mutations[sm] = total_mutations[sm]
        else:
            no_hit.append(sm)
    sample_name = vcf.split("/")[-1][0:-4] if len(vcf.split("/")) > 1 else vcf[0:-4]
    print(work)
    os.system("mkdir {}/results".format(work))
    get_report( "{}/results".format(work) , sample_name ,unique_mutations, shared_mutations, no_hit)

    unique_variants = list(set(sum(unique_mutations.values(), [])))

    if len(unique_variants) == 1 :
        u_counter = 0 
        s_counter = 0
        for mvs in shared_mutations:
            if unique_variants[0] in shared_mutations[mvs]:
                u_counter += 1
            else:   ##### Unfortunately BA.2 does not have unique mutations so check if its in the list which the other variants is not
                if "BA.2" in shared_mutations[mvs]:
                    s_counter += 1
                else:
                    print(shared_mutations[mvs])
        if u_counter == len(shared_mutations) :
            print("Found unique mutations of Variant {}".format(unique_variants[0]))
        elif u_counter + s_counter == len(shared_mutations) :
            print("Found unique mutations of Variant {} and BA.2".format(unique_variants[0]))  
    else:
        u_counter = 0 
        s_counter = 0   
        for mvs in shared_mutations:
            if not set(unique_variants).isdisjoint(set(shared_mutations[mvs])):
                u_counter += 1
            elif set(unique_variants).isdisjoint(set(shared_mutations[mvs])) and "BA.2" in shared_mutations[mvs] :
                s_counter += 1
        if u_counter == len(shared_mutations) :
            print("Found unique mutations of Variant {}".format(",".join(unique_variants)))
        elif u_counter + s_counter == len(shared_mutations) :
            print("Found unique mutations of Variant {} and BA.2".format(",".join(unique_variants)))
    

##################### ------------your arguments-----------------------
parser = argparse.ArgumentParser(description='Wanna derive SARS-COV-2 variants from your .bam or .vcf file?This is the code for you!',  epilog = "author: Maria Malliarou <maria.malliarou.ger@gmail.com> v1.1" )

parser.add_argument('--input_data',  type = str, required = True, help = "Please provide your .vcf / .bam file or respective direcories for multiple recognision.If you select .bam you also need to give argument --preprocess as True" )  ###θα παίρνει ένα ή πολλα vcf αρχείο
###parser.add_argument('--action')  #
parser.add_argument('--preprocess',default  = False, type = bool, help = "If select True, don't provide vcf but your initial .bam file which is trimmned, sorted and SARS-COV-2 aligned")  ###
#parser.add_argument('--bam',  type = str, help = "Please provide your bam file or vcf direcory for multiple recognision if you have selected preprocess True")  ###θα παίρνει ένα ή πολλα vcf αρχείο
parser.add_argument('--reference',  type = str, help = "Please provide your reference fasta file used if you choose preprocess") 
parser.add_argument('--alignment_name', nargs = 1, default  = "2019-nCoV")
parser.add_argument('--report', default = True, type = bool, help = "The name of the report" )

args = parser.parse_args()


######################-----------Parse database files -------------------
database = all_mutations("all_mutations.txt")  ###### get all of the mutations reported manually
gene_sequence = { i.split(":")[0] : str(f2dict("ncbi_dataset/data/gene.fna")[i].seq) for i in f2dict("ncbi_dataset/data/gene.fna")}  #####get correct positions with corresponding sequence of all SARS-COV-2 genes

######################-----------Core Code-------------------------------

if os.path.isfile(args.input_data):
    sample = args.input_data.split("/")[-1][0:-4] if len(args.input_data.split("/")) > 1 else args.input_data[0:-4]
    os.system("mkdir CoVarCaller_{}_results".format(sample))
    if args.preprocess == True:
        os.system("mkdir CoVarCaller_{}_results/vcf_files".format(sample))
        core_function(preprocess_file(args.bam, args.reference), gene_sequence, database, "CoVarCaller_{}_results".format(sample) )
    else:
        core_function(args.input_data, gene_sequence, database, "CoVarCaller_{}_results".format(sample) )
        
elif os.path.isdir(args.input_data):
    files = os.listdir(args.input_data)
    files = [q for q in files if q.endswith(".bam")]
    #print(files)
    os.system("mkdir {}all_results".format(args.input_data))
    path = "{}all_results".format(args.input_data)

    for file in files:
        print("Processing file {}".format(file) )
        sample = file.split("/")[-1][0:-4] if len(file.split("/")) > 1 else file[0:-4]
        os.system("mkdir {}/CoVarCaller_{}_results".format(path ,sample))
        if args.preprocess == True:
            
            os.system("mkdir {}/CoVarCaller_{}_results/vcf_results".format(path , sample))
            new_vcf = preprocess_file("{}/{}".format(args.input_data,file), args.reference, "{}/CoVarCaller_{}_results/vcf_results".format(path , sample))
            print(new_vcf)
            core_function( new_vcf, gene_sequence, database, "{}/CoVarCaller_{}_results".format(path,sample) )

        else:
            core_function("{}/{}".format(path,file), gene_sequence, database, "{}/CoVarCaller_{}_results".format(path,sample) )





 

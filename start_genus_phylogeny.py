import Bio
from Bio import SeqIO
from Bio import SeqFeature
import sys
import os
import csv
import subprocess

... Genus phylogeny
this is composed by two scripts 1) start_genus_phylogeny.pyhton and 2) continue_genus_phylogeny.sh
both scripts need to be in the same folder and only this script needs to be modified

before starting you need to create 2 folders:
one folder called prorteinortho, where you add the proteinortho result .tsv file
one folder called prokka_annotation where you place all genomes annotation from prokka (we just need .faa and .gbk for this)
...

##before you start type below the name of your .tsv file#
name_of_pr_ortho_file="nameofmy.proteinortho.tsv"
input_path = "proteinortho/"

#then change the following paramenter: how many genomes are you using?
z = 41

#name of the path to the prokka annotation
input_path_prokka="/path/to/prokka_annotation/"


#now you don't need to change anything now#
#create all the folders
output_path_RAW= "coregenes_locustags/"
final_ouput_path = "coregenes/"

os.mkdir(output_path_RAW)
os.mkdir(final_ouput_path)
os.mkdir("coregenes_alinment")
os.mkdir("coregenes_alinment_cleaned")



proteinortho_results = csv.reader(open(input_path + name_of_pr_ortho_file, 'r'), delimiter='\t')
x = 3
y = x + z

r = range(3, y)


#now we need to create a list with names of the genomes used in proteinortho
ICElist=[]
for row in proteinortho_results:
    for i in range(x, y):
        
        ICElist.append(row[i])
    break

ICElist = [s.replace(".faa", "") for s in ICElist]

#we now extract the locus tag of each single copy conserved gene
for line in proteinortho_results:
    
    if int(line[0])==z and int(line[1])==z:
        #print('extracting locus tags from	 ' + line[3])
        out1 = open(output_path_RAW + line[3] + '_singlegene.txt', 'w')
        for i in r:
            out1.write(line[x] + ',')
            x +=1
        x = 3 
        out1.close()    


#yes we did it already
proteinortho_results = csv.reader(open(input_path + name_of_pr_ortho_file, 'r'), delimiter='\t')
header = next(proteinortho_results)
	
print('all core genes locus tags list have been created, now I will extract the sequences' )

num = 1


for nome in os.listdir(output_path_RAW):
    ofile = open(final_ouput_path + nome + "_coregenes.fasta", "w")
    mylist=open(output_path_RAW + nome , 'r')
    coregenes = mylist.read().split(',')
    
    print('I am now extracting from all the genomes the gene number   ' + str(num))
    num +=1
    for element in ICElist:
        print('parsing   ' + element)
        #mylist=open(output_path_RAW + nome , 'r')
        #coregenes = mylist.read().split(',')
        
        ofile.write(">" + str(element) + "\n")
    

###this loop here is made so the genes are in the same order in all cores, does not follow the order the genes are found on the genome###
        
        gbk_input = SeqIO.parse(input_path_prokka + element + ".gbk", "genbank")
        print("i am parsing  " + element)
        for genome in gbk_input:
            for gene in genome.features:
                if(gene.type =="CDS"):
                    locustag=gene.qualifiers['locus_tag'][0]
                    prseq=gene.qualifiers['translation'][0]
                    DNAseq=gene.extract(genome.seq)
                    if locustag in coregenes:
                        ofile.write(str(DNAseq) + "\n")
                        continue
                    continue
                continue
            continue
        continue

    ofile.close()

#now we call the other script which is in bash, not elegant but does the job
subprocess.call("./continue_genus_phylogeny.sh", shell=True)

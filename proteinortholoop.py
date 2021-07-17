import os
import sys
import csv

#with this script we run proteinortho in a loop by increasing protein identity threshold at each loop.
#it starts from 20 and increases by 5 at each step, you can change than by changing 'number' at line 11 and line 35

#the folder prokka_annotation should contain all the genomes (.faa) you want to run proteinortho on


number=20
number_str=str(number)
out1 = open('summary.csv', 'w')
out1.write('pridentity, core, pan' + '\n')
out1.close()

for i in range (16):
    outputfile='41_genomes_singles_' + number_str + '_pridentity'
    command='proteinortho -clean -singles -identity={} /prokka_annotation/*.faa -project={}'.format(number,outputfile)
    os.system(command)  
    proteinortho_results = csv.reader(open('41_genomes_singles_' + number_str + '_pridentity.proteinortho.tsv'), delimiter='\t')
    header = next(proteinortho_results)
    core=0
    pan=0
    for lines in proteinortho_results:
        if int(lines[0])==41:
            core+=1
        else:
            pan+=1
    numbcore=str(core)
    numbpan=str(pan)
    out1 = open('iterative_proteinortho_summary.csv', 'a')
    out1.write(number_str + ',' + numbcore + ',' + numbpan + '\n')
    out1.close()
    number+=5
    number_str=str(number)
	
#you can use the summary file to draw a graph with you favourite program

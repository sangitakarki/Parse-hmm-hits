import os
import sys
import re
import shutil
from Bio import SeqIO
from collections import defaultdict
import random
import operator
import subprocess
import shlex
from Bio import Phylo

summary = "rawout.txt"
euk_genomes = "protein_file"
final_folder = "output_folder"

cog2proteins = defaultdict(list)
protein_identity = dict()
#Open summary file
with open(summary) as summary:
	for line in summary:
		line = line.rstrip()
		if line.startswith("protein"):
			pass
		else:
			tabs = line.split("\t")
			genome_protein = tabs[0]
			genome_protein = genome_protein.split(".")
			if len(genome_protein)>2:
				genome = genome_protein[0]+".faa"
				genome_protein.pop(0)
				protein = '.'.join(str(v) for v in genome_protein)
				
			else:
				genome = genome_protein[0]+".faa"
				protein = genome_protein[1]

			score = float(tabs [6])
			hit = tabs [2]
			if score>=200:
				cog2proteins[hit].append(protein)
				protein_identity[protein] = hit
summary.close()


#Go through each default dict and get sequences
sequences_dict = defaultdict(list)
for key, value in cog2proteins.items():
	for filegenomes in os.listdir(euk_genomes):
		if filegenomes.endswith(".faa"):
			protein_file = os.path.join(euk_genomes, filegenomes)
			for record in SeqIO.parse(protein_file, "fasta"):
				record_name = record.id.split(".")
				record_name.pop(0)
				record_name = '.'.join(str(v) for v in record_name)
				# print(record_name)
				if record_name in value:
					protein_id = protein_identity[record_name]
					print(protein_id, record_name)
					sequences_dict[protein_id].append(record)

###Getting final files
for key, value in sequences_dict.items():
	output = os.path.join(final_folder, key+".faa")
	print(output)
	SeqIO.write(value, output, "fasta")

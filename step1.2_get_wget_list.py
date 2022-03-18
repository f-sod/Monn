#!/usr/bin/python3

# =================================== protein-ligand complex =======================
def get_pdbid_list():
	pdbid_list = []
	with open('./pdbbind_index/INDEX_general_PL.2020') as f:
		for line in f.readlines():
			if line[0] != '#':
				pdbid_list.append(line.strip().split()[0])
	print(f"pdbid_list,{len(pdbid_list)}")
	return pdbid_list
pdbid_list = get_pdbid_list()

with open('out1.2_pdbid_list.txt','w') as fw:
	for pdbid in pdbid_list:
		fw.write(pdbid+'\n')


# creation d'un fichie de sortie des liens pour obtenir les complexes des pdb ID
with open('out1.2_pdbbind_wget_complex.txt','w') as fw:
	for pdbid in pdbid_list:
		fw.write('https://files.rcsb.org/download/'+pdbid+'.pdb\n')

# ================================= ligand pdb file =================================
def get_pdbid_to_ligand():
	pdbid_to_ligand = {}
	with open('./pdbbind_index/INDEX_general_PL.2020') as f:
		for line in f.readlines():
			if line[0] != '#':
				ligand = line.strip().split('(')[1].split(')')[0]
				if '-mer' in ligand:
					continue
				elif '/' in ligand:
					ligand = ligand.split('/')[0]
				if len(ligand) != 3:
					continue
				pdbid_to_ligand[line[:4]] = ligand
	print(f"pdbid_to_ligand,{len(pdbid_to_ligand)}")
	return pdbid_to_ligand
pdbid_to_ligand = get_pdbid_to_ligand()

# creation d'un fichier de sortie des liens des ligands retrouve dans les complexes
with open('out1.2_pdbbind_wget_ligand.txt','w') as fw : 
	for pdbid in pdbid_list:
		if pdbid in pdbid_to_ligand:
			ligand = pdbid_to_ligand[pdbid]
			fw.write('https://files.rcsb.org/ligands/download/'+ligand+'_ideal.pdb\n')

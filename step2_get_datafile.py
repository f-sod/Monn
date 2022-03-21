#!/usr/bin/python3

__auteurs__ = ""
__update_by__="ODJE Floriane"
__data__ = "2022-03"

from rdkit import Chem
import numpy as np
import os
from rdkit.Chem import AllChem
import pickle

def get_fasta_dict():
	"""
	Fill up a dictionary wile reading 
	out1.6_pdbbind_seqs.fasta and out1.6_uniprot_uniprot_mapping.tab
	Uniprot Id are keys , sequences are values
	
	Return
	------
	Dictionnary: uniprot_dict of 5256 elements
	"""
	uniprot_dict = {}
	name,seq = '',''
	with open('out1.6_pdbbind_seqs.fasta') as f:
		for line in f.readlines():
			if line[0] == '>':
				if name != '':
					uniprot_dict[name] = seq
				name = line.split('|')[1]
				seq = ''
			else:
				seq += line.strip()
		uniprot_dict[name] = seq
	print(f"uniprot_dict step1,{len(uniprot_dict)}")
	
	with open('out1.6_uniprot_uniprot_mapping.tab') as f:
		for line in f.readlines()[1:]:
			lines = line.split('\t')
			for name in lines[0].split(','):
				if name not in uniprot_dict and lines[1] in uniprot_dict:
					uniprot_dict[name] = uniprot_dict[lines[1]]
	print('uniprot_dict step2',len(uniprot_dict))
	return uniprot_dict

def get_pdbid_to_uniprotid(seq_dict):
	"""
	Fill up a dictionary wile reading 
	INDEX_general_PL_name.2020 and out1.4_pdb_uniprot_mapping.tab
	PDB IDs are keys , Uniprot IDs are values
	
	Return
	------
	Dictionnary: pdbid_to_uniprotid  of 19302 elements
	"""
	pdbid_set = set()
	pdbid_without_official_uniprot_id  = []
	pdbbind_mapping_dict = {}
	with open('./pdbbind_index/INDEX_general_PL_name.2020') as f:
		for line in f.readlines():
			if line[0] != '#':
				lines = line.strip().split('  ')
				if lines[2] != '------':
					pdbbind_mapping_dict[lines[0]] = lines[2] 
					pdbid_set.add(lines[0])
				else:
					pdbid_without_official_uniprot_id.append(lines[0])
	print('pdbbind_mapping_dict',len(pdbbind_mapping_dict))
	print('pdbid without matching uniport id', len(pdbid_without_official_uniprot_id ))
	
	uniprot_mapping_dict = {}
	with open('out1.4_pdb_uniprot_mapping.tab') as f:
		for line in f.readlines()[1:]:
			lines = line.split('\t')
			for pdbid in lines[0].split(','):
				pdbid_set.add(pdbid)
				if pdbid in uniprot_mapping_dict:
					uniprot_mapping_dict[pdbid].append(lines[1])
				else:
					uniprot_mapping_dict[pdbid] = [lines[1]]
	print('uniprot_mapping_dict',len(uniprot_mapping_dict))
	
	"""
	Some pdb id did not get official / matching UniProt id in the file
	Index general name but got several in out1.4_pdb_uniprot_mapping.tab. 
	Those elements are retrieved in the list inter_list. For those elements,
	We check if at least a sequence (in seq_dict) is associated with the UniProt id,
	Those having a sequence in the dictionary sequences are compared (by length)
	The longest is associated with the PDB id in uniprot_mapping_dict.  
	"""
	
	keys_with_multiple_value=[]
	for pdbid in uniprot_mapping_dict:
		if len(uniprot_mapping_dict[pdbid]) != 1:
			keys_with_multiple_value.append(pdbid)
	inter_list=list(set(pdbid_without_official_uniprot_id ) & set(keys_with_multiple_value))
	print(len(inter_list), inter_list)
	for pdbid in inter_list:
		referenced_uniprot_id=[uniprot_id for uniprot_id in uniprot_mapping_dict[pdbid] if uniprot_id in seq_dict ]
		print((referenced_uniprot_id))
		for idx in range(len(referenced_uniprot_id)-1):
			if len(seq_dict[referenced_uniprot_id[idx]]) > len(seq_dict[referenced_uniprot_id[idx + 1]]):
				uniprot_mapping_dict[pdbid]=[referenced_uniprot_id[idx]]
			else:
				uniprot_mapping_dict[pdbid]=[referenced_uniprot_id[idx + 1]]
	
	pdbid_to_uniprotid = {}
	count = 0
	for pdbid in pdbid_set:
		if pdbid in uniprot_mapping_dict and len(uniprot_mapping_dict[pdbid]) == 1:
			pdbid_to_uniprotid[pdbid] = uniprot_mapping_dict[pdbid][0] # index 0 to remove the list form an take the first uniprot id 
		else:
			pdbid_to_uniprotid[pdbid] = pdbbind_mapping_dict[pdbid]
			if pdbbind_mapping_dict[pdbid] not in seq_dict:
				count += 1
	print('final pdbid_to_uniprotid', len(pdbid_to_uniprotid))
	print ('no sequence', count)
	print (pdbid_to_uniprotid)
	return pdbid_to_uniprotid

def get_pdbid_to_ligand():
	""" 
	Extract ligand IDs of all protein-ligand complex
	While parsing INDEX_general_PL.2020 

	Fill up a dictionary 
	PDB IDs are keys , ligands are values
	
	Return
	------
	Dictionnary: pdbid_to_ligand of 16726 elements
	"""
	pdbid_to_ligand = {}
	with open('./pdbbind_index/INDEX_general_PL.2020') as f:
		count_error = 0
		for line in f.readlines():
			if line[0] != '#':
				ligand = line.strip().split('(')[1].split(')')[0]
				if '-mer' in ligand:
					count_error += 1
					continue
				elif '/' in ligand:
					ligand = ligand.split('/')[0]
				elif len(ligand) != 3:
					count_error += 1
					continue
				pdbid_to_ligand[line[:4]] = ligand
	print('pdbid_to_ligand',len(pdbid_to_ligand))
	print('count_error: no pdb ligand ID', count_error)
	return pdbid_to_ligand

def get_pdbid_to_affinity():
	""" 
	Parse the initial file to retrieve the values of the affinity constants.
	Set the constants to the same scale. 
	Count the number of data items for which there is no correct affinity value

	Return
	------

	Dictionnary : pdbid_to_measure
	Dictionnary : pdbid_to_value
	"""
	
	pdbid_to_measure, pdbid_to_value = {}, {}   # value: -log [M]
	with open('./pdbbind_index/INDEX_general_PL.2020') as f:
		count_error = 0
		for line in f.readlines():
			if line[0] != '#':
				lines = line.split('/')[0].strip().split('  ')
				pdbid = lines[0]
				if '<' in lines[3] or '>' in lines[3] or '~' in lines[3]:
					#print lines[3]
					count_error += 1
					continue
				measure = lines[3].split('=')[0]
				value = float(lines[3].split('=')[1][:-2])
				unit = lines[3].split('=')[1][-2:]
				if unit == 'nM':
					pvalue = -np.log10(value)+9
				elif unit == 'uM':
					pvalue = -np.log10(value)+6
				elif unit == 'mM':
					pvalue = -np.log10(value)+3
				elif unit == 'pM':
					pvalue = -np.log10(value)+12
				elif unit == 'fM':
					pvalue = -np.log10(value)+15
				else:
					print(unit)
				pdbid_to_measure[pdbid] = measure
				pdbid_to_value[pdbid] = pvalue
	print ('count_error not = measurement', count_error)
	return pdbid_to_measure, pdbid_to_value

def get_mol_dict():
	"""
	Using RDKIT module 
	While parsing ‘Components-pub.sdf’ SDM type of object are obtained
	The proper name of the SDM object is added in the dictionary
	Names are keys , rdkit object are values
	Return
	------

	Dictionnary : mol_dict of 36308 elements
	"""
	mol_dict = {}
	mols = Chem.SDMolSupplier('Components-pub.sdf')
	for m in mols:
		if m is None:
			continue
		name = m.GetProp("_Name")
		mol_dict[name] = m
	print('mol_dict',len(mol_dict))
	return mol_dict

if __name__ == "__main__":

	# Function call
	mol_dict = get_mol_dict()
	uniprot_dict = get_fasta_dict()
	pdbid_to_uniprotid = get_pdbid_to_uniprotid(uniprot_dict)
	pdbid_to_ligand = get_pdbid_to_ligand()
	pdbid_to_measure, pdbid_to_value = get_pdbid_to_affinity()
	print('pdbid_to_measure', len(pdbid_to_measure), 'pdbid_to_value', len(pdbid_to_value))

	# Pickle dictionnary to re-use later
	with open('pdbid_to_ligand_dict','wb') as dw:
		pickle.dump(pdbid_to_ligand,dw,protocol=0)

	count_success = 0
	fw = open('out2_pdbbind_all_datafile.tsv', 'w')
	error_step1, error_step2, error_step3, error_step4 = 0, 0, 0, 0
	for pdbid in pdbid_to_uniprotid:
		if pdbid not in pdbid_to_ligand:
			error_step1 += 1
			continue
		if pdbid not in pdbid_to_measure:
			error_step2 += 1
			continue
		ligand = pdbid_to_ligand[pdbid]
		if ligand not in mol_dict:
			print ('ligand', ligand)
			error_step3 += 1
			continue
		inchi = Chem.MolToInchi(mol_dict[ligand])

		uniprotid = pdbid_to_uniprotid[pdbid]
		if uniprotid in uniprot_dict:
			seq = uniprot_dict[uniprotid]
		else:
			print ('uniprotid', uniprotid)
			error_step4 += 1
			continue
		
		measure = pdbid_to_measure[pdbid]
		value = pdbid_to_value[pdbid]
		fw.write(pdbid+'\t'+uniprotid+'\t'+ligand+'\t'+inchi+'\t'+seq+'\t'+measure+'\t'+str(value)+'\n')
		count_success += 1
	fw.close()
	print ('count_success',count_success )

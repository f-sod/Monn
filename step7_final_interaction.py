#!/usr/bin/python3

__auteurs__ = ""
__update_by__="ODJE Floriane"
__data__ = "2022-03"


import numpy as np
import pickle


def get_pdb_to_uniprot_map(result_dict):
	"""
	Generate a dict which has as a key the pdb id 
	and as value the chain of the protein sequence 
	and instead of the amino sequence gets the indexes by calling seq_with_gap_to_index
	for both query and target function
	
	Return 
	------
	
	dict : pdb_to_uniprot_map_dict{'pdbid':{'chain id': 'pdb_to_uniprot_idx'}}
	"""
	pdb_to_uniprot_map_dict = {}
	for name in result_dict:
		pdbid, chain = name.split('_')
		seq_target, seq_query, align, target_start, query_start = result_dict[name]
		ratio = float(align.count('|'))/float(len(seq_target.replace('-','')))
		if ratio < 0.9:
			continue
		
		target_idx_list = seq_with_gap_to_idx(seq_target)
		query_idx_list = seq_with_gap_to_idx(seq_query)
		pdb_to_uniprot_idx = get_target_idx(target_idx_list, query_idx_list, align, target_start, query_start)
		if pdbid in pdb_to_uniprot_map_dict:
			pdb_to_uniprot_map_dict[pdbid][chain] = pdb_to_uniprot_idx
		else:
			pdb_to_uniprot_map_dict[pdbid] = {}
			pdb_to_uniprot_map_dict[pdbid][chain] = pdb_to_uniprot_idx
	return pdb_to_uniprot_map_dict

def get_result_dict():
	"""
	Read and Parsed the sequence alignement file
	For a query name and query target the seq is stripped
	and added as a value in the dict as well as the alignement 
	
	Return
	------
	result_dict : 
	"""
	result_dict = {}
	with open('./smith-waterman-src/out6.3_pdb_align.txt','r') as f:
		i = -1
		seq_target, seq_query, align = '', '', ''
		pdb_ratio_dict = {}
		for line in f.readlines():
			i += 1
			if i%4 == 0:
				if 'target_name' in line:
					if len(seq_target) != 0:
						result_dict[target_name] = (seq_target, seq_query, align, target_start, query_start)
					target_name = line.strip().split(' ')[-1]
					print(f"target_name,{target_name}")
					seq_target, seq_query, align = '', '', ''
				else:
					seq_target += line.split('\t')[1]
					print(f"seq_target,{seq_target}")
			elif i%4 == 1:
				if 'query_name' in line:
					query_name = line.strip().split(' ')[-1]
					print(f"query_name,{query_name}")
				else:
					align += line.strip('\n').split('\t')[1]
					#print(f"align,{align}")
			elif i%4 == 2:
				if 'optimal_alignment_score' in line:
					for item in line.strip().split('\t'):
						if item.split(' ')[0] == 'target_begin:':
							target_start = int(item.split(' ')[1])
						elif item.split(' ')[0] == 'query_begin:':
							query_start = int(item.split(' ')[1])
				else:
					seq_query += line.split('\t')[1]
		return result_dict

def read_fasta():
	"""
	Read the fasta file generated at the previous step
	Collect while parsing lines, pdb and uniprot ID
	The sequences is saved in the variable seq
	In order to return a dictionnary 
	
	Return
	-------
	Dict: uniprot_seq_dict {'pdbid': (seq, uniprot id)}
	
	"""
	uniprot_seq_dict = {}
	with open('out6.1_target_uniprot_pdb.fasta','r') as f:
		for line in f.readlines():
			if line[0] == '>':
				pdbid = line.split('_')[0][1:]
				name = line.strip().split('_')[-1]
			else:
				seq = line.strip()
				uniprot_seq_dict[pdbid] = (seq,name)
		return uniprot_seq_dict

def seq_with_gap_to_idx(seq):
	"""
	Create a list of index for each sequences, 
	if a gap occurs with the seq, the value -1 is affected 
	else the index is affected 
	
	Return
	------
	list : idx_list
	"""
	idx_list = []
	i = 0
	for aa in seq:
		if aa == '-':
			idx_list.append(-1)
		else:
			idx_list.append(i)
			i += 1
	return idx_list

def get_target_idx(target_idx_list, query_idx_list, align, target_start, query_start):
	"""
	Convert the alignement-gap sequence of pdb to exact position (uniprot) 
	With -1 instead of gap and value for | 
	
	Return 
	------
	List : pdb_to_uniprot_idx 
	"""
	pdb_to_uniprot_idx = []
	for i in range(target_start-1):
		pdb_to_uniprot_idx.append(-1)
	for i in range(len(target_idx_list)):
		if target_idx_list[i] != -1:
			if align[i]  == '|' and query_idx_list[i] != -1:
				pdb_to_uniprot_idx.append(query_idx_list[i] + query_start-1)
			else:
				pdb_to_uniprot_idx.append(-1)
	return pdb_to_uniprot_idx

def get_interact_in_uniprot_seq(pdb_to_uniprot, uniprot_seq, seq_dict, residue_interact):
	"""
	Switch the position number of the amino-acids with the corresponding interaction index 
	Concatenante the interacting residue within a list and their corresponding index within an other
	
	Return
	-----
	list : interact_in_uniprot_seq_list : index of the intereaction amino acids 
	list : residue_records : 1 letter amino acids codes list of the interacting amino acids
	list : residue_bond_type : similar to residue_interact with the interaction index instead of seq_pos nb 
	
	"""
	interact_in_uniprot_seq_list = []
	residue_bond_type = []
	interact_in_uniprot_seq_set = set()
	residue_record = '' 
	for item in residue_interact:
		chain, idx = item[0][0], int(item[0][1:])    # chain, idx (pdb) of interact residue
		if chain not in  pdb_to_uniprot:
			continue
		sequence, idx_list = seq_dict[chain]       # pdb seuqnce, idx of chain
		if idx_list.count(idx) != 1:             # some positions in pdb sequence may have the same idx
			print('idx_list.count(idx) != 1')
		seq_pos = idx_list.index(idx)   # position along pdb sequence
		if seq_pos >= len(pdb_to_uniprot[chain]):
			continue
		if pdb_to_uniprot[chain][seq_pos] == -1:
			continue
		interact_idx = pdb_to_uniprot[chain][seq_pos]
		#if interact_idx not in interact_in_uniprot_seq_set:
		interact_in_uniprot_seq_set.add(interact_idx)
		interact_in_uniprot_seq_list.append(interact_idx)
		residue_bond_type.append((interact_idx,item[2]))
		residue_record += item[1]     # to  check
	return interact_in_uniprot_seq_list, residue_record, residue_bond_type

if __name__ == "__main__":
	
	with open('out1.2_pdbid_list.txt') as f:
		pdbid_list = [line.strip() for line in f.readlines()]

	with open('out4_interaction_dict','rb') as f:
		interaction_dict_from_pdb = pickle.load(f)
	
	pdbid_to_ligand = pickle.load(open("pdbid_to_ligand_dict","rb"))
	result_dict = get_result_dict()
	print('result_dict',len(result_dict))
	pdb_to_uniprot_map_dict = get_pdb_to_uniprot_map(result_dict)
	
	i = 0
	count_no_uniprot_map = 0
	count_not_in_uniprot_seq_dict = 0
	uniprot_seq_dict = read_fasta()
	interaction_dict_from_pdb_final = {}

	for item in interaction_dict_from_pdb:
		pdbid, ligand = item.split('_')
		i += 1
		if pdbid in ['4p4s', '5ayf']:  #sequences are not same
			continue
		if pdbid not in uniprot_seq_dict:
			count_not_in_uniprot_seq_dict += 1
			continue
		if pdbid not in pdb_to_uniprot_map_dict:
			count_no_uniprot_map += 1
			continue
		ligand = pdbid_to_ligand[pdbid]       # get ligand id
		uniprot_seq, uniprot_id = uniprot_seq_dict[pdbid]    # get uniprot sequence

		seq_dict = interaction_dict_from_pdb[pdbid+'_'+ligand]['sequence']      # get pdb seq_dict 
		residue_interact = interaction_dict_from_pdb[pdbid+'_'+ligand]['residue_interact']   # get pdb residue_interact
		
		assert residue_interact is not None

		pdb_to_uniprot = pdb_to_uniprot_map_dict[pdbid]
		interact_in_uniprot_seq_list, residue_record, residue_bond_type \
		= get_interact_in_uniprot_seq(pdb_to_uniprot, uniprot_seq, seq_dict, residue_interact)
		
		assert residue_record == ''.join(np.array([aa for aa in uniprot_seq])[interact_in_uniprot_seq_list].tolist()), pdbid #check if uniprot aa is the same as contact aa

		bond = interaction_dict_from_pdb[pdbid+'_'+ligand]['bond']
		atom_idx = interaction_dict_from_pdb[pdbid+'_'+ligand]['atom_idx']
		atom_name = interaction_dict_from_pdb[pdbid+'_'+ligand]['atom_name']
		atom_element = interaction_dict_from_pdb[pdbid+'_'+ligand]['atom_element']
		atom_interact = interaction_dict_from_pdb[pdbid+'_'+ligand]['atom_interact']
		atom_bond_type = interaction_dict_from_pdb[pdbid+'_'+ligand]['atom_bond_type']
		
		assert len(atom_idx) != 0
		
		interaction_dict_from_pdb_final[pdbid] = {}
		interaction_dict_from_pdb_final[pdbid]['ligand'] = ligand
		interaction_dict_from_pdb_final[pdbid]['atom_idx'] = atom_idx
		interaction_dict_from_pdb_final[pdbid]['atom_name'] = atom_name
		interaction_dict_from_pdb_final[pdbid]['atom_element'] = atom_element
		interaction_dict_from_pdb_final[pdbid]['atom_interact'] = atom_interact
		interaction_dict_from_pdb_final[pdbid]['atom_bond_type'] = atom_bond_type
		interaction_dict_from_pdb_final[pdbid]['uniprot_id'] = uniprot_id
		interaction_dict_from_pdb_final[pdbid]['uniprot_seq'] = uniprot_seq
		interaction_dict_from_pdb_final[pdbid]['interact_in_uniprot_seq'] = interact_in_uniprot_seq_list
		interaction_dict_from_pdb_final[pdbid]['residue_bond_type'] = residue_bond_type
		
	print(f"interaction_dict_from_pdb_final",{len(interaction_dict_from_pdb_final)})

	print(f"count_no_uniprot_map, {count_no_uniprot_map}")
	print(f"count_not_in_uniprot_seq_dict,{count_not_in_uniprot_seq_dict}")

	with open('out7_final_pairwise_interaction_dict','wb') as f:
		pickle.dump(interaction_dict_from_pdb_final,f,protocol=0)
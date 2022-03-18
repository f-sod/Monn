#!/usr/bin/python3
# ============================= for retrival uniprot sequences ===============
# first get pdbbind_all_mapping.tab from uniprot using the pdbid list from the above step
def get_uniprotid():
	uniprotid_set = set()
	with open('./pdbbind_index/INDEX_general_PL_name.2020') as f:
		for line in f.readlines():
			if line[0] != '#':
				lines = line.strip().split('  ')
				if lines[2] != '------': #some uniprot id are unknown
					uniprotid_set.add(lines[2])
	print(f"uniprotid_set step1,{len(uniprotid_set)}")
	
	with open('out1.4_pdb_uniprot_mapping.tab') as f:
		for line in f.readlines()[1:]:
			lines = line.split('\t')
			uniprotid_set.add(lines[1])
	print(f"uniprotid_set step2,{len(uniprotid_set)}")
	return uniprotid_set

uniprotid_set = get_uniprotid()
uniprotid_list = list(uniprotid_set)
with open('out1.5_uniprotid_list.txt','w') as fw:
	for uniprotid in uniprotid_list:
		fw.write(uniprotid+'\n')


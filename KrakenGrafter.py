## the purpose of this script is to graft on the contents of a .fasta file to the nodes.dmp and names.dmp used by Kraken2
	## in effect "grafting" on a new branch to the tree of life inside the nodes.dmp and names.dmp files
	## this allows the user to then use the "kraken2-build --add-to-library" functionality during database building
	## thereby making sure their own custom sequences are included

## the script takes in three input files:
	## nodes.dmp to modify
	## names.dmp to modify
	## .fasta of sequences to graft
		## NOTE: the seqIDs need to contain no whitespace and no pipe characters ("|")
		## the entire seqID will be grafted into the nodes.dmp and names.dmp files so maybe keep it concise

## the script then generates three output files:
	## new_nodes.dmp (or specified name) of the modified nodes.dmp
	## new_names.dmp (or specified name) of the modified names.dmp
	## K2.fasta (or specified name) of the modified .fasta file

## the script can be used in two different ways -
## in both cases, brand new taxon IDs are generated that have not been seen in the nodes.dmp and names.dmp files:
	## 1) a given .fasta file's sequence(s) is/are inserted directly beneath (taxonomically speaking) a specified "root" node
		##
		## example:
		##
		## KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta input.fasta -root 32630
		##
	## 2) the user declares a brand new "parent" node that is inserted beneath the "root" node (at a user-specified taxonomic depth)
	## and then the given .fasta sequence(s) is/are inserted directly beneath this new "parent" node
	## this functionality is enabled by declaring both -parent_taxon and -parent_rank
		##
		## example:
		##
		## KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta input.fasta -root 32630 -parent_taxon new_genus -parent_rank genus
		##

## NOTE: i_nodes, i_names, root, o_nodes, o_names, o_fasta, and debug are all optional variables to declare
	## the default "root" node is 32630 - which is the NCBI's "synthetic construct" node - a safe habour for new sequences
	## otherwise, one can either go on the NCBI taxonomy browser or use grep 'taxon' nopes.dmp to try and work out what taxon ID to use for "root"

## written by INZ - 04/21/22
## Stanford Unversity
## provided with no acceptance of liability or promise of functionality
## version 0.1.0

## ===============================================================================================================================================

## import libraries:

import copy
import argparse

def main(input_nodes_name, input_names_name, input_fasta_name, root_taxID, new_parent_taxon, new_parent_rank, output_nodes_name, output_names_name, output_fasta_name, debug):
	
	## input variables:
	
	##input_nodes_name = 'nodes.dmp'
	##input_names_name = 'names.dmp'
	##input_fasta_name = 'HSVd.fasta'
	##root_taxID = '147262'
	##new_parent_taxon = ''
	##new_parent_rank = ''
	##output_nodes_name = 'new_nodes.dmp'
	##output_names_name = 'new_names.dmp'
	##output_fasta_name = 'HSVd_K2.fasta'
	##debug = 1

	version = '0.1.0'
	
	## necessary fixed strings/lists:
	
	blank_node_line = 'taxid\t|\tparent\t|\trank\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|'
	blank_name_line = 'taxid\t|\tvar\t|\t\t|\tscientific name\t|'
	kraken_seqID_portion = '|kraken:taxid|replace  '
	linnean_taxonomy_list = ['kingdom','phylum','class','order','family','genus','species','subspecies']
	last_linnean_rank = linnean_taxonomy_list[-1]
	
	## check that the ranks make sense
	
	if new_parent_taxon:
		if not new_parent_rank:
			print("error: because a new parent taxon has been specified, a taxonomic rank for it needs to be chosen")
			print("please choose a rank above " + last_linnean_rank)
			print("quitting")
			quit()
	##
	
	if (new_parent_rank == last_linnean_rank):
		print("error: the new parent taxonomic rank cannot be " + last_linnean_rank)
		print("please choose a higher rank")
		print("quitting")
		quit()
	##
	
	## open files:
	
	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry
		##
		return open_list
	##
	
	input_nodes_list = opener(input_nodes_name)
	input_names_list = opener(input_names_name)
	input_fasta_list = opener(input_fasta_name)
	
	## convert fasta to singleline (based on https://stackoverflow.com/a/50856787):
	
	def singleline(in_fasta_list):
		out_fasta_list = []
		new_seq_line = []
		for line in in_fasta_list:
			if line.startswith('>'):
				if new_seq_line:
					out_fasta_list.append(''.join(new_seq_line))
					new_seq_line = []
				out_fasta_list.append(line)
			else:
				new_seq_line.append(line.strip())
			##
		if new_seq_line:
			out_fasta_list.append(''.join(new_seq_line))
		##
		return out_fasta_list
	##
	
	input_fasta_list = singleline(input_fasta_list)
	
	## determine the rank of the sequences to be inserted:
		## one lower than either the root_taxID or the new_parent_rank
	
	if not new_parent_taxon:
		for node_ind, node in enumerate(input_nodes_list):
			node_taxid = node.split('\t')[0]
			if (node_taxid == root_taxID):
				root_rank = node.split('\t')[4]
				break
		if (root_rank == last_linnean_rank):
			print("error: selected root node is too low in taxonomic rank (no available ranks beneath)")
			print("root rank: " + root_rank)
			print("quitting")
			quit()
		else:
			input_seq_rank = linnean_taxonomy_list[(linnean_taxonomy_list.index(root_rank)+1)]
			if debug:
				print("new sequences will be inserted at the rank: " + input_seq_rank)
	else:
		input_seq_rank = linnean_taxonomy_list[(linnean_taxonomy_list.index(new_parent_rank)+1)]
		if debug:
			print("new sequences will be inserted at the rank: " + input_seq_rank)
	##
	
	## rename the input .fasta seqIDs to have the kraken2 apropriate names with unique taxIDs
		## make sure to reserve the first unique taxID for a new parent taxon if specified
	
	if new_parent_taxon:
		new_parent_taxID = int(input_nodes_list[-1].split('\t')[0]) + 1
		first_child_taxID = new_parent_taxID + 1
	else:
		first_child_taxID = int(input_nodes_list[-1].split('\t')[0]) + 1
	##
	
	def renamer(in_fasta_list, init_taxID):
		out_fasta_list = []
		out_taxid_list = []
		out_seqID_list = []
		for line_ind, line in enumerate(in_fasta_list):
			if line.startswith('>'):
				new_seqID = line + kraken_seqID_portion + line[1:]
				new_seqID = new_seqID.replace('replace', str(init_taxID))
				out_fasta_list.append(new_seqID)
				out_fasta_list.append(in_fasta_list[line_ind+1])
				out_taxid_list.append(init_taxID)
				out_seqID_list.append(line[1:])
				init_taxID += 1
		##
		return out_fasta_list, out_taxid_list, out_seqID_list
	##
	
	renamed_fasta_list, appended_taxID_list, appended_seqID_list = renamer(input_fasta_list, first_child_taxID)
	
	## graft on the nodes to the input_nodes_list
	
	def nodegrafter(in_node_list, in_taxID_list, in_parent_taxID, in_child_rank):
		out_nodes_list = copy.deepcopy(in_node_list)
		for in_child_taxID in in_taxID_list:
			new_node_line = blank_node_line.replace('taxid', str(in_child_taxID))
			new_node_line = new_node_line.replace('parent', str(in_parent_taxID))
			new_node_line = new_node_line.replace('rank', str(in_child_rank))
			out_nodes_list.append(new_node_line)
		##
		return out_nodes_list
	##
	
	if new_parent_taxon:
		output_nodes_list = nodegrafter(input_nodes_list, [new_parent_taxID], root_taxID, new_parent_rank)
		output_nodes_list = nodegrafter(output_nodes_list, appended_taxID_list, new_parent_taxID, input_seq_rank)
	else:
		output_nodes_list = nodegrafter(input_nodes_list, appended_taxID_list, root_taxID, input_seq_rank)
	##
	
	## graft on the names to the input_names_list
	
	def namegrafter(in_name_list, in_taxID_list, in_seqID_list):
		out_names_list = copy.deepcopy(in_name_list)
		for in_child_taxID_ind, in_child_taxID in enumerate(in_taxID_list):
			new_name_line = blank_name_line.replace('taxid', str(in_child_taxID))
			new_name_line = new_name_line.replace('var', str(in_seqID_list[in_child_taxID_ind]))
			out_names_list.append(new_name_line)
		##
		return out_names_list
	##
	
	if new_parent_taxon:
		output_names_list = namegrafter(input_names_list, [new_parent_taxID], [new_parent_taxon])
		output_names_list = namegrafter(output_names_list, appended_taxID_list, appended_seqID_list)
	else:
		output_names_list = namegrafter(input_names_list, appended_taxID_list, appended_seqID_list)
	##
	
	## report some stats if debug is enabled
	
	if debug:
		print("KrakenGrafter.py version: " + version)
		print("added " + input_fasta_name + " to nodes.dmp / names.dmp")
		if new_parent_taxon:
			print("expecting 1 (new parent taxon) + " + str(len(appended_seqID_list)) + " = " + str(1+len(appended_seqID_list)) + " new lines added")
		else:
			print("expecting: " + str(len(appended_seqID_list)) + " new lines added")
		##
		print("added: " + str(len(output_nodes_list)-len(input_nodes_list)) + " new nodes")
		print("added: " + str(len(output_names_list)-len(input_names_list)) + " new names")
		print("saving")
	##	
	
	## save the new files
	
	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##
	
	saver(output_nodes_name, output_nodes_list)
	saver(output_names_name, output_names_list)
	saver(output_fasta_name, renamed_fasta_list)
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-i_nodes', type=str, default = 'nodes.dmp', help='input nodes.dmp file - default = "nodes.dmp"')
	parser.add_argument('-i_names', type=str, default = 'names.dmp', help='input names.dmp file - default = "names.dmp"')
	parser.add_argument('-i_fasta', type=str, help='input .fasta file - no whitespace in seqIDs!')
	parser.add_argument('-root', type=str, default = '32630', help='taxID of the root node under which new sequences will be grafted - default = 32630 ("synthetic construct")')
	parser.add_argument('-parent_taxon', type=str, default = '', help='optionally declare the name of new node to be grafted beneath the root node')
	parser.add_argument('-parent_rank', type=str, default = '', help='if a new node is to be grafted, a taxonomic rank for it must be declared - cannot be "subspecies"!')
	parser.add_argument('-o_nodes', type=str, default = 'new_nodes.dmp', help='output nodes.dmp filename - default = "new_nodes.dmp"')
	parser.add_argument('-o_names', type=str, default = 'new_names.dmp', help='output names.dmp filename - default = "new_names.dmp"')
	parser.add_argument('-o_fasta', type=str, default = 'K2.fasta', help='output .fasta filename - default = "K2.fasta"')
	parser.add_argument('-debug', type=int, default = 1, help='print more detailed results - default = 1')
	
	args = parser.parse_args()
	main(args.i_nodes, args.i_names, args.i_fasta, args.root, args.parent_taxon, args.parent_rank, args.o_nodes, args.o_names, args.o_fasta, args.debug)
##
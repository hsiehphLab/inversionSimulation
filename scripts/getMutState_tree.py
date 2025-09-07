import sys, os
import pandas as pd
from ete3 import Tree
from Bio import Phylo
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer

#sys.argv = ["","chr8-7301025-INV-5297356.treefile", "mapping_hap_INV.txt", "temp"]
#sys.argv = ["","12001024-12051024.treefile", "mapping_hap_INV.txt", "temp"]
#sys.argv = ["","8201024-8301024.treefile", "mapping_hap_INV.txt", "temp"]
#sys.argv = ["","output/chr2-152052016-152053604-IGC/iqtree/chr2-152052016-152053604-IGC.treefile", "CMP_CMP_Clint_1,CMP_CMP_Clint_2", "output/chr2-152052016-152053604-IGC/data/mapping_hap_SV.txt", "output/chr2-152052016-152053604-IGC/data/chr2-152052016-152053604-IGC.masked.filter.fa.hapIDs", "tmp"]

#sys.argv = ["","output/chr1-25199846-25212045-IGC/iqtree/chr1-25199846-25212045-IGC.treefile", "CMP_CMP_Clint_1,CMP_CMP_Clint_2", "output/chr1-25199846-25212045-IGC/data/mapping_hap_SV.txt", "output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.masked.filter.fa.hapIDs", "tmp"]

#sys.argv = ["", "output/chr6-161870034-161890814-IGC/iqtree/chr6-161870034-161890814-IGC.treefile","CMP_CMP_Clint_1,CMP_CMP_Clint_2","output/chr6-161870034-161890814-IGC/data/mapping_hap_SV.txt","output/chr6-161870034-161890814-IGC/data/chr6-161870034-161890814-IGC.masked.filter.fa.hapIDs","output/chr6-161870034-161890814-IGC/iqtree/chr6-161870034-161890814-IGC.treefile.wMut"]

fname_tree = sys.argv[1]
outG = sys.argv[2]
tree = Phylo.read(fname_tree, "newick")
tips_tree = [n.name for n in tree.get_terminals()]
for g in outG.split(","):
	if g in tips_tree:
		tree.collapse(g)

fname_mapping_hap_SV = sys.argv[3]
fout_prefix = sys.argv[4]
df_mapping_hap_SV = pd.read_csv(fname_mapping_hap_SV, delimiter="\t", na_values=['.'], index_col="orig_hapID")


list_traits = [1 if x==True else 0 for x in df_mapping_hap_SV.loc[:,"SV"]]
list_haps = [x for x in df_mapping_hap_SV.index]
list_trait_aln = []
list_hapTrait = zip(list_haps, list_traits)
for pair in list_hapTrait:
	if pair[1] == 0:
		list_trait_aln.append(SeqRecord(Seq("A"), id=pair[0]))
	elif pair[1] == 1:
		list_trait_aln.append(SeqRecord(Seq("T"), id=pair[0]))

trait_aln = MultipleSeqAlignment(list_trait_aln)
terms = tree.get_terminals()
terms.sort(key=lambda term: term.name)
trait_aln.sort()
column_i = trait_aln[:, 0]
column_i == len(column_i) * column_i[0]
clade_states = dict(zip(terms, [{c} for c in column_i]))

score_i = 0
count = 0

for clade in tree.get_nonterminals(order="postorder"):
	addmut = False
	clade_childs = clade.clades
	try:
		left_state = clade_states[clade_childs[0]]
		right_state = clade_states[clade_childs[1]]
	except IndexError:
		continue
	state = left_state & right_state
#	print(state, left_state, right_state)
	if not state:
		state = left_state | right_state
		score_i = score_i + 1
		addmut = True
	if addmut:
		clade_states[clade] = state
		clade.name =  "Inner%s_" % count + "".join(list(state)) + ":"
		count += 1
	else:
		clade_states[clade] = state
		clade.name =  "Inner%s_" % count + "".join(list(state)) + "notMut:"
		count += 1

print ("NumberMutationEvents:%s" % score_i)

Phylo.write(tree,"%s.nwk" % fout_prefix  ,"newick", format_branch_length='%1.8f')

#Phylo.write(tree,"temp.xml","phyloxml")
#Phylo.convert("temp.nwk", "newick", "temp.nex","nexus")
#Phylo.convert("temp.nwk", "newick", "temp.xml","phyloxml")





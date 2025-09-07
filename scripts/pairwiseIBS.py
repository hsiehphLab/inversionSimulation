import numpy as np
import pandas as pd
import sys, os
from sklearn.metrics.pairwise import pairwise_distances
import itertools, math


if __name__ == "__main__":

	tped = pd.read_csv(sys.argv[1], header=None, delimiter=" ", na_values=["."])
	tfam = pd.read_csv(sys.argv[2], header=None, delimiter=" ", na_values=["."])
	# read IGC genotypes
	IGC_mapping = pd.read_csv(sys.argv[3], delimiter="\t", na_values=["."], index_col="orig_hapID")
	fout = sys.argv[4]
	
	# make a list of hapIDs in the order as they appear in the tfam file
	l_hapIDs_tfam = [ x+"_"+str(y) for x in tfam.iloc[:,1] for y in range(1,3)]
	
	# build a matrix for SNV genotypes and reindex the matrix by l_hapIDs_tfam
	pd_haplo_matrix = tped.iloc[:,4:].transpose().reset_index(drop=True)
	pd_haplo_matrix.index = l_hapIDs_tfam

	# reorder the dataframe IGC_mapping as the tfam file (l_hapIDs_tfam)
	IGC_mapping = IGC_mapping.reindex(l_hapIDs_tfam)

	# find rows/hapIDs that have NaN in the SV column of IGC_mapping
	missingG_hapIDs = IGC_mapping.loc[pd.isna(IGC_mapping["SV"]), :].index

	# remove rows of IGC_mapping and pd_haplo_matrix listed in missingG_hapIDs
	pd_haplo_matrix.drop(missingG_hapIDs, inplace=True)

	
	IGC_mapping.drop(missingG_hapIDs, inplace=True)

	# compute element-wise difference for pairs of rows in pd_haplo_matrix
	l_comb = []
	for c in itertools.combinations(pd_haplo_matrix.index,2): 
		try:
			np_pair_haplo_matrix = pd_haplo_matrix.reindex(c).to_numpy()

			# make a similarity array btw the two haplotypes
			array_similarity = np_pair_haplo_matrix[0,] == np_pair_haplo_matrix[1,]

			# set up a masked array w/ the shaple of the expected similarity array
			dim = array_similarity.shape
			array_mask = np.full((dim),0)

			for i in range(2):
				array_mask[np.where(np.isnan(np_pair_haplo_matrix[i,]))] = 1

			new_array_similarity = np.ma.masked_array(array_similarity, mask=array_mask)
			IBS = np.nansum(new_array_similarity)/(array_mask.shape[0] - sum(array_mask==1))
			l_comb.append((c[0], c[1], IBS))
			# store the similarity array
			ll = np.array([0 if x == True else 1 for x in array_similarity])
			l_pair_Diff = [math.nan if v == 1 else ll[i] for i,v in enumerate(array_mask)]
			
			l_pair_hapID = [y for y in IGC_mapping.reindex(c).index.tolist()]
			if np.nan in l_pair_hapID:
				print (l_pair_Diff)
				print (l_pair_hapID)
				print (c)
				exit("Error!")
		except KeyError:
			print (c)
			exit("Error!")
		try:
			l_out = [int(x) for x in IGC_mapping.loc[l_pair_hapID,"SV"].values.tolist()]
		except ValueError:
			print ("#####")
			print( inv.loc[l_pair_hapID,"SV"].values.tolist())
			print (l_pair_Diff)
			print (l_pair_hapID)
			print (c)
			print( inv.loc["EUR_CEU_NA12878_hap1","SV"])
			exit("Error!")
			
		l_out.extend(l_pair_hapID)
		l_out.extend(l_pair_Diff)
		sys.stdout.write("\t".join([str(x) for x in l_out]) + "\n")
	

	# compute pairwise IBS
	# matrix to pairwise
	df_pw_IBS = pd.DataFrame(l_comb)
	df_pw_IBS.columns = ["h1","h2","IBS"]

	df_pw_IBS = pd.merge(df_pw_IBS, IGC_mapping.loc[:,["SV"]], left_on="h1", right_index=True)
	df_pw_IBS = pd.merge(df_pw_IBS, IGC_mapping.loc[:,["SV"]], left_on="h2", right_index=True)
	

	df_pw_IBS.columns = ["h1","h2","IBS","h1_SVgeno","h2_SVgeno"]

	df_pw_IBS.to_csv(fout, sep="\t", index=False)




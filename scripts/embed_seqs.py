"""Code to embed a set of sequences in an embedding space using a trained protvec model and an embedding set of sequences. Creates a .csv file of the embedded sequences. Also returns file of sequences which could not be successfully embedded""" 

import pickle 
import numpy as np
from Bio import SeqIO
import pandas as pd
from sklearn.preprocessing import minmax_scale
import warnings
from random import shuffle	
warnings.filterwarnings(action='once')

def seq_3mers(sequence): 
	"""Takes a sequence to overlapping 3-mers"""
	
	seq_size = len(sequence) 
	seq_3mers = list() #intialise list 
	
	#iterate through sequence to obtain 3mers 
	for i in range (1,seq_size-1):
		seq_3mers.append(sequence[i-1]+sequence[i]+sequence[i+1])
		
	return seq_3mers

def seq_vector(seq_3mers, vec_dict, d, seq): 
	"""Converts a list of 3-mers to a vector, vec_dict is the model vectors and d is the number of dimensions in the embedding"""
	
	#intialise empty array 
	seq_arr = np.zeros((len(seq_3mers), 100))
	
	#populate the array
	for i in range (0,len(seq_arr)): 
		if seq_3mers[i] in vec_dict: 
			seq_arr[i] = vec_dict.get(seq_3mers[i])
		else: 
			
			print('\n3-mer not in model!!!!')
			print('BAD 3MER: '+seq_3mers[i])
			print(seq)
			
	#take the array to a vector 
	sequence_vec = pd.DataFrame(seq_arr, index = seq_3mers).sum(axis = 0, skipna = True)
	
	return sequence_vec

def normalise_vec(embedded_seqs, d): 
	"""Normalises embedded sequences in d dimensions. Got this from Phage-Prot https://github.com/mikejhuang/PhageProtVec/blob/master/protvec.py"""
	
	features = []
	for i in range(d):
		tempvec = [vec[i] for vec in embedded_seqs]
		print(tempvec)
		mean = np.mean(tempvec)
		var = np.var(tempvec)
		features.append([(vec[i]-mean)/var for vec in embedded_seqs])
	features = np.array(features).reshape(len(embedded_seqs),len(features))
	return features

def standardise_vec(embedded_seqs,d): 
	"""Standardisation for embedded sequences in d dimensions"""
	
	#intialise an empty matrix 
	stand_embed = np.zeros((len(embedded_seqs),d))
	
	#normalise each vector 
	for i in range(0,len(embedded_seqs)):
		x = embedded_seqs[i]
		z = (x-np.mean(x))/np.std(x)
		stand_embed[i] = z
	
	return stand_embed 

def missing_3mers(three_mers, vec_dict): 
	"""Checks that the 3mers are in the model"""
	
	missing = [x for x in three_mers if x not in list(vec_dict.keys())]
	
	return missing 
	
def embed_sequences(sequences,vec_dict, d): 
	"""Embeds sequences in d dimensions using some protvec model vec_dict"""
	
	#intialise an empty matrix for sequence vectors 
	embedded_seqs = np.zeros((len(sequences), d))

	#intialise an empty matrix for the id of sequences embedded 
	embedded_seqs_keys = list()
 
	#intialise a list of dictionaries for sequences with 3mers which are not in the model 
	bad_seqs = list()
	
	#get the keys of the sequences 
	seq_keys = list(sequences.keys()) #list of the sequence ids 

	#reshuffle the sequence keys 
	shuffle(seq_keys) 

	#iterate through the sequence to obtain vectors
	for i in range(0, len(sequences)):	

		#get 3mers in the sequence 
		three_mers = seq_3mers(str(sequences.get(seq_keys[i]).seq))
		
		#check the 3-mers are in the model (if not exclude these sequences)
		missing = missing_3mers(three_mers, vec_dict)
		
		if len(missing) == 0 and len(sequences.get(seq_keys[i]).seq) <= 1024:
 
			#get the vector of the sequence
			vec = seq_vector(three_mers, vec_dict, d, sequences.get(seq_keys[i]).seq)
			#add vector to array 
			embedded_seqs[i] = vec
			#add corresponding sequence key to the list
			embedded_seqs_keys.append(seq_keys[i])

		else: 
			bad_seqs.append({'seq_ID': seq_keys[i], 'sequence':str(sequences.get(seq_keys[i]).seq), 'bad_3mers': missing})
   
	#drop rows with zeroes (where there was a sequence could not be embedded) 
	embedded_seqs = embedded_seqs[~np.all(embedded_seqs == 0, axis=1)]
 
	#standardise the vectors 
	standard_embedding = standardise_vec(embedded_seqs,d)
  
	return standard_embedding, embedded_seqs_keys, bad_seqs

#import the pickled vector dictionary 
with open('../protvec_models/bacillus_3mervectors.pkl', 'rb') as f: 
	vec_dict = pickle.load(f) 

#load in real sequences to embed (see what the speed is like when I run it locally)
seqs = SeqIO.index("../sequences/bacillus_embeddingset.fa", 'fasta')

#embed the sequences  
print('EMBEDDING SEQUENCES')
embedding,seqs_keys, missing = embed_sequences(seqs, vec_dict, 100)  
print('sequences embedded') 

#Assemble the embedding in a dataframe and drop the empty rows 
embedding_df = pd.DataFrame(embedding)
embedding_df.index = seqs_keys

#save the embedding and missing 3mers so we can evaluate them locally 
print('Saving the embedding') 
embedding_df.to_csv('../embedded_sequences/bacillus_filtered_embedded.csv', sep = '\t')
(pd.DataFrame(missing)).to_csv('../embedded_sequences/bacillus_filtered_notembedded.csv', sep = '\t')
print('Embedding saved!') 




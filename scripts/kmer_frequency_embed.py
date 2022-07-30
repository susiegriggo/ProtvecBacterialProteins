""" 
Script which takes sequences amino acid sequences based on the frequency of 3-mers 
Requires biopython version 1.77 as Alphabet was deprecated in vesion 1.78 (September 2020)
""" 

#imports 
import numpy as np 
import pandas as pd
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import Reduced
import itertools
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from collections import Counter
import seaborn as sns
from sklearn.utils.extmath import randomized_svd 
from scipy.spatial import distance 
import random

def seq_3mers(sequence): 
    """Takes a sequence to overlapping 3-mers"""
    
    seq_size = len(sequence) 
    seq_3mers = list() #intialise list 
    
    #iterate through sequence to obtain 3mers 
    for i in range (1,seq_size-1):
        seq_3mers.append(sequence[i-1]+sequence[i]+sequence[i+1])
        
    return seq_3mers

def murphy10(seq_str): 
    """Takes an amino acid sequence using the standard 20 amino acid code to reduced 10 letter alphabet. 
        This funcseqs_swiss_keystion requires biopython version 1.77 or lower. Input is a a string of amino acids"""
    
    #turn starting sequence into a sequence object 
    intial_seq = Seq(seq_str, Alphabet.ProteinAlphabet())
    
    #intialise sequence object 
    new_seq = Seq('', Alphabet.Reduced.Murphy10())
    
    #iterate through the letters in the sequence and convert to murphy10 
    for aa in intial_seq: 
        new_seq += Alphabet.Reduced.murphy_10_tab[aa]
    
    return str(new_seq)

def seq_vector(seq, embedding): 
    """Embeds a sequence as a kmer frequency embedding"""
    
    #break the seq into kmers 
    seq_kmers = seq_3mers(seq)
    
    #intialise a vector for the sequence to be embedded 
    seq_vec = np.zeros(len(embedding))
    
    #iterate through the kmers in the sequence 
    for kmer in seq_kmers: 
        
        #add the kmer vector to make the sequence vector 
        seq_vec += embedding[kmer]
        
    #divide the sequence by the number of kmers (number of kmer counts) (NOT SURE IF THIS IS CORRECT - PLAY AROUND WITH this)
    seq_vec = seq_vec/len(seq_kmers)
    
    return seq_vec  

def embedkmers_seqs(seqs, embedding): 
    """Embed a list of sequences with a dataframe of kmer frequency"""
    
    #intialise an array to hold the embeddings
    embed_kmerfreq = np.zeros((len(seqs), len(embedding)))
    
    #iterate through the sequences 
    for i in range(len(seqs)): 
        #get the sequence 
        seq = seqs[i]
        
        #get the vector
        seq_vec= seq_vector(seq, embedding)
        
        #add the sequnce vector to the embeddings matrix 
        embed_kmerfreq[i] = seq_vec 
        
    return embed_kmerfreq 
    
#import the embedding sequences 
embed_seqs_dict = SeqIO.index("../Annotating_Bacterial_Function_Space_DATA/bacillus_carbohydratemetabolism_embeddingset.fa", 'fasta')
embed_seqs_keys = list(embed_seqs_dict.keys()) #gives md5 hashes of the sequences 
embed_seqs = [str(embed_seqs_dict.get(key).seq) for key in embed_seqs_keys]

#get a random subset of 16763 sequences to embed (this was the number of sequences embedded for bacteroides)
randint = random.sample(range(len(embed_seqs)), 16763) 
embed_seqs_keys = embed_seqs_keys

#determine which sequences contain the invalid character 'X' and remove them from the set of sequences to embed 
embed_seqs_containsX = ['X' not in seqs for seqs in embed_seqs]
keys_containsX = [embed_seqs_keys[i] for i in range(len(embed_seqs_keys)) if embed_seqs_containsX[i] == True]
embed_seqs = [str(embed_seqs_dict.get(key).seq) for key in keys_containsX]
embed_seqs_keys = keys_containsX 

#generate a list of all possible kmeres for the murphy10 alphabet 
murphy10_sub = Alphabet.Reduced.murphy_10_tab 
murphy10_l = set([d[1] for d in list(murphy10_sub.items())]) #list of letters in the murphy10 alphabet  
k = 3 #intialise the length of the kmer 
kmers = [''.join(kmer) for kmer in list(itertools.product(murphy10_l, repeat = k))]

#intialise idnetity matrix size of kmers to represent the kmer embedding (each 1 is denotes a different kmer)
kmerfreq = np.identity(len(kmers))
#represent as a dataframe 
kmerfreq_df = pd.DataFrame(kmerfreq) 
kmerfreq_df.columns = kmers
kmerfreq_df.index = kmers

#convert the embedded sequences to murphy 10 
embed_seqs_murphy10 = [murphy10(seq) for seq in embed_seqs]
#embed the sequences 
embed_kmerfreq = embedkmers_seqs(embed_seqs_murphy10, kmerfreq_df)

#save the embedding 
embed_kmerfreqDf = pd.DataFrame(embed_kmerfreq, index = embed_seqs_keys)
embed_kmerfreqDf.columns = kmers
embed_kmerfreqDf.to_csv('../Annotating_Bacterial_Function_Space_DATA/bacillus_carbohydratemetabolism_kmer_frequency.csv')

""""
Train Protvec model using Word2vec algorithm. 
Takes a filtered training set of sequences and returns a pickled dictionary of vectors for each possible 3-mer in the training set.
Code modified from biovec https://github.com/kyu999/biovec
""" 

#imports
from Bio import SeqIO
from gensim.models import word2vec
import numpy as np 
import sys
import gzip 
import os
import pickle
from pyfasta import Fasta
from tqdm import tqdm  


def split_kmers(seq, k):
    """Splits a sequence into overlapping kmers"""
    a, b, c = zip(*[iter(seq)]*k), zip(*[iter(seq[1:])]*k), zip(*[iter(seq[2:])]*k)
    str_kmers = []
    for kmers in [a,b,c]:
        x = []
        for kmer in kmers:
            x.append("".join(kmer))
        str_kmers.append(x)
    return str_kmers


def load_protvec(model_fname): 
    """Loads an already made protvec model"""
    return word2vec.Word2Vec.load(model_fname)


def generate_corpusfile(fasta_fname, n, corpus_fname):
    """
    Args:
        fasta_fname: corpus file name
        n: the number of chunks to split. In other words, "n" for "n-gram"
        corpus_fname: corpus_fnameput corpus file path
    Description:
        Protvec uses word2vec inside, and it requires to load corpus file
        to generate corpus.
    """
    f = open(corpus_fname, "w")
    fasta = Fasta(fasta_fname)
    for record_id in tqdm(fasta.keys(), desc='corpus generation progress'):
        r = fasta[record_id]
        seq = str(r)
        ngram_patterns = split_kmers(seq, n)
        for ngram_pattern in ngram_patterns:
            f.write(" ".join(ngram_pattern) + "\n")
    f.close()


class ProtVec(word2vec.Word2Vec):
    """
    Either fame or corpus is required.
	corpus_fname: fasta file for corpus
    corpus: corpus object implemented by gensim
    k: k of kmer
    out: corpus output file path
    min_count: least appearance count in corpus. if the n-gram appear k times which is below min_count, the model does not remember the n-gram
    """
    def __init__(self,k, size, min_count=1, corpus_fname=None, corpus=None,
                 out="corpus.txt",  sg=1, window=25, workers=64):

        self.k= k
        self.size = size 
        self.corpus_fname = corpus_fname


        if corpus is None and corpus_fname is None:
            raise Exception("Either corpus_fname or corpus is needed!")

        if corpus_fname is not None:
            print ('Generate Corpus file from fasta file')
            generate_corpusfile(corpus_fname, k, out)
            corpus = word2vec.Text8Corpus(out)
            print ("\n... OK\n")

        #create the model
        print('creating the model')
        word2vec.Word2Vec.__init__(self, corpus, size=self.size, sg=sg, window=window, min_count=min_count, workers=workers)
        print('model created!') 

    def pickleDictionary(self, dict_fname):
        """Gets the vectors in the model and takes them to a pickled dictionary"""
        vectors = self.wv
  
        
        print('Creating a dictionary of vectors')
        
        vec_dict = ({})
        for idx, key in enumerate(vectors.vocab): 
            vec_dict[key] = vectors[key]
            
        print('pickling the dictionary')
        pickle.dump(vec_dict, open(dict_fname, 'wb'))

#set variables for the embeding  - these could be changed
k = 3 #kmer length 
dim = 100 #number of dimensions 
min_count = 1
print('reading in files')
c_fname ="../Annotating_Bacterial_Function_Space_DATA/bacillus_carbohydratemetabolism_trainingset.fa"  #fasta file of sequences used to train the model
dict_fname ='../Annotating_Bacterial_Function_Space_DATA/bacillus_carbohydratemetabolism_3mervectors.pkl'  #3mer dictionary name 
kmer_corpus ="../Annotating_Bacterial_Function_Space_DATA/bacillus_carbohydratemetabolism_corps.txt" #output kmer corpus 
print('files read!') 

pv = ProtVec(k, dim, min_count, corpus_fname=c_fname, 
                 out=kmer_corpus)
pv.pickleDictionary(dict_fname)

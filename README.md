# Embedding Bacterial Proteins with Protvec 

Here we demonstrate the use of Protvec sequence embeddings for evaluating bacterial protein annotations. 

This repository contains the code to train and embed sequences using Protvec models with _Bacillus_ carbohydrate metabolism sequences and unannotated _Bacillus_ sequences. We also supply code to analyse the resulting models and embeddings.  

Training and embedding sequenes for the _Bacillus_ carbohydrate metabolism model are contained in this repository. Sequences to train the unknkown Protvec model can be downloaded [here](https://cloudstor.aarnet.edu.au/plus/remote.php/webdav/bacillus_unknown_trainingset.fa.gz). 

## Usage 
- Training Protvec models requires the python [gensim](https://pypi.org/project/gensim/) module 
- Clone the repository and unzip the sequences 
```
git clone https://github.com/susiegriggo/ProtvecBacterialProteins.git
gunzip data/*  
```

## Protvec Models 

- Protvec models are trained (using code adapted from [biovec](https://github.com/kyu999/biovec)) with the Protvec algorithm [1]. A vector size of 100 and a context size of 25 is used. We provide the training script `train_protvec.py`. 
- The Protvec model trained on _Bacillus_ carbohydrate metabolism seqeuences was compared with BLOSUM62 matrix. This analysis is shown in notebooks, `BLOSUM_comparison.ipynb`
- Alternatively trained Protvec models can be downloaded: 
  - Protvec model trained with [_Bacillus_ carbohydrate metabolism sequences](https://doi.org/10.25451/flinders.19770379)  
  - Protvec model trained with [unannotated _Bacillus_ sequences](https://doi.org/10.25451/flinders.19770742)  

### Protvec Sequence Embedding 
- Sequences are embedded using the provided `embed_seqs.py` script.
- In this study sequences were embedded using the _Bacillus_ carbohydrate metabolism Protvec model as a [Protvec model trained with 324,018 sequences from the Swiss-Prot database](http://dx.doi.org/10.7910/DVN/JMFHTN) [1]. 

### K-mer frequency Embedding 
- Sequences were also embedded as vectors based on the frequeny of each 3-mer using the scipt `kmer_frequency_embed.py`. 

### Visualisation 
- Sequence embeddings are visualised using PCA in the notebook `vis_embeddings.ipynb`. 

## Cluster Prediction 
- The number of clusters present within the embeddings is estimated using the Calinski-Harabasz index. The procedure for identifying clusters is included in the R script `get_CHcriterion.R`. 
- Calinski Harabasz index for different embedding techniques is compared in the notebook `num_clusters.ipynb`  

### Comapring Protvec Clustering with Subsystems 

- The clusters of sequences which arise from the Protvec sequence embedding are compared with their subsystems annotations by building a tanglegram using the R  dendextend package [2] as shown in the notebook `tanglegram.ipynb` 

## Unknown Protein Sequences 

- We trained a Protvec model with _Bacillus_ sequences with unknown function to embed and cluster _Bacillus_ sequences with unknown function. Clusters were formed using k-means and visualised using _t_-SNE `cluster_unknowns.ipynb`


## References 

[1] Asgari, E. and Mofrad, M.R., 2015. Continuous distributed representation of biological sequences for deep proteomics and genomics. PloS one, 10(11), p.e0141287.
[2] Galili, T., 2015. dendextend: an R package for visualizing, adjusting and comparing trees of hierarchical clustering. Bioinformatics, 31(22), pp.3718-3720.

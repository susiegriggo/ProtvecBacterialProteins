# Embedding Bacterial Proteins with Protvec 

**NEW:** This work has been published in BMC Bioinformatics. You can see the publication [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04930-5). 

Here we demonstrate the use of Protvec sequence embeddings for evaluating bacterial protein annotations. 

This repository contains the code to train and embed sequences using Protvec models with _Bacillus_ carbohydrate metabolism sequences and unannotated _Bacillus_ sequences. We also supply code to analyse the resulting models and embeddings.  

All training sequences, embedding sequences, and annotations can be downloaded [here](https://cloudstor.aarnet.edu.au/plus/s/wxRIVLKBejzAutc/download). 

## Usage 

Clone the repository and add sequences and data to the data directory. 

```
git clone https://github.com/susiegriggo/ProtvecBacterialProteins.git
cd ProtvecBacterialProteins/data 
wget https://cloudstor.aarnet.edu.au/plus/s/wxRIVLKBejzAutc/download
unzip download 
gunzip Annotating_Bacterial_Function_Space/*
mv Annotating_Bacterial_Function_Space/* ./
```

## Protvec Models 
- Training Protvec models requires the [gensim](https://pypi.org/project/gensim/) module 
- Protvec models are trained (using code adapted from [biovec](https://github.com/kyu999/biovec)) with the Protvec algorithm [1]. A vector size of 100 and a context size of 25 is used. We provide the training script `train_protvec.py`. This script can be modified for different training data sets. 
- The Protvec model trained on _Bacillus_ carbohydrate metabolism seqeuences was compared with BLOSUM62 matrix. This analysis is shown in the notebook, `BLOSUM_comparison.ipynb`
- Alternatively, pre-trained Protvec models can be downloaded: 
  - Protvec model trained with [_Bacillus_ carbohydrate metabolism sequences](https://doi.org/10.25451/flinders.19770379)  
  - Protvec model trained with [unannotated _Bacillus_ sequences](https://doi.org/10.25451/flinders.19770742)  

### Protvec Sequence Embedding 
- Sequences are embedded using the provided `embed_seqs.py` script. This script can be modified for different training data sets and Protvec models. 
- In this study sequences were embedded using the _Bacillus_ carbohydrate metabolism Protvec model as a [Protvec model trained with 324,018 sequences from the Swiss-Prot database](http://dx.doi.org/10.7910/DVN/JMFHTN) [1]. 

### K-mer frequency Embedding 
- Sequences were also embedded as vectors based on the frequeny of each 3-mer using the scipt `kmer_frequency_embed.py`. 

### Visualisation 
- Sequence embeddings are visualised using PCA in the notebook `vis_embeddings.ipynb`. 

## Cluster Prediction 
- The number of clusters present within the embeddings is estimated using the Calinski-Harabasz index. The procedure for identifying clusters is included in the R script `get_CHcriterion.R`. This script can be modified for different sets of embedded sequences. 
- Calinski Harabasz index for different embedding techniques is compared in the notebook `num_clusters.ipynb`  

### Comapring Protvec Clustering with Subsystems 

- The clusters of sequences which arise from the Protvec sequence embedding are compared with their subsystems annotations by building a tanglegram using the R  dendextend package [2] as shown in the notebook `tanglegram.ipynb` 

## Unknown Protein Sequences 

- We trained a Protvec model with _Bacillus_ sequences with unknown function to embed and cluster _Bacillus_ sequences with unknown function. Clusters were formed using k-means and visualised using _t_-SNE `cluster_unknowns.ipynb`
- The sequence similarity of the 100 sequences closest to the centroid for each cluster is compared in the notebook `unknown_cluster_sequence_similarity.ipynb` 

## Citation 

- If you use this work please cite it! 

Grigson, S.R., McKerral, J.C., Mitchell, J.G. and Edwards, R.A., 2022. Organizing the bacterial annotation space with amino acid sequence embeddings. BMC bioinformatics, 23(1), pp.1-14.

## References 

[1] Asgari, E. and Mofrad, M.R., 2015. Continuous distributed representation of biological sequences for deep proteomics and genomics. PloS one, 10(11), p.e0141287.
[2] Galili, T., 2015. dendextend: an R package for visualizing, adjusting and comparing trees of hierarchical clustering. Bioinformatics, 31(22), pp.3718-3720.

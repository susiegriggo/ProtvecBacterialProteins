# ProtvecBacterialAnnotations

Using Protvec sequence embeddings [1] for evaluating bacterial protein sequence annotations. 

## Protvec Models 

- Protvec models are trained using code adapted from biovec https://github.com/kyu999/biovec using the Protvec algorithm [1] with a vector size of 100 and a context size of 25. We provide the training script `train_protvec.py`. 

- The Protvec model trained on _Bacillus_ carbohydrate metabolism seqeuences was compared with BLOSUM62 matrix. This analysis is shown in notebooks, `BLOSUM_comparison.ipynb`

- As well as the models trained in this work, sequences are embedded using a [Protvec model trained with 324,018 sequences from the Swiss-Prot database](http://dx.doi.org/10.7910/DVN/JMFHTN). 

- Trained Protvec models are available at #TODO - add figshare/roads link - include all and the Bacillus model 

## Sequence Embedding 

### K-mer frequency Embedding 

Sequences are embedded as vectors based on the frequeny of _3_-mers using the scipt `kmer_frequency_embed.py`. 

### Protvec Embedding 

A set of 25,000 _Bacillus_ carbohydrate metabolism sequences are embedded using the _Bacillus_ Protvec model and the Swissprot Protvec model using `embed_seqs.py`. 

### Visualisation 
Sequence embeddings are visualised using PCA in the notebook `vis_embeddings.ipynb`. 


## Cluster Prediction 
The number of clusters present within the embeddings is estimated using the Calinski-Harabasz index. The procedure for identifying clusters is included in the R script `get_CHcriterion.R` 

The predicted number of clusters is visualised in the notebook `num_clusters.ipynb`  


### Comapring Protvec Clustering with Subsystems 

The clusters of sequences which arise from the Protvec sequence embedding are compared with their subsystems annotations by building a tanglegram using the R  dendextend package [2] as shown in the notebook `tanglegram.ipynb` 

## Unknown Protein Sequences 

We trained a Protvec model with _Bacillus_ sequences with unknown function to embed and cluster _Bacillus_ sequences with unknown function. Clusters were formed using k-means and visualised using _t_-SNE. #TODO


## Citing this work 
Please cite this work if you use it!
#TODO add citation - do I need a zendodo link or something 

## References 

[1] Asgari, E. and Mofrad, M.R., 2015. Continuous distributed representation of biological sequences for deep proteomics and genomics. PloS one, 10(11), p.e0141287.
[2] Galili, T., 2015. dendextend: an R package for visualizing, adjusting and comparing trees of hierarchical clustering. Bioinformatics, 31(22), pp.3718-3720.

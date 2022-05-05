# ProtvecBacterialAnnotations

Using Protvec sequence embeddings [1] for evaluating bacterial protein sequence annotations. 

### Protvec Models 

- Protvec models are trained using code adapted from biovec https://github.com/kyu999/biovec using the Protvec algorithm [1] with a vector size of 100 and a context size of 25. We provide the training script `train_protvec.py`. 

- The Protvec model trained on _Bacillus carbohydrate metabolism_ seqeuences was compared with BLOSUM62 matrix. This analysis is shown in notebooks, 

- As well as the models trained in this work, sequences are embedded using a [Protvec model trained with 324,018 sequences from the Swiss-Prot database](http://dx.doi.org/10.7910/DVN/JMFHTN). 

### Protvec Embedding 

Visualisation of Protvec embeddings using PCA included in notebooks 

### K-mer frequency embedding 

Visualisation of K-mer frequency is included in 

### Cluster Prediction 
The number of clusters present within the embeddings is estimated using the Calinski-Harabasz index. R script 

Visualisation of the Calinski Harabasz index is included in 

### Comapring Protvec Clustering with Subsystems 

Using tanglegram. Example is provided in notebooks

### Unknown Protein Sequences 
Trained model. Visualisation is included in notebooks 

## Citing this work 
Please cite this work if you use it!
#TODO add citation

## Dependecies



## References 

[1] Asgari, Ehsaneddin, 2015, "Replication Data for: Continuous Distributed Representation of Biological Sequences for Deep Proteomics and Genomics", https://doi.org/10.7910/DVN/JMFHTN

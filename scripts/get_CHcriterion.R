"""
Script to generate the Calinski Harabasz Index and sum of square errors
Completed with 500 bootstraps 
Works on a  subsample of 5000 sequences

Written in R to use R cluster library functions hclust, agnes and diana
""" 

#imports 
library(dplyr)
library(cluster)
library(tidyr)
library(ggplot2)
library(parallelDist)

Distance <- function(cluster)
{
    # the center of the cluster, mean of all the points
    center <- colMeans(cluster)

    # calculate the summed squared error between every point and 
    # the center of that cluster 
    distance <- apply( cluster, 1, function(row)
    {
        sum( ( row - center )^2 )
    }) %>% sum()

    return(distance)
}

WSS <- function( data, groups )
{
    k <- max(groups)

    # loop through each groups (clusters) and obtain its 
    # within sum squared error 
    total <- lapply( 1:k, function(k)
    {
    # extract the data point within the cluster
        cluster <- subset( data, groups == k )

        distance <- Distance(cluster)
        return(distance)
    }) %>% unlist()

    return( sum(total) )
}


CHCriterion <- function( data, kmax, clustermethod, metric, ...  )
{
    if( !clustermethod %in% c("hclust", 'agnes', 'diana' ) )
        stop( "method must be one of 'hclust', 'agnes' or 'diana'" )

    # total sum squared error (independent with the number of cluster k)
    tss <- Distance( cluster = data )

    
    # initialize a numeric vector storing the score
    wss <- numeric(kmax)
    
  
    #get the distance matrix 
     d <- parDist( data, method = metric )
   
    
    # k starts from 2, cluster 1 is meaningless
    start_time <- Sys.time()
    if( clustermethod == "agnes" )
    {
        clustering <- agnes( d, diss = TRUE, keep.dis = FALSE, keep.data = FALSE, ... )
    
    }else if( clustermethod == "diana" )
    {
        clustering <- diana( d, diss = TRUE, keep.dis = FALSE, keep.data = FALSE )

    }else # "hclust"
    {
        clustering <- hclust( d, ... )
    }
    
    for( k in 2:kmax )
        {
            groups <- cutree( clustering, k )
            wss[k] <- WSS( data = data, groups =  groups )
        }
    end_time <- Sys.time()
    print(end_time - start_time)
    
    # between sum of square
    bss <- tss - wss[-1]

    # cluster count start from 2! 
    
    numerator <- bss / ( 1:(kmax-1) )
    denominator <- wss[-1] / ( nrow(data) - 2:kmax )
    
        
    criteria <- cbind(numerator / denominator, wss[-1])
    
}


#read in bacillus data embedded with bacillus model
bacil_embedding <- read.delim(file = 'bacillus_carbohydratemetabolism_embedded.tsv', sep = '\t')
#set the sequence md5s as the index
rownames(bacil_embedding) <- bacil_embedding$X
bacil_embedding$X <- NULL
#standardise the embedding
bacil_stand <- scale(bacil_embedding)

#drop any rows (sequences) containing nan values
bacil_stand <- bacil_stand[complete.cases(bacil_stand), ]

#generate a random sample of 5,000 sequences for each embedding
#find the number of embedded bacillus sequences
samples <- dim(bacil_stand)[1]
#get indexes of the random sequences
subset <- sample.int(samples-1,5000)
#get a sample of the embedded sequences
bacil_sample <- bacil_stand[subset,] #sample of 5,000 sequence


#create loop for getting these metrics with bootstrapping 
B <- 500  #number of bootstraps 
k <- 150  #maximum number of clusters to consider
n <- nrow(bacil_sample)

#intialise matrices to hold data generated 
ch_data <- matrix(ncol = (k-1), nrow = B)
wss_data <- matrix(ncol = (k-1), nrow = B)


for (b in 1:B) {
  cat("Iteration ", b, "\n") #state which bootstrap up to 
    
  #generate a bootstrap sample 
  bb <- bacil_sample[sample(1:n, n, replace=TRUE), ]
    
  #generate indicies 
  tmp <- CHCriterion(data = bb, kmax = k, clustermethod = "hclust", metric = 'euclidean',method = 'ward.D')
    
  #add to output matrixward
  ch_data[b,] <- tmp[,1]
  wss_data[b,] <- tmp[,2]
    }

#save and give labels to the output matrices 
i <- 1:B
b_labels <- paste0('B',i)
j <- 2:k
k_labels <- paste0('k', j)

#create dataframes
ch_data <- data.frame(ch_data)
colnames(ch_data) <- k_labels
rownames(ch_data) <- b_labels

wss_data <- data.frame(wss_data)
colnames(wss_data) <- k_labels
rownames(wss_data) <- b_labels

#save output 
f <- file('calinski_bacil.csv', "wb")
write.csv(ch_data, file=f, eol="\n")
close(f)
f <- file('wss_bacil.csv', "wb")
write.csv(wss_data, file=f, eol="\n")
close(f)



library(RMySQL)
library(car)
library(RColorBrewer)
library(cluster)
library(tm)
library(slam)

colMax <- function(data){
  #http://stackoverflow.com/questions/24212739/how-to-find-the-highest-value-of-a-column-in-a-data-frame-in-r
  apply(data, 2, max, na.rm = TRUE)
}

colSort <- function(data, ...) sapply(data, sort, ...)

dtm.filter <- function(dtm, filter.start, filter.end){
  if(is.null(filter)){
    cat("Dimension of filtered matrix: ", dim(dtm))
    return(dtm)
  }else{
    wc<-data.frame(colMax(dtm))
    allowed_words_idx <- which(wc[,1]>=min(c(filter.start, filter.end)) & wc[,1]<= max(c(filter.start, filter.end)))
    cat("Dimension of filtered matrix: ", dim(dtm[,allowed_words_idx]))
    return(dtm[,allowed_words_idx])
  }
}

dtm.matrix.tfidf.cluster.kmeans <- function(dtm.matrix.tfidf, DTM, count.cluster){
  #http://michael.hahsler.net/SMU/CSE7337/install/tm.R
  
  ### k-means (this uses euclidean distance)
  #dtm.matrix.tfidf <- dtm.matrix.tfidf[rowSums(dtm.matrix.tfidf)!=0, ]
  rownames(dtm.matrix.tfidf) <- 1:nrow(dtm.matrix.tfidf)
  
  ### don't forget to normalize the vectors so Euclidean makes sense
  norm_eucl <- function(m){
    norm.m <- apply(m, MARGIN=1, FUN=function(x) sum(x^2)^.5)
    norm.m[norm.m==0] = 1
    m/norm.m
  }
  dtm.matrix.tfidf.norm <- norm_eucl(dtm.matrix.tfidf)
  
  ### cluster into 10 clusters
  dtm.matrix.tfidf.cluster <- kmeans(dtm.matrix.tfidf.norm, count.cluster)
  print("Cluster distribution:")
  print(table(dtm.matrix.tfidf.cluster$cluster))
  
  ### show clusters using the first 2 principal components
  plot(prcomp(dtm.matrix.tfidf.norm)$x, col=dtm.matrix.tfidf.cluster$cluster)
  
  dtm.matrix.tfidf.cluster$terms.frequent <- rep("", length(dtm.matrix.tfidf.cluster$cluster))
  for(cluster.i in 1:count.cluster){
    docs.cluster.index <- which(dtm.matrix.tfidf.cluster$cluster==cluster.i)
    terms.frequency <- sort(col_sums(DTM[docs.cluster.index,]), decreasing = T)
    dtm.matrix.tfidf.cluster$terms.frequent[docs.cluster.index] <- paste(names(terms.frequency)[1:min(5,length(terms.frequency))], collapse=",")
    print(paste0(c("Cluster: ", cluster.i, " | Terms: ", dtm.matrix.tfidf.cluster$terms.frequent[docs.cluster.index[1]] , "\n"), collapse=""))
  }
  
  return(dtm.matrix.tfidf.cluster)
}

rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

tf <- function(dtm.matrix){
  count.docs <- nrow(dtm.matrix)
  count.terms <- ncol(dtm.matrix)

  tf <- dtm.matrix %/% rep.col(rowSums(dtm.matrix), count.terms)
  tf[is.nan(tf)] <- 0
  
  return(tf)
}

idf <- function(dtm.matrix){
  count.docs <- nrow(dtm.matrix)
  count.terms <- ncol(dtm.matrix)
  
  dtm.matrix.flag <- as.matrix(dtm.matrix>=1)
  idf <- log(rep(count.docs, count.terms) %/% colSums(dtm.matrix.flag))
  idf <- rep.row(idf, count.docs)
  colnames(idf) <- colnames(dtm.matrix)
  rownames(idf) <- rownames(dtm.matrix)
  
  return(idf)
}

tfidf.image.write <- function(m, file){
  png(file)
  #par(mar = rep(0, 4))
  image(m, axes = FALSE, col = grey(seq(0, 1, length = 256)))
  dev.off()
}


tfidf.important.terms <- function(dtm.matrix.tfidf, start, end){
  index <- data.frame(which(dtm.matrix.tfidf>=start & dtm.matrix.tfidf<=end, arr.ind = T))
  return(colnames(dtm.matrix.tfidf)[unique(index$Terms)])
}

tfidf.important.terms.range <- function(dtm.matrix.tfidf, dtm.matrix){
  tfidf.important.terms.df = data.frame(term = c(), start = c(), end = c(), word.count = c())
  dtm.matrix.term.freq = data.frame(colSums(dtm.matrix))
  step = 0.05
  
  for(start in seq(0,max(dtm.matrix.tfidf), step)){
    terms <- tfidf.important.terms(dtm.matrix.tfidf, start = start, end = start + step)
    
    if(length(terms)>0){
      df <- data.frame(term=terms, start = start, end = start + step, word.count = dtm.matrix.term.freq[match(terms, colnames(dtm.matrix)),1])
      tfidf.important.terms.df = rbind(tfidf.important.terms.df, df)
      print(df)
    }else{
      print(paste0("No terms found | Start:", start, " | End: ", start + step))
    }
  }
  
  return(tfidf.important.terms.df)
}

tfidf.important.docs <- function(dtm.matrix.tfidf, doc_list, start, end){
  index <- data.frame(which(dtm.matrix.tfidf>=start & dtm.matrix.tfidf<=end, arr.ind = T))
  #return(rownames(dtm.matrix.tfidf)[unique(index$Docs)])
  results <- doc_list[sort(unique(index$Docs))]
  names(results) <- sort(unique(index$Docs))
  return(results)
}

tfidf.dist.hcluster <- function(dtm.matrix.cluster){
  tfidf.dist <- dist(dtm.matrix.cluster)
}

pca.visualize.3d <- function(dtm.filter.result, groups){
  #remove zero vector points
  dtm.filter.result.clean <- dtm.filter.result[rowSums(dtm.filter.result) != 0, ]
  groups.clean <- groups[rowSums(dtm.filter.result) != 0]
  
  cat("Cleaned zero dimensions: ", dim(dtm.filter.result), " to ", dim(dtm.filter.result.clean), " | groups.clean: ", length(groups.clean))
  
  dtm.matrix.pca <- prcomp(dtm.filter.result.clean, center = TRUE, scale. = TRUE)
  raw <- dtm.matrix.pca$x
  colors <- colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(groups.clean)))
  car::scatter3d(raw[,1], raw[,2], raw[,3],surface=F,groups = as.factor(groups.clean),
                 point.col = colors[as.numeric(as.factor(groups.clean))], 
                 surface.col = colors[as.numeric(as.factor(groups.clean))])
  # use below to see the color palette
  #pie(rep(1, length(groups.clean)),labels = "", border =  "transparent", col = colors[as.numeric(as.factor(groups.clean))])
}

cluster.best.fit <- function(dtm.filter.result){
  wss <- (nrow(dtm.filter.result)-1)*sum(apply(dtm.filter.result,2,var))
  
  for (i in 2:30){
    wss[i] <- sum(kmeans(dtm.filter.result, centers=i)$withinss)
  }
  
  plot(1:30, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
  return(wss)
}

mydb = dbConnect(MySQL(), user='xxxxx', password='xxxxx', dbname='JMLR', host='xxxxx.amazonaws.com')
rs = dbSendQuery(mydb, "select * from volume")
volume = fetch(rs, n=-1)

View(volume)

rs = dbSendQuery(mydb, "select * from abstract")
abstract = fetch(rs, n=-1)
abstract <- abstract[nchar(trimws(abstract$abstract))!=0, ]

View(abstract)

myCorpus <- Corpus(VectorSource(as.character(abstract$abstract)))
getTransformations()

myCorpus <- tm_map(myCorpus, content_transformer(tolower))
myCorpusCopy <- myCorpus

myCorpus <- tm_map(myCorpus, removePunctuation)
myCorpus <- tm_map(myCorpus, removeNumbers)
myStopwords <- c(stopwords('english'), "available", "via", "can", "learning", "algorithms", "algorithm",
                 "problems", "data", "problem", "model", "method")
myCorpus <- tm_map(myCorpus, removeWords, myStopwords)
#myCorpus <- tm_map(myCorpus, stemDocument)
#myCorpus <- tm_map(myCorpus, stemCompletion, dictionary=myCorpusCopy)

DTM <- DocumentTermMatrix(myCorpus, control = list(minWordLength = 1,  sparse=TRUE))
TDM <- TermDocumentMatrix(myCorpus, control = list(minWordLength = 1))

dtm.matrix <- as.matrix(DocumentTermMatrix(myCorpus, control = list(minWordLength = 1)))
dtm.matrix.tfidf <- as.matrix(DocumentTermMatrix(myCorpus, control = list(weighting = weightTfIdf)))

#write.csv(dtm.matrix.tfidf, 'dtm.matrix.tfidf.csv')

count.docs <- nrow(dtm.matrix)
count.terms <- ncol(dtm.matrix)

# Get important terms and documents
terms.important.range <- tfidf.important.terms.range(dtm.matrix.tfidf, dtm.matrix)
docs.important <-tfidf.important.docs(dtm.matrix.tfidf, abstract$abstract, start = 0.8, end = 1)
docs.length <- data.frame(length=unlist(lapply(abstract$abstract,nchar)))

#PCA Visualize
dtm.filter.result <- dtm.filter(dtm.matrix.tfidf, 0.3, 0.35)
pca.visualize.3d(dtm.filter.result, abstract$volume_id)

#Clustering
dtm.filter.result <- dtm.filter(dtm.matrix.tfidf, 0.2, 10)
cluster.best.fit.series <- cluster.best.fit(dtm.filter.result)

count.cluster = 20

dtm.cluster.result <- dtm.matrix.tfidf.cluster.kmeans(dtm.filter.result, DTM,count.cluster = count.cluster)
abstract$cluster <- dtm.cluster.result$cluster
abstract$terms.frequent <- dtm.cluster.result$terms.frequent
#save output------------------------------------------------------------------------
write.csv(abstract, 'result.abstract.csv')
write.csv(terms.important.range, 'terms.important.range.csv')

#print image
tfidf <- tf(dtm.matrix) * idf(dtm.matrix)
tfidf.image <- floor(tfidf*(255/max(tfidf)))
tfidf.image.write(tfidf, "tfidf.png")

#print image
dtm.matrix.tfidf.image <- floor(dtm.matrix.tfidf*(255/max(dtm.matrix.tfidf)))
tfidf.image.write(dtm.matrix.tfidf.image, "dtm.matrix.tfidf.png")











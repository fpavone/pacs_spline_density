rm(list=ls())
data <- read.csv("topsecret/Dati_PSD.csv", sep=";")
data <- data[,-c(1:3)]

delete <- NULL
count <- 1

for(i in 1:nrow(data)){
  for(j in 1:ncol(data)){
    if(data[i,j] > 100) data[i,j] <- 100
    if(j>1) data[i,j] <- data[i,j] - sum(data[i,c(1:j-1)])
    if(data[i,j] < 0){
      delete[count] <- i
      count <- count + 1
    }
  }
}

if(count > 1){
  data <- data[-delete,]
}

names <- colnames(data)
for(it in 1:length(names)) names[it] <- substring(names[it],3,nchar(names[it]))
colnames(data) <- names
  
write.table(data, file="data.txt", row.names = FALSE, col.names=TRUE, quote = FALSE)









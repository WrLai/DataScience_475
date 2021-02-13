data<-read.csv("dataSpikeAllSequences.csv")

string_similarity<-function(a,b){
  n=nchar(a)
  m=nchar(b)
  S=matrix(0,n+1,m+1)
  S[1,]=c(0:m)
  S[,1]=c(0:n)
  a=substring(a,1:n,1:n)
  b=substring(b,1:m,1:m)
  for (i in 2:(n+1)){
    for (j in 2:(m+1)){
      if (a[i-1]==b[j-1]){
        temp=0
      }
      else{
        temp=1
      }
      S[i,j]=min(S[i-1,j]+1,S[i,j-1]+1,S[i-1,j-1]+temp)
    }
  }
  return (1-S[n+1,m+1]/max(n,m))
}


library(dplyr)

# China
data_China<-data%>%
  filter(Location=="China") %>%
  arrange(collection.Date) %>%
  group_by(collection.Date) %>%
  slice(1)
data_China<-data_China[-29,]
data_China$collection.Date<-gsub("'","",data_China$collection.Date)
data_China$collection.Date<-as.Date(data_China$collection.Date,format="%Y-%m-%d")

similarity=c()
sequences=data_China$sequences
for (i in 1:(length(sequences)-1)){
  similarity=c(similarity,string_similarity(sequences[i],sequences[i+1]))
}
similarity_China<-data.frame(date=data_China$collection.Date[2:28],similarity=similarity)
library(ggplot2)
ggplot(similarity_China,aes(x=date,y=similarity))+geom_line()+
  labs(x="date",y="similarity",title="similarity of RNA sequences in China")+
  scale_x_date(date_labels="%Y-%m-%d",breaks="10 day")

library(utils)
similarity=c()
sequences=data_China$sequences
for (i in 1:(length(sequences)-1)){
  similarity=c(similarity,adist(sequences[i],sequences[i+1])[1])
}
similarity1_China<-data.frame(date=data_China$collection.Date[2:28],similarity=similarity)
library(ggplot2)
ggplot(similarity1_China,aes(x=date,y=similarity))+geom_line()+
  labs(x="date",y="similarity",title="distance of RNA sequences in China (using adist)")+
  scale_x_date(date_labels="%Y-%m-%d",breaks="10 day")

library(stringdist)
similarity=c()
sequences=data_China$sequences
for (i in 1:(length(sequences)-1)){
  similarity=c(similarity,stringdist(sequences[i],sequences[i+1]))
}
similarity2_China<-data.frame(date=data_China$collection.Date[2:28],similarity=similarity)
library(ggplot2)
ggplot(similarity2_China,aes(x=date,y=similarity))+geom_line()+
  labs(x="date",y="similarity",title="distance of RNA sequences in China (using stringdist)")+
  scale_x_date(date_labels="%Y-%m-%d",breaks="10 day")

# USA
data_USA<-data%>%
  filter(Location=="USA") %>%
  arrange(collection.Date) %>%
  group_by(collection.Date) %>%
  slice(1)
data_USA<-data_USA[2:227,]
data_USA$collection.Date<-gsub("'","",data_USA$collection.Date)
data_USA$collection.Date<-as.Date(data_USA$collection.Date,format="%Y-%m-%d")

similarity=c()
sequences=data_USA$sequences
for (i in 1:(length(sequences)-1)){
  similarity=c(similarity,string_similarity(sequences[i],sequences[i+1]))
}
similarity_USA<-data.frame(date=data_USA$collection.Date[2:226],similarity=similarity)
library(ggplot2)
ggplot(similarity_USA,aes(x=date,y=similarity))+geom_line()+
  labs(x="date",y="similarity",title="similarity of RNA sequences in USA")+
  scale_x_date(date_labels="%Y-%m-%d",breaks="1 month")

#Austrilia
data_Austrilia<-data%>%
  filter(Location=="Australia") %>%
  arrange(collection.Date) %>%
  group_by(collection.Date) %>%
  slice(1)
data_Austrilia<-data_USA[2:128,]
data_Austrilia$collection.Date<-gsub("'","",data_Austrilia$collection.Date)
data_Austrilia$collection.Date<-as.Date(data_Austrilia$collection.Date,format="%Y-%m-%d")

similarity=c()
sequences=data_Austrilia$sequences
for (i in 1:(length(sequences)-1)){
  similarity=c(similarity,1-adist(sequences[i],sequences[i+1])[1]/max(nchar(sequences[i]),nchar(sequences[i+1])))
}
similarity_Austrilia<-data.frame(date=data_Austrilia$collection.Date[2:127],similarity=similarity)
library(ggplot2)
ggplot(similarity_Austrilia,aes(x=date,y=similarity))+geom_line()+
  labs(x="date",y="similarity",title="similarity of RNA sequences in Austrilia")+
  scale_x_date(date_labels="%Y-%m-%d",breaks="1 month")

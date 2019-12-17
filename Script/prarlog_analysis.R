####get the paralog annotation from "biomaRt"
Get_ensembl_paralog <- function(Ensemble_ID_DAPPL,outputfile="temp.csv")
{
  library("biomaRt")
  ensembl<- useMart("ensembl")
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

  #all Ensemble ID in DAPPL
  Group <- rep("",length(Ensemble_ID_DAPPL))
  index <- 1
  for(i in 1147:length(Ensemble_ID_DAPPL))
  {
    print(i)
    if(Group[i] == "")
    {
      t <- getBM(attributes = "hsapiens_paralog_ensembl_gene",filters="ensembl_gene_id",value=Ensemble_ID_DAPPL[i],mart = ensembl)
      if(nrow(t)==0)
      {
        Group[i] <- paste0("No_paralog")
      }else
      {
        DAPPL_in_list_ID<- which(Ensemble_ID_DAPPL %in% t[,1])
        if(length(DAPPL_in_list_ID) == 0)
        {
          Group[i] <- paste0("No_paralog_inDAPPL")
        }else
        {
          Group[i] <- paste0("Group_",index)
          Group[which(Ensemble_ID_DAPPL %in% t[,1])] <- paste0("Group_",index)
          index <- index+1
          print(Group[i] )
        }
      }
    }
  }
  Temp_Group <- cbind(Ensemble_ID_DAPPL,Group)
  write.csv(Temp_Group,file = outputfile)
}
##########main##################
library("readxl")
annotation_paralog <- read_xlsx("annotation/Paralog_TFs.xlsx",sheet = "Paralog_TFs")
#####here we only get the all file names from DAPPL,just get the file names
Files_name <- as.character(annotation_paralog$Files_name)
TFs <- unlist(lapply(Files_name, function(x){return(toupper(unlist(strsplit(x,split = "-"))[3]))}))
####get the groups,Ensemble ID and paralog
TF_ensembleID_Groups <- read_xlsx(path = "annotation/Paralog_TFs.xlsx",sheet = "Paralog_TFs_Ensamble")
TF_Groups <- cbind(annotation_paralog ,TF_ensembleID_Groups$Group[match(TFs,as.character(TF_ensembleID_Groups$Gene_DAPPL))])
####get the grounp and plateID to get the index#################
temp_plate <- unlist(lapply(TF_Groups[,1],function(x){
  temp_ID <- unlist(strsplit(x,split = "A-"))[1]
  temp_ID <- as.integer(unlist(strsplit(temp_ID,split = "p"))[2])
  return( min(8,ceiling(temp_ID/2)))
}))
temp <- cbind(temp,temp_plate)
###sequences count
Sum_seq_N20_N20me <- read.csv("~/Documents/DAPPL_review/annotation/Sum_seq_N20_N20me.csv")
temp <- cbind(temp,Sum_seq_N20_N20me[match(temp[,1],Sum_seq_N20_N20me[,1]),c(2,3)])
write.csv(temp,file = "annotation/Paralog_TFs_temp.csv",quote = F,row.names = F)
###select the dupilicate data
temp <- temp[temp[,2]!="No_paralog_inDAPPL" & temp[,2]!="No_paralog" & temp[,2]!="Unknown",]#1286
index <- paste0(temp[,3],"_",temp[,2])
temp <- temp[index != "NA_NA",]
index <- paste0(temp[,3],"_",temp[,2])
###select the dupilicate data
count <- unlist(tapply(index, as.factor(index), function(x){return(length(x))}))
count_index <- unlist(tapply(index, as.factor(index), function(x){return(x[1])}))
count_index <- count_index[count > 1]#[count >=5 ]#
temp <- temp[index %in% count_index, ]
temp <- temp[order(index),]
write.csv(temp,file = "annotation/Paralog_TFs_temp.csv",quote = F,row.names = F)
###picture
index <- paste0(temp[,3],"_",temp[,2])
temp <- temp[,c(1:5)]
temp <- data.frame(cbind(temp,index))
ratio_max_min <- tapply(as.integer((temp$Sum_seq_N20)), index, function(x){return(max(log(x,10))/min(log(x,10)))})
hist(ratio_max_min,xlab= "log10(Max)/log10(min)")
max_min <- tapply(as.integer((temp$Sum_seq_N20)), index, function(x){return(max(x)/min(x))})
hist(max_min,xlab= "Max/min")


#######picture#######
index <- paste0(temp[,3],"_",temp[,2])
temp <- temp[,c(1:5)]
temp <- data.frame(cbind(temp,index))
TFs <- unlist(lapply(as.character(temp[,1]),function(x){
  return(toupper(unlist(strsplit(x,split = "-"))[3]))
}))

TFs_duplicated <- TFs
TFs[which(duplicated(TFs))] <- paste0(TFs[which(duplicated(TFs))],"_2")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_3")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_4")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_5")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_6")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_7")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_8")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_9")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_10")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_11")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_12")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_13")
TFs[which(duplicated(TFs))] <- paste0(TFs_duplicated[which(duplicated(TFs))],"_14")
TFs <- factor(TFs,levels = TFs)
temp <- data.frame(cbind(temp,index,TFs))
temp <- data.frame(temp, x=seq(1:nrow(temp)))

library(ggplot2)#figure only for count >= 5 and here are only 46 groups the bigest group with count 20
P1 <- ggplot(temp[1:10,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))
P2 <- ggplot(temp[11:21,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
   geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
   theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))
P3 <- ggplot(temp[22:28,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))

P4 <- ggplot(temp[29:43,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))

P5 <- ggplot(temp[44:56,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))

P6 <- ggplot(temp[57:70,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))

P7 <- ggplot(temp[71:91,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))

P8 <- ggplot(temp[92:115,], aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",angle = 45,hjust =1))+guides(fill=guide_legend(title ="Paralog group"))

ggplot(temp, aes(x=TFs, y=log(as.numeric(Sum_seq_N20),10), fill=index)) +
  geom_bar(stat="identity") + ylab("Count of sequences (log10)")+xlab("")+
  theme_bw()+theme(axis.text = element_text(color = "black",hjust =1))+guides(fill=guide_legend(title ="Paralog group"))+
  coord_flip()+scale_fill_manual(values = rep(c("cyan", "gold1", "darkseagreen","brown","darkolivegreen","burlywood","coral","aquamarine",
    "cyan1", "gold", "darkseagreen1","brown1","darkolivegreen1","burlywood1","coral1","aquamarine1",
    "cyan2", "gold2", "darkseagreen2","brown2","darkolivegreen2","burlywood2","coral2","aquamarine2",
    "cyan3", "gold3", "darkseagreen3","brown3","darkolivegreen3","burlywood3","coral3","aquamarine3",
    "cyan4", "gold4", "darkseagreen4","brown4","darkolivegreen4","burlywood4","aquamarine4",
    "darkslategray2",  "darkslategray1","darkorange","bisque","darkorange1","coral4","darkorange2")))+
  theme(panel.grid=element_line(colour=NA),panel.border = element_blank(),axis.line = element_line(color = "black"))+ggsave(filename = "paralog.pdf",height = 35,width = 10)

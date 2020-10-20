N16.alignment.matching <- function(fastq.filename=NULL, alignment.filename=NULL)
{
  library("Biostrings");
  read.seq <- readLines(fastq.filename);
  read.seq <- read.seq[seq(2,length(read.seq),4)]
  # matching matrix
  mat <- matrix(-1,5,5);
  mat[1,1] <- 1;
  mat[2,2] <- 1;
  mat[3,3] <- 1;
  mat[4,4] <- 1;
  mat[5,5] <- 1;
  colnames(mat) <- c("A", "C", "G", "T", "N");
  rownames(mat) <- c("A", "C", "G", "T", "N");
  # left sequence "GGGAGAAGGTCATCAAGAGG" local alignment
  left.seq <- "GGGAGAAGGTCATCAAGAGG";
  left.align <- pairwiseAlignment(read.seq, left.seq, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -1);
  left.align.start <- start(pattern(left.align));
  left.align.end <- end(pattern(left.align));	
  left.align.nchar <- nchar(left.align);
  left.align.score <- score(left.align);
  left.align.last4letter <- substr(as.character(pattern(left.align)), left.align.nchar-3, left.align.nchar);
  
  cat("left alignment Done!\n");
  # right sequence "CCTACGT" local alignment
  right.seq <- "CCTACGT";
  sub.read.seq <- substr(read.seq, 81, nchar(read.seq));
  right.align <- pairwiseAlignment(sub.read.seq, right.seq, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -1);
  right.align.start <- 80+start(pattern(right.align));
  right.align.end <- 80+end(pattern(right.align));	
  right.align.nchar <- nchar(right.align);
  right.align.score <- score(right.align);
  right.align.first4letter <- substr(as.character(pattern(right.align)), 1, 4);
  cat("right alignment Done!\n");

  # middle sequence "GGCTAGCAGCCACTATAAGCTTCGAAGACTGG" local alignment
  middle.seq <- "GGCTAGCAGCCACTATAAGCTTCGAAGACTGG";
  middle.align <- pairwiseAlignment(read.seq, middle.seq, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -0.2);
  middle.align.start <- start(pattern(middle.align));
  middle.align.end <- end(pattern(middle.align));	
  middle.align.nchar <- nchar(middle.align);
  middle.align.score <- score(middle.align)
  middle.align.first5letter <- substr(as.character(pattern(middle.align)), 1, 5);
  middle.align.last4letter <- substr(as.character(pattern(middle.align)), middle.align.nchar-3, middle.align.nchar);
  cat("middle alignment Done!\n");
 
  # middle.align.first5letter,  
  middle.seq_2 <- "GCAGCCACTATAAGCTT";
  middle.align_2 <- pairwiseAlignment(read.seq, middle.seq_2, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -0.2);
  middle.align.start_2 <- start(pattern(middle.align_2));
  middle.align.end_2 <- end(pattern(middle.align_2));	
  middle.align.nchar_2 <- nchar(middle.align_2);
  middle.align.score_2 <- score(middle.align_2);
  middle.align.first4letter <- substr(as.character(pattern(middle.align_2)), 1, 4);
  middle.align.before2letter <- substr(read.seq, middle.align.start_2-2, middle.align.start_2-1);
  
  
  DNA.align <- data.frame(left.align.score, left.align.nchar, left.align.start, left.align.end, left.align.last4letter, right.align.score, right.align.nchar, right.align.start, right.align.end, right.align.first4letter, middle.align.score, middle.align.nchar, middle.align.start, 
                          middle.align.end, middle.align.first5letter, middle.align.last4letter,middle.align.first4letter,middle.align.before2letter);
  write.csv(DNA.align, file=alignment.filename, row.names=FALSE);
}
call.N16.alignment.matching <- function()
{
  All_split_files <- list.files(paste0("Raw_Data/Split/"));
  fastq.file <- All_split_files[grep("A_D",All_split_files)]
  for(i in 1:length(fastq.file))
  {
    fastq.filename <- paste0("Raw_Data/Split/", fastq.file[i]);
    alignment.filename <-paste0("Raw_Data/Alignment/", fastq.file[i], ".N16.alignment");
    N16.alignment.matching(fastq.filename=fastq.filename, alignment.filename=alignment.filename);
      
   cat(paste(fastq.filename, "\tDone!", "\n"));
  }  
}
# call.N16.alignment.matching()

N16.matching.score.statistics <- function(index="A_D")
{
  fastq.file <- list.files(paste0("Raw_data/Split/"));
  fastq.file <- fastq.file[grep(index,fastq.file)]
  
  # mapping score statistics
  matching.sta <- NULL;
  for(i in 1:20)
  {
    alignment.filename <- paste0("Raw_Data/Alignment/", fastq.file[i], ".N16.alignment");
    alignment <- read.csv(alignment.filename);
    matching.sta <- rbind(matching.sta, rbind(alignment[,c(1,6,11)]));
    cat(alignment.filename, "\tDone!", "\n");
  }
  # left matching score histogram
  hist.left <- hist(matching.sta[,1], plot=FALSE, breaks=20);
  png(paste0("Raw_Data/matching.sta/",index,".N16.left.align.score.png"))
  plot(hist.left, main="left sequence alignment", ylab="read number", xlab="matching score", col="lightblue", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5); 
  dev.off();
  # right matching score histogram
  hist.right <- hist(matching.sta[,2], plot=FALSE, breaks=20);
  png(paste0("Raw_Data/matching.sta/",index,".N16.right.align.score.png"))
  plot(hist.right, main="right sequence alignment", ylab="read number", xlab="matching score", col="lightblue", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5); 
  dev.off();
  # middle matching score histogram
  hist.middle <- hist(matching.sta[,3], plot=FALSE, breaks=20);
  png(paste0("Raw_Data/matching.sta/",index,".N16.middle.align.score.png"))
  plot(hist.middle, main="middle sequence alignment", ylab="read number", xlab="matching score", col="lightgreen", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5);  
  dev.off();
}
# N16.matching.score.statistics("A_D")
# N16.matching.score.statistics("E_H")

library("Biostrings");
N16.protein.barcode.and.bseq.sequence<- function( index="A_D",annotation.file="annotation/annotation_ETS.csv",
                                                  Split_Protein_dir="Raw_Data/Split_Protein",control_lib="N16",plateID=c("OneCyc_rep1","OneCyc_rep2","OneCyc_rep3"),
                                                  plateID_indexL=c("AA","AC","T"),plateID_indexR=c("CT","CT","CT"))
{
  system(paste0("mkdir ",Split_Protein_dir))
  system(paste0("mkdir ",Split_Protein_dir,"/",control_lib,"/"))
  
  fastq.file <- list.files(paste0("Raw_data/Split/"));
  fastq.file <- fastq.file[grep(index,fastq.file)]
  
  library("Biostrings");
  TF.annotation<- read.csv(annotation.file,stringsAsFactors = FALSE)
  TF.annotation.Barcode <- gsub(" ", "",toupper(as.character(TF.annotation[,"Barcode"])))
  TF.file <- gsub(" ", "",toupper(as.character(TF.annotation[,4])))
  
  all.tf.barcode <- as.character(TF.annotation[,"Barcode"])
  RevCom <- c()
  for(i in 1:length(all.tf.barcode))
  {
    temp <- DNAString(x = as.character(all.tf.barcode[i]))
    temp <- reverseComplement(temp) 
    RevCom <- c(RevCom,as.character(temp))
  }
  all.tf.barcode <-RevCom
  
  total.ALL <- 0
  total.num <- c(0,0,0,0);
  left.num <- c(0,0,0,0);
  middle.num <- c(0,0,0,0);
  right.num <- c(0,0,0,0);
  bseq20.num <- c(0,0,0,0);
  A.num <- c(0,0,0,0);
  protein.num <- c(0,0,0,0);
  P12.protein.num <- c(0,0,0,0);
  
  for(i in 1:length(fastq.file))
  {
    fastq.filename <- paste0("Raw_data/Split/", fastq.file[i]);
    alignment.filename <- paste0("Raw_data/Alignment/", fastq.file[i], ".N16.alignment");
    
    
    read.seq <- readLines(fastq.filename);
    read.seq_ALL <- read.seq[seq(2,length(read.seq),4)]
    alignment_ALL <- read.csv(alignment.filename);
    total.ALL <- total.ALL + nrow(alignment_ALL)
    
    Library_index_L <- substr(read.seq_ALL, as.numeric(as.character(alignment_ALL[,"left.align.start"]))-2, as.numeric(as.character(alignment_ALL[,"left.align.start"]))-1);
    for(j in 1:length(plateID))
    {
      barcode.bseq.filename <- paste0("Raw_data/Barcode.bseq/",plateID[j],"_", fastq.file[i], ".N16.barcode.bseq");
      
      read.seq <- read.seq_ALL[Library_index_L==plateID_indexL[j]]
      alignment <- alignment_ALL[Library_index_L==plateID_indexL[j],]
      total.num[j] <- total.num[j]+nrow(alignment)
      
      read.left <- read.seq[alignment[,"left.align.last4letter"]=="GAGG" & alignment[,"left.align.score"]>=16];
      align.left <- alignment[alignment[,"left.align.last4letter"]=="GAGG" & alignment[,"left.align.score"]>=16,];
      left.num[j] <- left.num[j] + nrow(align.left);#GGCAT
      
      read.left.middle <- read.left[align.left[,"middle.align.first5letter"]=="GGCTA" & align.left[,"middle.align.last4letter"]=="CTGG" & align.left[,"middle.align.score"]>=20];
      align.left.middle <- align.left[align.left[,"middle.align.first5letter"]=="GGCTA" & align.left[,"middle.align.last4letter"]=="CTGG" & align.left[,"middle.align.score"]>=20, ];
      middle.num[j] <- middle.num[j] + nrow(align.left.middle);
      
      read.left.middle.right <- read.left.middle[align.left.middle[,"right.align.first4letter"]=="CCTA"];
      align.left.middle.right <- align.left.middle[align.left.middle[,"right.align.first4letter"]=="CCTA", ];
      right.num[j] <- right.num[j]+nrow(align.left.middle.right);
      
      read.left.middle.right.bseq <- read.left.middle.right[(align.left.middle.right[,"middle.align.start"]-align.left.middle.right[,"left.align.end"]-1)==16];
      align.left.middle.right.bseq <- align.left.middle.right[(align.left.middle.right[,"middle.align.start"]-align.left.middle.right[,"left.align.end"]-1)==16, ];
      bseq20.num[j] <- bseq20.num[j]+nrow(align.left.middle.right.bseq);
      barcode <- substr(read.left.middle.right.bseq, as.numeric(align.left.middle.right.bseq[,"middle.align.end"])+1, as.numeric(align.left.middle.right.bseq[,"right.align.start"])-1);

      bseq <- substr(read.left.middle.right.bseq, as.numeric(align.left.middle.right.bseq[,"left.align.end"])+2, as.numeric(align.left.middle.right.bseq[,"middle.align.start"])-1);
      UMI <- substr(barcode, 1, 8);
      A <- substr(barcode, 9, 9);
      barcode <- substr(barcode, 9, nchar(barcode));
      barcode.bseq <- data.frame(bseq, UMI, A, barcode);
      barcode.bseq <- barcode.bseq[barcode.bseq[,"A"]=="A" & nchar(as.character(barcode.bseq[,"UMI"]))==8, ];
      A.num[j] <- nrow(barcode.bseq)+A.num[j];
      barcode.bseq <- barcode.bseq[!is.na(match(barcode.bseq[,"barcode"], all.tf.barcode)), ];
      protein.num[j] <- nrow(barcode.bseq)+protein.num[j];
      #write.csv(barcode.bseq, file=barcode.bseq.filename, row.names=FALSE);
      P12.barcode.bseq <- barcode.bseq[!is.na(match(barcode.bseq[,"barcode"], all.tf.barcode)), ];
      P12.protein.num[j] <- nrow(P12.barcode.bseq)+P12.protein.num[j];
      cat(alignment.filename, "\tDone!", "\n");
    }
  }
  cat(total.ALL)
  total.ALL <- c(total.ALL,total.ALL,total.ALL,total.ALL)
  all.num <- cbind(total.ALL ,total.num, left.num, middle.num, right.num, bseq20.num, A.num, protein.num, P12.protein.num);
  write.csv(all.num, paste0("Raw_data/matching.sta/",index,"all.num.N16.csv"), row.names=FALSE);
}
# N16.protein.barcode.and.bseq.sequence()
# N16.protein.barcode.and.bseq.sequence("E_H",annotation.file="annotation/annotation_ETS.csv",Split_Protein_dir="Raw_Data/Split_Protein",control_lib="N16",plateID=c("twoCyc_rep1","twoCyc_rep2"),
#                                       plateID_indexL=c("","AA"),plateID_indexR=c("CT","CT"))

##########################################################
#
# individual protein barcode bseq relationship
# 
##########################################################
annotation.file <- "annotation/annotation_ETS.csv"
N16.individual.protein.barcode.bseq <- function(index="OneCyc_rep1")
{
  system(paste0("mkdir Raw_data/protein/",index,"/"))
  TF.annotation <- read.csv(annotation.file,stringsAsFactors = FALSE)
  TF.annotation.Barcode <- gsub(" ", "",toupper(as.character(TF.annotation[,"Barcode"])))
  TF.file <- gsub(" ", "",toupper(as.character(TF.annotation[,4])))
  
  all.tf.barcode <- as.character(TF.annotation[,"Barcode"])
  RevCom <- c()
  for(i in 1:length(all.tf.barcode))
  {
    temp <- DNAString(x = as.character(all.tf.barcode[i]))
    temp <- reverseComplement(temp) 
    RevCom <- c(RevCom,as.character(temp))
  }
  all.tf.barcode <-RevCom
  
  #barcode.bseq.filename <- paste0("Raw_data/Barcode.bseq/", fastq.file[i], ".N16.barcode.bseq");
  barcode.bseq.file <- list.files(paste0("Raw_data/Barcode.bseq/"));
  barcode.bseq.file <- barcode.bseq.file[grep(index,barcode.bseq.file)]
  ##################################
  for(i in 1:length(barcode.bseq.file))
  {
    barcode.bseq.filename <- paste0("Raw_data/Barcode.bseq/", barcode.bseq.file[i]);
    barcode.bseq <- read.csv(barcode.bseq.filename);
    for(j in 1:length(all.tf.barcode))
    {
      tf.barcode.bseq <- barcode.bseq[barcode.bseq[,"barcode"]==all.tf.barcode[j], c(1,2,4)];
      tf.barcode.bseq.filename <- paste0("Raw_data/protein/",index,"/",as.character(TF.file[j]),".txt");
      if(nrow(tf.barcode.bseq)>0)
      {
        print(tf.barcode.bseq.filename)
        write.table(tf.barcode.bseq, file=tf.barcode.bseq.filename, sep="\t", append=TRUE, row.names=FALSE, col.names=FALSE);
      }   
    }
    cat(barcode.bseq.filename, "\tDone!", "\n");
  }
}
# barcode.bseq.file <- list.files(paste0("Raw_data/Barcode.bseq/"));
# rep1.file <- barcode.bseq.file[grep("twoCyc_rep1",barcode.bseq.file)]
# rep2.file <- barcode.bseq.file[grep("twoCyc_rep2",barcode.bseq.file)]
# for(i in 1:length(rep1.file))
# {
#   rep1 <- read.csv(paste0("Raw_data/Barcode.bseq/",rep1.file[i]),header =TRUE  );
#   rep2 <- read.csv(paste0("Raw_data/Barcode.bseq/",rep2.file[i]),header =TRUE  );
#   barcode.bseq <- rbind(rep1,rep2)
#   barcode.bseq.filename <- paste0("Raw_data/Barcode.bseq/Merge",rep1.file[i])
#   write.csv(barcode.bseq, file=barcode.bseq.filename, row.names=FALSE);
# }
# barcode.bseq.file <- list.files(paste0("Raw_data/Barcode.bseq/"));

# N16.individual.protein.barcode.bseq("OneCyc_rep1")
# N16.individual.protein.barcode.bseq("OneCyc_rep2")
# N16.individual.protein.barcode.bseq("OneCyc_rep3")
# N16.individual.protein.barcode.bseq("MergetwoCyc_rep1")

N16.individual.protein.barcode.bseq.frequency <- function(index="OneCyc_rep1")
{ 
  
  library("Biostrings");
  all.barcode.bseq.file <- list.files(paste0("Raw_data/protein/",index));
  
  TF.annotation <- read.csv(annotation.file,stringsAsFactors = FALSE)
  all.tf <- read.csv(annotation.file,stringsAsFactors = FALSE)
  TF.annotation.Barcode <- gsub(" ", "",toupper(as.character(TF.annotation[,"Barcode"])))
  TF.file <- gsub(" ", "",toupper(as.character(TF.annotation[,4])))
  
  all.tf.barcode <- as.character(TF.annotation[,"Barcode"])
  RevCom <- c()
  for(i in 1:length(all.tf.barcode))
  {
    temp <- DNAString(x = as.character(all.tf.barcode[i]))
    temp <- reverseComplement(temp) 
    RevCom <- c(RevCom,as.character(temp))
  }
  all.tf.barcode <-RevCom
  all.tf.barcode <- substr(all.tf.barcode, 2, nchar(all.tf.barcode))
  
  system(paste0("mkdir Raw_data/protein/",index,"/protein.fqc/"))
  for(i in 1:length(all.tf.barcode))
  { 
    tf.barcode.bseq.filename <- paste0("Raw_data/protein/",index,"/",as.character(TF.file[i]),".txt");
    if(!is.na(match(tf.barcode.bseq.filename, paste0("Raw_data/protein/",index,"/", all.barcode.bseq.file))))
    {
      tf.barcode.bseq <- read.table(tf.barcode.bseq.filename, sep="\t");
      # bseq and UMI frequency 
      paste.bseq.UMI <- paste0(tf.barcode.bseq[,1], "-",tf.barcode.bseq[,2]);
      bseq.UMI.fqc <- as.data.frame(table(paste.bseq.UMI));
      split.bseq.UMI <- data.frame(do.call('rbind', strsplit(as.character(bseq.UMI.fqc[,1]),'-',fixed=TRUE)));
      bseq.UMI.fqc <- cbind(bseq.UMI.fqc, split.bseq.UMI);
      bseq.UMI.fqc <- bseq.UMI.fqc[,c(3,4,2)];
      colnames(bseq.UMI.fqc) <- c("bseq", "UMI", "frequency");
      bseq.UMI.fqc <- bseq.UMI.fqc[order(bseq.UMI.fqc[,1], decreasing=TRUE), ];
      bseq.UMI.fqc <- bseq.UMI.fqc[order(bseq.UMI.fqc[,3], decreasing=TRUE), ];
      
      # bseq frequency
      bseq.fqc <- as.data.frame(table(tf.barcode.bseq[,1]));
      colnames(bseq.fqc) <- c("bseq", "frequency");
      bseq.fqc <- bseq.fqc[order(bseq.fqc[,2], decreasing=TRUE), ];
      bseq.match <- match(bseq.UMI.fqc[,"bseq"], bseq.fqc[,"bseq"]);
      bseq.match.fqc <- as.data.frame(table(bseq.match));
      bseq.fqc <- cbind(bseq.fqc, bseq.match.fqc[,2]);
      colnames(bseq.fqc) <- c("bseq", "frequency", "UMI.num");
      tf.bseq.fqc.filename <- paste0("Raw_data/protein/",index,"/protein.fqc/",as.character(TF.file[i]),".csv");
      write.csv(bseq.fqc, file=tf.bseq.fqc.filename, row.names=FALSE);
      
      cat(i, "\t", as.character(TF.file[i]), "\tDone!", "\n");
    }
  }
}
# N16.individual.protein.barcode.bseq.frequency("OneCyc_rep1")
# N16.individual.protein.barcode.bseq.frequency("OneCyc_rep2")
# N16.individual.protein.barcode.bseq.frequency("OneCyc_rep3")
# N16.individual.protein.barcode.bseq.frequency("MergetwoCyc_rep1")


N16.individual.protein.6mer.frequency.4096 <- function(index="OneCyc_rep1")
{
  
  library("Biostrings");
  # all 6mer
  voc <- c("A","C","G","T");
  x.mer <- voc; 
  for(i in 2:6)
  {
    cur.rep.times <- rep(length(x.mer), length(voc));
    x.mer <- rep(x.mer, length(voc));
    add.mer <- rep(voc, cur.rep.times);
    x.mer <- paste0(add.mer, x.mer);
  }
  
  # calculate 6mer frequency---DNA20
  system(paste0("mkdir Raw_data/protein/",index,"/protein.fqc.6mer.4096/"))
  all.tf.bseq.fqc.DNA20.filename <- list.files(paste0("Raw_data/protein/",index,"/protein.fqc/"));
  for(i in 1:length(all.tf.bseq.fqc.DNA20.filename))
  {
    tf.bseq.fqc.DNA20.filename <- paste0("Raw_data/protein/",index,"/protein.fqc/", all.tf.bseq.fqc.DNA20.filename[i]);
    bseq.DNA20 <- read.csv(tf.bseq.fqc.DNA20.filename);
    bseq <- as.character(bseq.DNA20[,"bseq"]);
    fqc.6mer <- array(0, dim=length(x.mer));
    for(j in 1:length(x.mer))
    {
      bseq.6mer <- grep(x.mer[j], bseq);
      idv.6mer.fqc <- array(0, dim=length(bseq));
      idv.6mer.fqc[bseq.6mer] <- 1;
      fqc.6mer[j] <- sum(idv.6mer.fqc*bseq.DNA20[,"UMI.num"]);
      if(floor(j/100)==j/100)
      {
        cat(j, "\n");
      }
    }
    tf.6mer.fqc <- data.frame(x.mer, fqc.6mer);
    colnames(tf.6mer.fqc) <- c("bseq", "fqc");
    tf.6mer.fqc <- tf.6mer.fqc[order(tf.6mer.fqc[,2], decreasing=TRUE), ];
    tf.6mer.fqc.DNA20.filename <- paste0("Raw_data/protein/",index,"/protein.fqc.6mer.4096/", all.tf.bseq.fqc.DNA20.filename[i]);
    write.csv(tf.6mer.fqc, tf.6mer.fqc.DNA20.filename, row.names=FALSE);
    cat(i, "\t", as.character(all.tf.bseq.fqc.DNA20.filename[i]), "\tDone!", "\n");
  }
}
# N16.individual.protein.6mer.frequency.4096("OneCyc_rep1")
# N16.individual.protein.6mer.frequency.4096("OneCyc_rep2")
# N16.individual.protein.6mer.frequency.4096("OneCyc_rep3")
# N16.individual.protein.6mer.frequency.4096("MergetwoCyc_rep1")



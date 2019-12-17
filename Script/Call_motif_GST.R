##################################
#
# get the special sequences than call motifs with the GST background
#
##################################
protein.GST.sequences.two.steps <- function(protein.file,
                                            protein.dir,
                                            full_sequence.dir,
                                            GST.file,
                                            Scanplot_GST.dir,
                                            Seq_GST_fasta.dir,
                                            modif.lib,
                                            Solpe.guarantee=1,
                                            reads.cutoff)
{
  ##@protein.file,6mer file name for TF
  ##@protein.dir,6mer dir
  ##@full_sequence.dir,sequence dir
  ##@GST.file,6mer file name for GST
    
  library("seqinr")
  library(showtext)
  GST <- read.csv(GST.file) 
  print(protein.file)
  full_sequence <- read.csv(paste0(full_sequence.dir,"/",protein.file),stringsAsFactors = F)
  if(sum(as.numeric(full_sequence[,3])) < reads.cutoff)
  {
    print(paste0("!!!!!!!!sequence count is less then ",reads.cutoff,"!!!!!!!!"))
    return()
  }
  ###6mer
  modif.file <- read.csv(paste0(protein.dir,"/",protein.file),stringsAsFactors = F) 
  modif.GST <- cbind(GST,modif.file[match(GST[,1],modif.file[,1]),2])
  modif.GST[,3] <- as.numeric(modif.GST[,3])*sum(as.numeric(modif.GST[,2]),na.rm = T)/sum(as.numeric(modif.GST[,3]),na.rm = T)
  
  y.modif.lib = modif.lib
  ##############xlab fold for  picture#####################
  xy_max = max(c(500,as.numeric(modif.GST[,2]),as.numeric(modif.GST[,3])),na.rm = T)
  xy_max_temp <- xy_max
  while (xy_max_temp > 10)
  {
    xy_max_temp = xy_max_temp/10
  }
  xlim_fold <- log(xy_max/xy_max_temp,base = 10)
  xy_max <- ceiling(xy_max_temp) * xy_max/xy_max_temp
  
  pdf(file = paste0(Scanplot_GST.dir,"/",unlist(strsplit(protein,split = ".csv"))[1],".",Solpe.guarantee,".pdf"),height = 5,width = 5)
  par(mar = c(4,6,5,2))
  showtext_begin()
  plot(modif.GST[, c(2, 3)], xlab = eval(parse(text = paste0("expression('GST (X 10'^", xlim_fold, "*')')"))), ylab = eval(parse(text = paste0("expression('", toupper(unlist(strsplit(protein, "-"))[3]), " (X 10'^", xlim_fold, "*')')"))), xlim = c(0, xy_max), 
       ylim = c(0, xy_max),  cex = 0.5, cex.main = 2, font.main = 2, font.lab = 1.6, cex.lab = 1.6, cex.axis = 1.6, mgp = c(3, 0.1, 0), col = "black", xaxt = "n", yaxt = "n", family = "Arial")
  if (xy_max_temp < 2) {
    axis(side = 1, at = seq(0, ceiling(xy_max_temp), by = 0.5) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 0.5), 
         cex.axis = 1.6, family = "Arial")
    axis(side = 2, at = seq(0, ceiling(xy_max_temp), by = 0.5) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 0.5), 
         las = 1, cex.axis = 1.6, family = "Arial")
  } else {
    if (xy_max_temp <= 4) {
       axis(side = 1, at = seq(0, ceiling(xy_max_temp), by = 1) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 1), 
           cex.axis = 1.6, family = "Arial")
      axis(side = 2, at = seq(0, ceiling(xy_max_temp), by = 1) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 1), 
           las = 1, cex.axis = 1.6, family = "Arial")
    } else 
      {
       axis(side = 1, at = seq(0, ceiling(xy_max_temp), by = 2) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 2), 
           cex.axis = 1.6, family = "Arial")
      axis(side = 2, at = seq(0, ceiling(xy_max_temp), by = 2) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 2), 
           las = 1, cex.axis = 1.6, family = "Arial")
    }
  }
  sequencematrix <- modif.GST[order(modif.GST[, 3], decreasing = T), c(1, 2, 3)]
  sequencematrix <- data.frame(as.character(sequencematrix[, 1]), as.numeric(sequencematrix[, 2]), as.numeric(sequencematrix[,3]))
  sequencematrix <- sequencematrix[which(as.numeric(sequencematrix[, 3]) >= 1 &  as.numeric(sequencematrix[, 2]) >= 1),]
  sequencematrix <- data.frame(sequencematrix, as.numeric(sequencematrix[, 3])/as.numeric(sequencematrix[, 2]))
  colnames(sequencematrix) <- c("basq", "GST.count", modif.lib, paste0(modif.lib, "_GST.ratio"))
  Slope_cutoff <- max(summary(as.numeric(sequencematrix[, 4]))[5], Solpe.guarantee)
  addlineID1 <- which(as.numeric(sequencematrix[, 4]) > Slope_cutoff )
  addlineID2 <- which(as.numeric(sequencematrix[, 3]) > median(as.numeric(sequencematrix[, 3])))
  addlineID <- addlineID2[which(addlineID2 %in%  addlineID1 )]
  points(sequencematrix[addlineID, c(2, 3)], col = "red", cex = 0.5)  
  dev.off()
  
  ############
  if(length(addlineID) <= 3)
  {
    addlineID <- NULL
  }

  if(length(addlineID) > 0)
  {
    specific_mer <-  as.character(sequencematrix$basq[addlineID])
    has.specific_mer.sequence <- c()
    for(i in 1:length(specific_mer)){
      has.specific_mer.index <- grep(specific_mer[i],as.character(full_sequence[,1]))
      if(length(has.specific_mer.index) > 0){
        has.specific_mer.sequence <- c(has.specific_mer.sequence, 
                                       rep(as.character(full_sequence$bseq[has.specific_mer.index]),as.integer(full_sequence$UMI.num[has.specific_mer.index])))
        full_sequence <- full_sequence[-has.specific_mer.index,]
      }
    }
    write.fasta(sequences = as.list(has.specific_mer.sequence), names = as.character(seq(1:length(has.specific_mer.sequence))), paste0(Seq_GST_fasta.dir,"/",unlist(strsplit(protein,split = ".csv"))[1],".",Solpe.guarantee,".fasta"), nbchar =60, as.string = TRUE,open = "w")
  }else
  {
    pdf(file = paste0(Scanplot_GST.dir,"/",unlist(strsplit(protein,split = ".csv"))[1],".",Solpe.guarantee,".pdf"),height = 5,width = 5)
    par(mar = c(4, 6, 5, 2))
    showtext_begin()
    plot(modif.GST[, c(2, 3)], xlab = eval(parse(text = paste0("expression('GST (X 10'^", xlim_fold, "*')')"))), ylab = eval(parse(text = paste0("expression('", toupper(unlist(strsplit(protein, "-"))[3]), " (X 10'^", xlim_fold, "*')')"))), xlim = c(0, xy_max), 
         ylim = c(0, xy_max),  cex = 0.5, cex.main = 2, font.main = 2, font.lab = 1.6, cex.lab = 1.6, cex.axis = 1.6, mgp = c(3, 0.1, 0), col = "black", xaxt = "n", yaxt = "n", family = "Arial")
    if (xy_max_temp < 2) {
      axis(side = 1, at = seq(0, ceiling(xy_max_temp), by = 0.5) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 0.5), 
           cex.axis = 1.6, family = "Arial")
      axis(side = 2, at = seq(0, ceiling(xy_max_temp), by = 0.5) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), by = 0.5), 
           las = 1, cex.axis = 1.6, family = "Arial")
    } else {
      if (xy_max_temp <= 4) { 
        axis(side = 1, at = seq(0, ceiling(xy_max_temp), by = 1) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), 
                                                                                              by = 1), cex.axis = 1.6, family = "Arial")
        axis(side = 2, at = seq(0, ceiling(xy_max_temp), by = 1) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), 
                                                                                              by = 1), las = 1, cex.axis = 1.6, family = "Arial")
      } else
      {
        axis(side = 1, at = seq(0, ceiling(xy_max_temp), by = 2) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), 
                                                                                              by = 2), cex.axis = 1.6, family = "Arial")
        axis(side = 2, at = seq(0, ceiling(xy_max_temp), by = 2) * 10^xlim_fold, labels = seq(0, ceiling(xy_max_temp), 
                                                                                              by = 2), las = 1, cex.axis = 1.6, family = "Arial")
      }
    }
    dev.off()
  }
  return()
}
call_protein.GST.sequences.two.steps<- function(index= "DAPPL")
{
  protein.dir <- paste0("data/protein.fqc.6mer.4096/",index,"_protein.fqc.6mer.4096/")
  protein.file.list <- list.files(protein.dir)
  protein.file.list <- protein.file.list[-grep("p99A",protein.file.list)]#p99 is the GST index
  
  full_sequence.dir <- paste0("data/protein.fqc/",index,"_protein.fqc/")
  GST.file <- paste0("data/protein.fqc.6mer.4096/",index,"_protein.fqc.6mer.4096/p99A-H99-GST total-ATAGTCTCT.csv")
  
  Scanplot_GST.dir <- paste0("Result/",index,"_Scanplot")
  system(paste0("mkdir ",Scanplot_GST.dir))
  
  Seq_GST_fasta.dir <- paste0("Result/",index,"_logo.seq")
  system(paste0("mkdir ",Seq_GST_fasta.dir))
  
  modif.lib <- "N20"
  Solpe.guarantee <- 1
  y.modif.lib <- "N20"
  for(i in 1:length(protein.file.list))
  {
    protein.GST.sequences.two.steps(protein.file.list[i], protein.dir, full_sequence.dir, GST.file, Scanplot_GST.dir, Seq_GST_fasta.dir, modif.lib, Solpe.guarantee)
  }
}
#call_protein.GST.sequences.two.steps()

##################################
#
# mask the shell script and ues the homer to call motifs
#
##################################
make_sh <- function(index= "DAPPL",GST_file=P1112.N20.GST.fasta,output="P1112.N20.sh")
{
  
  protein.dir <- paste0("Result/",index,"_logo.seq/")
  Homer.dir <- paste0("Result/",index,"_Homer/")
  system(paste0("mkdir ",Homer.dir))
  protein.file.list <- list.files(protein.dir)
  system.line <- c()
  for(i in 1:length(protein.file.list))
  {
    system.line <- c(system.line,paste0("findMotifs.pl ",protein.dir,protein.file.list[i]," fasta ",Homer.dir,
                                    unlist(strsplit(protein.file.list[i],split = ".fasta"))[1],
                                    " -fasta data/protein.fqc/",GST_file," -noknown -S 2 -len 6,8,10"))
  }
  write.table(system.line,file = output,col.names = F,row.names = F,quote = F)
}
#make_sh()


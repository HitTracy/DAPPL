######compare the TF croess the modification###############
get.GST.correct<- function(scantplot_dir="/Volumes/Run_data/DAPPL_review/Result/Compare_modified/Scant_plot_fold/",mer_inputdir="/Volumes/Run_data/DAPPL_Before_review/DAPPL_Terminator_Modifed_compare4/DAPPL.2080/",mer_outputdir="/Volumes/Run_data/DAPPL_review/Result/Compare_modified/6mer_for_logo/",Ksd=6)
{ 
  library(xlsx)
  library(showtext) 
  for(symorhemi in c("sym","hemi"))
  {
    ###read the summary about the TF's motif, here just compare the TF's binding which has motif in the one libraray.
    Has.GST.log <- as.matrix(read_xlsx(path  = "annotation/GST.summary.xlsx",sheet = "Sym.Hemi"))
    Has.GST.log[which(is.na(Has.GST.log))] <- "-"
    lib1.name <- paste0(symorhemi,".CGT")
    lib2.name.ALL <- paste0(symorhemi,c(".CGT",".AGC",".CTA",".AAC",".ATA"))
    lib1.xlab <- paste0(symorhemi,".CG")
    lib2.ylab.ALL <- paste0(symorhemi,c(".CG",".mCG",".hmCG",".fCG",".caCG"))
    ### for the figure ylab
    if(symorhemi=="hemi")
    {
      lib2.ylab.ALL_alt <- c("C/C","mC/C","hmC/C","fC/C","caC/C")
    }else
    {
      lib2.ylab.ALL_alt <- c("C/C","mc/mc","hmC/hmC","fC/fC","caC/caC")
    }
    Has.GST.log <- Has.GST.log[order(Has.GST.log[,1],decreasing = F),]
    Has.GST.log <- Has.GST.log[which(Has.GST.log[,lib2.ylab.ALL[2]]=="+"|
                                       Has.GST.log[,lib2.ylab.ALL[3]]=="+"|
                                       Has.GST.log[,lib2.ylab.ALL[4]]=="+"|
                                       Has.GST.log[,lib2.ylab.ALL[5]]=="+"|
                                       Has.GST.log[,lib2.ylab.ALL[1]]=="+"),]
    Plate.All <- c("P12","P34","P56","P78","P910","P1112","P1314","P1517")
    Border_plate <- NULL
    Four_pattern_ALL <- NULL    
    for(plate in Plate.All)
    {
      protein.basename <-  as.character(Has.GST.log[,"File_ID"])
      plateID <-unlist(lapply(protein.basename, function(x){return(min(8,ceiling(as.numeric(unlist(strsplit(unlist(strsplit( x ,split = "A-"))[1],split = "p"))[2])/2)))}))
      TF_all <- as.character(Has.GST.log[which(Plate.All[plateID]==plate),1])
      for(protein.basename_TF in TF_all)
      {
        j <- which(as.character(Has.GST.log[,1])==protein.basename_TF)
        print(protein.basename_TF)
        Four_pattern <- NULL
        mean_distance <- NULL
        sd_distance <- NULL
        for(libID in c(2:5))
        {
          print(lib2.name.ALL[libID])
          lib2.name <- lib2.name.ALL[libID]
          if(!file.exists(paste0(scantplot_dir,symorhemi,"/",lib2.ylab.ALL[libID])))
          {
            system(paste0("mkdir ",scantplot_dir,symorhemi,"/",lib2.ylab.ALL[libID]))
          }
          if(!file.exists(paste0("/Volumes/Run_data/DAPPL_review/Result/Compare_modified/Scant_plot_fold/",symorhemi,"/",lib2.ylab.ALL[libID])))
          {
            system(paste0("mkdir /Volumes/Run_data/DAPPL_review/Result/Compare_modified/Scant_plot_fold/",symorhemi,"/",lib2.ylab.ALL[libID]))
          }
          if(!file.exists(paste0(mer_outputdir,symorhemi,"/",lib2.name)))
          {
            system(paste0("mkdir ",mer_outputdir,symorhemi,"/",lib2.name))
          }
          
          modifID <- unlist(strsplit(lib2.name.ALL[libID],split = "\\."))[2]
          hmCG.TF <- read.csv(paste0(mer_inputdir,plate,".",symorhemi,".2080/",modifID,"/",protein.basename_TF,".csv"))
          CG.TF <- read.csv(paste0(mer_inputdir,plate,".",symorhemi,".2080/CGT/",protein.basename_TF,".csv"))
          CG_hmCG.TF <- cbind(CG.TF,hmCG.TF[match(as.character(CG.TF[,1]),as.character(hmCG.TF[,1])),3])
          modif.control_TF <- cbind(as.numeric(as.character(CG_hmCG.TF[,3])),as.numeric(as.character(CG_hmCG.TF[,4])))
           
          rownames(modif.control_TF) <- as.character(CG_hmCG.TF[,1])
          colnames(modif.control_TF) <- c(lib1.name,lib2.name)#CG, modified
          
          #################read GST ##########################
          hmCG.GST <- read.csv(paste0(mer_inputdir,plate,".",symorhemi,".2080/",modifID,"/p99A-H99-GST total-ATAGTCTCT.csv"))
          CG.GST <- read.csv(paste0(mer_inputdir,plate,".",symorhemi,".2080/CGT/p99A-H99-GST total-ATAGTCTCT.csv"))
          CG_hmCG.GST <- cbind(CG.GST,hmCG.GST[match(as.character(CG.GST[,1]),as.character(hmCG.GST[,1])),3])
          modif.control_GST <- cbind(as.numeric(as.character(CG_hmCG.GST[,3])),as.numeric(as.character(CG_hmCG.GST[,4])))
      
          row.names(modif.control_GST) <- as.character(CG_hmCG.GST[,1])
          colnames( modif.control_GST) <- c(lib1.name,lib2.name)
          
          ##########lowess GST############################
          x <- modif.control_GST[,1]
          y <- modif.control_GST[,2]
          index <- which(x >1 & y>1)
          x <- x[index]
          y <- y[index]
          modif.control_GST <- modif.control_GST[index,]

          fold <- max(modif.control_TF[,1])/max(x)
          x <- x*fold
          y <- y*fold
          M <- log2(x/y)
          A <- 1/2*log2(x*y)
          out_lowess <- lowess(A,M,f = 0.2,iter = 2)
          xy_max <- max(c(x,y))
 
          y_new <- rep(0,length(y))
          for(i in 1:length(y))
          {
            y_new[i] <-y[i]*2^out_lowess$y[out_lowess$x==A[i]]
          }
          
          modif.control_GST_adjust <- cbind(x,y_new)
          rownames(modif.control_GST_adjust) <- rownames(modif.control_GST)
          rm(modif.control_GST)
          distance_GST <- abs(modif.control_GST_adjust[,2]-modif.control_GST_adjust[,1])/sqrt(2) #distance to y=x

          modif.control_TF <- modif.control_TF[row.names(modif.control_GST_adjust),]
          
          
          #################lowess ############################
          x_TF <- modif.control_TF[,1]
          y_TF <- modif.control_TF[,2]
          M_TF <- log2(x_TF/y_TF)
          A_TF <- 1/2*log2(x_TF*y_TF)
          
          MA_GST_model <- smooth.spline(out_lowess$x,out_lowess$y)
          M_TF_cent <- predict(MA_GST_model,x = A_TF)$y
          y_TF_new <-y_TF*2^(M_TF_cent)
          
          modif.control_TF_adjust <- cbind(x_TF,y_TF_new)
          row.names(modif.control_TF_adjust)<- row.names(modif.control_TF)
          distance_TF <- abs(modif.control_TF_adjust[,2]-modif.control_TF_adjust[,1])/sqrt(2)# distance y=x!!!!!!!!!!!!!

          ##########compare##########
          Intercept_GST <- (mean(distance_GST)*fold+Ksd*sd(distance_GST)*fold) 

          index_red <- which(distance_TF > Intercept_GST & (modif.control_TF_adjust[,2]+1)/(modif.control_TF_adjust[,1]+1) >1)
          index_blue <- which(distance_TF > Intercept_GST & (modif.control_TF_adjust[,2]+1)/(modif.control_TF_adjust[,1]+1) < 1)
          
          pdf(file  = paste0(scantplot_dir,symorhemi,"/",lib2.ylab.ALL[libID],"/",plate,"_",protein.basename_TF,"_Corrected_distance.pdf"),height=5, width=5)
          par(mar = c(4, 6, 5, 2))
          showtext_begin()
          modif.control_GST_adjust_rev <- modif.control_GST_adjust
          xy_max <- max(na.rm = T,c(500,modif.control_GST_adjust_rev[,1],modif.control_GST_adjust_rev[,2],modif.control_TF_adjust[,2],modif.control_TF_adjust[,1]))
          xy_max_temp <- xy_max
          while (xy_max_temp > 10) {xy_max_temp=(xy_max_temp/10)}
          xlim_fold <- log(xy_max/xy_max_temp,base = 10)
          xy_max <- ceiling(xy_max_temp)*xy_max/xy_max_temp
          plot(modif.control_GST_adjust_rev,xlab=eval(parse(text=paste0("expression('C/C (X 10'^",xlim_fold,"*')')"))),ylab=eval(parse(text=paste0("expression('",lib2.ylab.ALL_alt[libID]," (X 10'^",xlim_fold,"*')')"))),xlim=c(0,xy_max),ylim = c(0,xy_max),main=paste0("Corrected ",toupper(unlist(strsplit(protein.basename_TF,"-"))[3])),
                 cex=0.5,cex.main=2,font.main=2,font.lab=1.6,cex.lab=1.6,cex.axis=1.6,mgp=c(3,0.1,0),col="grey",xaxt="n",yaxt="n",family = 'Arial')
          if(xy_max_temp<=2)
           {
            axis(side=1,at=seq(0,ceiling(xy_max_temp),by=0.5)*10^xlim_fold,labels =seq(0,ceiling(xy_max_temp),by=0.5),cex.axis=1.6,family = 'Arial')
            axis(side=2,at=seq(0,ceiling(xy_max_temp),by=0.5)*10^xlim_fold,labels =seq(0,ceiling(xy_max_temp),by=0.5),las = 1,cex.axis=1.6,family = 'Arial')
           }else
          {
            if(xy_max_temp<=4)
            {
              axis(side=1,at=seq(0,ceiling(xy_max_temp),by=1)*10^xlim_fold,labels =seq(0,ceiling(xy_max_temp),by=1),cex.axis=1.6,family = 'Arial')
              axis(side=2,at=seq(0,ceiling(xy_max_temp),by=1)*10^xlim_fold,labels =seq(0,ceiling(xy_max_temp),by=1),las = 1,cex.axis=1.6,family = 'Arial')

            }else
            {
              axis(side=1,at=seq(0,ceiling(xy_max_temp),by=2)*10^xlim_fold,labels =seq(0,ceiling(xy_max_temp),by=2),cex.axis=1.6,family = 'Arial')
              axis(side=2,at=seq(0,ceiling(xy_max_temp),by=2)*10^xlim_fold,labels =seq(0,ceiling(xy_max_temp),by=2),las = 1,cex.axis=1.6,family = 'Arial')
            }
          }
          points(modif.control_TF_adjust,col="black",cex=0.5)
          
          ### After select the point we have to check the GST motif
          if(Has.GST.log[j,lib2.ylab.ALL[libID]]!="+")
          {
            index_red <- NULL
          }
          if(Has.GST.log[j,lib1.xlab]!="+")
          {
            index_blue <- NULL
          }

          count_point_cutoff = 3
          if(length(index_blue) < count_point_cutoff & length(index_red) < count_point_cutoff)
          {
          }else
          {
            if(length(index_blue) >= count_point_cutoff & length(index_red) >= count_point_cutoff )#"Bilateral"
            {
              redlm <- lm(modif.control_TF_adjust[index_red,2]~modif.control_TF_adjust[index_red,1])
              redlm_2<- lm(modif.control_TF_adjust[index_red,1]~modif.control_TF_adjust[index_red,2])
              if(redlm$coefficients[2]> 1/redlm_2$coefficients[2])
              {
              }else
              {
                redlm$coefficients[2] <- 1/redlm_2$coefficients[2]
                redlm$coefficients[1] <- -redlm_2$coefficients[1]*(1/redlm_2$coefficients[2])
              }
              abline(redlm,col="red")
              points(modif.control_TF_adjust[index_red,],col="red",cex=0.5)
              bluelm <- lm(modif.control_TF_adjust[index_blue,2]~modif.control_TF_adjust[index_blue,1])
              abline(bluelm,col="blue")
              points(modif.control_TF_adjust[index_blue,],col="blue",cex=0.5)
            }else
            {
              if(length(index_blue) >= count_point_cutoff )#"Suppress"
              {
                bluelm <- lm(modif.control_TF_adjust[index_blue,2]~modif.control_TF_adjust[index_blue,1])
                abline(bluelm,col="blue")
                points(modif.control_TF_adjust[index_blue,],col="blue",cex=0.5)
              }else#"Enhance"
              {
                redlm <- lm(modif.control_TF_adjust[index_red,2]~modif.control_TF_adjust[index_red,1])
                redlm_2<- lm(modif.control_TF_adjust[index_red,1]~modif.control_TF_adjust[index_red,2])
                if(redlm$coefficients[2]> 1/redlm_2$coefficients[2] & redlm$coefficients[2] > 0)
                {

                }else
                {
                  redlm$coefficients[2] <- 1/redlm_2$coefficients[2]
                  redlm$coefficients[1] <- -redlm_2$coefficients[1]*(1/redlm_2$coefficients[2])
                }
                abline(redlm,col="red")
                points(modif.control_TF_adjust[index_red,],col="red",cex=0.5)   
              }
            }
          }
          showtext_end()
          dev.off()

          if(length(index_blue) < count_point_cutoff & length(index_red) < count_point_cutoff)
          {
            index_blue <- NULL
            index_red <- NULL
            if(Has.GST.log[j,lib1.xlab]=="+" & Has.GST.log[j,lib2.ylab.ALL[libID]]=="+" )
            {
              Four_pattern <- c(Four_pattern,"No_preference")#"No_preference"
            }else
            {
              Four_pattern <- c(Four_pattern,"No_data")#"No_data"
            }

          }else
          {
            if(length(index_blue) >= count_point_cutoff & length(index_red) >= count_point_cutoff )#"Bilateral"
            {

              Four_pattern <- c(Four_pattern,"Bilateral")
              red <- modif.control_TF_adjust[index_red,]
              red <- red[order(red[,2],decreasing = T),]
              colnames(red ) <- c("CG",lib2.ylab.ALL[libID])
              write.csv(red,file = paste0(mer_outputdir,symorhemi,"/",lib2.name,"/",plate,"_",protein.basename_TF,".red.csv"))

              blue <- modif.control_TF_adjust[index_blue,]
              blue <- blue[order(blue[,1],decreasing = T),]
              colnames(blue )<- c("CG",lib2.ylab.ALL[libID])
              write.csv(blue,file = paste0(mer_outputdir,symorhemi,"/",lib2.name,"/",plate,"_",protein.basename_TF,".blue.csv"))
            }else
            {
              if(length(index_blue) >= count_point_cutoff )#"Suppress"
              {
                index_red <- NULL
                Four_pattern <- c(Four_pattern,"Suppress")
                blue <- modif.control_TF_adjust[index_blue,]
                blue <- blue[order(blue[,1],decreasing = T),]
                colnames(blue )<- c("CG",lib2.ylab.ALL[libID])
                write.csv(blue,file = paste0(mer_outputdir,symorhemi,"/",lib2.name,"/",plate,"_",protein.basename_TF,".blue.csv"))
              }else#"Enhance"
              {
                index_blue <- NULL
                Four_pattern <- c(Four_pattern,"Enhance")
                red <- modif.control_TF_adjust[index_red,]
                red <- red[order(red[,2],decreasing = T),]
                colnames(red ) <- c("CG",lib2.ylab.ALL[libID])
                write.csv(red,file = paste0(mer_outputdir,symorhemi,"/",lib2.name,"/",plate,"_",protein.basename_TF,".red.csv"))
              }
            }
          }
         }
        Four_pattern_ALL <- rbind(Four_pattern_ALL,c(protein.basename_TF,Four_pattern))
       }
     }#########这里还没整体！
    colnames(Four_pattern_ALL) <-c("ID", paste0(symorhemi,c(".mCG",".hmCG",".fCG",".caCG")))
    write.csv(Four_pattern_ALL,file = paste0(scantplot_dir,symorhemi,"_four_pattern_summary.csv"),row.names = F)
    Four_pattern_ALL_plot <- Four_pattern_ALL[-which(Four_pattern_ALL[,2]=="No_data" & Four_pattern_ALL[,3]=="No_data"&
                                                       Four_pattern_ALL[,4]=="No_data" & Four_pattern_ALL[,5]=="No_data"),]
    write.csv(Four_pattern_ALL_plot,file = paste0(scantplot_dir,symorhemi,"_four_pattern_summary_select.csv"),row.names = F)
  }
}
#########main##############
get.GST.correct()
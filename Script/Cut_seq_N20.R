##########################################################
#
#	DNA20 sequence alignment
#	
##########################################################
DNA20.alignment.matching <- function(fastq.filename=NULL, alignment.filename=NULL)
{
	library("Biostrings");
	read.seq <- readLines(fastq.filename);
	read.seq <- read.seq[seq(2,length(read.seq),4)]
	print(length(read.seq))
	### matching matrix
	mat <- matrix(-1,5,5);
	mat[1,1] <- 1;
	mat[2,2] <- 1;
	mat[3,3] <- 1;
	mat[4,4] <- 1;
	mat[5,5] <- 1;
	colnames(mat) <- c("A", "C", "G", "T", "N");
	rownames(mat) <- c("A", "C", "G", "T", "N");

	###left sequence "GGGAGAAGGTCATCAAGAGG" local alignment
	left.seq <- "GGGAGAAGGTCATCAAGAGG";
	print("strat left.align")
	left.align <- pairwiseAlignment(read.seq, left.seq, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -1);
	left.align.start <- start(pattern(left.align));
	left.align.end <- end(pattern(left.align));
	left.align.nchar <- nchar(left.align);
	left.align.score <- score(left.align);
	left.align.last4letter <- substr(as.character(pattern(left.align)), left.align.nchar-3, left.align.nchar);
	cat("left alignment Done!\n");
	### right sequence "CCTACGT" local alignment
	right.seq <- "CCTACGT";
	sub.read.seq <- substr(read.seq, 81, nchar(read.seq));
	right.align <- pairwiseAlignment(sub.read.seq, right.seq, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -1);
	right.align.start <- 80+start(pattern(right.align));
	right.align.end <- 80+end(pattern(right.align));
	right.align.nchar <- nchar(right.align);
	right.align.score <- score(right.align);
	right.align.first4letter <- substr(as.character(pattern(right.align)), 1, 4);
	cat("right alignment Done!\n");
	### middle sequence "GGCTAGCAGCCACTATAAGCTTCGAAGACTGG" local alignment
	middle.seq <- "GGCTAGCAGCCACTATAAGCTTCGAAGACTGG";
	middle.align <- pairwiseAlignment(read.seq, middle.seq, type = "local", substitutionMatrix = mat, gapOpening = -1, gapExtension = -0.2);
	middle.align.start <- start(pattern(middle.align));
	middle.align.end <- end(pattern(middle.align));
	middle.align.nchar <- nchar(middle.align);
	middle.align.score <- score(middle.align);
	middle.align.first5letter <- substr(as.character(pattern(middle.align)), 1, 5);
	middle.align.last4letter <- substr(as.character(pattern(middle.align)), middle.align.nchar-3, middle.align.nchar);
	cat("middle alignment Done!\n");

	DNA.align <- data.frame(left.align.score, left.align.nchar, left.align.start, left.align.end, left.align.last4letter, right.align.score, right.align.nchar, right.align.start, right.align.end, right.align.first4letter, middle.align.score, middle.align.nchar, middle.align.start, middle.align.end, middle.align.first5letter, middle.align.last4letter);
	write.csv(DNA.align, file=alignment.filename, row.names=FALSE);
}
call.DNA20.alignment.matching <- function()
{
  fastq.file <- list.files(paste0("data/Raw_data/",index,"_split/"));
  for(i in 1:length(fastq.file))
  {
	  fastq.filename <- paste0("data/Raw_data/",index,"_split/", fastq.file[i]);
	  alignment.filename <- paste0("data/alignment/", fastq.file[i], ".DNA20.alignment");
	  DNA20.alignment.matching(fastq.filename=fastq.filename, alignment.filename=alignment.filename);
	
	  alignment.filename <- paste0("data/alignment/", fastq.file[i], ".DNA20m.alignment");
	  DNA20m.alignment.matching(fastq.filename=fastq.filename, alignment.filename=alignment.filename);
	
	  cat(paste(fastq.filename, "\tDone!", "\n"));
  }
}

##########################################################
#
#	sequence alignment statistics
#	
##########################################################
DNA20.matching.score.statistics <- function()
{
	fastq.file <- list.files(paste0("data/Raw_data/",index,"_split/"));
	# mapping score statistics
	matching.sta <- NULL;
	for(i in 1:length(fastq.file))
	{
		alignment.filename <- paste0("data/alignment/", fastq.file[i], ".DNA20.alignment");
		alignment <- read.csv(alignment.filename);
		matching.sta <- rbind(matching.sta, rbind(alignment[,c(1,6,11)]));
		cat(alignment.filename, "\tDone!", "\n");
	}
	# left matching score histogram
	hist.left <- hist(matching.sta[,1], plot=FALSE, breaks=20);
	png(paste0("data/matching.sta/",index,"DNA20.left.align.score.png"))
	plot(hist.left, main="left sequence alignment", ylab="read number", xlab="matching score", col="lightblue", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5);	
	dev.off();
	# right matching score histogram
	hist.right <- hist(matching.sta[,2], plot=FALSE, breaks=20);
	png(paste0("data/matching.sta/",index,"DNA20.right.align.score.png"))
	plot(hist.right, main="right sequence alignment", ylab="read number", xlab="matching score", col="lightblue", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5);	
	dev.off();
	# middle matching score histogram
	hist.middle <- hist(matching.sta[,3], plot=FALSE, breaks=20);
	png(paste0("data/matching.sta/",index,"DNA20.middle.align.score.png"))
	plot(hist.middle, main="middle sequence alignment", ylab="read number", xlab="matching score", col="lightgreen", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5);	
	dev.off();
}

##########################################################
#
#	retrieve protein barcode and bseq sequence
#	
##########################################################
DNA20.protein.barcode.and.bseq.sequence <- function()
{
  fastq.file <- list.files(paste0("data/Raw_data/",index,"_split/"));
	all.tf <- read.csv(Barcode_annotation_files);
  all.tf.barcode <- as.character(all.tf[,"Barcode"])
	RevCom <- c()
	for(i in 1:length(all.tf.barcode))
	{
	  temp <- DNAString(x = as.character(all.tf.barcode[i]))
	  temp <- reverseComplement(temp) 
	  RevCom <- c(RevCom,as.character(temp))
	}
	all.tf.barcode <-RevCom
	
	total.num <- 0;
	left.num <- 0;
	middle.num <- 0;
	right.num <- 0;
	bseq20.num <- 0;
	A.num <- 0;
	protein.num <- 0;
	P12.protein.num <- 0;
	for(i in 1:length(fastq.file))
	{
		fastq.filename <- paste0("data/Raw_data/",index,"_split/", fastq.file[i]);
		alignment.filename <- paste0("data/alignment/", fastq.file[i], ".DNA20.alignment");
		barcode.bseq.filename <- paste0("data/barcode.bseq/", fastq.file[i], ".DNA20.barcode.bseq");
		read.seq <- readLines(fastq.filename);
		read.seq <- read.seq[seq(2,length(read.seq),4)]
		alignment <- read.csv(alignment.filename);
		total.num <- total.num + nrow(alignment);
		read.left <- read.seq[alignment[,"left.align.last4letter"]=="GAGG" & alignment[,"left.align.score"]>=16];
		align.left <- alignment[alignment[,"left.align.last4letter"]=="GAGG" & alignment[,"left.align.score"]>=16, ];
		left.num <- left.num + nrow(align.left);
		read.left.middle <- read.left[align.left[,"middle.align.first5letter"]=="GGCTA" & align.left[,"middle.align.last4letter"]=="CTGG" & align.left[,"middle.align.score"]>=20];
		align.left.middle <- align.left[align.left[,"middle.align.first5letter"]=="GGCTA" & align.left[,"middle.align.last4letter"]=="CTGG" & align.left[,"middle.align.score"]>=20, ];
		middle.num <- middle.num + nrow(align.left.middle);
		read.left.middle.right <- read.left.middle[align.left.middle[,"right.align.first4letter"]=="CCTA"];
		align.left.middle.right <- align.left.middle[align.left.middle[,"right.align.first4letter"]=="CCTA", ];
		right.num <- right.num+nrow(align.left.middle.right);
		read.left.middle.right.bseq <- read.left.middle.right[(align.left.middle.right[,"middle.align.start"]-align.left.middle.right[,"left.align.end"]-1)==20];
		align.left.middle.right.bseq <- align.left.middle.right[(align.left.middle.right[,"middle.align.start"]-align.left.middle.right[,"left.align.end"]-1)==20, ];
		bseq20.num <- bseq20.num+nrow(align.left.middle.right.bseq);
		barcode <- substr(read.left.middle.right.bseq, as.numeric(align.left.middle.right.bseq[,"middle.align.end"])+1, as.numeric(align.left.middle.right.bseq[,"right.align.start"])-1);
		bseq <- substr(read.left.middle.right.bseq, as.numeric(align.left.middle.right.bseq[,"left.align.end"])+1, as.numeric(align.left.middle.right.bseq[,"middle.align.start"])-1);
		UMI <- substr(barcode, 1, 8);
		A <- substr(barcode, 9, 9);
		barcode <- substr(barcode, 9, nchar(barcode));
		barcode.bseq <- data.frame(bseq, UMI, A, barcode);
		barcode.bseq <- barcode.bseq[barcode.bseq[,"A"]=="A" & nchar(as.character(barcode.bseq[,"UMI"]))==8, ];
		A.num <- nrow(barcode.bseq)+A.num;
	  barcode.bseq <- barcode.bseq[!is.na(match(barcode.bseq[,"barcode"], all.tf.barcode)), ];
	  protein.num <- nrow(barcode.bseq)+protein.num;
		write.csv(barcode.bseq, file=barcode.bseq.filename, row.names=FALSE);
		P12.barcode.bseq <- barcode.bseq[!is.na(match(barcode.bseq[,"barcode"], all.tf.barcode[1:204])), ];
		P12.protein.num <- nrow(P12.barcode.bseq)+P12.protein.num;
		cat(alignment.filename, "\tDone!", "\n");
	}
	all.num <- data.frame(total.num, left.num, middle.num, right.num, bseq20.num, A.num, protein.num, P12.protein.num);
	write.csv(all.num,paste0("data/matching.sta/",index,"all.num.DNA20.csv"), row.names=FALSE);
}

##########################################################
#
#	individual protein barcode bseq relationship
#	
##########################################################
DNA20.individual.protein.barcode.bseq <- function()
{
  system(paste0("mkdir data/protein/",index,"_protein/"))
  
  fastq.file <- list.files(paste0("data/Raw_data/",index,"_split/"));
  
	all.tf <- read.csv(Barcode_annotation_files);
  all.tf.barcode <- as.character(all.tf[,"Barcode"])
	RevCom <- c()
	for(i in 1:length(all.tf.barcode))
	{
	  temp <- DNAString(x = as.character(all.tf.barcode[i]))
	  temp <- reverseComplement(temp) 
	  RevCom <- c(RevCom,as.character(temp))
	}
	all.tf.barcode <-RevCom
	##################################
	for(i in 1:length(fastq.file))
	{
		barcode.bseq.filename <- paste0("data/barcode.bseq/", fastq.file[i], ".DNA20.barcode.bseq");
		barcode.bseq <- read.csv(barcode.bseq.filename);
		for(j in 1:length(all.tf.barcode))
		{
			tf.barcode.bseq <- barcode.bseq[barcode.bseq[,"barcode"]==all.tf.barcode[j], c(1,2,4)];
			tf.barcode.bseq.filename <- paste0("data/protein/",index,"_protein/",as.character(all.tf[j,"proteinID"]),".txt");
			if(nrow(tf.barcode.bseq)>0)
			{
			  print(tf.barcode.bseq.filename)
				write.table(tf.barcode.bseq, file=tf.barcode.bseq.filename, sep="\t", append=TRUE, row.names=FALSE, col.names=FALSE);
			}		
		}
		cat(barcode.bseq.filename, "\tDone!", "\n");
	}
}
 
##########################################################
#
#	calculate barcode and bseq pair frequency in matrix
#	
##########################################################
DNA20.individual.protein.barcode.bseq.frequency <- function()
{ 
  system(paste0("mkdir data/protein.fqc/",index,"_protein.fqc/"))
	library("Biostrings");
	all.barcode.bseq.file <- list.files(paste0("data/protein/",index,"_protein"));
	all.tf <- read.csv(Barcode_annotation_files);
	#all.tf.barcode <- substr(as.character(all.tf[,"TFbarcodes"]), 2, nchar(as.character(all.tf[,"TFbarcodes"])));
	all.tf.barcode <- as.character(all.tf[,"Barcode"])
	RevCom <- c()
	for(i in 1:length(all.tf.barcode))
	{
	  temp <- DNAString(x = as.character(all.tf.barcode[i]))
	  temp <- reverseComplement(temp) 
	  RevCom <- c(RevCom,as.character(temp))
	}
	all.tf.barcode <-RevCom
	all.tf.barcode <- substr(all.tf.barcode, 2, nchar(all.tf.barcode))
	
	for(i in 1:length(all.tf.barcode))
	{	
	  tf.barcode.bseq.filename <- paste0("data/protein/",index,"_protein/",as.character(all.tf[i,"proteinID"]),".txt");
		if(!is.na(match(tf.barcode.bseq.filename, paste0("data/protein/",index,"_protein/", all.barcode.bseq.file))))
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
			tf.bseq.fqc.filename <- paste0("data/protein.fqc/",index,"_protein.fqc/",as.character(all.tf[i,"proteinID"]),".csv");
			write.csv(bseq.fqc, file=tf.bseq.fqc.filename, row.names=FALSE);

			cat(i, "\t", as.character(all.tf[i,"proteinID"]), "\tDone!", "\n");
		}
	}
}

##########################################################
#	
#	DNA20 6mer frequency 4096
#	
##########################################################
DNA20.individual.protein.6mer.frequency.4096 <- function()
{
  system("mkdir data/protein.fqc.6mer.4096/")
  system(paste0("mkdir data/protein.fqc.6mer.4096/",index,"_protein.fqc.6mer.4096/"))
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
	all.tf.bseq.fqc.DNA20.filename <- list.files(paste0("data/protein.fqc/",index,"_protein.fqc/"));
	for(i in 1:length(all.tf.bseq.fqc.DNA20.filename))
	{
		tf.bseq.fqc.DNA20.filename <- paste0("data/protein.fqc/",index,"_protein.fqc/", all.tf.bseq.fqc.DNA20.filename[i]);
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
		tf.6mer.fqc.DNA20.filename <- paste0("data/protein.fqc.6mer.4096/",index,"_protein.fqc.6mer.4096/", all.tf.bseq.fqc.DNA20.filename[i]);
		write.csv(tf.6mer.fqc, tf.6mer.fqc.DNA20.filename, row.names=FALSE);
		cat(i, "\t", as.character(all.tf.bseq.fqc.DNA20.filename[i]), "\tDone!", "\n");
	}
}

##########################################################
#	
#	DNA20 6mer frequency 2080
#	
##########################################################
DNA20.individual.protein.6mer.frequency.2080 <- function()
{
  system("mkdir protein.fqc.6mer.2080/")
  system(paste0("mkdir protein.fqc.6mer.2080/",index,"_protein.fqc.6mer.2080/"))
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
	rev.x.mer <- x.mer;
	for(i in 1:length(rev.x.mer))
	{
		rev.x.mer[i] <- as.character(reverseComplement(DNAString(x.mer[i])));
	}
	x.mer <- x.mer[x.mer<=rev.x.mer];
	rev.x.mer <- x.mer;
	for(i in 1:length(rev.x.mer))
	{
		rev.x.mer[i] <- as.character(reverseComplement(DNAString(x.mer[i])));
	}
	
	# calculate 6mer frequency---DNA20
	all.tf.bseq.fqc.DNA20.filename <- list.files(paste0("data/protein.fqc/",index,"_protein.fqc"));
	for(i in 1:length(all.tf.bseq.fqc.DNA20.filename))
	{
		tf.bseq.fqc.DNA20.filename <- paste0("data/protein.fqc/",index,"_protein.fqc/", all.tf.bseq.fqc.DNA20.filename[i]);
		bseq.DNA20 <- read.csv(tf.bseq.fqc.DNA20.filename);
		bseq <- as.character(bseq.DNA20[,"bseq"]);
		fqc.6mer <- array(0, dim=length(x.mer));
		for(j in 1:length(x.mer))
		{
			bseq.6mer <- grep(x.mer[j], bseq);
			rev.bseq.6mer <- grep(rev.x.mer[j], bseq);
			idv.6mer.fqc <- array(0, dim=length(bseq));
			idv.6mer.fqc[unique(c(bseq.6mer, rev.bseq.6mer))] <- 1;
			fqc.6mer[j] <- sum(idv.6mer.fqc*bseq.DNA20[,"UMI.num"]);
			if(floor(j/100)==j/100)
			{
				cat(j, "\n");
			}
		}
		tf.6mer.fqc <- data.frame(x.mer, rev.x.mer, fqc.6mer);
		colnames(tf.6mer.fqc) <- c("bseq", "reverse", "fqc");
		tf.6mer.fqc <- tf.6mer.fqc[order(tf.6mer.fqc[,3], decreasing=TRUE), ];
		tf.6mer.fqc.DNA20.filename <- paste0("data/protein.fqc.6mer.2080/",index,"_protein.fqc.6mer.2080/", all.tf.bseq.fqc.DNA20.filename[i]);
		write.csv(tf.6mer.fqc, tf.6mer.fqc.DNA20.filename, row.names=FALSE);
		cat(i, "\t", as.character(all.tf.bseq.fqc.DNA20.filename[i]), "\tDone!", "\n");
	}
}





##############main##################
# the index is used to mark the sequenceing.
# the sequecen library(.fq) should be split by row counts with prefix = "index_split/index"
index <- "DAPPL"
Barcode_annotation_files <- "annotation/annotation_sample.csv"
call.DNA20.alignment.matching()
DNA20.matching.score.statistics()
DNA20.protein.barcode.and.bseq.sequence()
DNA20.individual.protein.barcode.bseq()
DNA20.individual.protein.barcode.bseq.frequency()
DNA20.individual.protein.6mer.frequency.2080()
DNA20.individual.protein.6mer.frequency.4096()

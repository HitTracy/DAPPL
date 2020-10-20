matrix.motif <- function(motif = NULL) {
    matrix = matrix(0, nrow = 4, ncol = nchar(motif))
    if (gregexpr("[aA]", motif)[[1]][1] > 0) {
        matrix[1, gregexpr("[aA]", motif)[[1]]] = 1
    }
    if (gregexpr("[cC]", motif)[[1]][1] > 0) {
        matrix[2, gregexpr("[cC]", motif)[[1]]] = 1
    }
    if (gregexpr("[gG]", motif)[[1]][1] > 0) {
        matrix[3, gregexpr("[gG]", motif)[[1]]] = 1
    }
    if (gregexpr("[tT]", motif)[[1]][1] > 0) {
        matrix[4, gregexpr("[tT]", motif)[[1]]] = 1
    }
    matrix
}

N16_draw_logo <- function(motif.file = "Seq_GST_Homer/CTCF/homerResults/motif1.motif", output_raw_logo) {
    if (!("package:Biostrings" %in% search())) {
        library(Biostrings)
    }
    if (!("package:stringr" %in% search())) {
        library(stringr)
    }
    library(ggplot2)
    library(ggseqlogo)
    
    motif <- read.delim(motif.file)[, c(1:4)]  #ACGT
    motif <- t(motif)
    rownames(motif) <- c("A", "C", "G", "T")
    
    ################ write raw motif###########
    New.motif <- motif
    csl2 <- make_col_scheme(chars = c("A", "C", "G", "T"), group = c("1", "2", "3", "4"), cols = c("#08d61d", "#0d09ed", "#edb409", 
        "#ed0909"))
    rownames(New.motif) <- c("A", "C", "G", "T")
    gp <- ggplot() + geom_logo(New.motif, method = "prob", col_scheme = csl2, show_guide = FALSE) + labs(y = "") + theme(panel.background = element_rect(fill = "transparent", 
        color = NA)) + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank())
    ggsave(output_raw_logo, gp, height = 4, width = 20)
}

match.logo.seq.modified_logo <- function(motif.file = "Seq_GST_Homer/CTCF/homerResults/motif1.motif", fasta.sequence.file = "logo.seq/CTCF.fasta", 
    modification.index = 9, output_raw_logo, output_modifed_logo) {
    # modification.index is the modified CG site in the N-mer,here is 18-mer.
    if (!("package:Biostrings" %in% search())) {
        library(Biostrings)
    }
    if (!("package:stringr" %in% search())) {
        library(stringr)
    }
    library(ggplot2)
    library(ggseqlogo)
    
    fasta.sequence <- readDNAStringSet(fasta.sequence.file)
    if (length(fasta.sequence) > 30000) {
        asta.sequence <- fasta.sequence[floor(runif(30000, 1, length(fasta.sequence)))]
    }
    
    motif <- read.delim(motif.file)[, c(1:4)]  #ACGT
    motif <- t(motif)
    rownames(motif) <- c("A", "C", "G", "T")
    
    ################ write raw motif###########
    New.motif <- motif
    csl2 <- make_col_scheme(chars = c("A", "C", "G", "T"), group = c("1", "2", "3", "4"), cols = c("#08d61d", "#0d09ed", "#edb409", 
        "#ed0909"))
    rownames(New.motif) <- c("A", "C", "G", "T")
    gp <- ggplot() + geom_logo(New.motif, method = "prob", col_scheme = csl2, show_guide = FALSE) + labs(y = "") + theme(panel.background = element_rect(fill = "transparent", 
        color = NA)) + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank())
    ggsave(output_raw_logo, gp, height = 4, width = 20)
    
    ################ match modified motif ###########
    SQ <- as.matrix(motif)
    NTperc = matrix(c(0.295, 0.205, 0.205, 0.295))  # percentage of nucleotide in A, C, G, and T
    pseudo = NTperc %*% sqrt(colSums(SQ))
    PWM = pseudo + SQ
    PSSM = log2(PWM/matrix(rep(colSums(PWM), 4), 4, byrow = TRUE)/matrix(rep(NTperc, dim(PWM)[2]), 4))
    
    sequence <- as.character(fasta.sequence)
    first.score <- array(-100, dim = length(sequence))
    first.position <- array(-100, dim = length(sequence))
    first.sequence <- array("", dim = length(sequence))
    first.strand <- array("", dim = length(sequence))
    len.pssm <- ncol(PSSM)
    
    seq.match.ALL <- NULL
    for (i in 1:length(sequence)) {
        seq.match <- NULL
        seq <- as.character(sequence[i])
        reverse.seq <- as.character(reverseComplement(DNAString(sequence[i])))
        m.length <- nchar(seq) - len.pssm + 1
        seq.length <- nchar(seq)
        matrix.score <- matrix(-100, nrow = 2, ncol = m.length)
        
        if (m.length > 0) {
            seq.matrix <- matrix.motif(seq)
            reverse.seq.matrix <- matrix.motif(as.character(reverse.seq))
            for (m in 1:m.length) {
                motif.matrix <- seq.matrix[, m:(m + len.pssm - 1)]
                reverse.motif.matrix <- reverse.seq.matrix[, (seq.length - len.pssm - m + 2):(seq.length - m + 1)]
                matrix.score[1, m] <- round(sum(motif.matrix * PSSM), 2)
                matrix.score[2, m] <- round(sum(reverse.motif.matrix * PSSM), 2)
            }
            
            position <- order(matrix.score, decreasing = TRUE)[1]
            
            if (position%%2 != 0) {
                first.score[i] <- matrix.score[position]
                first.position[i] <- floor(position/2) + 1
                first.sequence[i] <- substr(seq, first.position[i], first.position[i] + len.pssm - 1)
                first.strand[i] <- "+"
                
                if (modification.index - first.position[i] >= 0 & modification.index - first.position[i] < len.pssm) {
                  if (modification.index - first.position[i] == 0) {
                    seq.match <- paste0("E", str_sub(first.sequence[i], 2))
                  } else {
                    if (modification.index - first.position[i] + 2 > nchar(first.sequence[i])) {
                      seq.match <- paste0(str_sub(first.sequence[i], 1, modification.index - first.position[i]), "E")
                    } else {
                      seq.match <- paste0(str_sub(first.sequence[i], 1, modification.index - first.position[i]), "E", str_sub(first.sequence[i], 
                        modification.index - first.position[i] + 2))
                    }
                  }
                }
            } else {
                first.score[i] <- matrix.score[position]
                first.position[i] <- floor(position/2)
                first.sequence[i] <- substr(reverse.seq, seq.length - len.pssm - first.position[i] + 2, seq.length - first.position[i] + 
                  1)
                first.strand[i] <- "-"
                
                rev.first.position <- as.integer(nchar(sequence[i])) - first.position[i] - len.pssm + 2
                rev.modification.index <- as.integer(nchar(sequence[i])) - modification.index
                if (rev.modification.index - rev.first.position >= 0 & rev.modification.index - rev.first.position < len.pssm) {
                  if (rev.modification.index - rev.first.position == 0) {
                    seq.match <- paste0("E", str_sub(first.sequence[i], 2))
                  } else {
                    if (rev.modification.index - rev.first.position + 2 > nchar(first.sequence[i])) {
                      seq.match <- paste0(str_sub(first.sequence[i], 1, rev.modification.index - rev.first.position), "E")
                    } else {
                      seq.match <- paste0(str_sub(first.sequence[i], 1, rev.modification.index - rev.first.position), "E", str_sub(first.sequence[i], 
                        rev.modification.index - rev.first.position + 2))
                    }
                  }
                }
            }
        }
        
        if (length(seq.match) != 0) {
            seq.match.ALL <- c(seq.match.ALL, seq.match)
        }
    }
    
    New.motif <- matrix(nrow = nchar(seq.match.ALL[1]), ncol = 5)  #ACGTE
    for (i in 1:nchar(seq.match.ALL[1])) {
        temp.char <- unlist(lapply(seq.match.ALL, function(x) {
            return(str_sub(x, i, i))
        }))
        temp.char <- summary(as.factor(temp.char))/length(temp.char)
        New.motif[i, ] <- as.array(temp.char[c("A", "C", "G", "T", "E")])
    }
    
    New.motif <- t(New.motif)
    rownames(New.motif) <- c("A", "C", "G", "T", "E")
    ################## write modified motif######################
    New.motif <- t(New.motif)
    csl2 <- make_col_scheme(chars = c("A", "C", "G", "T", "E"), group = c("1", "2", "3", "4", "5"), cols = c("#08d61d", "#0d09ed", 
        "#edb409", "#ed0909", "#cc00ff"))
    gp <- ggplot() + geom_logo(New.motif, method = "prob", col_scheme = csl2, show_guide = FALSE) + labs(y = "") + theme(panel.background = element_rect(fill = "transparent", 
        color = NA)) + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank())
    ggsave(output_modifed_logo, gp, height = 4, width = 20)
}

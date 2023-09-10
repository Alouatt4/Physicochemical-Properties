####################################################################
########### Physicochemical Properties of FASTA sequence ###########
####################################################################

install.packages("seqinr")
install.packages("ggplot2")
install.packages("tidyverse")
library("tidyverse")
library("ggplot2")
library("seqinr")

#read FASTA seq
dna1 <- read.fasta(file="BTER_mrjp.fa")  
dna2 <- read.fasta(file="BTER_y-e.fa")  
dna3 <- read.fasta(file="BTER_y-e3.fa") 
dna4 <- read.fasta(file="BTER_y-g.fa")  
dna5 <- read.fasta(file="BTER_y-g2.fa") 
dna6 <- read.fasta(file="BTER_y-h.fa")  

#vector FASTA seq
dna1.seq <- dna1[[1]]                           
dna2.seq <- dna2[[1]]                           
dna3.seq <- dna3[[1]]                           
dna4.seq <- dna4[[1]]                           
dna5.seq <- dna5[[1]]                           
dna6.seq <- dna6[[1]]                           

#translate
amino1 <- translate(dna1.seq, sens = F, ambiguous = F, frame = 0) 
amino1                                                
amino2 <- translate(dna2.seq, sens = F, ambiguous = F, frame = 0) 
amino2                                                
amino3 <- translate(dna3.seq, sens = F, ambiguous = F, frame = 0) 
amino3                                                
amino4 <- translate(dna4.seq, sens = F, ambiguous = F, frame = 0) 
amino4                                                
amino5 <- translate(dna5.seq, sens = F, ambiguous = F, frame = 0) 
amino5                                                
amino6 <- translate(dna6.seq, sens = F, ambiguous = F, frame = 0) 
amino6                                                


amino.s <- aaa(amino1)                                 
amino.s
summary(amino.s)

# amino.3 <- count(dna1.seq, start = 0, 3)              
# amino.3
# amino.3 <- count(dna1.seq, start = 1, 3)              
# amino.3
# amino.3 <- count(dna1.seq, start = 2, 3)              
# amino.3
# nucleo <- count(dna1.seq, 2)                          
# nucleo

x11()

#Physicolchemical properties
AAstat(amino1, plot = T)                       
txt1 <- AAstat(amino1)
as.data.frame(txt1)
txt1 <- write.table(txt1,"BTER_mrjp.txt", sep="\t", row.names=T)

AAstat(amino2, plot = T)                       
txt2 <- AAstat(amino2)
as.data.frame(txt2)
txt2 <- write.table(txt2,"BTER_y-e.fa.txt", sep="\t", row.names=T)

AAstat(amino3, plot = T)                       
txt3 <- AAstat(amino3)
as.data.frame(txt3)
txt3 <- write.table(txt3,"BTER_y-e3.fa.txt", sep="\t", row.names=T)

AAstat(amino4, plot = T)                       
txt4 <- AAstat(amino4)
as.data.frame(txt4)
txt4 <- write.table(txt4,"BTER_y-g.fa.txt", sep="\t", row.names=T)

AAstat(amino5, plot = T)                       
txt5 <- AAstat(amino5)
as.data.frame(txt5)
txt5 <- write.table(txt5,"BTER_y-g2.fa.txt", sep="\t", row.names=T)

AAstat(amino6, plot = T)                       
txt6 <- AAstat(amino6)
as.data.frame(txt6)
txt6 <- write.table(txt6,"BTER_y-h.fa.txt", sep="\t", row.names=T)


###########################################
################### end ###################
###########################################





# ggsave("BTER_mrjp.tiff", units = "in", width=10, height=9, dpi=300, compression = 'lzw')
# ggsave("BTER_y-e.tiff", units="in", width=1300, height=700, compression = 'lzw')
# ggsave("BTER_y-e3.tiff", units="in", width=10, height=9, dpi=300, compression = 'lzw')
# ggsave("BTER_y-g.tiff", units="in", width=10, height=9, dpi=300, compression = 'lzw')
# ggsave("BTER_y-g2.tiff", units="in", width=10, height=9, dpi=300, compression = 'lzw')
# ggsave("BTER_y-h.tiff", units="in", width=10, height=9, dpi=300, compression = 'lzw')
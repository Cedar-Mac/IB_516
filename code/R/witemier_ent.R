
#program to calculate entropy of the different codon positions of the MOTUs 
#prepared to read directories with varying values of alpha
#mock directories for two alpha values (2 and 5) provided

library(entropy)
library(stringr)
library(here)

ent<-function(dada,v)
{
  print(paste("processing alpha",alpha[g],"MOTU ",str_remove(fastas[i],".denoise_summary.csv")))
  #assumes the first column contains an ID, customize if necessary
  reads<-rowSums(dada[,3:(dim(dada)[2]-2)])
  seqs<-dada$sequences
  long<-nchar(seqs[1])
  seqqs<-matrix("character",length(seqs),long)
  for (j in 1:length(seqs)) seqqs[j,]<-substring(seqs[j], seq(1, nchar(seqs[j]), 1), seq(1, nchar(seqs[j]), 1))
  
  entro<-vector("numeric",long)
  for (j in 1:long) entro[j]<-entropy(tapply(reads,seqqs[,j],sum))
  
  return(entro)
  
}

data_directory_header <- c("witemier_input")
analysis_directory_header <- c("code")

primer_df <- read.csv("witemier_input/PrimerCombine_Codon_start_positions.csv",
                      stringsAsFactors = FALSE)

primer_list <- primer_df$Primer[primer_df$Coding]

### main program ###
for (current_primer in primer_list){
  regex_primer <- str_replace(current_primer, "\\+", "\\\\+")
  codon_start <- primer_df$CodonStart[primer_df$Primer == current_primer]
  main_outfile <- here("witemier_input", paste("means_entropy_",current_primer,".csv",sep=""))
  
  #alpha should be ordered from higher to lower, so that the
  #last value is the more stringent, i.e., the one having
  #less MOTUs in the corresponding directory
  
  
  alpha<-c(1, 3, 5, 7, 9, 11, 13) #variable alpha
  
  seq_length <- primer_df$Length[primer_df$Primer == current_primer] #length of sequences, needs customization
  
  #generate a file with results, that will be placed in the same directory as the script
  write(paste("alpha","nMOTUs","nseqs","entropy1","entropy2","entropy3","entropy2/3",sep=","),
        file=main_outfile)
  
  #the first nucleotide of the codon is given in PrimerComine_Codon_start_positions.csv
  pos1<-seq(codon_start,seq_length,by=3)
  pos2<-seq((codon_start%%3)+1,seq_length,by=3)
  pos3<-seq(((codon_start+1)%%3)+1,seq_length,by=3)
  
  ### main iteration ###
  
  #To perform the right comparison, we need to consider only the MOTUs present at all alpha values
  #To do so, we consider only MOTUs present in the directory with the lowest alpha, the most stringent
  #(i.e, the directory corresponding to the last value of the variable)
  
  #if you prefer not to follow this advice and want to consider all MOTUs, then comment the following lines
 #  g<-length(alpha)
 # stringent_directori<-here("witemier_input",
  #                          paste("Denoise_summary_alpha_",as.character(alpha[g]),sep=""))#customize directory name as needed
 # fastasreduced<- list.files(path=stringent_directori, pattern=paste("*",regex_primer,".denoise_summary.csv",sep=""))
  #end of lines to comment
  
  for (g in 1:length(alpha))
  {
    
    message("processing alpha ",alpha[g])
    
    alpha_dir_name <- paste("maxee0_5_alpha_", as.character(alpha[g]), sep="")
    alpha_entropy_dir_name <- paste("Alpha_",alpha[g],"_entropy_summaries",sep="")
    alpha_directori<-here("witemier_input", alpha_dir_name) #customize name of directory as needed
    
    #generate output file in anaysis output directory
    alpha_primer_outfile <- here("witemier_input",
                                 alpha_entropy_dir_name, paste(current_primer, "_entropy.csv", sep=""))
    write(paste("id","nseqs","entropy1","entropy2","entropy3",sep=","),alpha_primer_outfile)
    
    #reading files
    fastas <- list.files(path=alpha_directori)
    #select those MOTUs present in the most stringent directory (comment this line if you want to use all MOTUs)
    #fastas<-fastas[fastas%in%fastasreduced]
    
    n_seqs<-vector("integer",length(fastas))
    entropy1<-vector("numeric",length(fastas))
    entropy2<-vector("numeric",length(fastas))
    entropy3<-vector("numeric",length(fastas))
    
    for (i in 1:length(fastas))
    {
      dades<-read.csv(here("witemier_input", alpha_dir_name, fastas[i]),stringsAsFactors=F)
      if(dim(dades)[1] > 0){
        result<-ent(dades,i)
        
        n_seqs[i]<-dim(dades)[1]
        
        entropy1[i]<-mean(result[pos1])
        entropy2[i]<-mean(result[pos2])
        entropy3[i]<-mean(result[pos3])
        
        #fill output file
        write(paste(str_remove(fastas[i],".denoise_summary.csv"),n_seqs[i],entropy1[i],entropy2[i],entropy3[i],sep=","),alpha_primer_outfile,append=T)
      } #end file analysis >0 rows
    } #end i
    
    #fill results file
    
    secs<-sum(n_seqs)
    
    meanent1<-mean(entropy1)
    meanent2<-mean(entropy2)
    meanent3<-mean(entropy3)
    entropy_ratio<-meanent2/meanent3
    
    write(paste(alpha[g],length(fastas),secs,meanent1,meanent2,meanent3,entropy_ratio,sep=","),main_outfile,append=T)
    
  } #end g
  
  message("done")
}



#More details can be found in the README.md file
#The function takes takes the following arguments:
#g: name of the genetype file in plink format (PED/MAP). Only the root name of the file should be specified, and both the .ped and .map 
#files should be in the wrking directory
#tp: an intiger specifying the the number of top principal components (PC) to include analysis
#ts: an intiger spcifying the number of top SNP to pick
#simultaneous: this is a logical(T/F) indicating if simultaneous correlation with top PC (specified by ts) or just the one pc specified
#by'ts' should be used to rank SNP. 
#Dependencies: the function makes use of PLINK1.9, so make sure plink.exe (version for windows commandline) is in your working directory.
#It helps to know how PLINK works and how output files are structured to undertand why some parts of this function are the way they are.

topsnp<- function (g, tp = 1, ts, simultaneous = F) {
  #perform PCA on data, after doing some quality control
  system(paste("plink","-cow", "-file", g, "-pca","-maf 0.1", "-geno 0.1", "-mind 0.1", "-out", "pca")) 
  ev<-read.table("pca.eigenvec", header = F) #read eigenvectors into ev
  eva<-read.table("pca.eigenval", header = F) #read eigenvalues into eva
  
  #if 'simultaneous' is false, take the 'ts', do GWA with it as a response variable and save output in p which has two colums: animal IDs and -log(p-value)
    if (!simultaneous){
    write.table(ev[, c(1, 2, 2+tp)], "pheno.txt", row.names = F, col.names = F, 
                quote = F, sep = " ") # writing a phenotype file to input into plink
    tname<-paste("asso", ts, sep="") # making a temporary name for plink output file
    system(paste("plink","-cow", "-file", "250k", "-assoc", "-pheno", "pheno.txt", 
                 "-allow-no-sex", "-out", tname)) #performing GWA with selected PC as a response and putting result in tname
    p<- read.table(paste(tname, ".qassoc", sep = ""), header = T)[,c(2,9)] # reading relevant parts of GWA result (SNP id & p-values) 
    p[, 2]<- log10(p[, 2]) * (-1)  #converting p-values into -logs
    } 
    else 
    {
    #if 'simultaneous' is true Perform GWA for 'ts' number of top PC and put results in the dataframe p whose 1st colmn is SNP ID, and subsequent columns are -log(p-value) for each PC
    marker_lists<-lapply(1:tp, function(x) {
    write.table(ev[, c(1, 2, 2+x)], "pheno.txt", row.names = F, col.names = F, 
                quote = F, sep = " ") # writing a phenotype ("pheno") file to input into plink
    tname<-paste("asso", x, sep = "") # making a temporary name for plink output file
    system(paste("plink","-cow", "-file", g, "-assoc", "-pheno", "pheno.txt", 
                 "-allow-no-sex", "-out", tname)) #performing GWA with selected PC as a response and putting result in tname
    return(read.table(paste(tname, ".qassoc", sep = ""), header = T)[,c(2,9)]) # relevant parts of GWA result (SNP id & p-values) 
      })
    p<- as.data.frame(marker_lists) #changing to dataframe
    p<- p[,-(seq(3, ncol(p), by = 2))] #get rid of redundant columns (snp IDs)
    p[, 2:(tp+1)]<- log10(p[,2:(tp+1)]) * (-1) #change raw p-values into -logs
    }
    
   #if 'simultaneous' is true, use the GWA result for the top 'ts' PC  and calculate linnear combination of -log(p) and eigenvalues. else use -log(p-val) for 1 PC only to rank SNP
   if (simultaneous)
   {
   y<- sapply(1:nrow(p), function (x) {  
   sum((p[x, 2:(tp + 1)]) * eva[,1][1:tp])
   })   # make 'y', a value to rank SNP based on -log(p-val) weighted by eigenvalue (eva) for corresponding PC
   p<- data.frame(p[,1], y)   #z is a data frame with marker names and their associated y value  
   }
  
   p<- p[order(p[,2], decreasing = T), ] #reordering rows of 'p' based on values of y which is the -log(p-value) derivation
   
  topsnp<- p[1:ts, 1]    # extract 'ts' amount of top SNP ranked based on y
  write.table(topsnp, file = "topsnp.txt", 
              row.names = F, col.names = F, sep = " ", quote = F)   #write top selected markers into a file for plink to read
  
  return(topsnp)
}

The function returns a vector with a specified number of top SNP selected based on correlation with 
the specified number of principal components which is determided by performing GWA with the PC as response variable. 
When using only one PC, the negative log of the p-value is used for ranking. When correlation with more than PC is preffered,  
the sum of the the -log (p-value) of the PC weighted by their eigenvalue is used for ranking (equation shown below).
The function also writes them into a text file that can be read by plink.

It takes the following arguments:

1. g: name of the genetype file in plink format (PED/MAP). Only the root name of the file should be specified, and both
the .ped and .map files should be in the wrking directory

2. tp: an intiger specifying the the number of top SNP to consider 

3. ts: an intiger spcifying the number of top SNP to pick 

4. simultaneous: this is a logical indicating if simultaneous correlation with top PC (specified by ts) or just the one pc 
specified by 'ts' should be used to rank SNP. 
Dependencies: the function makes use of PLINK1.9, so make sure plink.exe (for windows commandline) is in your working directory
It helps to know how PLINK works and how output files are structured to undertand why some parts of this function are the way they are.


The following equation was used to derive the value used to rank the SNP if correlation with more than PC is preffered:

xj = [-log(pc_p1 ) * ev1] + [-log(pc_p2 ) * ev2] … + [-log(pc_pi) * evi]
  pc_pi is the p-value from GWA for  ith PC and jth SNP
  evi   is the eigenvalue of ith PC. The idea is to find SNP with relation to the top PC while accounting for eigenvalue (amount of variation or population structure explained).

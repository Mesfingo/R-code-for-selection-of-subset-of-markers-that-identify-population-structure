# The following equation was used to derive the value used to rank the SNP:
xj = [-log(pc_p1 ) * ev1] + [-log(pc_p2 ) * ev2] … + [-log(pc_pi) * evi]
  pc_pi is the p-value from GWA for  ith PC and jth SNP
  evi   is the eigenvalue of ith PC

#The idea is to find SNP with relation to the top PC while accounting for eigenvalue.

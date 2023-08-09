library(SNPassoc)
genotype_data = read.csv('C:/Stas/LabWork/Bioinformatics/Projects/Ch4_GWAS/epigen/sim/Pure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02_SNPassoc.csv', header = TRUE, sep='\t')
genotype_data_s = setupSNP(data=genotype_data, colSNPs=3:ncol(genotype_data), sep="")



i = 0
for (x in colnames(genotype_data_s)){
  if(i > 2){
    print(x)
    ans = association(Phenotype ~ genotype_data_s[,c(x)], data = genotype_data_s)
    #print(ans)
  }
  i = i + 1
}
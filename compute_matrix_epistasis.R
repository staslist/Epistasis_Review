library("MatrixEpistasis")

compute_matrix_epistasis <- function(fname) {
  fname1 = paste(fname, sep = '', '.csv')
  fname2 = paste(fname, sep = '', '.pheno')
  fname3 = paste(fname, sep = '', '_MatrixEpistasis.txt')
  
  data = read.csv(fname1)
  snpA = as.matrix(data)
  snpB = as.matrix(data)
  trait = read.csv(fname2, header = FALSE)
  MatrixEpistasis_main(snpA, snpB, trait=trait, pvalThreshold=5e-8, outputFileName=fname3)
}
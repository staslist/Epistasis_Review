#include <mkl.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <time.h>
#include <omp.h>

int remma_epi_c(long num_ID, long num_SNP, long start_SNP,
long end_SNP, long total_part, long which_part, long num_test_each, char* kin_file,
char* output_prefix, double chi_value){
	
	time_t t;
	
	//////
	printf("Check the paramters from python to C.\n");
	printf("Individuals' number: %ld\n", num_ID);
	printf("SNP number: %ld\n", num_SNP);
	printf("Start and end for the first SNP: %ld, %ld\n", start_SNP, end_SNP);
	printf("The number of test in memory: %ld\n", num_test_each);
	printf("Kinship file prefix: %s\n", kin_file);
	printf("The prefix for output_file: %s\n", output_prefix);
	printf("Threshold chi value: %g\n", chi_value);
	
	
	
	//////
	printf("Prepare the SNP information.\n");
	
	long i = 0;
    char **chro = (char**) calloc(num_SNP, sizeof(char*));
    for(i = 0; i < num_SNP; i++){
    	chro[i] = (char*) calloc(100, sizeof(char));
    }
    
    
    char **SNP_ID = (char**) calloc(num_SNP, sizeof(char*));
    for(i = 0; i < num_SNP; i++){
    	SNP_ID[i] = (char*) calloc(100, sizeof(char));
    }
	
	  
    char **pos = (char**) calloc(num_SNP, sizeof(char*));
    for(i = 0; i < num_SNP; i++){
    	pos[i] = (char*) calloc(100, sizeof(char));
    }
    
    char bed_file[256];
    strcpy(bed_file, kin_file);
    strcat(bed_file, ".bim");
    
    FILE *in_bed = fopen(bed_file, "r");
    if(in_bed == NULL){
    	printf("Fail to open the bed file: %s.\n",bed_file);
    	exit(1);
    }
    
    i = 0;
    while(fscanf(in_bed,"%s%s%*s%s%*s%*s",chro[i],SNP_ID[i],pos[i])==3){
    	i++;
    }
    fclose(in_bed);
    in_bed = NULL;
	
		
	//////
	printf("Read temp file from python.\n");
	
	//eigen vector
	double *eigen_vec = (double*) calloc(num_ID*num_ID, sizeof(double));
	if(eigen_vec == NULL){
		printf("Not enough memory for eigen vector.\n");
		exit(1);
	}
	
	char eigen_vec_file[256];
	strcpy(eigen_vec_file, output_prefix);
	strcat(eigen_vec_file, ".eigen_vec");
	
	FILE *in_eigen_vec = fopen(eigen_vec_file, "r");
	if(in_eigen_vec==NULL){
		printf("Fail to open the eigen vector file.\n");
		exit(1);
	}
	
	i = 0;
	while(fscanf(in_eigen_vec,"%lf",&eigen_vec[i])==1){
		i++;
	}
	
	fclose(in_eigen_vec);
	
	//eigen value
	double *eigen_val = (double*) calloc(num_ID, sizeof(double));
	if(eigen_val == NULL){
		printf("Not enough memory for eigen value.\n");
		exit(1);
	}
	
	char eigen_val_file[256];
	strcpy(eigen_val_file, output_prefix);
	strcat(eigen_val_file, ".eigen_val");
	
	FILE *in_eigen_val = fopen(eigen_val_file, "r");
	if(in_eigen_val == NULL){
		printf("Fail to open the eigen value file.\n");
		exit(1);
	}
	
	i = 0;
	while(fscanf(in_eigen_val,"%lf",&eigen_val[i])==1){
		i++;
	}
	
	fclose(in_eigen_val);
		
	//temp effect
	double *temp_effect = (double*) calloc(num_ID, sizeof(double));
	if(temp_effect == NULL){
		printf("Not enough memory for temp effect.\n");
		exit(1);
	}
	
	char temp_effect_file[256];
	strcpy(temp_effect_file, output_prefix);
	strcat(temp_effect_file, ".eff");
	
	FILE *in_temp_effect = fopen(temp_effect_file, "r");
	if(in_temp_effect == NULL){
		printf("Fail to open the temp effect file.\n");
		exit(1);
	}
	
	i = 0;
	while(fscanf(in_temp_effect,"%lf",&temp_effect[i])==1){
		i++;
	}
	
	fclose(in_temp_effect);
	
	
	//////
	printf("Read the normalized marker matrix.\n");
	
	double *norm_SNP_mat = (double*) calloc(num_ID*num_SNP, sizeof(double));
	if(norm_SNP_mat == NULL){
		printf("Not enough memory for normlized SNP matrix.\n");
		exit(1);
	}
	
	char norm_SNP_file[256];
	strcpy(norm_SNP_file, kin_file);
	strcat(norm_SNP_file, ".addMarkerMat.bin");
	
	FILE *in_norm_SNP = fopen(norm_SNP_file, "rb");
	if(in_norm_SNP == NULL){
		printf("Fail to open the normalized matrix file.\n");
		exit(1);
	}
	
	i = 0;
	while(fread(&norm_SNP_mat[i], sizeof(double), 1, in_norm_SNP) == 1){
		i++; 
	}
	
	fclose(in_norm_SNP);
	in_norm_SNP = NULL;
	
	
	//////
	printf("Start epistatic test.\n");
	double *temp_epi_mat = (double*) calloc(num_ID * num_test_each, sizeof(double));
	if(temp_epi_mat == NULL){
		printf("Not enough memory for temp epistatic marker matrix, please use smaller value for num_test_each.\n");
		exit(1);
	}
	
	double *temp_epi_mat2 = (double*) calloc(num_ID * num_test_each, sizeof(double));
	if(temp_epi_mat2 == NULL){
		printf("Not enough memory for temp epistatic marker matrix, please use smaller value for num_test_each.\n");
		exit(1);
	}
	
	
	double *epi_effect = (double*) calloc(num_test_each, sizeof(double));
	if(epi_effect == NULL){
		printf("Not enough memory for epistatic effect.\n");
		exit(1);
	}

	double *epi_effect_var = (double*) calloc(num_test_each, sizeof(double));
	if(epi_effect_var == NULL){
		printf("Not enough memory for variance of epistatic effect.\n");
		exit(1);
	}	
	
	long j = 0, k = 0, m = 0, n = 0, p = 0, q = 0;
	long nth_set_test = 0;
	long SNP1 = start_SNP, SNP2 = start_SNP;
	
	double SNP_chi = 0.0;
	double alpha = 1.0, beta = 0.0;
	int incx = 1, incy = 1;
	
	
	char output_file[256];
	char suffix[256];
	sprintf(suffix, ".epires.%ld-%ld", total_part, which_part);
	strcpy(output_file, output_prefix);
	strcat(output_file, suffix);
	FILE *out_res = fopen(output_file, "w");
	if(out_res==NULL){
		printf("Fail to build the output file.\n");
		exit(1);
	}
	
	fprintf(out_res,"order1 order2 chro1 SNP_ID1 pos1 chro2 SNP_ID2 pos2 epi_effect chi_value\n");
	
	//long x = 0;
	printf("%ld %ld %ld \n", start_SNP, end_SNP, num_SNP);
	for(i = start_SNP; i < end_SNP; i++){
		for(j = i + 1; j < num_SNP; j++){
			//x++;
			//printf("The %ld th test.\n", x);
			
			for(n = 0; n < num_ID; n++){
				temp_epi_mat[m*num_ID + n] = norm_SNP_mat[i*num_ID + n] * norm_SNP_mat[j*num_ID + n];
			}
			m++;
			
			if(m >= num_test_each){
				printf("%ld\n", m);
				m = 0;
				nth_set_test++;
				time(&t);
				printf("##%s##\n",ctime(&t));
				printf("The %ld th subset SNPs for test.\n", nth_set_test);
				
				//epi effect
				cblas_dgemv(CblasRowMajor, CblasNoTrans,
				num_test_each, num_ID, alpha, temp_epi_mat, num_ID, temp_effect,
				incx, beta, epi_effect, incy);
				
				//var of epi effect
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_test_each,
				num_ID, num_ID, alpha, temp_epi_mat, num_ID, eigen_vec, num_ID,
				beta, temp_epi_mat2, num_ID);
				
			    #pragma omp parallel for
				for(p = 0; p < num_test_each; p++){
					epi_effect_var[p] = 0.0;
					for(q = 0; q < num_ID; q++){
						epi_effect_var[p] += temp_epi_mat2[p*num_ID + q]*eigen_val[q]*temp_epi_mat2[p*num_ID + q];
					}
					
				}
				
				for(p = 0; p < num_test_each; p++){
					
					SNP2++;
					if(SNP2 > num_SNP - 1){
						SNP1++;
						SNP2 = SNP1 + 1;
					}
					
					SNP_chi = epi_effect[p]*epi_effect[p]/epi_effect_var[p];
					
					if(SNP_chi > chi_value){
						fprintf(out_res,"%ld %ld %s %s %s %s %s %s %g %g\n",SNP1+1, SNP2+1, chro[SNP1],
						SNP_ID[SNP1],pos[SNP1],chro[SNP2],SNP_ID[SNP2],pos[SNP2],
						epi_effect[p],SNP_chi);
					}
				}
				
			}
			
		}
		
	}
	free(norm_SNP_mat);
	norm_SNP_mat = NULL;
	
    time(&t);
    printf("##%s##\n",ctime(&t)); 
	printf("Test for the remaining epistatic SNPs.\n");
	long mod = 0;
	
	mod = ((num_SNP - start_SNP - 1 + num_SNP - end_SNP)*(end_SNP - start_SNP)/2)%num_test_each;
	
	if(mod != 0){
		//epi effect
		cblas_dgemv(CblasRowMajor, CblasNoTrans,
		mod, num_ID, alpha, temp_epi_mat, num_ID, temp_effect,
		incx, beta, epi_effect, incy);
		
		//var of epi effect
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mod,
		num_ID, num_ID, alpha, temp_epi_mat, num_ID, eigen_vec, num_ID,
		beta, temp_epi_mat2, num_ID);
		
		
		#pragma omp parallel for
		for(p = 0; p < mod; p++){
			epi_effect_var[p] = 0.0;		
			for(q = 0; q < num_ID; q++){
				epi_effect_var[p] += temp_epi_mat2[p*num_ID + q]*eigen_val[q]*temp_epi_mat2[p*num_ID + q];
			}
					
		}
				
		for(p = 0; p < mod; p++){
					
			SNP2++;
			if(SNP2 > num_SNP - 1){
				SNP1++;
				SNP2 = SNP1 + 1;
			}
					
			SNP_chi = epi_effect[p]*epi_effect[p]/epi_effect_var[p];
					
			if(SNP_chi > chi_value){
				fprintf(out_res,"%ld %ld %s %s %s %s %s %s %g %g\n", SNP1+1, SNP2+1, chro[SNP1],
				SNP_ID[SNP1],pos[SNP1],chro[SNP2],SNP_ID[SNP2],pos[SNP2],
				epi_effect[p],SNP_chi);
			}
		}

		
		
	}
	
	free(epi_effect);
	epi_effect = NULL;
	free(epi_effect_var);
	epi_effect_var = NULL;
	free(temp_epi_mat);
	temp_epi_mat = NULL;
	free(temp_epi_mat2);
	temp_epi_mat2 = NULL;	
	free(temp_effect);
	temp_effect = NULL;
	fclose(out_res);
	out_res = NULL;
	
    time(&t);
    printf("##%s##\n",ctime(&t));  
	printf("Finished.\n");  
	
	return 0;		
		
}
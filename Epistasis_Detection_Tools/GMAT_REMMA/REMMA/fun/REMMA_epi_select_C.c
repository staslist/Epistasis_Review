#include <mkl.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <time.h>
#include <omp.h>

int remma_epi_select_c(long num_ID, long num_SNP, long num_SNP_select, long num_test_each, 
char* plink_prefix, char* output_prefix, char* SNP_select_file,
double chi_value){
	
	time_t t;
	
	//////
	printf("Check the paramters from python to C.\n");
	printf("Individuals' number: %ld\n", num_ID);
	printf("SNP number: %ld\n", num_SNP);
	printf("Selected SNP pair number: %ld.\n", num_SNP_select);
	printf("The number of test in memory: %ld\n", num_test_each);
	printf("The plink file prefix: %s\n", plink_prefix);
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
    strcpy(bed_file, plink_prefix);
    strcat(bed_file, ".bim");
    
    FILE *in_bed = fopen(bed_file, "r");
    if(in_bed == NULL){
    	printf("Fail to open the bed file: %s.\n",bed_file);
    	exit(1);
    }
    
    i = 0;
    while(fscanf(in_bed,"%s%s%*s%s%*s%*s",chro[i],SNP_ID[i],pos[i]) == 3){
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
	printf("Read the SNP select file.\n");
	long *SNP_select_pair1 = (long*) calloc(num_SNP_select, sizeof(long));
	if(SNP_select_pair1 == NULL){
		printf("Not enough memory for selected SNP pair.\n");
		exit(1);
	}

	long *SNP_select_pair2 = (long*) calloc(num_SNP_select, sizeof(long));
	if(SNP_select_pair2 == NULL){
		printf("Not enough memory for selected SNP pair.\n");
		exit(1);
	}
	
	//char SNP_select_file[256];
	//strcpy(SNP_select_file, output_prefix);
	//strcat(SNP_select_file, ".random_SNP");
	
	FILE *in_SNP_select = fopen(SNP_select_file, "r");
	if(in_SNP_select == NULL){
		printf("Fail to open the SNP select pair file.\n");
		exit(1);
	}
	
	fscanf(in_SNP_select, "%*[^\n]%*c"); //header line
	i = 0;
	while(fscanf(in_SNP_select, "%ld%ld", &SNP_select_pair1[i], &SNP_select_pair2[i]) == 2){
		fscanf(in_SNP_select, "%*[^\n]%*c");
		i++;
	}
	
	//////
	printf("Read the normalized marker matrix.\n");
	
	double *norm_SNP_mat = (double*) calloc(num_ID*num_SNP, sizeof(double));
	if(norm_SNP_mat == NULL){
		printf("Not enough memory for normlized SNP matrix.\n");
		exit(1);
	}
	
	char norm_SNP_file[256];
	strcpy(norm_SNP_file, plink_prefix);
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
	long pair1 = 0, pair2 = 0;
	long SNP1 = -1, SNP2 = -1;
	long SNP1_info = 0, SNP2_info = 0;
	
	double SNP_chi = 0.0;
	double alpha = 1.0, beta = 0.0;
	int incx = 1, incy = 1;
	
	char output_file[256];
	strcpy(output_file, SNP_select_file);
	strcat(output_file, "_res");
	FILE *out_res = fopen(output_file, "w");
	if(out_res==NULL){
		printf("Fail to build the output file.\n");
		exit(1);
	}
	
	fprintf(out_res,"order1 order2 chro1 SNP_ID1 pos1 chro2 SNP_ID2 pos2 epi_effect chi_value\n");
	
	for(i = 0; i < num_SNP_select; i++){
			
			pair1 = SNP_select_pair1[i] - 1;
			pair2 = SNP_select_pair2[i] - 1;
			for(n = 0; n < num_ID; n++){
				temp_epi_mat[m*num_ID + n] = norm_SNP_mat[pair1*num_ID + n] * norm_SNP_mat[pair2*num_ID + n];
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
					
					SNP1++;
					SNP2++;
					
					SNP1_info = SNP_select_pair1[SNP1] - 1;
					SNP2_info = SNP_select_pair2[SNP2] - 1;
					
					SNP_chi = epi_effect[p]*epi_effect[p]/epi_effect_var[p];
					
					if(SNP_chi > chi_value){
						fprintf(out_res,"%ld %ld %s %s %s %s %s %s %g %g\n",SNP1_info+1, SNP2_info+1, chro[SNP1_info],
						SNP_ID[SNP1_info],pos[SNP1_info],chro[SNP2_info],SNP_ID[SNP2_info],pos[SNP2_info],
						epi_effect[p],SNP_chi);
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
	
	mod = num_SNP_select%num_test_each;
	
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
					
			SNP1++;
			SNP2++;
					
			SNP1_info = SNP_select_pair1[SNP1] - 1;
			SNP2_info = SNP_select_pair2[SNP2] - 1;
					
			SNP_chi = epi_effect[p]*epi_effect[p]/epi_effect_var[p];
					
			if(SNP_chi > chi_value){
				fprintf(out_res,"%ld %ld %s %s %s %s %s %s %g %g\n", SNP1_info+1, SNP2_info+1, chro[SNP1_info],
				SNP_ID[SNP1_info],pos[SNP1_info],chro[SNP2_info],SNP_ID[SNP2_info],pos[SNP2_info],
				epi_effect[p],SNP_chi);
			}
		}	
	}
	
	free(SNP_select_pair1);
	SNP_select_pair1 = NULL;
	free(SNP_select_pair2);
	SNP_select_pair2 = NULL;	
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
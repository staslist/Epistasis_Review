#include <mkl.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <time.h>
#include <omp.h>

int remma_epi_scan_c(long num_ID, long num_SNP, long start_SNP,
long end_SNP, long total_part, long which_part, char* kin_file, 
char* output_prefix, double epi_effect_threshold){
	
	time_t t;
	
	//////
	printf("Check the paramters from python to C.\n");
	printf("Individuals' number: %ld\n", num_ID);
	printf("SNP number: %ld\n", num_SNP);
	printf("Start and end for the first SNP: %ld, %ld\n", start_SNP, end_SNP);
	printf("Kinship file prefix: %s\n", kin_file);
	printf("The prefix for output file: %s\n", output_prefix);
	printf("Threshold chi value: %g\n", epi_effect_threshold);
	
	
	
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
		
	long j = 0, k = 0, m = 0, n = 0, p = 0, q = 0;
	long nth_set_test = 0;
	long SNP1 = start_SNP, SNP2 = start_SNP;
	
	double SNP_chi = 0.0;
	double alpha = 1.0, beta = 0.0;
	int incx = 1, incy = 1;
	
	char output_file[256];
	char suffix[256];
	sprintf(suffix, ".epieff.%ld-%ld", total_part, which_part);
	strcpy(output_file, output_prefix);
	strcat(output_file, suffix);	
	
	FILE *out_res = fopen(output_file, "w");
	if(out_res==NULL){
		printf("Fail to build the output file.\n");
		exit(1);
	}
	
	fprintf(out_res,"order1 order2 chro1 SNP_ID1 pos1 chro2 SNP_ID2 pos2 epi_effect\n");
	
	
	long median = (long)((start_SNP + end_SNP - 1)/2);
	if((start_SNP + end_SNP - 1)%2 == 0){
		median = median - 1;
	}
	
	
    int printOut(long i, long start_SNP, long end_SNP, long num_SNP, long num_ID,
	double *norm_SNP_mat, double *temp_effect, double epi_effect_threshold,
	FILE *out_res, char **chro, char **SNP_ID,char **pos);
	
	#pragma omp parallel for
    for(i = start_SNP; i <= median; i++){
		
		printOut(i, start_SNP, end_SNP, num_SNP, num_ID,
		norm_SNP_mat, temp_effect, epi_effect_threshold,
		out_res, chro, SNP_ID, pos)	;	    
	    
    }
    
    
    double epi_effect = 0;
    i = median + 1;
    if((end_SNP + start_SNP	- 1)%2 == 0){
    	
    	for(j = i + 1; j < num_SNP; j++){
	    	epi_effect = 0;
	    	for(k =0 ; k < num_ID; k++){
	    		epi_effect += norm_SNP_mat[i*num_ID + k] * norm_SNP_mat[j*num_ID + k] * temp_effect[k];
    		}
	    	
            if(epi_effect >= epi_effect_threshold || epi_effect	<= -epi_effect_threshold){
                fprintf(out_res, "%ld %ld %s %s %s %s %s %s %g\n", i+1, j+1, 
				chro[i], SNP_ID[i], pos[i], chro[j], SNP_ID[j], pos[j], epi_effect);
            }				   	
	    } 
		   	
    }
    		
}



int printOut(long i, long start_SNP, long end_SNP, long num_SNP, long num_ID,
double *norm_SNP_mat, double *temp_effect, double epi_effect_threshold,
FILE *out_res, char **chro, char **SNP_ID, char **pos){
		
		double epi_effect = 0;
		long j, p, q, k;
		for(j = i + 1; j < num_SNP; j++){
			
	    	epi_effect = 0;	    	
	    	for(k = 0; k < num_ID; k++){
	    		epi_effect += norm_SNP_mat[i*num_ID + k] * norm_SNP_mat[j*num_ID + k] * temp_effect[k];
	    	}
	    	
            if(epi_effect >= epi_effect_threshold || epi_effect	<= -epi_effect_threshold){
                fprintf(out_res,"%ld %ld %s %s %s %s %s %s %g\n",i + 1, j + 1, chro[i], SNP_ID[i], 
				pos[i], chro[j], SNP_ID[j], pos[j], epi_effect);
            }
							   	
	    }
	    
	    p = end_SNP + start_SNP - 1 - i;
    	for(q = end_SNP	- i + start_SNP; q < num_SNP; q++){
    		
	    	epi_effect = 0;
	    	for(k = 0; k < num_ID; k++){
	    		epi_effect += norm_SNP_mat[p*num_ID + k] * norm_SNP_mat[q*num_ID + k] * temp_effect[k];
	    	}
	    	
            if(epi_effect >= epi_effect_threshold || epi_effect	<= -epi_effect_threshold){
                fprintf(out_res,"%ld %ld %s %s %s %s %s %s %g\n",p + 1, q + 1, chro[p], SNP_ID[p],
				pos[p], chro[q], SNP_ID[q], pos[q], epi_effect);
            }
						   	
	    }
	    
	    return 0;
	
}

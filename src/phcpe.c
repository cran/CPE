/*
Qianxing Mo, moq@mskcc.org
Department of Epidemiology and Biostatistics
Memorial Sloan-Kettering Cancer Center, NY 10021

These C functions are used to compute the concordance probability for the Cox proportional hazard model
(see Gonen and Heller, Biometrika(2005), 92, 4, pp. 965-970)

The algorithm used to compute the probability is based on the R codes written by  Ennapadam. S. Venkatraman at MSKCC, with further optimization by Qianxing Mo.  

/* Last updated 8/21/08 */

#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>

/* x = x - y*/
void vector_sub(int size, double * x, const double *y) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = x[i] - y[i];
  }
}
/* x = x + y */
void vector_add(int size, double * x, const double *y) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = x[i] + y[i];
  }
}
/* x = x*alpha */
void vector_scale(int size, double * x, const double alpha) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = alpha*x[i];
  }
}

/* Function to compute cpe and its variance */
void coxcpe(const int *ROW, const int *COL, const double *bandwidth, const double *xbeta, const double * Design, const double *varbeta, double * result) {
  
  int i, j, loc;
  int incx, incy;
  double ONE,ZERO;
  char *trans = "N";
  double CPE=0, CPEapp=0, varterm1=0, varterm2=0, UU2=0, totalrowsum2=0, tempCPE, tempCPEapp; 
  double bxjxi,bxixj,denomji,denomij,cdfbxij,cdfbxji,pdfbxij,pdfbxji,Uji,Uij,Scale1, Scale2; 
  double *xji, *xij, *tempv, *tempxij, *tempxji, *sumdervec, *rowsum2;
  double ** design;

  incx = 1;
  incy = 1;
  ONE = 1.0;
  ZERO = 0;

  xji = (double *) malloc((*COL)*sizeof(double));
  xij = (double *) malloc((*COL)*sizeof(double));
  tempv = (double *) malloc((*COL)*sizeof(double));
  tempxij = (double *) malloc((*COL)*sizeof(double));
  tempxji = (double *) malloc((*COL)*sizeof(double));
  sumdervec = (double *) calloc((*COL), sizeof(double));
  rowsum2 = (double *) calloc((*ROW), sizeof(double));

  if(xji==NULL || xij==NULL || tempv==NULL || tempxij==NULL || tempxji==NULL || sumdervec==NULL || rowsum2==NULL){
    printf("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
    exit(1);
  }

  design = (double **) malloc((*ROW)*sizeof(double *));
  if(design==NULL){
    printf("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
    exit(1);
  }
  
  for(i=0; i<(*ROW); i++) {
    design[i] = (double *) malloc((*COL)*sizeof(double));
    if(design[i] == NULL){
      printf("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
      exit(1);
    }
    for(j=0; j<(*COL); j++){
      loc = i*(*COL) + j;
      design[i][j] = Design[loc];
    }
  }

  Scale1 = 1.0/(*ROW);
  Scale2 = 2.0/((*ROW)*((*ROW)-1));

  for(i=0; i<((*ROW)-1); i++) {
    UU2 = 0.0;
    tempCPE = 0;
    tempCPEapp = 0;
    for(j=(i+1); j<(*ROW); j++) {
      bxjxi = xbeta[j] - xbeta[i];
      bxixj = 0 - bxjxi;
      denomji = 2 + expm1(bxjxi);
      denomij = 2 + expm1(bxixj);
      cdfbxij = pnorm(bxixj/(*bandwidth), 0, 1, 1, 0);
      cdfbxji = pnorm(bxjxi/(*bandwidth), 0, 1, 1, 0);
      pdfbxij = dnorm(bxixj/(*bandwidth), 0, 1, 0);
      pdfbxji = dnorm(bxjxi/(*bandwidth), 0, 1, 0);

      Uji = cdfbxij/denomji;
      Uij = cdfbxji/denomij;

      rowsum2[i] += Uji + Uij;
      rowsum2[j] += Uji + Uij;

      tempCPE += 1.0*(bxjxi <= 0)/denomji + 1.0*(bxixj < 0)/denomij;  /* [1] Kn(beta.hat) */
      tempCPEapp += Uji + Uij;                  /* smooth approximation for Kn(beta.hat)*/
      UU2 += (Uji + Uij) * (Uji + Uij);
      
      /* The following is for the asymptotic variance, form [5] */
 
      F77_CALL(dcopy)(COL,design[j],&incx, xji,&incy); /* copy design[j] to xji */
      F77_CALL(dcopy)(COL,design[i],&incx, tempv,&incy); /* copy design[i] to tempv */
      vector_sub(*COL, xji, tempv);
      F77_CALL(dcopy)(COL, xji, &incx, xij, &incy); 
      vector_scale(*COL, xij,-1.0);
  
      F77_CALL(dcopy)(COL, xij, &incx, tempxij, &incy); /* copy xij to tempxij */
      vector_scale(*COL, tempxij, pdfbxij/(denomji*(*bandwidth)));
      vector_scale(*COL, xij, cdfbxij*(denomji - 1.0)/(denomji*denomji));
      F77_CALL(dcopy)(COL, xji, &incx, tempxji, &incy);
      vector_scale(*COL, tempxji, pdfbxji/(denomij*(*bandwidth)));
      vector_scale(*COL, xji, cdfbxji*(denomij - 1)/(denomij*denomij));

      vector_add(*COL, tempxij, xij); /* tempxij = tempxij + xij */
      vector_add(*COL, tempxij, tempxji);
      vector_add(*COL, tempxij, xji);
      vector_scale(*COL, tempxij, Scale1);
      vector_add(*COL, sumdervec, tempxij);
     }		
    totalrowsum2 += Scale2*((rowsum2[i] + 0.5) * (rowsum2[i] + 0.5) - 2.0*UU2); /* for i < j */
    CPE += Scale1 * tempCPE;
    CPEapp += Scale1 * tempCPEapp;
  }
  
  /* the last row sum, i = (*ROW) -1 */
  totalrowsum2 += Scale2*((rowsum2[(*ROW) - 1] + 0.5) * (rowsum2[(*ROW) - 1] + 0.5));
  vector_scale(*COL, sumdervec, 2.0/((*ROW)-1));
  F77_CALL(dgemv)(trans, COL, COL, &ONE, varbeta, COL, sumdervec, &incx, &ZERO, tempv, &incy); 
  varterm2 = F77_CALL(ddot)(COL, sumdervec, &incx, tempv, &incy);
  varterm1 = totalrowsum2 - Scale2*(0.25*(*ROW) + 4.0*CPEapp*(2.0*(*ROW)*CPEapp + 0.5*(*ROW)) - 4.0*(*ROW)*(*ROW)*CPEapp*CPEapp/((*ROW)-1));
  varterm1 = 2.0*varterm1/((*ROW)*((*ROW)-1));

  CPE = 2.0*CPE/((*ROW) - 1);
  CPEapp =  2.0*CPEapp/((*ROW) - 1);

  result[0] = CPE;
  result[1] = CPEapp;
  result[2] = sqrt(varterm1 + varterm2);
  printf("Result --- varterm1 = %f, varterm2 = %f\n",varterm1, varterm2);

  for(i = 0; i < (*ROW); i++){ free(design[i]); }
  free(design);
  free(rowsum2);
  free(xij);
  free(xji);
  free(tempv);
  free(tempxij);
  free(tempxji);
  free(sumdervec);
}

/* Function to calculate the CPE only */ 
void coxcpeOnly(const int *ROW, const double *xbeta, double * result){
  int i, j;
  double CPE, tempCPE, bxjxi,bxixj,denomji,denomij,Scale1; 
  
  CPE = 0;
  Scale1 = 1.0/(*ROW);

  for(i=0; i<((*ROW)-1); i++) {
    tempCPE = 0;
    for(j=(i+1); j<(*ROW); j++) {
      bxjxi = xbeta[j] - xbeta[i];
      bxixj = 0 - bxjxi;
      denomji = 2 + expm1(bxjxi);
      denomij = 2 + expm1(bxixj);
      tempCPE += 1.0*(bxjxi <= 0)/denomji + 1.0*(bxixj < 0)/denomij;  /* [1] Kn(beta.hat) */
     }		
    CPE += Scale1 * tempCPE;
  }
  *result = 2.0*CPE/((*ROW) - 1);
}

/* The following two functions added on 8/18/08 to exclude data with ties */

/* Function to calculate the CPE only */ 
void cpeOnlyNoTies(const int *ROW, const double *xbeta, double * result){
  int i,j,N;
  double CPE,tempCPE,bxjxi,bxixj,denomji,denomij,Scale1; 
  
  CPE = 0;
  Scale1 = 1.0/(*ROW);
  /* actual size of no ties of covariates */
  N = 0;
  for(i=0; i<((*ROW)-1); i++) {
    tempCPE = 0;
    for(j=(i+1); j<(*ROW); j++) {
      if(xbeta[j] != xbeta[i]){
	N = N + 1;
	bxjxi = xbeta[j] - xbeta[i];
	bxixj = 0 - bxjxi;
	denomji = 2 + expm1(bxjxi);
	denomij = 2 + expm1(bxixj);
	tempCPE += 1.0*(bxjxi < 0)/denomji + 1.0*(bxixj < 0)/denomij;  /* [1] Kn(beta.hat) */
      }	
    }	
    CPE += Scale1 * tempCPE;
  }
  *result = (*ROW)*CPE/N;
}

/* Function to compute cpe and its variance, tied  */
void cpeNoTies(const int *ROW, const int *COL, const double *bandwidth, const double *xbeta, const double * Design, const double *varbeta, double * result) {
  
  int i, j, loc,N;
  int incx, incy;
  double ONE,ZERO;
  char *trans = "N";
  double CPE=0,CPEapp=0,varterm1=0,varterm2=0,UUKSUM=0,UUKij=0,UUKsumj,UUK,tempCPE,tempCPEapp; 
  double bxjxi,bxixj,denomji,denomij,cdfbxij,cdfbxji,pdfbxij,pdfbxji,Uji,Uij,Scale1; 
  double *xji, *xij, *tempv, *tempxij, *tempxji, *sumdervec;
  double ** design;

  incx = 1;
  incy = 1;
  ONE = 1.0;
  ZERO = 0;
  xji = (double *) malloc((*COL)*sizeof(double));
  xij = (double *) malloc((*COL)*sizeof(double));
  tempv = (double *) malloc((*COL)*sizeof(double));
  tempxij = (double *) malloc((*COL)*sizeof(double));
  tempxji = (double *) malloc((*COL)*sizeof(double));
  sumdervec = (double *) calloc((*COL), sizeof(double));
  if(xji==NULL || xij==NULL || tempv==NULL || tempxij==NULL || tempxji==NULL || sumdervec==NULL){
    printf("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
    exit(1);
  }
  design = (double **) malloc((*ROW)*sizeof(double *));
  if(design==NULL){
    printf("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
    exit(1);
  }
  for(i=0; i<(*ROW); i++) {
    design[i] = (double *) malloc((*COL)*sizeof(double));
    if(design[i] == NULL){
      printf("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
      exit(1);
    }
    for(j=0; j<(*COL); j++){
      loc = i*(*COL) + j;
      design[i][j] = Design[loc];
    }
  }
  
  N = 0;
  Scale1 = 1.0/(*ROW);
  for(i=0; i<((*ROW)-1); i++) {
    tempCPE = 0;
    tempCPEapp = 0;
    for(j=(i+1); j<(*ROW); j++) {
      if(xbeta[j] != xbeta[i]){
	N = N + 1;
	bxjxi = xbeta[j] - xbeta[i];
	bxixj = 0 - bxjxi;
	denomji = 2 + expm1(bxjxi);
	denomij = 2 + expm1(bxixj);
	cdfbxij = pnorm(bxixj/(*bandwidth), 0, 1, 1, 0);
	cdfbxji = pnorm(bxjxi/(*bandwidth), 0, 1, 1, 0);
	pdfbxij = dnorm(bxixj/(*bandwidth), 0, 1, 0);
	pdfbxji = dnorm(bxjxi/(*bandwidth), 0, 1, 0);
	Uji = cdfbxij/denomji;
	Uij = cdfbxji/denomij;
	tempCPE += 1.0*(bxjxi <= 0)/denomji + 1.0*(bxixj < 0)/denomij;  /* [1] Kn(beta.hat) */
	tempCPEapp += Uji + Uij;                  /* smooth approximation for Kn(beta.hat)*/
	/* The following is for the asymptotic variance, form [5] */
	F77_CALL(dcopy)(COL,design[j],&incx, xji,&incy); /* copy design[j] to xji */
	F77_CALL(dcopy)(COL,design[i],&incx, tempv,&incy); /* copy design[i] to tempv */
	vector_sub(*COL, xji, tempv);
	F77_CALL(dcopy)(COL, xji, &incx, xij, &incy); 
	vector_scale(*COL, xij,-1.0);
  
	F77_CALL(dcopy)(COL, xij, &incx, tempxij, &incy); /* copy xij to tempxij */
	vector_scale(*COL, tempxij, pdfbxij/(denomji*(*bandwidth)));
	vector_scale(*COL, xij, cdfbxij*(denomji - 1.0)/(denomji*denomji));
	F77_CALL(dcopy)(COL, xji, &incx, tempxji, &incy);
	vector_scale(*COL, tempxji, pdfbxji/(denomij*(*bandwidth)));
	vector_scale(*COL, xji, cdfbxji*(denomij - 1)/(denomij*denomij));

	vector_add(*COL, tempxij, xij); /* tempxij = tempxij + xij */
	vector_add(*COL, tempxij, tempxji);
	vector_add(*COL, tempxij, xji);
	vector_scale(*COL, tempxij, Scale1);
	vector_add(*COL, sumdervec, tempxij);
      }	
    }	
    CPE += Scale1 * tempCPE;
    CPEapp += Scale1 * tempCPEapp;
  }
  
  vector_scale(*COL, sumdervec, 1.0*(*ROW)/N); /* ROW = 1/scale1, attention (*ROW)/N == 0*/
  F77_CALL(dgemv)(trans, COL, COL, &ONE, varbeta, COL, sumdervec, &incx, &ZERO, tempv, &incy); 
  varterm2 = F77_CALL(ddot)(COL, sumdervec, &incx, tempv, &incy);

  CPE = (*ROW)*CPE/N;
  CPEapp = (*ROW)*CPEapp/N;

  /* now to calculate varterm1 */
  UUKij = 0;
  for(i=0; i<(*ROW); i++) {
    UUKsumj = 0.0;
    for(j=0; j<(*ROW); j++) {
      if(xbeta[j] != xbeta[i]){ 
	bxjxi = xbeta[j] - xbeta[i];
	bxixj = 0 - bxjxi;
	denomji = 2 + expm1(bxjxi);
	denomij = 2 + expm1(bxixj);
	cdfbxij = pnorm(bxixj/(*bandwidth), 0, 1, 1, 0);
	cdfbxji = pnorm(bxjxi/(*bandwidth), 0, 1, 1, 0);
	UUK = cdfbxij/denomji + cdfbxji/denomij - CPEapp;
	UUKsumj += UUK;
	UUKij += Scale1 * UUK*UUK;
      }
    }
    UUKSUM += Scale1 * UUKsumj*UUKsumj;
  }

  varterm1 = (UUKSUM - UUKij)*(*ROW)/N/N; /*note if divided by (N*N) will produce overflow problems */
  result[0] = CPE;
  result[1] = CPEapp;
  result[2] = sqrt(varterm1 + varterm2);
  printf("Result - %d, %d,%f, %f, %f,%f, %f\n",N,*ROW,UUKSUM,UUKij,varterm1,varterm2,result[2]);

  for(i = 0; i < (*ROW); i++){ free(design[i]); }
  free(design);
  free(xij);
  free(xji);
  free(tempv);
  free(tempxij);
  free(tempxji);
  free(sumdervec);
}

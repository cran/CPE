/*
Qianxing Mo, qianxing.mo@moffitt.org
Department of Biostatistics & Bioinformatics
H. Lee Moffitt Cancer Center & Research Institute
12902 Magnolia Drive, Tampa, FL 33612 

These C functions are used to compute the concordance probability for the Cox proportional hazard model
(see Gonen and Heller, Biometrika(2005), 92, 4, pp. 965-970)

The algorithm used to compute the probability is based on the R codes written by  Ennapadam. S. Venkatraman at MSKCC, with further optimization by Qianxing Mo.  

/* Last updated 6/14/2012; updated cpeNoTies function for the new algorithm to calculate the variance when covariates are tied*/

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
    error("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
    /*   error(1); */
  }

  design = (double **) malloc((*ROW)*sizeof(double *));
  if(design==NULL){
    error("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
  }
  
  for(i=0; i<(*ROW); i++) {
    design[i] = (double *) malloc((*COL)*sizeof(double));
    if(design[i] == NULL){
      error("Error: Fail to allocate memory space. Your computer may not have enough memory. \n");
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
  /*  Rprintf("Result --- varterm1 = %f, varterm2 = %f\n",varterm1, varterm2); */

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
  *result = CPE/N*(*ROW);
}


/* Function to compute cpe and its variance, not tied; last updated 7/8/2012  */
void cpeNoTies(const int *ROW, const int *COL, const double *bandwidth, const double *xbeta, const double * Design, const double *varbeta, double * result) {
  
  int i, j, k, loc;
  int incx, incy, TWO;
  double ONE,ZERO,dROW;
  char *trans = "N";
  double CPE=0,sumu2=0,varterm1=0,varterm2=0,rij1=0,rij2=0,sumrij1=0,sumrij2=0,tempCPE,kappa1,kappa2; 
  double bxjxi,bxixj,denomji,denomij,uij2,uji2,pdfbxij,pdfbxji,uji1,uij1,Scale1; 
  double *xji, *xij, *tempv, *tempxij, *tempxji, *sumdu1, *duji1, *duij1,*duji2, *duij2;
  double **design;
  double *dtran, *V1, *V2, *tempv2;

  incx = 1;
  incy = 1;
  TWO = 2;
  ONE = 1.0;
  ZERO = 0;
  dROW = (*ROW)*1.0;

  xji = (double *) malloc((*COL)*sizeof(double));
  xij = (double *) malloc((*COL)*sizeof(double));
  tempv = (double *) malloc((*COL)*sizeof(double));
  tempxij = (double *) malloc((*COL)*sizeof(double));
  tempxji = (double *) malloc((*COL)*sizeof(double));
  sumdu1 = (double *) calloc((*COL), sizeof(double));
  duji1 = (double *) malloc((*COL)*sizeof(double));
  duij1 = (double *) malloc((*COL)*sizeof(double));
  duji2 = (double *) malloc((*COL)*sizeof(double));
  duij2 = (double *) malloc((*COL)*sizeof(double));
  dtran = (double *) calloc(2,sizeof(double));
  V1 = (double *) calloc(4,sizeof(double));
  V2 = (double *) calloc(4,sizeof(double));
  tempv2 = (double *) calloc(2,sizeof(double));


  if(xji==NULL || xij==NULL || tempv==NULL || tempxij==NULL || tempxji==NULL || sumdu1==NULL || duji1==NULL || duij1==NULL || duji2==NULL || duij2==NULL || dtran==NULL || V1==NULL || V2==NULL || tempv2==NULL){
    error("Error: Fail to allocate memory space. \n");
  }
  design = (double **) malloc((*ROW)*sizeof(double *));
  if(design==NULL){
    error("Error: Fail to allocate memory space. \n");
  }
  for(i=0; i<(*ROW); i++) {
    design[i] = (double *) malloc((*COL)*sizeof(double));
    if(design[i] == NULL){
      error("Error: Fail to allocate memory space. \n");
    }
    for(j=0; j<(*COL); j++){
      loc = i*(*COL) + j;
      design[i][j] = Design[loc];
    }
  }
  

  kappa2 = 0;
  Scale1 = 1.0/(*ROW);
  kappa1 = 0;
  for(i=0; i<((*ROW)-1); i++) {
    tempCPE = 0;
    for(j=(i+1); j<(*ROW); j++) {
      if(xbeta[j] != xbeta[i]){
        kappa2 = kappa2 + 1;
        bxjxi = xbeta[j] - xbeta[i];   /* lcji */
        bxixj = 0 - bxjxi;
        denomji = 2 + expm1(bxjxi);   /* denominator of u_ji1*/
        denomij = 2 + expm1(bxixj);
        uij2 = pnorm(bxjxi/(*bandwidth), 0, 1, 1, 0);
        /*      uji2 = pnorm(bxixj/(*bandwidth), 0, 1, 1, 0); */
        uji2 = 1- uij2;              /* because uij2 + uji2 = 1 */
        pdfbxij = dnorm(bxjxi/(*bandwidth), 0, 1, 0);
        /*      pdfbxji = dnorm(bxixj/(*bandwidth), 0, 1, 0); */
        pdfbxji = pdfbxij;         /* because of symmetry */
        uji1 = uji2/denomji;        /* u_ji1 */
        uij1 = uij2/denomij;        /* u_ij1 */
        tempCPE += 1.0*(bxjxi < 0)/denomji + 1.0*(bxixj < 0)/denomij;  /* [1] Kn(beta.hat) */
        kappa1 += uji1 + uij1;                  /* smooth approximation for Kn(beta.hat)*/

        /* The following is for the asymptotic variance term 2 */
        F77_CALL(dcopy)(COL,design[j],&incx, xji,&incy);   /* copy design[j] to xji */
        F77_CALL(dcopy)(COL,design[i],&incx, tempv,&incy); /* copy design[i] to tempv */
        vector_sub(*COL, xji, tempv);                      /* X[j,] - X[i,]  */
        F77_CALL(dcopy)(COL, xji, &incx, xij, &incy);     /* xij = -xji */
        vector_scale(*COL, xij,-1.0);

        F77_CALL(dcopy)(COL, xji,  &incx, duji2, &incy); /* copy xji to duji2 */
        F77_CALL(dcopy)(COL, xij, &incx, duij2, &incy);
        vector_scale(*COL, duji2,-pdfbxji/(*bandwidth));   /* duji2 = -Xji/h*pdfbxji */
        vector_scale(*COL, duij2,-pdfbxij/(*bandwidth));   /* duij2 = -Xij/h*pdfbxij */

        F77_CALL(dcopy)(COL, duij2, &incx, duij1, &incy); /* copy duij2 to duij1 */
        vector_scale(*COL, duij1, denomij);               /* duij1 = denomij * duij2 */
        vector_scale(*COL, xij, uij2*(denomij-1.0));      /* xij = xij*uij2*(denomij-1) */
        vector_sub(*COL,duij1, xij);
        vector_scale(*COL, duij1, Scale1/denomij/denomij);

        F77_CALL(dcopy)(COL, duji2, &incx, duji1, &incy); /* copy duji2 to duji1 */
        vector_scale(*COL, duji1, denomji);
        vector_scale(*COL, xji, uji2*(denomji-1.0));
        vector_sub(*COL,duji1, xji);
        vector_scale(*COL, duji1, Scale1/denomji/denomji);

        vector_add(*COL, sumdu1, duij1);
        vector_add(*COL, sumdu1, duji1);
        sumu2 += uji2 + uij2;
      }
    }
    CPE += Scale1 * tempCPE;
  }

  vector_scale(*COL, sumdu1, 1.0/sumu2); /* now sumdu1 = sumdu1/sumu2 */
  vector_scale(*COL, sumdu1, 1.0/Scale1); /* now sumdu1 = sumdu1/sumu2 */
  /* Rprintf("sumdu1 = %f, sumdu2 = %f, sumdu3 = %f \n",sumdu1[0],sumdu1[1],sumdu1[2]);                                                                   
     Rprintf("sumu2 = %f, kappa1 = %f, kappa2 = %f \n",sumu2,kappa1,kappa2); */
  F77_CALL(dgemv)(trans, COL, COL, &dROW, varbeta, COL, sumdu1, &incx, &ZERO, tempv, &incy);
  varterm2 = F77_CALL(ddot)(COL, sumdu1, &incx, tempv, &incy);

  CPE = CPE/kappa2*(*ROW);                 /* scale1 = 1/ROW   */
  kappa1 = kappa1/(*ROW)/(*ROW-1)*2.0;
  kappa2 = kappa2/(*ROW)/(*ROW-1)*2.0;
  dtran[0] = 1.0/kappa2;
  dtran[1] = 0 - kappa1/kappa2/kappa2;

  /* now to calculate varterm1 */
  for(i=0; i<(*ROW); i++) {
    sumrij1 = 0;
    sumrij2 = 0;
    for(j=0; j<(*ROW); j++) {
      if(xbeta[j] != xbeta[i]){
	bxjxi = xbeta[j] - xbeta[i];
	bxixj = 0 - bxjxi;
	denomji = 2 + expm1(bxjxi);   /* denominator of u_ji1*/
	denomij = 2 + expm1(bxixj);
	uij2 = pnorm(bxjxi/(*bandwidth), 0, 1, 1, 0);
	uji2 = 1- uij2;              /* because uij2 + uji2 = 1 */
	uji1 = uji2/denomji;        /* u_ji1 */
	uij1 = uij2/denomij;        /* u_ij1 */
      
	rij1 = uij1+uji1 - kappa1;
	rij2 = uij2+uji2 - kappa2;
      }else{
	rij1 = 0 - kappa1;
	rij2 = 0 - kappa2;
      }
      sumrij1 += rij1;
      sumrij2 += rij2;
      /* V2 is for the 3rd loops when j = k */
      V2[0] = V2[0] + Scale1*rij1*rij1;  /* V[0][0] */
      V2[2] = V2[2] + Scale1*rij1*rij2;  /* V[0][1] */
      V2[3] = V2[3] + Scale1*rij2*rij2;  /* V[1][1] */
    }
    V1[0] = V1[0] + Scale1*sumrij1*sumrij1;  /* V[0][0] */
    V1[2] = V1[2] + Scale1*sumrij1*sumrij2;  /* V[0][1] */
    V1[3] = V1[3] + Scale1*sumrij2*sumrij2;  /* V[1][1] */
  }

  for(i=0; i<4; i++){
    V1[i] = V1[i] - V2[i];
  }
  V1[1] = V1[2];

  for(i=0; i<4; i++){
    V1[i] = V1[i]/(*ROW-1)/(*ROW-1)*4.0;
  }

  /* printf("Good before dgemv; %f, %f, %f \n",dtran[0],dtran[1],varterm2); */
  F77_CALL(dgemv)(trans, &TWO, &TWO, &ONE, V1, &TWO, dtran, &incx, &ZERO, tempv2, &incy); 
  varterm1 = F77_CALL(ddot)(&TWO, dtran, &incx, tempv2, &incy);

  result[0] = CPE;
  result[1] = 0;
  result[2] = sqrt((varterm1 + varterm2)/(*ROW));
  /* Rprintf("Result - %f,%f, %f\n",varterm1,varterm2,result[2]); */

  for(i = 0; i < (*ROW); i++){ free(design[i]); }
  free(design);
  free(xij);
  free(xji);
  free(tempv);
  free(tempxij);
  free(tempxji);
  free(sumdu1);
  free(duij1);
  free(duji1);
  free(duij2);
  free(duji2);
  free(dtran);
  free(V1);
  free(V2);
  free(tempv2);
}


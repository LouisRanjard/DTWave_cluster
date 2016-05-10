/*
* =============================================================
* DTWave.c
*
* Mex C function to calculate the dynamic time warp distance between
* two arrays. The arrays may be of different length. Also computes a 
* weighted average vector sequence.
*
*
* =============================================================
*/

#include <mex.h>
/* #include <matrix.h> */
#include <math.h>
/* #include <mcheck.h> for debugging */
/* #include <stdio.h> for debugging */

/* PROTOTYPES */
double warpav(double *traceb, double *d, int numrows, int numcols, double *indelc);
double warpav_ce(double *traceb, double *d, int numrows, int numcols, double *indelc);
double *averseq(int *duree, double *traceb, double *mat1, double *mat2, double weight, int numrows, int numcols, int numv);
int mincost(double cost[], int a, int b);
double louiround(double x);
void euclidist(double *d, double *mat1, double *mat2, double *indelc, int numrows, int numcols, int numv, double *cow);

double warpav_ce(double *traceb, double *d, int numrows, int numcols, double *indelc)
{
    int i,j,k = 0 ; /* i rows ; j cols ; k traceback */
    double cost[5] ;
    double *dsf,dist ;

    dsf = (double *)mxMalloc((numrows+1) * (numcols+1) * sizeof(double));

    /* first cell (1,1) */
    *dsf = *d * 2;
    *traceb = 0;

    /* first column */
   for (i=1;i<=numcols;i++) {
       *(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
       *(traceb+((numrows+1)*i)) = 1;
   }

    /* first row */
   for (j=1;j<=numrows;j++) {
       *(dsf+j) = *(dsf+(j-1)) + *indelc;
       *(traceb+j) = 3;
   }

    /* cell (2,2) */
   cost[1] = *(dsf+1) + *indelc;
   cost[2] = *dsf + *d;
   cost[3] = *(dsf+numrows+1) + *indelc;
   k = mincost(cost,1,3);
   *(traceb+numrows+2) = k;
   *(dsf+numrows+2) = cost[k];

    /* second column */
   j=1;
   for (i=2;i<=numcols;i++) {
       cost[0] = *(dsf+((numrows+1)*(i-2))+(j-1)) + *(d+((numrows)*(i-2)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-1)))*0.5;
       cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
       cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
       cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
       k = mincost(cost,0,3);
       *(traceb+((numrows+1)*i)+1) = k;
       *(dsf+((numrows+1)*i)+1) = cost[k];
   }

    /* second row */
   i=1;
   for (j=2;j<=numrows;j++) {
       cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
       cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
       cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
       cost[4] = *(dsf+((numrows+1)*(i-1))+(j-2)) + *(d+((numrows)*(i-1)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-2)))*0.5;
       k = mincost(cost,1,4);
       *(traceb+((numrows+1))+j) = k;
       *(dsf+(numrows+1)+j) = cost[k];
   }

    /* rest of the matrix */
   for (i=2;i<=numcols;i++) {
       for (j=2;j<=numrows;j++) {
           cost[0] = *(dsf+((numrows+1)*(i-2))+(j-1)) + *(d+((numrows)*(i-2)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-1)))*0.5;
           cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
           cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
           cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
           cost[4] = *(dsf+((numrows+1)*(i-1))+(j-2)) + *(d+((numrows)*(i-1)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-2)))*0.5;
           k = mincost(cost,0,4);
           *(traceb+((numrows+1)*i)+j) = k;
           *(dsf+((numrows+1)*i)+j) = cost[k];
       }
   }
/* print matrices
for (i=0;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       printf(" %f", *(traceb+((numrows+1)*(i))+j));
   }
   printf("\n");
}
printf("\n");
for (i=0;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       printf(" %f", *(dsf+((numrows+1)*(i))+j));
   }
   printf("\n");
} */
   /* compute distance */
   /* last cell */
   dist = (*(dsf+((numrows+1)*(numcols+1))-1)) / (double)(numrows+numcols);
   /* end free space
   double mini=*(dsf+((numrows+1)*(numcols+1))-1), minj=*(dsf+((numrows+1)*(numcols+1))-1) ;
   for (i=0;i<=numcols;i++) {if ((*(dsf+((numrows+1)*(i))+numrows))<mini) {mini=*(dsf+((numrows+1)*(i))+numrows);} }
   for (j=0;j<=numrows;j++) {if ((*(dsf+((numrows+1)*(numcols))+j))<minj) {minj=*(dsf+((numrows+1)*(numcols))+j);} }
   if (minj<mini) {dist = minj;} else {dist = mini;} */

   mxFree(dsf);
   return dist ;
}

double warpav(double *traceb, double *d, int numrows, int numcols, double *indelc)
{
    int i,j,k = 0 ; /* i rows ; j cols ; k traceback */
    double cost[4] ; /* need to keep 4 cells for compatibility with warpav_ce */
    double *dsf,dist ;

    dsf = (double *)mxMalloc((numrows+1) * (numcols+1) * sizeof(double));

    /* first cell (1,1) */
    *dsf = *d * 2;
    *traceb = 0;

    /* first column */
   for (i=1;i<=numcols;i++) {
       *(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
       *(traceb+((numrows+1)*i)) = 1;
   }

    /* first row */
   for (j=1;j<=numrows;j++) {
       *(dsf+j) = *(dsf+(j-1)) + *indelc;
       *(traceb+j) = 3;
   }

    /* rest of the matrix */
   for (i=1;i<=numcols;i++) {
       for (j=1;j<=numrows;j++) {
           cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
           cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
           cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
           k = mincost(cost,1,3);
           *(traceb+((numrows+1)*i)+j) = k;
           *(dsf+((numrows+1)*i)+j) = cost[k];
       }
   }
/* print matrices
for (i=0;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       printf(" %f", *(traceb+((numrows+1)*(i))+j));
   }
   printf("\n");
}
printf("\n");
for (i=0;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       printf(" %f", *(dsf+((numrows+1)*(i))+j));
   }
   printf("\n");
} */
   /* compute distance */
   /* last cell */
   dist = (*(dsf+((numrows+1)*(numcols+1))-1)) / (double)(numrows+numcols);
   /* end free space
   double mini=*(dsf+((numrows+1)*(numcols+1))-1), minj=*(dsf+((numrows+1)*(numcols+1))-1) ;
   for (i=0;i<=numcols;i++) {if ((*(dsf+((numrows+1)*(i))+numrows))<mini) {mini=*(dsf+((numrows+1)*(i))+numrows);} }
   for (j=0;j<=numrows;j++) {if ((*(dsf+((numrows+1)*(numcols))+j))<minj) {minj=*(dsf+((numrows+1)*(numcols))+j);} }
   if (minj<mini) {dist = minj;} else {dist = mini;} */

   mxFree(dsf);
   return dist ;
}

/* compute the average sequence given the alignment */
double *averseq(int *duree, double *traceb, double *mat1, double *mat2, double weight, int numrows, int numcols, int numv)
{
    int i,j,n=0;
    int l=0;
    int maxt=0; /* length of the matrix */
    int ta,tb;
    double *mat4,*t,*mata;

    mat4 = (double *)mxMalloc((numcols+numrows)*numv*sizeof(double)) ;
    t = (double *)mxMalloc((numcols+numrows)*sizeof(double)) ;

    i = numcols;
    j = numrows;

    /* trace back the alignment and implement the average vector sequence */
    while ( i>0 || j>0 ) {
       *(t+l) = weight*j + (1-weight)*i;           /*printf(" %f->%f\n",*(t+l),louiround(*(t+l)) );*/
       if ( (int)(louiround(*(t+l)))>maxt ) { maxt=(int)(louiround(*(t+l))); } 
       switch ( (int) *(traceb+((numrows+1)*(i))+(j)) ) {
           case 0:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = weight*(*(mat1+n+(numv*(j-1)))) + 0.5*(1-weight)*(*(mat2+n+(numv*(i-1)))) + 0.5*(1-weight)*(*(mat2+n+(numv*(i-2)))) ; }
               i-=2;
               j--;
               break;
           case 1:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = (*(mat2+n+(numv*(i-1)))); }
               i--;
               break;
           case 2:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = weight*(*(mat1+n+(numv*(j-1)))) + (1-weight)*(*(mat2+n+(numv*(i-1)))); }
               i--;
               j--;
               break;
           case 3:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = (*(mat1+n+(numv*(j-1)))); }
               j--;
               break;
           case 4:
               for (n=0;n<numv;n++) {
                    *(mat4+n+(numv*l)) = 0.5*weight*(*(mat1+n+(numv*(j-1)))) + 0.5*weight*(*(mat1+n+(numv*(j-2)))) + (1-weight)*(*(mat2+n+(numv*(i-1)))); }
               i--;
               j-=2;
               break;
           default:
               printf("error: traceback no match\n");
               break;
       }
       l++;
    }
    mata = (double *)mxMalloc(maxt*numv*sizeof(double)) ;

    /* compute average with correct time direction and number of time points */
    for ( ta=1 ; ta<=maxt ; ta++ ) {
       for ( tb=l-1 ; tb>=0 ; tb-- ) {
           if ( *(t+tb)==(float)ta ) {
               for ( n=0 ; n<numv ; n++ ) {
                    *(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
               break;
           } else if ( *(t+tb)>(float)ta ) {
               if ( tb==l-1 ) {
                    for ( n=0 ; n<numv ; n++ ) {
                            *(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
               } else {
                    for ( n=0 ; n<numv ; n++ ) {
                            *(mata+n+(numv*(ta-1))) = 0.5*(*(mat4+n+(numv*(tb+1)))) + 0.5*(*(mat4+n+(numv*(tb)))) ; }
               }
               break;
           } else if ( tb==0 ) {
               for ( n=0 ; n<numv ; n++ ) { 
                    *(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
           }
       }
    }

    *duree = maxt ;
    mxFree(mat4);
    mxFree(t);
    return mata ;
}

/* find the minimum value in an array between indexes a and b
int mincost(double cost[], int a, int b)
{
    int j=a;
    while (a<b+1) {
        if (cost[a]<cost[j]) {j=a;}
        a++;
    }
    return j;
}*/
/* find the minimum value in an array between indexes a and b
and returns its index
if different values are equals, randomly chose one using time, quite slower */
int mincost(double cost[], int a, int b)
{
   int j=a;
   while (a<=b) {
       if (cost[a]<cost[j]) {j=a;}
       else if (cost[a]==cost[j]) {
            if (time(NULL)%2==0) {j=a;} }
       a++;
   }
   return j;
}

/* no round() function in math.h so here it is, SOMETIMES THIS FUNCTION IS BUILT IN */
double louiround(double x)
{
    double intpart;
    if(modf(x,&intpart)>=0.5){
      return ceil(x);}
    else{
      return floor(x);}
}

/* compute Euclidean distance between vectors of 2 matrices 
   the insertion deletion cost is set as (half) the average distance between vectors */
void euclidist(double *d, double *mat1, double *mat2, double *indelc, int numrows, int numcols, int numv, double *cow)
{
    int i,j,n ;

    for (i=0;i<numcols;i++) {
       for (j=0;j<numrows;j++) {
           *(d+((numrows)*i)+j)=0;
           for ( n=0 ; n<numv ; n++ ) {
               *(d+((numrows)*i)+j) = *(d+((numrows)*i)+j) + *(cow+n)*pow( ((*(mat2+n+(numv*(i))))-(*(mat1+n+(numv*(j))))), 2 ) ; }
           *(d+((numrows)*i)+j) = sqrt( *(d+((numrows)*i)+j) ) ;
           *indelc += *(d+(numrows*i)+j) ;
       }
    }
    *indelc = *indelc/(numrows*numcols) ;
    return ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *indelc, *d, *mat1, *mat2, *traceb, *mataverage, *distance, *cow ;
    int numrows, numcols, numv, *duree, ce ;
    /* mwSize numrows, numcols, numv, *duree ; */
    double weight ;

    mat1 = mxGetPr(prhs[0]);
    mat2 = mxGetPr(prhs[1]);
    cow = mxGetPr(prhs[2]);
    weight = mxGetScalar(prhs[3]);
    ce = mxGetScalar(prhs[4]);

    numrows = mxGetN(prhs[0]); /* Number of columns in mxArray */
    numcols = mxGetN(prhs[1]);
    numv = mxGetM(prhs[0]); /* Number of rows in mxArray */

    /* memory allocation */
    d = (double *)mxMalloc((numrows) * (numcols) * sizeof(double));
    traceb = (double *)mxMalloc((numrows+1) * (numcols+1) * sizeof(double));
    indelc = (double *)mxMalloc(sizeof(double));
    duree = (int *)mxMalloc(sizeof(int));

    /* Create variable for return values and get pointers */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL); /* DISTANCE */
    distance = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL); /* AVERAGE VECTOR SEQUENCE */

    /* compute pairwise vector distance */
    *indelc = 0 ;
    euclidist(d,mat1,mat2,indelc,numrows,numcols,numv,cow);
    
    /* Call the warpav subroutine */
    if (ce){
        *distance = warpav_ce(traceb,d,numrows,numcols,indelc) ;
    }else{
        *distance = warpav(traceb,d,numrows,numcols,indelc) ;
    }
    
    /* Call the traceback subroutine to get the average sequence */
    if ( weight!=0 ) {
        mataverage = averseq(duree,traceb,mat1,mat2,weight,numrows,numcols,numv);
        /* set the pointer for return matrix */
        mxSetPr( plhs[1], mataverage );
        mxSetM( plhs[1], numv );
        mxSetN( plhs[1], *duree );
        mexMakeMemoryPersistent(mataverage); /* memory deallocation bug with octave */
    }

    mxFree(duree);
    mxFree(d);
    mxFree(traceb);
    mxFree(indelc);
}

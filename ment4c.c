#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h> 

#include <stddef.h> 
#include <string.h> 

#define PI 3.141592653589793238462643383279 

/*
%--------------------------------------------------------
%
%    To reconstruct image from sinogram
%    using the maximum entropy method described in
%
%    C. T. Mottershead, "Maximum entropyy beam diagnostic tomography,"
%    IEEE Transactions on Nuclear Science, Vol. NS-32, No. 5, October 1985.
%
%    Kai Hock 28 December 2009
%
%
%    Input:  ProjectFile - name of file for projection data
%    Output: ReconFile   - name of file for reconstructed image
%    tol:    tolerance for Gauss-Seidel calculations (try 0.01 or 0.001)
%
%--------------------------------------------------------------------

function ment(ProjectFile, ReconFile, tol)

load(ProjectFile);

%--------------------------------------------------------------------
%
% npos      - number of positions
% positions - matrix [nproj x npos] for positions of the projections
% nproj     - number of projections;
% angles    - column vector [nproj x 1] of intervals, i.e. angles
% centre    - column vector [nproj x 1], 
%             the position of the centre of rotation for each projection
% sinogram  - matrix of projection data [nproj x npos], 
%             where npos = number of positions across the projections
%
%--------------------------------------------------------------------
*/

#define NR_END 1
#define FREE_ARG char*


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
   


double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
   
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}




double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
   
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
   
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
   
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   
	/* return pointer to array of pointers to rows */
	return m;
}





void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
   

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
   
int npos;          /*- number of positions                                     */
int nproj;         /*- number of projections;                                  */
double *weights;
double *angles;    /*- column vector [nproj x 1] of intervals, i.e. angles     */
double *centre;   /*- column vector [nproj x 1],                               */ 
                 /* the position of the centre of rotation for each projection */
double **positions;/*- matrix [nproj x npos] for positions of the projections  */
double **sinogram;/*- matrix of projection data [nproj x npos],                */
                 /* where npos = number of positions across the projections    */

void fread_sino(char *ProjectFile)
/* projection data loaded into global arrays declared above: 
 npos      - number of positions
 positions - matrix [nproj x npos] for positions of the projections
 nproj     - number of projections;
 angles    - column vector [nproj x 1] of intervals, i.e. angles
 centre    - column vector [nproj x 1], 
             the position of the centre of rotation for each projection
 sinogram  - matrix of projection data [nproj x npos], 
             where npos = number of positions across the projections
*/
{
  int i;
  FILE *fid;

  fid = fopen(ProjectFile, "rb");
  if (fid == NULL) nrerror("Error reading ProjectFile");

  /* read array sizes */
  fread(&nproj, sizeof(int), 1, fid);
  fread(&npos, sizeof(int), 1, fid);

  /* define arrays */
  angles = dvector(0, nproj-1);
  weights = dvector(0, nproj-1);
  centre = dvector(0, nproj-1);
  
  positions = dmatrix(0, nproj-1, 0, npos-1);
  sinogram = dmatrix(0, nproj-1, 0, npos-1);
  
  /* read arrays */  
  fread(angles, nproj*sizeof(double), 1, fid);
  fread(weights, nproj*sizeof(double), 1, fid);
  fread(centre, nproj*sizeof(double), 1, fid);

  for (i=0; i<nproj; i++)
    fread(positions[i], npos*sizeof(double), 1, fid);
  
  for (i=0; i<nproj; i++)
    fread(sinogram[i], npos*sizeof(double), 1, fid);
  
  fclose(fid);
}

double max2(double **A, int m, int n)
{
  int i, j;
  double a;
 
  a = 0;
  for (i=0; i<m; i++)
    for (j=0; j<n; j++) {

      if (a < A[i][j])  a = A[i][j];
  }

  return a;
}

void sino_txt(char *filename)
{
  int i, j;
  double a;
  FILE *u;
  u = fopen(filename,"w");

/*  a = max2(sinogram, nproj, npos);*/
  fprintf(u,"# nproj %d, npos %d   \n", nproj, npos);

  for (i=0; i<nproj; i++) {
  
    fprintf(u,"\n\n# %6.3f deg  \n", angles[i]/PI*180);

    for (j=0; j<npos; j++)
      fprintf(u,"%9.6f   %6.3f\n", positions[i][j], sinogram[i][j]);

    fprintf(u,"\n");
  }
  
  fclose(u);

}



void h_txt(double **h)
{
  int i, j;
  double a;
  FILE *u;
  u = fopen("h.txt","w");

/*  a = max2(sinogram, nproj, npos);*/
  fprintf(u,"# nproj %d, npos %d   \n", nproj, npos);

  for (i=0; i<nproj; i++) {
  
    fprintf(u,"\n\n# %6.3f deg  \n", angles[i]/PI*180);

    for (j=0; j<npos; j++)
      fprintf(u,"%9.6f   %6.3f\n", positions[i][j], h[i][j]);

    fprintf(u,"\n");
  }
  
  fclose(u);

}



void image_txt(char *file, double *x, double *y, double **A, int m, int n)
{
  int i, j;
  double a;
  FILE *u;
  u = fopen(file,"w");

  a = max2(A, m, n);

  for (i=0; i<m; i++) {
  
    for (j=0; j<n; j++)
      fprintf(u,"%9.6f   %9.6f   %9.6f\n", x[i], y[n-j-1], A[i][n-j-1]/a);

    fprintf(u,"\n");
  }
  
  fclose(u);

}





double interp1(double *x, double *y, double xi, int n)
/* 
   linear interpolation over x, y 
   if outside range of x return 0
*/
{
  int i;
  double yi;

  if ((xi < x[0]) || (xi >= x[n-1]))  return 0;

  for (i=0; i<n-1; i++) if ((xi >= x[i]) && (xi < x[i+1])) break;
  
  yi = y[i] + (xi - x[i])/(x[i+1] - x[i])*(y[i+1] - y[i]);

  return yi;
}           


double check_diff(double **h, double **h0, int nproj, int npos)
{
  int i, j;
  double d, a;

  a = 0;
  d = 0;

  for (i=0; i<nproj; i++)
    for (j=0; j<npos; j++) {

      if (d < fabs(h[i][j] -  h0[i][j]))  d = fabs(h[i][j] -  h0[i][j]);
      if (a < fabs(h0[i][j]))           a = fabs(h0[i][j]);
    }

  if (a == 0) nrerror("double check_diff() - Error: a = 0.");
  d = d/a;

  return d;
}



void convergence(int ii, double *d, int ni)
{
  int i;
  printf("-------------------------------------------\n");
  
  for (i=0; i<ii; i++)
    printf("Iteration %d of %d, error = %f\n", i, ni, d[i]);
}

void matrix_eq(double **A, double **B, int m, int n)
{
  int i, j;

  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      A[i][j] = B[i][j];
}


double trapezium_rule(double *y, int n, double h)
{
  int i;
  double a;
  a = 0;
  for (i=0; i<n-1; i++)
    a = a + (y[i] + y[i+1])*h/2;

  return a;
}


double integrate(double x1, double x2, double *y, double *x, int n)
/*
  integrate y wrt x from x1 to x2, using the trapezium rule

  x1   lower limit
  x2   upper limit
  y    y array
  x    x array
  n    size of both arrays
*/
{
  int i, i1, i2;
  double y1, y2, sum;

/* locate x1 and x2 in x array */
  i1 = -1;
  i2 = -1;
  for (i=0; i<n-1; i++) {
    if ((x1>=x[i]) && (x1<=x[i+1])) {
      i1 = i;
      break;
    }
  }
  for (i=0; i<n-1; i++) {
    if ((x2>=x[i]) && (x2<=x[i+1])) {
      i2 = i;
      break;
    }
  }
  if (x2 < x1) nrerror("upper limit < lower limit");
  if (i1 == -1) { printf("x1 = %lf\n", x1); nrerror("x1 not in range of x array"); }
  if (i2 == -1) nrerror("x2 not in range of x array");

/* integrate from x1 to x2 */

  /* at lower limit */
  y1 = interp1(x, y, x1, n);
  sum = 0.5*(x[i1+1] - x1)*(y1 + y[i1+1]);

  /* in between */
  for (i=i1; i<i2-1; i++)
    sum = sum + 0.5*(x[i+1] - x[i])*(y[i] + y[i+1]);

  /* at upper limit */
  y2 = interp1(x, y, x2, n);
  sum = sum + 0.5*(x2 - x[i2])*(y[i2] + y2);

  return sum;
}



double integrate_1(double x1, double x2, double *y, double *x, int n)
/*
  integrate y wrt x from x1 to x2, using the trapezium rule

  x1   lower limit
  x2   upper limit
  y    y array
  x    x array
  n    size of both arrays
  
  If x is partly outside of x1-x2, just integrate the part that is inside.
*/
{
  int i, i1, i2;
  double y1, y2, sum;

/* locate x1 and x2 in x array */
  i1 = -1;
  i2 = -1;
  if (x1<x[0]) i1=0;
  for (i=0; i<n-1; i++) {
    if ((x1>=x[i]) && (x1<=x[i+1])) {
      i1 = i;
      break;
    }
  }
  if (x2>x[n-1]) i2=n-1;
  for (i=0; i<n-1; i++) {
    if ((x2>=x[i]) && (x2<=x[i+1])) {
      i2 = i;
      break;
    }
  }
  if (x2 < x1) nrerror("upper limit < lower limit");
  if (i1 == -1) { printf("x1 = %lf\n", x1); nrerror("x1 not in range of x array"); }
  if (i2 == -1) nrerror("x2 not in range of x array");

/* integrate from x1 to x2 */

  /* at lower limit */
  y1 = interp1(x, y, x1, n);
  sum = 0.5*(x[i1+1] - x1)*(y1 + y[i1+1]);

  /* in between */
  for (i=i1; i<i2-1; i++)
    sum = sum + 0.5*(x[i+1] - x[i])*(y[i] + y[i+1]);

  /* at upper limit */
  y2 = interp1(x, y, x2, n);
  sum = sum + 0.5*(x2 - x[i2])*(y[i2] + y2);

  return sum;
}



double TrapeziumRule(double *y, double *x, int n)
/*
  integrate y wrt x using the trapezium rule

  y    y array
  x    x array
  n    size of both arrays (must be > 1)
  
*/
{
  int i;
  double h, sum;

  if (n <= 1) nrerror("TrapeziumRule(): n is <= 1.  n must be > 1.");
  
  h = (x[n-1] - x[0]) / (n-1); 
  sum = 0;

  for (i=1; i<n-1; i++)
    sum = sum + y[i];

  sum = y[0] + 2*sum + y[n-1];
	
  return h*sum/2;
}


double maxrange()
{
  int i;
  double xmax, a, b;
  
  xmax = 0;
  
  for (i=0; i<nproj; i++) {

    a = fabs(positions[i][0]   - centre[i]);
    b = fabs(positions[i][npos] - centre[i]);
  
    if (xmax < a) xmax = a;
    if (xmax < b) xmax = b;
  }
  
  return xmax;
}


double maxrange_1()
{
  int i;
  double xmax, a, b;
  
  xmax = 0;
  
  for (i=0; i<nproj; i++) {

    a = fabs(positions[i][0]   - centre[i]);
    b = fabs(positions[i][npos] - centre[i]);
  
    if (xmax < a) xmax = a;
    if (xmax < b) xmax = b;
  }
  
  xmax = 2*xmax;
  
  return xmax;
}

void ComputeH(double **h, double **rho, double R, double dt, int nr, int N1)
/*
 As defined in Mottershead (1985),  h_n(s) = exp(lambda_n(s) - 1/N),
 where lambda_n(s) is the Lagrange multiplier for 
 the nth angle and at position s, and N is the total number of angles.
*/
{
  int it, k, ipos, iproj;
  double s, t, t1, sk, hk, integral;
  double *xn, *yn, *tn, *Pk;	
  xn = dvector(0, nr*N1-1);   /* for integration along each ray */
  yn = dvector(0, nr*N1-1);
  tn = dvector(0, nr*N1-1);
  Pk = dvector(0, nr*N1-1);   /* product of h_k(s), k != n */
  
  printf("     ComputeH():\n");
  
    for (iproj=0; iproj<nproj; iproj++)  {
    /*------------projection loop start ------------*/
      printf("               proj %d\n", iproj);

      for (ipos=0; ipos<npos; ipos++)  {
      /*------------position loop start ------------*/
        /* to find the integral */       
        
        /* 1. determine xn, yn along t (ray) */
        s = rho[iproj][ipos];
 	    t1 = sqrt(R*R - s*s);     /* integration limit */

        for (it=0; it<nr*N1; it++) { 

          t = -R + (double)it * dt;                 /* position along ray */  
          tn[it] = t;
          xn[it] = s * cos(angles[iproj]) - t * sin(angles[iproj]);
          yn[it] = s * sin(angles[iproj]) + t * cos(angles[iproj]);
        }

        /* 2. for each xn, yn                            */
        /*        - for each angle k (except n), find s_k  */
        /*        - then find h_k(s_k)                    */

        for (it=0; it<nr*N1; it++) {    /* for each position along ray */
        /*--------------ray loop start --------------------*/

		  if (fabs(tn[it]) >= t1) Pk[it] = 0; 

		  if (fabs(tn[it]) < t1) {
		  
          Pk[it] = 1;           /* initialise product of h_k(s_k) */
          for (k=0; k<nproj; k++) {  /* for each angle */
          /*--------------product loop start--------------------*/

            if (k != iproj) {        /*  except iproj (n) */

              sk = xn[it] * cos(angles[k]) + yn[it] * sin(angles[k]); /* position on ray k */
              hk = interp1(rho[k], h[k], sk, npos); /* find projection at corresponding s */

              /* 3. multiply all h_k(s_k) (except for k=n) and store in array Pk  */
              Pk[it] = Pk[it] * hk;
            }
          /*--------------product loop end--------------------*/
          }
		  }
        /*--------------ray loop end--------------------*/
        }

       /*  4.  integrate Pk using the trapezium rule, within circle of radius R */
       integral = TrapeziumRule(Pk, tn, nr*N1);

       /* Gauss Seidel */
       if (integral == 0) h[iproj][ipos] = 0;
       else h[iproj][ipos] = sinogram[iproj][ipos] / integral;  

      /*------------position loop end ------------*/
      }
    /*------------projection loop end ------------*/
    }

  free_dvector(xn, 0, nr*N1-1);  
  free_dvector(yn, 0, nr*N1-1);
  free_dvector(tn, 0, nr*N1-1);
  free_dvector(Pk, 0, nr*N1-1);  
}

	
void ComputeProjections(double **p, double **h, double **rho, 
                        double R, double dt, int nr, int N1)
{
  int it, k, ipos, iproj;
  double s, t, t1, sk, hk, integral;
  double *xn, *yn, *tn, *Pk;	
  xn = dvector(0, nr*N1-1);   /* for integration along each ray */
  yn = dvector(0, nr*N1-1);
  tn = dvector(0, nr*N1-1);
  Pk = dvector(0, nr*N1-1);   /* product of h_k(s), k != n */
  
  printf("     ComputeProjections():\n");
  
    for (iproj=0; iproj<nproj; iproj++)  {
    /*------------projection loop start ------------*/
      printf("                        proj %d\n", iproj);

      for (ipos=0; ipos<npos; ipos++)  {
      /*------------position loop start ------------*/
/*        printf("iter %d, proj %d, pos %d\n", ii, iproj, ipos);*/

        /* to find the integral */       
        
        /* 1. determine xn, yn along t (ray) */
        s = rho[iproj][ipos];
        t1 = sqrt(R*R - s*s);     /* integration limit */

        for (it=0; it<nr*N1; it++) { 

          t = -R + (double)it * dt;                 /* position along ray */  
          t = t;
          tn[it] = t;
          xn[it] = s * cos(angles[iproj]) - t * sin(angles[iproj]);
          yn[it] = s * sin(angles[iproj]) + t * cos(angles[iproj]);

        }

        /* 2. for each xn, yn                            */
        /*        - for each angle k, find s_k  */
        /*        - then find h_k(s_k)                    */

        for (it=0; it<nr*N1; it++) {    /* for each position along ray */
        /*--------------ray loop start --------------------*/

		  if (fabs(tn[it]) >= t1) Pk[it] = 0; 

		  if (fabs(tn[it]) < t1) {
		  
          Pk[it] = 1;           /* initialise product of h_k(s_k) */
          for (k=0; k<nproj; k++) {  /* for each angle */
          /*--------------product loop start--------------------*/

              sk = xn[it] * cos(angles[k]) + yn[it] * sin(angles[k]); /* position on ray k */
              hk = interp1(rho[k], h[k], sk, npos); /* find projection at corresponding s */

              /* 3. multiply all h_k(s_k) and store in array Pk  */
              Pk[it] = Pk[it] * hk;

          }
          /*--------------product loop end--------------------*/
          }
        /*--------------ray loop end--------------------*/
        }

       /*  4.  integrate Pk using the trapezium rule, within circle of radius R */
/*       integral = integrate_1(-t1, t1, Pk, tn, nr*N1);*/
       integral = TrapeziumRule(Pk, tn, nr*N1);

       p[iproj][ipos] = integral;

      /*------------position loop end ------------*/
      }
    /*------------projection loop end ------------*/
    }
	
  free_dvector(xn, 0, nr*N1-1);  
  free_dvector(yn, 0, nr*N1-1);
  free_dvector(tn, 0, nr*N1-1);
  free_dvector(Pk, 0, nr*N1-1);  
}

void ReconstructImage(double **fdata, double *x1, double *y1, 
                      double **h, double **rho, double a1, double b1, int N1)
/* -- reconstruct image ----------    */
/*  f(x, y) = product of all h_n(s_n(x, y)) */
{
  int i, j, n;
  double x2, y2, s, hn;
  
  for (i=0; i<N1; i++) 
    for (j=0; j<N1; j++) {
 
      fdata[i][j] = 1;

      for (n=0; n<nproj; n++) {

	  /* transform from normalised to real coordinates 
	    x2 =  x1[i]*sqrt(b1);
		y2 = -a1*x1[i]/sqrt(b1) + y1[j]/sqrt(b1);*/

	  /* transform from real to normalised coordinates */
	    x2 =  x1[i]/sqrt(b1);
		y2 =  a1*x1[i]/sqrt(b1) + y1[j]*sqrt(b1);
		
/*        s = x1[i] * cos(angles[n]) + y1[j] * sin(angles[n]);*/
        s = x2 * cos(angles[n]) + y2 * sin(angles[n]);
        hn = interp1(rho[n], h[n], s, npos);

        /* 3. multiply all h_n(s_n) and store in array f */
        fdata[i][j] = fdata[i][j] * hn;
      }
    }
}
	
/*----------------------------------------------*/
/*
  MENT code.
  Transform from normalised to real space during reconstruction.
*/
int main(int argc, char *argv[])
{
  char time1[50];
  time_t stime; 
  stime = time(NULL); 
  sprintf(time1, "%s", asctime(localtime(&stime))); 

/* . . . . . . . . . . . . . . . . . . . . . .  */
  int ni, nr, N1, ii, n, i, it, k, j, i1, j1, iproj, ipos;
  double tol, xmax, dx, dt, R, a1, b1, a, threshold;
  double s, t, t1, hk, integral, hn, sk, x, x2, y2;
  double *d, *x1, *y1, **p, **p0, **h, **h0, **rho, **fdata;
  char ProjectFile[100], ReconFile[100];

/* Twiss parameters for transforming to normalised phase space */
/* for no transformation */	
  a1 = 0;    
  b1 = 1;        

/* N1 x N1 is no. of pixels of image */
/*  N1 = 200;        
  xmax = 0.02; */

  if (argc != 3) nrerror("use: ment4c <xmax>  <N1>");  
  N1 = atoi(argv[2]);
  xmax = atof(argv[1]);
  printf("N1 = %d, xmax = %f\n", N1, xmax);
  
  tol = .01;       /* tolerance of Gauss-Seidel iterations - try 0.01   */
  ni = 30;         /* maximum no. of Gauss-Seidel iterations - try 10   */
  nr = 2;          /* no. of divisions of dx - for ray interval         */
  threshold = 0.01;/*set sinogram to 0 at this fraction of maximum value*/

  sprintf(ProjectFile, "sinogram.bin");
  sprintf(ReconFile, "recon_MENT.txt"); 
  
  fread_sino(ProjectFile);  

/*  xmax = maxrange_1();*/
  R  = xmax / 2;          /* radius of image circle */
  dx = xmax/(double)N1;   /* pixel size of reconstructed image */
  dt = dx/((double)nr);   /* interval along each ray, for integration purpose */

/*--------------------------------*/
  x1 = dvector(0, N1-1);      /* for reconstruction grid */
  y1 = dvector(0, N1-1);
  d = dvector(0, ni-1);        /* measure of convergence error */

  p = dmatrix(0, nproj-1, 0, npos-1);
  p0 = dmatrix(0, nproj-1, 0, npos-1);
  h = dmatrix(0, nproj-1, 0, npos-1);      /*  h_n(s), for all n */
  h0 = dmatrix(0, nproj-1, 0, npos-1);     /*  h_n(s) of previous iteration */     
  rho = dmatrix(0, nproj-1, 0, npos-1);    /* positions[][] - centre[] */ 
  fdata = dmatrix(0, N1-1, 0, N1-1);       /* reconstructed image */ 

  a = max2(sinogram, nproj, npos);
  
  for (i=0; i<nproj; i++)
    for (j=0; j<npos; j++) { 
      p[i][j] = 1;                             
      p0[i][j] = 1;

      h[i][j] = 1;                              /* initialising  h_n(s) */
      h0[i][j] = 1;
      rho[i][j] = positions[i][j] - centre[i];  /* centering each projection */ 
      
      /* threshold to avoid divergence of h at edge of image */
      if (sinogram[i][j] < threshold*a) sinogram[i][j] = 0;	  
    }

  sino_txt("sino_ment.txt");  /* save sinogram as text for checking with gnuplot */
  
/*---- Gauss-Seidel ----------------------------*/
  for (ii=0; ii<ni; ii++)  {     
    printf("iter %d, \n", ii);

	ComputeH(h, rho, R, dt, nr, N1);
	ComputeProjections(p, h, rho, R, dt, nr, N1);
	
    /* check convergence */
    d[ii] = check_diff(p, p0, nproj, npos);
    printf("                               error %lf\n", d[ii]);

    matrix_eq(p0, p, nproj, npos);	
    if (d[ii] < tol) break;
  }

  convergence(ii, d, ni);
  h_txt(h);

  for (i=0; i<N1; i++) {

    x1[i] = (-xmax/2 + (double)i * dx) * 1;
    y1[i] = x1[i] ;
  }

  ReconstructImage(fdata, x1, y1, h, rho, a1, b1, N1);
  image_txt(ReconFile, x1, y1, fdata, N1, N1);

  printf("xmax = %f m\n", xmax);

/* . . . . . . . . . . . . . . . . . . . . . .  */
  printf(" ----------------------------------- \n");
  stime = time(NULL); 
  printf(" Start %s End   %s \n", time1, asctime(localtime(&stime)));

  return 0;
}

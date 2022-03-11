#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>

int LMAXX = 10, NMAXX = 10;
int MAX_MESH_POINTS = 500;
int nmax;
	
double V(double r)
{
	double v;
	v = 0.5*r*r;
	return v;
}	

void Hsolve(int n, int L, int N, double hbar, double rmax, double E[N], double u[N][n])
{
	int i,j;
	double d, e1, e2;
	double h, r, vv;
	
	gsl_vector *eval = gsl_vector_alloc (n);
	gsl_matrix *evec = gsl_matrix_alloc (n, n);
	gsl_matrix *A = gsl_matrix_alloc (n, n);
	
	h = rmax/n;
	
	gsl_matrix_set_zero(A);
	
	d = hbar*2./(h*h);
	e1 = -hbar/(h*h);
	
	for(i=1;i<n;i++)
	{
		
		r = (i+1)*h;
		vv = V(r);
		gsl_matrix_set(A,i,i,d + vv + hbar*L*(L+1)/(r*r));
		if(i<n-1) gsl_matrix_set(A,i,i+1,e1);
		if(i>0) gsl_matrix_set(A,i,i-1,e1);
	};
	
//	gsl_matrix_fprintf(stdout,A,"%f"); 
	
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv(A, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	
	for(i=0;i<N;i++)
	{
		E[i] = gsl_vector_get(eval,i);
		for(j=0;j<n-1;j++)
		{
			u[i][j] = gsl_matrix_get(evec,j,i);
		}
	}
	return;
}

int main(void)
{
	double pi,hbar,rmax,h,r,unorm,uv;
	int Lmax,L,N,n;
	int i,j,k,NL;
	FILE *out;
	out=fopen("u.out","w+");
	
	printf("Enter max angular momentum: ");
	scanf("%i",&Lmax);
	printf("Enter max principal number: ");
	scanf("%i",&N);
	printf("Enter hbar: ");
	scanf("%lf",&hbar);
	printf("Enter rmax: ");
	scanf("%lf",&rmax);
	printf("Enter the length of the mesh: ");
	scanf("%i",&n);
	NL = N*Lmax;
	double E[N],u[N][n];
	double UT[Lmax*N][n+1],EET[NL];
	double U[Lmax*N][n+1],EE[NL];
	int LL[NL],NN[NL];
	
	h=rmax/n;
	pi=4.*atan(1.);
	
	k=0;
	for(L=0;L<Lmax;L++)
	{
		Hsolve(n,L,N,hbar,rmax,E,u);
		for(i=0;i<N;i++)
		{
			EET[k]=E[i];
			if(L==0)
			{
				UT[k][0]=u[i][0]/h;
			}
			else
			{
				UT[k][0]=0.;
			}
			for(j=0;j<n;j++)
			{
				r = (j+1)*h;
				UT[k][j+1] = u[i][j]/r;
			}
			LL[k] = L;
			NN[k] = i+1;
			k++;
		}
	}
	
	gsl_permutation *p = gsl_permutation_alloc(NL);
	gsl_vector *eall = gsl_vector_alloc(NL);
	for(i=0;i<NL;i++)
	{
		gsl_vector_set(eall,i,EET[i]);
	}
	gsl_sort_vector_index(p,eall);
	for(i=0;i<NL;i++)
	{
		EE[i] = gsl_vector_get(eall,p->data[i]);
		printf("%3i  N:%i   L:%i   E:%f\n",i+1,NN[p->data[i]],LL[p->data[i]],EE[i]);
		unorm = 0.;
		for(j=0;j<n+1;j++)
		{
			r=j*h;
			uv = UT[p->data[i]][j];
			U[i][j] = uv/h;
			unorm += uv*uv*r*r/h;
			fprintf(out,"%10.5f		%10.5g\n",r,U[i][j]);
		}
		printf("%3i: NORM: %f\n",i,unorm*h);
		fprintf(out,"\n");
	}
	return 0;
}
		
	


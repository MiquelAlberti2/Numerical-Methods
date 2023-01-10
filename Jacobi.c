/*Ejercicio práctica 1, Mètodes Numèrics 2, Jacobi*/
/*MIQUEL ALBERTÍ BINIMELIS*/
/*AYLA AYADI ANTÓN*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-8
#define MAX_ITER 10000


int jacobi(int,double**,double**, double*, double*);
double maxim(double, double);
int A_matrix(int, int, int);
double g_function(int, int, double, int);
double f_function(int, int, double);
double u_function(int, int, int);
void crear_documento(double**, int);

int main(void){
	int i,j,n;
	double **x, **xprev;
	double *radi_espectral, *error_sistema;

	radi_espectral=(double*)malloc(sizeof(double));
	error_sistema=(double*)malloc(sizeof(double));

	printf("Introdueix la n:\n");
	scanf("%d",&n);

	x=(double**)calloc(n,sizeof(double));
	if(x==NULL){
		exit(2);
	}
	for(i=0;i<n;i++){
		x[i]=(double*)calloc(n,sizeof(double));
		if(x[i]==NULL) {
			exit(2);
		}
	}
	
    xprev=(double**)calloc(n,sizeof(double));
	if(xprev==NULL){
		exit(2);
	}
	for(i=0;i<n;i++){
		xprev[i]=(double*)calloc(n,sizeof(double));
		if(xprev[i]==NULL) {
			exit(2);
		}
	}


	int iter = jacobi(n,x,xprev, radi_espectral, error_sistema);

	printf("\nSOLUCIÓ OBTINGUDA (x[i]) \tVALOR REAL (u(x,y))\n");

	double max, dif;

	max=0;

	for(j=0; j<n; j++){
		for(i=0; i<n; i++){
			printf("%15.7e \t\t %15.7e \n",x[i][j], u_function(i+1,j+1,n));

			/*
			APROFITO PER CALCULAR LA NORMA DEL MÀXIM DE LA DIFERÈNCIA
			*/
	        dif = fabs( x[i][j] - u_function(i+1,j+1,n) );
			//actualitzem si trobem un valor major
	        max = maxim(max,dif); //max(max, dif)
		    
		}
	}
	printf("\n");

	printf("\nNombre d'iterats per la convergència de Jacobi : %d \n",iter);
	printf("Error entre la solució obtinguda i la real : %15.7e \n", max);
	printf("Error del mètode iteratiu : %15.7e \n", *error_sistema);
	printf("Aproximació del radi espectral : %15.7e \n\n", *radi_espectral);

	crear_documento(x, n);

	//alliberem memòria dinàmica
	for(i=0;i<n;i++){
		free(x[i]);
	}
	free(x);
    
    for(i=0;i<n;i++){
		free(xprev[i]);
	}
	free(xprev);

	return 0;
}


int jacobi(int n,double **x, double **xprev, double *radi_espectral, double *error_sistema){
	int i,j, i2, j2;
	double **b;
	int n_quadrat = n*n;

	/*
	DEFINIM VECTOR DE TERMES INDEP
	*/

	b=(double**)malloc(n*sizeof(double*));

	if(b==NULL){
		exit(3);
	}

	for(i=0;i<n;i++){
		b[i]=(double*)malloc(n*sizeof(double));
		if(b[i]==NULL) {
			exit(2);
		}
	}


	double h;
	h=1./(n+1);


	//OMPLIM EL "VECTOR" DE TERMES INDEPENDENTS
	for(j=0; j<n; j++){
		for(i=0; i<n; i++){
			//ALERTA: els índex de l'enunciat comencen per 1
			b[i][j]=h*h*f_function(i+1,j+1,h);

			//afegim les falses incògnites al terme independent
	        if(i==0){
	        	b[i][j] += g_function(0,j+1,h,n);
	        }
	        else if(i==n-1){
	        	b[i][j] += g_function(n+1,j+1,h,n);
	        }

	        if(j==0){
	        	b[i][j] += g_function(i+1,0,h,n);
	        }
	        else if(j==n-1){
	        	b[i][j] += g_function(i+1,n+1,h,n);
	        }
		}
	}


	/*
	CALCUL
	*/

	int iter = 0;
	double sum, max;
	double dif, error_ant;
	

    max = 2*TOL; //valor simbòlic per entrar al bucle 

	while(max>TOL && iter<MAX_ITER){
		error_ant=max; //guardo l'error anterior per calcular el quocient entre aquests
		max=0;
		
		for(i=0;i<n_quadrat;i++){
			/* ATENCIO
			els índex "i","j" són els de l'algorisme,
			els índex "i2", "j2", son els índex "i","j" de l'enunciat -1 pq estan dins una matriu
			*/

			//CALCULEM EL SUMATORI DE L'ALGORISME
			sum=0;
			for(j=0;j<n_quadrat;j++){
				if (j!=i){
					//calculem els índex de l'enunciat
			        i2= j%n;
			        j2= j/n;
					sum += A_matrix(i,j,n)*xprev[i2][j2];
				}
	        }
	        
	        //calculem els índex de l'enunciat
	        i2= i%n;
	        j2= i/n;

	        //CALCULEM LA x DE LA NOVA ITERACIÓ APLICANT L'ALGORISME
	        
	        x[i2][j2] = (b[i2][j2]-sum) / 4; //a[i][i]=4
	        
	         //CALCULEM LA NORMA DEL MÀXIM DE LA DIFERÈNCIA PEL CRITERI D'ATURADA
	        dif = fabs(x[i2][j2] - xprev[i2][j2]);
	       
			//actualitzem si trobem un valor major
	        max = maxim(max,dif); //max(max, dif)

		}

		//copiem a xprev el contingut de x per calcular la següent iteració
		for(j=0; j<n; j++){
			for(i=0; i<n; i++){
				xprev[i][j] = x[i][j];
			}
		}

		printf("\nk=%d l'error usant la norma infinit és: %15.9e \n", iter, max);
		if(iter>0){ //en la primera iteració no tenim error amb el que comparar
			printf("quocient de ( error actual / error anterior ) : %15.9e \n", max/error_ant);
		}
		iter++;

	}

	(*radi_espectral)=max/error_ant;
	(*error_sistema)=max;
	
    for(i=0;i<n;i++){
        free(b[i]);
    }
	free(b);

	return iter;
}

int A_matrix(int i, int j, int n){
	//Diagonal
	if (i==j){
		return 4;
	}
	//les bandes més externes diagonals
	else if(i==j+n || i==j-n){
		return -1;
	}
	//les bandes contigües a la diagonal
	else if(i==j+1 || i==j-1){
		if(((i%n==0)  && (j<i)) || ((j%n==0) && (j>i))){
			return 0;
		}
		else{
			return -1;
		}
	}
	else{
		return 0;
	}

}

double maxim(double x, double y){
	if (x > y){
		return x;
	}
	else{
		return y;
	}
}

double u_function(int i, int j, int n){
	double h;
	double x,y;
	h=1./(n+1);
	x=i*h;
    y=j*h;

    return x*x*sin(y) - 2*y*y*cos(x);
}

double f_function(int i, int j, double h){

	double x, y;
	x=i*h;
	y=j*h;
	return (-2+x*x)*sin(y) + (-2*y*y+4)*cos(x);
}

double g_function(int i, int j, double h, int n){
	double x,y;

	x = i*h;
	y = j*h;

	if(i==0){ //g(0,y) 
		return -2*y*y;
	}
	else if(j==0){ //g(x,0)
		return 0;		
	}
	else if(i==n+1){ //g(n+1,y) 
		return sin(y)-2*y*y*cos(1);
	}
	else if(j==n+1){ //g(x,n+1) 
		return x*x*sin(1) - 2*cos(x);
	}

	//per retornar alguna cosa i no peti en cas que no s'usi correctamnet
	return -1;
}

void crear_documento(double **x, int n){
	FILE *file=fopen("resultadosJacobi", "w");

	double xi,yi, h;
	h=1./(n+1);

	int i,j;

	//imprimim frontera
	for(i=1;i<n+1;i++){
		xi=i*h;
		fprintf(file, "%lf \t %lf \t %lf \n", (double)0, xi, g_function(0,i,h,n));
		fprintf(file, "%lf \t %lf \t %lf \n", (double)1, xi, g_function(n+1,i,h,n));
		fprintf(file, "%lf \t %lf \t %lf \n", xi, (double)0, g_function(i,0,h,n));
		fprintf(file, "%lf \t %lf \t %lf \n", xi, (double)1, g_function(i,n+1,h,n));
	}
	fprintf(file, "%lf \t %lf \t %lf \n", (double)0, (double)0, g_function(0,0,h,n));
	fprintf(file, "%lf \t %lf \t %lf \n", (double)1, (double)0, g_function(n+1,0,h,n));
	fprintf(file, "%lf \t %lf \t %lf \n", (double)0, (double)1, g_function(0,n+1,h,n));
	fprintf(file, "%lf \t %lf \t %lf \n", (double)1, (double)1, g_function(n+1,n+1,h,n));

	//imprimim interior
	for(j=0;j<n;j++){
		for(i=0;i<n;i++){
			xi=(i+1)*h;
			yi=(j+1)*h;
			fprintf(file, "%lf \t %lf \t %lf \n", xi, yi, x[i][j]);
		}
	}

	fclose(file);
}

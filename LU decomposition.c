#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void calculMatrixULP(double** u, double** l, int* p, int n);
double* resolucioTrianInf(double** m, double* b, int n);
double* resolucioTrianSup(double** m, double* b, int n);

void permutacioFila(double** matriu, int i, int j, int n);
void permutacioVectorInt(int* vector, int i, int j);
void permutacioVectorDouble(double* vector, int i, int j);

int main(void){
	int n, *p;
	double **a, *b, *pb, *x, *y, **u, **l;
	
	//printf("Numero de ecuaciones: ");
	//scanf("%d", &n);
	n=3;

	a=(double**) malloc (n*sizeof(double*));
	u=(double**) malloc (n*sizeof(double*));
	l=(double**) malloc (n*sizeof(double*));

	if (a==NULL || u==NULL || l==NULL){
		printf("No hi ha prou memoria\n");
		exit(1);
	}

	for (int i=0;i<n;i++){
		a[i]=(double*) malloc (n*sizeof(double));
		u[i]=(double*) malloc (n*sizeof(double));
		l[i]=(double*) malloc (n*sizeof(double));
		if (a[i] == NULL || u[i] == NULL || l[i] == NULL){
			printf("No hi ha prou memoria\n");
			exit(2);
		}
	}

	b=(double*) malloc (n*sizeof(double));
	pb=(double*) malloc (n*sizeof(double));
	x=(double*) malloc (n*sizeof(double));
	y=(double*) malloc (n*sizeof(double));
	p=(int*) malloc (n*sizeof(int));

	if (b==NULL || pb==NULL || x==NULL || p==NULL || y==NULL){
		printf ("No hi ha prou memoria\n");
		exit(3);
	}

	for(int i=0; i<n; i++){
		p[i]=i+1;
	}
	
/* PER DEMANAR PER TECLAT LA MATRIU
	printf("\nMatriz del sistema:\n");
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			printf("a(%d, %d): ", i+1, j+1);
			scanf("%lf", &a[i][j]);
			u[i][j]=a[i][j];
		}
	}
	
	printf("\n\nTerminos independientes:\n");
	for (int i=0; i<n; i++){
		printf("b(%d): ", i+1);
		scanf("%lf", &b[i]);
	}
*/
	
//*/PER HARDCODEAR	
	a[0][0]=1 ;
	a[0][1]=4 ;
	a[0][2]=9 ;
	a[1][0]=4 ;
	a[1][1]=14 ;
	a[1][2]=30 ;
	a[2][0]=-2;
	a[2][1]=-2 ;
	a[2][2]=-3;

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			u[i][j]=a[i][j];
		}
	}

	b[0]=17;
	b[1]=58;
	b[2]=-2;
//*/

	//Fins aquí la preparació de les matrius a, b, x, p, u
	//La matriu u ara és una copia d'a, i sobre ella faré les operacions per no perdre l'original
	

	calculMatrixULP(u, l, p, n);

	//calcul Pb
	for(int i=0; i<n; i++){
		pb[i]=b[p[i]-1];
	}
		
	y=resolucioTrianInf(l, pb, n);

/*
	for(int i=0; i<n; i++){
		printf("%f ", y[i]);
	}
	printf(" \n ");
*/
	
	x=resolucioTrianSup(u, y, n);

/*
	printf("MATRIZ U: \n");
	for (int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			printf(" %f", u[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("MATRIZ L: \n");
	for (int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			printf(" %f", l[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
*/
	printf("SOLUCIONS DEL SISTEMA: \n");
	for (int i=0; i<n; i++){
		printf("X%d = %f \n", i+1, x[i]);
	}
	

	//Alliberam la memòria dinàmica
	for(int i=0; i<n; i++){
		free(a[i]);
		free(u[i]);
		free(l[i]);
	}

	free(a);
	free(u);
	free(l);
	free(b);
	free(x);
	free(y);
	free(p);

	
	return 0;
}


void calculMatrixULP(double** u, double** l, int* p, int n){
	double pivot;
	int filaPiv;
	for(int i=0; i<n-1; i++){ //OJO, FAIG N-1 PASSOS
		pivot=u[i][i];
		filaPiv=i;
		for(int j=i+1; j<n; j++){
			if(fabs(pivot)<fabs(u[j][i])){ //ojo, en valor absolut
				pivot=u[j][i];
				filaPiv=j;
			}
		}//ara ja tinc que pivot és l'element màxim de la columna
		
//OJO NO POSAR EL ZERO, USAR LA TOLERANCIA

		if(pivot==0){ 
			printf("Error pivot==0\n");
			exit(4);
		}

		//faig la permutació si cal
		if(filaPiv!=i){
			permutacioFila(u, filaPiv, i, n);
			permutacioVectorInt(p, filaPiv, i);
		}	
			
		//resto a cada fila un múltiple de la fila i
		for(int j=i+1; j<n; j++){
			u[j][i]= u[j][i]/pivot; //guardo en la mateixa matriu u els multiplicadors(ja considero que son 0)
			for(int k=i+1; k<n; k++){
				u[j][k]=u[j][k] - u[j][i]*u[i][k];
			}
		}	
	}

	//Escric be la L i la U
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(j==i){
				l[i][j]=1;
			}
			else if(j<i){
				l[i][j]=u[i][j];
				u[i][j]=0;
			}
			else{
				l[i][j]=0;
			}
		}
	}
}

double* resolucioTrianInf(double** m, double* b, int n){
	double sum;
	double *x;
	x=(double*) malloc (n*sizeof(double));
	
	//No cal comprovar pivot != 0 pq ja m'ho garanteix la factorització LU
	for (int i=0; i<n; i++){		

		sum=0.0;
		for(int j=0; j<i; j++){
			sum+=m[i][j]*x[j];
		}
		x[i]=(b[i]-sum)/m[i][i];
	}
	return x;
}

double* resolucioTrianSup(double** m, double* b, int n){
	double sum;
	double *x;
	x=(double*) malloc (n*sizeof(double));
	
	//No cal comprovar pivot != 0 pq ja m'ho garanteix la factorització LU
	for (int i=n-1; i>=0; i--){		

		sum=0.0;
		for(int j=i+1; j<n; j++){
			sum+=m[i][j]*x[j];
		}
		x[i]=(b[i]-sum)/m[i][i];
	}
	return x;
}

void permutacioFila(double** matriu, int i, int j, int n){
	double *aux;
	aux=(double*) malloc (n*sizeof(double)); //li assigno memòria a aux
	
	for(int z=0; z<n; z++){ //copio el contingut
		aux[z]=matriu[i][z];
	}

	free(matriu[i]);

	matriu[i]=matriu[j]; //aquí només moc referències
	matriu[j]=aux;
}

void permutacioVectorInt(int* vector, int i, int j){
	int aux;

	aux=vector[i];
	vector[i]=vector[j];
	vector[j]=aux;
	
}

void permutacioVectorDouble(double* vector, int i, int j){
	double aux;

	aux=vector[i];
	vector[i]=vector[j];
	vector[j]=aux;
	
}

/*Ejercicio práctica 2, Mètodes Numèrics 2*/
/*MIQUEL ALBERTÍ BINIMELIS*/
/*AYLA AYADI ANTÓN*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A -0.8
#define B 0.5
#define C -1.2
#define D 0.3

#define PRE 1e-8
#define h 1e-2
#define ITER 5
#define TOL 1e-6
#define n 100

double funcion(double x, double y){
    double i = 16*pow(x,4) + pow(y,4) + A*x*y*y - 16*x*x*y - 1;
    double j = x*x + (y-1)*(y-1) + B*x*y + C;

	return i*j+D;
}

double funcion_dx(double x, double y){
    double i = 16*pow(x,4) + pow(y,4) + A*x*y*y - 16*x*x*y - 1;
    double j = x*x + (y-1)*(y-1) + B*x*y + C;

    double di = 64*pow(x,3) + A*y*y - 32*x*y;
    double dj = 2*x + B*y;

	return di*j + i*dj;
}

double funcion_dy(double x, double y){
	double i = 16*pow(x,4) + pow(y,4) + A*x*y*y - 16*x*x*y - 1;
    double j = x*x + (y-1)*(y-1) + B*x*y + C;

    double di = 4*pow(y,3) + 2*A*x*y - 16*x*x;
    double dj = 2*(y-1) + B*x;

	return di*j + i*dj;
}

void aprox_inicial_zeros_gnuplot();


typedef struct Punt{
    double x;
    double y;
} punt;


punt prediccio(punt, int);
punt correccio(punt, punt);
punt correccio_inicial(punt);
punt inversa_producto(punt, punt);
punt corba(double, punt);




int main(void){
    FILE *file=fopen("resultados_metode_cont", "w");
    
    punt aprox_inicial;

    //pts que hem trobat amb la funció aprox_inicial_zeros_gnuplot()
    aprox_inicial.x = 0.17  ;
    aprox_inicial.y = -0.862  ;
    
    printf("f(punt inicial)= %lf\n", funcion(aprox_inicial.x, aprox_inicial.y));

    aprox_inicial = correccio_inicial(aprox_inicial);
    printf("punt inicial després de la correccio: (%lf, %lf)\n", aprox_inicial.x, aprox_inicial.y);
    printf("f(punt inicial corregit)= %lf\n", funcion(aprox_inicial.x, aprox_inicial.y));

    int i;
    
    punt pt = aprox_inicial;
    punt pt_anterior;

    //calculem n pts en una direcció
    for(i=0; i<n; i++){
        pt_anterior = pt;
        if(fabs(sqrt((funcion_dx(pt.x, pt.y)*funcion_dx(pt.x, pt.y)+funcion_dy(pt.x, pt.y)*funcion_dy(pt.x, pt.y)))) < TOL){
            printf("Punt singular\n");
            break;
        }

        pt=prediccio(pt, 1);
        printf("\n\nf(prediccio)= %lf\n", funcion(pt.x, pt.y));
        pt=correccio(pt, pt_anterior);
        printf("f(correccio)= %lf\n", funcion(pt.x, pt.y));

        
        fprintf(file, "%lf \t %lf \t %lf \n", pt.x,pt.y,funcion(pt.x, pt.y));
    }
    
    fprintf(file, "\n"); //escrivim una línia en blanc pq el plot quedi bé
    pt = aprox_inicial;

    //calculem n pts en l'altra direcció
    for(i=0; i<n; i++){
        pt_anterior = pt;
        if(fabs(sqrt((funcion_dx(pt.x, pt.y)*funcion_dx(pt.x, pt.y)+funcion_dy(pt.x, pt.y)*funcion_dy(pt.x, pt.y)))) < TOL){
            printf("Punt singular\n");
            break;
        }

        pt=prediccio(pt, -1);
        printf("\n\nf(prediccio)= %lf\n", funcion(pt.x, pt.y));
        pt=correccio(pt, pt_anterior);
        printf("f(correccio)= %lf\n", funcion(pt.x, pt.y));

        
        fprintf(file, "%lf \t %lf \t %lf \n", pt.x,pt.y,funcion(pt.x, pt.y));
    }
    
    fclose(file);
    
    return 0;
}



punt prediccio(punt inici, int direccio){
    punt pred;
    double norma;

    norma = sqrt(funcion_dy(inici.x,inici.y)*funcion_dy(inici.x,inici.y) + funcion_dx(inici.x,inici.y)*funcion_dx(inici.x,inici.y));
    //ja hem verificat abans que aquesta es diferent de 0
    
    pred.x = inici.x + direccio*h* (-1*funcion_dy(inici.x,inici.y)/norma);
    pred.y = inici.y + direccio*h* funcion_dx(inici.x,inici.y)/norma;
      
    return pred;
}

punt correccio(punt prediccio, punt inici){
    int i=0;
    punt aux;
    while(i<ITER && fabs(funcion(prediccio.x, prediccio.y))>=PRE){
        aux = inversa_producto(inici, prediccio);

        if (fabs(aux.x)<TOL && fabs(aux.y)<TOL){
            //vol dir que hem tingut un denominador = 0 o que la iteració de Newton ens ha retornat el mateix pt
            //en ambdós casos cal aturar Newton
            break;
        }

        prediccio.x=prediccio.x-aux.x;
        prediccio.y=prediccio.y-aux.y;
        
        i++;
    }
    
    if(i==ITER){
        printf("Newton no ha convergit\n");
    }
    
    return prediccio;
}

punt inversa_producto(punt inici, punt prediccio){
    double a = funcion_dx(prediccio.x, prediccio.y);
    double b = funcion_dy(prediccio.x, prediccio.y);
    double c = 2*(prediccio.x - inici.x);
    double d = 2*(prediccio.y - inici.y);
    double e = funcion(prediccio.x, prediccio.y);
    double f = (prediccio.x - inici.x)*(prediccio.x - inici.x) + (prediccio.y - inici.y)*(prediccio.y - inici.y) - h*h;

    punt result;

    double den = a*d - b*c;

    if(fabs(den) >= TOL){
        result.x = (d*e - b*f) / den;
        result.y = (-c*e + a*f) / den;
    }
    else{ //ho indiquem així perquè aturi Newton
        result.x=0;
        result.y=0;
    }
    
    return result;
}


punt correccio_inicial(punt inici){
    double t=0;
    double num, den;
    punt pt = inici;
    
    int i=0;
    while(i<ITER && fabs(funcion(pt.x, pt.y)) >= PRE){
        
        
        num = funcion(pt.x, pt.y);
        den = funcion_dx(pt.x, pt.y)*funcion_dx(inici.x, inici.y) + funcion_dy(pt.x, pt.y)*funcion_dy(inici.x, inici.y);

        if(fabs(den)<TOL){
            printf("Error en la correccio inicial\n");
            break;
        }
        
        t = t - (num/den);
        
        pt=corba(t, inici);
        i++;
    }
    
    if(i==ITER){
        printf("Newton pel punt inicial no ha convergit\n");
    }
    
    return pt;
}

punt corba(double t, punt inici){
    punt resultat;
    resultat.x = inici.x + t*funcion_dx(inici.x, inici.y);
    resultat.y = inici.y + t*funcion_dy(inici.x, inici.y);
    
    return resultat;
}


/*
 * PROGRAMA PER IMPRIMIR UNA PRIMERA APROXIMACIÓ DELS ZEROS

void aprox_inicial_zeros_gnuplot(){
    FILE *file=fopen("resultados_aprox_inicial", "w");
    
    double i,j;
    for(i=-10 ; i<=10 ; i+=0.002){
        for(j=-10 ; j<=10 ; j+=0.002){
            if(fabs(funcion(i,j)) < 1e-2){
                fprintf(file, "%lf \t %lf \t %lf \n", i,j,funcion(i,j));
            }
        }
    }
    
    fclose(file);
} */

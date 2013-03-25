#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define GRAV_CONSTANT 66.7E-12 /* (meter^3)/(kilogram*second^2) */ 
#define MSTOMGAL 1E5 /* converts meter per squared second to miliGal */

#define DIMENSAO_STRING 100

double **aloca_matriz_double (FILE *, int, int);

double **libera_matriz_double (int, double **);

double *aloca_vetor_double (FILE *, int);

double *libera_vetor_double (double *);

void leitura_parametros (FILE *relatorio, int *M, double *Z_layer, char arquivo_camada[DIMENSAO_STRING], double *declination_earth, double *inclination_earth, double *declination_body, double *inclination_body, int *continuacao_interpolacao, char arquivo_continuacao_interpolacao[DIMENSAO_STRING], int *tensor, char arquivo_tensor[DIMENSAO_STRING]);

void leitura_camada(FILE *relatorio, int M, char arquivo_camada[DIMENSAO_STRING], double *Y, double *X, double *p);

void gz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gz);

void gxx_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gxx);

void gyy_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gyy);

void gzz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gzz);

void gxy_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gxy);

void gxz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gxz);

void gyz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gyz);

double gz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double gxx_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double gyy_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double gzz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double gxy_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double gxz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double gyz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

void impressao_componente(FILE *arquivo, int N, double *y, double *x, double *z, double *componente);

void leitura_arquivo_dados_calculados1(FILE *relatorio, FILE *entrada, char arquivo_continuacao_interpolacao[DIMENSAO_STRING], int N, int *ny, int *nx, int *nz, double *ymin, double *xmin, double *zmin, double *dy, double *dx, double *dz);

void leitura_arquivo_dados_calculados2(FILE *relatorio, FILE *entrada, char arquivo_continuacao_interpolacao[DIMENSAO_STRING], int N, double *yp, double *xp, double *zp);

void constroi_grid(int ny, int nx, int nz, double ymin, double xmin, double zmin, double dy, double dx, double dz, double *yp, double *xp, double *zp);

double dipolo_mag (int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer);

double dipolo_mag2 (int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer, double ly, double lx, double lz, double Ly, double Lx, double Lz);

void mag_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *p, double *gz);

void main () {

	int N, M;
	int i, j, k;
	int nx, ny, nz;
	int continuacao_interpolacao, tensor;
	double *xp, *yp, *zp, *dados_calculados, dx, dy, dz, xmin, ymin, zmin;
	double declination_earth, inclination_earth, Lx, Ly, Lz, declination_body, inclination_body, lx, ly, lz;
	double *X, *Y, Z_layer, *p;
	double aux0, aux1, aux2;
	double stdev_dados;
	char arquivo_camada[DIMENSAO_STRING], arquivo_continuacao_interpolacao[DIMENSAO_STRING], arquivo_tensor[DIMENSAO_STRING];

	time_t start_global, end_global, start_local, end_local;

	FILE *entrada, *saida, *relatorio;
	
	time (&start_global);

	/*

	Sistema de coordenadas
	**********************

	x aponta para norte (metros)
	y aponta para leste (metros)
	z aponta para baixo (metros)

	Descrição das variáveis
	***********************

	M = Número total de fontes equivalentes.
	X = Vetor que armazena as coordenadas X das fontes equivalentes.
	Y = Vetor que armazena as coordenadas Y das fontes equivalentes.
	Z_layer = Coordenada Z_layer das fontes equivalentes.
	p = Vetor que armazena a densidade de cada fonte equivalente.

	N = Número total de dados calculados.
	*xp = Vetor que armazena as coordenadas x dos dados calculados.
	*yp = Vetor que armazena as coordenadas y dos dados calculados.
	*zp = Vetor que armazena as coordenadas z dos dados calculados.
	*dados_calculados = Vetor que amazena os dados calculados.
	nx = Número de pontos, na direção x de um grid, em que os dados serão calculados.
	ny = Número de pontos, na direção y de um grid, em que os dados serão calculados.
	nz = Número de pontos, na direção z de um grid, em que os dados serão calculados.
	xmin = Coordenada mais a sul do grid de dados calculados.
	ymin = Coordenada mais a oeste do grid de dados calculados.
	zmin = Coordenada mais rasa do grid de dados calculados.
	dx = Intervalo entre os pontos, na direção x, do grid de dados calculados.
	dy = Intervalo entre os pontos, na direção y, do grid de dados calculados.
	dz = Intervalo entre os pontos, na direção z, do grid de dados calculados.

	continuacao_interpolacao = Label que define os parâmetros da continuação/interpolação:
		(0) Não faz continuação/interpolação;
		(1) Faz continuação/interpolação em um grid;
		(2) Faz continuação/interpolação em um conjunto de pontos existente em um arquivo;
	tensor = Label que define os parâmetros do cálculo das componentes do tensor:
		(0) Não calcula as componentes do tensor;
		(1) Calcula as componentes do tensor em um grid;
		(2) Calcula as componentes do tensor em um conjunto de pontos existente em um arquivo;

	Descrição breve do programa
	***************************

	Este programa lê uma camada equivalente e calcula alguma transformação.
	Esta transformação pode ser: interpolação em um conjunto de pontos estruturado (grid)
	ou não, continuação para cima ou para baixo ou cálculo das componentes do tensor.

	*/

	relatorio = fopen("relatorio.txt", "w");

	leitura_parametros(relatorio, &M, &Z_layer, arquivo_camada, &declination_earth, &inclination_earth, &declination_body, &inclination_body, &continuacao_interpolacao, arquivo_continuacao_interpolacao, &tensor, arquivo_tensor);

	declination_earth = (double)((3.14159265358979323846*declination_earth)/180.0);
	inclination_earth = (double)((3.14159265358979323846*inclination_earth)/180.0);
	declination_body  = (double)((3.14159265358979323846*declination_body)/180.0);
	inclination_body  = (double)((3.14159265358979323846*inclination_body)/180.0);
	
	Lx = cos(inclination_earth)*cos(declination_earth);
	Ly = cos(inclination_earth)*sin(declination_earth);
	Lz = sin(inclination_earth);

	lx = cos(inclination_body)*cos(declination_body);
	ly = cos(inclination_body)*sin(declination_body);
	lz = sin(inclination_body);

	X = aloca_vetor_double(relatorio, M);
	Y = aloca_vetor_double(relatorio, M);
	p = aloca_vetor_double(relatorio, M);

	leitura_camada(relatorio, M, arquivo_camada, Y, X, p);

	if ((continuacao_interpolacao != 0) && ((continuacao_interpolacao == 1) || (continuacao_interpolacao == 2))) {

		if (fopen (arquivo_continuacao_interpolacao, "r") == NULL) {
		
			fprintf (relatorio, "Arquivo %s nao encontrado!!\n\n", arquivo_continuacao_interpolacao);
			
			fclose (relatorio);

			abort();
		
		}
		
		entrada = fopen (arquivo_continuacao_interpolacao, "r");
	
		fscanf (entrada, "%d", &N);

		yp = aloca_vetor_double(relatorio, N);
		xp = aloca_vetor_double(relatorio, N);
		zp = aloca_vetor_double(relatorio, N);

		if (continuacao_interpolacao == 1) {
		
			leitura_arquivo_dados_calculados1(relatorio, entrada, arquivo_continuacao_interpolacao, N, &ny, &nx, &nz, &ymin, &xmin, &zmin, &dy, &dx, &dz);

			constroi_grid(ny, nx, nz, ymin, xmin, zmin, dy, dx, dz, yp, xp, zp);
			
		}
		
		if (continuacao_interpolacao == 2) {

			leitura_arquivo_dados_calculados2(relatorio, entrada, arquivo_continuacao_interpolacao, N, yp, xp, zp);
			
		}
		
		fclose (entrada);

		dados_calculados = aloca_vetor_double(relatorio, N);

		/* Cálculo da continuação/interpolação */
		//gz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		mag_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, Ly, Lx, Lz, ly, lx, lz, p, dados_calculados);
		
		/* Impressão do arquivo de saida */
		saida = fopen("t_calculado.txt", "w");
		impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		fclose(saida);
		
		dados_calculados = libera_vetor_double(dados_calculados);

		yp = libera_vetor_double(yp);
		xp = libera_vetor_double(xp);
		zp = libera_vetor_double(zp);

	}

	if ((tensor != 0) && ((tensor == 1) || (tensor == 2))) {
	
		//i = strcmp(arquivo_continuacao_interpolacao, arquivo_tensor);
	
		//if (i != 0) {
	
			if (fopen (arquivo_tensor, "r") == NULL) {
			
				fprintf (relatorio, "Arquivo %s nao encontrado!!\n\n", arquivo_tensor);
				
				fclose (relatorio);

				abort();
			
			}
			
			entrada = fopen (arquivo_tensor, "r");
		
			fscanf (entrada, "%d", &N);

			yp = aloca_vetor_double(relatorio, N);
			xp = aloca_vetor_double(relatorio, N);
			zp = aloca_vetor_double(relatorio, N);

			if (tensor == 1) {
			
				leitura_arquivo_dados_calculados1(relatorio, entrada, arquivo_tensor, N, &ny, &nx, &nz, &ymin, &xmin, &zmin, &dy, &dx, &dz);

				constroi_grid(ny, nx, nz, ymin, xmin, zmin, dy, dx, dz, yp, xp, zp);
				
			}
			
			if (tensor == 2) {

				leitura_arquivo_dados_calculados2(relatorio, entrada, arquivo_tensor, N, yp, xp, zp);
				
			}
			
			fclose (entrada);
		
			dados_calculados = aloca_vetor_double(relatorio, N);

			/* Cálculo das componentes do tensor */
			gxx_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
			
			/* Impressão do arquivo de saida */
			saida = fopen("gxx_calculado.txt", "w");
			impressao_componente(saida, N, yp, xp, zp, dados_calculados);
			fclose(saida);

			/* Cálculo das componentes do tensor */
			gxy_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
			
			/* Impressão do arquivo de saida */
			saida = fopen("gxy_calculado.txt", "w");
			impressao_componente(saida, N, yp, xp, zp, dados_calculados);
			fclose(saida);

			/* Cálculo das componentes do tensor */
			gxz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
			
			/* Impressão do arquivo de saida */
			saida = fopen("gxz_calculado.txt", "w");
			impressao_componente(saida, N, yp, xp, zp, dados_calculados);
			fclose(saida);

			/* Cálculo das componentes do tensor */
			gyy_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
			
			/* Impressão do arquivo de saida */
			saida = fopen("gyy_calculado.txt", "w");
			impressao_componente(saida, N, yp, xp, zp, dados_calculados);
			fclose(saida);

			/* Cálculo das componentes do tensor */
			gyz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
			
			/* Impressão do arquivo de saida */
			saida = fopen("gyz_calculado.txt", "w");
			impressao_componente(saida, N, yp, xp, zp, dados_calculados);
			fclose(saida);

			/* Cálculo das componentes do tensor */
			gzz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
			
			/* Impressão do arquivo de saida */
			saida = fopen("gzz_calculado.txt", "w");
			impressao_componente(saida, N, yp, xp, zp, dados_calculados);
			fclose(saida);

			dados_calculados = libera_vetor_double(dados_calculados);
		
		//}
		//
		//else {
		//
		//	dados_calculados = aloca_vetor_double(relatorio, N);

		//	/* Cálculo das componentes do tensor */
		//	gxx_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		//	
		//	/* Impressão do arquivo de saida */
		//	saida = fopen("gxx_calculado.txt", "w");
		//	impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		//	fclose(saida);

		//	/* Cálculo das componentes do tensor */
		//	gyy_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		//	
		//	/* Impressão do arquivo de saida */
		//	saida = fopen("gyy_calculado.txt", "w");
		//	impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		//	fclose(saida);

		//	/* Cálculo das componentes do tensor */
		//	gzz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		//	
		//	/* Impressão do arquivo de saida */
		//	saida = fopen("gzz_calculado.txt", "w");
		//	impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		//	fclose(saida);

		//	/* Cálculo das componentes do tensor */
		//	gxy_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		//	
		//	/* Impressão do arquivo de saida */
		//	saida = fopen("gxy_calculado.txt", "w");
		//	impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		//	fclose(saida);

		//	/* Cálculo das componentes do tensor */
		//	gxz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		//	
		//	/* Impressão do arquivo de saida */
		//	saida = fopen("gxz_calculado.txt", "w");
		//	impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		//	fclose(saida);

		//	/* Cálculo das componentes do tensor */
		//	gyz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, dados_calculados);
		//	
		//	/* Impressão do arquivo de saida */
		//	saida = fopen("gyz_calculado.txt", "w");
		//	impressao_componente(saida, N, yp, xp, zp, dados_calculados);
		//	fclose(saida);
		//	
		//	dados_calculados = libera_vetor_double(dados_calculados);
		//
		//}
	
	}

	time (&end_global);
	
	aux0 = difftime (end_global, start_global);

	printf("\nPrograma finalizado com sucesso em %.3lf segundo(s)!\n\n", aux0);
	
	fprintf(relatorio, "\nPrograma finalizado com sucesso em %.3lf segundo(s)!\n\n", aux0);
	
	fclose (relatorio);

	system ("PAUSE");

}

double **aloca_matriz_double (FILE *arq, int linha, int coluna) {

    double **m;  /* ponteiro para a matriz */
    int i;

    /* aloca as linhas da matriz */

    m = (double **)calloc(linha, sizeof(double *));

    if (m == NULL) {

        fprintf (arq, "Memoria Insuficiente (linhas)!\n\n");

        fclose (arq);

		abort();

    }

    /* aloca as colunas da matriz */

    for (i = 0; i < linha; i++ ) {

        m[i] = (double *)calloc(coluna, sizeof(double));

        if (m[i] == NULL) {

            fprintf (arq, "Memoria Insuficiente (colunas)!\n\n");

            fclose (arq);

            abort();

        }

    }

    return (m); /* retorna o ponteiro para a matriz */

}


double **libera_matriz_double (int linha, double **m) {
      
    int  i;  /* variavel auxiliar */
    
    if (m == NULL) { 
          
        return (NULL);
        
    }
    
    for (i = 0; i < linha; i++) { 
        
        free (m[i]); /* libera as linhas da matriz */
        
    }
    
    free (m); /* libera a matriz */
        
    return (NULL); /* retorna um ponteiro nulo */
    
}

double *aloca_vetor_double (FILE *arq, int tamanho) {
       
    double *v; /* ponteiro para o vetor */
    
    v = (double *)calloc(tamanho, sizeof(double));
        
    if (v == NULL) { /*** verifica se há memória suficiente ***/
          
        fprintf (arq, "Memoria Insuficiente!\n\n");
        
        fclose (arq);
        
        abort();
        
    }
     
    return (v); /* retorna o ponteiro para o vetor */
       			           
}

double *libera_vetor_double (double *v) {
      
    if (v == NULL) {
          
        return (NULL);
        
    }
    
	free(v); /* libera o vetor */
    
	return (NULL); /* retorna o ponteiro */

}

void leitura_parametros (FILE *relatorio, int *M, double *Z_layer, char arquivo_camada[DIMENSAO_STRING], double *declination_earth, double *inclination_earth, double *declination_body, double *inclination_body, int *continuacao_interpolacao, char arquivo_continuacao_interpolacao[DIMENSAO_STRING], int *tensor, char arquivo_tensor[DIMENSAO_STRING]) {

	char str[20];
	
	FILE *entrada;

	sprintf (str, "parametros.txt");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		abort();

	}

	entrada = fopen(str, "r");

	if (fscanf (entrada, "%d %lf %s %lf %lf %lf %lf %d %s %d %s", M, Z_layer, arquivo_camada, declination_earth, inclination_earth, declination_body, inclination_body, continuacao_interpolacao, arquivo_continuacao_interpolacao, tensor, arquivo_tensor) != 11) {
	
		fprintf (relatorio, "Erro de leitura no arquivo %s!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		abort();
	
	}

	fclose (entrada);

}

void leitura_camada(FILE *relatorio, int M, char arquivo_camada[DIMENSAO_STRING], double *Y, double *X, double *p) {

	int i;
	
	FILE *entrada;

	if (fopen(arquivo_camada, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", arquivo_camada);

		fclose (relatorio);

		printf ("Erro!\n\n");

		abort();

	}

	entrada = fopen(arquivo_camada, "r");

	for (i = 0; i < M; i++) {

		if (fscanf(entrada, "%lf %lf %lf", &Y[i], &X[i], &p[i]) != 3) {

			fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", arquivo_camada);

			fclose (relatorio);

			printf ("Erro!\n\n");

			abort();

		}

	}
	
	fclose (entrada);

}

void gz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gz) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gz[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gz[i] += gz_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) { 

	double aux0, aux1, aux2, aux3, aux4, aux5, aux6;

	aux1 = (xp[i] - X[j]); /* x */
	aux2 = (yp[i] - Y[j]); /* y */
	aux3 = (zp[i] - Z_layer); /* z */
				
	aux4 = pow(aux1, 2) + pow(aux2, 2) + pow(aux3, 2); /* r^2 */
	aux5 = pow(aux4, 1.5); /* r^3 */

	/* -z/r^3 */			
	aux0 = (double)((-aux3)/aux5);

	return aux0*GRAV_CONSTANT*MSTOMGAL;

}

void gxx_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gxx) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gxx[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gxx[i] += gxx_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gxx_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) { 

	double aux0, aux1, aux2, aux3, aux4, aux5, aux6;

	aux1 = pow((xp[i] - X[j]), 2); /* x^2 */
	aux2 = pow((yp[i] - Y[j]), 2); /* y^2 */
	aux3 = pow((zp[i] - Z_layer), 2); /* z^2 */
			
	aux4 = aux1 + aux2 + aux3; /* r^2 */
	aux5 = pow(aux4, 1.5);     /* r^3 */
	aux6 = pow(aux4, 2.5);     /* r^5 */

	/* (3*x^2)/r^5 - (1/r^3)  */ 			
	aux0 = ((double)((3*aux1)/aux6)) - ((double)(1.0/aux5));

	return aux0*GRAV_CONSTANT*1E9;

}

void gyy_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gyy) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gyy[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gyy[i] += gyy_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gyy_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) { 

	double aux0, aux1, aux2, aux3, aux4, aux5, aux6;

	aux1 = pow((xp[i] - X[j]), 2); /* x^2 */
	aux2 = pow((yp[i] - Y[j]), 2); /* y^2 */
	aux3 = pow((zp[i] - Z_layer), 2); /* z^2 */
			
	aux4 = aux1 + aux2 + aux3; /* r^2 */
	aux5 = pow(aux4, 1.5);     /* r^3 */
	aux6 = pow(aux4, 2.5);     /* r^5 */

	/* (3*y^2)/r^5 - (1/r^3)  */ 			
	aux0 = ((double)((3*aux2)/aux6)) - ((double)(1.0/aux5));

	return aux0*GRAV_CONSTANT*1E9;

}

void gzz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gzz) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gzz[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gzz[i] += gzz_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gzz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) { 

	double aux0, aux1, aux2, aux3, aux4, aux5, aux6;

	aux1 = pow((xp[i] - X[j]), 2); /* x^2 */
	aux2 = pow((yp[i] - Y[j]), 2); /* y^2 */
	aux3 = pow((zp[i] - Z_layer), 2); /* z^2 */
			
	aux4 = aux1 + aux2 + aux3; /* r^2 */
	aux5 = pow(aux4, 1.5);     /* r^3 */
	aux6 = pow(aux4, 2.5);     /* r^5 */

	/* (3*z^2)/r^5 - (1/r^3)  */ 			
	aux0 = ((double)((3*aux3)/aux6)) - ((double)(1.0/aux5));

	return aux0*GRAV_CONSTANT*1E9;

}

void gxy_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gxy) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gxy[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gxy[i] += gxy_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gxy_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) { 

	double aux0, aux1, aux2, aux3, aux4, aux5;

	aux1 = (xp[i] - X[j]); /* x^2 */
	aux2 = (yp[i] - Y[j]); /* y^2 */
	aux3 = (zp[i] - Z_layer); /* z^2 */
			
	aux4 = pow(aux1, 2) + pow(aux2, 2) + pow(aux3, 2); /* r^2 */
	aux5 = pow(aux4, 2.5);     /* r^5 */

	/* (3*x*y)/r^5 */
	aux0 = ((double)((3*aux1*aux2)/aux5));

	return aux0*GRAV_CONSTANT*1E9;

}

void gxz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gxz) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gxz[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gxz[i] += gxz_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gxz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) { 

	double aux0, aux1, aux2, aux3, aux4, aux5;

	aux1 = (xp[i] - X[j]); /* x^2 */
	aux2 = (yp[i] - Y[j]); /* y^2 */
	aux3 = (zp[i] - Z_layer); /* z^2 */
			
	aux4 = pow(aux1, 2) + pow(aux2, 2) + pow(aux3, 2); /* r^2 */
	aux5 = pow(aux4, 2.5);     /* r^5 */

	/* (3*x*z)/r^5 */
	aux0 = ((double)((3*aux1*aux3)/aux5));

	return aux0*GRAV_CONSTANT*1E9;

}

void gyz_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double *p, double *gyz) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gyz[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			gyz[i] += gyz_monopolo(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
		
		}
	
	}

}

double gyz_monopolo(int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) {

	double aux0, aux1, aux2, aux3, aux4, aux5;

	aux1 = (xp[i] - X[j]); /* x^2 */
	aux2 = (yp[i] - Y[j]); /* y^2 */
	aux3 = (zp[i] - Z_layer); /* z^2 */
			
	aux4 = pow(aux1, 2) + pow(aux2, 2) + pow(aux3, 2); /* r^2 */
	aux5 = pow(aux4, 2.5);     /* r^5 */

	/* (3*y*z)/r^5 */
	aux0 = ((double)((3*aux2*aux3)/aux5));

	return aux0*GRAV_CONSTANT*1E9;

}

void impressao_componente(FILE *arquivo, int N, double *y, double *x, double *z, double *componente) {

	int i;

	fprintf (arquivo , "              y              x              z\n");

	for (i = 0; i < N; i++) {

		fprintf (arquivo, "%15.3lf %15.3lf %15.3lf %15.3E\n", y[i], x[i], z[i], componente[i]);

	}

}

void leitura_arquivo_dados_calculados1(FILE *relatorio, FILE *entrada, char arquivo_continuacao_interpolacao[DIMENSAO_STRING], int N, int *ny, int *nx, int *nz, double *ymin, double *xmin, double *zmin, double *dy, double *dx, double *dz) {

	if (fscanf(entrada, "%d %d %d %lf %lf %lf %lf %lf %lf", ny, nx, nz, ymin, xmin, zmin, dy, dx, dz) != 9) {

		fprintf (relatorio, "Erro na leitura no arquivo %s!!!\n\n", arquivo_continuacao_interpolacao);
		
		fclose (relatorio);
		
		abort();

	}

	if (((*ny)*(*nx)*(*nz)) != N) {

		fprintf (relatorio, "Erro na leitura no arquivo %s!!!\n\n", arquivo_continuacao_interpolacao);
		
		fclose (relatorio);
		
		abort();

	}

}

void constroi_grid(int ny, int nx, int nz, double ymin, double xmin, double zmin, double dy, double dx, double dz, double *yp, double *xp, double *zp) {

	int i, j, k, l;
	double y, x, z;
	
	for (l = 0, z = zmin, i = 0; i < nz; i++, z += dz) {
	
		for (x = xmin, j = 0; j < nx; j++, x += dx) {
		
			for (y = ymin, k = 0; k < ny; k++, y += dy, l++) {

				yp[l] = y;
				xp[l] = x;
				zp[l] = z;
			
			}
		
		}
	
	}

}

void leitura_arquivo_dados_calculados2(FILE *relatorio, FILE *entrada, char arquivo_continuacao_interpolacao[DIMENSAO_STRING], int N, double *yp, double *xp, double *zp) {

int i;

	for (i = 0; i < N; i++) {
	
		if (fscanf(entrada, "%lf %lf %lf", &yp[i], &xp[i], &zp[i]) != 3) {
		
			fprintf (relatorio, "Erro na leitura do arquivo %s!!!\n\n", arquivo_continuacao_interpolacao);
			
			fclose (relatorio);
			
			abort();
		
		}
	
	}

}

double dipolo_mag (int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer) {

	double aux0, aux1, aux2, aux3, aux4, aux5, aux6;

	aux1 = (xp[i] - X[j]); /* x */
	aux2 = (yp[i] - Y[j]); /* y */
	aux3 = (zp[i] - Z_layer); /* z */
	aux4 = pow(aux3, 2);  /* z^2 */
				
	aux5 = pow(aux1, 2) + pow(aux2, 2) + aux4; /* r^2 */
	aux6 = pow(aux5, 2.5); /* r^5 */

	/* (3z^2 - r^2)/r^5) */			
	aux0 = (3.0*aux4) - aux5;
	aux0 = (double)(aux0/aux6);

	return aux0*418.87902047863909846;

}

double dipolo_mag2 (int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer, double ly, double lx, double lz, double Ly, double Lx, double Lz) {

	double aux0, aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8;
	double V1, V2, V3, V4, V5, V6;

	aux1 = (xp[i] - X[j]); /* x */
	aux2 = (yp[i] - Y[j]); /* y */
	aux3 = (zp[i] - Z_layer); /* z */
	aux4 = pow(aux1, 2);  /* x^2 */
	aux5 = pow(aux2, 2);  /* y^2 */
	aux6 = pow(aux3, 2);  /* z^2 */
				
	aux7 = aux4 + aux5 + aux6; /* r^2 */
	aux8 = pow(aux7, 2.5); /* r^5 */
	
	V1 = (double)(((3*aux4) - aux7)/aux8);
	V4 = (double)(((3*aux5) - aux7)/aux8);
	V6 = (double)(((3*aux6) - aux7)/aux8);
		
	V3 = (double)((3*aux1*aux3)/aux8);
	V2 = (double)((3*aux1*aux2)/aux8);
	V5 = (double)((3*aux2*aux3)/aux8);

	/* lx*(V1*Lx + V2*Ly + V3*Lz) + ly*(V2*Lx + V4*Ly + V5*Lz) + lz*(V3*Lx + V5*Ly + V6*Lz) */
	aux0  = lx*((V1*Lx) + (V2*Ly) + (V3*Lz));
	aux0 += ly*((V2*Lx) + (V4*Ly) + (V5*Lz));
	aux0 += lz*((V3*Lx) + (V5*Ly) + (V6*Lz));

	return (aux0*100.0); /* 100 = Cm*(10^-9) ==> Cm = (10^-7) henry/meter, (10^-9) transforms Tesla to nano Tesla */

}

void mag_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *p, double *gz) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		gz[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			//gz[i] += dipolo_mag(i, yp, xp, zp, j, Y, X, Z_layer)*p[j];
			gz[i] += dipolo_mag2(i, yp, xp, zp, j, Y, X, Z_layer, Ly, Lx, Lz, ly, lx, lz)*p[j];
		
		}
	
	}

}
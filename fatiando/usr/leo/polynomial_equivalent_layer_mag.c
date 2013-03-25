#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Funções retiradas do Numerical Recipes, utilizadas para contaminar os dados com rúido
aleatório e com distribuição Gaussiana ==> */

	/* variável utilizada nas funções gasdev e ran1, abaixo */

	long idum;

	/* Função gasdev do capítulo sete (Random Numbers) do Numerical Recipes */

	/* Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
	as the source of uniform deviates. */

	float gasdev(long *);

	/* Função ran1 do capítulo sete (Random Numbers) do Numerical Recipes */

	/* "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
	safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
	values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
	successive deviates in a sequence. RNMX should approximate the largest floating value that is
	less than 1. */

	#define IA 16807
	#define IM 2147483647
	#define AM (1.0/IM)
	#define IQ 127773
	#define IR 2836
	#define NTAB 32
	#define NDIV (1+(IM-1)/NTAB)
	#define EPS 1.2e-7
	#define RNMX (1.0-EPS)

	float ran1(long *);

/* <== Funções retiradas do Numerical Recipes, utilizadas para contaminar os dados com rúido
aleatório e com distribuição Gaussiana */


double **aloca_matriz_double (FILE *, int, int);

double **libera_matriz_double (int, double **);

double *aloca_vetor_double (FILE *, int);

double *libera_vetor_double (double *);

double sist_lin_LU (FILE *relatorio, int M, double **GtG, double *Gtd, double *p, double uridge);

void choldc(FILE *relatorio, double **a, int n, double *p);

void cholsl(double **a, int n, double *p, double *b, double *x);

void leitura_parametros (FILE *relatorio, int *N, double *uridge, double *usuavidade, double *Z_layer, int *grau, double *stdev_dados, int *semente, double *declination_earth, double *inclination_earth, double *declination_body, double *inclination_body, int *Npolx, int *Npoly, int *Nspolx, int *Nspoly, double *lambidax, double *lambiday, double *uBtAtAB);

void leitura_dados(FILE *relatorio, int N, double *Xp, double *Yp, double *Zp, double *gobs, double *xmin, double *xmax, double *ymin, double *ymax, double stdev_dados);

void coordenadas_fontes(int Npolx, int Npoly, int Nspolx, int Nspoly, double dx, double dy, double X1, double Y1, double *X, double *Y);

void calcula_BtAtAB (int N, int Q, int Npol, int Nspol, int grau, int M, int H, double *Y, double *X, double Z_layer, double *yp, double *xp, double *zp, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *a, double **B, double *g, double **GtG, double *Gttobs, double *tobs, double *fator_normalizacao);

void calcula_B (int Npol, int Nspol, int grau, double *Y, double *X, double **B, double Y1, double Y2, double X1, double X2);

void calcula_BtRtRB (int H, int Q, int Npolx, int Npoly, int Nspol, int Nsy, int Nspolx, int Nspoly, double **B, double *g, double **GtG, double fator_normalizacao, double usuavidade);

void incorpora_ridge (int H, double **GtG, double uridge, double fator_normalizacao);

void calcula_magnetizacoes(int M, int H, double *c, double *p, double **B);

void linha_a (int i, int M, double *Y, double *X, double Z_layer, double *yp, double *xp, double *zp, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *a);

void mag_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *p, double *t);

double dipolo_mag2 (int i, double *yp, double *xp, double *zp, int j, double *Y, double *X, double Z_layer, double ly, double lx, double lz, double Ly, double Lx, double Lz);

void impressao_componente(FILE *arquivo, int N, double *y, double *x, double *z, double *tobs, double *componente);

void impressao_camada(FILE *arquivo, int M, double *Y, double *X, double *p);

void main () {

	int N, M, Npol, Nspol, Npoly, Npolx, Nsy, Nsx, Nspoly, Nspolx, grau, Q, H;
	int i, j, k;
	int semente;
	double *tobs, *xp, *yp, *zp, xmin, xmax, ymin, ymax;
	double Y1, X1, Y2, X2, lambiday, lambidax, dy, dx;
	double *t, *X, *Y, Z_layer, *c, *p;
	double declination_earth, inclination_earth, Lx, Ly, Lz, declination_body, inclination_body, lx, ly, lz;
	double *a, **B, *g;
	double aux0, aux1, aux2;
	double uridge, usuavidade, fator_normalizacao, **BtRtRB, **BtAtAB, **GtG, *Gttobs, uBtAtAB;
	double stdev_dados;
	char str[100];

	time_t start1, end1, start2, end2;

	FILE *entrada, *saida, *relatorio;
	
	time (&start1);

	/*

	Sistema de coordenadas
	**********************
	x aponta para norte (metros)
	y aponta para leste (metros)
	z aponta para baixo (metros)

	Descrição das variáveis
	***********************
	N = Número de dados observados.
	M = Número total de fontes equivalentes.
	Nsx = Número total de fontes equivalentes na direção x.
	Nsy = Número total de fontes equivalentes na direção y.
	Npolx = Número de polinômios na direção x.
	Npoly = Número de polinômios na direção y.
	Nspolx = Número de fontes equivalentes, na direção x, em cada polinômio.
	Nspoly = Número de fontes equivalentes, na direção y, em cada polinômio.
	Nspol = Número de fontes equivalentes em cada polinômio.
	dx = Distância entre as fontes equivalentes na direção x.
	dy = Distância entre as fontes equivalentes na direção y.
	tobs = Vetor que armazena t observado.
	xp = Vetor que armazena as coordenadas x do t observado.
	yp = Vetor que armazena as coordenadas y do t observado.
	zp = Vetor que armazena as coordenadas z do t observado.
	xmin = Coordenada mais a sul do t.
	xmax = Coordenada mais a norte do t.
	ymin = Coordenada mais a oeste do t.
	ymax = Coordenada mais a leste do t.
	X1 = Coordenada mais a sul da camada equivalente.
	X2 = Coordenada mais a norte da camada equivalente.
	Y1 = Coordenada mais a oeste da camada equivalente.
	Y2 = Coordenada mais a leste da camada equivalente.
	lambidax = extensão da camada equivalente na direção x para além dos limites dos dados observados.
	lambiday = extensão da camada equivalente na direção y para além dos limites dos dados observados.
	t = Vetor que armazena t predito pela camada equivalente. 
	X = vetor que armazena as coordenadas X das fontes equivalentes.
	Y = vetor que armazena as coordenadas Y das fontes equivalentes.
	Z_layer = coordenada Z_layer das fontes equivalentes.
	c = vetor que armazena os coeficientes dos polinômios que descrevem as densidades das fontes equivalentes.
	grau = Variável que armazena o grau dos polinomios que descrevem as densidades das fontes equivalentes.
	Q = Número de coeficientes de cada um dos polinomios que descrevem as densidades das fontes equivalentes.
	H = Número total de coeficientes dos polinomios que descrevem as densidades das fontes equivalentes.
	uridge = parâmetro de regularização Ridge Regression (Tikhonov de ordem zero) para o cálculo da distribuição
		de densidades na camada equivalente.
	usuavidade = parâmetro de regularização Suavidade (Tikhonov de ordem um) para o cálculo da distribuição
		de densidades na camada equivalente.

	Descrição breve do programa
	***************************

	Seja um vetor N-dimensional de dados observados (gravidade) tobs, é possível calcular uma
	camada equivalente que ajusta esses dados. Essa camada é formada por monopolos com densidades
	diferentes. O cálculo da distribuição de densidades nos monopolos é feita resolvendo-se o
	seguinte sistema linear:

	Ap = tobs,												(1)

	sendo p um vetor M-dimensional com as densidades dos monopolos e A uma matriz N x M, cujo
	elemento ij representa o efeito gravitacional, na posição xp[i], yp[i] e zp[i], de um monopolo
	com densidade unitária e localizado na posição X[j], Y[j] e Z_layer.
	
	A camada equivalente pode ser dividida regularmente nas direções x e y, de forma que a distribuição
	de densidades dentro de cada setor retangular seja descrita por um polinômio de grau "grau". Essa
	distribuição de densidades polinomial pode ser representada da seguinte forma:

	Bc = p,													(2)

	em que c é um vetor H-dimensional, cujos elementos são os coeficientes dos polinômios que descrevem
	a distribuição de densidades na camada equivalente e B é uma matriz M x H, cujo elemento ij é a
	derivada do conjunto de polinômios, na posição X[i], Y[i] e Z_layer, em relação ao j-ésimo coeficiente.

	Como a camada é dividida em setores retangulares e cada setor tem uma distribuição de densidades
	descrita por um polinômio, a derivada em relação a um coeficiente j de um polinômio que pertence
	a um setor k, calculada em uma posição X[j], Y[j] e Z_layer pertencente qualquer setor l diferente de k
	possui valor nulo. Essa derivada não é nula apenas quando calculada em uma posição X[j], Y[j] e Z_layer
	pentencento ao mesmo setor k do polinômio.

	Substituindo 2 em 1 chega-se ao seguinte sistema linear

	ABc = tobs												(3)

	e o problema de estimar a densidade em cada monopolo da camada equivalente passa a ser estimar
	os valores dos coeficientes que descrevem a distribuição de densidades em cada setor da camada
	equivalente. Esse procedimento pode diminuir drásticamente o esforço computacional envolvido no
	cálculo da camada equivalente, visto que é possível descrever uma grande quantidade de monopolos
	com um polinômio de grau baixo.

	Para que seja possível, por exemplo, fazer continuação para cima ou calcular o tensor de gravidade
	por meio da camada equivalente, é necessário que a distribuição de densidades estimada seja suave.
	Sendo assim, essa informação foi imposta por meio do regularizador de Tikhonov de ordem um

	Rp = 0,													(4)

	sendo R uma matriz de diferenças finitas, e o problema de estimar a distribuição de densidades na
	camada equivalente passa a ser estimar o vetor c que minimiza a expressão

	(ABc - tobs)t(ABc - tobs) + u(RBc)t(RBc)				(5)

	sendo u	o parâmetro de regularização. A solução da expressão 5 é dada em termos da solução do
	seguinte sistema linear

	(BtAtAB + uBtRtRB)c = BtAttobs.						(6)

	*/

	relatorio = fopen("relatorio.txt", "w");

	leitura_parametros (relatorio, &N, &uridge, &usuavidade, &Z_layer, &grau, &stdev_dados, &semente, &declination_earth, &inclination_earth, &declination_body, &inclination_body, &Npolx, &Npoly, &Nspolx, &Nspoly, &lambidax, &lambiday, &uBtAtAB);

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

	xp = aloca_vetor_double(relatorio, N);
	yp = aloca_vetor_double(relatorio, N);
	zp = aloca_vetor_double(relatorio, N);
	tobs = aloca_vetor_double(relatorio, N);

	leitura_dados(relatorio, N, xp, yp, zp, tobs, &xmin, &xmax, &ymin, &ymax, stdev_dados);

	/* Cálculo dos limites da camada equivalente ==> */
	X1 = xmin - (lambidax*(xmax - xmin));
	X2 = xmax + (lambidax*(xmax - xmin));
	Y1 = ymin - (lambiday*(ymax - ymin));
	Y2 = ymax + (lambiday*(ymax - ymin));
	/* <== Cálculo dos limites da camada equivalente */

	/* Cálculo das coordenadas X e Y das fontes equivalentes ==> */
	printf ("calculo das coordenadas das fontes\n\n");

	Npol = Npolx*Npoly;

	Nsx = Npolx*Nspolx;
	Nsy = Npoly*Nspoly;

	Nspol = Nspolx*Nspoly;
	
	dx = (double)((X2 - X1)/(Nsx - 1.0));
	dy = (double)((Y2 - Y1)/(Nsy - 1.0));
	
	M = Nsx*Nsy;

	/* Cálculo do número de coeficientes dos polinômios ==> */
	for (Q = 1, i = 2; i <= (grau+1); i++) {
	
		Q += i; 
	
	}
	/* <== Cálculo do número de coeficientes dos polinômios */

	H = Npol*Q;

	X = aloca_vetor_double(relatorio, M);
	Y = aloca_vetor_double(relatorio, M);

	coordenadas_fontes(Npolx, Npoly, Nspolx, Nspoly, dx, dy, X1, Y1, X, Y);
	/* <== Cálculo das coordenadas X e Y das fontes equivalentes */

	/* Cálculo da matriz B ==> */
	time (&start2);
	
	B = aloca_matriz_double (relatorio, M, H);

	calcula_B (Npol, Nspol, grau, Y, X, B, Y1, Y2, X1, X2);
	
	time (&end2);
	
	aux0 = difftime (end2, start2);
	
	fprintf (relatorio, "Cálculo da matriz B (%.5lf segundos)\n\n", aux0);
	
	/* <== Cálculo da matriz B */

	GtG = aloca_matriz_double (relatorio, H, H);
	BtRtRB = aloca_matriz_double (relatorio, H, H);
	BtAtAB = aloca_matriz_double (relatorio, H, H);
	Gttobs = aloca_vetor_double (relatorio, H);
	g = aloca_vetor_double (relatorio, H);

	/* Cálculo da matriz GtG = BtAtAB e do vetor Gttobs ==> */
	if (uBtAtAB == 1.0) {

		time (&start2);
	
		a = aloca_vetor_double (relatorio, M);

		printf ("matriz BtAtAB\n\n");
		calcula_BtAtAB (N, Q, Npol, Nspol, grau, M, H, Y, X, Z_layer, yp, xp, zp, Ly, Lx, Lz, ly, lx, lz, a, B, g, BtAtAB, Gttobs, tobs, &fator_normalizacao);
	
		a = libera_vetor_double (a);
	
		time (&end2);
	
	}
	
	else {
	
		time (&start2);

		entrada = fopen ("BtAtAB.txt", "r");
		
		for (i = 0; i < H; i++) {
		
			for (j = i; j < H; j++) {
			
				fscanf(entrada, "%lf", &BtAtAB[i][j]);
			
			}
			
			fscanf(entrada, "%lf", &Gttobs[i]);
		
		}
		
		fclose (entrada);

		for (aux0 = 0.0, i = 0; i < H; i++) {

			aux0 += BtAtAB[i][i];

		}

		aux0 = (double)(aux0/H);

		fator_normalizacao = aux0;
	
		time (&end2);
	
	}
	
	aux0 = difftime (end2, start2);
	
	fprintf (relatorio, "Cálculo da matriz BtAtAB (%.5lf segundos)\n\n", aux0);

	/* <== Cálculo da matriz GtG = BtAtAB e do vetor Gttobs */

	/* Incorporação do vínculo de suavidade. Cálculo da matriz BtRtRB ==> */
	time (&start2);
	
	printf ("matriz BtRtRB\n\n");

	calcula_BtRtRB (H, Q, Npolx, Npoly, Nspol, Nsy, Nspolx, Nspoly, B, g, BtRtRB, fator_normalizacao, usuavidade);
	
	time (&end2);
	
	aux0 = difftime (end2, start2);
	
	fprintf (relatorio, "Cálculo da matriz BtRtRB (%.5lf segundos)\n\n", aux0);
	
	/* <== Incorporação do vínculo de suavidade. Cálculo da matriz BtRtRB */
	
	/* Calcula GtG = (BtRtRB + BtAtAB) ==> */
	for (i = 0; i < H; i++) {
	
		for (j = i; j < H; j++) {
		
			GtG[i][j] = (BtRtRB[i][j] + BtAtAB[i][j]);
		
		}
	
	}
	/* <== Calcula GtG = (BtRtRB + BtAtAB) */

	/* incorporação do vínculo Ridge */
	incorpora_ridge (H, GtG, uridge, fator_normalizacao);
	
	/* Estimativa dos coeficientes dos polinômios que descrevem as densidades 
	na camada equivalente ==> */
	c = aloca_vetor_double (relatorio, H);

	time (&start2);
	
	printf ("\nResolucao do sistema linear\n");
	/* Resolução do sistema GtGc = Gttobs */
	/*sist_lin_LU (relatorio, H, GtG, Gttobs, c, 0.0);*/
	
	choldc(relatorio, GtG, H, g); /* Decomposição de Cholesky */
	
	cholsl(GtG, H, g, Gttobs, c); /* Resolução do sistema linear */

	time (&end2);

	aux0 = difftime (end2, start2);
	
	fprintf (relatorio, "Resolucao do sistema linear (%.5lf segundos)\n\n", aux0);

	GtG = libera_matriz_double(H, GtG);
	BtRtRB = libera_matriz_double(H, BtRtRB);
	BtAtAB = libera_matriz_double(H, BtAtAB);
	Gttobs = libera_vetor_double(Gttobs);
	g = libera_vetor_double (g);
	/* <== Estimativa dos coeficientes dos polinômios que descrevem as densidades 
	na camada equivalente */

	/* Cálculo das intensidades de magnetização na camada equivalente ==> */
	printf ("\nCalculo das intensidades de magnetizacao na camada\n");
	p = aloca_vetor_double (relatorio, M);
	
	calcula_magnetizacoes(M, H, c, p, B);
	
	B = libera_matriz_double (M, B);
	/* <== Cálculo das densidades na camada equivalente */

	/* Cálculo de t e das componentes do tensor ==> */
	printf ("\nCalculo de t\n\n");
	
	t = aloca_vetor_double (relatorio, N);

	/*t_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, t);*/
	mag_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, Ly, Lx, Lz, ly, lx, lz, p, t);

	/* <== Cálculo de t e das componentes do tensor */

	/* impressão dos arquivos de saida ==> */

	saida = fopen ("ajuste.txt", "w");
	impressao_componente (saida, N, yp, xp, zp, tobs, t);
	fclose (saida);

	saida = fopen ("camada.txt", "w");
	impressao_camada (saida, M, Y, X, p);
	fclose (saida);

	/*tz_equivalente(N, M, yp, xp, zp, Y, X, Z_layer, p, t);

	saida = fopen ("tz_calculado.txt", "w");
	impressao_componente (saida, N, yp, xp, zp, t);
	fclose (saida);*/
	
	/* <== impressão dos arquivos de saida */

	time (&end1);
	
	aux0 = difftime (end1, start1);

	printf("\nPrograma finalizado com sucesso em %.3lf segundos!\n\n", aux0);
	
	fprintf(relatorio, "\nPrograma finalizado com sucesso em %.3lf segundos!\n\n", aux0);
	
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

        system ("PAUSE");
        
        return (NULL);

    }

    /* aloca as colunas da matriz */

    for (i = 0; i < linha; i++ ) {

        m[i] = (double *)calloc(coluna, sizeof(double));

        if (m[i] == NULL) {

            fprintf (arq, "Memoria Insuficiente (colunas)!\n\n");

            fclose (arq);

            system ("PAUSE");
            
            return (NULL);

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
        
        system ("PAUSE");
        
		return (NULL);
        
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

double sist_lin_LU (FILE *relatorio, int M, double **GtG, double *Gtd, double *p, double uridge) {

	int i, j, k, l;
	double *b, **L;
	double aux0, aux1;

	/************* decomposição LU da matriz (GtG + uridgeI) ==> ******************/
	printf ("   Decomposicao LU da matriz GtG\n");

	L = aloca_matriz_double (relatorio, M, M);

	/**** primeira linha de U, primeira coluna de L e diagonal principal de L ****/

	aux0 = 0.01*GtG[0][0];
	aux1 = uridge*aux0;

	L[0][0] = GtG[0][0] + aux1;

	for (i = 1; i < M; i++) {

		L[0][i] = GtG[0][i];
		L[i][0] = (double)(GtG[i][0]/L[0][0]);

	}

	/**** varre as linhas de U e as colunas de L ****/
	for (i = 1; i < M; i++) {

		/**** i-ésima linha de U ****/
		j = i;

			for (k = 0; k < i; k++) {

				L[i][j] -= L[i][k]*L[k][j];

			}

			aux0 = 0.01*GtG[i][j];
			aux1 = uridge*aux0;

			L[i][j] += GtG[i][j] + aux1;

		for (j = i+1; j < M; j++) {

			for (k = 0; k < i; k++) {

				L[i][j] -= L[i][k]*L[k][j];

			}

			L[i][j] += GtG[i][j];

		}

		/**** i-ésima coluna de L ****/
		for (j = i+1; j < M; j++) {

			for (k = 0; k < i; k++) {

				L[j][i] -= L[j][k]*L[k][i];

			}

			L[j][i] += GtG[j][i];

			L[j][i] = (double)(L[j][i]/L[i][i]);

		}

	}

	/*********** <== decomposição LU da matriz (GtG + uridgeI) **************/

	/******* resolução do sistema linear (GtG uridgeI)p = Gtd ==> ***********/

	b = aloca_vetor_double(relatorio, M);

	/********** (GtG uridgeI)p = Gtd, LUp = Gtd, Lb = Gtd, Up = b ***********/

	/**** resolução do sistema Lb = Gtd ***/
	printf ("   Resolucao do sistema Lb = Gtd\n");
	
	for (k = 0; k < M; k++) {

		b[k] = 0.0;

		for (l = 0; l < k; l++) {

			b[k] -= L[k][l]*b[l];

		}

		b[k] += Gtd[k];

	}

	/*** resolução do sistema Up = b ***/
	printf ("   Resolucao do sistema Up = b\n\n");
	
	for (j = M-1; j >= 0; j--) {

		p[j] = 0.0;

		for (k = j+1; k < M; k++) {

			p[j] -= L[j][k]*p[k];

		}

		p[j] += b[j];

		p[j] = (double)(p[j]/L[j][j]);

	}

	/******* <== resolução do sistema linear (GtG uridgeI)p = Gtd ***********/

	b = libera_vetor_double(b);
	L = libera_matriz_double (M, L);

	return 0 ;

}

/* Função gasdev do capítulo sete (Random Numbers) do Numerical Recipes */

/* Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
as the source of uniform deviates. */

float gasdev(long *idum) {

	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	/* We don't have an extra deviate handy, so pick two uniform numbers in the
	square extending from -1 to +1 in each direction, see if they are in the unit
	circle, and if they are not, try again. */

	if (iset == 0) {

		do {

			v1=2.0*ran1(idum)-1.0;

			v2=2.0*ran1(idum)-1.0;

			rsq=v1*v1+v2*v2;

		} while (rsq >= 1.0 || rsq == 0.0);

		fac=sqrt(-2.0*log(rsq)/rsq);

		/* Now make the Box-Muller transformation to get two normal deviates. Return one and
		save the other for next time. */

		gset=v1*fac;

		iset=1; /* Set flag. */

		return v2*fac;

	}

	/* We have an extra deviate handy, so unset the flag, and return it. */

	else {

		iset=0;

		return gset;

	}

}

/* Função ran1 do capítulo sete (Random Numbers) do Numerical Recipes */

/* "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest floating value that is
less than 1. */

float ran1(long *idum) {

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) { /* Initialize. */

		if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0. */

		else *idum = -(*idum);

		for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups). */

			k=(*idum)/IQ;

			*idum=IA*(*idum-k*IQ)-IR*k;

			if (*idum < 0) *idum += IM;

			if (j < NTAB) iv[j] = *idum;

		}

		iy=iv[0];

	}

	/* Start here when not initializing. */
	k=(*idum)/IQ;

	/* Compute idum=(IA*idum) % IM without overflows by Schrage's method. */
	*idum=IA*(*idum-k*IQ)-IR*k;

	if (*idum < 0) *idum += IM;

	/* Will be in the range 0..NTAB-1. */
	j=iy/NDIV;

	/* Output previously stored value and refill the shuffe table. */
	iy=iv[j];

	iv[j] = *idum;

	/* Because users don't expect endpoint values. */
	if ((temp=AM*iy) > RNMX) return RNMX;

	else return temp;

}

void leitura_parametros (FILE *relatorio, int *N, double *uridge, double *usuavidade, double *Z_layer, int *grau, double *stdev_dados, int *semente, double *declination_earth, double *inclination_earth, double *declination_body, double *inclination_body, int *Npolx, int *Npoly, int *Nspolx, int *Nspoly, double *lambidax, double *lambiday, double *uBtAtAB) {

	char str[20];
	
	FILE *entrada;

	sprintf (str, "parametros.txt");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	entrada = fopen(str, "r");

	fscanf (entrada, "%d %lf %lf %lf %d", N, uridge, usuavidade, Z_layer, grau);

	fscanf(entrada, "%lf %d", stdev_dados, semente);
	
	fscanf(entrada, "%lf %lf %lf %lf", declination_earth, inclination_earth, declination_body, inclination_body);
	
	fscanf(entrada, "%d %d %d %d", Npolx, Npoly, Nspolx, Nspoly);
	
	fscanf(entrada, "%lf %lf", lambidax, lambiday);
	
	fscanf(entrada, "%lf", uBtAtAB);
	
	idum = semente;

	fclose (entrada);

}

void leitura_dados(FILE *relatorio, int N, double *Xp, double *Yp, double *Zp, double *tobs, double *xmin, double *xmax, double *ymin, double *ymax, double stdev_dados) {

	int i;
	double aux;
	char str[20];
	
	FILE *entrada;

	sprintf (str, "tobs.txt");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	entrada = fopen(str, "r");

	if (fscanf(entrada, "%lf %lf %lf %lf", &Yp[0], &Xp[0], &Zp[0], &tobs[0]) != 4) {

		fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	tobs[0] += stdev_dados*gasdev(&idum);

	(*xmin) = Xp[0];
	(*xmax) = Xp[0];
	(*ymin) = Yp[0];
	(*ymax) = Yp[0];

	for (i = 1; i < N; i++) {

		if (fscanf(entrada, "%lf %lf %lf %lf", &Yp[i], &Xp[i], &Zp[i], &tobs[i]) != 4) {

			fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", str);

			fclose (relatorio);

			printf ("Erro!\n\n");

			system ("PAUSE");

			return 0;

		}

		tobs[i] += stdev_dados*gasdev(&idum);

		/* Determinação dos limites máximo e mínimo em X e Y */
		if ((*xmin) > Xp[i]) {

			(*xmin) = Xp[i];

		}
		if ((*xmax) < Xp[i]) {

			(*xmax) = Xp[i];

		}
		if ((*ymin) > Yp[i]) {

			(*ymin) = Yp[i];

		}
		if ((*ymax) < Yp[i]) {

			(*ymax) = Yp[i];

		}

	}
	
	fclose (entrada);

}

void linha_a (int i, int M, double *Y, double *X, double Z_layer, double *yp, double *xp, double *zp, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *a) { 

	int j;
	double aux1, aux2, aux3, aux4, aux5;

	for (j = 0; j < M; j++) {

		/*a[j] = t_monopolo(i, yp, xp, zp, j, Y, X, Z_layer);*/
		a[j] = dipolo_mag2(i, yp, xp, zp, j, Y, X, Z_layer, Ly, Lx, Lz, ly, lx, lz);
	
	}

}

void calcula_BtAtAB (int N, int Q, int Npol, int Nspol, int grau, int M, int H, double *Y, double *X, double Z_layer, double *yp, double *xp, double *zp, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *a, double **B, double *g, double **GtG, double *Gttobs, double *tobs, double *fator_normalizacao) {

	double aux0;
	int i, j, k, k1, k2, l, m, n;
	
	FILE *arquivo;

	/* percorre as N linhas da matriz A */
	for (i = 0; i < N; i++) { /* for i */

		/* Cálculo da i-ésima linha da matriz A */
		linha_a (i, M, Y, X, Z_layer, yp, xp, zp, Ly, Lx, Lz, ly, lx, lz, a);

		for (l = 0, j = 0; j < Npol; j++) { /* for j */
		
			k1 = (j*Nspol);
			k2 = k1 + Nspol;

			for (k = 0; k < Q; k++, l++) { /* for k */
			
				/* Produto escalar entre a i-ésima linha de A
				e a l-ésima coluna de B ==> */
				g[l] = a[k1]*B[k1][l];

				for (m = (k1+1); m < k2; m++) {
						
					g[l] += a[m]*B[m][l];
						
				}					
				/* <== Produto escalar entre a i-ésima linha de A
				e a l-ésima coluna de B */

				/* Cômputo de uma parcela do l-ésimo
				elemento do vetor Gttobs */
				Gttobs[l] += g[l]*tobs[i];

				for (n = 0; n <= l; n++) {

					/* Cômputo de uma parcela do
					elemento nl da matriz GtG */
					/* GtG[l][n] += g[l]*g[n]; preenche a triangular inferior de GtG */
					GtG[n][l] += g[l]*g[n]; /* preenche a triangular superior de GtG */

				}
						
				printf ("i = %5d N = %5d l = %5d H = %5d\n", i, N, l, H);
				
			} /* for k */

		} /* for j */
		
	} /* for i */

	arquivo = fopen ("BtAtAB.txt", "w");
	
	for (i = 0; i < H; i++) {
	
		for (j = i; j < H; j++) {
	
			fprintf (arquivo, "%.15lf ", GtG[i][j]);
	
		}
		
		fprintf (arquivo, "%.15lf\n", Gttobs[i]);
	
	}
	
	fclose (arquivo);
	
	for (aux0 = 0.0, i = 0; i < H; i++) {

		aux0 += GtG[i][i];

	}

	aux0 = (double)(aux0/H);
	
	(*fator_normalizacao) = aux0;

}

void calcula_BtRtRB (int H, int Q, int Npolx, int Npoly, int Nspol, int Nsy, int Nspolx, int Nspoly, double **B, double *g, double **GtG, double fator_normalizacao, double usuavidade) {

	int i, j, k, l, m, n, o, p, q;
	int aux1, aux2;
	double aux0;

	aux1 = Nsy*Nspolx;
	aux2 = ((Nspolx - 1)*Nspoly) + 1;
	
	/* Suavidade na direção y ==> */
	for (i = 0; i < Npolx; i++) { /* for i */

		/*l = Nspoly - 1 + (i*Nsy*Nspolx);*/
		l = Nspoly - 1 + (i*aux1);

		for (j = 0; j < (Npoly - 1); j++) { /* for j */

			/* a impressão abaixo é para avaliar a velocidade do programa, de maneira visual */
			printf ("i = %5d Npolx = %5d j = %5d (Npoly - 1) = %5d\n", (i+1), Npolx, (j+1), (Npoly - 1));

			for (k = 0; k < Nspolx; k++, l += Nspoly) { /* for k */

				/*for (m = 0; m < H; m++) { * for m * */
				o = (l/Nspol);
				o *= Q;
				p = o + (2*Q);
				
				for (m = o; m < p; m++) {

					/* Produto escalar entre a l-ésima linha de R e a m-ésima coluna de B */
					/*g[m] = B[l][m] - B[(l + ((Nspolx - 1)*Nspoly) + 1)][m];*/
					g[m] = B[l][m] - B[(l + aux2)][m];

					/*for (n = 0; n <= m; n++) { * for n * */
					for (n = o; n <= m; n++) {

						/* Cômputo de uma parcela do
						elemento nm da matriz GtG */
						/*GtG[m][n] += g[m]*g[n]; preenche a triangular inferior de GtG */
						GtG[n][m] += g[m]*g[n]; /* preenche a triangular superior de GtG */

					} /* for n */

				} /* for m */
				
			} /* for k */

		} /* for j */

	} /* for i */
	/* <== Suavidade na direção y */
	
	aux1 = Nspoly*(Nspolx - 1);
	aux2 = ((Npoly - 1)*Nspol) + Nspoly;
	
	q = (aux2/Q);
	q *= Q;
	q ++;

	/* Suavidade na direção x ==> */
	for (l = 0, i = 0; i < Npoly; i++) { /* for i */

		for (j = 0; j < (Npolx - 1); j++) { /* for j */

		/* a impressão abaixo é para avaliar a velocidade do programa, de maneira visual */
		printf ("i = %5d Npoly = %5d j = %5d (Npolx - 1) = %5d\n", (i+1), Npoly, (j+1), (Npolx - 1));

			/*l += Nspoly*(Nspolx - 1);*/
			l += aux1;

			for (k = 0; k < Nspoly; k++, l++) { /* for k */

				/*for (m = 0; m < H; m++) { * for m * */

				o = (l/Nspol);
				o *= Q;
				p = o + q;
				
				if (p > H) {
				
					p = H;
					
				}
				
				for (m = o; m < p; m++) {

					/* Produto escalar entre a l-ésima linha de R e a m-ésima coluna de B */
					/*g[m] = B[l][m] - B[(l + ((Npoly - 1)*Nspol) + Nspoly)][m];*/
					g[m] = B[l][m] - B[(l + aux2)][m];

					/*for (n = 0; n <= m; n++) { * for n * */
					for (n = o; n <= m; n++) {

						/* Cômputo de uma parcela do
						elemento nm da matriz GtG */
						/* GtG[m][n] += g[m]*g[n]; preenche a triangular inferior de GtG */
						GtG[n][m] += g[m]*g[n]; /* preenche a triangular superior de GtG */

					} /* for n */

				} /* for m */

			} /* for k */

		} /* for j */

	} /* for i */
	/* <== Suavidade na direção x */

	for (aux0 = 0.0, i = 0; i < H; i++) {

		aux0 += GtG[i][i];

	}

	aux0 = (double)(aux0/H);
	aux0 = (double)(fator_normalizacao/aux0);

	usuavidade *= aux0;

	/* Multiplica a parte de GtG acima da diagonal principal
	pelo parâmetro de regularização usuavidade */
	for (i = 0; i < H; i++) {

		for (j = i; j < H; j++) {

			GtG[i][j] *= usuavidade;

		}

	}

}

void coordenadas_fontes(int Npolx, int Npoly, int Nspolx, int Nspoly, double dx, double dy, double X1, double Y1, double *X, double *Y) {

	int i, j, k, l, m;
	double x, y, x0, y0;

	for (m = 0, i = 0; i < Npolx; i++) {
	
		x0 = X1 + (i*Nspolx*dx);
		
		for (j = 0; j < Npoly; j++) {
		
			y0 = Y1 + (j*Nspoly*dy);
			
			for (k = 0, x = x0; k < Nspolx; k++, x += dx) {
			
				for (l = 0, y = y0; l < Nspoly; l++, y += dy, m++) {
			
					Y[m] = y;
					X[m] = x;
			
				}
			
			}
		
		}
	
	}

}

void impressao_componente(FILE *arquivo, int N, double *y, double *x, double *z, double *tobs, double *tcalc) {

	int i;

	fprintf (arquivo , "              y              x              z           tobs          tcalc       diferenca\n");

	for (i = 0; i < N; i++) {

		fprintf (arquivo, "%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3E\n", y[i], x[i], z[i], tobs[i], tcalc[i], (tcalc[i] - tobs[i]));

	}

}

void calcula_magnetizacoes(int M, int H, double *c, double *p, double **B) {

	int i, j;

	for (i = 0; i < M; i++) {
	
		p[i] = 0.0;
		
		for (j = 0; j < H; j++) {
		
			p[i] += B[i][j]*c[j];
		
		}
	
	
	}

}

void impressao_camada(FILE *arquivo, int M, double *Y, double *X, double *p) {

	int i;

	fprintf (arquivo , "              Y               X               p\n");

	for (i = 0; i < M; i++) {

		fprintf (arquivo, "%15.3lf %15.3lf %15.5E\n", Y[i], X[i], p[i]);

	}

}

void coluna_b (int M, int j, int Nspol, int l, int k, double *X, double *Y, double *b) {

	int i, i_inicial, i_final;

	i_inicial = (j*Nspol);
	i_final = i_inicial + Nspol;

	for (i = 0; i < M; i++) {

		b[i] = 0.0;

	}

	for (i = i_inicial; i < i_final; i++) {	

		b[i] = pow(Y[i], l)*pow(X[i], (k-l));
	
	}
	
}


void calcula_B (int Npol, int Nspol, int grau, double *Y, double *X, double **B, double Y1, double Y2, double X1, double X2) {

	int i, i_inicial, i_final;
	int j, k, l, m;
	double Y_medio, X_medio, constante_ponderacao;
	
	Y_medio = (Y2 - Y1)*0.5;
	X_medio = (X2 - X1)*0.5;

	constante_ponderacao = (Y_medio + X_medio)*0.5;
	constante_ponderacao = (double)(1.0/constante_ponderacao);

	/* os próximos 3 laços percorrem as H colunas da matriz B */
	for (m = 0, j = 0; j < Npol; j++) { /* for j */
			
		i_inicial = (j*Nspol);
		i_final = i_inicial + Nspol;

		for (k = 0; k <= grau; k++) { /* for k */
			
			for (l = 0; l <= k; l++, m++) { /* for l */

				/* Cálculo da m-ésima coluna da matriz B ==> */
				for (i = i_inicial; i < i_final; i++) {

					/*B[i][m] = constante_ponderacao*pow(Y[i], l)*pow(X[i], (k-l));*/
					B[i][m] = pow(Y[i], l)*pow(X[i], (k-l));
					
				}
				/* <== Cálculo da m-ésima coluna da matriz B */

			} /* for l */

		} /* for k */

	} /* for j */

}

void choldc(FILE *relatorio, double **a, int n, double *p) {
/*Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
decomposition, A = LLt . On input, only the upper triangle of a need be given; it is not
modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
elements which are returned in p[1..n].*/

	int i,j,k;
	double sum;
	
	for (i = 0; i < n; i++) {
	
		j = i;
		
			for (sum = a[i][j], k = i-1; k >= 0; k--) {
			
				sum -= a[i][k]*a[j][k];
				
			}

			if (sum <= 0.0) { /* a, with rounding errors, is not positive definite */ 

				fprintf (relatorio, "\nCholesky falhou!!\n\n");
				
				fclose (relatorio);
				
				abort();
					
			}
			
			p[i] = sqrt(sum);

		for (j = (i+1); j < n; j++) {
		
			for (sum = a[i][j], k = i-1; k >= 0; k--) {
			
				sum -= a[i][k]*a[j][k];
				
			}

			a[j][i] = (double)(sum/p[i]);
		
		}

	}
	
}

void incorpora_ridge (int H, double **GtG, double uridge, double fator_normalizacao) {

	int i;
	double aux0, aux1;
	
	if (uridge != 0.0) {
	
		aux1 = uridge*fator_normalizacao;

		for (i = 0; i < H; i++) {

			GtG[i][i] += aux1;
	
		}

	}

}

void cholsl(double **a, int n, double *p, double *b, double *x) {

/*Solves the set of n linear equations Ax = b, where a is a positive-definite symmetric matrix.
a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. Only the lower
triangle of a is accessed. b[1..n] is input as the right-hand side vector. The solution vector is
returned in x[1..n]. a, n, and p are not modified and can be left in place for successive calls
with diferent right-hand sides b. b is not modified unless you identify b and x in the calling
sequence, which is allowed. */

	int i,k;
	double sum;

	for (i = 0; i < n; i++) { /* Solve Ly = b, storing y in x. */

		for (sum = b[i], k = i-1; k >= 0; k--) {
		
			sum -= a[i][k]*x[k];
			
		}
		
		x[i] = (double)(sum/p[i]);

	}

	for (i = (n - 1); i >= 0; i--) { /* Solve Ltx = y. */

		for (sum = x[i], k = i+1; k < n; k++) {
		
			sum -= a[k][i]*x[k];
			
		}
		
		x[i] = (double)(sum/p[i]);

	}
	
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

void mag_equivalente(int N, int M, double *yp, double *xp, double *zp, double *Y, double *X, double Z_layer, double Ly, double Lx, double Lz, double ly, double lx, double lz, double *p, double *t) {

	int i, j;

	for (i = 0; i < N; i++) {
	
		t[i] = 0.0;
	
		for (j = 0; j < M; j++) {
		
			t[i] += dipolo_mag2(i, yp, xp, zp, j, Y, X, Z_layer, Ly, Lx, Lz, ly, lx, lz)*p[j];
		
		}
	
	}

}

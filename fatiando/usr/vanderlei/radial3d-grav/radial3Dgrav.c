/* *****************************************************************************
 Copyright 2010 The Fatiando a Terra Development Team

 This file is part of Fatiando a Terra.

 Fatiando a Terra is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Fatiando a Terra is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Fatiando a Terra.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************** */

/*
Author: Vanderlei C. Oliveira Jr.
Date: 24 Jun 2011
Contact: vandscoelho@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Divide_macro(a, b) ((a)/((b) + (1E-15)))
#define Divide_macro2(a, b) (((a) + (1E-15))/((b) + (1E-15)))
#define PI 3.14159265358979323846

double GRAV = 6.67E-20; /* (kilômetro cúbico)/((kilograma)*(segundo ao quadrado)) */

double **aloca_matriz_double (FILE *, int, int);

double **libera_matriz_double (int, double **);

double *aloca_vetor_double (FILE *, int);

double *libera_vetor_double (double *);

int *aloca_vetor_int (FILE *, int);

int *libera_vetor_int (int *);

double mod_direto (int N, int M, int P, int nvertices, double *Xp, double *Yp, double *Zp, double *x0, double *y0, double *rho, double *z1, double *z2, double **raio, double teta, double **x, double **y, double *gcalc, double *gobs);

double sist_lin_LU (FILE *, int, int, int, double **, double *, double **, double);

double derivada_ajuste_raio (int ponto, int poligono, int vertice, double *Xp, double *Yp, double *Zp, double *x0, double *y0, double *rho, double *z1, double *z2, int nvertices, double **raio, double teta, double **x, double **y, double dr);

double derivada_ajuste(FILE *relatorio, int N, int M, int nvertices, int P, double **H_ajuste, double *grad_ajuste, double *gobs, double *gcalc, double **raio, double dr, double teta, double *rho, double *z1, double *z2, double **x, double **y, double *x0, double *y0, double *Xp, double *Yp, double *Zp);

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

int paraview (int, int, int, double, double **, double *, double *, double *, double *, int);

int paraview_error (int, int, double, double **, double **, double *, double *, double *, double *, double *, int);

double derivada_ajuste_x0 (int ponto, int prisma, int nvertices, double *Xp, double *Yp, double *Zp, double *x0, double *y0, double *rho, double *z1, double *z2, double **raio, double teta, double **x, double **y, double dr);

double derivada_ajuste_y0 (int ponto, int prisma, int nvertices, double *Xp, double *Yp, double *Zp, double *x0, double *y0, double *rho, double *z1, double *z2, double **raio, double teta, double **x, double **y, double dr);

/* Contribuição dos vínculos na função objetivo ==> */

double fobj_ridge (double uridge, double **raio, int M, int nvertices);

double fobj_absolut_equality_raios (double uae_raios, double **raio, double *pae, int nvertices);

double fobj_absolut_equality_origem (double uae_x0, double uae_y0, double x0_afloramento, double y0_afloramento, double *x0, double *y0, int M);

double fobj_flatness_raios (double uflatness_rad, double uflatness_vert, double **raio, int M, int nvertices);

double	fobj_flatness_origens (double uflatness_x0, double uflatness_y0, double *x0, double *y0, int M);

/* <== Contribuição dos vínculos na função objetivo */

int main() {

	int i, j, k, l, m, n;
	int N, M, Q, P;

	double **x, **y, teta, **raio, **raio_inicial, **raio_trans, **raio_medio, **raio_std, **dp, *x0_inicial, *y0_inicial, *x0, *y0, *x0_trans, *y0_trans, *x0_medio, *y0_medio, *x0_std, *y0_std, *rho, *z1, *z2;
	int nvertices;
    double *Xp, *Yp, *Zp, *gobs, *gobs_original, *gcalc;

	int ITMAX, ITMAX_marq, iteracao, iteracao_marq, ITMAX_repetibilidade;
	double lambida_inicial, lambida, dlambida, fobj1, fobj2, variacao_relativa, epsilon, dr, dz, z0;
	double *grad, *grad_ajuste, **H, **H_ajuste;
	double uridge, rmax, rmin, x0max, x0min, y0max, y0min, uflatness_rad, uflatness_vert;
	double uflatness_x0, uflatness_y0;

	int semente;
	double stdev_dados, stdev_perturbacao, r0;
	double aux0, aux1, aux2, aux3;

	int L;
	int *npontos;
	double uae_raios, uae_x0, uae_y0, x0_afloramento, y0_afloramento;
	double **contorno, *pae;

	int contador_repetibilidade;

	double massa_pre, massa_pre_media = 0.0, massa_pre_std = 0.0;
	double ajuste, ajuste_medio = 0.0, ajuste_std = 0.0, s, s_medio = 0.0, s_std = 0.0;
	
	time_t start, end;

    /*

    N: número de pontos onde se deseja calcular a atração gravitacional.
    M: número de polígonos.
    Q: número total de vértices (raios).
	P: Q mais as coordenadas x0 e y0 de cada placa, exceto da primeira (P = Q + 2M).
    **x: Matriz que armazena as coordenadas x dos vértices dos polígonos. A posição
        ij armazena a coordenada x do j-ésimo vértice do i-ésimo polígono. Esta coordenada
        está na direção norte-sul.
    **y: Matriz que armazena as coordenadas y dos vértices dos polígonos. A posição
        ij armazena a coordenada y do j-ésimo vértice do i-ésimo polígono. Esta coordenada
        está na direção leste-oeste.
    teta: Ângulo teta entre os raios dos polígono.
    **raio_inicial: Matriz que armazena a estimativa inicial das coordenadas dos vértices
		dos polígonos. A posição ij armazena a coordenada raio do j-ésimo vértice do i-ésimo
		polígono.
	r0: armazena o valor númérico da estimativa inicial das coordenadas dos vétices dos polígonos.
    **raio: Matriz que armazena as coordenadas raio dos vértices dos polígonos. A posição
        ij armazena a coordenada raio do j-ésimo vértice do i-ésimo polígono.
    *x0_inicial: Vetor que armazena a coordenada inicial x da origem de cada polígono.
    *y0_inicial: Vetor que armazena a coordenada inicial y da origem de cada polígono.
    *x0: Vetor que armazena a coordenada x da origem de cada polígono.
    *y0: Vetor que armazena a coordenada y da origem de cada polígono.
    *rho: Vetor que armazena o contraste de densidade de cada polígono.
    nvertices: Número de vértices dos polígonos.
    *z1, *z2: Vetores que armazenam a coordenada z do topo e da base de cada polígono,
        respectivamente.
    *Xp, *Yp, *Zp: Vetores que armazenam as coordenadas x, y e z, respectivamente, dos
        pontos onde se deseja calcular a atração gravitacional. Xp e Yp estão
        nas direções norte-sul e leste-oeste, respectivamente.
    *gobs: Vetor que armazena os dados observados.
    *gpre: Vetor que armazena os dados preditos pelo modelo.
	**H_ajuste: matriz Hessiana do funcional ajuste.
	*grad_ajuste: vetor gradiente do funçional ajuste.
	**H: Matriz Hessiana da funçãoobjetivo.
	*grad: Vetor gradiente da função objetivo.
	dr: perturbação nos parâmetros para o cálculo da matriz Jacobiana do modelo direto.
	lambida: parâmetros de Marquardt.
	dlambida: constante que multiplica ou divide o parâmetro de Marquardt.
	ITMAX: número máximo de iterações do looping externo.
	iteracao: contador do looping externo.
	ITMAX_marq: número máximo de iterações do algoritmo de Marquardt.
	iteracao_marq: contador do algoritmo de Marquardt.
	variacao_relativa: é uma razão que mede a variação da função objetivo na iteração k+1
	    em relação a k.
	epsilon: valor que deterina o critério de parada.
	uridge: parâmetro de regularização Ridge Regression.
	stdev_dados: estimativa do erro dos dados.
	semente: número inteiro negativo que inicializa o gerador de números aleatórios.
	rmax: limite superior para os parâmetros.
	rmin: limite inferior para os parâmetros.
	uflatness_rad: parâmetro de regularização Flatness entre os raios de um mesmo prisma.
	uflatness_vert: parâmetro de regularização Flatness entre os raios de prismas adjacentes.
	utv_x0: parâmetro de regularização Total Variation entre as coordenadas x0 de prismas adjacentes.
	utv_y0: parâmetro de regularização Total Variation entre as coordenadas y0 de prismas adjacentes.
	beta: constante real positiva utilizada no funcional Total Variation.
	ajuste: valor do funcional ajuste.
	dz: espessura dos prismas.
	z0: topo do primeiro polígono.
	beta: parâmetro que pondera o Ridge.
   	uae_raios: parâmetro de regularização Absolut Equality nos raios da primeira placa.
   	uae_x0: parâmetro de regularização Absolut Equality do x0 da primeira placa.
   	uae_y0: parâmetro de regularização Absolut Equality do y0 da primeira placa.
   	contorno: matriz que contem as coordenadas em Y, X e raios, tetas dos pontos que descrevem o
   	    contorno do topo da parte aflorante do corpo.
   	pae: vetor que contem os raios que descrevem o contorno da parte aflorante do corpo.
	npontos: vetor cujo i-ésimo elemento contem o número de pontos do contorno da parte aflorante do corpo
	    que são utilizados para computar o valor do i-ésimo elemento de pae.
   	
    */

	char str[100];

	FILE *relatorio, *entrada, *saida;

	time (&start);

	relatorio = fopen ("relatorio_inicial.txt", "w");

	/******** leitura do arquivo da estimativa inicial ==> **************/

	sprintf (str, "input.txt");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	entrada = fopen(str, "r");

	fscanf(entrada, "%lf %lf %lf %lf", &lambida_inicial, &dlambida, &dr, &epsilon);
	fscanf(entrada, "%d %d %d", &ITMAX, &ITMAX_marq, &ITMAX_repetibilidade);
	fscanf(entrada, "%lf %lf %d", &stdev_dados, &stdev_perturbacao, &semente);

	idum = semente;
	
	/* "Esquenta" gerador de números aleatórios ==> */
	for (i = 0; i < 1000; i++) {
	
		aux0 = gasdev(&idum);
	
	}
	/* <== "Esquenta" gerador de números aleatórios */

	fscanf (entrada, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &uridge, &uflatness_rad, &uflatness_vert, &rmin, &rmax, &y0min, &y0max, &x0min, &x0max, &uflatness_y0, &uflatness_x0);

	fscanf(entrada, "%d %d %lf %lf", &M, &nvertices, &dz, &z0);

	x0_inicial = aloca_vetor_double (relatorio, M);
	y0_inicial = aloca_vetor_double (relatorio, M);
	x0 = aloca_vetor_double (relatorio, M);
	y0 = aloca_vetor_double (relatorio, M);
	x0_trans = aloca_vetor_double (relatorio, M);
	y0_trans = aloca_vetor_double (relatorio, M);
	rho = aloca_vetor_double (relatorio, M);
	z1 = aloca_vetor_double (relatorio, M);
	z2 = aloca_vetor_double (relatorio, M);

	raio_inicial = aloca_matriz_double (relatorio, M, nvertices);
	raio = aloca_matriz_double (relatorio, M, nvertices);
	raio_trans = aloca_matriz_double (relatorio, M, nvertices);
	x = aloca_matriz_double (relatorio, M, nvertices);
	y = aloca_matriz_double (relatorio, M, nvertices);

	Q = M*nvertices;
	P = Q + (2*M);
	teta = (double)((2*PI)/nvertices);

	for (i = 0; i < M; i++) {

        if(fscanf(entrada, "%lf %lf %lf %lf", &rho[i], &y0_inicial[i], &x0_inicial[i], &r0) != 4) {

			fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", str);

			fclose (relatorio);

			printf ("Erro!\n\n");

			system ("PAUSE");

			return 0;

        }

        for (j = 0; j < nvertices; j++) {

			raio_inicial[i][j] = r0;

        }

		z1[i] = z0 + (i*dz);
		z2[i] = z1[i] + dz;

	}

	fclose (entrada);

	/* impressao do modelo inicial */
	contador_repetibilidade = 9999;
	paraview (M, Q, nvertices, teta, raio_inicial, z1, z2, x0_inicial, y0_inicial, contador_repetibilidade);

	/******** <== leitura do arquivo da estimativa inicial **************/

	/******** leitura do arquivo de dados observados ==> **************/

	sprintf (str, "gobs.txt");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	entrada = fopen(str, "r");

	fscanf(entrada, "%d", &N);

	Xp = aloca_vetor_double(relatorio, N);
	Yp = aloca_vetor_double(relatorio, N);
	Zp = aloca_vetor_double(relatorio, N);
	gobs_original = aloca_vetor_double(relatorio, N);

	for (i = 0; i < N; i++) {

		if (fscanf(entrada, "%lf %lf %lf %lf", &Yp[i], &Xp[i], &Zp[i], &gobs_original[i]) != 4) {

			fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", str);

			fclose (relatorio);

			printf ("Erro!\n\n");

			system ("PAUSE");

			return 0;

		}

		gobs_original[i] += stdev_dados*gasdev(&idum);

	}

	fclose (entrada);

	/******** <== leitura do arquivo de dados observados **************/

	/******* leitura do contorno da parte aflorante do corpo e construção do vetor pae ==> ********/

	sprintf (str, "contorno.txt", "r");

	if (fopen(str, "r") == NULL) {

		fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

		fclose (relatorio);

		printf ("Erro!\n\n");

		system ("PAUSE");

		return 0;

	}

	entrada = fopen(str, "r");

	fscanf (entrada, "%d %lf %lf %lf %lf %lf", &L, &uae_raios, &y0_afloramento, &x0_afloramento, &uae_x0, &uae_y0);

	pae = aloca_vetor_double (relatorio, nvertices);

	if (L != 0) {

       	contorno = aloca_matriz_double (relatorio, L, 4);
		npontos = aloca_vetor_int (relatorio, nvertices);

		for (i = 0; i < L; i++) {

			if (fscanf(entrada, "%lf %lf", &contorno[i][0], &contorno[i][1]) != 2) {

				fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", str);

				fclose (relatorio);

				printf ("Erro!\n\n");

				system ("PAUSE");

				return 0;

			}

			/*aux0 = contorno[i][0] - y0_inicial[0];*/
			aux0 = contorno[i][0] - y0_afloramento;
			/*aux1 = contorno[i][1] - x0_inicial[0];*/
			aux1 = contorno[i][1] - x0_afloramento;
			contorno[i][2] = pow((pow(aux0, 2) + pow(aux1, 2)), 0.5);
			contorno[i][3] = -atan2(aux1,aux0);

			if (contorno[i][3] < 0.0) {

				contorno[i][3] += 2*PI;

			}

			aux2 = (double)(contorno[i][3]/teta);

			j = (int)(aux2);

			if (aux2 < ((j + 0.5)*teta)) {

				npontos[j] ++;
				pae[j] += contorno[i][2];

			}
			else {

				if ((j+1) == nvertices) {

					npontos[0] ++;
					pae[0] += contorno[i][2];

				}
				else {

					npontos[j+1] ++;
					pae[j+1] += contorno[i][2];

				}

			}

		}

		for (i = 0; i < nvertices; i++) {

			pae[i] *= (double)(1.0/npontos[i]);

		}

		fclose (entrada);

		npontos = libera_vetor_int(npontos);
		contorno = libera_matriz_double(L, contorno);

	}

	/******* <== leitura do contorno da parte aflorante do corpo e construção do vetor pae ********/

	gcalc = aloca_vetor_double (relatorio, N);
	gobs = aloca_vetor_double (relatorio, N);

	x0_medio = aloca_vetor_double (relatorio, M);
	y0_medio = aloca_vetor_double (relatorio, M);
	raio_medio = aloca_matriz_double (relatorio, M, nvertices);

	fclose (relatorio);

	for (contador_repetibilidade = 1; contador_repetibilidade <= ITMAX_repetibilidade; contador_repetibilidade++) { /* <== looping repetibilidade */

	printf("\n\n*************** Teste %3d ********************\n\n", contador_repetibilidade);

	/* Perturbação dos dados */
	for (i = 0; i < N; i++) {

		gobs[i] = gobs_original[i] + stdev_perturbacao*gasdev(&idum);

	}

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

			raio[i][j] = raio_inicial[i][j];

		}
		
		x0[i] = x0_inicial[i];
		y0[i] = y0_inicial[i];
		
	}

	sprintf(str, "relatorio%d.txt", contador_repetibilidade);

	relatorio = fopen(str, "w");

	/******** inversão ==> ********************************************/
	dp = aloca_matriz_double (relatorio, M, (nvertices + 2));
	grad = aloca_vetor_double (relatorio, P);
	H = aloca_matriz_double (relatorio, P, P);
	grad_ajuste = aloca_vetor_double (relatorio, P);
	H_ajuste = aloca_matriz_double (relatorio, P, P);

	/* gcalc avaliado em raio_inicial */
	/* calculo do vetor residuo */
	/* funcao objetiva avaliada em raio_inicial */
	ajuste = mod_direto (N, M, P, nvertices, Xp, Yp, Zp, x0, y0, rho, z1, z2, raio, teta, x, y, gcalc, gobs);

	fobj2 = ajuste;

	/* Contribuição do Ridge dos raios */
	fobj2 += fobj_ridge(uridge, raio, M, nvertices);

	/* Contribuição do Absolut Equality dos raios */
	fobj2 += fobj_absolut_equality_raios (uae_raios, raio, pae, nvertices);
	
	/* Contribuição do Absolut Equality das coordenadas x0 e y0 */
	fobj2 += fobj_absolut_equality_origem (uae_x0, uae_y0, x0_afloramento, y0_afloramento, x0, y0, M);

	/* Contribuição do Flatness dos raios */
	fobj2 += fobj_flatness_raios (uflatness_rad, uflatness_vert, raio, M, nvertices);
	
	/* Contribuição do Flatness das coordenadas x0 e y0 */
	fobj2 += fobj_flatness_origens (uflatness_x0, uflatness_y0, x0, y0, M);
	
	iteracao = 0;
	
	lambida = lambida_inicial;

	fprintf (relatorio, "%5d %15.5E %10.3E      --    \n", iteracao, fobj2, ajuste);
	printf ("\n%5d %10.3E %10.3E      --    \n", iteracao, fobj2, ajuste);

	iteracao = 1;

	do {

		fobj1 = fobj2;

		/* calculo do vetor grad_ajuste */
		/* calculo da matriz H_ajuste */
		derivada_ajuste(relatorio, N, M, nvertices, P, H_ajuste, grad_ajuste, gobs, gcalc, raio, dr, teta, rho, z1, z2, x, y, x0, y0, Xp, Yp, Zp);

		/* cálculo do vetor grad e da matriz H ==> */

		/* incorporação da parcela proveniente do funcional do ajuste ==> */
		for (i = 0; i < P; i++) {

			grad[i] = grad_ajuste[i];

			for (j = 0; j < P; j++) {

                H[i][j] = H_ajuste[i][j];
                
			}

		}
		/* <== incorporação da parcela proveniente do funcional do ajuste */

		/* incorporação do vínculo Ridge nos raios ==> */
		for (i = 0, k = 0; i < M; i++) {
		
			for (j = 0; j < nvertices; j++, k++) {

				grad[k] += (uridge*raio[i][j]);

				/* diagonal principal */
				H[k][k] += uridge;

			}

			k += 2; /* para pular x0 e y0 */

		}
        /* <== incorporação do vínculo Ridge nos raios */

        /* incorporação do vínculo Absolut Equality nos raios do prisma da superfície ==> */
		if (uae_raios != 0.0) {

			for (i = 0; i < nvertices; i++) {

				grad[i] += (uae_raios*(raio[0][i] - pae[i]));

				H[i][i] += uae_raios;

			}

		}
		/* <== incorporação do vínculo Absolut Equality nos raios do prisma da superfície */

        /* incorporação do vínculo Absolut Equality nas coordenadas x0 e y0 do prisma da superfície ==> */
		if ((uae_x0 != 0.0) && (uae_y0 != 0.0)) {

			grad[nvertices] += uae_x0*(x0[0] - x0_afloramento);

			grad[nvertices + 1] += uae_y0*(y0[0] - y0_afloramento);

			H[nvertices][nvertices] += uae_x0;

			H[nvertices + 1][nvertices + 1] += uae_y0;

		}
		/* <== incorporação do vínculo Absolut Equality nas coordenadas x0 e y0 do prisma da superfície */

        /* incorporação do vínculo Flatness nos raios intraplaca ==> */
        if (uflatness_rad != 0.0) {

			for (i = 0, k = 0; i < M; i++) {

				j = 0;

				grad[k] += (uflatness_rad*((2*raio[i][j]) + (-raio[i][j + 1]) + (-raio[i][nvertices - 1])));

				H[k][k] += (2*uflatness_rad);
				H[k][k + 1] += -uflatness_rad;
				H[k][k + nvertices - 1] += -uflatness_rad;

				k++;

				for (j = 1; j < (nvertices - 1); j++, k++) {

	                grad[k] += (uflatness_rad*((2*raio[i][j]) + (-raio[i][j - 1]) + (-raio[i][j + 1])));

					H[k][k] += (2*uflatness_rad);
					H[k][k - 1] += -uflatness_rad;
					H[k][k + 1] += -uflatness_rad;


				}

				grad[k] += (uflatness_rad*((2*raio[i][j]) + (-raio[i][j - 1]) + (-raio[i][0])));

				H[k][k] += (2*uflatness_rad);
				H[k][k - 1] += -uflatness_rad;
				H[k][k - nvertices + 1] += -uflatness_rad;

				k++;

                k += 2; /* pula x0 e y0 */

			}

		}

		/* <== incorporação do vínculo Flatness nos raios intraplaca */

        /* incorporação do vínculo Flatness nos raios entre placas e nas
		coordenadas x0 e y0 entre placas ==> */
		if ((uflatness_vert != 0.0) && (M != 1)) {

			i = 0;
			k = 0;

				/* vínculo Flatness nos raios entre placas */
				for (j = 0; j < nvertices; j++, k++) {

					grad[k] += (uflatness_vert*(raio[i][j] - raio[i + 1][j]));

					H[k][k] += uflatness_vert;
					H[k][k + nvertices + 2] += -uflatness_vert;

				}

				/* vínculo Flatness entre as coordenadas x0 e y0 */
				grad[k] += (uflatness_x0*(x0[i] - x0[i + 1]));

				H[k][k] += uflatness_x0;
				H[k][k + nvertices + 2] += -uflatness_x0;

				k++;

				grad[k] += (uflatness_y0*(y0[i] - y0[i + 1]));

				H[k][k] += uflatness_y0;
				H[k][k + nvertices + 2] += -uflatness_y0;

				k++;

			for (i = 1; i < (M - 1); i++) {

                /* vínculo Flatness nos raios entre placas */
				for (j = 0; j < nvertices; j++, k++) {

					grad[k] += (uflatness_vert*(-raio[i - 1][j] + 2*raio[i][j] - raio[i + 1][j]));

					H[k][k - nvertices - 2] += -uflatness_vert;
					H[k][k] += (2*uflatness_vert);
					H[k][k + nvertices + 2] += -uflatness_vert;

				}

				/* vínculo Flatness entre as coordenadas x0 e y0 */
				grad[k] += (uflatness_x0*(-x0[i - 1] + 2*x0[i] - x0[i + 1]));

				H[k][k - nvertices - 2] += -uflatness_x0;
				H[k][k] += (2*uflatness_x0);
				H[k][k + nvertices + 2] += -uflatness_x0;

                k++;

				grad[k] += (uflatness_y0*(-y0[i - 1] + 2*y0[i] - y0[i + 1]));

				H[k][k - nvertices - 2] += -uflatness_y0;
				H[k][k] += (2*uflatness_y0);
				H[k][k + nvertices + 2] += -uflatness_y0;

                k++;

			}

            /* vínculo Flatness nos raios entre placas */
			for (j = 0; j < nvertices; j++, k++) {

				grad[k] += (uflatness_vert*(-raio[i - 1][j] + raio[i][j]));

				H[k][k - nvertices - 2] += -uflatness_vert;
				H[k][k] += uflatness_vert;

    		}

			/* vínculo Flatness entre as coordenadas x0 e y0 */
			grad[k] += (uflatness_x0*(-x0[i - 1] + x0[i]));

			H[k][k - nvertices - 2] += -uflatness_x0;
			H[k][k] += uflatness_x0;

			k++;

			grad[k] += (uflatness_y0*(-y0[i - 1] + y0[i]));

			H[k][k - nvertices - 2] += -uflatness_y0;
			H[k][k] += uflatness_y0;

			k++;

		}
		/* <== incorporação do vínculo Flatness nos raios entre placas e nas
		coordenadas x0 e y0 entre placas*/

		/* <== cálculo do vetor grad e da matriz H */
		
		/* incorporação do vínculo de positividade ==> */
		
		/* cálculo dos parâmetros transformados */
		for (i = 0; i < M; i++) {
		
			x0_trans[i] = -log((double)((x0max - x0[i])/(x0[i] - x0min)));
			y0_trans[i] = -log((double)((y0max - y0[i])/(y0[i] - y0min)));

			for (j = 0; j < nvertices; j++) {
			
				raio_trans[i][j] = -log((double)((rmax - raio[i][j])/(raio[i][j] - rmin)));
			
			}
		
		}
		
		/* cálculo da hessiana transformada */
		for (k = 0, i = 0; i < M; i++) {
		
			for (j = 0; j < nvertices; j++, k++) {
			
				aux0 = (double)(((raio[i][j] - rmin + 1E-5)*(rmax + 1E-5 - raio[i][j]))/(rmax - rmin));
				
				for (l = 0; l < P; l++) {
				
					H[l][k] *= aux0;
				
				}
			
			}
			
			aux0 = (double)(((x0[i] - x0min + 1E-5)*(x0max + 1E-5 - x0[i]))/(x0max - x0min));

			for (l = 0; l < P; l++) {
			
				H[l][k] *= aux0;
			
			} 			
			
			k++;
			
			aux0 = (double)(((y0[i] - y0min + 1E-5)*(y0max + 1E-5 - y0[i]))/(y0max - y0min));

			for (l = 0; l < P; l++) {
			
				H[l][k] *= aux0;
			
			} 			
			
			k++;
			
		}
		/* <== incorporação do vínculo de positividade */

		/* diminuir lambida porque este é aumentado na primeira iteração do looping */
		lambida = (double)(lambida/dlambida);

		iteracao_marq = 0;

		do {

			/* aumenta lambida Marquardt */
			lambida *= dlambida;

			/* resolve sistema */
			/* calcula raio k+1 */
			sist_lin_LU (relatorio, M, P, nvertices, H, grad, dp, lambida);

			/* cálculo de raio_trans + dp ==> */
			for (i = 0; i < M; i++) {
			
				x0_trans[i] += dp[i][nvertices];
				
				x0_trans[i] = x0min + ((double)((x0max - x0min)/(1 + exp(-x0_trans[i]))));

				y0_trans[i] += dp[i][nvertices+1];
				
				y0_trans[i] = y0min + ((double)((y0max - y0min)/(1 + exp(-y0_trans[i]))));

				for (j = 0; j < nvertices; j++) {
				
					raio_trans[i][j] += dp[i][j];
					
					raio_trans[i][j] = rmin + ((double)((rmax - rmin)/(1 + exp(-raio_trans[i][j]))));
				
				}
			
			}
			/* <== cálculo de raio_trans + dp */

			/* funcao objetiva avaliada em raio k+1 */
			ajuste = mod_direto (N, M, P, nvertices, Xp, Yp, Zp, x0_trans, y0_trans, rho, z1, z2, raio_trans, teta, x, y, gcalc, gobs);

			fobj2 = ajuste;

			/* Contribuição do Ridge dos raios */
			fobj2 += fobj_ridge (uridge, raio_trans, M, nvertices);

			/* Contribuição do Absolut Equality dos raios */
			fobj2 += fobj_absolut_equality_raios (uae_raios, raio_trans, pae, nvertices);
	
			/* Contribuição do Absolut Equality das coordenadas x0 e y0 */
			fobj2 += fobj_absolut_equality_origem (uae_x0, uae_y0, x0_afloramento, y0_afloramento, x0_trans, y0_trans, M);
			
			/* Contribuição do Flatness dos raios */
			fobj2 += fobj_flatness_raios (uflatness_rad, uflatness_vert, raio_trans, M, nvertices);

			/* Contribuição do Flatness das coordenadas x0 e y0 */
			fobj2 += fobj_flatness_origens (uflatness_x0, uflatness_y0, x0_trans, y0_trans, M);

			iteracao_marq ++;

            fprintf (relatorio, "%5d %10.3E %10.3E %10.3E %5d\n", iteracao, fobj2, ajuste, lambida, iteracao_marq);
            printf ("%5d %10.3E %10.3E %10.3E %5d\n", iteracao, fobj2, ajuste, lambida, iteracao_marq);

		} while ((fobj2 >= fobj1) && (iteracao_marq < ITMAX_marq));

		/* diminuir lambida */
		lambida = (double)(lambida/dlambida);

		if (lambida < 1E-10) {

			lambida = 1E-10;

		}
		
		if (lambida > 1E10) {
		
			lambida = 1E10;
		
		}

		variacao_relativa = (double)((fobj2 - fobj1)/fobj1);
		variacao_relativa = fabs(variacao_relativa);

		iteracao ++;
		
		if (fobj2 < fobj1) {
		
			/* atualização dos parâmetros ==> */
			for(i = 0; i < M; i++) {

				x0[i] = x0_trans[i];
				y0[i] = y0_trans[i];

				for(j = 0; j < nvertices; j++) {

					raio[i][j] = raio_trans[i][j];

				}

			}
			/* <== atualização dos parâmetros */

		}

	} while (((fobj2 < fobj1) && (variacao_relativa >= epsilon) && (iteracao <= ITMAX)) || (iteracao < 6));

	if (fobj2 > fobj1) {

		printf ("\nfobj aumentou!\n\n");
		fprintf (relatorio, "\nfobj aumentou!\n\n");

	}
	if (fobj2 == fobj1) {

		printf ("\nfobj igual!\n\n");
		fprintf (relatorio, "\nfobj igual!\n\n");

	}

	if (variacao_relativa < epsilon) {

		printf ("\nvariacao_relativa < epsilon!\n\n");
		fprintf (relatorio, "\nvariacao_relativa < epsilon!\n\n");

	}

	/******** <== inversão ********************************************/

	ajuste = mod_direto (N, M, P, nvertices, Xp, Yp, Zp, x0, y0, rho, z1, z2, raio, teta, x, y, gcalc, gobs);	
	
	massa_pre = 0.0;

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

			aux0 = teta*j;

			x[i][j] = -(raio[i][j]*sin(aux0)); /* considera o ponto x0[i], y0[i] como origem */
			y[i][j] =  (raio[i][j]*cos(aux0)); /* considera o ponto x0[i], y0[i] como origem */

		}
		
		aux0 = 0.0;
		
		for (j = 0; j < (nvertices - 1); j++) {

			aux0 += (x[i][j]*y[i][j+1]) - (x[i][j+1]*y[i][j]);

		}

		j = (nvertices - 1);

		aux0 += (x[i][j]*y[i][0]) - (x[i][0]*y[i][j]);

		massa_pre += aux0*rho[i];

	}

	massa_pre *= (dz*0.5);
	massa_pre_media += massa_pre;

	ajuste_medio += ajuste;

	for (s = 0, i = 0; i < N; i++) {
	
		s += fabs((gobs[i] - gcalc[i]));
	
	}
	
	s = (double)(s/N);
	
	s_medio += s;

	paraview (M, Q, nvertices, teta, raio, z1, z2, x0, y0, contador_repetibilidade);

	dp = libera_matriz_double (M, dp);
	grad = libera_vetor_double (grad);
	H = libera_matriz_double (P, H);
	grad_ajuste = libera_vetor_double (grad_ajuste);
	H_ajuste = libera_matriz_double (P, H_ajuste);

	fclose(relatorio);

	sprintf(str, "ajuste%d.txt", contador_repetibilidade);

	saida = fopen(str, "w");

	for (i = 0; i < N; i++) {

		fprintf (saida, "%15.3lf %15.3lf %15.3lf %15.5lf %15.5lf\n", Yp[i], Xp[i], Zp[i], gobs[i], gcalc[i]);

	}

	fprintf (saida, "\n%15.5E\n", ajuste);
	fprintf (saida, "%15.5E\n", s);

	fclose (saida);

	sprintf(str, "model%d.txt", contador_repetibilidade);

	saida = fopen(str, "w");

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

			fprintf (saida, "%15.5E\n", raio[i][j]);

			raio_medio[i][j] += raio[i][j];

		}

		fprintf (saida, "\n%15.5E\n", x0[i]);
		fprintf (saida, "%15.5E\n\n", y0[i]);

		x0_medio[i] += x0[i];
		y0_medio[i] += y0[i];

	}

	fprintf (saida, "%15.5E\n", massa_pre);

	fclose(saida);

	} /* <== looping repetibilidade */

	/* Cálculo da média dos parâmetros ==> */
	aux0 = (double)(1.0/ITMAX_repetibilidade);

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

			raio_medio[i][j] *= aux0;

		}

		x0_medio[i] *= aux0;
		y0_medio[i] *= aux0;

	}

	massa_pre_media *= aux0;

	ajuste_medio *= aux0;

	s_medio *= aux0;

	/* <== Cálculo da média dos parâmetros */

	/* Cálculo do desvio padrão dos parâmetros ==> */
	relatorio = fopen ("relatorio_final.txt", "w");
	
	raio_std = aloca_matriz_double (relatorio, M, nvertices);
	x0_std = aloca_vetor_double (relatorio, M);
	y0_std = aloca_vetor_double (relatorio, M);

	for (i = 1; i <= ITMAX_repetibilidade; i++) {

		sprintf(str, "model%d.txt", i);
		
		entrada = fopen(str, "r");
		
		for (j = 0; j < M; j++) {

			for (k = 0; k < nvertices; k++) {

				fscanf (entrada, "%lf", &raio[j][k]);

				raio_std[j][k] += pow ((raio[j][k] - raio_medio[j][k]), 2);

			}
			
			fscanf (entrada, "%lf", &x0[j]);

			x0_std[j] += pow((x0[j] - x0_medio[j]), 2);
			
			fscanf (entrada, "%lf", &y0[j]);

			y0_std[j] += pow((y0[j] - y0_medio[j]), 2);

		}

		fscanf (entrada, "%lf", &aux0);
		fscanf (entrada, "%lf", &massa_pre);

		massa_pre_std += pow((massa_pre - massa_pre_media), 2);

		fclose(entrada);

	}

	for (i = 1; i <= ITMAX_repetibilidade; i++) {

		sprintf(str, "ajuste%d.txt", i);
		
		entrada = fopen(str, "r");
		
		for (j = 0; j < N; j++) {

			fscanf (entrada, "%lf %lf %lf %lf", &aux0, &aux1, &aux2, &aux3);

		}

		fscanf (entrada, "%lf", &ajuste);
		fscanf (entrada, "%lf", &s);

		ajuste_std += pow((ajuste - ajuste_medio), 2);

		s_std += pow((s - s_medio), 2);

		fclose(entrada);

	}


	aux0 = (double)(1.0/(ITMAX_repetibilidade - 1));

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

			raio_std[i][j] = pow((aux0*raio_std[i][j]), 0.5);
			
		}

		x0_std[i] = pow((aux0*x0_std[i]), 0.5);
		y0_std[i] = pow((aux0*y0_std[i]), 0.5);

	}
	
	fclose (relatorio);

	massa_pre_std = pow((aux0*massa_pre_std), 0.5);
	ajuste_std = pow((aux0*ajuste_std), 0.5);
	s_std = pow((aux0*s_std), 0.5);

	/* <== Cálculo do desvio padrão dos parâmetros */

	sprintf(str, "ajuste0.txt");

	saida = fopen(str, "w");

	ajuste = mod_direto (N, M, P, nvertices, Xp, Yp, Zp, x0_medio, y0_medio, rho, z1, z2, raio_medio, teta, x, y, gcalc, gobs_original);

	s = 0.0;
	
	for (i = 0; i < N; i++) {

		fprintf (saida, "%15.3lf %15.3lf %15.3lf %15.5lf %15.5lf\n", Yp[i], Xp[i], Zp[i], gobs_original[i], gcalc[i]);

		s += fabs((gobs[i] - gcalc[i]));

	}
	
	s = (double)(s/N);

	massa_pre = 0.0;

	for (i = 0; i < M; i++) {

		aux0 = 0.0;

		for (j = 0; j < (nvertices - 1); j++) {

			aux0 += (raio_medio[i][j]*raio_medio[i][j+1]);

		}

		aux0 += (raio_medio[i][j]*raio_medio[i][0]);

		massa_pre += aux0*rho[i];

	}

	massa_pre *= (dz*sin(teta)*0.5);

	fprintf(saida, "\n");
	fprintf(saida, "%15.5lf\n", ajuste);
	fprintf(saida, "%15.5lf  %8.5lf\n\n", ajuste_medio, ajuste_std);
	fprintf(saida, "%15.5lf\n", s);
	fprintf(saida, "%15.5lf  %8.5lf\n\n", s_medio, s_std);
    fprintf(saida, "%15.5lf\n", massa_pre);
	fprintf(saida, "%15.5lf  %8.5lf\n\n", massa_pre_media, massa_pre_std);

	fprintf(saida, "ajuste modelo medio\n");
	fprintf(saida, "ajuste medio (std)\n\n");
	fprintf(saida, "s modelo medio\n");
	fprintf(saida, "s medio (std)\n\n");
    fprintf(saida, "massa predita modelo medio\n");
	fprintf(saida, "massa predita media (std)\n\n");

	fclose(saida);

	contador_repetibilidade = 0;

	paraview (M, Q, nvertices, teta, raio_medio, z1, z2, x0_medio, y0_medio, contador_repetibilidade);

	paraview_error (M, nvertices, teta, raio_medio, raio_std, z1, x0_medio, x0_std, y0_medio, y0_std, contador_repetibilidade);

	fclose (saida);

	time (&end);
	
	aux0 = difftime (end, start);

	printf("\nPrograma finalizado com sucesso em %.3lf segundos!\n\n", aux0);

	system("PAUSE");

	return 0;

}

double **aloca_matriz_double (FILE *arq, int linha, int coluna) {

    double **m;  /* ponteiro para a matriz */
    int   i, j;

    /* aloca as linhas da matriz */

    m = (double **)calloc(linha, sizeof(double *));

    if (m == NULL) {

        fprintf (arq, "Memoria Insuficiente (linhas)!\n\n");
        printf ("\n\nMemoria Insuficiente (linhas)!\n\n");

        fclose (arq);

        system ("PAUSE");
        return 0;

        return (NULL);

    }

    /* aloca as colunas da matriz */

    for ( i = 0; i < linha; i++ ) {

        m[i] = (double *)calloc(coluna, sizeof(double));

        if (m[i] == NULL) {

            fprintf (arq, "Memoria Insuficiente (colunas)!\n\n");

            fclose (arq);

            system ("PAUSE");
            return 0;

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

        return (NULL);

        fclose (arq);

        system ("PAUSE");
        return 0;

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

int *aloca_vetor_int (FILE *arq, int tamanho) {

    int *v; /* ponteiro para o vetor */

    v = (int *)calloc(tamanho, sizeof(int));

    if (v == NULL) { /*** verifica se há memória suficiente ***/

        fprintf (arq, "Memoria Insuficiente!\n\n");

        return (NULL);

        fclose (arq);

        system ("PAUSE");
        return 0;

    }

    return (v); /* retorna o ponteiro para o vetor */

}

int *libera_vetor_int (int *v) {

    if (v == NULL) {

        return (NULL);

    }

    free(v); /* libera o vetor */

    return (NULL); /* retorna o ponteiro */

}

double mod_direto (int N, int M, int P, int nvertices, double *Xp, double *Yp, double *Zp,
					double *x0, double *y0, double *rho, double *z1, double *z2,
					double **raio, double teta, double **x, double **y,
					double *gcalc, double *gobs) {

	int i, j, k;
	double gaux, aux, auxk1, auxk2, aux1k1, aux1k2, aux2k1, aux2k2;
	double Ak1, Ak2, Bk1, Bk2, Ck1, Ck2, Dk1, Dk2, E1k1, E1k2, E2k1, E2k2;
	double Xk1, Xk2, Yk1, Yk2, Z1, Z1_quadrado, Z2, Z2_quadrado, Qk1, Qk2;
	double R1k1, R1k2, R2k1, R2k2, p, p_quadrado, tetak;
	double ajuste;

	ajuste = 0.0;

	/******** Cálculo da atração gravitacional ==> **********/

	/* O cálculo das coordenadas cartesianas referente aos vétices dos polígonos
	é feito durante o cálculo da atração gravitacional do primeiro ponto, i = 0.
	Essas coordenas novas são utilizadas nos demais pontos. */

	/*for (i = 1; i < N; i++) { <== for dos pontos */
	i = 0;

		gcalc[i] = 0.0;

		for (j = 0; j < M; j++) { /* for dos prismas */

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */
			
			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];
			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			for (k = 0; k < (nvertices-1); k++) { /* <== for dos vértices */

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
       			tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

                tetak = teta*(k+1);

				x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
				y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

				Xk1 = x[j][k] - Xp[i];
				Xk2 = x[j][k+1] - Xp[i];
				Yk1 = y[j][k] - Yp[i];
				Yk2 = y[j][k+1] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			} /* <== for dos vértices */

			/***** última aresta ==> *******/

			Xk1 = x[j][nvertices-1] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][nvertices-1] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== última aresta *******/

			gaux *= aux;

			gcalc[i] += gaux;

		} /* for dos prismas */

		ajuste += pow((gobs[i] - gcalc[i]), 2);
		
	/*} <== for dos pontos */

	for (i = 1; i < N; i++) { /* <== for dos pontos */

		gcalc[i] = 0.0;

		for (j = 0; j < M; j++) { /* for dos prismas */

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];
			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			for (k = 0; k < (nvertices-1); k++) { /* <== for dos vértices */

				Xk1 = x[j][k] - Xp[i];
				Xk2 = x[j][k+1] - Xp[i];
				Yk1 = y[j][k] - Yp[i];
				Yk2 = y[j][k+1] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			} /* <== for dos vértices */

			/***** última aresta ==> *******/

			Xk1 = x[j][nvertices-1] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][nvertices-1] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== última aresta *******/

			gaux *= aux;

			gcalc[i] += gaux;

		} /* for dos prismas */

		ajuste += pow((gobs[i] - gcalc[i]), 2);
		
	} /* <== for dos pontos */

	/******** <== Cálculo da atração gravitacional **********/

	ajuste *= fabs((double)(1.0/(N-P)));
	/*ajuste *= ((double)(1.0/N));*/

	return ajuste;

}

double sist_lin_LU (FILE *relatorio, int M, int P, int nvertices, double **H, double *grad, double **dp, double lambida) {

	int i, j, k, l;
	double *b, **L, *solucao;

	/*************** decomposição LU da matriz H ==> ********************/
	L = aloca_matriz_double (relatorio, P, P);

	/**** primeira linha de U, primeira coluna de L e diagonal principal de L ****/

	L[0][0] = H[0][0] + lambida;

	for (i = 1; i < P; i++) {

		L[0][i] = H[0][i];
		L[i][0] = (double)(H[i][0]/L[0][0]);

	}

	/**** varre as linhas de U e as colunas de L ****/
	for (i = 1; i < P; i++) {

		/**** i-ésima linha de U ****/
		j = i;

			for (k = 0; k < i; k++) {

				L[i][j] -= L[i][k]*L[k][j];

			}

			L[i][j] += (H[i][j] + lambida);

		for (j = i+1; j < P; j++) {

			for (k = 0; k < i; k++) {

				L[i][j] -= L[i][k]*L[k][j];

			}

			L[i][j] += H[i][j];

		}

		/**** i-ésima coluna de L ****/
		for (j = i+1; j < P; j++) {

			for (k = 0; k < i; k++) {

				L[j][i] -= L[j][k]*L[k][i];

			}

			L[j][i] += H[j][i];

			L[j][i] = (double)(L[j][i]/L[i][i]);

		}

	}

	/*************** <== decomposição LU da matriz H ********************/

	/************ resolução do sistema linear Hdp = -grad ==> **************/

	b = aloca_vetor_double(relatorio, P);

	/**** Hdp = -grad, LUdp = -grad, Lb = -grad, Udp = b ****/

	/**** resolução do sistema Lb = -grad ***/
	for (k = 0; k < P; k++) {

			b[k] = 0.0;

			for (l = 0; l < k; l++) {

				b[k] -= L[k][l]*b[l];

			}

			b[k] -= grad[k]; /* esse é o menos que precede o grad */

	}

	/*** resolução do sistema Udraio = b ***/

	solucao = aloca_vetor_double (relatorio, P);

	for (j = P-1; j >= 0; j--) {

		solucao[j] = 0;

		for (k = j+1; k < P; k++) {

			solucao[j] -= L[j][k]*solucao[k];

		}

		solucao[j] += b[j];

		solucao[j] = (double)(solucao[j]/L[j][j]);

	}

	for (i = 0, k = 0; i < M; i++) {

		for (j = 0; j < (nvertices + 2); j++, k++) {

			dp[i][j] = solucao[k];

		}

	}

	/************ <== resolução do sistema linear Hdraio = -grad **************/

	b = libera_vetor_double(b);
	L = libera_matriz_double (P, L);
	solucao = libera_vetor_double(solucao);

	/************ <== resolução do sistema linear Hdraio = -grad **************/

	return 0 ;

}

double derivada_ajuste_raio (int ponto, int poligono, int vertice, double *Xp, double *Yp, double *Zp,
							double *x0, double *y0, double *rho, double *z1, double *z2,
							int nvertices, double **raio, double teta, double **x, double **y,
							double dr) {

	int i, j, k;
	double gaux, aux, auxk1, auxk2, aux1k1, aux1k2, aux2k1, aux2k2;
	double Ak1, Ak2, Bk1, Bk2, Ck1, Ck2, Dk1, Dk2, E1k1, E1k2, E2k1, E2k2;
	double Xk1, Xk2, Yk1, Yk2, Z1, Z1_quadrado, Z2, Z2_quadrado, Qk1, Qk2;
	double R1k1, R1k2, R2k1, R2k2, p, p_quadrado, tetak;

	double gcalc_perturb, gcalc, derivada;

	if (vertice == (nvertices - 1)) {

		/*for (i = 1; i < N; i++) { <== for dos pontos */
		i = ponto;

		gcalc_perturb = 0.0;

		gcalc = 0.0;

		/*for (j = 0; j < M; j++) {*/ /* for dos prismas */
		j = poligono;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];

			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			/*for (k = 0; k < (nvertices-1); k++) {*/ /* <== for dos vértices */
			k = vertice;

			/* Cálculo do modelo direto perturbado ==> */

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
       			tetak = teta*(k-1);

				x[j][k-1] = x0[j] - (raio[j][k-1]*sin(tetak));
				y[j][k-1] = y0[j] + (raio[j][k-1]*cos(tetak));

                tetak = teta*k;

				x[j][k] = x0[j] - ((raio[j][k] + dr)*sin(tetak));
				y[j][k] = y0[j] + ((raio[j][k] + dr)*cos(tetak));

				Xk1 = x[j][k-1] - Xp[i];
				Xk2 = x[j][k] - Xp[i];
				Yk1 = y[j][k-1] - Yp[i];
				Yk2 = y[j][k] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/*}*/ /* <== for dos vértices */

			/***** aresta seguinte ==> *******/

			/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
			tetak = teta*k;

			x[j][k] = x0[j] - ((raio[j][k] + dr)*sin(tetak));
			y[j][k] = y0[j] + ((raio[j][k] + dr)*cos(tetak));

			/*x[j][0] = x0[j] - (raio[j][0]*sin(0.0));*/
			x[j][0] = x0[j];

            /*y[j][0] = y0[j] + (raio[j][0]*cos(0.0));*/
            y[j][0] = y0[j] + raio[j][0];

			Xk1 = x[j][k] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][k] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== aresta seguinte *******/

			gaux *= aux;

			gcalc_perturb += gaux;

			/* Cálculo do modelo direto perturbado ==> */

			/* Cálculo do modelo direto não perturbado ==> */

				gaux = 0.0;

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
				tetak = teta*(k-1);

				x[j][k-1] = x0[j] - (raio[j][k-1]*sin(tetak));
				y[j][k-1] = y0[j] + (raio[j][k-1]*cos(tetak));

                tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

				Xk1 = x[j][k-1] - Xp[i];
				Xk2 = x[j][k] - Xp[i];
				Yk1 = y[j][k-1] - Yp[i];
				Yk2 = y[j][k] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/*}*/ /* <== for dos vértices */

			/***** aresta seguinte ==> *******/

			/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
            tetak = teta*k;

			x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
			y[j][k] = y0[j] + (raio[j][k]*cos(tetak));


			/*x[j][0] = x0[j] - (raio[j][0]*sin(0.0));*/
			x[j][0] = x0[j];

            /*y[j][0] = y0[j] + (raio[j][0]*cos(0.0));*/
            y[j][0] = y0[j] + raio[j][0];

			Xk1 = x[j][k] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][k] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== aresta seguinte *******/

			gaux *= aux;

			gcalc += gaux;

			/* Cálculo do modelo direto não perturbado ==> */

			derivada = (double)((gcalc_perturb - gcalc)/dr);

	}

	if (vertice == 0) {
	
		/*for (i = 1; i < N; i++) { <== for dos pontos */
		i = ponto;

		gcalc_perturb = 0.0;

		gcalc = 0.0;

		/*for (j = 0; j < M; j++) {*/ /* for dos prismas */
		j = poligono;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];

			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			/*for (k = 0; k < (nvertices-1); k++) {*/ /* <== for dos vértices */
			k = vertice;

			/* Cálculo do modelo direto perturbado ==> */

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
       			tetak = teta*(nvertices-1);

				x[j][nvertices-1] = x0[j] - (raio[j][nvertices-1]*sin(tetak));
				y[j][nvertices-1] = y0[j] + (raio[j][nvertices-1]*cos(tetak));

                tetak = teta*k;

				x[j][k] = x0[j] - ((raio[j][k] + dr)*sin(tetak));
				y[j][k] = y0[j] + ((raio[j][k] + dr)*cos(tetak));

				Xk1 = x[j][nvertices-1] - Xp[i];
				Xk2 = x[j][k] - Xp[i];
				Yk1 = y[j][nvertices-1] - Yp[i];
				Yk2 = y[j][k] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/*}*/ /* <== for dos vértices */

			/***** aresta seguinte ==> *******/

			/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
			tetak = teta*k;

			x[j][k] = x0[j] - ((raio[j][k] + dr)*sin(tetak));
			y[j][k] = y0[j] + ((raio[j][k] + dr)*cos(tetak));

            tetak = teta*(k+1);

			x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
			y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

			Xk1 = x[j][k] - Xp[i];
			Xk2 = x[j][k+1] - Xp[i];
			Yk1 = y[j][k] - Yp[i];
			Yk2 = y[j][k+1] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== aresta seguinte *******/

			gaux *= aux;

			gcalc_perturb += gaux;

			/* <== Cálculo do modelo direto perturbado */

			/* Cálculo do modelo direto não perturbado ==> */

				gaux = 0.0;

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
				tetak = teta*(nvertices-1);

				x[j][nvertices-1] = x0[j] - (raio[j][nvertices-1]*sin(tetak));
				y[j][nvertices-1] = y0[j] + (raio[j][nvertices-1]*cos(tetak));

                tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

				Xk1 = x[j][nvertices-1] - Xp[i];
				Xk2 = x[j][k] - Xp[i];
				Yk1 = y[j][nvertices-1] - Yp[i];
				Yk2 = y[j][k] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/*}*/ /* <== for dos vértices */

			/***** aresta seguinte ==> *******/

			/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
			tetak = teta*k;

			x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
			y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

            tetak = teta*(k+1);

			x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
			y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

			Xk1 = x[j][k] - Xp[i];
			Xk2 = x[j][k+1] - Xp[i];
			Yk1 = y[j][k] - Yp[i];
			Yk2 = y[j][k+1] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== aresta seguinte *******/

			gaux *= aux;

			gcalc += gaux;

			/* Cálculo do modelo direto não perturbado ==> */

			derivada = (double)((gcalc_perturb - gcalc)/dr);
	
	}

	if ((vertice != (nvertices - 1)) && (vertice != 0)) {

		/*for (i = 1; i < N; i++) { <== for dos pontos */
		i = ponto;

		gcalc_perturb = 0.0;

		gcalc = 0.0;

		/*for (j = 0; j < M; j++) {*/ /* for dos prismas */
		j = poligono;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];

			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			/*for (k = 0; k < (nvertices-1); k++) {*/ /* <== for dos vértices */
			k = vertice;

			/* Cálculo do modelo direto perturbado ==> */

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
       			tetak = teta*(k-1);

				x[j][k-1] = x0[j] - (raio[j][k-1]*sin(tetak));
				y[j][k-1] = y0[j] + (raio[j][k-1]*cos(tetak));

                tetak = teta*k;

				x[j][k] = x0[j] - ((raio[j][k] + dr)*sin(tetak));
				y[j][k] = y0[j] + ((raio[j][k] + dr)*cos(tetak));

				Xk1 = x[j][k-1] - Xp[i];
				Xk2 = x[j][k] - Xp[i];
				Yk1 = y[j][k-1] - Yp[i];
				Yk2 = y[j][k] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/*}*/ /* <== for dos vértices */

			/***** aresta seguinte ==> *******/

			/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
			tetak = teta*k;

			x[j][k] = x0[j] - ((raio[j][k] + dr)*sin(tetak));
			y[j][k] = y0[j] + ((raio[j][k] + dr)*cos(tetak));

            tetak = teta*(k+1);

			x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
			y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

			Xk1 = x[j][k] - Xp[i];
			Xk2 = x[j][k+1] - Xp[i];
			Yk1 = y[j][k] - Yp[i];
			Yk2 = y[j][k+1] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== aresta seguinte *******/

			gaux *= aux;

			gcalc_perturb += gaux;

			/* <== Cálculo do modelo direto perturbado */

			/* Cálculo do modelo direto não perturbado ==> */

				gaux = 0.0;

				/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
				tetak = teta*(k-1);

				x[j][k-1] = x0[j] - (raio[j][k-1]*sin(tetak));
				y[j][k-1] = y0[j] + (raio[j][k-1]*cos(tetak));

                tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

				Xk1 = x[j][k-1] - Xp[i];
				Xk2 = x[j][k] - Xp[i];
				Yk1 = y[j][k-1] - Yp[i];
				Yk2 = y[j][k] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/*}*/ /* <== for dos vértices */

			/***** aresta seguinte ==> *******/

			/* Cálculo das coordenadas cartesianas dos vértices dos prismas */
			tetak = teta*k;

			x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
			y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

            tetak = teta*(k+1);

			x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
			y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

			Xk1 = x[j][k] - Xp[i];
			Xk2 = x[j][k+1] - Xp[i];
			Yk1 = y[j][k] - Yp[i];
			Yk2 = y[j][k+1] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== aresta seguinte *******/

			gaux *= aux;

			gcalc += gaux;

			/* Cálculo do modelo direto não perturbado ==> */

			derivada = (double)((gcalc_perturb - gcalc)/dr);

	}

	return derivada;

}

double derivada_ajuste(FILE *relatorio, int N, int M, int nvertices, int P, double **H_ajuste,
				double *grad_ajuste, double *gobs, double *gcalc, double **raio, double dr,
				double teta, double *rho, double *z1, double *z2, double **x, double **y,
				double *x0, double *y0, double *Xp, double *Yp, double *Zp) {

	int i, j, k, l, m;
	double *vaux, aux;

	vaux = aloca_vetor_double (relatorio, P);

	for (i = 0; i < P; i++) {

		grad_ajuste[i] = 0.0;

		for(j = 0; j < P; j++) {

			H_ajuste[i][j] = 0.0;

		}

	}

	for (i = 0; i < N; i++) {

		for (j = 0, l = 0; j < M; j++) {

			/* derivadas em relação aos raios ==> */
			for (k = 0; k < nvertices; k++, l++) {

				vaux[l] = derivada_ajuste_raio (i, j, k, Xp, Yp, Zp, x0, y0, rho, z1, z2, nvertices, raio, teta, x, y, dr);

				grad_ajuste[l] -= vaux[l]*(gobs[i] - gcalc[i]); /* esse é o menos que precede o gradiente do ajuste */

				for (m = 0; m <= l; m++) {

					H_ajuste[m][l] += vaux[m]*vaux[l];

				}

			}
			/* <== derivadas em relação aos raios */

			/* derivadas em relação a x0 e y0 ==> */

				vaux[l] = derivada_ajuste_x0 (i, j, nvertices, Xp, Yp, Zp, x0, y0, rho, z1, z2, raio, teta, x, y, dr);

				grad_ajuste[l] -= vaux[l]*(gobs[i] - gcalc[i]); /* esse é o menos que precede o gradiente do ajuste */

				for (m = 0; m <= l; m++) {

					H_ajuste[m][l] += vaux[m]*vaux[l];

				}

				l++;
				
				vaux[l] = derivada_ajuste_y0 (i, j, nvertices, Xp, Yp, Zp, x0, y0, rho, z1, z2, raio, teta, x, y, dr);

				grad_ajuste[l] -= vaux[l]*(gobs[i] - gcalc[i]); /* esse é o menos que precede o gradiente do ajuste */

				for (m = 0; m <= l; m++) {

					H_ajuste[m][l] += vaux[m]*vaux[l];

				}

				l++;

			/* <== derivadas em relação a x0 e y0 */

		}

	}

	aux = (double)(1.0/(N-P));
	aux = fabs(aux);

	for (i = 0; i < P; i++) {

		grad_ajuste[i] *= aux;

		H_ajuste[i][i] *= aux;

		for (j = i+1; j < P; j++) {

			H_ajuste[i][j] *= aux;

			H_ajuste[j][i] = H_ajuste[i][j];

		}

	}

	vaux = libera_vetor_double (vaux);

	return 0;

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

int paraview (int M, int Q, int nvertices, double teta, double **raio, double *z1, double *z2, double *x0, double *y0, int contador_adaptativo) {

	int i, j, k, l;
	double aux0, aux1, aux2;
	double xmax, xmin, ymax, ymin, zmin, zmax, raiomax, x0min, x0max, y0min, y0max;
	char str[100];

	FILE *arquivo;

	raiomax = 0.0;
	x0min = x0[0];
	x0max = x0[0];
	y0min = y0[0];
	y0max = y0[0];

	sprintf(str, "vert_fit_laterais%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", (2*Q));

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

			/* cálculo dos limites da caixa */
			if (raio[i][j] > raiomax) {

				raiomax = raio[i][j];

			}
			if (x0min > x0[i]) {

				x0min = x0[i];

			}
			if (x0max < x0[i]) {

				x0max = x0[i];

			}
			if (y0min > y0[i]) {

				y0min = y0[i];

			}
			if (y0max < y0[i]) {

				y0max = y0[i];

			}

		}

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			/*fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, (z1[i] + (0.9*(z2[i] - z1[i]))));*/
			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z2[i]);

		}


	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "POLYGONS %d %d\n\n", (Q), (Q*5));

   	for (i = 0, k = 0; i < M; i++) {

		/* faces laterais do i-ésimo prisma */
		for (j = k; j <  (k + nvertices - 1); j++) {

			fprintf (arquivo, "  4 ");
			fprintf (arquivo, "%3d ", j);
			fprintf (arquivo, "%3d ", (j + nvertices));
			fprintf (arquivo, "%3d ", (j + 1 + nvertices));
			fprintf (arquivo, "%3d\n", (j + 1));

		}

			fprintf (arquivo, "  4 ");
			fprintf (arquivo, "%3d ", j);
			fprintf (arquivo, "%3d ", (j + nvertices));
			fprintf (arquivo, "%3d ", (k + nvertices));
			fprintf (arquivo, "%3d\n", k);

		k += 2*nvertices;

	}

	fclose(arquivo);

	sprintf(str, "vert_fit_tampas%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", (2*Q + 2*M));

	for (i = 0; i < M; i++) {

		/* face superior */
        fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", x0[i], y0[i], z1[i]);

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

		}

		/* face inferior */
		/*fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", x0[i], y0[i], (z1[i] + (0.98*(z2[i] - z1[i]))));*/
        fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", x0[i], y0[i], z2[i]);

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			/*fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, (z1[i] + (0.98*(z2[i] - z1[i]))));*/
            fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z2[i]);

			/* cálculo dos limites da caixa */
			if (raio[i][j] > raiomax) {

				raiomax = raio[i][j];

			}
			if (x0min > x0[i]) {

				x0min = x0[i];

			}
			if (x0max < x0[i]) {

				x0max = x0[i];

			}
			if (y0min > y0[i]) {

				y0min = y0[i];

			}
			if (y0max < y0[i]) {

				y0max = y0[i];

			}

		}

	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "POLYGONS %d %d\n\n", (2*Q), (Q*2*4));

	for (i = 0, l = 0, k = (l + 1); i < M; i++) {

		/* face de cima do i-ésimo prisma */
		for (j = 0; j < (nvertices - 1); j++, k++) {

			fprintf (arquivo, "%6d ", 3);
			fprintf (arquivo, "%6d ", l);
			fprintf (arquivo, "%6d ", k);
			fprintf (arquivo, "%6d ", (k+1));
            fprintf (arquivo, "\n");

		}

		fprintf (arquivo, "%6d ", 3);
		fprintf (arquivo, "%6d ", l);
		fprintf (arquivo, "%6d ", k);
		fprintf (arquivo, "%6d ", (k - nvertices + 1));

		fprintf (arquivo, "\n");

		l += (nvertices + 1);
		k = (l + 1);

		/* face de baixo do i-ésimo prisma */
		for (j = 0; j < (nvertices - 1); j++, k++) {

			fprintf (arquivo, "%6d ", 3);
			fprintf (arquivo, "%6d ", l);
			fprintf (arquivo, "%6d ", k);
			fprintf (arquivo, "%6d  ", (k+1));
			fprintf (arquivo, "\n");

		}

		fprintf (arquivo, "%6d ", 3);
		fprintf (arquivo, "%6d ", l);
		fprintf (arquivo, "%6d ", k);
		fprintf (arquivo, "%6d ", (k - nvertices + 1));

		fprintf (arquivo, "\n");

		l += (nvertices + 1);
		k = (l + 1);

	}

	fclose (arquivo);

	sprintf(str, "caixa%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

    fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Caixa\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS 8 float\n\n");

	zmin = z1[0];
	zmax = 1.1*z2[M-1];
	xmin = x0min - (1.1*raiomax);
	xmax = x0max + (1.1*raiomax);
	ymin = y0min - (1.1*raiomax);
	ymax = y0max + (1.1*raiomax);

	fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmin, ymin, zmin);
	fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmax, ymin, zmin);
	fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmax, ymax, zmin);
    fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmin, ymax, zmin);
	fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmin, ymin, zmax);
	fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmax, ymin, zmax);
	fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n", xmax, ymax, zmax);
    fprintf (arquivo, "%10.3lf %10.3lf %10.3lf\n\n", xmin, ymax, zmax);

	fprintf (arquivo, "POLYGONS 2 10\n\n");
	fprintf (arquivo, "4 4 5 6 7\n");
	fprintf (arquivo, "4 2 3 7 6");

	fclose(arquivo);

	sprintf(str, "vert_fit_contorno%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", (2*Q));

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

		}

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			/*fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, (z1[i] + (0.9*(z2[i] - z1[i]))));*/
            fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z2[i]);

		}

	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "LINES %d %d\n\n", (Q*2), (Q*2*3));

	for (i = 0, k = 0; i < M; i++) {

		for (j = 0; j < (nvertices - 1); j++, k++) {

			fprintf (arquivo, "2 ");
			fprintf (arquivo, "%6d ", k);
			fprintf (arquivo, "%6d ", (k+1));
			fprintf (arquivo, "\n");

		}

		fprintf (arquivo, "2 ");
		fprintf (arquivo, "%6d ", k);
		fprintf (arquivo, "%6d ", (k - nvertices + 1));
		fprintf (arquivo, "\n");

		k++;

		for (j = 0; j < (nvertices - 1); j++, k++) {

			fprintf (arquivo, "2 ");
			fprintf (arquivo, "%6d ", k);
			fprintf (arquivo, "%6d ", (k+1));
            fprintf (arquivo, "\n");

		}

		fprintf (arquivo, "2 ");
		fprintf (arquivo, "%6d ", k);
		fprintf (arquivo, "%6d ", (k - nvertices + 1));
		fprintf (arquivo, "\n");

		k++;

	}

	fclose(arquivo);

	sprintf(str, "vert_fit_origins%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", M);

	aux0 = z1[0] + ((double)((z2[0] - z1[0])/2.0));
	aux1 = z2[0] - z1[0];

	for (i = 0; i < M; i++) {

		fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", x0[i], y0[i], aux0);

		aux0 += aux1;

	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "LINES %d %d\n\n", (M-1), (3*(M-1)));

	for (i = 0; i < (M-1); i++) {

		fprintf (arquivo, "2 ");
		fprintf (arquivo, "%6d ", i);
		fprintf (arquivo, "%6d ", (i+1));
		fprintf (arquivo, "\n");

	}

	fclose(arquivo);
	
	sprintf(str, "esqueleto_fit%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", (2*Q));

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

		}

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			/*fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, (z1[i] + (0.9*(z2[i] - z1[i]))));*/
            fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z2[i]);

		}

	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "LINES %d %d\n\n", Q, (Q*3));

	for (i = 0, k = 0; i < M; i++) {

		for (j = 0; j < (nvertices - 1); j++, k++) {

			fprintf (arquivo, "2 ");
			fprintf (arquivo, "%6d ", k);
			fprintf (arquivo, "%6d ", (k+1));
			fprintf (arquivo, "\n");

		}

		fprintf (arquivo, "2 ");
		fprintf (arquivo, "%6d ", k);
		fprintf (arquivo, "%6d ", (k - nvertices + 1));
		fprintf (arquivo, "\n");

		k++;

		for (j = 0; j < (nvertices - 1); j++, k++) {

			/*fprintf (arquivo, "2 ");
			fprintf (arquivo, "%6d ", k);
			fprintf (arquivo, "%6d ", (k+1));
            fprintf (arquivo, "\n");*/

		}

		/*fprintf (arquivo, "2 ");
		fprintf (arquivo, "%6d ", k);
		fprintf (arquivo, "%6d ", (k - nvertices + 1));
		fprintf (arquivo, "\n");*/

		k++;

	}

	fclose(arquivo);

	sprintf(str, "vert_fit_casca%d.vtk", contador_adaptativo);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", (2*Q));

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

		}

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - (raio[i][j]*sin(aux0));
			aux2 = y0[i] + (raio[i][j]*cos(aux0));

            fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z2[i]);

		}

	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "POLYGONS %d %d\n\n", ((((2*M) - 1)*nvertices) + 2), ((5*(((2*M) - 1)*nvertices)) + (2*(nvertices + 1))));

	/* face superior */
	fprintf (arquivo, "%6d ", nvertices);
	
	for (j = 0; j < nvertices; j++) {

        fprintf (arquivo, "%6d ", j);

	}

    fprintf(arquivo, "\n");

   	for (i = 0, k = 0; i < ((2*M) - 1); i++) {

		for (j = 0; j < (nvertices-1); j++, k++) {

			fprintf (arquivo, "  4 ");
			fprintf (arquivo, "%3d ", k);
			fprintf (arquivo, "%3d ", (k + nvertices));
			fprintf (arquivo, "%3d ", (k + nvertices + 1));
			fprintf (arquivo, "%3d\n", (k + 1));

		}

			fprintf (arquivo, "  4 ");
			fprintf (arquivo, "%3d ", k);
			fprintf (arquivo, "%3d ", (k + nvertices));
			fprintf (arquivo, "%3d ", (k + 1));
			fprintf (arquivo, "%3d\n", k - nvertices + 1);

		k ++;

	}

	/* face inferior */
	fprintf (arquivo, "%6d ", nvertices);
	
	for (j = 0; j < nvertices; j++, k++) {

        fprintf (arquivo, "%6d ", k);

	}

    fprintf(arquivo, "\n");

	fclose(arquivo);
	
	return 0;

}

double derivada_ajuste_x0 (int ponto, int prisma, int nvertices, double *Xp, double *Yp, double *Zp,
					double *x0, double *y0, double *rho, double *z1, double *z2,
					double **raio, double teta, double **x, double **y, double dr) {

	int i, j, k;
	double gaux, aux, auxk1, auxk2, aux1k1, aux1k2, aux2k1, aux2k2;
	double Ak1, Ak2, Bk1, Bk2, Ck1, Ck2, Dk1, Dk2, E1k1, E1k2, E2k1, E2k2;
	double Xk1, Xk2, Yk1, Yk2, Z1, Z1_quadrado, Z2, Z2_quadrado, Qk1, Qk2;
	double R1k1, R1k2, R2k1, R2k2, p, p_quadrado, tetak;

	double gcalc, gcalc_perturb, derivada;

	/******** Cálculo da derivada em relação a coordenada x0 do centro de um prisma **********/

	/* Cálculo da atração gravitacional do prisma com centro na posição x0, y0 ==> */

	/*for (i = 1; i < N; i++) { <== for dos pontos */
	i = ponto;

		gcalc = 0.0;

		/*for (j = 0; j < M; j++) { for dos prismas */
		j = prisma;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];
			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			for (k = 0; k < (nvertices-1); k++) { /* <== for dos vértices */

				/* Cálculo das coordenadas cartesianas dos vértices do prisma */
       			tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

                tetak = teta*(k+1);

				x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
				y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

				Xk1 = x[j][k] - Xp[i];
				Xk2 = x[j][k+1] - Xp[i];
				Yk1 = y[j][k] - Yp[i];
				Yk2 = y[j][k+1] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			} /* <== for dos vértices */

			/***** última aresta ==> *******/

			Xk1 = x[j][nvertices-1] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][nvertices-1] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== última aresta *******/

			gaux *= aux;

			gcalc += gaux;

		/*} for dos prismas */

	/*} <== for dos pontos */

   	/* <== Cálculo da atração gravitacional do prisma com centro na posição x0, y0 */

   	/* Cálculo da atração gravitacional do prisma com centro na posição x0 + dr, y0 ==> */

	/*for (i = 1; i < N; i++) { <== for dos pontos */
	i = ponto;

		gcalc_perturb = 0.0;

		/*for (j = 0; j < M; j++) { for dos prismas */
		j = prisma;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];
			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			for (k = 0; k < (nvertices-1); k++) { /* <== for dos vértices */

				/* Cálculo das coordenadas cartesianas dos vértices do prisma */
       			tetak = teta*k;

				x[j][k] = (x0[j] + dr) - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

                tetak = teta*(k+1);

				x[j][k+1] = (x0[j] + dr) - (raio[j][k+1]*sin(tetak));
				y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

				Xk1 = x[j][k] - Xp[i];
				Xk2 = x[j][k+1] - Xp[i];
				Yk1 = y[j][k] - Yp[i];
				Yk2 = y[j][k+1] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			} /* <== for dos vértices */

			/***** última aresta ==> *******/

			Xk1 = x[j][nvertices-1] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][nvertices-1] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== última aresta *******/

			gaux *= aux;

			gcalc_perturb += gaux;

		/*} for dos prismas */

	/*} <== for dos pontos */

   	/* <== Cálculo da atração gravitacional do prisma com centro na posição x0 + dr, y0 */

	derivada = (double)((gcalc_perturb - gcalc)/dr);

	return derivada;

}

double derivada_ajuste_y0 (int ponto, int prisma, int nvertices, double *Xp, double *Yp, double *Zp, 
					double *x0, double *y0, double *rho, double *z1, double *z2,
					double **raio, double teta, double **x, double **y, double dr) {

	int i, j, k;
	double gaux, aux, auxk1, auxk2, aux1k1, aux1k2, aux2k1, aux2k2;
	double Ak1, Ak2, Bk1, Bk2, Ck1, Ck2, Dk1, Dk2, E1k1, E1k2, E2k1, E2k2;
	double Xk1, Xk2, Yk1, Yk2, Z1, Z1_quadrado, Z2, Z2_quadrado, Qk1, Qk2;
	double R1k1, R1k2, R2k1, R2k2, p, p_quadrado, tetak;

	double gcalc, gcalc_perturb, derivada;

	/******** Cálculo da derivada em relação a coordenada y0 do centro de um prisma **********/

	/* Cálculo da atração gravitacional do prisma com centro na posição x0, y0 ==> */

	/*for (i = 1; i < N; i++) { <== for dos pontos */
	i = ponto;

		gcalc = 0.0;

		/*for (j = 0; j < M; j++) { for dos prismas */
		j = prisma;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];
			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			for (k = 0; k < (nvertices-1); k++) { /* <== for dos vértices */

				/* Cálculo das coordenadas cartesianas dos vértices do prisma */
       			tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = y0[j] + (raio[j][k]*cos(tetak));

                tetak = teta*(k+1);

				x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
				y[j][k+1] = y0[j] + (raio[j][k+1]*cos(tetak));

				Xk1 = x[j][k] - Xp[i];
				Xk2 = x[j][k+1] - Xp[i];
				Yk1 = y[j][k] - Yp[i];
				Yk2 = y[j][k+1] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			} /* <== for dos vértices */

			/***** última aresta ==> *******/

			Xk1 = x[j][nvertices-1] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][nvertices-1] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== última aresta *******/

			gaux *= aux;

			gcalc += gaux;

		/*} for dos prismas */

	/*} <== for dos pontos */

   	/* <== Cálculo da atração gravitacional do prisma com centro na posição x0, y0 */

   	/* Cálculo da atração gravitacional do prisma com centro na posição x0, y0 + dr ==> */

	/*for (i = 1; i < N; i++) { <== for dos pontos */
	i = ponto;

		gcalc_perturb = 0.0;

		/*for (j = 0; j < M; j++) { for dos prismas */
		j = prisma;

			aux = GRAV*rho[j]*1E12*1E8;
			/* o 1E12 transforma a densidade de g/cm^3 para kg/km^3 */
			/* o 1E8 transforma a aceleração gravitacional de kg/s^2 para mGal */

			gaux = 0.0;

			Z1 = z1[j] - Zp[i];
			Z2 = z2[j] - Zp[i];
			Z1_quadrado = pow(Z1, 2);
			Z2_quadrado = pow(Z2, 2);

			for (k = 0; k < (nvertices-1); k++) { /* <== for dos vértices */

				/* Cálculo das coordenadas cartesianas dos vértices do prisma */
       			tetak = teta*k;

				x[j][k] = x0[j] - (raio[j][k]*sin(tetak));
				y[j][k] = (y0[j] + dr) + (raio[j][k]*cos(tetak));

                tetak = teta*(k+1);

				x[j][k+1] = x0[j] - (raio[j][k+1]*sin(tetak));
				y[j][k+1] = (y0[j] + dr) + (raio[j][k+1]*cos(tetak));

				Xk1 = x[j][k] - Xp[i];
				Xk2 = x[j][k+1] - Xp[i];
				Yk1 = y[j][k] - Yp[i];
				Yk2 = y[j][k+1] - Yp[i];

				p = (Xk1*Yk2) - (Xk2*Yk1);
				p_quadrado = pow(p, 2);
				Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
				Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

				Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
				Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

				R1k1 = Ak1 + Z1_quadrado;
				R1k1 = pow(R1k1, 0.5);
				R1k2 = Ak2 + Z1_quadrado;
				R1k2 = pow(R1k2, 0.5);
				R2k1 = Ak1 + Z2_quadrado;
				R2k1 = pow(R2k1, 0.5);
				R2k2 = Ak2 + Z2_quadrado;
				R2k2 = pow(R2k2, 0.5);

				Ak1 = pow(Ak1, 0.5);
				Ak2 = pow(Ak2, 0.5);

				Bk1 = pow(Qk1, 2) + p_quadrado;
				Bk1 = pow(Bk1, 0.5);
				Bk2 = pow(Qk2, 2) + p_quadrado;
				Bk2 = pow(Bk2, 0.5);

				Ck1 = Qk1*Ak1;
				Ck2 = Qk2*Ak2;

				Dk1 = Divide_macro(p, 2.0);
				Dk2 = Dk1;
				Dk1 *= Divide_macro(Ak1, Bk1);
				Dk2 *= Divide_macro(Ak2, Bk2);

				E1k1 = R1k1*Bk1;
				E1k2 = R1k2*Bk2;
				E2k1 = R2k1*Bk1;
				E2k2 = R2k2*Bk2;

				auxk1 = Divide_macro(Qk1, p);
				auxk2 = Divide_macro(Qk2, p);
				aux1k1 = Divide_macro(Z1, R1k1);
				aux1k2 = Divide_macro(Z1, R1k2);
				aux2k1 = Divide_macro(Z2, R2k1);
				aux2k2 = Divide_macro(Z2, R2k2);

				gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
				gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
				gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

				gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
				gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			} /* <== for dos vértices */

			/***** última aresta ==> *******/

			Xk1 = x[j][nvertices-1] - Xp[i];
			Xk2 = x[j][0] - Xp[i];
			Yk1 = y[j][nvertices-1] - Yp[i];
			Yk2 = y[j][0] - Yp[i];

			p = (Xk1*Yk2) - (Xk2*Yk1);
			p_quadrado = pow(p, 2);
			Qk1 = ((Yk2 - Yk1)*Yk1) + ((Xk2 - Xk1)*Xk1);
			Qk2 = ((Yk2 - Yk1)*Yk2) + ((Xk2 - Xk1)*Xk2);

			Ak1 = pow(Xk1, 2) + pow(Yk1, 2);
			Ak2 = pow(Xk2, 2) + pow(Yk2, 2);

			R1k1 = Ak1 + Z1_quadrado;
			R1k1 = pow(R1k1, 0.5);
			R1k2 = Ak2 + Z1_quadrado;
			R1k2 = pow(R1k2, 0.5);
			R2k1 = Ak1 + Z2_quadrado;
			R2k1 = pow(R2k1, 0.5);
			R2k2 = Ak2 + Z2_quadrado;
			R2k2 = pow(R2k2, 0.5);

			Ak1 = pow(Ak1, 0.5);
			Ak2 = pow(Ak2, 0.5);

			Bk1 = pow(Qk1, 2) + p_quadrado;
			Bk1 = pow(Bk1, 0.5);
			Bk2 = pow(Qk2, 2) + p_quadrado;
			Bk2 = pow(Bk2, 0.5);

			Ck1 = Qk1*Ak1;
			Ck2 = Qk2*Ak2;

			Dk1 = Divide_macro(p, 2.0);
			Dk2 = Dk1;
			Dk1 *= Divide_macro(Ak1, Bk1);
			Dk2 *= Divide_macro(Ak2, Bk2);

			E1k1 = R1k1*Bk1;
			E1k2 = R1k2*Bk2;
			E2k1 = R2k1*Bk1;
			E2k2 = R2k2*Bk2;

			auxk1 = Divide_macro(Qk1, p);
			auxk2 = Divide_macro(Qk2, p);
			aux1k1 = Divide_macro(Z1, R1k1);
			aux1k2 = Divide_macro(Z1, R1k2);
			aux2k1 = Divide_macro(Z2, R2k1);
			aux2k2 = Divide_macro(Z2, R2k2);

			gaux += (Z2 - Z1)*(atan(auxk2) - atan(auxk1));
			gaux += Z2*(atan(aux2k1*auxk1) - atan(aux2k2*auxk2));
			gaux += Z1*(atan(aux1k2*auxk2) - atan(aux1k1*auxk1));

			gaux += Dk1*(log(Divide_macro2((E1k1 - Ck1), (E1k1 + Ck1))) - log(Divide_macro2((E2k1 - Ck1), (E2k1 + Ck1))));
			gaux += Dk2*(log(Divide_macro2((E2k2 - Ck2), (E2k2 + Ck2))) - log(Divide_macro2((E1k2 - Ck2), (E1k2 + Ck2))));

			/***** <== última aresta *******/

			gaux *= aux;

			gcalc_perturb += gaux;

		/*} for dos prismas */

	/*} <== for dos pontos */

   	/* <== Cálculo da atração gravitacional do prisma com centro na posição x0, y0 + dr */

	derivada = (double)((gcalc_perturb - gcalc)/dr);

	return derivada;

}

int paraview_error (int M, int nvertices, double teta, double **raio, double **raio_std, double *z1, double *x0, double *x0_std, double *y0, double *y0_std, int contador) {

	int i, j, k, l;
	double aux0, aux1, aux2;
	char str[100];

	FILE *arquivo;

	sprintf(str, "vert_fit_erro%d.vtk", contador);

	arquivo = fopen(str, "w");

	fprintf (arquivo, "# vtk DataFile Version 2.0\n");
	fprintf (arquivo, "Modelo Inverso\n");
	fprintf (arquivo, "ASCII\n");
	fprintf (arquivo, "DATASET POLYDATA\n");
	fprintf (arquivo, "POINTS %d float\n\n", (2*M*(nvertices+2)));

	for (i = 0; i < M; i++) {

		for (j = 0; j < nvertices; j++) {

            aux0 = teta*j;

			aux1 = x0[i] - ((raio[i][j] + raio_std[i][j])*sin(aux0));
			aux2 = y0[i] + ((raio[i][j] + raio_std[i][j])*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

			aux1 = x0[i] - ((raio[i][j] - raio_std[i][j])*sin(aux0));
			aux2 = y0[i] + ((raio[i][j] - raio_std[i][j])*cos(aux0));

			fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", aux1, aux2, z1[i]);

		}

		fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", (x0[i] + x0_std[i]), y0[i], z1[i]);
		fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", (x0[i] - x0_std[i]), y0[i], z1[i]);

		fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", x0[i], (y0[i] + y0_std[i]), z1[i]);
		fprintf (arquivo, "%15.5lf %15.5lf %15.5lf\n", x0[i], (y0[i] - y0_std[i]), z1[i]);

	}

    fprintf(arquivo, "\n");

	fprintf (arquivo, "LINES %d %d\n\n", (M*(nvertices+2)), (3*M*(nvertices+2)));

	for (i = 0; i < ((2*M*(nvertices+2)) - 1); i += 2) {

		fprintf (arquivo, "2 %7d %7d\n", i, (i+1));

	}
	
	fprintf (arquivo, "2 %7d %7d\n", i, (i+1));

	fclose(arquivo);
	
	return 0;

}

double fobj_ridge (double uridge, double **raio, int M, int nvertices) {

	int i, j;
	double ridge;

	ridge = 0.0;

	if (uridge != 0.0) {

		for (i = 0; i < M; i++) {
		
			for (j = 0; j < nvertices; j++) {

				ridge += raio[i][j]*raio[i][j];

			}

		}

		ridge *= uridge;

	}

	return ridge;
	
}

double fobj_absolut_equality_raios (double uae_raios, double **raio, double *pae, int nvertices) {

	int i;
	double ae_raios;

	ae_raios = 0.0;

	if (uae_raios != 0.0) {

		for (i = 0; i < nvertices; i++) {

			ae_raios += pow((raio[0][i] - pae[i]), 2);

		}

        ae_raios *= uae_raios;

	}

	return ae_raios;

}

double fobj_absolut_equality_origem (double uae_x0, double uae_y0, double x0_afloramento, double y0_afloramento, double *x0, double *y0, int M) {

	double ae_x0, ae_y0, ae_origem;
	
	ae_x0 = 0.0;
	ae_y0 = 0.0;
	ae_origem = 0.0;

	if ((uae_x0 != 0.0) && (uae_y0 != 0.0)) {

		ae_x0 = pow((x0[0] - x0_afloramento), 2);

        ae_x0 *= uae_x0;

		ae_y0 = pow((y0[0] - y0_afloramento), 2);

        ae_y0 *= uae_y0;

	}
	
	ae_origem = ae_x0 + ae_y0;

	return ae_origem;

}

double fobj_flatness_raios (double uflatness_rad, double uflatness_vert, double **raio, int M, int nvertices) {

	int i, j;
	double flatness_rad, flatness_vert, flatness_raios;

	flatness_rad = 0.0;
	
	if (uflatness_rad != 0.0) {

		for (i = 0; i < M ; i++) {

			for (j = 0; j < (nvertices - 1); j++) {

				flatness_rad += (raio[i][j] - raio[i][j + 1])*(raio[i][j] - raio[i][j + 1]);

			}

            flatness_rad += (raio[i][j] - raio[i][0])*(raio[i][j] - raio[i][0]);

		}


		flatness_rad *= uflatness_rad;

	}

	flatness_vert = 0.0;

	if ((uflatness_vert != 0.0) && (M != 1)) {

		for (i = 0; i < (M - 1); i++) {

			for (j = 0; j < nvertices; j++) {

				flatness_vert += pow((raio[i][j] - raio[i + 1][j]), 2);

			}

		}


		flatness_vert *= uflatness_vert;

	}

	flatness_raios = flatness_rad + flatness_vert;
	
	return flatness_raios;

}

double	fobj_flatness_origens (double uflatness_x0, double uflatness_y0, double *x0, double *y0, int M) {
	
	int i;
	double flatness_x0, flatness_y0, flatness_origens;

	flatness_x0 = 0.0;
	flatness_y0 = 0.0;

	if ((uflatness_x0 != 0.0) && (uflatness_y0 != 0) && (M != 1)) {

		for (i = 0; i < (M - 1); i++) {

			flatness_x0 += pow((x0[i] - x0[i + 1]), 2);

			flatness_y0 += pow((y0[i] - y0[i + 1]), 2);

		}

		flatness_x0 *= uflatness_x0;
		flatness_y0 *= uflatness_y0;

	}

	flatness_origens = flatness_x0 + flatness_y0;
	
	return flatness_origens;

}
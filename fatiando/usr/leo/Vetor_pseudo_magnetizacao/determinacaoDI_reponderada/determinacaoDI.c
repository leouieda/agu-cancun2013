#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M_PI		3.14159265358979323846

double **aloca_matriz_double (FILE *, int, int);

double **libera_matriz_double (int, double **);

double *aloca_vetor_double (FILE *, int);

double *libera_vetor_double (double *);

void sist_lin_LU (FILE *relatorio, int P, double **H, double *grad, double *p, double lambida);

void incoporacao_ajuste(int N, int M, double Hx, double Hy, double Hz, double **B, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, double *g, double *Tobs, double *T, double *grad, double **H, double *residuo);

void mod_direto(int N, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, double Hx, double Hy, double Hz, double F, double **B, double *p, double *dp, double *T);

double fobj_ajuste (int N, double *Tobs, double *T, double *residuo);

double fobj_ridge(int M, double *p, double *dp, double uridge);

void incorporacao_ridge(int M, double **H, double *grad, double *p, double *dp, double uridge);

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

int main () {

	int i, j;
	int N, M;
	int semente, iteracao, ITMAX, iteracao_marq, ITMAX_marq;
	double *gxx, *gxy, *gxz, *gyy, *gyz, *gzz, *Tobs, *T, **B;
	double *Y, *X, *Z;
	double *p, *p_anterior, *p_inicial, *dp, inclinacao, declinacao, intensidade;
	double *grad, **H, *g, *residuo;
	double lambida_inicial, lambida, dlambida, fobj1, fobj2, variacao_relativa, epsilon, tau, tau_min;
	double stdev_mag, stdev_tensor;
	double D, I, F, Hx, Hy, Hz;
	double uridge;
	double aux0, aux1, aux2, aux3;
	char str[100];

	FILE *relatorio, *entrada, *saida;

	relatorio = fopen ("relatorio.txt", "w");

	/* Descrição das variáveis ==>
	
	<== Descrição das variáveis */

	/* leitura do arquivo da aproximação inicial e parâmetros da inversão ==> */

		sprintf (str, "input.txt");

		if (fopen(str, "r") == NULL) {

			fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

			fclose (relatorio);

			printf ("Erro!\n\n");

			system ("PAUSE");

			return 0;

		}

		entrada = fopen(str, "r");

		fscanf(entrada, "%lf %lf %lf %lf", &lambida_inicial, &dlambida, &epsilon, &tau_min);
		fscanf(entrada, "%d %d", &ITMAX, &ITMAX_marq);
		fscanf(entrada, "%lf %lf %d", &stdev_mag, &stdev_tensor, &semente);

		idum = semente;
		
		/* esquentar o gerador de números aleatórios ==> */
		for (i = 0; i < 5000; i++) {
		
			aux0 = gasdev(&idum);
		
		}		
		/* <== esquentar o gerador de números aleatórios */
		
		fscanf (entrada, "%lf", &uridge);

		fscanf (entrada, "%lf %lf %lf", &D, &I, &F);
		
		D = (double)(3.14159265358979323846*D/180.0);
		I = (double)(3.14159265358979323846*I/180.0);

		/* componente horizontal do campo geomagnético */
		aux0 = F*cos(I);

		Hx = aux0*cos(D);
		Hy = aux0*sin(D);
		/* componente vertical do campo geomagnético */
		Hz = F*sin(I);
		
		M = 3;

		p_inicial = aloca_vetor_double (relatorio, M);
		
		for (i = 0; i < M; i++) {
		
			fscanf (entrada, "%lf", &p_inicial[i]);

		}

		aux0 = (double)(3.14159265358979323846*p_inicial[0]/180.0);
		aux1 = (double)(3.14159265358979323846*p_inicial[1]/180.0);
		p_inicial[2] *= ((double)(100/66.7));

		/* componente horizontal do vetor magnetização inicial */
		aux2 = p_inicial[2]*cos(aux1);

		p_inicial[0] = aux2*cos(aux0);
		p_inicial[1] = aux2*sin(aux0);
		
		/* componente vertical do vetor magnetização inicial */
		p_inicial[2] = p_inicial[2]*sin(aux1);

	/* <== leitura do arquivo da aproximação inicial e parâmetros da inversão */

	/******** leitura do arquivo de dados observados ==> **************/

		sprintf (str, "dados.txt");

		if (fopen(str, "r") == NULL) {

			fprintf (relatorio, "Arquivo %s nao encontrado!\n\n", str);

			fclose (relatorio);

			printf ("Erro!\n\n");

			system ("PAUSE");

			return 0;

		}

		entrada = fopen(str, "r");
		
		fscanf(entrada, "%d", &N);

		X = aloca_vetor_double(relatorio, N);
		Y = aloca_vetor_double(relatorio, N);
		Z = aloca_vetor_double(relatorio, N);
		Tobs = aloca_vetor_double(relatorio, N);
		gxx = aloca_vetor_double(relatorio, N);
		gxy = aloca_vetor_double(relatorio, N);
		gxz = aloca_vetor_double(relatorio, N);
		gyy = aloca_vetor_double(relatorio, N);
		gyz = aloca_vetor_double(relatorio, N);
		gzz = aloca_vetor_double(relatorio, N);

		for (i = 0; i < N; i++) {

			if (fscanf(entrada, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Y[i], &X[i], &Z[i], &Tobs[i], &gxx[i], &gxy[i], &gxz[i], &gyy[i], &gyz[i], &gzz[i]) != 10) {

				fprintf(relatorio, "Erro na leitura do arquivo %s!\n\n", str);

				fclose (relatorio);

				printf ("Erro!\n\n");

				system ("PAUSE");

				return 0;

			}
			
			Tobs[i] += stdev_mag*gasdev(&idum);
			
			gxx[i] += stdev_tensor*gasdev(&idum);
			gxy[i] += stdev_tensor*gasdev(&idum);
			gxz[i] += stdev_tensor*gasdev(&idum);
			gyy[i] += stdev_tensor*gasdev(&idum);
			gyz[i] += stdev_tensor*gasdev(&idum);
			gzz[i] += stdev_tensor*gasdev(&idum);

		}

		fclose (entrada);

	/******** <== leitura do arquivo de dados observados **************/

	/******** inversão ==> ********************************************/

	p = aloca_vetor_double (relatorio, M);
	p_anterior = aloca_vetor_double (relatorio, M);

	dp = aloca_vetor_double (relatorio, M);

	for (i = 0; i < M; i++) {
	
		p[i] = p_inicial[i];
		
	}

	T = aloca_vetor_double (relatorio, N);
	B = aloca_matriz_double (relatorio, N, M);

	residuo = aloca_vetor_double (relatorio, N);
	
	for (i = 0; i < N; i++) {
	
		residuo[i] = 1.0;
	
	}
	
	mod_direto(N, gxx, gxy, gxz, gyy, gyz, gzz, Hx, Hy, Hz, F, B, p, dp, T);

	saida = fopen ("ajuste_inicial.txt", "w");

	for (i = 0; i < N; i++) {

		fprintf (saida, "%15.3lf %15.3lf %15.3lf %15.5lf %15.5lf\n", Y[i], X[i], Z[i], Tobs[i], T[i]);

	}
	
	fclose (saida);

	/* função objetivo avaliada avaliada em p_inicial */
	fobj2 = fobj_ajuste(N, Tobs, T, residuo);
	fobj2 += fobj_ridge(M, p, dp, uridge);

	grad = aloca_vetor_double (relatorio, M);
	H = aloca_matriz_double (relatorio, M, M);
	
	g = aloca_vetor_double (relatorio, M);

	/* Esquema de reponderação ==> */
	do {
		
		iteracao = 0;

		for (i = 0; i < M; i++) {
		
			dp[i] = 0.0;
		
		}
		
		for (i = 0; i < M; i++) {
		
			p_anterior[i] = p[i];
		
		}
		
		lambida = lambida_inicial;

		fprintf (relatorio, "%5d %10.3E ---------- -----\n", iteracao, fobj2);
		printf ("\n%5d %10.3E ---------- -----\n", iteracao, fobj2);

		do {

			fobj1 = fobj2;

			for (i = 0; i < M; i++) {

				p_anterior[i] = p[i];

			}

			iteracao ++;

			/* atualização do p */
			for (i = 0; i < M; i++) {

				p[i] += dp[i];

			}

			/* calculo do gradiente e da hessiana ==> */
			for (i = 0; i < M; i++) {
			
				grad[i] = 0.0;
			
				for (j = 0; j < M; j++) {
				
					H[i][j] = 0.0;
				
				}		
			
			}

			incoporacao_ajuste(N, M, Hx, Hy, Hz, B, gxx, gxy, gxz, gyy, gyz, gzz, g, Tobs, T, grad, H, residuo);
			
			incorporacao_ridge(M, H, grad, p, dp, uridge);
			
			for (i = 0; i < M; i++) {
			
				for (j = (i+1); j < M; j++) {
				
					H[j][i] = H[i][j];
				
				}
			
			}
			/* <== calculo do gradiente e da hessiana */

			/* diminuir lambida porque este é aumentado na primeira iteração do looping */
			lambida = (double)(lambida/dlambida);

			iteracao_marq = 0;

			do {

				/* aumenta lambida Marquardt */
				lambida *= dlambida;

				if (lambida > 1E20) {
			
					lambida = 1E20;
		
				}

				/* resolve sistema */
				/* calcula p k+1 */
				sist_lin_LU (relatorio, M, H, grad, dp, lambida);

				mod_direto(N, gxx, gxy, gxz, gyy, gyz, gzz, Hx, Hy, Hz, F, B, p, dp, T);

				/* funcao objetiva avaliada em p k+1 */			
				fobj2 = fobj_ajuste(N, Tobs, T, residuo);
				fobj2 += fobj_ridge(M, p, dp, uridge);			
				
				iteracao_marq ++;

				fprintf (relatorio, "%5d %10.3E %10.3E %5d\n", iteracao, fobj2, lambida, iteracao_marq);
				printf ("%5d %10.3E %10.3E %5d\n", iteracao, fobj2, lambida, iteracao_marq);

			} while ((fobj2 >= fobj1) && (iteracao_marq <= ITMAX_marq));

			/* diminuir lambida */
			lambida = (double)(lambida/dlambida);
			
			if (lambida < 1E-20) {
			
				lambida = 1E-20;
			
			}
			
			variacao_relativa = (double)((fobj2 - fobj1)/fobj2);
			variacao_relativa = fabs(variacao_relativa);
			
			/*printf ("%10.3E\n", variacao_relativa);*/

		} while ((fobj2 < fobj1) && (variacao_relativa >= epsilon) && (iteracao <= ITMAX));

		if (fobj2 == fobj1) {

			printf ("\nfobj igual!\n");
			fprintf (relatorio, "\nfobj igual!\n");

		}

		if (fobj2 > fobj1) {

			printf ("\nfobj aumentou!\n");
			fprintf (relatorio, "\nfobj aumentou!\n");

		}

		if (variacao_relativa < epsilon) {

			printf ("\nvariacao relativa menor que epsilon!\n");
			fprintf (relatorio, "\nvariacao relativa menor que epsilon!\n");

		}

		if (fobj2 > fobj1) {

			for (i = 0; i < M; i++) {

				dp[i] = 0.0;
				
			}

			mod_direto(N, gxx, gxy, gxz, gyy, gyz, gzz, Hx, Hy, Hz, F, B, p, dp, T);

			/* funcao objetiva avaliada em p k+1 */			
			fobj2 = fobj_ajuste(N, Tobs, T, residuo);
			fobj2 += fobj_ridge(M, p, dp, uridge);

		}
		
		for (i = 0; i < N; i++) {
		
			residuo[i] = fabs(Tobs[i] - T[i]) + 1E-10;
		
		}

		aux0 = 0.0;
		aux1 = 0.0;

		for (i = 0; i < M; i++) {
		
			aux0 += pow((p_anterior[i] - p[i]), 2);
			aux1 += pow(p[i], 2);
		
		}
		
		aux0 = pow(aux0, 0.5);
		aux1 = pow(aux1, 0.5);

		tau = (double)(aux0/(1.0 + aux1));

		fprintf (relatorio, "\n%10.3E %10.3E\n\n", tau, tau_min);
		printf ("\n%10.3E %10.3E\n\n", tau, tau_min);
		
	} while (tau >= tau_min);
	/* <== Esquema de reponderação */

		if (fobj2 > fobj1) {

			for (i = 0; i < M; i++) {

				dp[i] = 0.0;
				
				p[i] = p_anterior[i];
				
			}

			mod_direto(N, gxx, gxy, gxz, gyy, gyz, gzz, Hx, Hy, Hz, F, B, p, dp, T);

			for (i = 0; i < N; i++) {
		
				residuo[i] = fabs(Tobs[i] - T[i]) + 1E-10;
		
			}

			/* funcao objetiva avaliada em p k+1 */			
			fobj2 = fobj_ajuste(N, Tobs, T, residuo);
			fobj2 += fobj_ridge(M, p, dp, uridge);

		}

	/******** <== inversão ********************************************/

	/* impressão dos arquivos de saída ==> */

		sprintf(str, "ajuste.txt");

		saida = fopen (str, "w");

		for (i = 0; i < N; i++) {

			fprintf (saida, "%15.3lf %15.3lf %15.3lf %15.5lf %15.5lf\n", Y[i], X[i], Z[i], Tobs[i], T[i]);

		}

		fprintf (saida, "\nfobj %15.5lf\n\n", fobj2);

		fclose (saida);
			
		aux0 = pow(p[0], 2);
		aux1 = pow(p[1], 2);
		aux2 = pow(p[2], 2);
		aux3 = pow((aux0 + aux1), 0.5);

		intensidade = pow((aux0 + aux1 + aux2), 0.5);
		inclinacao = atan2(p[2], aux3);
		declinacao = acos(((double)(p[0]/aux3)));

		if (p[1] < 0) {
		
			declinacao = -declinacao;
		
		}

		inclinacao = (double)(180.0*inclinacao/3.14159265358979323846);
		declinacao = (double)(180.0*declinacao/3.14159265358979323846);
		
		sprintf(str, "model.txt");

		saida = fopen (str, "w");

		fprintf (saida, "        p_inicial               p\n");

		for (i = 0; i < M; i++) {
		
			fprintf (saida, "%d %15.3lf %15.3lf\n", i, p_inicial[i], p[i]);

		}

		fprintf (saida, "\ndeclinacao  %15.3lf\n", declinacao);
		fprintf (saida, "inclinacao  %15.3lf\n", inclinacao);
		fprintf (saida, "intensidade %15.3lf\n", intensidade);
		
		fclose (saida);
		
		saida = fopen ("reducao_polo.txt", "w");		
		
		p[0] = 0.0;
		p[1] = 0.0;
		p[2] = intensidade;

		for (i = 0; i < M; i++) {

			dp[i] = 0.0;
			
		}
		
		Hx = 0.0;
		Hy = 0.0;
		Hz = F;

		mod_direto(N, gxx, gxy, gxz, gyy, gyz, gzz, Hx, Hy, Hz, F, B, p, dp, T);

		for (i = 0; i < N; i++) {

			fprintf (saida, "%15.3lf %15.3lf %15.3lf %15.5lf\n", Y[i], X[i], Z[i], T[i]);

		}

		fclose (saida);

	/* <== impressão dos arquivos de saída */

	printf ("\n");

	system ("PAUSE");

	return 0;

}

double **aloca_matriz_double (FILE *arq, int linha, int coluna) {

    double **m;  /* ponteiro para a matriz */
    int   i, j;

    /* aloca as linhas da matriz */

    m = (double **)calloc(linha, sizeof(double *));

    if (m == NULL) {

        fprintf (arq, "Memoria Insuficiente (linhas)!\n\n");

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

void sist_lin_LU (FILE *relatorio, int P, double **H, double *grad, double *p, double lambida) {

	int i, j, k, l;
	double *b, **L;

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

	for (j = P-1; j >= 0; j--) {

		p[j] = 0;

		for (k = j+1; k < P; k++) {

			p[j] -= L[j][k]*p[k];

		}

		p[j] += b[j];

		p[j] = (double)(p[j]/L[j][j]);

	}

	/************ <== resolução do sistema linear Hdraio = -grad **************/

	b = libera_vetor_double(b);
	L = libera_matriz_double (P, L);

	/************ <== resolução do sistema linear Hdraio = -grad **************/

}

void incoporacao_ajuste(int N, int M, double Hx, double Hy, double Hz, double **B, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, double *g, double *Tobs, double *T, double *grad, double **H, double *residuo) {

	int i, j, k;
	double aux0, aux1, aux2, aux3;

	for (i = 0; i < N; i++) { /* for i */

		/* i-ésima linha de G ==> */
		aux0 = (Hx + B[i][0]);
		aux1 = (Hy + B[i][1]);
		aux2 = (Hz + B[i][2]);
		
		aux3 = pow(aux0, 2) + pow(aux1, 2) + pow(aux2, 2);
		aux3 = pow(aux3, 0.5);
		
		g[0] = (double)(((gxx[i]*aux0) + (gxy[i]*aux1) + (gxz[i]*aux2))/aux3);
		g[1] = (double)(((gxy[i]*aux0) + (gyy[i]*aux1) + (gyz[i]*aux2))/aux3);
		g[2] = (double)(((gxz[i]*aux0) + (gyz[i]*aux1) + (gzz[i]*aux2))/aux3);
		/* <== i-ésima linha de G */

		for (j = 0; j < M; j++) { /* for j */
		
			grad[j] -= g[j]*((double)((Tobs[i] - T[i])/residuo[i]));
			
			for (k = 0; k <= j; k++) { /* for k */
			
				H[k][j] += ((double)((g[k]*g[j])/residuo[i]));
			
			} /* for k */
				
		} /* for j */

	} /* for i */
	
	for (aux0 = 0.0, i = 0; i < M; i++) {
	
		aux0 += H[i][i];
	
	}
	
	aux0 = (double)(M/aux0);

	for (i = 0; i < M; i++) {
	
		grad[i] *= aux0;
		
		for (j = 0; j <= i; j++) {
		
			H[j][i] *= aux0;
		
		}
	
	}

}

void mod_direto(int N, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, double Hx, double Hy, double Hz, double F, double **B, double *p, double *dp, double *T) {

	int i;
	double aux0; 

	for (i = 0; i < N; i++) {
	
		B[i][0] = ((p[0] + dp[0])*gxx[i]) + ((p[1] + dp[1])*gxy[i]) + ((p[2] + dp[2])*gxz[i]);
		B[i][1] = ((p[0] + dp[0])*gxy[i]) + ((p[1] + dp[1])*gyy[i]) + ((p[2] + dp[2])*gyz[i]);
		B[i][2] = ((p[0] + dp[0])*gxz[i]) + ((p[1] + dp[1])*gyz[i]) + ((p[2] + dp[2])*gzz[i]);

		aux0 = pow((Hx + B[i][0]), 2) + pow((Hy + B[i][1]), 2) + pow((Hz + B[i][2]), 2);
		
		T[i] = pow(aux0, 0.5) - F;
	
	}

}

double fobj_ajuste (int N, double *Tobs, double *T, double *residuo) {

	int i;
	double resultado;

	resultado = 0.0;

	for (i = 0; i < N; i++) {

		/*resultado += (double)(pow((Tobs[i] - T[i]), 2)/residuo[i]);*/
		/*resultado += (double)(fabs(Tobs[i] - T[i])/residuo[i]);*/
		resultado += fabs(Tobs[i] - T[i]);

	}
	
	resultado = (double)(resultado/N);

	return resultado;

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

double fobj_ridge(int M, double *p, double *dp, double uridge) {

	int i;
	double ridge;
	
	ridge = 0.0;

	for (i = 0; i < M; i++) {
	
		ridge += pow((p[i] + dp[i]), 2);

	}
	
	ridge *= uridge;

	return ridge;

}

void incorporacao_ridge(int M, double **H, double *grad, double *p, double *dp, double uridge) {

	int i;
	double aux0, aux1;
	
	for (i = 0; i < M; i++) {

		/* Gradiente */	
		grad[i] += uridge*(p[i] + dp[i]);

		/* Hessiana */
		H[i][i] += uridge;

	}

}
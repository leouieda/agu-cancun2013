#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 681

void main () {

	int i;
	double gxx[N], gxy[N], gxz[N], gyy[N], gyz[N], gzz[N];
	double Bx[N], By[N], Bz[N], T[N];
	double inclinacao, declinacao, intensidade, I, D, F;
	double Mx, My, Mz, Hx, Hy, Hz;
	double Y[N], X[N], Z[N];
	double G, rho;
	double a;
	double aux0, aux1, aux2;
	
	FILE *entrada, *saida;;

	G = 6.67E-20;
	G *= (1E12*1E9);
	rho = 1.0;

	F = 23000.0;
	D = 20.0;
	I = 5.0;
	
	intensidade = 2.0;
	intensidade *= 1E-7;
	intensidade *= ((double)(1.0/(G*rho)));
	declinacao = 0.0;
	inclinacao = -30.0;

	a = (double)(intensidade/(G*rho));

	inclinacao = (double)(3.14159265358979323846*inclinacao/180.0);
	declinacao = (double)(3.14159265358979323846*declinacao/180.0);

	/* componente horizontal do vetor magnetização */
	aux0 = intensidade*cos(inclinacao);
	Mx = aux0*cos(declinacao);
	My = aux0*sin(declinacao);

	/* componente vertical do campo geomagnético */
	Mz = intensidade*sin(inclinacao);

	I = (double)(3.14159265358979323846*I/180.0);
	D = (double)(3.14159265358979323846*D/180.0);

	/* componente horizontal do campo geomagnético */
	aux0 = F*cos(I);
	Hx = aux0*cos(D);
	Hy = aux0*sin(D);

	/* componente vertical do vetor magnetização */
	Hz = F*sin(I);

	entrada = fopen ("tensor_obs.txt", "r");

	for (i = 0; i < N; i++) {
	
		fscanf (entrada, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &Y[i], &X[i], &Z[i], &gxx[i], &gxy[i], &gxz[i], &gyy[i], &gyz[i], &gzz[i]);
	
	}
	
	fclose (entrada);
	
	saida = fopen("saida.txt", "w");
	
	for (i = 0; i < N; i++) {
	
		/*Bx[i] = ((double)(1.0/(G*rho)))*((gxx[i]*Mx)+(gxy[i]*My)+(gxz[i]*Mz));
		By[i] = ((double)(1.0/(G*rho)))*((gxy[i]*Mx)+(gyy[i]*My)+(gyz[i]*Mz));
		Bz[i] = ((double)(1.0/(G*rho)))*((gxz[i]*Mx)+(gyz[i]*My)+(gzz[i]*Mz));*/

		Bx[i] = ((gxx[i]*Mx)+(gxy[i]*My)+(gxz[i]*Mz));
		By[i] = ((gxy[i]*Mx)+(gyy[i]*My)+(gyz[i]*Mz));
		Bz[i] = ((gxz[i]*Mx)+(gyz[i]*My)+(gzz[i]*Mz));
		
		Bx[i] *= 1E9;
		By[i] *= 1E9;
		Bz[i] *= 1E9;

		aux0 = pow((Hx + Bx[i]), 2);
		aux1 = pow((Hy + By[i]), 2);
		aux2 = pow((Hz + Bz[i]), 2);

		T[i] = pow((aux0 + aux1 + aux2), 0.5) - F;
		
		fprintf (saida, "%15.5lf %15.5lf %15.5lf %15.5lf\n", Y[i], X[i], Z[i], T[i]);

	}

	fclose (saida);
	
	system ("PAUSE"); 

}
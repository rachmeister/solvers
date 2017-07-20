#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-14

double X, dx;
double A;
char fname[100];

double yprime(double x, double y) {
	return A*y;
}
double absval(double x) {
	if( x < 0. ) return -x;
	else return x;
}
double newton(double xn, double yn) {
	
	double y_new=yn, y_old;
	int cnt=0;

	do {
		y_old = y_new;
		y_new = y_old - ((1.-A*dx)*y_old-yn)/(1-A*dx);
		cnt++;
		if( y_new == NAN || cnt > 10e6 ) {
			printf("\tERROR: Converge failed -- please try again.\n");
			return NAN;
		}
	} while(absval((y_new-y_old)/y_old) > EPS);

	return y_new;
}

void exact(double x0, double y0) {
	double x = x0+dx;
	double y = y0;
	FILE * fp = fopen("exact.csv","w");

	//printf("Exact\t|%5s%lf%4s","",y0,"");
	fprintf(fp,"%.3lf\n%.17g\n", x0, y0);
	while( x <= X ) {
		y = y0*exp(A*x);
		x+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x, y);
	}
	//printf("\t|%5s%lf\n","",y);
	fclose(fp);
	
}
void eulers(double x0, double y0) {

	double x = x0;
	double y = y0;
	FILE *fp;
	
	sprintf(fname,"eulers_%.3f.csv",dx);
	fp=fopen(fname,"w");

	//printf("Eulers\t|%5s%lf%6s","",y0,"");
	fprintf(fp,"%.3lf\n%.17g\n", x, y);

	while( x <= X ) {
		y = y + dx*yprime(x,y);
		x+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x, y);
	}
	//printf("\t|%5s%lf\n","",y);
	fclose(fp);
}
void improved_eulers(double x0, double y0) {
	double x = x0;
	double yp, y = y0;
	FILE *fp;
	
	sprintf(fname,"improved_eulers_%.3f.csv",dx);
	//sprintf(fname,"improved_eulers.csv");
	fp=fopen(fname,"w");

	fprintf(fp,"%.3lf\n%.17g\n", x, y);

	while( x <= X) {
		//yp = yprime(x,y);
		//y = y + dx/2.*(yprime(x+dx,y+dx*yp)+yp);
		y = ((2.+A*dx)/(2.-A*dx))*y;
		//y = pow(((2.+A*dx)/(2.-A*dx)),(int(x/dx)));
		x+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x, y);
	}
	fclose(fp);
}
void implicit_eulers(double x0, double y0) {
	double x = x0;
	double yp, y = y0;
	FILE *fp;

	sprintf(fname,"implicit_eulers_%.3f.csv",dx);
	fp=fopen(fname,"w");

	//printf("Eulers3\t|%5s%lf%6s","",y,"");
	fprintf(fp,"%.3lf\n%.17g\n", x, y);
	while( x <= X) {
		y = newton(x,y);
		x+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x, y);
	}
	//printf("\t|%5s%lf\n","",y);
	fclose(fp);

}
void leapfrog(double x0, double y0) {
	int i;
	double x[3] = {x0,dx,0.};
	double y[3] = {y0,y0*exp(A*(x0+dx)),0.};
	FILE *fp;

	sprintf(fname,"leapfrog_%.3f.csv",dx);
	fp=fopen(fname,"w");

	//printf("Lpfrog\t|%5s%lf%6s","",y0,"");
	fprintf(fp,"%.3lf\n%.17g\n", x[0], y[0]);
	fprintf(fp,"%.3lf\n%.17g\n", x[1], y[1]);
	
	while( x[1] <= X ) {
		y[2] = y[0] + 2*dx*yprime(x[1],y[1]);
		y[0] = y[1]; y[1] = y[2];

		for(i=0;i<=2;i++) x[i]+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x[1], y[1]);
	}
	//printf("\t|%5s%lf\n","",y[1]);
	fclose(fp);
}
void rk2(double x0, double y0) {
	double x = x0;
	double y = y0;
	FILE *fp;

	sprintf(fname,"rk2_%.3f.csv",dx);
	fp=fopen(fname,"w");

	fprintf(fp,"%.3lf\n%.17g\n", x0, y0);

	while( x <= X ) {
		y = y + dx*yprime(x+0.5*dx, y+0.5*dx*yprime(x,y));

		x+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x, y);
	}
	fclose(fp);
}
void milnes(double x0, double y0) {
	int i;
	double x[3] = {x0,dx,0.};
	double y[3] = {y0,y0*exp(A*(x0+dx)),0.};
	FILE *fp;

	sprintf(fname,"milnes_%.3f.csv",dx);
	//sprintf(fname,"milnes.csv");
	fp=fopen(fname,"w");

	//printf("Milne\t|%5s%lf%6s","",y0,"");
	fprintf(fp,"%.3lf\n%.17g\n", x[0], y[0]);
	fprintf(fp,"%.3lf\n%.17g\n", x[1], y[1]);
	
	while( x[1] <= X ) {
		y[2] = 3./(3.-A*dx)*(4./3.*A*dx*y[1]+(1.+A*dx/3.)*y[0]);
		y[0] = y[1]; y[1] = y[2];

		for(i=0;i<=2;i++) x[i]+=dx;
		fprintf(fp,"%.3lf\n%.17g\n", x[1], y[1]);
	}
	//printf("\t|%5s%lf\n","",y[1]);
	fclose(fp);

}

int main(int argc, char *argv[]) {

	double x;
	double y0;
	int method;
	char *methods[]={"Exact","Eulerian","Improved Eulerian","Implicit Eulerian","Leapfrog","RK2","Milne"};

	X = 10.0;
	dx = 0.1;
	y0 = 1.0;
	A = 1.;

	if(argc == 5) {
		sscanf(argv[1],"%d",&method);
		sscanf(argv[2],"%lf",&X);
		sscanf(argv[3],"%lf",&dx);
		sscanf(argv[4],"%lf",&A);
	}

	printf("\n=================== Ordinary Differential Equation Solver ======================\n");

	printf("Equation: y'=Ay with y(0)=%.2lf\n",y0);
	printf("Over [0,%.2lf] with dx = %.2lf\n",X,dx);
	printf("Using method %s\n", methods[method]);

	//printf("x value\t|");
	//for(x=0;x<=X;x+=X) printf("%8s%.1lf%12s|","",x,"");
	//printf("\n");
	//printf("-----------------------------------------------------------------\n");

	switch (method) {
	case 0 :
		exact(0,y0);
		break;
	case 1 :
		eulers(0,y0);
		break;
	case 2 :
		improved_eulers(0,y0);
		break;
	case 3 :
		implicit_eulers(0,y0);
		break;
	case 4 :
		leapfrog(0,y0);
		break;
	case 5 :
		rk2(0,y0);
		break;
	case 6 :
		milnes(0,y0);
		break;
	}

	return 0;

}

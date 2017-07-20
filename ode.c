#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double X, dx;

double yprime(double x, double y) {
	return x+y;
}
void exact(double x0, double y0) {
	double x = x0+dx;
	double y = y0;

	printf("Exact\t|%5s%lf%6s","",y0,"");
	while( x <= X ) {
		y = exp(x) - x - 1;
		x+=dx;
		printf("%6s%lf%6s","",y,"");
	}
	printf("\n");
	
}
void eulers(double x0, double y0) {

	double x = x0;
	double y = y0;

	printf("Eulers\t|%5s%lf%6s","",y0,"");

	while( x < X ) {
		y = y + dx*yprime(x,y);
		x+=dx;
		printf("%6s%lf%6s","",y,"");
	}

	printf("\n");
}
void rk4(double x0, double y0) {
	double x = x0;
	double y = y0;

	double k1, k2, k3, k4;

	printf("RK-4\t|%5s%lf%6s","",y0,"");

	while( x < X) {
		k1 = dx*yprime(x,y);
		k2 = dx*yprime(x+dx/2.,y+k1/2.);
		k3 = dx*yprime(x+dx/2.,y+k2/2.);
		k4 = dx*yprime(x+dx,y+k3);

		y = y + (k1 + 2*k2 + 2*k3 + k4)/6.;
		x += dx;

		printf("%6s%lf%6s","",y,"");
	}

	printf("\n");
}
void adamsBashforth2step(double x0, double y0) {
	double x = x0+dx;
	double y[3];

	y[0] = y0;
	y[1] = x*x/2.;

	printf("AB-2\t|%5s%lf%6s","",y[0],"");
	printf("%6s%lf%6s","",y[1],"");

	while(x < X) {
		y[2] = y[1] + dx/2.*(3*yprime(x,y[1])-yprime(x-dx,y[0]));

		y[0] = y[1];
		y[1] = y[2];

		x += dx;

		printf("%6s%lf%6s","",y[1],"");
	}
	printf("\n");
}
void adamsMoulton2step(double x0, double y0) {
	double x = x0+dx;
	double y[3];

	y[0] = y0;
	y[1] = x*x/2.;
	y[2] = (x+dx)*(x+dx)/2.;

	printf("AM-2\t|%5s%lf%6s","",y[0],"");
	printf("%6s%lf%6s","",y[1],"");

	while(x < X) {
		y[2] = y[1] +dx/12.*(5.*yprime(x+dx,y[2]) + 8*yprime(x,y[1]) - yprime(x-dx,y[0]));
		
		y[0] = y[1];
		y[1] = y[2];
		
		x+=dx;
		printf("%6s%lf%6s","",y[1],"");
	}
	printf("\n");

}
void adamsMoulton2step_ex(double x0, double y0) {
	double x = x0+dx;
	double y[3];

	y[0] = y0;
	y[1] = x*x/2.;

	printf("AM-2\t|%5s%lf%6s","",y[0],"");
	printf("%6s%lf%6s","",y[1],"");

	while(x < X) {
		y[2] = (y[0]+128.*y[1]+6.*0.1+12.*x)/115.;
		
		y[0] = y[1];
		y[1] = y[2];

		x+=dx;
		printf("%6s%lf%6s","",y[1],"");
	}
	printf("\n");


}
void improved_eulers(double x0, double y0) {
	double x = x0;
	double y = y0;

	printf("Eulers2\t|%5s%lf%6s","",y,"");

	while( x < X) {
		y = y + dx*yprime(x+0.5*dx,y+0.5*dx*yprime(x,y));

		x+=dx;
		printf("%6s%lf%6s","",y,"");
	}
	printf("\n");

}
void preCorr(double x0, double y0) {
	double x = x0;
	double y[6];

	y[0] = y0;
	y[1] = (x0+dx)*(x0+dx)/2.;
	y[2] = (x0+2*dx)*(x0+2*dx)/2.;
	y[3] = (x0+3*dx)*(x0+3*dx)/2.;
	y[4] = (x0+4*dx)*(x0+4*dx)/2.;

	x = x0+4*dx;

	y[5] = y[4] + dx/720.*(1901.*yprime(x,y[4])-2774.*yprime(x-dx,y[3])+2616.*yprime(x-2*dx,y[2])-1274.*yprime(x-3*dx,y[1])+251.*yprime(x-4*dx,y[0]));

	y[5] = y[4] + dx/24.*(9.*yprime(x+dx,y[5])+19.*yprime(x,y[4])-5.*yprime(x-dx,y[3])+yprime(x-2*dx,y[2]));

	printf("PreCorr\t|%5s%lf%6s","",y[0],"");
	x = x0;
	while( x < X) {
		x += dx;
		printf("%6s%lf%6s","",y[(int)(x/dx)],"");
	}
	printf("\n");
}

int main(int argc, char *argv[]) {

	double x;

	X = 0.5;
	dx = 0.1;

	printf("\n=================== Ordinary Differential Equation Solver ======================\n");

	printf("Equation: y'=x+y with y(0)=0\n");
	printf("Over [0,0.5] with dx = 0.1\n\n\n");

	printf("x value\t|");
	for(x=0;x<=X;x+=dx) printf("%8s%.1lf%8s|","",x,"");
	printf("\n");
	printf("---------------------------------------------------------------------------------------------------------------------------------\n");

	exact(0,0);
	eulers(0,0);
	rk4(0,0);
	adamsBashforth2step(0,0);
	//adamsMoulton2step(0,0);
	adamsMoulton2step_ex(0,0);
	improved_eulers(0,0);
	preCorr(0,0);

	return 0;

}

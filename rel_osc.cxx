// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
// edited by: markus, date: 2015-12-09
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 2;
	double dx = 0.05,x=0;
	const double L = 20;
  double y0[dim] = {2.0, 0};
	double yn[dim];
	double k1[dim], k2[dim], k3[dim], k4[dim];
	double p0 = 0;

	for(int i = 1; i<=10; i++){
	x=0;
	p0 += 0.5;
	y0[0] = p0;
	y0[1] = 0;

  //out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx, k1, k2, k3, k4);

		if(yn[1]<0 && y0[1]>0) break;

	for(int i=0; i<dim; i++) y0[i] = yn[i];
		//out << x << "\t" << y0[0] << "\t" << y0[1] << endl;


	}

	double b[3];
	double theta;
	double l = 0.0;
	double r = 1.0;
	double w = y0[1];

	while(abs(y0[1])>1E-8){
	theta= (r+l)/2.0;

	b[0] = theta - 3.0*pow(theta, 2.0)/2.0 + 2.0*pow(theta,3.0)/3.0;
	b[1] = pow(theta, 2.0) - 2.0*pow(theta,3.0)/3.0;
	b[2] = - pow(theta, 2.0)/2.0 + 2.0*pow(theta,3.0)/3.0;

	y0[1] = w + dx*b[0]*k1[1] + dx*b[1]*k2[1] + dx*b[1]*k3[1] + dx*b[2]*k4[1];

	if(y0[1]<0) r=theta;
	if(y0[1]>0) l=theta;

	}
	cout << p0 << "\t" << x+theta*dx << endl;
	//cout << x+theta*dx << endl;
	}
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;

	for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
	f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

	for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// relativistic motion of a particle in a strong electrostatic field E(x)
void f(double* const y0, const double x)
{
	double y[2] = { y0[0], y0[1]};

  	y0[0] = y[1];
	y0[1] = -y[0]/sqrt(1+y[0]*y[0]);
}

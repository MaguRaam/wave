#include "plot.h"

//wave 1D 2nd order stencil:
BZ_DECLARE_STENCIL4(wave12, P1, P2, P3, CFL)
P3(0, 0) = CFL * CFL * (P2(1, 0) - 2.0 * P2(0, 0) + P2(-1, 0)) + 2.0 * P2(0, 0) - P1(0, 0);
BZ_END_STENCIL


void wave(double tf = 2)
{
	cout.precision(6);
	auto print = [](const auto& x){std::cout<<x<<std::endl;};
	
	//grid:
	int N = 100;
	double xo = -1.0, xl = 1.0, h = (xl - xo)/double(N-1); 
	Array<double,1> x(N);
	x = xo + tensor::i*h;
	
	//wave speed:
 	double c = 1.0; 
	
	//cfl:
	double cfl = 0.5;
	
	//time-step:
	double dt = cfl*h/c;
	
	//initialize arrays:
	Array<double,1> P1,P2,P3,CFL;
	allocateArrays(shape(N), P1, P2, P3, CFL);
	
	//set cfl:
	CFL = cfl;

	//initial shape of the wave po(x):
	double sigma = 0.2;
	auto po = [sigma](double x){return exp( (-0.5*x*x)/(sigma*sigma) );};
	
	//exact solution
	double t;
	auto pexact = [&t,po,c](auto x){return 0.5*po(c*t - x) + 0.5*po(c*t + x);}; 
	
	
	//ics at t = 0:
	int nt = 0;
	t = 0.0;
	std::transform(x.begin(),x.end(),P1.begin(),pexact); 
	write(P1,x,t,nt);
	
	//ics at t = dt:
	nt++;
	t+=dt;
	std::transform(x.begin(),x.end(),P2.begin(),pexact); 
	write(P2,x,t,nt);
	
	while (t < tf){
		nt++;
		t+=dt;
		
		//update pressure:
		applyStencil(wave12(), P1, P2, P3, CFL);
		
		if (nt%100 == 0){
			write(P3,x,t,nt);
		}
		print(t);
		 
		cycleArrays(P1,P2,P3);
	}
	
}

int main()
{
	wave();
	return 0;
}


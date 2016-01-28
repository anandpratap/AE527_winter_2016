#include <cmath>
#include <cstdio>
#include "stdio.h"
#include "time.h"
#include <iostream>
#include "omp.h"

class VortexSolver{
public:
	int nx, ny, n;
	double Lx, Ly;
	double *x, *y, *alpha, *omega, *u, *v;
	double *k1u, *k1v;
	double *k2u, *k2v;
	double *k3u, *k3v;
	double *k4u, *k4v;
	double *xrk, *yrk;
	double dx, dy;
	double t, dt;
	int write_freq, integrator;
	VortexSolver(int, int, double);
	~VortexSolver(){
		delete[] x;
		delete[] y;
		delete[] alpha;
		delete[] omega;
		delete[] u;
		delete[] v;
		delete[] k1u;
		delete[] k1v;
		delete[] k2u;
		delete[] k2v;
		delete[] k3u;
		delete[] k3v;
		delete[] k4u;
		delete[] k4v;
		delete[] xrk;
		delete[] yrk;
	};
	
	void initialize(void);
	void calc_velocity(double*, double*, double*, double*, double*, double*, double*);
	void integrate(void);
	void run(double);
	void write_solution(int step);
};

VortexSolver::VortexSolver(int nxi, int nyi, double dti){
	nx = nxi;
	ny = nyi;
	dt = dti;
	
	n = nx*ny;
	
	x = new double[n];
	y = new double[n];

	xrk = new double[n];
	yrk = new double[n];

	alpha = new double[n];
	omega = new double[n];
	u = new double[n];
	v = new double[n];

	k1u = new double[n];
	k1v = new double[n];
	k2u = new double[n];
	k2v = new double[n];
	k3u = new double[n];
	k3v = new double[n];
	k4u = new double[n];
	k4v = new double[n];

	write_freq = int(1.0/dt);
	integrator = 4;
	initialize();
	
}

void VortexSolver::initialize(){
	Lx = 10.0;
	Ly = 5.0;
	dx = 2*Lx/(nx - 1);
	dy = 2*Ly/(ny - 1);
	double A = dx*dy;

	for(int i=0; i<n; i++){
		x[i] = -Lx + dx*(i%nx);
		y[i] = -Ly + dy*(i/nx);
	}
	
	for(int i=0; i<n; i++){
		alpha[i] = 0.0;
		omega[i] = 0.0;
		u[i] = 0.0;
		v[i] = 0.0;
	}

	//
	double tau, rc2;
	double xc, yc;

		
	auto vorticity = [](auto tau, auto rc2, auto r2){return tau/M_PI/rc2*exp(-r2/rc2);};


	tau = 1.0;
	rc2 = pow(1.0, 2);
	xc = -5.0;
	yc = 0.0;
	for(int i=0; i<n; i++){
		double r2 = pow(x[i] - xc, 2) + pow(y[i] - yc, 2);
		omega[i] += vorticity(tau, rc2, r2);
	}
	
	xc = 5.0;
	yc = 0.0;
	for(int i=0; i<n; i++){
		double r2 = pow(x[i] - xc, 2) + pow(y[i] - yc, 2);
		omega[i] -= vorticity(tau, rc2, r2);
	}


	tau = 0.1;
	rc2 = pow(0.5, 2);
	xc = -3.0;
	yc = 0.0;
	for(int i=0; i<n; i++){
		double r2 = pow(x[i] - xc, 2) + pow(y[i] - yc, 2);
		omega[i] -= vorticity(tau, rc2, r2);
	}
	
	xc = 3.0;
	yc = 0.0;
	for(int i=0; i<n; i++){
		double r2 = pow(x[i] - xc, 2) + pow(y[i] - yc, 2);
		omega[i] += vorticity(tau, rc2, r2);
	}

	for(int i=0; i<n; i++){
		alpha[i] = omega[i]*A;
	}
}

void VortexSolver::write_solution(int step){
	char filename[100];
	sprintf(filename, "data/step_%06d.dat", step);
	printf("Writing file %s\n", filename);

	FILE *fp;
	if ((fp=fopen(filename,"wt"))==NULL ){
		printf("Error: cannot open file %s\n", filename);
		throw(-1);
	}
	fprintf(fp, "Variables=\"X\",\"Y\",\"OMEGA\",\"ALPHA\", \"U\", \"V\"\nZone I=     %d,J=     %d, F=POINT\n", nx, ny);
	for(int i=0; i<n; i++){
		fprintf(fp, "%0.12lf %0.12lf %0.12lf %0.12lf %0.12lf %0.12lf\n", 
				x[i], y[i], omega[i], alpha[i], u[i], v[i]);
	}
	fclose(fp);
}

void VortexSolver::calc_velocity(double *xtarget, double *ytarget, double *xsource, double *ysource, double *alpha_source, double *ui, double *vi){
	double delta = 1.0;
	double delta2 = pow(delta, 2);
	
#pragma omp parallel for
	for(int i=0; i<n; i++){
		double dxij, dyij, r2, rt2, u_i, v_i;
		u_i = 0.0;
		v_i = 0.0;
		for(int j=0; j<n; j++){
			dxij = xtarget[i] - xsource[j];
			dyij = ytarget[i] - ysource[j];
			r2 = pow(dxij, 2) + pow(dyij, 2);
			rt2 = r2/delta2;
			double tmp = alpha_source[j]/(2.0*M_PI*(r2 + 1e-16))*(1.0 - exp(-rt2));
			u_i += dyij*tmp;
			v_i -= dxij*tmp;
		}
		ui[i] = u_i;
		vi[i] = v_i;
	}
}

void VortexSolver::integrate(void){
	for(int i=0; i<n; i++){
		xrk[i] = x[i];
		yrk[i] = y[i];
	}
	calc_velocity(xrk, yrk, xrk, yrk, alpha, k1u, k1v);
	for(int i=0; i<n; i++){
		xrk[i] = x[i] + k1u[i]*dt*0.5;
		yrk[i] = y[i] + k1v[i]*dt*0.5;
	}

	calc_velocity(xrk, yrk, xrk, yrk, alpha, k2u, k2v);

	for(int i=0; i<n; i++){
		xrk[i] = x[i] + k2u[i]*dt*0.5;
		yrk[i] = y[i] + k2v[i]*dt*0.5;
	}

	calc_velocity(xrk, yrk, xrk, yrk, alpha, k3u, k3v);

	for(int i=0; i<n; i++){
		xrk[i] = x[i] + k3u[i]*dt;
		yrk[i] = y[i] + k3v[i]*dt;
	}

	calc_velocity(xrk, yrk, xrk, yrk, alpha, k4u, k4v);
	for(int i=0; i<n; i++){
		u[i] = (k1u[i] + 2.0*k2u[i] + 2.0*k3u[i] + k4u[i])/6.0;
		v[i] = (k1v[i] + 2.0*k2v[i] + 2.0*k3v[i] + k4v[i])/6.0;
		x[i] += u[i]*dt;
		y[i] += v[i]*dt;
	}
}

void VortexSolver::run(double tf){
	int step = 0;
	
	struct timespec start, finish;
	double elapsed;

	while(1){
		// timing start
		clock_gettime(CLOCK_MONOTONIC, &start);
	
		integrate();

		// timings calculations
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		
		if(step%write_freq == 0) write_solution(step);

		printf("TIME: %0.8e CPU TIME: %0.8fs\n", t, elapsed);
		t += dt;
		step += 1;
		
		if(t > tf) break;
	}
}

int main(int argc, char** argv){
	VortexSolver vs(50, 25, 0.1);
	vs.run(100.0);
	return 0;
}

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>

#define NSYSTEMS 5

#define T 3

#define UI 5
#define UJ 1

struct grain {
	//Euler angles in radians
	double phi1, phi, phi2;
	gsl_matrix *rot_mat;
	std::vector<double> crss;
	grain() {}
	grain(double p1, double p, double p2);
	~grain();
};

struct vect {
	double x, y, z;
};

struct slip_system {
	std::string label;
	//planes
	vect n;
	//directions
	vect b;
	double crss;
	//deformation matrix
	gsl_matrix* def_mat=0;

	//symetry part of deformation matrix
	gsl_matrix* symp_def=0;

	//non symetry part of deformation matrix
	gsl_matrix* asymp_def=0;

	std::vector<double> activity;

	slip_system(const std::string &l, const vect &in, const vect &ib, double r, int max_def);
	~slip_system();
};

struct slip_system_matrix {
    //m(ij)=bi*nj

    //slip systems matrix
    gsl_matrix *mat=0;

    //inverse slip systems matrix
    gsl_matrix *matinv=0;

    gsl_matrix *shear_stress=0;

    //combination of slip systems
    std::vector<int> comb;
    slip_system_matrix();
    ~slip_system_matrix();
};

struct handler {
    gsl_matrix *hand;
	handler();
	~handler();

};

struct U_handler {
    gsl_matrix *U_hand;
	U_handler();
	~U_handler();
};

#endif // STRUCTURES_H

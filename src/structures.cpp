#include <string>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "includes/structures.h"

void norm(vect &vr);

grain::grain(double p1, double p, double p2) : phi1(p1), phi(p), phi2(p2) {
    rot_mat = gsl_matrix_alloc(3,3);
}

grain::~grain() {
    if(rot_mat!=0)
        gsl_matrix_free(rot_mat);
}

slip_system::slip_system(const std::string &l, const vect &in, const vect &ib, double r, int max_def) : label(l), crss(r) {

	n = {in.x, in.y, in.z};
	b = {ib.x, ib.y, ib.z};

	norm(n);
	norm(b);

    for(int i=0;i<max_def;i++) {
        activity.push_back(0);
    }

	gsl_matrix *bb = gsl_matrix_alloc(T,1);
	gsl_matrix_set(bb, 0, 0, b.x);
    gsl_matrix_set(bb, 1, 0, b.y);
    gsl_matrix_set(bb, 2, 0, b.z);

    gsl_matrix *nn = gsl_matrix_alloc(1,T);
    gsl_matrix_set(nn, 0, 0, n.x);
    gsl_matrix_set(nn, 0, 1, n.y);
    gsl_matrix_set(nn, 0, 2, n.z);

	def_mat = gsl_matrix_alloc(T,T);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, bb, nn, 0.0, def_mat);

	gsl_matrix* helper = gsl_matrix_alloc(T,T);
    gsl_matrix_transpose_memcpy(helper, def_mat);

    symp_def = gsl_matrix_alloc(T,T);
    gsl_matrix_memcpy(symp_def, def_mat);
    gsl_matrix_add(symp_def, helper);
    gsl_matrix_scale(symp_def, 0.5);

    asymp_def = gsl_matrix_alloc(T,T);
    gsl_matrix_memcpy(asymp_def,def_mat);
    gsl_matrix_sub(asymp_def, helper);
    gsl_matrix_scale(asymp_def, 0.5);

	gsl_matrix_free(nn);
	gsl_matrix_free(bb);
	gsl_matrix_free(helper);

}

slip_system::~slip_system() {

    if(def_mat!=0)
        gsl_matrix_free(def_mat);
    if(symp_def!=0)
        gsl_matrix_free(symp_def);
	if(asymp_def!=0)
        gsl_matrix_free(asymp_def);

}

slip_system_matrix::slip_system_matrix() {
    mat = gsl_matrix_alloc(NSYSTEMS,NSYSTEMS);
    matinv = gsl_matrix_alloc(NSYSTEMS, NSYSTEMS);
    shear_stress = gsl_matrix_alloc(NSYSTEMS,UJ);
}

slip_system_matrix::~slip_system_matrix() {
    if(mat!=0)
        gsl_matrix_free(mat);
    if(matinv!=0)
        gsl_matrix_free(matinv);
    if(shear_stress!=0)
        gsl_matrix_free(shear_stress);
}

handler::handler() {
    hand = gsl_matrix_alloc(T,T);
}

handler::~handler() {
    gsl_matrix_free(hand);
}

U_handler::U_handler() {
    U_hand = gsl_matrix_alloc(UI,UJ);
}
U_handler::~U_handler() {
    gsl_matrix_free(U_hand);
}

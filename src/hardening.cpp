#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>

#include <gsl/gsl_matrix.h>

#include "includes/hardening.h"

Hardening::Hardening(const int siza, const int sizb) {

    if(siza <=0 || sizb <=0)
        throw std::runtime_error("Wrong parameters for hardening model");

    hardening_mat = gsl_matrix_alloc(siza, sizb);

    if(hardening_mat == NULL)
        throw std::runtime_error("Cant initialize hardening matrix in hardening model");

    hardening_factor=1;
    crss_0=1;
    power=1;

}

Hardening::Hardening(const Hardening& h) {

    if(h.get_hardening_matrix() == NULL) {
        throw std::runtime_error("Copy constructor of Hardening, hardening matrix in parent hardening is null");
    }

    hardening_mat = gsl_matrix_alloc(h.get_hardening_matrix()->size1, h.get_hardening_matrix()->size2);

    if(hardening_mat == NULL)
        throw std::runtime_error("Copy constructor of Hardening, cant initialize hardening matrix in hardening model");

    if(gsl_matrix_memcpy(hardening_mat,h.get_hardening_matrix()) != 0) {
        throw std::runtime_error("Copy constructor of Hardening, catn copy matrices");
    }

    hardening_factor=h.get_factor();
    crss_0=h.get_crss_0();
    power=h.get_power();

}

Hardening& Hardening::operator=(const Hardening& h) {

    if(&h!=this) {

        if(h.get_hardening_matrix() == NULL) {
            throw std::runtime_error("Assignment operator of Hardening, hardening matrix in parent hardening is null");
        }

        gsl_matrix_free(hardening_mat);

        hardening_mat = gsl_matrix_alloc(h.get_hardening_matrix()->size1, h.get_hardening_matrix()->size2);

        if(hardening_mat == NULL)
            throw std::runtime_error("Assignment operator of Hardening, cant initialize hardening matrix in hardening model");

        if(gsl_matrix_memcpy(hardening_mat,h.get_hardening_matrix()) != 0) {
            throw std::runtime_error("Assignment operator of Hardening, catn copy matrices");
        }

        hardening_factor=h.get_factor();
        crss_0=h.get_crss_0();
        power=h.get_power();
    }

    return *this;

}

void Hardening::load_data(const std::string &hardening_file_name) {

    std::cout<<"Hardening model exp: loading hardening"<<std::endl;

  	std::ifstream hardening_file (hardening_file_name);

	if (hardening_file.is_open()) {

        double d;

        for(int i=0;i<hardening_mat->size1;i++) {
            for(int j=0;j<hardening_mat->size2;j++) {
                hardening_file >> d;
                gsl_matrix_set(hardening_mat,i,j,d);
            }
        }

        hardening_file >> hardening_factor;
        hardening_file >> crss_0;
        hardening_file >> power;

		hardening_file.close();

	} else throw std::runtime_error("No hardening file");

	std::cout<<"...done hardening model exp: loading hardening data"<<std::endl;

}

void Hardening::update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef) const {


    double result;
    int ims;
    int is;

    for(int i=0;i<comb.size();i++) {

        result = 0;
        ims = comb[i];

        for(int j=0;j<comb.size();j++) {
            is = comb[j];
            result += gsl_matrix_get(hardening_mat,ims,is)*hardening_factor*pow(1-g->crss[ims]/crss_0,power)*fabs(gsl_matrix_get(udef.U_hand,j,0));
        }
        g->crss[ims] += result;
    }

}



gsl_matrix* Hardening::get_hardening_matrix() const {
    return hardening_mat;
}

double Hardening::get_factor() const {
    return hardening_factor;
}
double Hardening::get_crss_0() const {
    return crss_0;
}
double Hardening::get_power() const {
    return power;
}



Hardening::~Hardening() {

    if(hardening_mat != NULL)
        gsl_matrix_free(hardening_mat);

}

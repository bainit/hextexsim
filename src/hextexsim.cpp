#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "includes/util.h"
#include "includes/functions.h"
#include "includes/structures.h"
#include "includes/hardening.h"

using namespace std;

int main(int argc, char *argv[]) {

	cout<<"hestexsim: software for texture simulations in metals"<<endl;

    if(!validate_inputs(argv)) {
        usage();
        exit(1);
	}

	Matrices_keeper mkeeper;

	//Identity matrix
	gsl_matrix *I = mkeeper.alloc_next(3,3);

	//deformation tensor
	gsl_matrix *deftensor = mkeeper.alloc_next(3,3);

	//constrains
	gsl_matrix *k1 = mkeeper.alloc_next(3,3);
	gsl_matrix *k2 = mkeeper.alloc_next(3,3);

	//gsl_matrix *hk1 = gsl_matrix_alloc(3,3);
	//gsl_matrix *hk2 = gsl_matrix_alloc(3,3);

	gsl_matrix_set_identity(I);


	double max_def;
	double def_step;
	double evM;

	string settings_file = string(argv[1]);

    load_settings(&max_def, &def_step, deftensor, k1, k2, settings_file);

	gsl_matrix_scale(deftensor,def_step);
	gsl_matrix_scale(k1,def_step);
	gsl_matrix_scale(k2,def_step);

	evM = calc_evM(deftensor);


	string grains_file = string(argv[2]);
	vector<grain*> *grains;
	grains = prepare_rot_matrices(grains_file);

	int total_grains_num = grains->size();

    string slip_systems_file = string(argv[3]);
	vector<slip_system*> slip_systems;

	load_slip_systems(slip_systems, slip_systems_file, max_def);

	fill_grain_crss(grains,slip_systems);

    string hardening_file = string(argv[4]);
    unique_ptr<Hardening> hard = unique_ptr<Hardening>(new Hardening(slip_systems.size(), slip_systems.size()));
    hard->load_data(hardening_file);

	vector<slip_system_matrix*> matrices;
	prepare_slip_systems_matrices(matrices, slip_systems);

	vector<double> total_taylor_factor;

	for(int i=0;i<max_def;i++) {

        cout<<"def: "<<i<<endl;
        double taylor_factor=0;

        #pragma omp parallel for
        for(int j=0;j<total_grains_num;j++) {

            handler def;
            U_handler udef;
            double internal_energy;

            handler handk1;
            U_handler uhandk1;

            handler handk2;
            U_handler uhandk2;

            Matrices_keeper keeper;

            gsl_matrix* shear = keeper.alloc_next(5,1);
            gsl_matrix* shear_min = keeper.alloc_next(5,1);
            gsl_matrix* k1helper = keeper.alloc_next(5,1);
            gsl_matrix* k2helper = keeper.alloc_next(5,1);

            double min_internal_energy=100;
            double max_energy=-1;
            double energy=0;

            slip_system_matrix* mat=0;

            grain *g = grains->at(j);

            gsl_matrix *rm = g->rot_mat;

            calc_handler(def, deftensor, rm);
            get_uhandler(udef, def);

            int k;

            for(k=0;k<matrices.size();k++) {

                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrices[k]->matinv, udef.U_hand, 0.0, shear);

                    internal_energy = calc_internal_energy(shear,matrices[k], g);

                    if(internal_energy<=min_internal_energy) {

                        energy=calc_energy(matrices[k]->mat, def.hand, matrices[k], g);

                        if(internal_energy<min_internal_energy || (internal_energy==min_internal_energy && energy>max_energy)) {
                            min_internal_energy=internal_energy;
                            max_energy=energy;
                            mat = matrices[k];
                            gsl_matrix_memcpy(shear_min, shear);
                        }
                    }
            }

            taylor_factor += min_internal_energy/evM;

            calc_handler(handk1,k1,rm);
            get_uhandler(uhandk1,handk1,1);

            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat->matinv, uhandk1.U_hand, 0.0, k1helper);

            calc_handler(handk2,k2,rm);
            get_uhandler(uhandk2,handk2,1);

            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat->matinv, uhandk2.U_hand, 0.0, k2helper);

            update_rot_matrix(mat->comb, shear_min, slip_systems, k1helper, k2helper, rm, i);
            hard->update_hardening(mat->comb, grains->at(j), udef);

        }

        total_taylor_factor.push_back(taylor_factor/total_grains_num);

	}

	cleanup_energy();

	string export_path = string(argv[5])+"/texture_result.csv";

    for(int i=0;i<grains->size();i++) {
        grains->at(i) = grain_from_rot_matrix(grains->at(i)->rot_mat);
    }

    export_grains_to_file(*grains,export_path," ");

	string taylor_factor_file = string(argv[5])+"/taylor_factor.csv";
    export_taylor_factor_to_file(total_taylor_factor, taylor_factor_file);

	string slip_system_activity_file = string(argv[5])+"/systems_activity.csv";
	export_slip_system_activity_to_file(slip_systems, slip_system_activity_file, total_grains_num, evM);

	free_grains(*grains);

    free_slip_systems_matrices(matrices);
	free_slip_systems(slip_systems);
}

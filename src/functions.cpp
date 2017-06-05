#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <cmath>

#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_linalg.h>

#include "includes/functions.h"
#include "includes/structures.h"

using namespace std;

double rad2deg(double rad) {
    return rad*180/M_PI;
}

double deg2rad(double deg) {
    return deg*M_PI/180;
}

void print_matrix(const gsl_matrix *m) {

	size_t i,j;
	for(i=0;i<m->size1;i++) {
		for(j=0;j<m->size2;j++) {
			printf("%.4f ", gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}

}

double vect_len(const vect &v) {

	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

}

void norm(vect &v) {
	double len = vect_len(v);

	v.x = v.x/len;
	v.y = v.y/len;
	v.z = v.z/len;

}

double double_dot_prod(const gsl_matrix* m1, const gsl_matrix* m2) {

    double c1 = gsl_matrix_get(m1,0,0)*gsl_matrix_get(m2,0,0) + gsl_matrix_get(m1,1,1)*gsl_matrix_get(m2,1,1) + gsl_matrix_get(m1,2,2)*gsl_matrix_get(m2,2,2);
    double c2 = gsl_matrix_get(m1,0,1)*gsl_matrix_get(m2,0,1) + gsl_matrix_get(m1,1,0)*gsl_matrix_get(m2,1,0);
    double c3 = gsl_matrix_get(m1,0,2)*gsl_matrix_get(m2,0,2) + gsl_matrix_get(m1,2,0)*gsl_matrix_get(m2,2,0);
    double c4 = gsl_matrix_get(m1,1,2)*gsl_matrix_get(m2,1,2) + gsl_matrix_get(m1,2,1)*gsl_matrix_get(m2,2,1);

    return c1+c2+c3+c4;

}

double calc_evM(const gsl_matrix* def) {

    //Equivalent von Misses strain
    double exx = 2.0/3.0*gsl_matrix_get(def, 0,0)-1.0/3.0*gsl_matrix_get(def, 1,1)-1.0/3.0*gsl_matrix_get(def,2,2);
    double eyy = -1.0/3.0*gsl_matrix_get(def,0,0)+2.0/3.0*gsl_matrix_get(def,1,1)-1.0/3.0*gsl_matrix_get(def,2,2);
    double ezz = -1.0/3.0*gsl_matrix_get(def, 0,0)-1.0/3.0*gsl_matrix_get(def,1,1)+2.0/3.0*gsl_matrix_get(def,2,2);
    double exy = 2*gsl_matrix_get(def,0,1);
    double eyz = 2*gsl_matrix_get(def,1,2);
    double ezx = 2*gsl_matrix_get(def,2,0);

    //evM = 2/3*sqrt(3/2*(exx*exx + eyy*eyy + ezz*ezz) + 3/4*(exy^2 + eyz^2 + ezx^2));
    //evM = sqrt(2/3*(exx*exx + eyy*eyy + ezz*ezz+exy*exy+eyz*eyz+ezx*ezx) );
    //evM = sqrt((2/9)*((exx-eyy)^2+(eyy-ezz)^2+(ezz-exx)^2)+(1/3)*(exy^2+eyz^2+ezx^2));
    double evM = 2.0/3.0*sqrt(3.0/2.0*(exx*exx + eyy*eyy + ezz*ezz) + 3.0/4.0*(exy*exy + eyz*eyz + ezx*ezx));
    return evM;

}

double calc_internal_energy(const gsl_matrix* s, const slip_system_matrix *smatrix, grain *g) {

    double result = 0;
    int ims = 0;

    for(int i=0;i<5;i++) {
        ims = smatrix->comb[i];
        //result += abs(gsl_matrix_get(s,i,0)*gsl_matrix_get(st,i,0));
        result += fabs(gsl_matrix_get(s,i,0)*g->crss.at(ims));
    }

    return result;

    //return abs(gsl_matrix_get(s, 0,0)) + abs(gsl_matrix_get(s, 1,0)) + abs(gsl_matrix_get(s, 2,0)) + abs(gsl_matrix_get(s, 3,0)) + abs(gsl_matrix_get(s, 4,0));

}

double calc_energy(const gsl_matrix* slpsysmat, const gsl_matrix* def_l, const slip_system_matrix *smatrix, grain *g, bool flag) {

    gsl_matrix* stress = gsl_matrix_alloc(5,1);
    gsl_matrix* shear_stress = gsl_matrix_alloc(5,1);

    if(flag) {
        gsl_matrix_free(stress);
        gsl_matrix_free(shear_stress);
        return 0;
    }

    int ims = 0;

    for(int i=0;i<5;i++) {

        ims = smatrix->comb[i];
        gsl_matrix_set(shear_stress,i,0,g->crss.at(ims));

    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, slpsysmat, shear_stress, 0.0, stress);

    int result = fabs((-1)*gsl_matrix_get(stress,1,0)*gsl_matrix_get(def_l,0,0) + ((-1)*gsl_matrix_get(stress,0,0)+gsl_matrix_get(stress,1,0))*gsl_matrix_get(def_l,1,1)+2*gsl_matrix_get(stress,2,0)*gsl_matrix_get(def_l,1,2) + 2*gsl_matrix_get(stress,3,0)*gsl_matrix_get(def_l,0,2) + 2*gsl_matrix_get(stress,4,0)*gsl_matrix_get(def_l,0,1));

    gsl_matrix_free(stress);
    gsl_matrix_free(shear_stress);

    return result;
}

void cleanup_energy() {

    calc_energy(0,0,0,0,true);

}

double update_rot_matrix(const vector<int>& comb, const gsl_matrix* shearmin, const vector<slip_system*>& systems, const gsl_matrix* uk1, const gsl_matrix* uk2, gsl_matrix* rmat, int def_step) {

    gsl_matrix* r = gsl_matrix_alloc(T,T);
    gsl_matrix* result = gsl_matrix_alloc(T,T);

    double tshear;

    for(int i=0;i<comb.size();i++) {

        slip_system *s = systems[comb[i]];

        gsl_matrix_memcpy(r,s->asymp_def);

        //Constrains
        tshear = gsl_matrix_get(shearmin,i,0) + fabs(gsl_matrix_get(uk1,i,0)) + fabs(gsl_matrix_get(uk2,i,0));

        //No constrains
        //tshear = gsl_matrix_get(shearmin,i,0);

        gsl_matrix_scale(r,tshear);

        gsl_matrix_add(result,r);

        s->activity[def_step] += fabs(gsl_matrix_get(shearmin,i,0));

    }

    gsl_matrix_set(result,0,0,gsl_matrix_get(result,0,0)+1);
    gsl_matrix_set(result,1,1,gsl_matrix_get(result,1,1)+1);
    gsl_matrix_set(result,2,2,gsl_matrix_get(result,2,2)+1);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, result, rmat, 0.0, r);

    gsl_matrix_memcpy(rmat, r);

    gsl_matrix_free(r);
    gsl_matrix_free(result);

}

//void update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef, gsl_matrix *hardening) {
//
//    double result;
//    int ims;
//    int is;
//
//    for(int i=0;i<comb.size();i++) {
//
//        result = 0;
//        ims = comb[i];
//
//        for(int j=0;j<comb.size();j++) {
//            is = comb[j];
//            result += gsl_matrix_get(hardening,ims,is)*9*pow(1-g->crss[ims]/6,3.25)*abs(gsl_matrix_get(udef.U_hand,j,0));
//        }
//        g->crss[ims] = result;
//    }
//
//}

//in Bunge
gsl_matrix* get_rot_matrix(const grain *g, bool free_m) {

	gsl_matrix *z1;
	gsl_matrix *x;
	gsl_matrix *z2;

	gsl_matrix *temp;

//	static bool flag = 1;

//	if(free_m) {
//		gsl_matrix_free(z1);
//		gsl_matrix_free(x);
//		gsl_matrix_free(z2);
//		gsl_matrix_free(temp);
//
//		flag = 1;
//
//		return 0;
//	}

	gsl_matrix *r = gsl_matrix_alloc(3,3);

	//if(flag) {

		z1 = gsl_matrix_alloc(3,3);
		x = gsl_matrix_alloc(3,3);
		z2 = gsl_matrix_alloc(3,3);

		gsl_matrix_set(z1,0,2,0);
		gsl_matrix_set(z1,1,2,0);
		gsl_matrix_set(z1,2,0,0);
		gsl_matrix_set(z1,2,1,0);
		gsl_matrix_set(z1,2,2,1);

		gsl_matrix_set(x,0,0,1);
		gsl_matrix_set(x,0,1,0);
		gsl_matrix_set(x,0,2,0);
		gsl_matrix_set(x,1,0,0);
		gsl_matrix_set(x,2,0,0);
		gsl_matrix_set(x,2,0,0);

		gsl_matrix_set(z2,0,2,0);
		gsl_matrix_set(z2,1,2,0);
		gsl_matrix_set(z2,2,0,0);
		gsl_matrix_set(z2,2,1,0);
		gsl_matrix_set(z2,2,2,1);

		temp = gsl_matrix_alloc(3,3);

//		flag = 0;
	//}

	gsl_matrix_set(z1,0,0,cos(g->phi1));
	gsl_matrix_set(z1,0,1,sin(g->phi1));
	gsl_matrix_set(z1,1,0,-sin(g->phi1));
	gsl_matrix_set(z1,1,1,cos(g->phi1));


	gsl_matrix_set(x,1,1,cos(g->phi));
	gsl_matrix_set(x,1,2,sin(g->phi));
	gsl_matrix_set(x,2,1,-sin(g->phi));
	gsl_matrix_set(x,2,2,cos(g->phi));


	gsl_matrix_set(z2,0,0,cos(g->phi2));
	gsl_matrix_set(z2,0,1,sin(g->phi2));
	gsl_matrix_set(z2,1,0,-sin(g->phi2));
	gsl_matrix_set(z2,1,1,cos(g->phi2));


	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, z2, x, 0.0, temp);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, temp, z1, 0.0, r);

	gsl_matrix_free(z1);
	gsl_matrix_free(x);
	gsl_matrix_free(z2);
	gsl_matrix_free(temp);

	return r;

}

grain* grain_from_rot_matrix(gsl_matrix *rot_mat) {

	double phi1,phi,phi2;

    	if(fabs(gsl_matrix_get(rot_mat,2,2)) < 1) {

       		phi = acos(gsl_matrix_get(rot_mat,2,2));
        	phi1 = atan2(gsl_matrix_get(rot_mat,2,0)/sin(M_PI),-gsl_matrix_get(rot_mat,2,1)/sin(M_PI));
        	phi2 = atan2(gsl_matrix_get(rot_mat,0,2)/sin(M_PI),gsl_matrix_get(rot_mat,1,2)/sin(M_PI));

	}

    	else {
        	phi = 0;
        	phi1 = atan2(gsl_matrix_get(rot_mat,0,1),gsl_matrix_get(rot_mat,0,0))/2;
        	phi2 = phi1;
	}

    	grain *g = new grain(phi1,phi,phi2);

	return g;
}


void load_settings(double *max_def, double *def_step, gsl_matrix *deformation, gsl_matrix *k1, gsl_matrix *k2, const string &settings_file_name) {

	cout<<"Loading settings"<<endl;

  	ifstream settings_file(settings_file_name);

	if (settings_file.is_open()) {
		settings_file >> *max_def;
		settings_file >> *def_step;

		size_t i,j;
		double v;

		for(i=0;i<deformation->size1;i++) {
			for(j=0;j<deformation->size2;j++) {
				settings_file >> v;
				gsl_matrix_set(deformation,i,j,v);
			}
		}

		for(i=0;i<k1->size1;i++) {
			for(j=0;j<k1->size2;j++) {
				settings_file >> v;
				gsl_matrix_set(k1,i,j,v);
			}
		}

		for(i=0;i<k2->size1;i++) {
			for(j=0;j<k2->size2;j++) {
				settings_file >> v;
				gsl_matrix_set(k2,i,j,v);
			}
		}

        settings_file.close();

	}else throw runtime_error("No settings file");

	cout<<"...done loading settigs"<<endl;

}

void load_grains(vector<grain*> &grains, const string &grains_file_name) {

	cout<<"Loading grains"<<endl;

  	ifstream grains_file (grains_file_name);

  	double d;

	if (grains_file.is_open()) {
		while(grains_file >> d) {

			grain *g = new grain;

			g->phi1 = d;
			grains_file >> g->phi;
			grains_file >> g->phi2;
			grains.push_back(g);
		}
		grains_file.close();
	}else throw runtime_error("No grains file");

	cout<<"...done loading grains"<<endl;

}

void load_slip_systems(vector<slip_system*> &slip_systems, const string &slip_systems_file_name, int max_def) {

	cout<<"Loading slip systems"<<endl;

  	ifstream slip_systems_file (slip_systems_file_name);

	if (slip_systems_file.is_open()) {
        string l;
        while(slip_systems_file >> l) {

            double n1,n2,n3;
            double b1,b2,b3;
            double crss;

                //slip_systems_file >> l >> n1 >> n2 >> n3 >> b1 >> b2 >> b3 >> crss;
            slip_systems_file >> n1 >> n2 >> n3 >> b1 >> b2 >> b3 >> crss;

            slip_system *s = new slip_system(l, {n1,n2,n3}, {b1,b2,b3}, crss, max_def);

            slip_systems.push_back(s);
        }
		slip_systems_file.close();
	} else throw runtime_error("No slip systems file");

	cout<<"...done loading slip systems"<<endl;

}

vector<grain*>* prepare_rot_matrices(const string &grains_file_name) {

    cout<<"Preparing rotation matrices"<<endl;
	//grains
	vector<grain*> *grains = new vector<grain*>();

	try{
        load_grains(*grains,grains_file_name);
    }catch(runtime_error& e) {
        delete grains;
        throw runtime_error(e.what());
    }

	for(int i=0;i<grains->size();i++) {

        grain* g = grains->at(i);

		g->rot_mat = get_rot_matrix(g);
	}
    //get_rot_matrix(0,true);

	cout<<"...done preparing rotation matrices"<<endl;

	return grains;
}

void fill_grain_crss(std::vector<grain*> *grains, const std::vector<slip_system*> &slip_systems) {

    for(int i=0;i<grains->size();i++) {
        grain *g = grains->at(i);
        for(int j=0;j<slip_systems.size();j++) {
            g->crss.push_back(slip_systems[j]->crss);
        }
    }
}

void prepare_slip_systems_matrices(std::vector<slip_system_matrix*> &matrices, std::vector<slip_system*> &slip_systems) {

    cout<<"Preparing slip system matrices"<<endl;
    gsl_combination *c;
    c = gsl_combination_calloc(slip_systems.size(), 5);

    gsl_matrix *mat_temp = gsl_matrix_alloc(NSYSTEMS,NSYSTEMS);
    gsl_matrix *sstress = gsl_matrix_alloc(NSYSTEMS,UJ);

    do{
       vector<int> comb(NSYSTEMS);

       gsl_matrix *mat = gsl_matrix_alloc(NSYSTEMS,NSYSTEMS);

        for(int i=0;i<NSYSTEMS;i++) {

            comb[i] = gsl_combination_get(c,i);
            slip_system *s = slip_systems[comb[i]];

            gsl_matrix_set(mat_temp, 0, i, s->b.y * s->n.y);
            gsl_matrix_set(mat_temp, 1, i, s->b.z * s->n.z);

            gsl_matrix_set(mat_temp, 2, i, (s->b.y * s->n.z + s->b.z * s->n.y));
            gsl_matrix_set(mat_temp, 3, i, (s->b.x * s->n.z + s->b.z * s->n.x));
            gsl_matrix_set(mat_temp, 4, i, (s->b.x * s->n.y + s->b.y * s->n.x));

            gsl_matrix_set(sstress, i, 0, s->crss);

        }

        gsl_matrix_memcpy(mat,mat_temp);

        int s;
        gsl_permutation *p = gsl_permutation_alloc(NSYSTEMS);
        gsl_linalg_LU_decomp (mat_temp, p, &s);

        double det = gsl_linalg_LU_det(mat_temp, s);

        if(det==0) {
            gsl_matrix_free(mat);
            gsl_permutation_free(p);
            continue;
        }

        gsl_matrix *matinv = gsl_matrix_alloc(NSYSTEMS,NSYSTEMS);

        gsl_linalg_LU_invert(mat_temp, p, matinv);

        slip_system_matrix *ssm = new slip_system_matrix;
        gsl_matrix_memcpy(ssm->mat, mat);
        gsl_matrix_memcpy(ssm->matinv, matinv);
        gsl_matrix_memcpy(ssm->shear_stress,sstress);
        ssm->comb = comb;

        matrices.push_back(ssm);

        gsl_matrix_free(mat);
        gsl_permutation_free(p);
        gsl_matrix_free(matinv);

    } while (gsl_combination_next(c) == GSL_SUCCESS);

    gsl_matrix_free(mat_temp);
    gsl_matrix_free(sstress);

    cout<<"...done preparing slip systems matrices"<<endl;

}

void convert2U(const gsl_matrix* mat, gsl_matrix* umat) {

    gsl_matrix_set(umat,0,0, gsl_matrix_get(mat,1,1));
    gsl_matrix_set(umat,1,0, gsl_matrix_get(mat,2,2));
    gsl_matrix_set(umat,2,0, gsl_matrix_get(mat,1,2));
    gsl_matrix_set(umat,3,0, gsl_matrix_get(mat,0,2));
    gsl_matrix_set(umat,4,0, gsl_matrix_get(mat,0,1));

}

void calc_handler(handler &h, const gsl_matrix *def, const gsl_matrix *rm) {

    gsl_matrix *helper = gsl_matrix_alloc(3,3);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rm, def, 0.0, helper);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, helper, rm, 0.0, h.hand);

    gsl_matrix_free(helper);

//    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rm, k1, 0.0, h.helper);
//    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, h.helper, rm, 0.0, h.constrains_k1);
//
//    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rm, k2, 0.0, h.helper);
//    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, h.helper, rm, 0.0, h.constrains_k2);

}

void get_uhandler(U_handler &uh, const handler &h, int s) {

    convert2U(h.hand,uh.U_hand);
//    convert2U(h.constrains_k1, uh.U_constrains_k1);
//    convert2U(h.constrains_k2, uh.U_constrains_k2);

    gsl_matrix_set(uh.U_hand, 2,0, s*gsl_matrix_get(uh.U_hand,2,0));
    gsl_matrix_set(uh.U_hand, 3,0, s*gsl_matrix_get(uh.U_hand,3,0));
    gsl_matrix_set(uh.U_hand, 4,0, s*gsl_matrix_get(uh.U_hand,4,0));

}

void write_grains(const std::vector<grain*> &grains, ofstream &file, const string &delimiter) {

    for(int i=0;i<grains.size();i++) {
			grain *g = grains[i];
			file << g->phi1 << delimiter << g->phi << delimiter << g->phi2 << delimiter << 1 <<endl;
		}
}


void export_grains_to_file(const std::vector<grain*> &grains, const string &export_file_name, const string &delimiter) {

    cout<<"Exporting grains to file"<<endl;
	ofstream export_file (export_file_name);

	cout<<export_file_name<<endl;

	if (export_file.is_open()) {

		write_grains(grains,export_file, delimiter);
		export_file.close();

	}else {

        cerr<<"Cant export grains to file: "<<export_file_name<<endl;
        cerr<<"Trying export grains to file: /tmp/texture_result.csv"<<endl;

        ofstream emergency_file("/tmp/texture_result.csv");

        if(emergency_file.is_open()) {

            write_grains(grains,emergency_file, delimiter);
            emergency_file.close();

        }else throw runtime_error("Cant export grains to file");

	}

    cout<<"...done exporting grains to file"<<endl;

}

void export_taylor_factor_to_file(const vector<double> &total_tf, const string &filename) {

    cout<<"Exporting Taylor factor to file"<<endl;
	ofstream export_file (filename);

	if (export_file.is_open()) {

		for(int i=0;i<total_tf.size();i++) {
			export_file << total_tf[i] << endl;
		}

		export_file.close();

	}else {

        cerr<<"Cant export Taylor factor to file: "<<filename<<endl;
        cerr<<"Trying export Taylor factor to file: /tmp/taylor_factor.csv"<<endl;

        ofstream emergency_file("/tmp/taylor_factor.csv");

        if(emergency_file.is_open()) {

            for(int i=0;i<total_tf.size();i++) {
                emergency_file << total_tf[i] << endl;
            }

            emergency_file.close();

        }else throw runtime_error("Cant export Taylor factor to file");

	}

    cout<<"...done exporting Taylor factor to file"<<endl;

}

void write_slip_system_activity(const std::vector<slip_system*> &slip_systems, std::ofstream &export_file, int grains_num, double evM){

    slip_system *st = slip_systems[0];
    int actvity_size = st->activity.size();

    for(int j=-1;j<actvity_size;j++) {
        for(int i=0;i<slip_systems.size();i++) {

            slip_system *s = slip_systems[i];

                if(j<0) {
                    if(i==0)
                        export_file << "- ";
                    else
                        export_file << s->label << " ";
                }
                else {
                    if(i==0)
                        export_file <<j*evM+evM<< " ";
                    else
                        export_file << s->activity[j]/grains_num << " ";
                }
            }
            export_file << endl;
    }

}

void export_slip_system_activity_to_file(const std::vector<slip_system*> &slip_systems, const std::string &filename, int grains_num, double evM) {


     cout<<"Exporting slip system activity to file"<<endl;
	 ofstream export_file (filename);

	if (export_file.is_open()) {

        write_slip_system_activity(slip_systems, export_file, grains_num, evM);
		export_file.close();

	}else {

        cerr<<"Cant export slip system activity to file: "<<filename<<endl;
        cerr<<"Trying export Taylor factor to file: /tmp/taylor_factor.csv"<<endl;

        ofstream emergency_file("/tmp/taylor_factor.csv");

        if(emergency_file.is_open()) {

            write_slip_system_activity(slip_systems, emergency_file, grains_num, evM);

            emergency_file.close();

        }else throw runtime_error("Cant export export slip system activity to file");
	}

    cout<<"...done exporting slip system activity to file"<<endl;

}

void free_grains(vector<grain*> &grains) {

    cout<<"Freeing grains"<<endl;
	for(int i;i<grains.size();i++) {
		delete grains[i];
	}
	cout<<"...done freeing grains"<<endl;

}

void free_slip_systems(vector<slip_system*> &slip_systems) {

    cout<<"Freeing slip systems"<<endl;
	for(int i;i<slip_systems.size();i++) {
		delete slip_systems[i];
	}
    cout<<"...done freeing slip systems"<<endl;
}

void free_rotation_matrices(std::vector<gsl_matrix*> rot_matrices) {

    cout<<"Freeing rotation matrices"<<endl;
    for(int i=0;i<rot_matrices.size();i++) {
        gsl_matrix_free(rot_matrices[i]);
    }
    cout<<"..done freeing rotation matrices"<<endl;

}

void free_slip_systems_matrices(vector<slip_system_matrix*> &matrices) {

    cout<<"Freeing slip systems matrices"<<endl;
	for(int i;i<matrices.size();i++) {
		delete matrices[i];
	}
    cout<<"...done freeing slip systems matrices"<<endl;

}

























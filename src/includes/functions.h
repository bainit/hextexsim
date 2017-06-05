#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <string>
#include <fstream>

#include <gsl/gsl_matrix.h>

#include "structures.h"

double rad2deg(double rad);
double deg2rad(double deg);
double vect_len(const vect &v);
void norm(vect &vr);


double double_dot_prod(const gsl_matrix* m1, const gsl_matrix* m2);
double calc_evM(const gsl_matrix* def);
double calc_internal_energy(const gsl_matrix* shear, const slip_system_matrix *smatrix, grain *g);
double calc_energy(const gsl_matrix* slpsysmat, const gsl_matrix* def_l, const slip_system_matrix *smatrix, grain *g, bool falg=false);
void cleanup_energy();

double update_rot_matrix(const std::vector<int>& comb, const gsl_matrix* shaermin, const std::vector<slip_system*>& systems, const gsl_matrix* uk1, const gsl_matrix* uk2, gsl_matrix* rmat, int def_step);


gsl_matrix* get_rot_matrix(const grain *g, bool free_m=false);
grain* grain_from_rot_matrix(gsl_matrix *rot_mat);



void print_matrix(const gsl_matrix *m);

void load_settings(double *max_def, double *def_step, gsl_matrix *deformation, gsl_matrix *k1, gsl_matrix *k2, const std::string &settings_file_name);
void load_grains(std::vector<grain*> &grains, const std::string &grains_file_name);
void load_slip_systems(std::vector<slip_system*> &slip_systems, const std::string &slip_system_file, int max_def);



std::vector<grain*>* prepare_rot_matrices(const std::string &grains_file_name);
void fill_grain_crss(std::vector<grain*> *grains, const std::vector<slip_system*> &slip_systems);
void prepare_slip_systems_matrices(std::vector<slip_system_matrix*> &matrices, std::vector<slip_system*> &slip_systems);


void convert2U(const gsl_matrix* mat, gsl_matrix* umat);
void calc_handler(handler &h, const gsl_matrix *def, const gsl_matrix *rm);
void get_uhandler(U_handler &uh, const handler &h, int s=2);


void write_grains(const std::vector<grain*> &grinas, std::ofstream &file, const std::string &delimiter);
void export_grains_to_file(const std::vector<grain*> &grains, const std::string &export_file, const std::string &delimiter);
void export_taylor_factor_to_file(const std::vector<double> &total_tf, const std::string& filename);
void write_slip_system_activity(const std::vector<slip_system*> &slip_systems, std::ofstream &export_file, int grains_num, double eVM);
void export_slip_system_activity_to_file(const std::vector<slip_system*> &slip_systems, const std::string& filename, int grains_num, double eVM);



void free_grains(std::vector<grain*> &grains);
void free_slip_systems(std::vector<slip_system*> &slip_systems);
void free_rotation_matrices(std::vector<gsl_matrix*> rot_matrices);
void free_slip_systems_matrices(std::vector<slip_system_matrix*> &matrices);


#endif	//FUNCTIONS_H

#ifndef HARDENIG_MODEL
#define HARDENING_MODEL

#include <string>
#include <vector>
#include "structures.h"

class hardening_model {

    public:
        virtual void load_data(const std::string &hardening_file)=0;
        virtual void update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef)=0;
        virtual ~hardening_model();
};

class hardening {

    private:
        hardening_model *model;
    public:
        hardening(hardening_model *m);
        hardening_model* get_model();
        ~hardening();
};

class hardening_exp_model : public hardening_model {

    private:
        gsl_matrix *hardening_mat;
        double hardening_factor;
        double crss_0;
        double power;

    public:
        hardening_exp_model(int siza, int sizeb);
        void load_data(const std::string &hardening_file);
        void update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef);
        ~hardening_exp_model();
        gsl_matrix* get_hardening_matrix();
        double get_factor();
        double get_crss_0();
        double get_power();

};

#endif // HARDENIG_MODEL

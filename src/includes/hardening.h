#ifndef HARDENIG
#define HARDENING

#include <string>
#include <vector>
#include "structures.h"

class Hardening {

    private:
        gsl_matrix *hardening_mat;
        double hardening_factor;
        double crss_0;
        double power;

    public:
        Hardening(const int siza, const int sizeb);

        Hardening(const Hardening& h);
        Hardening& operator=(const Hardening& h);

        void load_data(const std::string &hardening_file);

        void update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef) const;

        gsl_matrix* get_hardening_matrix() const;
        double get_factor() const;
        double get_crss_0() const;
        double get_power() const;

        ~Hardening();


};

#endif // HARDENIG


























//class Hardening_model {
//
//    public:
//        virtual void load_data(const std::string &hardening_file)=0;
//        virtual void update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef) const =0;
//        virtual ~Hardening_model();
//};
//
//class Hardening {
//
//    private:
//        Hardening_model *model=NULL;
//    public:
//        Hardening(Hardening_model *m);
//        Hardening(const Hardening& h);
//
//        Hardening& operator=(const Hardening& h);
//
//        Hardening_model* get_model() const;
//        ~Hardening();
//};
//
//class Hardening_exp_model : public Hardening_model {
//
//    private:
//        gsl_matrix *hardening_mat;
//        double hardening_factor;
//        double crss_0;
//        double power;
//
//    public:
//        Hardening_exp_model(const int siza, const int sizeb);
//
//        void load_data(const std::string &hardening_file);
//
//        void update_hardening(const std::vector<int> &comb, grain *g, const U_handler &udef) const;
//
//        gsl_matrix* get_hardening_matrix() const;
//        double get_factor() const;
//        double get_crss_0() const;
//        double get_power() const;
//
//        ~Hardening_exp_model();
//
//
//};



#ifndef HARDENING
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

#endif // HARDENING

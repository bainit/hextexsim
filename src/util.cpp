#include <iostream>
#include <stdexcept>

#include <sys/stat.h>

#include <gsl/gsl_matrix.h>

#include "includes/util.h"

using namespace std;

Matrices_keeper::Matrices_keeper(){

}

Matrices_keeper::Matrices_keeper(const Matrices_keeper& k) {

    gsl_matrix* m;
    gsl_matrix* d;
    int status;

    for(int i=0;i<k.matrices.size();i++) {

            m = k.matrices.at(i);
            if(m==NULL)
                throw runtime_error("Error in copy constructor of Matrices_keeper, gsl_matrix is null");

            d = this->alloc_next(m->size1, m->size2);
            if(d==NULL)
                throw runtime_error("Error in copy constructor of Matrices_keeper, cant alloc matrix");

            status = gsl_matrix_memcpy(d,m);
            if(status != 0)
                throw runtime_error("Error in copy constructor of Matrices_keeper, cant copy matrices");

    }

}

Matrices_keeper& Matrices_keeper::operator=(const Matrices_keeper& k) {

    if(this!=&k) {

        this->dalloc_all();

        gsl_matrix* m;
        gsl_matrix* d;
        int status;

        for(int i=0;i<k.matrices.size();i++) {

                m = k.matrices.at(i);
                if(m==NULL)
                    throw runtime_error("Error in assignment operator of Matrices_keeper, gsl_matrix is null");

                d = this->alloc_next(m->size1, m->size2);
                if(d==NULL)
                    throw runtime_error("Error in assignment operator of Matrices_keeper, cant alloc matrix");

                status = gsl_matrix_memcpy(d,m);
                if(status != 0)
                    throw runtime_error("Error in assignment operator of Matrices_keeper, cant copy matrices");

        }
    }

    return *this;

}

gsl_matrix* Matrices_keeper::alloc_next(int s1, int s2) {

    gsl_matrix *m = gsl_matrix_alloc(s1,s2);
    matrices.push_back(m);
    return m;
}

void Matrices_keeper::dalloc_all() {

    for(int i=0;i<matrices.size();i++) {
        gsl_matrix_free(matrices.at(i));
    }
    matrices.clear();
}

Matrices_keeper::~Matrices_keeper() {
    dalloc_all();
}

bool validate_inputs(char **args) {

    struct stat sb;

    if(args[1] == NULL) {
        cout<<"No setting file specified"<<endl;
        return false;
    }

    if (!(stat(args[1], &sb) == 0 && S_ISREG(sb.st_mode))){
        cout<<"Setting file does not exist"<<endl;
        return false;
    }

    if(args[2] == NULL) {
        cout<<"No grains file specified"<<endl;
        return false;
    }

    if (!(stat(args[2], &sb) == 0 && S_ISREG(sb.st_mode))){
        cout<<"Grains file does not exist"<<endl;
        return false;
    }

    if(args[3] == NULL) {
        cout<<"No slip system file specified"<<endl;
        return false;
    }

    if (!(stat(args[3], &sb) == 0 && S_ISREG(sb.st_mode))){
        cout<<"Slip system file does not exist"<<endl;
        return false;
    }

    if(args[4] == NULL) {
        cout<<"No hardening file specified"<<endl;
        return false;
    }

    if (!(stat(args[4], &sb) == 0 && S_ISREG(sb.st_mode))){
        cout<<"Hardening file does not exist"<<endl;
        return false;
    }

    if(args[5] == NULL) {
        cout<<"No output dir specified"<<endl;
        return false;
    }


    if (!(stat(args[5], &sb) == 0 && S_ISDIR(sb.st_mode))){
        cout<<"Output dir does not exist"<<endl;
        return false;
    }

    return true;

}

void usage() {

    cout<<"usage of hextexsim:"<<endl;
    cout<<"./hextexsim <settings_file> <grains_file> <slip_systems_file> <hardening_file> <output_dir>"<<endl;

}

void exit_on_status(int status, char *msg) {

    if(status != 0) {
        cout<<msg<<", status: "<<status<<endl;
        exit(1);
    }

}
void exit_on_null(void *p, char *msg) {

    if(p == NULL) {
        cout<<msg<<endl;
        exit(1);
    }

}

#ifndef UTIL_H
#define UTIL_H

#include <vector>

#include <gsl/gsl_matrix.h>

class Matrices_keeper {

    private:
        std::vector<gsl_matrix*> matrices;
    public:
        Matrices_keeper();
        Matrices_keeper(const Matrices_keeper& keeper);
        Matrices_keeper& operator=(const Matrices_keeper& keeper);

        gsl_matrix* alloc_next(int s1, int s2);
        void dalloc_all();
        ~Matrices_keeper();

};

bool validate_inputs(char **args);
void usage();

void exit_on_status(int status, char *msg);
void exit_on_null(void *p, char *msg);

#endif // UTIL_H

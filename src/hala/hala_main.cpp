#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>

#define GENCMA

#ifdef __cplusplus
extern "C" {
#endif


#include "cutest.h"
#include "cutest_routines.h"

#define HALA    hala
#define HALASPC halaspc
#define GETINFO getinfo

rp_ HALA(int, int, double*, int, double*, double*, double, double*, double*, double, double, double, int*);
void HALASPC(integer, char*);
void GETINFO(integer, integer, rp_*, rp_*, rp_*, rp_*, logical*, logical*, VarTypes*);

integer CUTEst_nvar;
integer CUTEst_ncon;
integer CUTEst_nnzj;
integer CUTEst_nnzh;

int MAINENTRY(void) {
    const char *fname = "OUTSDIF.d";
    integer funit = 42;
    integer iout = 6;
    integer io_buffer = 11;
    integer ierr;
    integer status;

    VarTypes vtypes;

    rp_ *x, *bl, *bu, *dummy1, *dummy2;
    rp_ *v = NULL, *cl = NULL, *cu = NULL;
    logical *equatn = NULL, *linear = NULL;
    char *pname, *vnames, *gnames, *cptr;
    char **Vnames, **Gnames;
    logical grad;
    integer e_order = 1, l_order = 0, v_order = 0;
    logical constrained = FALSE_;

    rp_ calls[7], cpu[4];
    integer nlin = 0, nbnds = 0, neq = 0;
    rp_ dummy;
    integer ExitCode;
    int i, j;

    ierr = 0;
    FORTRAN_open(&funit, fname, &ierr);
    if (ierr) {
        printf("Error opening file OUTSDIF.d.\nAborting.\n");
        exit(1);
    }

    CUTEST_cdimen_r(&status, &funit, &CUTEst_nvar, &CUTEst_ncon);
    if (status) {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    if (CUTEst_ncon) constrained = TRUE_;

    MALLOC(x, CUTEst_nvar, rp_);
    MALLOC(bl, CUTEst_nvar, rp_);
    MALLOC(bu, CUTEst_nvar, rp_);
    if (constrained) {
        MALLOC(equatn, CUTEst_ncon, logical);
        MALLOC(linear, CUTEst_ncon, logical);
        MALLOC(v, CUTEst_ncon, rp_);
        MALLOC(cl, CUTEst_ncon, rp_);
        MALLOC(cu, CUTEst_ncon, rp_);
        CUTEST_csetup_r(&status, &funit, &iout, &io_buffer,
                        &CUTEst_nvar, &CUTEst_ncon, x, bl, bu,
                        v, cl, cu, equatn, linear,
                        &e_order, &l_order, &v_order);
    } else {
        CUTEST_usetup_r(&status, &funit, &iout, &io_buffer,
                        &CUTEst_nvar, x, bl, bu);
    }
    if (status) {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    MALLOC(pname, FSTRING_LEN + 1, char);
    MALLOC(vnames, CUTEst_nvar * FSTRING_LEN, char);
    MALLOC(Vnames, CUTEst_nvar, char*);
    for (i = 0; i < CUTEst_nvar; i++)
        MALLOC(Vnames[i], FSTRING_LEN + 1, char);

    if (constrained) {
        MALLOC(gnames, CUTEst_ncon * FSTRING_LEN, char);
        MALLOC(Gnames, CUTEst_ncon, char*);
        for (i = 0; i < CUTEst_ncon; i++)
            MALLOC(Gnames[i], FSTRING_LEN + 1, char);
        CUTEST_cnames_r(&status, &CUTEst_nvar, &CUTEst_ncon,
                        pname, vnames, gnames);
    } else {
        CUTEST_unames_r(&status, &CUTEst_nvar, pname, vnames);
    }

    if (status) {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    pname[FSTRING_LEN] = '\0';
    printf(" Problem                 : %-s\n", pname);

    for (i = 0; i < CUTEst_nvar; i++) {
        cptr = vnames + i * FSTRING_LEN;
        for (j = 0; j < FSTRING_LEN; j++) {
            Vnames[i][j] = *cptr;
            cptr++;
        }
        Vnames[i][FSTRING_LEN] = '\0';
    }

    for (i = 0; i < CUTEst_ncon; i++) {
        cptr = gnames + i * FSTRING_LEN;
        for (j = 0; j < FSTRING_LEN; j++) {
            Gnames[i][j] = *cptr;
            cptr++;
        }
        Gnames[i][FSTRING_LEN] = '\0';
    }

    FREE(vnames);
    if (constrained) FREE(gnames);

    printf("Variable names:\n");
    for (i = 0; i < CUTEst_nvar; i++)
        printf("  %s\n", Vnames[i]);

    for (i = 0; i < CUTEst_nvar; i++) FREE(Vnames[i]);
    FREE(Vnames);

    if (constrained) printf("Constraint names:\n");
    for (i = 0; i < CUTEst_ncon; i++)
        printf("  %s\n", Gnames[i]);

    for (i = 0; i < CUTEst_ncon; i++) FREE(Gnames[i]);
    if (constrained) FREE(Gnames);

    GETINFO(CUTEst_nvar, CUTEst_ncon, bl, bu, cl, cu,
            equatn, linear, &vtypes);

    int maxIterations = 30;
    double tau_values[3] = {0.1, 10.0, 50.0};
    double tolerance_values[5] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    double lambda0_values[3] = {1.0, 10.0, 50.0};
    double alpha = 0.6;
    

    std::ofstream MyFile;
    MyFile.open("results.txt", std::ios_base::app);
    MyFile << "#" << pname << "\n";

    for (int t = 0; t < 3; t++) {
        for (int l = 0; l < 3; l++) {
            for (int tol = 0; tol < 5; tol++) {

                int currentIteration = 0;
                double tau = tau_values[t];
                double lambda0 = lambda0_values[l];
                double tolerance = tolerance_values[tol];

                clock_t start_time = clock();
                dummy = HALA(CUTEst_nvar, CUTEst_ncon, x, maxIterations, bl, bu, tau, cl, cu, tolerance, alpha, lambda0, &currentIteration);
                clock_t end_time = clock();
                double solve_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
                ExitCode = 0;

                CUTEST_creport_r(&status, calls, cpu);
                if (status) {
                    printf("** CUTEst error, status = %d, aborting\n", status);
                    exit(status);
                }
                /*
                printf("\n\n ************************ CUTEst statistics ************************\n\n");
                printf(" Code used               : HALA\n");
                printf(" Problem                 : %-s\n", pname);
                printf(" tau                     = %-10.2f\n", tau);
                printf(" lambda0                 = %-10.2f\n", lambda0);
                printf(" tolerance               = %-10.2e\n", tolerance);
                printf(" # variables             = %-10d\n", (int)CUTEst_nvar);
                printf(" # constraints           = %-10d\n", (int)CUTEst_ncon);
                printf(" # linear constraints    = %-10d\n", vtypes.nlin);
                printf(" # equality constraints  = %-10d\n", vtypes.neq);
                printf(" # inequality constraints= %-10d\n", vtypes.nineq);
                printf(" # bound constraints     = %-10d\n", vtypes.nbnds);
                printf(" # objective functions   = %-15.7g\n", calls[0]);
                printf(" # objective gradients   = %-15.7g\n", calls[1]);
                printf(" # objective Hessians    = %-15.7g\n", calls[2]);
                printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
                if (constrained) printf(" # Jacobians             = %-15.7g\n", calls[4]);
                if (constrained) printf(" # constraint gradients  = %-15.7g\n", calls[5]);
                if (constrained) printf(" # constraint Hessians   = %-15.7g\n", calls[6]);
                printf(" Final f                 = %-15.5e\n", dummy);
                printf(" # Iterations   = %-15d\n", currentIteration);
                printf(" Preparation time        = %-10.5f seconds\n", cpu[0]);
                printf(" Solve time              = %-10.5f seconds\n", solve_time);
                printf(" ******************************************************************\n\n");*/

                MyFile << "Tau: " << tau << "\n";
                MyFile << "Tolerance: " << tolerance << "\n";
                MyFile << "Lambda0: " << lambda0 << "\n";
                MyFile << "Final f: " << dummy << "\n";
                MyFile << "Iterations: " << currentIteration << "\n";
                MyFile << "Solve time: " << solve_time << "\n";
                MyFile << "\n";
            }
        }
    }
    MyFile.close();

    ierr = 0;
    FORTRAN_close(&funit, &ierr);
    if (ierr) {
        printf("Error closing %s on unit %d.\n", fname, (int)funit);
        printf("Trying not to abort.\n");
    }

    FREE(pname);
    FREE(x); FREE(bl); FREE(bu);
    FREE(v); FREE(cl); FREE(cu);
    FREE(equatn);
    FREE(linear);

    if (constrained)
        CUTEST_cterminate_r(&status);
    else
        CUTEST_uterminate_r(&status);

    return 0;
}
#ifdef __cplusplus
}
#endif

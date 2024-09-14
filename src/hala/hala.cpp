// runcutest -p hala -D HS3 -lstdc++
#include "stdio.h"
#include <iostream>
#include <Eigen/Core>
#include <LBFGSB.h>
#ifdef __cplusplus
extern "C" {    /*To prevent C++ compilers from mangling symbols*/ 
#endif

#include <math.h>
#include "cutest.h"


using Eigen::VectorXd;
using namespace LBFGSpp;
using namespace std;

//Function declaration
double CalculateHL(double*, double, int, double, double*, int*);
void CalculateHLGr(double*, double, int, int, double*, double*, double*, double*, int*);
void UpdateLagr(int, double, double*, double*, int*);
bool CheckFeasible(int, double, double*, double*);




class Evaluation{
    private:
        int n;
        int m;
        int* constraintType;
        double* tau;
        double* lambda;

    public:
        Evaluation(int n_, int m_, int* constraintType_, double* tau_, double* lambda_) : 
        n(n_), m(m_), constraintType(constraintType_), tau(tau_), lambda(lambda_){}
        //Function to evaluate Lagrangian and gradient
        double operator()(const VectorXd& x, VectorXd& grad){

            double constraints[m];
            double objectiveGradient[n];
            double constraintGradient[n*m];
            double gradResult[n];
            int status;
            double f;
            double x0[n];

            for (int j = 0; j < n; j++){
                x0[j] = x(j);
            }
            //Evaluate f(x) and g(x)
            CUTEST_cfn(&status, &n, &m, x0, &f, constraints);
            
            logical grglagf = FALSE_;
            logical jtrans = FALSE_;

            //Evaluate f'(x) and g'(x)
            CUTEST_cgr(&status, &n, &m, x0, NULL, &grglagf, objectiveGradient, &jtrans, &m, &n, constraintGradient);

            //Calculate Hyperbolic Lagrangian and its gradient
            double fx = CalculateHL(lambda, *tau, m, f, constraints, constraintType); 
            CalculateHLGr(lambda, *tau, m, n, objectiveGradient, constraintGradient, constraints, gradResult, constraintType);
            for (int j = 0; j < n; j++){
                grad[j] = gradResult[j];
            }


            return fx;
        }
};


rp_ hala(int n, int m, double* x, int maxIterations, double* bl, double* bu,
 double tau, double* cl, double* cu, double tolerance, double alpha, double lambda0, int* currentIteration, int* constraintType){

    //Memory allocation for arrays
    double f;
    double* newX = new double[n];
    double* lambda = new double[m];
    double* constraints = new double[m];
    double* objectiveGradient = new double[n];
    double* constraintGradient = new double[n*m];
    double* gradResult = new double[n];
    double* nullVector = new double[n];
    double fx;

    nullVector[n] = {0.0};

    //Set parameters for the LBFGS-B algorithm
    LBFGSBParam<double> param;
    param.max_iterations = 10000;
    LBFGSBSolver<double> solver(param);
    Evaluation fun(n, m, constraintType, &tau, lambda);

    //Set initial values for lambda
    for (int i = 0; i < m; i++){
        lambda[i] = lambda0;
    }
    VectorXd oldEigenX = VectorXd::Map(x, n);
    VectorXd newEigenX = VectorXd::Map(x, n);
    VectorXd lowerBound = VectorXd::Map(bl, n);
    VectorXd upperBound = VectorXd::Map(bu, n);
    
    int status;
    *currentIteration = 0;
    //Begin main loop
    for (int k = 0; k < maxIterations; k++){
        (*currentIteration)++;
        
        oldEigenX = newEigenX;

        int numberInnerIterations = solver.minimize(fun, newEigenX, fx, lowerBound, upperBound);

        copy(newEigenX.data(), newEigenX.data() + n, newX);

        CUTEST_cfn(&status, &n, &m, newX, &f, constraints);
        UpdateLagr(m, tau, lambda, constraints, constraintType);
        bool feasability = CheckFeasible(m, tolerance, lambda, constraints);
        bool xConvergence = oldEigenX.isApprox(newEigenX, tolerance);

        if (feasability && xConvergence) break;
        if (!feasability) tau = tau * alpha;

    }
    delete[] newX;
    delete[] lambda;
    delete[] constraints;
    delete[] objectiveGradient;
    delete[] constraintGradient;
    delete[] gradResult;

    return f;

}

double CalculateHL(double* lambda, double tau, int m, double f, double* constraints, int* constraintType){

    //Set hyperbolic penalty
    double Penalty = 0.0;
    for (int i = 0; i < m; i++){
        double gamma = lambda[i] * constraints[i];
        if (constraintType[i] == 0){
            Penalty += -gamma + sqrt((gamma * gamma) + (1/(tau * tau)));
    }
        else{
            Penalty += gamma + sqrt((gamma * gamma) + (1/(tau * tau)));
        }
    }
    return f + Penalty;
}


void CalculateHLGr(double* lambda, double tau, int m, int n, double* objectiveGradient, 
    double* constraintGradient, double* constraints, double* gradResult, int* constraintType){

    //Set penalty
    for (int i = 0; i < n; i++){
        double Penalty = 0.0;
        for (int j = 0; j < m; j++){
            double gamma = lambda[j] * constraints[j];
            if (constraintType[j] == 0){
                Penalty -= lambda[j] * (1 - (gamma / sqrt((gamma * gamma) + (1/(tau * tau))))) * constraintGradient[i * m + j];
        }
            else{
                Penalty += lambda[j] * (1 + (gamma / sqrt((gamma * gamma) + (1/(tau * tau))))) * constraintGradient[i * m + j];
            }
        }
    gradResult[i] = objectiveGradient[i] + Penalty;
    }
}


void UpdateLagr(int m, double tau, double* lambda, double* constraints, int* constraintType){

    for (int i = 0; i < m; i++){
        double gamma = lambda[i] * constraints[i];
        if (constraintType[i] == 0){
            lambda[i] = lambda[i] * (1 - ((gamma)/ sqrt((gamma * gamma) + (1/(tau * tau)))));
        }
        else{
            lambda[i] = lambda[i] * (1 + ((gamma)/ sqrt((gamma * gamma) + (1/(tau * tau)))));
        }
    }
}


/*bool CheckFeasible(int m, double tolerance, double* lambda, double* constraints, double* cl, double* cu){

    for (int i = 0; i < m; i++){
        if (constraints[i] < cl[i] - tolerance || constraints[i] > cu[i] + tolerance) return false;
    }
    return true;
}*/

bool CheckFeasible(int m, double tolerance, double* lambda, double* constraints){

    for (int i = 0; i < m; i++){
        if ((lambda[i] * constraints[i]) > tolerance || (lambda[i] * constraints[i]) < - tolerance) return false; 
    }
    return true;
    
}


void halaspc( integer funit, char *fname )
{

    integer ierr;

    /* This is a dummy routine to read a spec file.
       Possibly, this routine contains precision-dependent directives */

    /* Open relevant file */
    FORTRAN_open( &funit, fname, &ierr );
    if ( ierr )
    {
        printf( "Error opening spec file %s.\nAborting.\n", fname );
        exit(1);
    }

    /* ... Do something ... */

    FORTRAN_close( &funit, &ierr );
    return;

}

void getinfo( integer n, integer m, rp_ *bl, rp_ *bu,
              rp_ *cl, rp_ *cu, logical *equatn,
              logical *linear, VarTypes *vartypes )
{

    int i;

    vartypes->nlin = 0; vartypes->neq = 0; vartypes->nbnds = 0;
    vartypes->nrange = 0;
    vartypes->nlower = 0; vartypes->nupper = 0; vartypes->nineq = 0;
    vartypes->nineq_lin = 0; vartypes->nineq_nlin = 0;
    vartypes->neq_lin = 0; vartypes->neq_nlin = 0;

    for ( i = 0; i < n; i++ )
        if ( bl[i] > -CUTE_INF || bu[i] < CUTE_INF ) vartypes->nbnds++;
    for ( i = 0; i < m; i++ )
    {
        if ( linear[i] ) vartypes->nlin++;
        if ( equatn[i] )
        {
            vartypes->neq++;
            if ( linear[i] )
                vartypes->neq_lin++;
            else
                vartypes->neq_nlin++;
        }
        else
        {
            vartypes->nineq++;
            if ( cl[i] > -CUTE_INF )
            {
                if ( cu[i] < CUTE_INF )
                    vartypes->nrange++;
                else
                {
                    vartypes->nlower++;
                }
            }
            else
            {
                if ( cu[i] < CUTE_INF )
                {
                    vartypes->nupper++;
                }
            }
            if ( linear[i] )
            {
                vartypes->nineq_lin++;
            }
            else
            {
                vartypes->nineq_nlin++;
            }
        }
    }
    return;
}

#ifdef __cplusplus
}     /*Closing brace for  extern "C"  block*/ 
#endif

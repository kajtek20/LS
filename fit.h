#include<iostream>
#include<cmath>
#include<vector>



using namespace std;




class fit_sines
{
public:
    fit_sines();// set all parameters
    void lin_fit_A_and_ph(const vector<long double>  & date, const vector<long double>  & flux, const vector<long double>  & weights,\
                          int n_sines, int n_datapoints, vector<vector<long double> > & sine_parameters);
    void Levenberg_Marquardt_fit(const vector<long double> & t, const vector<long double> & y, const vector<long double> & w,\
                                    const int n_sines, const int n_datapoints,\
                                    const vector<vector<bool> > fit_control, vector<vector<long double> > & sine_parameters);
    
    //moze warto zrobic konstruktor dla tych parametrow
    //parameters for L-M
    long double lambda0; //initial value of L-M parameter 
    long double eps1; //convergence tolerance for gradient
    long double eps2; //convergence tolerance for parameters
    long double eps3; //convergence tolerance for red. Chi-square
    long double eps4; //determines acceptance of a L-M step
    long double eps5; // determines max absolute change of frequency
    int maxiter_multipler;  //maxiter=maxiter_multipler *  Npar
    
private:
    //solve system of linear equations (Gauss elimination). Matrix A is replaced by its LU decomposition
    //Ax=b
    void gauss(int unknows, vector<vector<long double> > &A, vector<long double> &b, vector<long double> &x);
    //replace A with its inverse matrix
    void gauss_inverse_matrix(int size, vector<vector<long double> > &A);
    void make_J_and_ymod(const int n_datapoints, const vector<long double> & t,\
                               const int n_sines, const vector<vector<long double> > & sine_parameters,\
                               const vector<vector<bool> > fit_control, vector<vector<long double> > &J, vector<long double> & ymod);
    
    void err_Levenberg_Marquardt(const vector<long double> & t, const vector<long double> & y, const vector<long double> & w,\
                                        const int n_sines, const int n_datapoints,\
                                        const vector<vector<bool> > fit_control, vector<vector<long double> > & sine_parameters,\
                                        const int fitted_parameters, vector<vector<long double> > & JTWJ);
    
    
};

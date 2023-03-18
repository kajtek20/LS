//g++ -Wall -O2 -finline-functions    fast_LS.cpp   -lnfft3 -lfftw3 -lm -lc



/*
 * Based on Leroy, A&A 545, A50 (2012)
 */

#include <fstream>
#include <complex>
#include <nfft3.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>


using namespace std;



class light_curve
{
public:
    vector<long double> date, flux, flux_err, flux_resid, weights;
    int data_points;
    int n_sines;
    long double t0; // t0 i set to the floor from first observation
    long double Rayleigh_resolution;
    vector<vector<bool> > fit_control; // control whether to fit a given parameter for a given frequency, in order as in sine_parameters
    
    light_curve(); //konstruktor
    vector<vector<long double> > sine_parameters; // in vector are stored parameters of nonlinear fitted sines (offset + sum_{i=1}^n A*sin(2PI(nu t+phase))), ie.
                                             //sine_parameters[0][0] -offset, sine_parameters[0][1] -offset error
                                             //sine_parameters[1][0] - first frequency, sine_parameters[1][1] - amplitude of the first frequency, 
                                             //sine_parameters[1][2] - phase of the first frequency, sine_parameters[1][3] - frequency error, 
                                             //sine_parameters[1][4] - amplitude error, sine_parameters[1][5] - phase error, 
                                             //sine_parameters[1][6] - S/N calculated in window sine_parameters[1][7] - S/N calculated from full spectrum
                                             //sine_parameters[2][0] - second frequency etc
    
    void read_data(string file_name);
    
    void prewithen_data();
    
    void print();
private:
    int count_columns(string file_name);
    int columns;
};




class fast_LS_periodogram
{
public:
    /* Computes the Lomb - Scargle normalised periodogram of a times - series .
    t the times , reduced to [ -1/2 ,1/2).    y the measurements , centred around <y >.
    npts the length of the times -series .     over the oversampling factor .
    hifac the highest frequency in units of " average " Nyquist frequency .
    This function returns the results in a structure , LS ( see text ).    */
    void calculate_LS(const vector<long double> & t, const vector<long double> & y, const int npts , const long double over=4 , const long double hifac=1, bool save_trf=true);
    void find_nu_max_amp(long double StoNlimit=4, bool noise_in_window=true, long double window_size=1);
    fast_LS_periodogram();
    
    vector<long double>  freqs; // (>0) frequencies
    vector<long double> Pn;  // periodogram ordinates (amplitude periodogram)
    vector<long double> Pn_window;
    vector<long double> Pn_noise;
    long double noise_full;
    int nfreqs ; // number of frequencies
    long double var, mean;
    long double nu_max_amp, SN_nu_max_amp, SN_nu_max_amp_full, amp_max;
    
private:
    void meanAndVariance(const int npts, const vector<long double> &y);
    void centerData(const int npts, const vector<long double> & y);
    void reduce(const vector<long double> & t, const int npts, const long double oversampling);
    void nfft ( const vector<long double> & t, const vector<long double> & y, int n, int m, vector< complex<long double> > &d, bool window_spectrum);
    void calculate_noise(bool noise_in_window, long double window_size);
    
    vector<long double> y_centered, t_reduced;
    vector< complex<long double> > sp, win;
    int output_file_number=0;
    int noise_window_points_div_2;
    bool calculate_nfft_win=true;
    bool cal_noise_window_points=true;
    bool save_trf_and_noise_files=false;
};





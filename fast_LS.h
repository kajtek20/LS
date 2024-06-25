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
#include "interface.h"

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
                                             
    vector<vector<int> > v_komb;  // 1 kolumna 0 -> to czestotliwosc niezalezna, jeśli 1  to harmonika, 2 to kombinacja 2 rodziców, 
                                  //3 to kombinacja 3 rodziców taka, ze:
                                  // nu_i = v_komb[i][1] * nu_{v_komb[i][2]} + v_komb[i][3] * nu_{v_komb[i][4]} + v_komb[i][5] * nu_{v_komb[i][6]}
                                  // indices corresponds to sine_parameters therefore start fill from 1
                                             
    vector<long double> freq_error_Czerny; //Schwarzenberg-Czerny 1991
    
    void read_data(interface & iface);
    
    void prewithen_data();
    
    void err_corr_Czerny();  //Schwarzenberg-Czerny 1991
    
    //two overloaded functions to search for harmonics and combinations
    void check_har_kom(int i,  interface & iface); // check i-th frequency; suggest possible harmonics and combinations; 
                                                   //user decides; depending on the response, it sets the frequency to be fitted
                                                   //or set exact combination or harmonic;
                                                   //parents can be all frequencies found up to i-th frequency
    void check_har_kom(interface & iface); // check after finding all significant frequencies;
                                           // only add information about possible frequency;
                                           // don't check frequencies marked as harmonic or combination for the other version of function
                                           // as parents are considered only frequencies with higher amplitude than possible combination 
    
    
    void print();
    void remove_close_freq(long double factor, ofstream & out);
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
    void calculate_LS(const vector<long double> & t, const vector<long double> & y, const int npts , const long double over=4 , const long double hifac=1,\
                      bool save_trf=true, const bool spec_mode_default=true, const long double spec_max_freq=25, const long double spec_resol=1e-4);
    void find_nu_max_amp(long double StoNlimit, bool noise_in_window=true, long double window_size=1, int n_highest_peaks=1);
    //void find_nu_max_amp(int n_period,  double Rayleigh_resolution, double factor, long double StoNlimit, bool noise_in_window=true, long double window_size=1, int n_highest_peaks=1);
    void set_Pn_mask(int n_period, int index_nu_max_amp_scalar, double nu_max_amp_scalar, double Rayleigh_resolution, double factor);
    fast_LS_periodogram();
    void fill_nu_user_max_amp(long double user_frequecy, int accept_sol);
    
    vector<long double>  freqs; // (>0) frequencies
    vector<long double> Pn;  // periodogram ordinates (amplitude periodogram)
    vector<long double> Pn_window;
    vector<long double> Pn_noise;
    long double noise_full;
    int nfreqs ; // number of frequencies
    long double var, mean;
    vector<long double> nu_max_amp, SN_nu_max_amp, SN_nu_max_amp_full, amp_max;
    vector<int> index_nu_max_amp;

    
private:
    void meanAndVariance(const int npts, const vector<long double> &y);
    void centerData(const int npts, const vector<long double> & y);
    void reduce(const vector<long double> & t, const int npts, const long double oversampling);
    void nfft ( const vector<long double> & t, const vector<long double> & y, int n, int m, vector< complex<long double> > &d, bool window_spectrum);
    void calculate_noise(bool noise_in_window, long double window_size);
    
    vector<long double> y_centered, t_reduced;
    vector< complex<long double> > sp, win;
    vector<bool> Pn_mask;

    int output_file_number=0;
    int noise_window_points_div_2;
    bool calculate_nfft_win=true;
    bool cal_noise_window_points=true;
    bool save_trf_and_noise_files=false;
};





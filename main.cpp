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
//#include <variant>
#include "interface.h"
#include "fast_LS.h"
#include "fit.h"


using namespace std;




int main(int argc, char* argv[])
{
    ofstream out;
    light_curve lc;
    fast_LS_periodogram LS;
    fit_sines fs;
    
    interface iface(argc, argv);

    lc.read_data(iface.cp["data_file"]);
    out.open(iface.cp["freq_file"]);
    
    
    
    time_t start, koniec;
    
    lc.sine_parameters.push_back(vector<long double>(2));//vector for shift and its uncertainty
    lc.fit_control.push_back(vector<bool>(1, true));
    
    
    for(int i=1; ; i++)
    {
        if(string_to_bool(iface.cp["spec_mode_default"]))
        {
            if(i==1)
                LS.calculate_LS(lc.date, lc.flux, lc.data_points , 4 , 1, string_to_bool(iface.cp["save_spectra"]));
            else
                LS.calculate_LS(lc.date, lc.flux_resid, lc.data_points , 4 , 1, string_to_bool(iface.cp["save_spectra"]));
        }
        else
        {
            if(i==1)
                LS.calculate_LS( lc.date, lc.flux, lc.data_points , 4 , 1, string_to_bool(iface.cp["save_spectra"]),\
                                 string_to_bool(iface.cp["spec_mode_default"]), stod(iface.cp["spec_high_freq"]), stod(iface.cp["spec_res"]) );
            else
                LS.calculate_LS( lc.date, lc.flux_resid, lc.data_points , 4 , 1, string_to_bool(iface.cp["save_spectra"]),\
                                 string_to_bool(iface.cp["spec_mode_default"]), stod(iface.cp["spec_high_freq"]), stod(iface.cp["spec_res"]) );
        }
        
        LS.find_nu_max_amp(stod(iface.cp["StoN_criterion"]), string_to_bool(iface.cp["noise_in_window"]), stod(iface.cp["window_size"]));
        cout<<"Z LS: nu-> "<<LS.nu_max_amp<<" amp-> "<<LS.amp_max<<" S/N-> "<<LS.SN_nu_max_amp<<" S/Nfull-> "<<LS.SN_nu_max_amp_full<<endl;
        
        if(LS.nu_max_amp>0)
        {
            lc.sine_parameters.push_back(vector<long double>(8));
            lc.fit_control.push_back(vector<bool>(3, true));
            
            if(string_to_bool(iface.cp["fit_all_freq"]) == false)
            {
                for(int kk=1; kk<=i- stoi(iface.cp["fit_n_last_freq"]); kk++)
                    lc.fit_control[kk][0]=false;
            }
            
            lc.sine_parameters[i][6]=LS.SN_nu_max_amp;
            lc.sine_parameters[i][7]=LS.SN_nu_max_amp_full;
            
            lc.sine_parameters[i][0]=LS.nu_max_amp;
            lc.n_sines=i;
            fs.lin_fit_A_and_ph(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.sine_parameters);
            
            time( & start );
            fs.Levenberg_Marquardt_fit(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.fit_control, lc.sine_parameters);
            time( & koniec );      
            cout<<"CZAS M-L  "<<difftime( koniec, start )<<endl<<endl;
            
            lc.prewithen_data();
            
            for(int j=1; j<=i; j++){
            out<<fixed<<setw(4)<<j<<setw(12)<<setprecision(8)<<lc.sine_parameters[j][0]<<setw(12)<<setprecision(8)<<lc.sine_parameters[j][3];
            out<<setw(12)<<setprecision(4)<<lc.sine_parameters[j][1]<<setw(12)<<setprecision(4)<<lc.sine_parameters[j][4];
            out<<setw(12)<<setprecision(4)<<lc.sine_parameters[j][2]<<setw(12)<<setprecision(4)<<lc.sine_parameters[j][5];
            out<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][6]<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][7]<<endl;}
            out<<endl<<endl;
        }
        else
            break;
    }
    out.close();
}



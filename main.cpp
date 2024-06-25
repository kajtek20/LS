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
    ofstream out, out2;
    light_curve lc;
    fast_LS_periodogram LS;
    fit_sines fs;
    string command;
    
    interface iface(argc, argv);

    lc.read_data(iface);
    out.open(iface.cp["freq_file"]+"_history");
    out2.open(iface.cp["freq_file"]);
    
    
    
    setenv("OMP_NUM_THREADS", iface.cp["threads"].c_str(),1);
    
    time_t start_f, koniec_f;
    
    time( & start_f );
    
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

        LS.find_nu_max_amp(stod(iface.cp["StoN_criterion"]), string_to_bool(iface.cp["noise_in_window"]), stod(iface.cp["window_size"]), stoi(iface.cp["try_n_highest_peaks"]));

        //if( string_to_bool(iface.cp["remove_too_close_hard"]) )
        //    LS.find_nu_max_amp(i,  lc.Rayleigh_resolution, stod(iface.cp["remove_too_close_fac"]), stod(iface.cp["StoN_criterion"]), string_to_bool(iface.cp["noise_in_window"]), stod(iface.cp["window_size"]));
        //else
        //    LS.find_nu_max_amp(stod(iface.cp["StoN_criterion"]), string_to_bool(iface.cp["noise_in_window"]), stod(iface.cp["window_size"]));


        if(stoi(iface.cp["set_try_n_highest_peaks_eq_1_after_found_n_freq"]) == i)
            iface.cp["try_n_highest_peaks"]="1";
        
        if(LS.nu_max_amp[0]>0)
        {
            int accept_sol=-1;
            double user_frequecy;
            if(stoi(iface.cp["try_n_highest_peaks"]) == 1)
            {
                accept_sol=1;
                cout<<"\E[1;34m" <<"####################################################################"<< "\E[0m"<<endl;
                cout<<"L-S: nu"<<i<<"-> " << "\E[1;31m" << LS.nu_max_amp[0] << "\E[0m"<<" amp-> "<<LS.amp_max[0]<<" S/N-> "<<LS.SN_nu_max_amp[0]<<" S/Nfull-> "<<LS.SN_nu_max_amp_full[0]<<endl;

                lc.sine_parameters.push_back(vector<long double>(8));
                lc.fit_control.push_back(vector<bool>(3, true));
            
                if(string_to_bool(iface.cp["fit_all_freq"]) == false)
                {
                    for(int kk=1; kk<=i- stoi(iface.cp["fit_n_last_freq"]); kk++)
                        lc.fit_control[kk][0]=false;
                }
            
                lc.sine_parameters[i][6]=LS.SN_nu_max_amp[0];
                lc.sine_parameters[i][7]=LS.SN_nu_max_amp_full[0];
            
                lc.sine_parameters[i][0]=LS.nu_max_amp[0];
                lc.n_sines=i;

                //check combinations and harmonics
                if( string_to_bool(iface.cp["set_exact_com_and_har"]) )
                {
                    lc.check_har_kom(i, iface); // suggest possible harmonics and combinations; user decides; depending on the response, it sets the frequency to be fitted or set exact combination or harmonic
                }
            
                fs.lin_fit_A_and_ph(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.sine_parameters);
            
                fs.Levenberg_Marquardt_fit(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.fit_control, lc.sine_parameters, lc.v_komb);
            }
            else
            {
                for(int try_n=0; try_n<stoi(iface.cp["try_n_highest_peaks"]); try_n++)
                {
                    cout<<"\E[1;34m" <<"####################################################################"<< "\E[0m"<<endl;
                    cout<< "\E[22;32m"<<"Trail solution nr "<<try_n+1<< "\E[0m"<<endl;
                    cout<<"L-S: nu"<<i<<"-> " << "\E[1;31m" << LS.nu_max_amp[try_n] << "\E[0m"<<" amp-> "<<LS.amp_max[try_n]<<" S/N-> "<<LS.SN_nu_max_amp[try_n]<<" S/Nfull-> "<<LS.SN_nu_max_amp_full[try_n]<<endl;

                    light_curve lc_trail;
                    lc_trail=lc;

                    lc_trail.sine_parameters.push_back(vector<long double>(8));
                    lc_trail.fit_control.push_back(vector<bool>(3, true));

                    if(string_to_bool(iface.cp["fit_all_freq"]) == false)
                    {
                        for(int kk=1; kk<=i- stoi(iface.cp["fit_n_last_freq"]); kk++)
                            lc_trail.fit_control[kk][0]=false;
                    }

                    lc_trail.sine_parameters[i][6]=LS.SN_nu_max_amp[try_n];
                    lc_trail.sine_parameters[i][7]=LS.SN_nu_max_amp_full[try_n];

                    lc_trail.sine_parameters[i][0]=LS.nu_max_amp[try_n];
                    lc_trail.n_sines=i;

                    //check combinations and harmonics
                    if( string_to_bool(iface.cp["set_exact_com_and_har"]) )
                    {
                        lc_trail.check_har_kom(i, iface); // suggest possible harmonics and combinations; user decides; depending on the response, it sets the frequency to be fitted or set exact combination or harmonic
                    }

                    fs.lin_fit_A_and_ph(lc_trail.date, lc_trail.flux, lc_trail.weights, lc_trail.n_sines, lc_trail.data_points, lc_trail.sine_parameters);
                    fs.Levenberg_Marquardt_fit(lc_trail.date, lc_trail.flux, lc_trail.weights, lc_trail.n_sines, lc_trail.data_points, lc_trail.fit_control, lc_trail.sine_parameters, lc_trail.v_komb);
                }
                cout<<endl<<"Provide the number of the accepted solution"<<endl;
                cout<<"or enter -100 to specify any frequency"<<endl;
                cin>>accept_sol;
                if(accept_sol==-100)
                {
                    accept_sol=1;
                    cout<<"Enter frequency to fit"<<endl;
                    cin>>user_frequecy;

                    //LS.SN_nu_max_amp[accept_sol-1]=0;
                    //LS.SN_nu_max_amp_full[accept_sol-1]=0;
                    //LS.nu_max_amp[accept_sol-1]=user_frequecy;
                    LS.fill_nu_user_max_amp(user_frequecy, accept_sol);
                }
                lc.sine_parameters.push_back(vector<long double>(8));
                lc.fit_control.push_back(vector<bool>(3, true));





                if(string_to_bool(iface.cp["fit_all_freq"]) == false)
                {
                    for(int kk=1; kk<=i- stoi(iface.cp["fit_n_last_freq"]); kk++)
                        lc.fit_control[kk][0]=false;
                }

                lc.sine_parameters[i][6]=LS.SN_nu_max_amp[accept_sol-1];
                lc.sine_parameters[i][7]=LS.SN_nu_max_amp_full[accept_sol-1];

                lc.sine_parameters[i][0]=LS.nu_max_amp[accept_sol-1];
                lc.n_sines=i;

                //check combinations and harmonics
                if( string_to_bool(iface.cp["set_exact_com_and_har"]) )
                {
                    lc.check_har_kom(i, iface); // suggest possible harmonics and combinations; user decides; depending on the response, it sets the frequency to be fitted or set exact combination or harmonic
                }

                fs.lin_fit_A_and_ph(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.sine_parameters);
                fs.Levenberg_Marquardt_fit(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.fit_control, lc.sine_parameters, lc.v_komb);
            }

            if( string_to_bool(iface.cp["remove_too_close_hard"]) )
            {
                LS.set_Pn_mask(i, LS.index_nu_max_amp[accept_sol-1], LS.nu_max_amp[accept_sol-1], lc.Rayleigh_resolution, stod(iface.cp["remove_too_close_fac"]));
            }
            cout<<endl<<endl<<endl;
             
            lc.prewithen_data();
            
            for(int j=1; j<=i; j++){
            out<<fixed<<setw(4)<<j<<setw(14)<<setprecision(8)<<lc.sine_parameters[j][0]<<setw(12)<<setprecision(8)<<lc.sine_parameters[j][3];
            out<<setw(12)<<setprecision(4)<<lc.sine_parameters[j][1]<<setw(12)<<setprecision(4)<<lc.sine_parameters[j][4];
            out<<setw(12)<<setprecision(6)<<lc.sine_parameters[j][2]<<setw(12)<<setprecision(6)<<lc.sine_parameters[j][5];
            out<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][6]<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][7]<<endl;}
            out<<endl<<endl;
        }
        else
            break;
    }
    
    if( string_to_bool(iface.cp["remove_too_close"]) )
    {
        ofstream out3;
        out3.open(iface.cp["freq_file"]+"_removed");
        lc.remove_close_freq(stod(iface.cp["remove_too_close_fac"]), out3);
        out3.close();
        fs.Levenberg_Marquardt_fit(lc.date, lc.flux, lc.weights, lc.n_sines, lc.data_points, lc.fit_control, lc.sine_parameters, lc.v_komb);
    }
    
    lc.check_har_kom(iface); //overloaded function, only check existance of possible harmonics or combinations
    
    if( string_to_bool(iface.cp["Czerny_err_corr"]) )
    {
        lc.err_corr_Czerny();
    
        out2<<"  ID        nu         nu_err nu_err_Czerny    A            A_err     phase[0-1]  phase_err S/N window S/N full    remarks"<<endl;
        for(int j=1; j<=lc.n_sines; j++){
            out2<<fixed<<setw(4)<<j<<setw(14)<<setprecision(8)<<lc.sine_parameters[j][0]<<setw(12)<<setprecision(8)<<lc.sine_parameters[j][3];
            out2<<setw(12)<<setprecision(8)<<lc.freq_error_Czerny[j];
            out2<<setw(13)<<setprecision(5)<<lc.sine_parameters[j][1]<<setw(13)<<setprecision(5)<<lc.sine_parameters[j][4];
            out2<<setw(12)<<setprecision(6)<<lc.sine_parameters[j][2]<<setw(12)<<setprecision(6)<<lc.sine_parameters[j][5];
            out2<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][6]<<setw(10)<<setprecision(2)<<lc.sine_parameters[j][7];
            
            if(lc.v_komb[j][0] == 0 || lc.v_komb[j][0] == 4)
                out2<<setw(10)<<setprecision(3)<<0<<" n * ";
            if(lc.v_komb[j][0] == 1 || lc.v_komb[j][0] == -1)
                out2<<setw(10)<<setprecision(3)<<lc.v_komb[j][0]<<" h * "<<lc.v_komb[j][1]<<"*f"<<lc.v_komb[j][2];
            if(lc.v_komb[j][0] == 2 || lc.v_komb[j][0] == -2)
                out2<<setw(10)<<setprecision(3)<<lc.v_komb[j][0]<<" k * "<<lc.v_komb[j][1]<<"*f"<<lc.v_komb[j][2]\
                    <<showpos<<lc.v_komb[j][3]<<noshowpos<<"*f"<<lc.v_komb[j][4];
            if(lc.v_komb[j][0] == 3 || lc.v_komb[j][0] == -3)
                out2<<setw(10)<<setprecision(3)<<lc.v_komb[j][0]<<" k * "<<lc.v_komb[j][1]<<"*f"<<lc.v_komb[j][2]\
                    <<showpos<<lc.v_komb[j][3]<<noshowpos<<"*f"<<lc.v_komb[j][4]\
                    <<showpos<<lc.v_komb[j][5]<<noshowpos<<"*f"<<lc.v_komb[j][6];
        
            
            out2<<endl;
        }
    }
    else
    {
        out2<<"  ID        nu         nu_err     A            A_err     phase[0-1]  phase_err S/N window S/N full    remarks"<<endl;
        for(int j=1; j<=lc.n_sines; j++){
            out2<<fixed<<setw(4)<<j<<setw(14)<<setprecision(8)<<lc.sine_parameters[j][0]<<setw(12)<<setprecision(8)<<lc.sine_parameters[j][3];
            out2<<setw(13)<<setprecision(5)<<lc.sine_parameters[j][1]<<setw(13)<<setprecision(5)<<lc.sine_parameters[j][4];
            out2<<setw(12)<<setprecision(6)<<lc.sine_parameters[j][2]<<setw(12)<<setprecision(6)<<lc.sine_parameters[j][5];
            out2<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][6]<<setw(8)<<setprecision(2)<<lc.sine_parameters[j][7];
            
            if(lc.v_komb[j][0] == 0)
                out2<<setw(10)<<setprecision(3)<<0<<" n * ";
            if(lc.v_komb[j][0] == 1 || lc.v_komb[j][0] == -1)
                out2<<setw(10)<<setprecision(3)<<lc.v_komb[j][0]<<" h * "<<lc.v_komb[j][1]<<"*f"<<lc.v_komb[j][2];
            if(lc.v_komb[j][0] == 2 || lc.v_komb[j][0] == -2)
                out2<<setw(10)<<setprecision(3)<<lc.v_komb[j][0]<<" k * "<<lc.v_komb[j][1]<<"*f"<<lc.v_komb[j][2]\
                    <<showpos<<lc.v_komb[j][3]<<noshowpos<<"*f"<<lc.v_komb[j][4];
            if(lc.v_komb[j][0] == 3 || lc.v_komb[j][0] == -3)
                out2<<setw(10)<<setprecision(3)<<lc.v_komb[j][0]<<" k * "<<lc.v_komb[j][1]<<"*f"<<lc.v_komb[j][2]\
                    <<showpos<<lc.v_komb[j][3]<<noshowpos<<"*f"<<lc.v_komb[j][4]\
                    <<showpos<<lc.v_komb[j][5]<<noshowpos<<"*f"<<lc.v_komb[j][6];
            
            out2<<endl;
        }
    }

    
    
    out.close(); out2.close();
    
    time( & koniec_f );
    cout<<"Calculation time:  "<<difftime( koniec_f, start_f )<<endl<<endl;
}



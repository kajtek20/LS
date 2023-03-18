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
#include "fast_LS.h"

using namespace std;


template <typename T> int sgn(T val)  // szablon funkcji signum
{
    return (T(0) < val) - (val < T(0));
}

inline long double sign(long double a, long double b)
{
  return ((b >= 0) ? 1 : -1) * abs(a);
}

inline long double square(long double a)
{
  return a * a; 
}






int light_curve::count_columns(string file_name)
{
    int col=0;
    string line, item;
    stringstream licznikStr;
    ifstream in(file_name);
    if(!in.good())
    {
        cout<<"Problem with "<<file_name<<" file"<<endl;
        exit(-1);
    }
    
    getline(in, line);
    licznikStr.str(line);
    while(licznikStr>>item) col++;
    in.close();
    return col;
}

void light_curve::print()
{
    cout<<"col "<<columns<<"  "<<data_points<<endl;
    
    for(int i=0; i<data_points; i++)
        cout<<fixed<<setw(20)<<setprecision(6)<<date[i]<<" "<<setw(20)<<setprecision(6)<<flux[i]<<endl;
}

void light_curve::read_data(string file_name)
{
    string tmp;
    long double t, f, ferr;
    long double mean=0;
    ifstream in(file_name);
    if(!in.good())
    {
        cout<<"Problem with "<<file_name<<" file"<<endl;
    }
    columns=count_columns(file_name);
    if(columns<2)
        cout<<"Ups! something wrong with file "<<file_name<<endl;
    
    data_points=0;
    
    
    if(columns>0)
    {
        while(true)
        {
            in>>t>>f;
            if(in.eof())
                break;
            
            mean+=f;
            if(columns>2)
            {
                in>>ferr;
                flux_err.push_back(ferr);
                weights.push_back(1/(ferr*ferr));
            }
            else
            {
                flux_err.push_back(1);
                weights.push_back(1);
            }
            
            date.push_back(t); flux.push_back(f);
            data_points++;
            
            for(int i=3; i<columns; i++)
                in>>tmp;
        }
    }
    
    t0=floor(date[0]);
    Rayleigh_resolution=1/(date[data_points-1]-date[0]);
    mean/=data_points;
    for(int i=0; i<data_points; i++)
    {
        flux[i]-=mean;
        date[i]-=t0;
    }

    cout<<"flux_mean "<<mean<<" "<<data_points<<"   T0 "<<date[0]<<endl;
        
    in.close();
}






void fast_LS_periodogram::meanAndVariance(const int npts, const vector<long double> &y)
{
    mean = 0;
    long double M2 = 0;

    int nn = 1;
    for (int i = 0; i < npts; i++, nn++)
    {
        long double delta = y[i] - mean;
        mean +=  delta / nn;
        M2 += delta * (y[i] - mean); // This expression uses the new value of mean
    }
    var = M2 / (npts - 1);
}

void fast_LS_periodogram::centerData(const int npts, const vector<long double> & y)
{
  meanAndVariance(npts, y);

  y_centered.clear();
  y_centered.resize(npts);
  for (int i = 0; i < npts; i++)
    y_centered[i]=(y[i] - mean);
}


/**
 * Reduces the times to the interval [-1/2, 1/2).
 * param[in,out] the times.
 * param[in] npts the number of time values.
 * param[in] oversampling the oversampling factor.
 */
void fast_LS_periodogram::reduce(const vector<long double> & t, const int npts, const long double oversampling)
{
  long double tmax = t[npts - 1];
  long double tmin = t[0];

  long double range = oversampling * (tmax - tmin);

  const long double eps_border = 1e-20;
  const long double a = 0.5 - eps_border;

  t_reduced.clear();
  t_reduced.resize(npts);
  // Reduce to [-1/2, 1/2[
  for (int i = 0; i < npts; i++)
    t_reduced[i] = 2 * a * (t[i] - tmin) / range - a;
}


void fast_LS_periodogram::nfft ( const vector<long double> & t, const vector<long double> & y, int n, int m, vector< complex<long double> > &d, bool window_spectrum)
{
    fftw_init_threads();
    
	d.resize(m+1);
    // Creates NFFT plan for 2*m Fourier coefficients ( positive and negative/ frequencies ) and n data samples .
	nfft_plan p;
    nfft_init_1d (&p, 2 * m, n);
	if( !window_spectrum ) // data spectrum
	{
		for (int i = 0; i < n; i++)
		{
			p.x[i] = t[i];
			p.f[i][0] = y[i];
			p.f[i][1] = 0;
		}
	}
	else  // window spectrum
	{
		for (int i = 0; i < n; i++)
		{
			p.x[i] = t[i];
            p.f[i][0] = 1.0;
			p.f[i][1] = 0;
		}
	}

   if (p.flags & PRE_ONE_PSI )   // Possibly optimises .
	   nfft_precompute_one_psi (&p);

   nfft_adjoint (&p);  // Computes the adjoint transform .

   for (int i = 0; i < m; i ++) // Outputs the positive frequency Fourier coefficients .
   {
	   d[i].real( p.f_hat[i+m][0] );
	   d[i].imag( p.f_hat[i+m][1] );
   }
   d[m].real(p.f_hat[0][0]); // = conj (p. f_hat [0]);
   d[m].imag(-p.f_hat[0][1]);
   nfft_finalize (&p);
}


void fast_LS_periodogram::calculate_LS(const vector<long double> & t, const vector<long double> & y, const int npts , const long double over , const long double hifac, bool save_trf)
{
    long double df = 1.0 / ( over * (t[npts - 1] - t[0]));
    // Index of the highest frequency in the positive frequency part of spectrum .
    int m = floor (0.5 * npts * over * hifac );
	
	// Centers the data.
	centerData(npts, y);
	// Reduces the times to [-1/2, 1/2).
    reduce(t, npts, over);
	
	freqs.resize(m);
    Pn.resize(m);
    nfreqs = m;
    
    
    //time_t start, koniec;
    //time( & start );
    // Unnormalised FTs of the data and window .
    //vector<complex<long double> > sp ;
	sp.resize(m+1);
    nfft(t_reduced, y_centered, npts , m, sp, false);
    
    
    //vector<complex<long double> > win;
    if(calculate_nfft_win)
    {
        calculate_nfft_win=false;
        win.resize(2*m+1);
	    vector<long double> empty_vector;
        nfft (t_reduced, empty_vector , npts , 2 * m, win, true);
    }
    //time( & koniec );      
    //cout<<"CZAS nfft  "<<difftime( koniec, start )<<endl;


    complex<long double>  z1, z2;
    long double absz2, hc2wt, hs2wt, cwt, swt, den, cterm, sterm;
    int jm1;
    // Computes the periodogram ordinates , and store the results in the LS structure .
	for (int j = 1; j <= m; j ++)
    //for (int j = 0; j <= m; j ++)
	{
		z1 = sp[j];
        // FT of data at \ omega
        z2 = win [2 * j]; // FT of window at 2\ omega
        absz2 = abs(z2);
        hc2wt = 0.5 * z2.imag() / absz2 ;
        hs2wt = 0.5 * z2.real() / absz2 ;
        cwt = sqrt (0.5 + hc2wt );
		swt = sign ( sqrt (0.5 - hc2wt), hs2wt );
        den = 0.5 * npts + hc2wt * z2.real() + hs2wt * z2.imag();
        cterm = square (cwt * z1.real() + swt * z1.imag( )) / den ;
        sterm = square (cwt * z1.imag() - swt * z1.real( )) / ( npts - den );
        jm1=j - 1;
        //jm1=j;
        freqs [jm1] = j * df;
        //cout<<j*df<<" "<<freqs [jm1]<<endl;
        Pn[jm1] = sqrt(2*( cterm + sterm) / (npts));
	}
	
	
	if(save_trf && output_file_number == 0)
    {
        Pn_window.resize(nfreqs);
        
//        for (int j = 2; j <= m; j ++)
            for (int j = 1; j <= m; j ++)
        {
            z1 = win[j];
            // FT of data at \ omega
            z2 = win [2 * j]; // FT of window at 2\ omega
            absz2 = abs(z2);
            hc2wt = 0.5 * z2.imag() / absz2 ;
            hs2wt = 0.5 * z2.real() / absz2 ;
            cwt = sqrt (0.5 + hc2wt );
		    swt = sign ( sqrt (0.5 - hc2wt), hs2wt );
            den = 0.5 * npts + hc2wt * z2.real() + hs2wt * z2.imag();
            cterm = square (cwt * z1.real() + swt * z1.imag( )) / den ;
            sterm = square (cwt * z1.imag() - swt * z1.real( )) / ( npts - den );
  //          jm1=j - 2;
                      jm1=j-1;
            freqs [jm1] = j * df;
            Pn_window[jm1] = sqrt(2*( cterm + sterm) / (npts));
        }
        
        
        ofstream out;
        out.clear();
        string trf_file_name="spec_win.trf";
        out.open(trf_file_name);
        if(!out.good())
            cout<<"Problem with output file"<<endl;
        
        //for some reasons Pn_window for the first frequency is not well calculated :(
        for(int i=1; i<nfreqs; i++)
            out<<fixed<<setw(12)<<setprecision(7)<<freqs[i]<<setw(14)<<setprecision(7)<<Pn_window[i]<<'\n';
        
        out.close();
    }
	
	if(save_trf)
    {
        save_trf_and_noise_files=true; // control for save noise files
        ofstream out;
        out.clear();
        string trf_file_name;
        stringstream licznikStr;
        output_file_number++;
        licznikStr.str("");
        licznikStr<<setw(4)<<setfill('0')<<output_file_number;
        trf_file_name="spec"+licznikStr.str()+".trf";
        out.open(trf_file_name);
        if(!out.good())
            cout<<"Problem with output file"<<endl;
        
        for(int i=0; i<nfreqs; i++)
            out<<fixed<<setw(12)<<setprecision(7)<<freqs[i]<<setw(14)<<setprecision(7)<<Pn[i]<<'\n';
        
        out.close();
    }
}


void fast_LS_periodogram::calculate_noise(bool noise_in_window, long double window_size)
{
    int noise_window_points=0;
    int noise_licznik=0;
    Pn_noise.resize(nfreqs);
    if(noise_in_window)
    {      
      if(cal_noise_window_points)
      {
          cal_noise_window_points=false; //it is enough calculate only once noise window points
          for(int i=0; i<nfreqs; i++)
          {
              if(freqs[i]-freqs[0] < window_size )
                  noise_window_points++;
              else
                  break;
          }
          noise_window_points_div_2=noise_window_points/2;
      }
      
      long double noise=0, noise_final=0;
      noise_licznik=0;
      for(int i=0; i<nfreqs; i++)
      {
          if(i==0)
          {
              for(int j=0; j<noise_window_points_div_2; j++)
              {
                  noise+=Pn[j];
              }
              noise_final=noise/noise_window_points_div_2;
              noise_licznik=noise_window_points_div_2;
          }
          else 
          {
              if(i<noise_window_points_div_2 && i<nfreqs-noise_window_points_div_2)
              {
                  noise+=Pn[i+noise_window_points_div_2-1];
                  noise_licznik++;
                  noise_final=noise/noise_licznik;
              }
              else if(i>=noise_window_points_div_2 && i<nfreqs-noise_window_points_div_2)
              {
                  noise+=Pn[i+noise_window_points_div_2-1];
                  noise-=Pn[i-noise_window_points_div_2];
                  
                  noise_final=noise/noise_licznik;
              }
              else
              {
                  noise-=Pn[i-noise_window_points_div_2];
                  noise_licznik--;
                  noise_final=noise/noise_licznik;
              }
          }
          Pn_noise[i]=noise_final;
          
          
      }
    }
    noise_licznik=0;
    noise_full=0;
    for(int i=0; i<nfreqs; i++)
    {
        noise_licznik++;
		noise_full+=Pn[i];
    }
    noise_full/=noise_licznik;
    if(save_trf_and_noise_files)
    {
        string trf_file_name;
        ofstream out;
        out.clear();
        string noise_file_name;
        stringstream licznikStr;
        licznikStr.str("");
        licznikStr<<setw(4)<<setfill('0')<<output_file_number;
        trf_file_name="noise"+licznikStr.str()+".trf";
        out.open(trf_file_name);
        if(!out.good())
            cout<<"Problem with output file"<<endl;
        
        for(int i=0; i<nfreqs; i++)
            out<<fixed<<setw(12)<<setprecision(7)<<freqs[i]<<setw(14)<<setprecision(7)<<Pn_noise[i]<<'\n';
        
        out.close();
    }
}



void fast_LS_periodogram::find_nu_max_amp(long double StoNlimit, bool noise_in_window, long double window_size)
{
    nu_max_amp=SN_nu_max_amp=amp_max=-1;
    
    calculate_noise(noise_in_window, window_size);
    for(int i=0; i<nfreqs; i++)
    {
        if( Pn[i] > amp_max && Pn[i] >= StoNlimit*Pn_noise[i] )
        {
            amp_max=Pn[i];
            nu_max_amp=freqs[i];
            SN_nu_max_amp=Pn[i]/Pn_noise[i];
            SN_nu_max_amp_full=Pn[i]/noise_full;
        }
    }
}



fast_LS_periodogram::fast_LS_periodogram()
{
    output_file_number=0;
    calculate_nfft_win=true;
    cal_noise_window_points=true;
    save_trf_and_noise_files=false;
}


//konstruktor
light_curve::light_curve()
{
    n_sines=0;
}



void light_curve::prewithen_data()
{
  long double arg, sum;
  long double dwapi=2*M_PI;
  flux_resid.resize(data_points);
  
  for(int i=0; i<data_points; i++)
  {
    sum=sine_parameters[0][0];
    for(int j=1; j<=n_sines; j++)
    {
      arg=dwapi*(sine_parameters[j][0]*date[i]+sine_parameters[j][2]);
      sum+=sine_parameters[j][1]*sin(arg);
    }
    flux_resid[i] = flux[i]-sum;
  }
}

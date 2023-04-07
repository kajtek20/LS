//Levenberg-Marquardt based on Gavin 2022


#include<iostream>
#include<cmath>
#include"fit.h"
#include<iomanip>
#include <ctime>

using namespace std;


void add_matrix(vector<vector<long double> > & m1, vector<vector<long double> > & m2)
{
    for(std::vector<std::vector<long double> >::size_type i=0; i<m1.size(); i++)
        for(std::vector<std::vector<long double> >::size_type j=0; j<m1.size(); j++)
            m1[i][j]+=m2[i][j];
}


const long double dwaPI=2*M_PI;


void fit_sines::lin_fit_A_and_ph(const vector<long double>  & date, const vector<long double>  & flux, const vector<long double>  & weights,\
                                 int n_sines, int n_datapoints, vector<vector<long double> > & sine_parameters)
{
    vector<vector<long double> > XTWX;
    vector<long double> XWy, beta, X_row;
    int parameters=n_sines*2+1;
    long double arg, X_row_j_wt;
    
    XTWX.resize(parameters, vector<long double>(parameters,0));
    XWy.resize(parameters, 0);
    beta.resize(parameters, 0);
    X_row.resize(parameters);
    
    X_row[0]=1;
    
    
    //cout<<"**********::**********"<<endl;
    #pragma omp declare reduction (add_m : vector<vector<long double> > : add_matrix(omp_out, omp_in)) initializer(omp_priv = vector<vector<long double> >(omp_orig.size(), vector<long double>(omp_orig.size())) )
            
    #pragma omp declare reduction(add_v : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                              initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
            
    #pragma omp parallel for reduction (add_m : XTWX) reduction (add_v : XWy) private (arg, X_row_j_wt) firstprivate (X_row)
    for(int i=0; i<n_datapoints; i++)
    {   
        for(int l=1; l<=n_sines; l++)
        {
            
            arg=dwaPI*sine_parameters[l][0]*date[i];
            X_row[l*2-1]=sin(arg);
            X_row[l*2]=cos(arg);
        }
        for(int j=0; j<parameters; j++)
        {
            X_row_j_wt=X_row[j] * weights[i];
            XWy[j] += X_row_j_wt * flux[i];
            for(int k=0; k<=j; k++)
            {
                XTWX[j][k] += X_row_j_wt * X_row[k];
            }
        }
    }
  
    

    for (int j=1; j<parameters; j++)    // Fill in above the diagonal from symmetry.
        for (int k=0; k<j; k++)
            XTWX[k][j]=XTWX[j][k];
        

    
 
   int l1, k1;
   long double fi;
   gauss(parameters, XTWX, XWy, beta);
   sine_parameters[0][0]=beta[0];
   for(int i=1; i<=n_sines; i++)
   {
       l1=2*i-1;
       k1=l1+1;
       sine_parameters[i][1]=sqrt(beta[l1]*beta[l1]+beta[k1]*beta[k1]);
       
       fi = atan(beta[k1]/beta[l1]);
       if (beta[l1] > 0 && beta[k1] < 0) fi += dwaPI;
       else if (beta[l1] < 0) fi += M_PI;
       sine_parameters[i][2]=fi;
       
       sine_parameters[i][2]/=(dwaPI);
   }
}



//solve system of linear equations (Gauss elimination). Matrix A is replaced by its LU decomposition
//Ax=b
void fit_sines::gauss(int unknows, vector<vector<long double> > &A, vector<long double> &b, vector<long double> &x)
{
    vector<long double>  s;
    vector<int> p;
    long double elem_max, z, suma;
    int tmp;
    
    p.resize(unknows);
    s.resize(unknows);
    x.resize(unknows);
    
    for(int i=0; i<unknows; i++)
    {
        p[i]=i;
        elem_max=0;
        for(int j=0; j<unknows; j++)
        {
            if(fabs(A[i][j])>elem_max)
            {
                elem_max=fabs(A[i][j]);
                s[i]=elem_max;
            }
        }
        //wywazanie wierszy
        for(int j=0; j<unknows; j++)
        {
            A[i][j]/=s[i];
        }
        b[i]/=s[i];
        //
    }
    
    
    
    for(int k=0; k<unknows-1; k++)
    {        
        for(int iter=k; iter<unknows; iter++)
        {
            if( fabs(A[p[iter]][k])/s[p[iter]] > fabs(A[p[k]][k])/s[p[k]] )
            {
                tmp=p[k];
                p[k]=p[iter];
                p[iter]=tmp;
                
            }
        }
        
        for(int i=k+1; i<unknows; i++)
        {
            z = A[p[i]][k] / A[p[k]][k];
            A[p[i]][k]=z;
            for(int j=k+1; j<unknows; j++)
            {
                A[p[i]][j]=A[p[i]][j]-z*A[p[k]][j];
            }
        }
    }
    
    for(int k=0; k<unknows-1; k++)
    {
        for(int i=k+1; i<unknows; i++)
        {
            b[p[i]]=b[p[i]]-A[p[i]][k]*b[p[k]];
        }
    }
    
    
    
    for(int i=unknows-1; i>=0; i--)
    {
        suma=0;
        for(int j=i+1; j<unknows; j++)
            suma+=A[p[i]][j]*x[j];
        
        x[i]=(b[p[i]]-suma)/A[p[i]][i];
    }
}


void fit_sines::gauss_inverse_matrix(int size, vector<vector<long double> > &A)
{
    vector<long double>  s;
    vector<int> p;
    long double elem_max, z, suma;
    int tmp;
    
    p.resize(size);
    s.resize(size);
    
    vector<vector<long double> > b;
    vector<vector<long double> > x;
    b.resize(size, vector<long double>(size, 0));
    x.resize(size, vector<long double>(size, 0));
    
    for(int i=0; i<size; i++)
        b[i][i]=1.0l;
    
    for(int i=0; i<size; i++)
    {
        p[i]=i;
        elem_max=0;
        for(int j=0; j<size; j++)
        {
            if(fabs(A[i][j])>elem_max)
            {
                elem_max=fabs(A[i][j]);
                s[i]=elem_max;
            }
        }
        //wywazanie wierszy
        for(int j=0; j<size; j++)
        {
            A[i][j]/=s[i];
        }
        b[i][i]/=s[i];
    }
    
    
    
    for(int k=0; k<size-1; k++)
    {        
        for(int iter=k; iter<size; iter++)
        {
            if( fabs(A[p[iter]][k])/s[p[iter]] > fabs(A[p[k]][k])/s[p[k]] )
            {
                tmp=p[k];
                p[k]=p[iter];
                p[iter]=tmp;
                
            }
        }
        
        for(int i=k+1; i<size; i++)
        {
            z = A[p[i]][k] / A[p[k]][k];
            A[p[i]][k]=z;
            for(int j=k+1; j<size; j++)
            {
                A[p[i]][j]=A[p[i]][j]-z*A[p[k]][j];
            }
        }
    }
    
    for(int l=0; l<size; l++)
    {
        for(int k=0; k<size-1; k++)
        {
            for(int i=k+1; i<size; i++)
            {
                b[p[i]][l]=b[p[i]][l]-A[p[i]][k]*b[p[k]][l];
            }
        }
    }
    
    
    
    for(int l=0; l<size; l++)
    {
        for(int i=size-1; i>=0; i--)
        {
            suma=0;
            for(int j=i+1; j<size; j++)
                suma+=A[p[i]][j]*x[j][l];
        
            x[i][l]=(b[p[i]][l]-suma)/A[p[i]][i];
        }
    }
    A=x;
}


void fit_sines::make_J_and_ymod(const int n_datapoints, const vector<long double> & t,\
                                      const int n_sines, const vector<vector<long double> > & sine_parameters,\
                                      const vector<vector<bool> > fit_control, vector<vector<long double> > &J, vector<long double> & ymod)
{
    long double arg;
    long double sinus, cosinus;
    int kk=0, i, j, k;
    
    if(fit_control[0][0])
    {
        kk=1;
        for(int i=0; i<n_datapoints; i++)
            J[i][0]=1;    //first parameter is constant offset
    }
            
    #pragma omp parallel for private (i, j, k, arg, sinus, cosinus)
    for(i=0; i<n_datapoints; i++)
    {
        ymod[i]=sine_parameters[0][0];
        for(j=1, k=kk; j<=n_sines; j++)
        {
            arg=dwaPI*(sine_parameters[j][0]*t[i]+sine_parameters[j][2]);
            sinus=sin(arg);
            cosinus=cos(arg);
            ymod[i] += sine_parameters[j][1]*sinus;
         
//            if(i<3)
  //              cout<<fixed<<setprecision(20)<<"aaaaa "<<sinus<<" "<<arg<<" "<<dwaPI<<" "<<sine_parameters[j][0]<<" "<<t[i]<<" "<<sine_parameters[j][1]<<" "<<sine_parameters[j][2]<<endl;

                        
            if(fit_control[j][0])
            {   
                J[i][k++] = dwaPI * sine_parameters[j][1] * t[i] * cosinus;
            }
            if(fit_control[j][1])
            {
                J[i][k++] = sinus;
            }
            if(fit_control[j][2])
                J[i][k++] = dwaPI * sine_parameters[j][1] * cosinus;
        }
    }
}



void fit_sines::Levenberg_Marquardt_fit(const vector<long double> & t, const vector<long double> & y, const vector<long double> & w,\
                                    const int n_sines, const int n_datapoints,\
                                    const vector<vector<bool> > fit_control, vector<vector<long double> > & sine_parameters,\
                                    const vector<vector<int> > & v_komb)
{
    int fitted_parameters=0, iter;
    vector<vector<long double> > J, JTWJ, A, AA, bb;
    vector<long double>  ymod, y_sub_ymod, JTWy, hlm, b;
    int max_iterations;
    long double wt;
    long double chi_p, chi_p_ad_h;
    long double ymod_p_ad_h;
    long double rho, rho_denominator, convergence_denominator=1e-3, max_par_actual;
    long double big_number=1e7;
    bool recalculate_eq_coeff=true;
    long double max_grad, max_par, chi_red, max_nu_change;
    
    
    for(std::vector<std::vector<bool> >::size_type i=0; i<fit_control.size(); i++)
    {
        for(std::vector<std::vector<bool> >::size_type j=0; j<fit_control[i].size(); j++)
        {
            if(fit_control[i][j])
                fitted_parameters++;
        }
    }
    J.resize(n_datapoints, vector<long double>(fitted_parameters));
    ymod.resize(n_datapoints);
    y_sub_ymod.resize(n_datapoints);
    
    JTWJ.resize(fitted_parameters, vector<long double>(fitted_parameters, 0));
    A.resize(fitted_parameters, vector<long double>(fitted_parameters));
    JTWy.resize(fitted_parameters, 0);
    hlm.resize(fitted_parameters);
    
    max_iterations=maxiter_multipler*fitted_parameters;
    
    //cout<<"max iter "<<max_iterations<<" "<<maxiter_multipler<<" "<<fitted_parameters<<" "<<eps1<<" "<<eps2<<" "<<eps5<<endl;
    
    for(iter=0; iter<max_iterations; iter++)
    {
        if(recalculate_eq_coeff)
        {
            for(int i=0; i<fitted_parameters; i++)
            {
                for(int j=0; j<=i; j++)  //fil symmetric matrix
                {
                    JTWJ[i][j]=0;
                }
                JTWy[i]=0;
            }

            
            if(v_komb.size() != 0)
            {
                for(std::vector<std::vector<long double> >::size_type i=1; i<sine_parameters.size(); i++)
                {
                    if(v_komb[i][0] == 1)
                        sine_parameters[i][0] = v_komb[i][1]*sine_parameters[v_komb[i][2]][0];
                    else if(v_komb[i][0] == 2)
                        sine_parameters[i][0] = v_komb[i][1]*sine_parameters[v_komb[i][2]][0]+v_komb[i][3]*sine_parameters[v_komb[i][4]][0];
                    else if(v_komb[i][0] == 3)
                        sine_parameters[i][0] = v_komb[i][1]*sine_parameters[v_komb[i][2]][0]+v_komb[i][3]*sine_parameters[v_komb[i][4]][0]\
                                               +v_komb[i][5]*sine_parameters[v_komb[i][6]][0];
                }
            }
            
            make_J_and_ymod(n_datapoints, t, n_sines, sine_parameters, fit_control, J, ymod);
            chi_p=0;


            
            /*
            for(int i=0; i<fitted_parameters; i++)
            {
                //y_sub_ymod[i]=y[i]-ymod[i];
                //chi_p+=y_sub_ymod[i]*y_sub_ymod[i]*w[i];
                        
                for(int j=0; j<=i; j++)
                {
                    for (int k=0; k<n_datapoints; k++)
                    {
                        if(i==0 && j==0)
                        {
                            y_sub_ymod[k]=y[k]-ymod[k];
                            chi_p+=y_sub_ymod[k]*y_sub_ymod[k]*w[k];
                        }
                        JTWJ[i][j] += J[k][i] * w[k] * J[k][j];

                        //cout<<i<<" "<<j<<" "<<k<<" "<<JTWJ[i][j]<<" "<<J[k][i]<<" "<<w[k]<<"  "<<J[k][j]<<endl;

                        if(j==0){
                            JTWy[i] += J[k][i] * w[k] * y_sub_ymod[k];
                            //if(i==1) cout<<fixed<<JTWy[i]<<endl;
                        }
                    }
                }
            }
            
           */
            
            
            //this version is faster
            long double dy;

            /*
            for (int i=0, j, l, k, m; i<n_datapoints; i++)    //Summation loop over all data.
            {
                //(*funcs)(x[i],a,&ymod,dyda,ma);
                 //sig2i=w[i];
                 dy=y[i]-ymod[i];
                 
                 for (j=-1, l=0; l<fitted_parameters; l++) 
                 {
                     { 
                         wt=J[i][l]*w[i];
                         for (j++, k=-1, m=0; m<=l; m++)
                         {
                            JTWJ[j][++k] += wt*J[i][m];
                         }
    	                 JTWy[j] += dy*wt;
        
                     } 
                }
                chi_p += dy*dy*w[i];   //And find χ2 .
            }
            */
            
            //this version is faster
            //export OMP_NUM_THREADS=8;
            
            
            
            int i, j, l, k, m;
            #pragma omp declare reduction (add_m : vector<vector<long double> > : add_matrix(omp_out, omp_in)) initializer(omp_priv = vector<vector<long double> >(omp_orig.size(), vector<long double>(omp_orig.size())) )
            
            #pragma omp declare reduction(add_v : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                              initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
            
            #pragma omp parallel for reduction (+ : chi_p) reduction (add_m : JTWJ) reduction (add_v : JTWy) private (i, j, l, k, m, dy, wt)
            for (i=0; i<n_datapoints; i++)    //Summation loop over all data.
            {
                 dy=y[i]-ymod[i];
                 
                 for (j=-1, l=0; l<fitted_parameters; l++) 
                 {
                     { 
                         wt=J[i][l]*w[i];
                         for (j++, k=-1, m=0; m<=l; m++)
                         {
                            JTWJ[j][++k] += wt*J[i][m];
                         }
    	                 JTWy[j] += dy*wt;
        
                     } 
                }
                chi_p += dy*dy*w[i];   //And find χ2 .
            }
            
            
        
            #pragma omp parallel for 
            for (int j=1;j<fitted_parameters; j++)   //Fill in the symmetric side.
                for (int k=0; k<j; k++) JTWJ[k][j]=JTWJ[j][k];
                
                
        }


        A=JTWJ;
        for(int i=0; i<fitted_parameters; i++)
            A[i][i] += lambda0 * A[i][i];
        
        b=JTWy;
        //gauss change A and b
        gauss(fitted_parameters, A, b, hlm);
        


        long double arg_p_ad_h;
        vector<vector<long double> > sine_parameters_add_h=sine_parameters;
        chi_p_ad_h=0;
        
        

        for(int i=0, hindex=0; i<=n_sines; i++)
        {
            if(i==0 && fit_control[i][0])
                sine_parameters_add_h[i][0] += hlm[hindex++];

            if(i>0 && fit_control[i][0])
                sine_parameters_add_h[i][0] += hlm[hindex++];
            if(i>0 && fit_control[i][1])
                sine_parameters_add_h[i][1] += hlm[hindex++];
            if(i>0 && fit_control[i][2])
                sine_parameters_add_h[i][2] += hlm[hindex++];
        }
        
        if(v_komb.size() != 0)
        {
            for(std::vector<std::vector<long double> >::size_type i=1; i<sine_parameters_add_h.size(); i++)
            {
                if(v_komb[i][0] == 1)
                    sine_parameters_add_h[i][0] = v_komb[i][1]*sine_parameters_add_h[v_komb[i][2]][0];
                else if(v_komb[i][0] == 2)
                    sine_parameters_add_h[i][0] = v_komb[i][1]*sine_parameters_add_h[v_komb[i][2]][0]\
                                                 +v_komb[i][3]*sine_parameters_add_h[v_komb[i][4]][0];
                else if(v_komb[i][0] == 3)
                    sine_parameters_add_h[i][0] = v_komb[i][1]*sine_parameters_add_h[v_komb[i][2]][0]\
                                                 +v_komb[i][3]*sine_parameters_add_h[v_komb[i][4]][0]\
                                                 +v_komb[i][5]*sine_parameters_add_h[v_komb[i][6]][0];
            }
        }

         ////////////////////////////////////
        #pragma omp parallel for reduction(+ : chi_p_ad_h) private (ymod_p_ad_h)
        for(int i=0; i<n_datapoints; i++)
        {
            ymod_p_ad_h=sine_parameters_add_h[0][0];
                        
            for(int j=1; j<=n_sines; j++)
            {
                arg_p_ad_h = dwaPI * (sine_parameters_add_h[j][0]*t[i] + sine_parameters_add_h[j][2]);
                ymod_p_ad_h += sine_parameters_add_h[j][1] * sin(arg_p_ad_h);
                
            }
            chi_p_ad_h += (y[i]-ymod_p_ad_h)*(y[i]-ymod_p_ad_h)*w[i];
        }
        
        /*
        for(int i=0; i<n_datapoints; i++)
        {
            ymod_p_ad_h=sine_parameters_add_h[0][0];
                        
            for(int j=1; j<=n_sines; j++)
            {
                arg_p_ad_h = dwaPI * (sine_parameters_add_h[j][0]*t[i] + sine_parameters_add_h[j][2]);
                ymod_p_ad_h += sine_parameters_add_h[j][1] * sin(arg_p_ad_h);
                
            }
            chi_p_ad_h += (y[i]-ymod_p_ad_h)*(y[i]-ymod_p_ad_h)*w[i];
        }
        */
        
        rho_denominator=0;
        for(std::vector<std::vector<bool> >::size_type i=0; i<fit_control.size(); i++)
        {
            rho_denominator+=hlm[i] * (lambda0*JTWJ[i][i]*hlm[i] + JTWy[i]);
        }
        rho=(chi_p-chi_p_ad_h)/fabs(rho_denominator);
        
        
        if(rho>eps4)
        {
            sine_parameters=sine_parameters_add_h;
            //lambda0 = max(lambda0*0.1, small_number);
            lambda0 *= 0.1;
            recalculate_eq_coeff=true;
        }
        else
        {
            lambda0 = min(lambda0*10, big_number);
            recalculate_eq_coeff=false;
        }
        
        //convergence test
        max_grad=max_nu_change=max_par=-1e6;
        
        for(int i=0; i<fitted_parameters; i++)
        {
            if(fabs(JTWy[i])>max_grad)
                max_grad = fabs(JTWy[i]);
        }
        
        for(int i=0, hindex=0; i<=n_sines; i++)
        {
            if(i==0 && fit_control[i][0])
            {
                fabs(sine_parameters[i][0])>1e-3 ? max_par_actual = fabs(hlm[hindex]/sine_parameters[i][0])\
                                                 : max_par_actual = fabs(hlm[hindex]/convergence_denominator);
                if(max_par_actual>max_par)
                    max_par=max_par_actual;
                hindex++;
            }
                
            if(i>0 && fit_control[i][0])
            {
                fabs(sine_parameters[i][0])>1e-3 ? max_par_actual = fabs(hlm[hindex]/sine_parameters[i][0])\
                                                 : max_par_actual = fabs(hlm[hindex]/convergence_denominator);
                if(max_par_actual>max_par)
                    max_par=max_par_actual;
                
                if(fabs(hlm[hindex])>max_nu_change)
                    max_nu_change=fabs(hlm[hindex]);
                hindex++;
            }
            
            if(i>0 && fit_control[i][1])
            {
                fabs(sine_parameters[i][1])>1e-3 ? max_par_actual = fabs(hlm[hindex]/sine_parameters[i][1])\
                                                 : max_par_actual = fabs(hlm[hindex]/convergence_denominator);
                if(max_par_actual>max_par)
                    max_par=max_par_actual;
                hindex++;
            }
            
            if(i>0 && fit_control[i][2])
            {
                fabs(sine_parameters[i][2])>1e-3 ? max_par_actual = fabs(hlm[hindex]/sine_parameters[i][2])\
                                                 : max_par_actual = fabs(hlm[hindex]/convergence_denominator);
                if(max_par_actual>max_par)
                    max_par=max_par_actual;
                hindex++;
            }
        }
        chi_red=chi_p/(n_datapoints-fitted_parameters+1);
        
//        cout<<"chi "<<chi_p<<" "<<chi_red<<"  "<<chi_p-chi_p_ad_h<<endl<<endl;
        
        if((max_par<eps2 || chi_red <eps3) && max_nu_change<eps5)
        {
            cout<<"converged after "<<iter<<" iterations"<<endl;
            break;
        }
    }
    cout<<"STOP: "<<"relatvie max parameter change: "<<max_par<<" -- relatvie max frequency change: "<<max_nu_change<<"  --- chi_red: "<<chi_red<<endl<<endl;
    if(iter == max_iterations)
        cout<<"ups ---- probably L-M not converged"<<endl;
    
    
    //calculate errors  
    err_Levenberg_Marquardt(t, y, w, n_sines, n_datapoints, fit_control, sine_parameters, fitted_parameters, JTWJ);
    
    
    
}

void fit_sines::err_Levenberg_Marquardt(const vector<long double> & t, const vector<long double> & y, const vector<long double> & w,\
                                        const int n_sines, const int n_datapoints,\
                                        const vector<vector<bool> > fit_control, vector<vector<long double> > & sine_parameters,\
                                        const int fitted_parameters, vector<vector<long double> > & JTWJ)

{
    gauss_inverse_matrix(fitted_parameters, JTWJ);
    
    long double sum, var, oc, sdev2, rat, swag;
  
    var = swag = 0;
    
    for(int i=0; i<n_datapoints; i++)
    {
        sum = sine_parameters[0][0];
        //sum=0;
        swag += w[i];
        for(int j=1; j<=n_sines; j++)
            sum +=  sine_parameters[j][1]*sin(dwaPI*(sine_parameters[j][0]*t[i] + sine_parameters[j][2]));
        
        oc = y[i] - sum;
        var = var + oc*oc;
    }
    //cout<<"var  "<<var<<endl;
    sdev2 = var/(n_datapoints - 1);
    rat = n_datapoints/(swag*sdev2);
    
    //cout<<" RAT  "<<rat<<endl;
    
    for(int i=0, hindex=0; i<=n_sines; i++)
    {
        if(i==0 && fit_control[i][0])
        {
            sine_parameters[i][1] = sqrt(JTWJ[hindex][hindex]/rat);
            hindex++;
        }
        if(i>0 && fit_control[i][0])
        {
            sine_parameters[i][3] = sqrt(JTWJ[hindex][hindex]/rat);
            hindex++;
        }
        if(i>0 && fit_control[i][1])
        {
            sine_parameters[i][4] = sqrt(JTWJ[hindex][hindex]/rat);
            hindex++;
        }
        if(i>0 && fit_control[i][2])
        {
            sine_parameters[i][5] = sqrt(JTWJ[hindex][hindex]/rat);
            hindex++;
        }
    }
}






fit_sines::fit_sines()
{
    //parameters for L-M
    lambda0=0.001; //initial value of L-M parameter 
    eps1=1e-3; //convergence tolerance for gradient
    eps2=1e-3; //convergence tolerance for parameters
    eps3=1e-1; //convergence tolerance for red. Chi-square
    eps4=1e-1; //determines acceptance of a L-M step
    eps5=1e-8; // determines max absolute change of frequency
    maxiter_multipler=50;  //maxiter=maxiter_multipler *  Npar
}

#include "interface.h"

using namespace std;


bool string_to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}


interface::interface(int _narguments, char* argv[]): narguments(_narguments)
{
    cout<< "\E[22;32m"<<"#############################################"<< "\E[0m"<<endl;
    cout<< "\E[1;32m"<<"LS and fit for Linux, version of May 13, 2023"<< "\E[0m"<<endl;
    cout<< "\E[1;32m"<<"by Wojtek Szewczuk"<< "\E[0m"<<endl;
    cout<< "\E[22;32m"<<"#############################################"<< "\E[0m"<<endl<<endl;

    for(int i=0; i<narguments; i++)
        call_parameters.push_back(argv[i]);

    if(narguments != 3 && narguments !=5)
    {
        cout<<"Minimal usage: LS -f <DataFile>"<<endl;
        cout<<"Other parameters are read from the default configuration file <Config>"<<endl;
        cout<<"Optional usage: LS -c <myConfigFile>"<<endl;
        cout<<"or            : LS -f <DataFile> -c <myConfigFile>"<<endl;
        exit(1);
    }


    for(int i=0; i<narguments; i++)
    {
        if(call_parameters[i]=="-f")
        {
            cp.insert(std::pair{"data_file", call_parameters[i+1]});
        }

        if(call_parameters[i]=="-c")
        {
            cp.insert(std::pair{"config_file", call_parameters[i+1]});
        }
    }
    cp.insert(std::pair{"config_file", "Config"}); //will add the key only if it has not been added before from command line
    make_default_cp();



    read_config();



    cout<<"#####################################";
    cout<<endl<<endl<<endl;
}



void interface::read_config()
{
    ifstream in(cp["config_file"]);
    string line;
    string mapkey, mapvalue;
    map<string,string>::iterator mapiter;

    if(in.good())
    {
        while(true)
        {
            getline(in, line);
            if(in.eof())
                break;
            config_content.push_back(line);
        }
    }
    else
    {
        cout<<"Lack of the configfile"<<endl;
        cout<<"Default configuration (perhaps inappropriate) is used"<<endl;
    }
    in.close();


    std::string::size_type n;
    for(std::vector<std::__cxx11::basic_string<char> >::size_type i=0; i<config_content.size(); i++)
    {
        n=config_content[i].find("//");
        if( n!= string::npos)
            config_content[i]=config_content[i].substr(0, n);//remove comments

        config_content[i].erase(std::remove_if(config_content[i].begin(), config_content[i].end(), ::isspace), config_content[i].end());//remove white spaces

        if(config_content[i].size() == 0 || std::all_of(config_content[i].cbegin(),config_content[i].cend(),[](char c) { return std::isspace(c); }) ) //remove empty lines
        {
            config_content.erase(config_content.begin()+i);
            i--;
        }
    }




    for(std::vector<std::__cxx11::basic_string<char> >::size_type i=0; i<config_content.size(); i++)
    {
        n = config_content[i].find("=");
        mapkey = config_content[i].substr(0, n);
        mapvalue = config_content[i].substr(n+1, config_content[i].back());
        mapiter = cp.find(mapkey);
        if( mapiter != cp.end() )
            cp.erase( mapiter );  //overwrite previous entries
        cp.insert(std::pair{mapkey,mapvalue});
    }


    cout<<"++++++++++++++++++++++++++++"<<endl;
    for (const auto& p : cp) {
        cout << "(" << p.first << ", " << p.second << ") "<<endl;
    }
    cout<<"++++++++++++++++++++++++++++"<<endl;


}



void interface::make_default_cp()
{
    cp.insert(std::pair{"freq_file", "freq.dat"});
    cp.insert(std::pair{"save_spectra", "true"});
    cp.insert(std::pair{"noise_in_window", "true"});
    cp.insert(std::pair{"window_size", "1"});
    cp.insert(std::pair{"StoN_criterion", "5"});
    cp.insert(std::pair{"fit_all_freq", "true"});
    cp.insert(std::pair{"fit_n_last_freq", "1"});

    cp.insert(std::pair{"spec_mode_default", "true"});
    cp.insert(std::pair{"spec_res", "0.0001"});
    cp.insert(std::pair{"spec_high_freq", "25"});

    cp.insert(std::pair{"remove_too_close", "false"});
    cp.insert(std::pair{"remove_too_close_hard", "false"});
    cp.insert(std::pair{"remove_too_close_fac", "2.5"});
    cp.insert(std::pair{"threads", "1"});

    cp.insert(std::pair{"Czerny_err_corr", "true"});

    cp.insert(std::pair{"set_exact_com_and_har", "false"});

    cp.insert(std::pair{"har_range", "10"});
    cp.insert(std::pair{"com_2par_range", "10"});
    cp.insert(std::pair{"com_3par_range", "3"});

    cp.insert(std::pair{"com_par_range_v2", "10"});
    cp.insert(std::pair{"har_range_v2", "10"});
    //cp.insert(std::pair{"T0", "1870"}); //no default value -> get floor from the first observation

    cp.insert(std::pair{"allow_3_parents_v2", "true"});
    cp.insert(std::pair{"try_n_highest_peaks", "1"});
    cp.insert(std::pair{"set_try_n_highest_peaks_eq_1_after_found_n_freq", "100000000000000000000"}); // set try_n_highest_peaks equal to 1 after preweithen on n frequency
}

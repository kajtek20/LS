#include "interface.h"

using namespace std;

interface::interface(int _narguments, char* argv[]): narguments(_narguments)
{
    cout<<"#####################################"<<endl;
    cout<<"LS for Linux, version of Mar 16, 2023"<<endl;
    cout<<"by Wojtek Szewczuk"<<endl;
    cout<<"#####################################"<<endl<<endl;
    
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
            //datafile=call_parameters[i+1];
            cp.insert(std::pair{"data_file", call_parameters[i+1]});
        }
        
        if(call_parameters[i]=="-c")
        {
            //configfile=call_parameters[i+1];
            cp.insert(std::pair{"config_file", call_parameters[i+1]});
        }
    }
    cp.insert(std::pair{"config_file", "Config"}); //will add the key only if it has not been added before from command line
    
    

    read_config();


    cp.insert(std::pair{"mass", "%%%%"});
    map<string,string>::iterator iter = cp.find("mass") ;
    if( iter != cp.end() )
    cp.erase( iter );
    cp.insert(std::pair{"mass", "#####"});



    cout<<"------------------ "<<cp["mass"]<<" "<<cp.count("mass0")<<endl;
    cout<<"------------------ "<<cp["data_file"]<<" "<<cp.count("mass0")<<endl;
    cout<<"------------------ "<<cp["config_file"]<<" "<<cp.count("mass0")<<endl;
    
    cout<<"#####################################";
    cout<<endl<<endl<<endl;
}



void interface::read_config()
{
    ifstream in(cp["config_file"]);
    string line;
    //int config_size;
    
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
    
    for(std::vector<std::__cxx11::basic_string<char> >::size_type i=0; i<config_content.size(); i++)
        cout<<config_content[i].size()<<" "<<config_content[i]<<endl;
    
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
    
    
    
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    for(std::vector<std::__cxx11::basic_string<char> >::size_type i=0; i<config_content.size(); i++)
        cout<<config_content[i].size()<<" "<<config_content[i]<<endl;

    //config_parameters_double dq;
    //config_parameters_double dw;
    //par_v.push_back(dq);
}




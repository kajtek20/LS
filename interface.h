#ifndef interface_h
#define interface_h

#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>

using namespace std;








class interface
{
public:
    interface(int _narguments, char* argv[]);
    
    int narguments;
    vector<string> call_parameters;
    

    std::map<std::string, std::string> cp;
    
private:
    void read_config();
    void make_default_cp();
    vector<string> config_content;
};




bool string_to_bool(std::string str);


#endif

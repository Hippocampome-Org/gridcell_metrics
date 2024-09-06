#include <cstdio>
#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <regex>

using namespace std;

int main() {
	ifstream input_file("config.txt");
	//ofstream output_file("run_gridcell_metrics.sh");
	FILE *output_file = fopen("run_gridcell_metrics.sh", "w");
    string str, str2;
    vector<vector<string>> config;
    vector<string> config_line;
    regex newlines_re("\n+");
    regex str_expr ("(.*)\n?");
    smatch sm;    

    while (getline(input_file, str))
    {
		stringstream ss(str);
		string substr;
		while(getline(ss, substr, ','))
		{
		    config_line.push_back(substr);
		}
		//printf("%s \n",config_line[1].c_str());
		config.push_back(config_line);
    	config_line.clear();
    }

    for(int i = 1; i < size(config); i++) {
    	str.clear();
    	str = config[i][1];
    	//printf("%s \n",str.c_str());
    	//printf("%s ",str.c_str());
  		//auto str2 = regex_replace(str, newlines_re, "");
  		//char* result2 = result.c_str(); 
  //		char* arr = new char[result.length()];
 //	 		strcpy(arr, result.c_str());
    	//fprintf(output_file,"%s ",arr);
   		//regex_match (str,sm,str_expr);
   		//for (unsigned i=0; i<sm.size(); ++i) {
   		//	string str2 = sm[1];
	    //	fprintf(output_file,"%s ",str2.c_str());
	   	//}
    	//fprintf(output_file,"%s ",str.c_str());
    	//str2.clear();
    	//char* str2 = str.c_str();
    	char str3[str.size()];
		strncpy(str3, str.c_str(), sizeof(str3));
    	for(int j = 0; j<str.length(); j++) {
    	//	cout << str.at(j);
    	}
				string token;    	
    	for (char * p = str3; p != str3 + sizeof(str3) / sizeof(str3[0]); ++p)
		{	
			if (*p != ' ' && *p != (char)10 && *p != (char)13) {
				token.push_back(*p);
		 		//cout << *p << endl;
			}
		}
		 		fprintf(output_file,"%s ",token.c_str());		
		 		token.clear();
    }
    fclose(output_file);
    input_file.close();
    //cout << "size: " << size(config);
    
	return 0;
}
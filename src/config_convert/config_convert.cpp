#include <cstdio>
#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

int main() {
	ifstream input_file("config.txt");
	FILE *output_file = fopen("../run_gridcell_metrics_with_config.sh", "w");
    string str, token;
    vector<vector<string>> config;
    vector<string> config_line;

    // collect tokens
    while (getline(input_file, str)) {
		stringstream ss(str);
		string substr;
		while(getline(ss, substr, ',')) {
		    config_line.push_back(substr);
		}
		config.push_back(config_line);
    	config_line.clear();
    }

    // write tokens to file
    token.append("./run_gridcell_metrics ");
    for(int i = 1; i < size(config); i++) {
    	str.clear();
    	str = config[i][1];
    	char cl_char[str.size()];
		strncpy(cl_char, str.c_str(), sizeof(cl_char));
    	for (char * p = cl_char; p != cl_char + sizeof(cl_char) / sizeof(cl_char[0]); ++p) {	
			if (*p != ' ' && *p != (char)10 && *p != (char)13) {
				token.push_back(*p);
			}
		}
 		fprintf(output_file,"%s ",token.c_str());		
 		token.clear();
    }
    fclose(output_file);
    input_file.close();
    
	return 0;
}
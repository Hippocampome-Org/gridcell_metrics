#include <cstdio>
#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
    if(argc > 3) {
        printf("Too many arguments supplied. Include arguments: output_filename compiled_filename\n");
        return 1;
    }
    else if (argc <= 2) {
        printf("Too few arguments supplied. Include arguments: output_filename compiled_filename\n");
        return 1;
    }

    string output_filename = argv[1];
    string compiled_filename = argv[2];

	ifstream input_file("rate_map_config.txt");
	FILE *output_file = fopen(output_filename.c_str(), "w");
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
    token.append(compiled_filename);
    token.append(" ");
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
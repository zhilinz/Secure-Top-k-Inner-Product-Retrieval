#ifndef BASE_H_INCLUDED
#define BASE_H_INCLUDED

//#include <iostream>
//#include <stdint.h>
#include <string>
//#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
//#include <queue>
//#include <cmath>
//#include <dirent.h>
//#include <map>
//#include <stdexcept>
//#include <algorithm>
//#include <cfloat>
//#include <cstdlib>
#include <sys/stat.h>

using namespace std;

typedef unsigned long long ull;

string getCurrentTime(const char *pattern = "%Y%m%d-%H%M%S");

string to_string(int n);

void split(const string &s, vector<string> &elems, char delim=',') ;


// the size of a file in K bytes

double evaluateFileSize(string filepath);


#endif //BASE_H_INCLUDED

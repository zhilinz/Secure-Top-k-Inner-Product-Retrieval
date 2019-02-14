#include "base.h"

string getCurrentTime(const char *pattern) {
    time_t cur_time = time(0);
    char buf[100];
    strftime(buf, 100, pattern, localtime(&cur_time)); //format date and time.
    return buf;
}

string to_string(int n) {
    ostringstream stm;
    stm << n;
    return stm.str();
}

void split(const string &s, vector<string> &elems, char delim) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}



// the size of a file in K bytes

double evaluateFileSize(string filepath)
{

    double fileSize;

    struct stat st;
    stat(filepath.c_str(), &st);

    fileSize=st.st_size;

    return fileSize/1000;
}



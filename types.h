#ifndef TYPES_H
#define TYPES_H

#include <cstring>

//typedef long double myFloat;
typedef float myFloat;

class My_string{
public:
    My_string(char str[256] ) {
        strcpy(this->str,str);
    }
    char str[256];
};




#endif // TYPES_H

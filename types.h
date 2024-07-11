#ifndef TYPES_H
#define TYPES_H

typedef long double myFloat;

class My_string{
public:
    My_string(char str[256] ) {
        strcpy(this->str,str);
    }
    char str[256];
};




#endif // TYPES_H

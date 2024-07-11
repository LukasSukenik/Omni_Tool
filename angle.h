#ifndef ANGLE_H
#define ANGLE_H

#include <vector>

using namespace std;

class Angle {
public:
    Angle() : N(-1), type(-1), at1(-1), at2(-1), at3(-1) {}
    Angle(int N, int type, int at1, int at2, int at3) : N(N), type(type), at1(at1), at2(at2), at3(at3) {}
    int N, type, at1, at2, at3;
};

class Angles : public vector< Angle >
{
public:
    void offset(int offs)
    {
        for(Angle& item : (*this))
        {
            item.N += offs;
            item.at1 += offs;
            item.at2 += offs;
            item.at3 += offs;
        }
    }
};

#endif // ANGLE_H

#ifndef TYPES_H
#define TYPES_H

#include <algorithm>>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <cstring>

using namespace std;

//typedef long double myFloat;
typedef float myFloat;




class My_string{
public:
    My_string(char str[256] ) {
        strcpy(this->str,str);
    }
    char str[256];
};




ostream& operator<<(ostream& os, const vector<int>& values)
{
    for(const auto& val : values)
    {
        os << val << ' ';
    }
    return os;
}




class TMD
{
public:
    TMD(){}

    int size;
    int proximal_n;
    int distal_n;

    void clear()
    {
        size = 0;
        proximal_n = 0;
        distal_n = 0;
    }
};




template <typename T>
unordered_set<T> operator+(unordered_set<T> a, unordered_set<T>& b)
{
    a.merge(b);
    return a;
}

template <typename T>
class Param_Dictionary : public unordered_map<string, T>
{
public:
    unordered_set<string> valid_keys;

    Param_Dictionary(unordered_set<string> valid_key_list) : valid_keys(valid_key_list) {}

    bool is_key_valid(string key)
    {
        return valid_keys.contains(key);
    }

    void validate_keyword(string keyword, string default_value)
    {
        if( !this->contains(keyword) )
        {
            cerr << "Missing keyword; " << keyword << ": " << default_value << endl;
            exit(-1);
        }
    }
};




/**
 * @brief Populate
 * Defines the population of particles in the simulation box
 */
class Population
{
public:
    Population(){}

    bool random=false;
    int count=0;

    bool empty()
    {
        if(count == 0)
            return true;
        return false;
    }

    void clear()
    {
        random=false;
        count=0;
    }

    friend std::ostream& operator<<(std::ostream& os, const Population& pop)
    {
        if(pop.random)
        {
            os << pop.count << " random";
        }
        else
        {
            os << "not defined";
        }
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Population& pop)
    {
        string str;
        is >> pop.count >> str;

        transform(str.begin(), str.end(), str.begin(), ::tolower);

        if (str == "random")
        {
            pop.random = true;
        }

        return is;
    }
};




enum class IO_Type { none, xyz, pdb, gro, lammps_full };

std::ostream& operator<<(std::ostream& os, const IO_Type out)
{
    switch (out) {
    case IO_Type::gro:
        os << "gro"; return os;
    case IO_Type::xyz:
        os << "xyz"; return os;
    case IO_Type::pdb:
        os << "pdb"; return os;
    case IO_Type::lammps_full:
        os << "lammps_full"; return os;
    default:
        os << "none"; return os;
    }
    return os;
}

class IO{
public:
    IO() {}

    IO_Type type = IO_Type::none;

    void clear()
    {
        type = IO_Type::none;
    }

    friend std::istream& operator>>(std::istream& is, IO& out)
    {
        string str;
        is >> str;

        transform(str.begin(), str.end(), str.begin(), ::tolower);

        out.type = IO_Type::none;
        if (str == "xyz")
        {
            out.type = IO_Type::xyz;
        }
        if (str == "gro")
        {
            out.type = IO_Type::gro;
        }
        if (str == "pdb")
        {
            out.type = IO_Type::pdb;
        }
        if (str == "lammps_full")
        {
            out.type = IO_Type::lammps_full;
        }

        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const IO& io)
    {
        os << io.type;
        return os;
    }
};




#endif // TYPES_H

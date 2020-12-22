#ifndef HH4bStruct_H
#define HH4bStruct_H

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <TROOT.h>
#include <boost/bind.hpp>

struct HH4bStruct
{
    float mass=0;
    TLorentzVector momentum;
    float Ht=0;
    std::vector<float> pt;
    int index;
};

#endif

#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <set>
#include <iostream>
#include <chrono>
#include <map>
#include <vector>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <numeric>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "pcg_random.hpp"
#include "pcg_uint128.hpp"
#include "pcg_extras.hpp"

pcg_extras::seed_seq_from<std::random_device> seed;
pcg32 rnEngine(seed);
// std::default_random_engine rnEngine;
std::uniform_real_distribution<double> distDouble0To1(0,nextafter(1,2));
std::uniform_int_distribution<int_fast64_t> distInt0To2(0,2);
std::uniform_int_distribution<int_fast64_t> distInt0To1(0,1);

typedef std::vector<int_fast64_t> vint;
typedef std::vector<double> vdbl;
typedef std::vector<vint> vvint;
typedef std::vector<vvint> vvvint;
typedef std::vector<std::pair<int_fast64_t,int_fast64_t>> vpint;
typedef std::vector<std::string> vstr;
typedef std::set <int_fast64_t> sint;

template <typename IntType>
std::vector<IntType> range(IntType start, IntType stop, IntType step)
{
    if (step == IntType(0))
    {
        throw std::invalid_argument("step for range must be non-zero");
    }

    std::vector<IntType> result;
    IntType i = start;
    while ((step > 0) ? (i < stop) : (i > stop))
    {
        result.push_back(i);
        i += step;
    }
    return result;
}

template <typename IntType>
std::vector<IntType> range(IntType start, IntType stop)
{
    return range(start, stop, IntType(1));
}

template <typename IntType>
std::vector<IntType> range(IntType stop)
{
    return range(IntType(0), stop, IntType(1));
}

template <typename Type>
bool Type_in_vType(std::vector<Type> vec,Type intg)
{
    return (std::find(vec.begin(),vec.end(),intg) != vec.end());
}

template <typename Type>
bool Type_notin_vType(std::vector<Type> vec,Type intg)
{
    return (std::find(vec.begin(),vec.end(),intg) == vec.end());
}

int_fast64_t fmod(const int_fast64_t input, const int_fast64_t ceil) 
{
    // apply the modulo operator only when needed
    // (i.e. when the input is greater than the ceiling)
    return input >= ceil ? input % ceil : input;
    // NB: the assumption here is that the numbers are positive
}

int_fast64_t mod(int_fast64_t a,int_fast64_t b)
{
    //we assume a<1000000*b
    if (a<0)
        return fmod(abs(a)*b+a,b);
    else 
        return fmod(a,b);
}

bool one_in_vector(vint vec)
{
    return (std::find(vec.begin(),vec.end(),1)!=vec.end());
}

int_fast64_t countone_in_vector(vint vec)
{
    int_fast64_t count=0;
    for(int_fast64_t i=0; i<vec.size();++i)
    {
        if (vec[i]==1)
        {
            count+=1; 
        }
    }
    return count;
}

int_fast64_t mod2dot_fun(vint logicalcoords,vint error)
{
    int_fast64_t count=0;
    for(int_fast64_t i=0;i<logicalcoords.size();++i)
    {
        if(error[logicalcoords[i]]==1) 
            count=fmod(count+1,2);
    }
    return count;
}

vint sum_arr(vint a,vint b)
{
    vint c;
    for(int_fast64_t i=0; i<a.size();++i)
    {
        c.push_back(a[i]+b[i]);
    }
    return c;
}


vint sum_arr_mod(int_fast64_t L,vint a,vint b)
{
    vint c;
    for(int_fast64_t i=0; i<a.size();++i)
    {
        c.push_back(fmod((a[i]+b[i]),2*L));
    }
    return c;
}

vint subtract_arr(vint a,vint b)
{
    vint c;
    for(int_fast64_t i=0; i<a.size();++i)
    {
        c.push_back(a[i]-b[i]);
    }
    return c;
}
vint mul_arr(vint a,vint b)
{
    vint c;
    for(int_fast64_t i=0; i<a.size();++i)
    {
        c.push_back(a[i]*b[i]);
    }
    return c;
}

vint vectimesc(int_fast64_t k,vint v)
{
	vint vk;
    for(int_fast64_t x:v)
    {
        vk.push_back(x*k);
    }
    return vk;
}

void updateset(sint& set,int position)
{
    auto it = set.find(position);
    if (it == set.end())
    {
        set.insert(position);
    }
    else
    {
        set.erase(it);
    }
}

void tests()
{
    vint a{1,0,0,1};
    vint b{0,1,0,1};
    vint c{0,0,0,0};
    
    assert(mod2dot_fun(a,b)==1);
    assert(one_in_vector(a)==true);
    assert(one_in_vector(b)==true);
    assert(one_in_vector(c)==false);
}
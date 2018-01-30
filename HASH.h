#ifndef HASH_H
#define HASH_H
#include <utility>
using namespace std;
struct HASH
{
    size_t operator()(const pair<int, int> &x) const
    {
        return (size_t)x.first * 37U + (size_t)x.second;
    }
};
#endif
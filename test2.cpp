#include<fstream>
#include<set>
#include<unordered_set>
#include<iostream>
using namespace std;

struct cmp {
    bool operator() (const int &l, const int &r) const {
        if(l%2 == r%2) return l>r;
        else return (l%2) > (r%2);
    }
};

void f(){
    set<int,cmp> s;
    for(int i=0;i<10;i++) s.insert(i);
    for(auto it:s) cout << it << '\n';
}

int main(){
    f();
    return 0;   
}

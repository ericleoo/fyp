#include<fstream>
#include<set>
#include<unordered_set>
#include<iostream>
#include "ProgressBar.h"
using namespace std;

void f(){
    unordered_set<int> s;
    for(int i=0;i<10000000;i++){
        s.insert(i);
    }
    s.clear();
    s.rehash(0);
    cout << s.bucket_count() << '\n';
}

int main(){
    f();
    return 0;
    
}

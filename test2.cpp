#include<fstream>
#include<set>
#include<unordered_set>
#include<iostream>
using namespace std;

void f(){
    set<int> s;
    
    cout << "One\n";
    for(int i=0;i<50000000;i++){
        s.insert(i);
    }
    s.clear();
    // s.rehash(0);

    cout << "Two\n";
    /*
    un./a.oordered_set<int> w;
    for(int i=0;i<100000000;i++){
        w.insert(i);
    }
    w.clear();
    */

    cout << "Done\n";
    getchar();
}

int main(){
    f();
    return 0;
    
}

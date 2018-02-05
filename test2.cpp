#include<fstream>
#include<set>
#include "ProgressBar.h"
using namespace std;

int main(){
    ProgressBar bar(5);
    
    bar.SetNIter(1000000);
    bar.Reset();
        
    set<int> s;
    for(int i=0;i<1000000;i++){
        bar.Update();
        s.insert(i);
    }
    
    fstream fout("hihi.txt",fstream::app);
    fout << "fgrg\n";
    fout << "fbsbs\n";
    fout.close();
        
    return 0;
}

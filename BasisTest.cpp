#include "Basis.h"
#include <cstdio>
#include <cassert>

using namespace std;

//a series of tests
int main(){
    WedgeBasis basis(2,5,7);
    auto iter = basis.getIter();
    int count = 0;
    do{
        count++;
    }while(iter.next());
    printf("Check %d==%d\n",count,basis.size());
    assert(count==basis.size());
    return 0;
}

#include "Basis.h"
#include <cstdio>
#include <cassert>

using namespace std;

//a series of tests
int main(){
    WedgeBasis basis(createReducedWedgeBasis(2,7,20));
    auto iter = basis.getIter();
    long long count = 0;
    do{
        printf("%lld\n",count);
        assert(count==basis.rank(iter.getCurr()));
        assert(iter.getCurr()==basis.unrank(count));
        count++;
    }while(iter.next());
    printf("Check %lld==%lld\n",count,basis.size());
    assert(count==basis.size());
    return 0;
}

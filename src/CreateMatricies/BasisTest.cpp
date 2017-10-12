#include "Basis.h"
#include <cstdio>
#include <cassert>

using namespace std;

//a series of tests
int main(){
    WedgeBasis basis(2,6,15);
    auto iter = basis.getIter();
    long long count = 0;
    do{
        printf("%d\n",count);
        assert(count==basis.rank(iter.getCurr()));
        assert(iter.getCurr()==basis.unrank(count));
        count++;
    }while(iter.next());
    printf("Check %lld==%lld\n",count,basis.size());
    assert(count==basis.size());
    return 0;
}

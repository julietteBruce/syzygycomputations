#include "Basis.h"

#include <vector>
#include <utility>
#include <cstdio>
#include <list>

using namespace std;

void printRow(const vector<long long>& row){
    for(size_t i=0;i<row.size();i++)
        printf("%lld ",row[i]);
    printf("\n");
}

bool isArtinianRedundent(int d, const basis_elem_t& elem,const WedgeBasis & basis){
    for(auto x : elem){
        for(auto y : basis.getBasicElements()[x]){
            if(y>=d)
                return true;
        }
    }
    return false;
}

vector<long long> computeImage(const basis_elem_t& elem,const WedgeBasis& codomain){
    vector<long long> imgList(elem.size());
    vector<monomial_t> imgElem = elem;
    monomial_t removed = imgElem[0];
    imgElem.erase(imgElem.begin());
    for(size_t i=0;i<elem.size();i++){
        imgList[i]=codomain.rank(imgElem);
        swap(imgElem[i],removed);
    }
    return imgList;
}


void outputMatrix(int n, int d, int p){
    WedgeBasis domainBasis(createWedgeBasis(n,d,p));
    WedgeBasis codomainBasis(createWedgeBasis(n,d,p-1));
    printf("%d %d %d\n",n,d,p);
    auto domainIter = domainBasis.getIter();
    do{
        basis_elem_t elem = domainIter.getCurr();
        printRow(computeImage(elem,codomainBasis));
    }while(domainIter.next());
}


int main(int argc, char** argv){
    int n,d,p;
    if(argc<4){
        printf("Please specify n, d, and p on the command line\n");
        return 1;
    }
    n=atoi(argv[1]);
    d=atoi(argv[2]);
    p=atoi(argv[3]);

    outputMatrix(n,d,p);
    return 0;
}

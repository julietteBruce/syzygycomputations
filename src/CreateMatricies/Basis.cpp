#include "Basis.h"
#include "Combinatorics.h"

#include <vector>
#include <utility>
#include <cstdio>
#include <map>
#include <list>
#include <algorithm>

using namespace std;

//note we only every care about small values, so this version is fine
template<typename T>
static T binom(T n, T k){
    T val = 1;
    for(int i=0;i<k;i++)
        val*=n-i;
    for(int i=1;i<k+1;i++)
        val/=i;
    return val;
}

/*
This gives the lexicographical rank of the basis element. Note that the basis elements are always strictly increasing
The basic idea is that at every step, we count how elements lie lexicographically between the current prefix
and the next prefix, which is given by the number of words starting with the current prefix, minus the number of
words starting with the next prefix.
*/
long long WedgeBasis::rank(const basis_elem_t& elem) const{
    long long ret = 0;
    int prev = -1;
    for(size_t i=0;i<elem.size();i++){
        //add the number of basis elements between the +1 case, and actuality
        ret += cachedSum(prev+1,elem[i],p-i-1);
        prev = elem[i];
    }
    return ret;
}

//Inefficient.
basis_elem_t WedgeBasis::unrank(long long r) const{
    basis_elem_t elem(p);
    int prev = -1;
    for(int k=0;k<p;k++){
        //find the cachedSum that makes this entry possible
        for(int j=prev+1;j<N;j++){
            long long s = cachedSum(prev+1,j,p-k-1);
            if(s>r){
                elem[k]=j-1;
                r -= cachedSum(prev+1,elem[k],p-k-1);
                break;
            }
            else if(s==r){
                elem[k]=j;
                r=0;
                break;
            }
        }
        prev = elem[k];
    }
    return elem;
}

vector<int> WedgeBasis::multidegree(const basis_elem_t& elem) const{
    vector<int> ret(m);
    for(monomial_t e : elem){
        for(int i=0;i<m;i++)
            ret[i]+=values[e][i];
    }
    return ret;
}

void WedgeBasis::multidegree(const basis_elem_t& elem, vector<int>& ret) const{
    for(int i=0;i<m;i++)
        ret[i]=0;
    for(monomial_t e : elem){
        for(int i=0;i<m;i++)
            ret[i]+=values[e][i];
    }
}

static vector<vector<long long> > createBinomialCache(int N){
    vector<vector<long long> > ret(N+1);
    for(int i=0;i<=N;i++){
        ret[i]=vector<long long>(i+1);
        ret[i][0]=1;
        ret[i][i]=1;
        for(int j=1;j<=i-1;j++)
            ret[i][j]=ret[i-1][j-1]+ret[i-1][j];
    }
    return ret;
}

WedgeBasis::WedgeBasis(int _m, vector<vector<int> > monomials, int _p) : m(_m),p(_p),N(monomials.size()),values(monomials),binomCache(createBinomialCache(N)){
}


long long WedgeBasis::cachedSum(int start, int stop, int k) const{
    return cachedBinom(N-start,k+1)-cachedBinom(N-stop,k+1);
}

SubBasis::SubBasis(const WedgeBasis &parent, const vector<int>& md){
    WedgeBasisIterator iter = parent.getIter();
    int counter = 0;
    int parent_counter = 0;
    do{
        if(is_below(parent.multidegree(iter.getCurr()),md)){
            mdMap[parent_counter]=counter;
            counter++;
        }
        parent_counter++;
    }while(iter.next());
}

long long SubBasis::convert_rank(long long old_rank) const{
    auto iter = mdMap.find(old_rank);
    if(iter!=mdMap.end())
        return iter->second;
    return -1;
}

WedgeBasis createWedgeBasis(int n, int d, int p){
    return WedgeBasis(n+1, createIntegerVectors(n,d),p);
}

WedgeBasis createProductWedgeBasis(int n1, int n2, int d1, int d2, int p){
    return WedgeBasis(n1+n2+2, createProductIntegerVectors(n1,n2,d1,d2),p);
}

WedgeBasis createReducedWedgeBasis(int n, int d, int p){
    vector<vector<int> > reducedIntegerVectors;
    for(vector<int> elem : createIntegerVectors(n,d)){
        if(all_of(elem.begin(),elem.end(),[d](int x){return x<d;})){
            reducedIntegerVectors.push_back(elem);
        }
    }
    return WedgeBasis(n+1,reducedIntegerVectors,p);
}

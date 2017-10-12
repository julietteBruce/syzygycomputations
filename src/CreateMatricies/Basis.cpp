#include "Basis.h"
#include "Combinatorics.h"

#include <vector>
#include <utility>
#include <cstdio>
#include <map>
#include <list>

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
This gives the lexagraphical rank of the basis element. Note that the basis elements are always strictly increasing
The basic idea is that at every step, we count how elements lie lexagraphically between the current prefix
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
    for(size_t k=0;k<p;k++){
        //find the cachedSum that makes this entry possible
        for(int j=prev+1;j<N;j++){
            int s = cachedSum(prev+1,j,p-k-1);
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

bool WedgeBasis::isArtinianTrivial(long long r) const{
    return trivialityCache[r];
}

vector<int> WedgeBasis::multidegree(const basis_elem_t& elem) const{
    vector<int> ret(n+1);
    for(int e : elem){
        for(int i=0;i<n+1;i++)
            ret[i]+=values[e][i];
    }
    return ret;
}

void WedgeBasis::multidegree(const basis_elem_t& elem, vector<int>& ret) const{
    for(int i=0;i<n+1;i++)
        ret[i]=0;
    for(int e : elem){
        for(int i=0;i<n+1;i++)
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

WedgeBasis::WedgeBasis(int _n, int _d, int _p) : n(_n), d(_d), p(_p), N(binom(n+d,n)), values(createIntegerVectors(n,d))/*,sumCache(N+1)*/,binomCache(createBinomialCache(N)),trivialityCache(binomCache[N][p]){
    WedgeBasisIterator iter = getIter();
    for(int i=0;i<size();i++){
        trivialityCache[i] = false;
        for(long long w : iter.getCurr()){
            for(int x : values[w]){
                if(x>=d)
                    trivialityCache[i] = true;
            }
        }
        iter.next();
    }
}

long long WedgeBasis::cachedSum(int start, int stop, int k) const{
    return cachedBinom(N-start,k+1)-cachedBinom(N-stop,k+1);
    //return sumCache[start][stop-start][k];
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

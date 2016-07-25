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

long long WedgeBasis::rank(const basis_elem_t& elem) const{
    long long ret = 0;
    int prev = -1;
    for(size_t i=0;i<elem.size();i++){
        ret += cachedSum(prev+1,elem[i],p-i-1);
        prev = elem[i];
    }
    return ret;
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

//the <long long> on the binom should really be automatic, but for some reason I get linking errors if I don't
//specify
WedgeBasis::WedgeBasis(int _n, int _d, int _p) : n(_n), d(_d), p(_p), N(binom<long long>(n+d,n)), values(createIntegerVectors(n,d)),sumCache(N+1),binomCache(N+1){

    //precompute a bunch of stuff
    for(int i=0;i<=N;i++){
        binomCache[i]=vector<long long>(i+1);
        binomCache[i][0]=1;
        binomCache[i][i]=1;
        for(int j=1;j<=i-1;j++)
            binomCache[i][j]=binomCache[i-1][j-1]+binomCache[i-1][j];
    }

    for(int start=0;start<N;start++){
        sumCache[start]=vector<vector<long long> >(N-start);
        for(int stop=start;stop<N;stop++){
            sumCache[start][stop-start]=vector<long long>(p);
            for(int k=0;k<p;k++){
                if(start==stop)
                    sumCache[start][stop-start][k]=0;
                else
                    sumCache[start][stop-start][k]=sumCache[start][stop-1-start][k]+cachedBinom(N-stop,k);
            }
        }
    }
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

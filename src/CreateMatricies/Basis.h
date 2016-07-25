
#include <vector>
#include <map>

typedef std::vector<long long> basis_elem_t;

long long binom(long long,long long);

//note that this iterator uses a different style from the c++ iterators, because the c++ iterators are needlessly
//complicated for our purposes
class WedgeBasisIterator{
public:
    WedgeBasisIterator(long long _N, int _k) : N(_N),k(_k), curr(_k) {
        for(int i=0;i<k;i++)
            curr[i]=i;
    }
    inline basis_elem_t getCurr() const {return curr;}
    inline bool next(){
        for(int i=1;i<=k;i++){
            if(curr[k-i]!=N-i){
                //something to advance
                curr[k-i]++;
                int val=curr[k-i]+(i-1);
                for(int j=1;j<i;j++,val--){
                    curr[k-j]=val;//curr[k-i]+(i-j);
                }
                return true;
            }
        }
        return false;
    }
    void reset(){
        for(int i=0;i<curr.size();i++)
            curr[i]=i;
    }

private:
    const long long N;
    const int k;
    std::vector<long long> curr;
};

class WedgeBasis{
public:
    WedgeBasis(int n, int d, int p);
    long long rank(const basis_elem_t& elem) const;
    std::vector<int> multidegree(const basis_elem_t& elem) const;
    void multidegree(const basis_elem_t& elem,std::vector<int>& ret) const;
    inline WedgeBasisIterator getIter() const{
        return WedgeBasisIterator(N,p);
    }
    inline long long size() const{
        return cachedBinom(N,p);
    }
    inline const std::vector<std::vector<int> > getBasicElements() const { return values; }
private:
    const long long n,d,p,N;
    const std::vector<std::vector<int> > values;
    std::vector<std::vector<std::vector<long long> > > sumCache;
    //a very weird function
    //we need this for all pairs 0<=start<=stop<N?
    //and all k<=p
    inline long long cachedSum(int start, int stop, int k) const{
        return sumCache[start][stop-start][k];
    }

    std::vector<std::vector<long long> > binomCache;
    inline long long cachedBinom(int n, int k) const{
        if(k>n || k<0)
            return 0;
        return binomCache[n][k];
    }
};

//TODO make this more intelligent
class SubBasis{
public:
    SubBasis(const WedgeBasis& parent,const std::vector<int>& md);
    long long convert_rank(long long old_rank) const;
    inline long long size() const{
        return mdMap.size();
    }
private:
    const std::vector<int> maxMultidegree;
    std::map<long long,long long> mdMap;
};




static inline bool is_below(const std::vector<int>& a, const std::vector<int>& b){
    for(size_t i = 0; i<a.size();i++){
        if(a[i]>b[i])
            return false;
    }
    return true;
}

#include <vector>
#include <utility>
#include <cstdio>
#include <map>
#include <list>

using namespace std;


typedef vector<int> basis_elem_t;

//note we only every care about small values, so this version is fine
int binom(int n, int k){
    int val = 1;
    for(int i=0;i<k;i++)
        val*=n-i;
    for(int i=1;i<k+1;i++)
        val/=i;
    return val;
}

//the number of ways to get an integer vector of n+1 things that sums to d is given by the stars and bars formula
//namely it's \binom{n+d}{n}
vector<vector<int> > createIntegerVectors(int n, int d){
    vector<vector<int> > ret(binom(n+d,n),vector<int>(n+1));
    //intialize values
    vector<int> partialSums(n+1);
    for(int i=0;i<n;i++)
        partialSums[i]=0;
    partialSums[n]=d;
    for(size_t i=0;i<ret.size();i++){
        //create a value out of a partial sum
        ret[i][0]=partialSums[0];
        for(int j=1;j<n+1;j++)
            ret[i][j]=partialSums[j]-partialSums[j-1];
        for(int j=n-1;j>=0;j--){
            //find an shiftable entry and shift to the right
            if(partialSums[j]!=partialSums[j+1]-1){
                partialSums[j]++;
                //reset everything past here
                for(int k=j+1;k<n;k++)
                    partialSums[k]=partialSums[j]+(k-j);
                break;
            }
        }
    }
    return ret;
}

//note that this iterator uses a different style from the c++ iterators, because the c++ iterators are needlessly
//complicated for our purposes
class WedgeBasisIterator{
public:
    WedgeBasisIterator(long long _N, int k) : N(_N), curr(k) {
        for(int i=0;i<k;i++)
            curr[i]=i;
    }
    inline vector<int> getCurr() const {return curr;}
    inline bool next(){
        int k = curr.size();
        for(int i=1;i<=k;i++){
            if(curr[k-i]!=N-i){
                //something to advance
                curr[k-i]++;
                for(int j=1;j<i;j++){
                    curr[k-j]=curr[k-i]+(i-j);
                }
                return true;
            }
        }
        return false;
    }
private:
    const long long N;
    vector<int> curr;
};

class WedgeBasis{
public:
    WedgeBasis(int n, int d, int p);
    long long rank(const basis_elem_t& elem) const;
    vector<int> multidegree(const basis_elem_t& elem) const;
    inline WedgeBasisIterator getIter() const{
        return WedgeBasisIterator(N,p);
    }
    inline long long size() const{
        return cachedBinom(N,p);
    }
private:
    const long long n,d,p,N;
    const vector<vector<int> > values;
    vector<vector<vector<long long> > > sumCache;
    //a very weird function
    //we need this for all pairs 0<=start<=stop<N?
    //and all k<=p
    inline long long cachedSum(int start, int stop, int k) const{
        return sumCache[start][stop-start][k];
    }

    vector<vector<long long> > binomCache;
    inline long long cachedBinom(int n, int k) const{
        if(k>n || k<0)
            return 0;
        return binomCache[n][k];
    }

};


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

WedgeBasis::WedgeBasis(int _n, int _d, int _p) : n(_n), d(_d), p(_p), N(binom(n+d,n)), values(createIntegerVectors(n,d)),sumCache(N+1),binomCache(N+1){

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


void printRow(const vector<int>& row){
    for(size_t i=0;i<row.size();i++)
        printf("%d,%d ",row[i],i%2 ? -1 : 1);
    printf("\n");
}

void printMultidegree(const vector<int>& md){
    printf("(");
    bool first = true;
    for(int curr : md){
        if(!first)
            printf(",");
        printf("%d",curr);
        first = false;
    }
    printf(")");
}

//space separated list
void printList(const list<int>& items){
    for(int i : items){
        printf("%d ", i);
    }
}

void printMDInfo(const map<vector<int>,list<int> >& mdInfo){
    printf("%d\n",mdInfo.size());
    for(pair<vector<int>,list<int> > md : mdInfo){
        printMultidegree(md.first);
        printf(": ");
        printList(md.second);
        printf("\n");
    }
}

vector<int> computeImage(const basis_elem_t& elem,const WedgeBasis& codomain){
    vector<int> imgList(elem.size());
    vector<int> imgElem = elem;
    int removed = imgElem[0];
    imgElem.erase(imgElem.begin());
    for(size_t i=0;i<elem.size();i++){
        imgList[i]=codomain.rank(imgElem);
        swap(imgElem[i],removed);
    }
    return imgList;
}


void outputMatrix(int n, int d, int p){
    WedgeBasis domainBasis(n,d,p);
    WedgeBasis codomainBasis(n,d,p-1);
    printf("%lld\n",domainBasis.size());
    auto domainIter = domainBasis.getIter();
    map<vector<int>,list<int> > domainMDs;
    int i=0;
    do{
        basis_elem_t elem = domainIter.getCurr();
        printRow(computeImage(elem,codomainBasis));
        domainMDs[domainBasis.multidegree(elem)].push_back(i);
        i++;
    }while(domainIter.next());
    printMDInfo(domainMDs);
    domainMDs.clear();
    map<vector<int>,list<int> > codomainMDs;
    auto codomainIter = codomainBasis.getIter();
    i=0;
    do{
        basis_elem_t elem = codomainIter.getCurr();
        codomainMDs[codomainBasis.multidegree(elem)].push_back(i);
        i++;
    }while(codomainIter.next());
    printMDInfo(codomainMDs);
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

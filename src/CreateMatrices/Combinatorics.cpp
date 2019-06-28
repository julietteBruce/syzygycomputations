#include "Combinatorics.h"

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
  An integer vector of length with the entries summing to d is equivalent to a
  choice of n numbers between 0 and d, namely the partial sums. In particular,
  to iterate over all integers vectors we instead iterate over the partial sums
  then to reconstruct the integer vector from the partial sums
 */
vector<vector<int> > createIntegerVectors(int n, int d){
    //the number of ways to get an integer vector of n+1 things that sums to d is given by
    //the stars and bars formula namely \binom{n+d}{n}
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
            if(partialSums[j]!=partialSums[j+1]){
                partialSums[j]++;
                //reset everything past here
                for(int k=j+1;k<n;k++)
                    partialSums[k]=partialSums[j];
                break;
            }
        }
    }
    return ret;
}

vector<vector<int> > createProductIntegerVectors(int n1, int n2, int d1, int d2){
  // the number of ways to get pairs of integer vectors of ni+1 things that sum to di, i=1,2
  vector<vector<int> > ret(binom(n1+d1,n1)*binom(n2+d2,n2),vector<int>(n1+n2+2));
  vector<vector<int> > vec1(createIntegerVectors(n1, d1));
  vector<vector<int> > vec2(createIntegerVectors(n2, d2));
  for (int i=0;i<vec1.size();i++){
    for (int j=0;j<vec2.size();j++){
      auto iter=copy(vec1[i].begin(), vec1[i].end(), ret[i*vec2.size()+j].begin());
      copy(vec2[j].begin(), vec2[j].end(), iter);
	}
  }
  return ret;
}

#include <functional>
#include <algorithm>

#include "ToricVariety.h"
#include "Combinatorics.h"

using namespace std;

vector<vector<int> > Pn::multidegrees(const vector<int>& degree, bool dedup) const {
    vector<vector<int> > initialList = createIntegerVectors(n,degree[0]);
    if(dedup){
        vector<vector<int> > ret;
        for(auto md : initialList){
            bool isIncreasing = true;
            int last = md[0];
            for(size_t i=1;i<md.size();i++){
                int curr = md[i];
                if(last>curr){
                    isIncreasing = false;
                    break;
                }
                last = curr;
            }
            if(isIncreasing)
                ret.push_back(md);
        }
        return ret;
    }
    else{
        return initialList;
    }
}

int sumBy(const vector<ToricVariety*>& values, function<int(ToricVariety*)> f){
    int ret = 0;
    for(ToricVariety* x : values){
        ret += f(x);
    }
    return ret;
}

int getDegreeSize(ToricVariety *v){
    return v->degreeSize;
}

int getMdSize(ToricVariety *v){
    return v->mdegSize;
}

Product::Product(const vector<ToricVariety* >& _tvs) :
    ToricVariety(sumBy(_tvs,getDegreeSize),sumBy(_tvs,getMdSize)),tvs(_tvs.size()){
    for(size_t i = 0;i<_tvs.size();i++)
        tvs[i] = unique_ptr<ToricVariety>(_tvs[i]);
}

vector<vector<int> > Product::multidegrees(const vector<int>& degree, bool dedup) const {
    vector<vector<vector<int> > > mdParts(tvs.size());
    size_t size = 0;
    auto degreeIter = degree.begin();
    for(size_t i=0;i<tvs.size();i++){
        int length = tvs[i]->degreeSize;
        vector<int> degreePiece(degreeIter,degreeIter+length);
        degreeIter+=length;
        mdParts[i] = tvs[i]->multidegrees(degreePiece, dedup);
        size += tvs[i]->mdegSize;
    }
    vector<vector<int> > ret;
    vector<size_t> indicies(tvs.size(),0);
    vector<int> currMd(size);
    while(true){
        //create the multidegree
        auto iter = currMd.begin();
        for(size_t i = 0;i<indicies.size();i++){
            vector<int> &currPiece = mdParts[i][indicies[i]];
            iter = copy(currPiece.begin(),currPiece.end(),iter);
        }
        ret.push_back(currMd);
        //increment the index
        for(size_t i = 0;i<indicies.size();i++){
            if(indicies[i]==mdParts[i].size()-1){
                indicies[i]=0;
                //if we are about the increment the last term, we are done
                if(i==indicies.size()-1)
                    return ret;
            }
            else{
                indicies[i]++;
                break;
            }
        }
    }
}

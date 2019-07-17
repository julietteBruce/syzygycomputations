#include <functional>
#include <algorithm>

#include "ToricVariety.h"
#include "Combinatorics.h"

using namespace std;

vector<vector<int> > Pn::multidegrees(const vector<int>& degree, bool dedup) const {
    if(degree[0]<0)
        return {};
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

static int sumBy(const vector<ToricVariety*>& values, function<int(ToricVariety*)> f){
    int ret = 0;
    for(ToricVariety* x : values){
        ret += f(x);
    }
    return ret;
}

static int getDegreeSize(ToricVariety *v){
    return v->degreeSize;
}

static int getMdSize(ToricVariety *v){
    return v->mdegSize;
}

static vector<unique_ptr<ToricVariety> > mkTvsVector(vector<ToricVariety*> tvs){
    vector<unique_ptr<ToricVariety> > ret (tvs.size());
    for(size_t i = 0;i<tvs.size();i++)
        ret[i] = unique_ptr<ToricVariety>(tvs[i]);
    return ret;
}

Product::Product(const vector<ToricVariety* >& _tvs) :
    ToricVariety(sumBy(_tvs,getDegreeSize),sumBy(_tvs,getMdSize)),tvs(mkTvsVector(_tvs)){
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

P1Bundle::P1Bundle(ToricVariety* _baseSpace, const std::vector<int>& _twistDegree)
    : ToricVariety(_baseSpace->degreeSize+1,_baseSpace->mdegSize+2),
      baseSpace(_baseSpace), twist(_twistDegree){
}

vector<vector<int> > P1Bundle::multidegrees(const vector<int>& degree, bool dedup) const{
    size_t baseDegreeSize = baseSpace->degreeSize;
    int d = degree[baseDegreeSize];
    vector<vector<vector<int> > > mds(d+1);
    vector<int> currentDegree(degree.begin(),degree.begin()+baseDegreeSize);
    for(int i=0;i<=d;i++){
        mds[i] = baseSpace->multidegrees(currentDegree,dedup);
        for(size_t j=0;j<baseDegreeSize;j++){
            currentDegree[j] -= twist[j];
        }
    }
    size_t totalMds = 0;
    for(auto& entry : mds)
        totalMds += entry.size();
    vector<vector<int> > ret(totalMds);
    size_t index = 0;
    for(int i=0;i<=d;i++){
        for(size_t j=0;j<mds[i].size();j++){
            ret[index].resize(baseSpace->mdegSize+2);
            auto iter = copy_n(mds[i][j].begin(),baseSpace->mdegSize,ret[index].begin());
            *(iter++)=i;
            *(iter++)=d-i;
            index++;
        }
    }
    return ret;
};

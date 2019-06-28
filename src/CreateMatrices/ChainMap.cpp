#include "ChainMap.h"
#include "Combinatorics.h"

#include <algorithm>
#include <sstream>

using namespace std;

static void printMatlabRow(FILE* fileInfo,long long rowNum, const vector<long long>& row){
    for(size_t i=0;i<row.size();i++)
        if(row[i]>=0)
            fprintf(fileInfo,"%lld %lld %d\n",rowNum,row[i]+1,i%2==0 ? 1 : -1);
}

static FILE *openOutputFile(const char* outputDir,const char* namePrefix,const vector<int>& md){
    stringstream name;
    name << outputDir << "/";
    name << namePrefix;
    for(int x : md)
        name  << "_" << x;
    name << ".dat";
    return fopen(name.str().c_str(),"w");
}

vector<long long> ChainMap::computeTrimmedImage(const basis_elem_t& elem,const vector<int>&){
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

void ChainMap::createMatrixFiles(const char* directory,
                                 const char* namePrefix,
                                 const vector<vector<int> >& mdInfo){
    if(domain.size()==0 || mdInfo.empty())
        return;
    vector<FILE*> files(mdInfo.size(),(FILE*)NULL);
    auto domainIter = domain.getIter();
    vector<int> rowMd(mdInfo[0].size());
    vector<long long> rowNums(mdInfo.size(),1LL);
    //files.assign(mdInfo.size(),NULL);
    //rowNums.assign(mdInfo.size(),1);
    do{
        basis_elem_t elem = domainIter.getCurr();
        domain.multidegree(elem,rowMd);
        for(size_t i=0;i<mdInfo.size();i++){
            if(is_below(rowMd,mdInfo[i])){
                vector<long long> trimmedImage = computeTrimmedImage(elem,mdInfo[i]);
                if(!all_of(trimmedImage.begin(),trimmedImage.end(),[](long long i){return i<0;})){
                    if(!files[i]){
                        if(!(files[i] = openOutputFile(directory,namePrefix,mdInfo[i]))){
                            fprintf(stderr,"Could not open output file\n");
                            exit(1);
                        }
                    }
                    printMatlabRow(files[i],rowNums[i]++,trimmedImage);
                }
            }
        }
    }while(domainIter.next());

}

vector<long long> PnChainMap::computeTrimmedImage(const basis_elem_t& elem, const vector<int>& md){
    vector<int> wedgePartMd = domain.multidegree(elem);
    //vector<int> tensorPartMd(wedgePartMd.size());
    //the largest value that the removed part can have so that the tensor part ends up nontrivial in the image
    vector<int> removedPartLimit(wedgePartMd.size());
    for(size_t i = 0;i<wedgePartMd.size();i++){
        removedPartLimit[i] = embeddingDegree-(md[i] - wedgePartMd[i])-1;
    }

    vector<long long> imgList(elem.size());
    vector<monomial_t> imgElem = elem;
    monomial_t removed = imgElem[0];
    imgElem.erase(imgElem.begin());

    for(size_t i=0;i<elem.size();i++){
        if(is_below(domain.getMonomials()[removed],removedPartLimit)){
            imgList[i]=codomain.rank(imgElem);
        }
        else{
            imgList[i]=-1;
        }
        swap(imgElem[i],removed);
    }
    return imgList;
}

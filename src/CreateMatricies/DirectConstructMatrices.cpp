#include "Basis.h"
#include "Combinatorics.h"

#include <cstdio>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <map>
#include <utility>

#include <sys/stat.h>

using namespace std;

void printMatlabRow(FILE* fileInfo,long long rowNum, const vector<long long>& row){
    for(size_t i=0;i<row.size();i++)
        if(row[i]>=0)
            fprintf(fileInfo,"%lld %lld %d\n",rowNum,row[i],i%2==0 ? 1 : -1);
}

bool isTensorPartTrivial(const basis_elem_t& elem,vector<int> md, int d,const WedgeBasis& basis){
    vector<int> elemMd = basis.multidegree(elem);
    for(int i=0;i<elemMd.size();i++)
        if(md[i]-elemMd[i] >=d)
            return true;
    return false;
}

vector<long long> computeImage(const basis_elem_t& elem,const WedgeBasis& codomain){
    vector<long long> imgList(elem.size());
    vector<long long> imgElem = elem;
    long long removed = imgElem[0];
    imgElem.erase(imgElem.begin());
    for(size_t i=0;i<elem.size();i++){

        //if(!isTrivial(imgElem,_,codomain)){
        imgList[i]=codomain.rank(imgElem);
            //}
        swap(imgElem[i],removed);
    }
    return imgList;
}


//returns true if there's still something in the image
bool trimImage(vector<long long>& img, vector<int> md, int d, const WedgeBasis& codomain){
    bool ret = false;
    for(int i=0;i<img.size();i++){
        long long term = img[i];
        basis_elem_t termElem = codomain.unrank(term);
        if(codomain.isArtinianTrivial(term) || isTensorPartTrivial(termElem,md,d,codomain)){
            img[i] = -1;
        }
        else{
            ret = true;
        }
    }
    return ret;
}

FILE *openOutputFile(const char* outputDir,const vector<int>& md){
    stringstream name;
    name << outputDir << "/";
    name << "multidegree";
    for(int x : md)
        name  << "_" << x;
    name << ".dat";
    return fopen(name.str().c_str(),"w");
}

void outputMatrix(int n, int d, int p,const map<vector<int>,FILE*>& mdInfo){
    WedgeBasis domainBasis(n,d,p);
    WedgeBasis codomainBasis(n,d,p-1);
    auto domainIter = domainBasis.getIter();
    vector<int> rowMd(n+1);
    map<vector<int>,long long> rowNums;
    for(const auto &info : mdInfo){
        rowNums[info.first] = 1;
    }
    do{
        basis_elem_t elem = domainIter.getCurr();
        domainBasis.multidegree(elem,rowMd);
        bool exists = false;
        vector <long long> img;
        for(const auto &info : mdInfo){
            if(is_below(rowMd,info.first)){
                if(!exists){
                    img = computeImage(elem,codomainBasis);
                    exists = true;
                }
                vector<long long> tmpImg = img;
                if(trimImage(tmpImg,info.first,d,codomainBasis)){
                    printMatlabRow(info.second,rowNums[info.first]++,tmpImg);
                }
            }
        }
    }while(domainIter.next());
}

int main(int argc, char ** argv){
    if(argc<6){
        fprintf(stderr,"Not enough arguments\n");
        fprintf(stderr,"usage: ./DirectConstructMatricies n d p md_list out_dir\n");
        return -1;
    }
    int n=atoi(argv[1]);
    int d=atoi(argv[2]);
    int p=atoi(argv[3]);

    map<vector<int>,FILE*>allMds;

    const char* outputDir = argv[5];
    ifstream mdListFile (argv[4]);
    mdListFile >> ws;
    while(!mdListFile.eof()){
        vector<int> md(n+1);
        for(int i=0;i<=n;i++){
            mdListFile >> md[i];
        }
        allMds[md]=openOutputFile(outputDir,md);
        mdListFile >> ws;
    }
    if(!allMds.empty())
        outputMatrix(n,d,p,allMds);
    return 0;
}

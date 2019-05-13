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

vector<long long> computeTrimmedImage(const basis_elem_t& elem, vector<int> md, int d ,const WedgeBasis& domain,const WedgeBasis& codomain){
    vector<int> wedgePartMd = domain.multidegree(elem);
    vector<int> tensorPartMd(wedgePartMd.size());
    //the largest value that the removed part can have so that the tensor part ends up nontrivial in the image
    vector<int> removedPartLimit(wedgePartMd.size());
    for(size_t i = 0;i<wedgePartMd.size();i++){
        tensorPartMd[i]=md[i]-wedgePartMd[i];
        removedPartLimit[i] = d-tensorPartMd[i]-1;
        printf("%d ",tensorPartMd[i]);
    }
    printf("\n");

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

FILE *openOutputFile(const char* outputDir,const vector<int>& md){
    stringstream name;
    name << outputDir << "/";
    name << "multidegree";
    for(int x : md)
        name  << "_" << x;
    name << ".dat";
    return fopen(name.str().c_str(),"w");
}

void outputMatrix(const char* outputDir,int n, int d, int p,const vector<vector<int> >& mdInfo){
    vector<FILE*> files;
    WedgeBasis domainBasis(createReducedWedgeBasis(n,d,p));
    WedgeBasis codomainBasis(createReducedWedgeBasis(n,d,p-1));
    auto domainIter = domainBasis.getIter();
    vector<int> rowMd(n+1);
    vector<long long> rowNums;
    files.assign(mdInfo.size(),NULL);
    rowNums.assign(mdInfo.size(),1);
    if(domainBasis.size()==0)
      return;
    do{
        basis_elem_t elem = domainIter.getCurr();
        domainBasis.multidegree(elem,rowMd);

        //bool exists = false;
        //vector <long long> img;
        for(size_t i=0;i<mdInfo.size();i++){
            if(is_below(rowMd,mdInfo[i])){
                vector<long long> trimmedImage = computeTrimmedImage(elem,mdInfo[i],d,domainBasis,codomainBasis);
                if(!all_of(trimmedImage.begin(),trimmedImage.end(),[](long long i){return i<0;})){
                    if(!files[i]){
                        if(!(files[i] = openOutputFile(outputDir,mdInfo[i]))){
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

int main(int argc, char ** argv){
    if(argc<6){
        fprintf(stderr,"Not enough arguments\n");
        fprintf(stderr,"usage: ./DirectConstructMatricies n d p md_list out_dir\n");
        return -1;
    }
    int n=atoi(argv[1]);
    int d=atoi(argv[2]);
    int p=atoi(argv[3]);

    vector<vector<int> >allMds;

    const char* outputDir = argv[5];
    ifstream mdListFile (argv[4]);
    mdListFile >> ws;
    while(!mdListFile.eof()){
        vector<int> md(n+1);
        for(int i=0;i<=n;i++){
            mdListFile >> md[i];
        }
        allMds.push_back(md);
        mdListFile >> ws;
    }
    if(!allMds.empty())
        outputMatrix(outputDir,n,d,p,allMds);
    return 0;
}

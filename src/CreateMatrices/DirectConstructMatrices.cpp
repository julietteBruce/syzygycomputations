#include "Basis.h"
#include "Combinatorics.h"
#include "ChainMap.h"
#include "ToricVariety.h"

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

int main(int argc, char ** argv){
    if(argc<6){
        fprintf(stderr,"Not enough arguments\n");
        fprintf(stderr,"usage: ./DirectConstructMatricies n d p out_dir md_list|(q k)\n");
        fprintf(stderr,"or\n");
        fprintf(stderr,"usage: ./DirectConstructMatricies --product n1 n2 d1 d2 p out_dir md_list|(q k1 k2)\n");
        return -1;
    }

    ToricVariety *var;
    int p;
    int q=-1;
    vector<int> d;
    vector<int> b;
    ChainMap *chainMap = NULL;
    const char *outputDir = NULL;
    const char *mdListFilename = NULL;
    if(strcmp(argv[1],"--product")==0){
        int n1=atoi(argv[2]);
        int n2=atoi(argv[3]);
        d = {atoi(argv[4]),
             atoi(argv[5])};
        p=atoi(argv[6]);
        outputDir = argv[7];
        if(argc>9){
            q = atoi(argv[8]);
            b = {atoi(argv[9]),
                 atoi(argv[10])};
        }
        else{
            mdListFilename = argv[8];
        }
        var = new Product({new Pn(n1),new Pn(n2)});
        chainMap = new ToricChainMap(*var,d,p);
    }
    else{
        int n=atoi(argv[1]);
        var = new Pn(n);
        d = {atoi(argv[2])};
        p = atoi(argv[3]);
        outputDir = argv[4];
        chainMap = new PnChainMap(n,d[0],p);
        if(argc>6){
            q = atoi(argv[5]);
            b = {atoi(argv[6])};
        }
        else{
            mdListFilename = argv[5];
        }
    }
    //Create the list of multidegrees to process
    vector<vector<int> > allMds;
    if(mdListFilename){
        ifstream mdListFile (mdListFilename);
        mdListFile >> ws;
        while(!mdListFile.eof()){
            vector<int> md(var->mdegSize);
            for(size_t i=0;i<var->mdegSize;i++){
                mdListFile >> md[i];
            }
            allMds.push_back(md);
            mdListFile >> ws;
        }
    }
    else{
        vector<int> totalDegree(var->degreeSize);
        for(size_t i=0;i<totalDegree.size();i++)
            totalDegree[i] = (p+q)*d[i]+b[i];
        allMds = var->multidegrees(totalDegree,true);
    }
    //create the matrices
    if(!allMds.empty())
        chainMap->createMatrixFiles(outputDir,"multidegree",allMds);
    delete chainMap;
    delete var;
    return 0;
}

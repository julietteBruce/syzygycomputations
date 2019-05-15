#include "Basis.h"
#include "Combinatorics.h"
#include "ChainMap.h"

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


    vector<vector<int> >allMds;
    ChainMap *chainMap = NULL;
    const char* outputDir = NULL;
    if(strcmp(argv[1],"--product")==0){
        int n1=atoi(argv[2]);
        int n2=atoi(argv[3]);
        int d1=atoi(argv[4]);
        int d2=atoi(argv[5]);
        int p=atoi(argv[6]);
        chainMap = new ProductChainMap(n1,n2,d1,d2,p);
        outputDir = argv[7];
        if(argc>9){
            int q = atoi(argv[8]);
            int k1 = atoi(argv[9]);
            int k2 = atoi(argv[10]);
            allMds = createProductIntegerVectors(n1,n2,p*d1+k1+q*d1,p*d2+k2+q*d2);
        }
        else{
            ifstream mdListFile (argv[8]);
            mdListFile >> ws;
            while(!mdListFile.eof()){
                vector<int> md(n1+n2+2);
                for(int i=0;i<=n1+n2+2;i++){
                    mdListFile >> md[i];
                }
                allMds.push_back(md);
                mdListFile >> ws;
            }
        }
    }
    else{
        int n=atoi(argv[1]);
        int d=atoi(argv[2]);
        int p=atoi(argv[3]);
        outputDir = argv[4];
        chainMap = new PnChainMap(n,d,p);
        if(argc>6){
            int q = atoi(argv[6]);
            int k = atoi(argv[7]);
            allMds = createIntegerVectors(n,p*d+k+q*d);
        }
        else{
            ifstream mdListFile (argv[5]);
            mdListFile >> ws;
            while(!mdListFile.eof()){
                vector<int> md(n+1);
                for(int i=0;i<=n;i++){
                    mdListFile >> md[i];
                }
                allMds.push_back(md);
                mdListFile >> ws;
            }
        }
    }

    if(!allMds.empty())
        chainMap->createMatrixFiles(outputDir,"multidegree",allMds);
    delete chainMap;
    return 0;
}

#include "Basis.h"
#include "Combinatorics.h"
#include "Util.h"
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
    if(argc<5){
        fprintf(stderr,"Not enough arguments\n");
        fprintf(stderr,"usage: ./DirectConstructMatricies n d p out_dir [md_list]\n");
        return -1;
    }
    int n=atoi(argv[1]);
    int d=atoi(argv[2]);
    int p=atoi(argv[3]);
    //For testing purposes
    int k=0;
    int q=1;

    vector<vector<int> >allMds;

    const char* outputDir = argv[4];
    if(argc<6){
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

    PnChainMap chainMap(n,d,p);
    if(!allMds.empty())
        chainMap.createMatrixFiles(outputDir,"multidegree",allMds);
    return 0;
}

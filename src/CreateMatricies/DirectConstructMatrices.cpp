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
    PnChainMap chainMap(n,d,p);

    if(!allMds.empty())
        chainMap.createMatrixFiles(outputDir,"multidegree",allMds);
    return 0;
}

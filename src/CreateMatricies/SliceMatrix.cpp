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

#include <sys/stat.h>

using namespace std;

/*
  Note In theory we could be more efficient by computing all the slices at once,
  but this causes issues with the limit on the number of open files, so instead
  we compute the slices sequentially.
*/


class RowSource{
public:
    RowSource (const char* file) : subBasis(NULL),md(NULL),f(file){
        f >> n >> d >> p;
        f >> ws;
        domainBasis = new WedgeBasis(n,d,p);
        codomainBasis = new WedgeBasis(n,d,p-1);
        rowIter = new WedgeBasisIterator(domainBasis->getIter());
    }
    ~RowSource(){
        delete rowIter;
        delete codomainBasis;
        delete domainBasis;
        if(subBasis)
            delete subBasis;
        if(md)
            delete md;
    }
    inline int getN() const {return n;};
    inline int getD() const {return d;};
    inline int getP() const {return p;};
    bool nextRow(vector<long long>& row){
        vector<int> rowMd;
        string str;
        bool ret = true;
        do{
            if(!ret)
                return ret;
            rowMd = domainBasis->multidegree(rowIter->getCurr());
            getline(f,str);
            if(!rowIter->next())
                ret = false;
        }while(md && !is_below(rowMd,*md));
        //parse a line, loads of fun
        for(int i=0;i<p;i++){
            size_t idx=0;
            row[i] = stoll(str,&idx);
            if(subBasis)
                row[i] = subBasis->convert_rank(row[i]);
            str = str.substr(idx+1);
        }
        return ret;
    }

    void setMultidegree(const vector<int>& _md){
        md = new vector<int>(_md);
        subBasis = new SubBasis(*codomainBasis,_md);
    }

    bool hasRows(){
        return !subBasis || subBasis->size()!=0;
    }

private:
    WedgeBasis *domainBasis,*codomainBasis;
    SubBasis *subBasis;
    vector<int> *md;
    WedgeBasisIterator *rowIter;
    int n,d,p;
    ifstream f;
};

//output something matlab can understand
void formatOutputMatlab(const char* filename, RowSource& source){
    long long rowNum = 1;//mathlab 1 indexes
    vector<long long> row(source.getP());
    FILE *out = NULL;
    bool hasRows = true;
    while(hasRows){
        hasRows = source.nextRow(row);
        if(!out){
            out = fopen(filename,"w");
            if(!out){
                printf("error opening output file %s\n",filename);
                return;
            }
        }
        for(int i=0;i<row.size();i++){
            fprintf(out,"%lld %lld\t%d\n",rowNum ,row[i]+1, i%2 ? -1 : 1);
        }
        ++rowNum;
    }
    if(out)
        fclose(out);
}

int main(int argc, char ** argv){
    if(argc<4){
        printf("Not enough arguments\n");
        return -1;
    }

    RowSource source(argv[1]);
    int q = atoi(argv[2]);
    string directory(argv[3]);
    mkdir(argv[3],S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
    int totalDegree=source.getD()*source.getP()+source.getD()*q;
    vector<vector<int> > allMds(createIntegerVectors(source.getN(),totalDegree));
    for(const vector<int>& md : allMds){
        bool skip = false;
        for(size_t i=0;i<md.size()-1;i++){
            if(md[i]>md[i+1]){
                skip = true;
                break;
            }
        }
        if(skip)
            continue;
        RowSource realSource(argv[1]);
        realSource.setMultidegree(md);
        if(realSource.hasRows()){
            stringstream name;
            name << directory << "/";
            name << "multidegree";
            for(int x : md)
                name  << "_" << x;
            name << ".dat";
            formatOutputMatlab(name.str().c_str(),realSource);
        }
    }


    return 0;
}

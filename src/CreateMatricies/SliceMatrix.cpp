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
#include <set>

#include <sys/stat.h>

using namespace std;

/*
  Note In theory we could be more efficient by computing all the slices at once,
  but this causes issues with the limit on the number of open files, so instead
  we compute the slices sequentially.
*/


class RowSource{
public:
    RowSource (const char* file) : done(false),md(NULL),f(file){
        f >> n >> d >> p;
        f >> ws;
        domainBasis = new WedgeBasis(n,d,p);
        rowIter = new WedgeBasisIterator(domainBasis->getIter());
        rowMd.resize(n+1);
    }
    ~RowSource(){
        delete rowIter;
        delete domainBasis;
        if(md)
            delete md;
    }
    inline int getN() const {return n;};
    inline int getD() const {return d;};
    inline int getP() const {return p;};
    /* Retrives the next row into the output parameter row.
     * Returns true if a row was successfully retrieved */
    bool nextRow(vector<long long>& row){
        //vector<int> rowMd(n+1);
        string str;
        bool first = true;
        do{
            if(done)
                return false;
            domainBasis->multidegree(rowIter->getCurr(),rowMd);
            if(first)
                first = false;
            else
                f.ignore(numeric_limits<streamsize>::max(),'\n');//getline(f,str);
            if(!rowIter->next())
                done = true;
        }while(md && !is_below(rowMd,*md));
        getline(f,str);
        //parse a line, loads of fun
        for(int i=0;i<p;i++){
            size_t idx=0;
            row[i] = stoll(str,&idx);
            str = str.substr(idx+1);
        }
        return true;
    }

    void setMultidegree(const vector<int>& _md){
        md = new vector<int>(_md);
    }

private:
    bool done;
    WedgeBasis *domainBasis;
    vector<int> *md;
    vector<int> rowMd;
    WedgeBasisIterator *rowIter;
    int n,d,p;
    ifstream f;
};

/* This function writes the rows in source to the file given by filename.
 * So that it doesn't create empty files, it doesn't open the file until we get at least one row
 * from source to write out */
void formatOutputMatlab(const char* filename, RowSource& source){
    long long rowNum = 1;//mathlab 1 indexes
    vector<long long> row(source.getP());
    FILE *out = NULL;
    while(source.nextRow(row)){
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
        fprintf(stderr,"Not enough arguments\n");
        fprintf(stderr,"usage: ./SliceMatrix inputFile q k outputDir\n");
        return -1;
    }

    RowSource source(argv[1]);
    int q = atoi(argv[2]);
    int k = atoi(argv[3]);
    string directory(argv[4]);
    mkdir(argv[3],S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
    set<vector<int> > allMds;
    WedgeBasis domainBasis(source.getN(),source.getD(),source.getP());
    WedgeBasisIterator domainIter = domainBasis.getIter();
    vector<vector<int> > otherMds(createIntegerVectors(source.getN(),source.getD()*q+k));
    do{
        vector<int> firstMd(domainBasis.multidegree(domainIter.getCurr()));
        for(const vector<int>& secondMd : otherMds){
            vector<int> md(firstMd);
            for(int i=0;i<secondMd.size();i++)
                md[i] += secondMd[i];
            bool skip = false;
            for(size_t i=0;i<md.size()-1;i++){
                if(md[i]>md[i+1]){
                    skip = true;
                    break;
                }
            }
            if(!skip)
                allMds.insert(md);
        }
    }while(domainIter.next());

    for(const vector<int>& md : allMds){
        RowSource realSource(argv[1]);
        realSource.setMultidegree(md);
        stringstream name;
        name << directory << "/";
        name << "multidegree";
        for(int x : md)
            name  << "_" << x;
        name << ".dat";
        formatOutputMatlab(name.str().c_str(),realSource);
    }


    return 0;
}

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

    /*
    void runListeners(){
        vector<function<void(long long, vector<long long>)> > real_listeners;
        vector<vector<int> > mds;
        //create listeners that listen to the unmodified line
        //using subbasis converters
        for(auto elem : listeners){
            auto md = elem.first;
            mds.push_back(md);
            auto callback = elem.second;
            auto sub_basis = SubBasis(*codomainBasis,md);
            real_listeners.push_back([sub_basis,callback](long long rowNum,vector<long long> line){
                    std::vector<long long> mapped_line(line.size());
                    transform(line.begin(),line.end(),mapped_line.begin(),
                              [&sub_basis](long long x) -> long long {return sub_basis.convert_rank(x);});
                    callback(rowNum,mapped_line);
                });
        }

        auto rowIter = domainBasis->getIter();
        vector<int> rowIndicies(real_listeners.size());

        for(size_t i=0;i<nrows;i++){
            vector<int> rowMd = domainBasis->multidegree(rowIter.getCurr());
            rowIter.next();
            string str;
            getline(f,str);
            //parse a line, loads of fun
            std::vector<long long> line(p);
            for(size_t i=0;i<p;i++){
                size_t idx=0;
                line[i] = stoll(str,&idx);
                str = str.substr(idx+1);
            }
            //now feed it to the listeners
            for(size_t j=0;j<real_listeners.size();j++){
                if(is_below(rowMd,mds[j]))
                    real_listeners[j](rowIndicies[j]++,line);
            }
        }
    }
    */
private:
    WedgeBasis *domainBasis,*codomainBasis;
    SubBasis *subBasis;
    vector<int> *md;
    WedgeBasisIterator *rowIter;
    int n,d,p;
    ifstream f;
};

//output something matlab can understand
void formatOutputMatlab(FILE *out, RowSource& source){
    long long rowNum = 1;//mathlab 1 indexes
    vector<long long> row(source.getP());
    while(source.nextRow(row)){
        for(int i=0;i<row.size();i++){
            fprintf(out,"%lld %lld\t%d\n",rowNum ,row[i]+1, i%2 ? -1 : 1);
        }
        ++rowNum;
    }
}

int main(int argc, char ** argv){
    if(argc<3){
        printf("Not enough arguments\n");
        return -1;
    }

    RowSource source(argv[1]);
    int q = atoi(argv[2]);
    char * end = strrchr(argv[1],'.');
    string directory(argv[1],end-argv[1]);
    directory += "_";
    directory += argv[2];
    mkdir(directory.c_str(),S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
    int totalDegree=source.getD()*source.getP()+source.getD()*q;
    vector<vector<int> > allMds(createIntegerVectors(source.getN(),totalDegree));
    for(auto md : allMds){
        RowSource realSource(argv[1]);
        realSource.setMultidegree(md);
        stringstream name;
        name << directory << "/";
        name << "multidegree";
        for(int x : md)
            name  << "_" << x;
        name << ".dat";
        FILE *output = fopen(name.str().c_str(),"w");
        formatOutputMatlab(output,realSource);
        fclose(output);
    }


    return 0;
}

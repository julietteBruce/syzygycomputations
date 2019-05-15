#ifndef SYZYGY_CHAIN_MAP_H_
#define SYZYGY_CHAIN_MAP_H_ 1

#include <vector>
#include <utility>

#include "Basis.h"

class ChainMap{
public:
    virtual ~ChainMap() = default;
    ChainMap(const WedgeBasis& _domain,const WedgeBasis& _codomain) : domain(_domain),codomain(_codomain){};
    virtual std::vector<long long> computeTrimmedImage(const basis_elem_t& elem,
                                                       const std::vector<int>& md);
    void createMatrixFiles(const char* directory,
                           const char* namePrefix,
                           const std::vector<std::vector<int> >& mdInfo);
protected:
    const WedgeBasis domain;
    const WedgeBasis codomain;
};

class PnChainMap : public ChainMap{
public:
    PnChainMap(int n, int d, int p) :
        ChainMap(createReducedWedgeBasis(n,d,p),
                 createReducedWedgeBasis(n,d,p-1)),
        embeddingDegree(d) {};
    std::vector<long long> computeTrimmedImage(const basis_elem_t& elem,
                                                const std::vector<int>& md);
protected:
    int embeddingDegree;
};

class ProductChainMap : public ChainMap{
public:
    ProductChainMap(int n1, int n2, int d1, int d2, int p) :
        ChainMap(createProductWedgeBasis(n1,n2,d1,d2,p),
                 createProductWedgeBasis(n1,n2,d1,d2,p-1)) {}
};

#endif

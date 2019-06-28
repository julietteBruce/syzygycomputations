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

class ToricChainMap : public ChainMap {
public:
    ToricChainMap(ToricVariety& tv, const std::vector<int>& degree, int p)
        : ChainMap(createToricWedgeBasis(tv,degree,p),
                   createToricWedgeBasis(tv,degree,p)) {}
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

#endif

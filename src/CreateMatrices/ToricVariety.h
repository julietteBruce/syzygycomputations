#ifndef SYZYGY_TORIC_VARIETY_H_
#define SYZYGY_TORIC_VARIETY_H_ 1

#include <vector>
#include <memory>

/**
 * This class represents the data that the matrix construction code needs
 * about a toric variety to do it's work. For each toric variety, we need
 * degreeSize = rank Cl(X) = rank Pic(X)
 * mdegSize = |\Sigma(1)| = number of rays in \Sigma
 * a function multidegrees that given a element of Cl(X) represented as an
 * integer vector, gives the corresponding vectors representing the appropriate
 *  multidegrees.
 * In theory, this is just about computing interior lattice points for the
 * polytope corresponding to the divisior, but that's hard, so instead
 * we just deal with a few special cases.
 */
class ToricVariety {
private:
public:
    ToricVariety(size_t _degreeSize,size_t _mdegSize)
        : degreeSize(_degreeSize),mdegSize(_mdegSize){};
    virtual ~ToricVariety() = default;
    const size_t degreeSize;
    const size_t mdegSize;
    virtual std::vector<std::vector<int> > multidegrees(const std::vector<int>& deg, bool dedup) const = 0;
};

class Pn : public ToricVariety{
protected:
    const int n;
public:
    Pn(int _n) : ToricVariety(1,_n+1),n(_n) {};
    virtual std::vector<std::vector<int> > multidegrees(const std::vector<int>& deg, bool dedup) const;
};

class Product : public ToricVariety {
protected:
    std::vector<std::unique_ptr<ToricVariety> > tvs;
public:
    Product(const std::vector<ToricVariety* >& tvs);
    virtual std::vector<std::vector<int> > multidegrees(const std::vector<int>& deg, bool dedup) const;
};

/** A P1 bundle over a toric variety
 */
class P1Bundle : public ToricVariety {
protected:
    const std::unique_ptr<ToricVariety> baseSpace;
    const std::vector<int> twist;//the cartier data of the divisor of the bundle
public:
    P1Bundle(ToricVariety* baseSpace, const std::vector<int>& twistDegree);
    virtual std::vector<std::vector<int> > multidegrees(const std::vector<int>& deg, bool dedup) const;
};

#endif

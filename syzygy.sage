
import operator

#we represent an element of S_d by a n-tuple of integers adding up to d
#compute the basis for S_d
def compute_basis(n,d):
    return map(tuple,IntegerVectors(d,n+1))

#a wedge is represented by taking subsets of the basis
def wedge_basis(p,basis):
    if p>=0:
        return Subsets(basis,p);
    else:
        return Set()

#a tensor is a cartesian product on the basis
def tensor_basis(basis1,basis2):
    return CartesianProduct(basis1,basis2);

#put everything together for a basis of an entry in the complex
def complex_basis(n,d,p,q):
    sd_basis = compute_basis(n,d)
    sqd_basis = compute_basis(n,q*d)
    return tensor_basis(wedge_basis(p,sd_basis),sqd_basis)

def tuple_sum(a,b):
    return tuple([x+y for (x,y) in zip(a,b)])

def multidegree(elem):
    (ms,f) = elem
    md = f
    for m in ms:
        md = tuple_sum(md,m)
    return md

def normalize_elem(elem):
    (ms,f) = elem
    ms_list = sorted(ms);
    return (tuple(ms_list),f)

def normalize_md(md):
    md_list = list(md)
    md_list.sort()
    return tuple(md_list)


def wedge_part(n,d,p):
    sd_basis = compute_basis(n,d);
    def normalize_wedge(elem):
        return tuple(sorted(elem))
    return map(normalize_wedge,wedge_basis(p,sd_basis)); #for various boring reasons, we're going to just use the list

#computes the image of a basis element, the returned value is a dict from basis elements in the codomain to a sign
def compute_image(elem,codomain_basis_ids):
    def without(lst,i):
        new_list = list(lst)
        del new_list[i]
        return tuple(new_list)
    return {codomain_basis_ids[without(elem,i)]:(1 if i%2==0 else -1) for i in [0..(len(elem)-1)]}


#for some n,d,q, it returns (as indicies into the wedge basis) the applicable elements of the wedge basis
def slice_multidegree(wedge_basis,n,d,q):
    def wedge_multidegree(elem):
        md = (0,)*(n+1)
        for m in elem:
            md = tuple_sum(md,m)
        return md
    ret = dict()
    for (i,elem) in enumerate(wedge_basis):
        wedge_md = wedge_multidegree(elem);
        for f in compute_basis(n,q*d):
            md = tuple_sum(wedge_md,f);
            if md in ret:
                ret[md].append(i);
            else:
                ret[md] = [i];
    return ret

#compute the rank of the \delta_p map in row q
def compute_rank(n,d,p,q):
    #first compute the matrix on the basis of wedges
    domain_basis = wedge_part(n,d,p)
    codomain_basis = wedge_part(n,d,p-1)
    codomain_basis_ids = {elem:n for (n,elem) in enumerate(codomain_basis)}
    rows = map(lambda x : compute_image(x,codomain_basis_ids),domain_basis)

    #for each multidegree, compute the rows/columns of the slice
    #We don't actually need the codomain_slices part, it won't change the rank
    #but using it reduces the size of the largest matricies
    domain_slices = slice_multidegree(domain_basis,n,d,q)
    codomain_slices = slice_multidegree(codomain_basis,n,d,q+1)

    counts = dict()
    ranks = dict()
    for (md,curr_dom_slice) in domain_slices.iteritems():
        normalized = normalize_md(md)
        if md not in codomain_slices:
            continue;
        curr_codom_slice = codomain_slices[md]

        if normalized in counts:
            counts[normalized]+=1
        else:
            counts[normalized]=1
        if normalized!=md:
            continue;
        def check_image(domainId,codomainId):
            realDomainId = curr_dom_slice[domainId];
            realCodomainId = curr_codom_slice[codomainId]
            if realCodomainId in rows[realDomainId]:
                return rows[realDomainId][realCodomainId]
            else:
                return 0

        ranks[md] = matrix(ZZ,len(curr_dom_slice),len(curr_codom_slice),check_image).rank(algorithm='linbox')

    r=0
    for md in counts:
        r+=counts[md]*ranks[md]
    return r

#computes a single betti number
def compute_betti(n,d,p,q):
    #if we can remove the call to .cardinality and compute the cardinality directly, we can remove the
    #complex_basis function entirely
    ker_rank = complex_basis(n,d,p,q).cardinality() - compute_rank(n,d,p,q)
    img_rank = compute_rank(n,d,(p+1),q-1)
    return ker_rank-img_rank

#test code
print compute_betti(2,2,2,0)
print compute_betti(2,2,3,1)
print compute_betti(2,2,2,1)
print compute_betti(2,2,1,1)
print compute_betti(2,2,0,0)
print compute_betti(2,2,0,1)

print compute_betti(2,3,0,0)
print compute_betti(2,3,1,1)
print compute_betti(2,3,2,1)
print compute_betti(2,3,3,1)
print compute_betti(2,3,4,1)
print compute_betti(2,3,5,1)
print compute_betti(2,3,6,1)
print compute_betti(2,3,7,1)

print compute_betti(2,3,1,2)
print compute_betti(2,3,2,2)
print compute_betti(2,3,3,2)
print compute_betti(2,3,4,2)
print compute_betti(2,3,5,2)
print compute_betti(2,3,6,2)
print compute_betti(2,3,7,2)

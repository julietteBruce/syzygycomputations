
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
    ms_list = list(ms)
    ms_list.sort()
    return (tuple(ms_list),f)

#sort the elements into bins by multidegree
def bin_multidegrees(basis):
    bins = dict()
    for elem in basis:
        md = multidegree(elem)
        if md in bins:
            bins[md].append(normalize_elem(elem))
        else:
            bins[md]=[normalize_elem(elem)]
    return bins

#computes the image of a basis element, the returned value is a dict from basis elements in the codomain to a sign
def compute_image(elem,codomain_basis_ids):
    (ms,f) = elem
    def without(lst,i):
        new_list = list(lst)
        del new_list[i];
        return tuple(new_list)
    return {codomain_basis_ids[(without(ms,i),tuple_sum(f,ms[i]))]:(1 if i%2==0 else -1) for i in [0..(len(ms)-1)]}

#some useful functions for displaying elements
def format_monomial(elem):
    var_names = ["x_{" + str(i) + "}" for i in [0..(len(elem)-1)]]
    parts = [v + "^" + "{" + str(i) + "}" for (v,i) in zip(var_names,elem) if i!=0]
    return ''.join(parts)

#pretty print an element of the complex basis
def format_basis(elem):
    (ms,f) = elem
    ms_list = list(ms)
    ms_list.sort()
    return '\wedge '.join(map(format_monomial,ms_list)) + "\otimes " + format_monomial(f)
    

def normalize_md(md):
    md_list = list(md)
    md_list.sort()
    return tuple(md_list)

#compute the rank of the \delta_p map in row q
def compute_rank(n,d,p,q):
    #an element of S_d is given by n integers (
    domain_basis = bin_multidegrees(complex_basis(n,d,p,q))
    codomain_basis = bin_multidegrees(complex_basis(n,d,p-1,q+1))
    #for each multidegree compute the matrix
    counts = dict();
    ranks = dict();
    for md in domain_basis:
        normalized = normalize_md(md)
        if normalized in counts:
            counts[normalized]+=1
        else:
            counts[normalized]=1
        if md!=normalized:
            continue;
        sub_domain_basis = domain_basis[md]
        #safe to skip those where the codomain doesn't contain any of that multidegree
        #this is needed to prevent crashes on certain edge cases
        if md not in codomain_basis:
            ranks[md]=0;
            continue;

        sub_codomain_basis = codomain_basis[md]
        sub_codomain_ids = {sub_codomain_basis[n]:n for n in [0..(len(sub_codomain_basis)-1)]}
        rows = map(lambda x : compute_image(x,sub_codomain_ids),sub_domain_basis)
        def check_image(domainId,codomainId):
            if codomainId in rows[domainId]:
                return rows[domainId][codomainId]
            else:
                return 0
        ranks[md] = matrix(ZZ,len(sub_domain_basis),len(sub_codomain_basis),check_image).rank();
    r=0
    for md in counts:
        r+=counts[md]*ranks[md]
    return r


def compute_betti(n,d,p,q):
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

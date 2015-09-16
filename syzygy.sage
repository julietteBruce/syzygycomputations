
import operator
from collections import defaultdict
import sys
import itertools

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

def normalize_md(md):
    md_list = list(md)
    md_list.sort()
    return tuple(md_list)

def normalize_wedge(elem):
    return tuple(sorted(elem))

def wedge_part(n,d,p):
    sd_basis = compute_basis(n,d);
    return OrderedSubsets(sd_basis,p)

#implements the set of subsets of size k, but instead of returning sets, return ordered tuples
class OrderedSubsets:
    def __init__(self,base_set,k):
        if isinstance(base_set,dict):
            self.base_list = base_set.values();
            self.base_dict = base_set.copy();
        else:
            self.base_list = list(base_set);
            self.base_dict = {elem:n for (n,elem) in enumerate(self.base_list)}
        self.k = k
        self.n = len(self.base_list)
    def __len__(self):
        return int(binomial(self.n,self.k))
    def rank(self,subset):
        indicies = map(lambda x : self.base_dict[x],subset)
        ret = 0
        prev = -1
        for j in range(0,self.k):
            for i in range(prev+1,indicies[j]):
                ret += binomial(self.n-i-1,self.k-j-1)
            prev = indicies[j]
        return ret
    def __as_elements(self,indicies):
        return tuple(map(lambda i : self.base_list[i],indicies))
    def __getitem__(self,idx):
        ret = [-1]*self.k;
        prev = -1
        curr_idx = 0
        for j in range(0,self.k):
            for i in range(prev+1,self.n):
                c = binomial(self.n-i-1,self.k-j-1)
                if curr_idx+c>idx:
                    ret[j]=i
                    prev = i
                    break
                curr_idx+=c
        return self.__as_elements(ret)
    def __iter__(self):
        final_elem = list(range(self.n-self.k,self.n))
        def next_elem(curr):
            for i in range(1,self.k+1):
                if curr[-i]!=self.n-i:
                    break;
            curr[-i] +=1; #advance an entry
            for j in range(1,i):
                curr[-j]=curr[-i]+(i-j)
        def gen():
            curr = list(range(0,self.k))
            yield self.__as_elements(curr)
            while curr!=final_elem:
                next_elem(curr)
                yield self.__as_elements(curr)
        return gen()

#computes the image of a basis element, the returned value is a dict from basis elements in the codomain to a sign
def compute_image(elem,codomain_basis):
    def without(lst,i):
        new_list = list(lst)
        del new_list[i]
        return tuple(new_list)
    return {codomain_basis.rank(without(elem,i)):(1 if i%2==0 else -1) for i in xrange(0,len(elem))}


#for some n,d,q, it returns (as indicies into the wedge basis) the applicable elements of the wedge basis
def slice_multidegree(wedge_multidegrees,n,d,q):
    ret = defaultdict(list)
    f_part = compute_basis(n,q*d)
    for (wedge_md,indicies) in wedge_multidegrees.iteritems():
        for f in f_part:
            md = tuple_sum(wedge_md,f);
            ret[md].extend(indicies);
    return ret

#returns a matrix in row dictionary form, along with with a map from mutlidegrees to slices (lists of rows and columns of the resulting matrix
def construct_matrices(n,d,p,q):
    #first compute the matrix on the basis of wedges
    domain_basis = wedge_part(n,d,p)
    codomain_basis = wedge_part(n,d,p-1)

    rows = map(lambda x : compute_image(x,codomain_basis),domain_basis)

    #instead of constructing slices directly, group the rows and columns by multidegree
    def wedge_multidegree(elem):
        md = (0,)*(n+1)
        for m in elem:
            md = tuple_sum(md,m)
        return md
    def group_basis(basis):
        ret = defaultdict(list)
        for (k,elem) in enumerate(basis):
            md = wedge_multidegree(elem)
            ret[md].append(k)
        return ret
    domain_mds = group_basis(domain_basis)
    codomain_mds = group_basis(codomain_basis)
    return (rows,domain_mds,codomain_mds)

#prints a space separated list
def print_list(outStream,lst):
    for val in lst:
        outStream.write(str(val) + " ");
    outStream.write("\n")

#the format used is as follows
#<nrows> <ndom_mds> <ncodom_mds>
#<rows>
#<domain multidegrees>
#<codomain multidegrees>
#
#where <rows> is represented as follows
#<col>,<val> <col>,<val> ....
#note, no effort will be made to ensure that the col are in increasing order
#
#where each of the multidegree parts are represented by
#<multidegree>: <row/col id> <row/col id> ...

def print_matrices(outStream,rows,domain_mds,codomain_mds):
    print_list(outStream,[len(rows),len(domain_mds),len(codomain_mds)])
    for row in rows:
        for (col,val) in row.iteritems():
            outStream.write('%d,%d ' % (col,val))
        outStream.write('\n')

    for (md,rows) in domain_mds.iteritems():
        outStream.write(str(md) + ": ")
        print_list(outStream,rows)
    for (md,cols) in codomain_mds.iteritems():
        outStream.write(str(md) + ": ")
        print_list(outStream,cols)

def read_list(inStream):
    return inStream.readline().split();

#undoes print_matricies
def read_matrices(inStream):
    nrows,ndom_mds,ncodom_mds = map(int,read_list(inStream))
    rows = [dict()]*nrows
    for i in range(0,nrows):
        entries = read_list(inStream)
        rows[i]=dict(map(lambda s : map(int,s.split(',')),entries))
    domain_mds = dict()
    codomain_mds = dict()
    def read_multidegree_list():
        line = inStream.readline();
        mdstr,idstr = line.split(':')
        md = tuple(map(int,mdstr.strip()[1:-1].split(',')))
        ids = map(int,idstr.split())
        return (md,ids)
    for i in range(0,ndom_mds):
        md,ids = read_multidegree_list();
        domain_mds[md]=ids
    for i in range(0,ncodom_mds):
        md,ids = read_multidegree_list();
        codomain_mds[md]=ids
    return (rows,domain_mds,codomain_mds)

#compute the rank using data in the format returned by construct_matricies
#TODO this function uses far too many arguments, in particular, n,d, and p (but not q) are fixed given
#the other data, but there's no easy way to infer any of them (except maybe n)
def compute_rank(n,d,p,q,rows,domain_mds,codomain_mds):
    domain_slices = slice_multidegree(domain_mds,n,d,q)
    codomain_slices = slice_multidegree(codomain_mds,n,d,q+1)
    counts = defaultdict(int)
    ranks = dict()
    for (md,curr_dom_slice) in domain_slices.iteritems():
        if md not in codomain_slices:
            continue
        curr_codom_slice = codomain_slices[md]
        normalized = normalize_md(md)
        counts[normalized]+=1

        if normalized!=md:
            continue;
        def check_image(domainId,codomainId):
            realDomainId = curr_dom_slice[domainId];
            realCodomainId = curr_codom_slice[codomainId]
            if realCodomainId in rows[realDomainId]:
                return rows[realDomainId][realCodomainId]
            else:
                return 0

        ranks[md] = matrix(ZZ,len(curr_dom_slice),len(curr_codom_slice),check_image).rank()

    r=0
    for md in counts:
        r+=counts[md]*ranks[md]
    return r

#computes a single betti number
def compute_betti(n,d,p,q):
    #if we can remove the call to .cardinality and compute the cardinality directly, we can remove the
    #complex_basis function entirely
    ker_rank = complex_basis(n,d,p,q).cardinality() - compute_rank(n,d,p,q,*construct_matrices(n,d,p,q))
    img_rank = compute_rank(n,d,p+1,q-1,*construct_matrices(n,d,(p+1),q-1))
    return ker_rank-img_rank


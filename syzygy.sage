
#we represent an element of S_d by a n-tuple of integers adding up to d
#compute the basis for S_d
def compute_basis(n,d):
    return map(tuple,IntegerVectors(d,n+1).list())

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

#computes the image of a basis element, the returned value is a list of basis elements in the image
#along with their signs
def compute_image(elem):
    (ms,f) = elem
    #wildly ineffcient way of doing things so that ordering works
    ms_list = list(ms)
    ms_list.sort()
    def without(lst,i):
        new_list = list(lst)
        del new_list[i];
        return new_list
    def tuple_sum(a,b):
        return tuple([x+y for (x,y) in zip(a,b)])
    return [((i%2),[Set(without(ms_list,i)),tuple_sum(f,ms_list[i])]) for i in [0..(len(ms_list)-1)]]

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
    

#compute the rank of the \delta_p map in row q
def compute_rank(n,d,p,q):
    #an element of S_d is given by n integers (
    domain_basis = complex_basis(n,d,p,q)
    codomain_basis = complex_basis(n,d,p-1,q+1)

    rows = map(compute_image,domain_basis)
    #converts a list of (sign,basis_element) to a vector for a row of the matrix
    def to_vector(elem):
        res = [0]*codomain_basis.cardinality()
        for (sign,b) in elem:
            i = codomain_basis.rank(b)
            res[i] = 1 if sign==0 else -1
        return res
    
    return matrix(map(to_vector,rows)).rank()

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

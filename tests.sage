load("syzygy.sage")

print list(OrderedSubsets([0,1,2,3,4,5],2))
test = list()
test2 = OrderedSubsets([0,1,2,3,4,5],2)
for i in range(0,len(test2)):
    test.append(test2[i])
print test

for (i,elem) in enumerate(test2):
    print i," ",test2.rank(elem)


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


with open("test.txt",'w') as f:
    print_matrices(f,*lazy_construct_matrices(2,3,5,1))

with open("test.txt",'r') as f:
    print compute_rank(2,3,5,1,*read_matrices(f))

print compute_rank(2,3,5,1,*construct_matrices(2,3,5,1))

#WARNING this test creates a ~40Mb file
with open("big_test.txt",'w') as f:
    print_matrices(f,*lazy_construct_matrices(2,5,10,1))

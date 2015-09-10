load("syzygy.sage")

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
    print_matrices(f,*construct_matrices(2,3,5,1))

with open("test.txt",'r') as f:
    print compute_rank(*read_matrices(f))

print compute_rank(*construct_matrices(2,3,5,1))

import sage.coding.linear_code
import sage.combinat.permutation


def gg(C):
    if C < 2 :
        return 0
    else :
        return 1

print gg(0)        


def minimum_distance(C):
    G = C.generator_matrix()
    n = C.length()
    a = range(1,n + 1)
    print(a)
    shuffle(a)
    print(a)
    a = Permutation(a)
    print(G) 
    print(" ")
    G.permute_columns(a)
    print(G)
    G = G.echelon_form()
    print(" ")
    print(G)    
    


C = codes.LinearCode(random_matrix(GF(2),3,10))
minimum_distance(C)
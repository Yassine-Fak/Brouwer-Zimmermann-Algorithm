import sage.coding.linear_code
import sage.combinat.permutation


def minimum_distance_brouwn(C):
    G = C.generator_matrix()
    n = C.length()
    k = C.dimension()
    a = range(1,n + 1)
    print(a)
    shuffle(a)
    print(a)
    a = Permutation(a)
    print(" ")
    print(G) 
    print(" ")
    G.permute_columns(a)
    print(G)
    G = G.echelon_form()
    print(" ")
    print(G)

    l = 0
    for i in range(n//k): 

        L = G.matrix_from_columns(range(i*k,(i+1)*k))
        if L.rank() == k:
            l = l + 1 
        else :
            break 
           

    print(l)

         
    


C = codes.LinearCode(random_matrix(GF(5),3,20))
minimum_distance_brouwn(C)


def information_set(M):
    n = M.nrows()
    k = M.ncols()
    # Je vais sauvegarder le resultat de chaque information_set dans une ligne de la matrice resultat
    res = matrix(n//k,k)
    M2 = copy(M)

    i = 0
    while M2.rank() == k:
        res[i] = M2.pivot()





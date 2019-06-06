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
            # il faut verifier que les Ai sont dijointe 

    print(l)
    
    # Mnt je vais trouver le maximum d'informations set que je peux avoir :
    # Et si la matrice G a une colonne aui est dupplique plusierus fois ? 
    # Comment trouver le corps sur lequelle les elments sont ecrit ?


         
    


C = codes.LinearCode(random_matrix(GF(5),3,20))
minimum_distance_brouwn(C)
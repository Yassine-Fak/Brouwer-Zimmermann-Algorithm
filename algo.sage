import sage.coding.linear_code
import sage.combinat.permutation

def systematique_form(G):
    M = copy(G)
    a = M.pivots()
    b = [ i+1 for i in a ]
    c = b + [ i for i in range(1, M.ncols()) if i not in b]
    c = Permutation(c)
    M.permute_columns(c)
    return M.echelon_form()

# La fonction suivante rernvoie les informations set, meme ceux aui ne sont pas dijoints
def greedy_information_set(M):
    k = M.nrows()
    n = M.ncols()
    Com = Combinations(range(n),k)
    res = []
    l = Com.list()
    num_info_set = 0
    for i in range(len(l)):
        if M.matrix_from_columns(l[i]).rank() == k:
            num_info_set = num_info_set + 1
            res = res + [l[i]]
    print(" The number of infomation Sets is : " , num_info_set )  
    print res  

def infomation_set(G):
    M = copy(G)
    k = G.nrows()
    n = G.ncols()
    res = matrix(k,n)
    num_info_set = 0
    while M.rank() == k :
        num_info_set = num_info_set + 1
        a = M.pivots()
        a = list(a)
        L = M.matrix_from_columns(a)
        for i in range










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
#minimum_distance_brouwn(C)

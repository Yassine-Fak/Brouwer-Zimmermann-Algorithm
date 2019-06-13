

def systematique_form(G):
    M = copy(G)
    a = M.pivots()
    b = [ i+1 for i in a ]
    c = b + [ i for i in range(1, M.ncols()) if i not in b]
    c = Permutation(c)
    M.permute_columns(c)
    return M.echelon_form()


def infomation_set(G):

    M = copy(G)
    k = G.nrows()
    num_info_set = 0
    L = M.matrix_from_columns(range(num_info_set*k,G.ncols()))

    while L.rank() == k:

        a = L.pivots()
        b = range(0,num_info_set*k) + [num_info_set*k + i for i in a]
        c = b
        for i in range(G.ncols()):
          if i in b:
            pass
          else:
            c = c + [i]
        d = [i + 1 for i in c]
        d = Permutation(d)
        M.permute_columns(d)
        num_info_set = num_info_set + 1
        L = M.matrix_from_columns(range(num_info_set*k,G.ncols()))

    return (M,num_info_set)

def infomation_set_brouwer(G, maxiter = 150):

    k = G.nrows()
    n = G.ncols()
    M, num_info_set = infomation_set(G)
    q = n//k
    i = 0

    while i < maxiter and num_info_set != q :

        anc_num_info_set = num_info_set
        anc_M = copy(M)
        M = copy(G)
        p = Permutations(n).random_element()
        M.permute_columns(p)
        M, num_info_set = infomation_set(M)

        if anc_num_info_set > num_info_set :
          M = copy(anc_M)
          num_info_set = anc_num_info_set
        
        i = i + 1

    return (M,num_info_set)



def infomation_set_brouwer_zimmer(G, maxiter = 100):

    k = G.nrows()
    n = G.ncols()
    M, num_info_set = infomation_set(G)
    q = n//k
    r = n%k
    
    for i in range(maxiter):

      anc_M = copy(M)
      anc_num_info_set = num_info_set
      M = copy(G)
      p = Permutations(n).random_element()
      M.permute_columns(p)
      M, num_info_set = infomation_set(M)

      if num_info_set < anc_num_info_set :
        M = copy(anc_M)
        num_info_set = anc_num_info_set
        continue

      if num_info_set > anc_num_info_set :
        continue
      
      if num_info_set == q and M.matrix_from_columns(range(num_info_set*k,n)).rank() == r :
        break 
      
      if num_info_set == anc_num_info_set :
        if M.matrix_from_columns(range(num_info_set*k,n)).rank() >= anc_M.matrix_from_columns(range(num_info_set*k,n)).rank() :
          continue
        
        else : 
          M = copy(anc_M)
          num_info_set = anc_num_info_set
        
    return (M,num_info_set)
    

def list_of_system_gen_mat(M,m):

    k = M.nrows()
    L = []

    for i in range(m):
      A = M.matrix_from_columns(range(i*k , i*k + k))
      L = L + [A.inverse()*M]
    
    return L


def minimum_distance_brouwer(C):

    G1 = C.generator_matrix()
    k = G1.nrows()
    n = G1.ncols()
    F = C.base_field()
    V = VectorSpace(F,k)
    G2, num_info_set = infomation_set_brouwer(G1)
    L = list_of_system_gen_mat(G2,num_info_set)
    ub = n - k + 1
    lb = num_info_set
    w = 1

    print("The number of disjoint information set is : {} ".format(num_info_set))

    while w <= k and lb < ub :

      elmts_of_vec_spa = [i for i in V if i.hamming_weight() == w]

      for j in range(0,num_info_set) :
        for x in elmts_of_vec_spa :
          ub = min(ub, (x*L[j]).hamming_weight())
          if ub <= lb :
            return ub
        lb = lb + 1
      w = w + 1
    return ub


C = codes.random_linear_code(GF(2),46,3)
#C = codes.random_linear_code(GF(2),61,5)
G = C.generator_matrix()


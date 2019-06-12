

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


def infomation_set_brouwn(G, maxiter = 50):

    M = copy(G)
    k = G.nrows()
    n = G.ncols()
    q = n//k
    M, num_info_set = infomation_set(M)
    i = 0

    while i < maxiter and num_info_set == q :

        anc_num_info_set = num_info_set
        p = Permutations(n).random_element()
        M.permute_columns(p)
        M, num_info_set = infomation_set(M)

        if anc_num_info_set > num_info_set : 
          M.permute_columns(p.inverse())
          num_info_set = anc_num_info_set

        i = i + 1
    
    return (M,num_info_set)








def minimum_distance_brouwn(C):

    G1 = C.generator_matrix()
    k = G1.nrows()
    n = G1.ncols()
    F = C.base_field()
    V = VectorSpace(F,k)
    G2, num_info_set = infomation_set(G1) 
    ub = n - k + 1
    lb = num_info_set 
    L = []
    w = 1
    
    print("The number of disjoint information set is : {} ".format(num_info_set))
    
    for i in range(num_info_set):

      A = G2.matrix_from_columns(range(i*k , i*k + k))
      L = L + [A.inverse()*G2]
      #print L[i]
      #print " " 
    #print L

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
          

C = codes.random_linear_code(GF(2),15,3)
G = C.generator_matrix()
Permutations(4).random_element()

# http://doc.sagemath.org/html/en/reference/modules/sage/modules/vector_mod2_dense.html
# http://doc.sagemath.org/html/en/reference/modules/sage/modules/free_module_element.html



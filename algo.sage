from sage.combinat.gray_codes import combinations


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
    L = M.matrix_from_columns(range(num_info_set*k,G.ncols())).echelon_form()

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
      num_info_set += 1
      L = M.matrix_from_columns(range(num_info_set*k,G.ncols())).echelon_form()

    return (M,num_info_set)
   

def infomation_set_brouwer(G, maxiter = 200):

    k = G.nrows()
    n = G.ncols()
    M, num_info_set = infomation_set(G)
    q = n//k
    i = 0
    
    while i < maxiter and num_info_set != q :
      
      M_inter = copy(G)
      p = Permutations(n).random_element()
      M_inter.permute_columns(p)
      M_inter, num_info_set_inter = infomation_set(M_inter)

      if num_info_set_inter > num_info_set :
        M = copy(M_inter)
        num_info_set = num_info_set_inter
      
      i += 1

    return (M,num_info_set)


def infomation_set_brouwer_zimmer(G, maxiter = 200):

    k = G.nrows()
    n = G.ncols()
    M, num_info_set = infomation_set(G)
    R = M.matrix_from_columns(range(num_info_set*k,n))
    q = n//k
    r = n%k
    i = 0
    
    while i < maxiter and num_info_set != q and R.rank() != r:

      M_inter = copy(G)
      p = Permutations(n).random_element()
      M_inter.permute_columns(p)
      M_inter, num_info_set_inter = infomation_set(M_inter)
      R_inter = M_inter.matrix_from_columns(range(num_info_set_inter*k,n))

      if num_info_set_inter > num_info_set :
        M = copy(M_inter)
        R = copy(R_inter)
        num_info_set = num_info_set_inter
        i += 1
        continue

      if num_info_set_inter < num_info_set :
        i += 1
        continue
      
      if num_info_set == num_info_set_inter :

        if R.rank() >= R_inter.rank() :
          i += 1
          continue

        else :
          M = copy(M_inter)
          R = copy(R_inter)
          i += 1
          continue
      
      if num_info_set_inter == q and R_inter.rank() == r :
        M = copy(M_inter)
        num_info_set = num_info_set_inter
        break

    return (M,num_info_set)


def list_of_system_gen_mat(M,m):

    k = M.nrows()
    L = []

    for i in range(m):
      A = M.matrix_from_columns(range(i*k , i*k + k))
      L = L + [A.inverse()*M]
    
    return L


def minimum_distance_brouwer_ancien(C):

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
        lb += 1
      w += 1
    return ub


C = codes.random_linear_code(GF(2),46,3)
#C = codes.random_linear_code(GF(2),61,5)
G = C.generator_matrix()

# http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/combination.html
# http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/gray_codes.html


def minimum_distance_brouwer(C):

    G1 = C.generator_matrix()
    n,k = C.length(), C.dimension()
    F = C.base_field()
    g = F.multiplicative_generator()
    q = F.cardinality() 

    G2, num_info_set = infomation_set_brouwer(G1)
    L = list_of_system_gen_mat(G2,num_info_set)
    ub = n - k + 1 
    lb = num_info_set
    w = 1

    print("The number of disjoint information set is : {} ".format(num_info_set))

    while w <= k and lb < ub :

      for m in range(0,num_info_set) : # pour calculer G22 = L[m]

        for x in [g**t for t in range(q-1)]: # pour enumerer les element du corps 
          X = zero_matrix(1,k)
          for e in range(w):
            X[0,e] = x
          X = X[0]
          X = list(X)
          ub = min(ub, (X*L[m]).hamming_weight())
          if ub <= lb :
            return ub
          for i,j in combinations(k,w): # il faut discuter le cas lorsque on met x de partout
            X[i] = 0
            X[j] = x
            ub = min(ub, (X*L[m]).hamming_weight())
            if ub <= lb :
              return ub
        # Ici je traite le cas lorsqu'on a uniquement un seul element dans la liste X 


        lb += 1     

      w += 1

    # le prb est aue ca genere uniauement tous les mots de pd = 1, et le pd = 2, 3 4 ???????

    return ub




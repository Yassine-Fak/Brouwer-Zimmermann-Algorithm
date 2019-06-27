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


def infomation_set_brouwer(G, maxiter = 3):

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


def infomation_set_brouwer_zimmer(G, maxiter = 3):

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


def list_of_system_gen_mat(M,m,k):

    L = []
    for i in xrange(m):
      A = M.matrix_from_columns(xrange(i*k , i*k + k))
      L = L + [A.inverse()*M]

    return L

def minimum_distance_brouwer(C):

    G1 = C.generator_matrix()
    n,k = C.length(), C.dimension()
    F = C.base_field()
    g = F.multiplicative_generator()
    q = F.cardinality()
    G2, num_info_set = infomation_set_brouwer(G1)
    L = list_of_system_gen_mat(G2,num_info_set,k)
    M = [g^i for i in xrange(q-1)]
    ub = n - k + 1
    lb = num_info_set
    w = 1
    Z = IntegerRing()
    print("The number of disjoint information set is : {} ".format(num_info_set))
    
    if F == GF(2) :
      while w <= k and lb < ub :
        for m in xrange(0,num_info_set) : # pour calculer G22 = L[m]
          A = L[m].row(0)
          for i in xrange(1,w):
            A += L[m].row(i)
          ub = min(ub, A.hamming_weight())
          if ub <= lb :
            return ub

          for i,j in combinations(k,w):
            A += L[m].row(i) + L[m].row(j)
            ub = min(ub, A.hamming_weight())
            if ub <= lb :
              return ub
          lb += 1
        w += 1
      return ub

    while w <= k and lb < ub :
      for m in xrange(num_info_set) : # pour calculer G22 = L[m]
        A = L[m].row(0)
        for i in xrange(1,w):
          A += L[m].row(i)
        a = [0]*w
        for v in xrange((q-1)^w):
          a_anc = copy(a)
          a = Z(v).digits(q-1,padto=w) 
          X = [M[a[w-1-i]] for i in xrange(w)] + [F.zero()]*(k-w)
          a_supp = []
          for i in xrange(w):
            if a[i] != a_anc[i]:
              a_supp += [i]
          
          for i in a_supp :
            A += (M[a[i]] - M[a_anc[i]])*L[m].row(w-1-i)

          ub = min(ub, A.hamming_weight())
          if ub <= lb :
            return ub
          
          A_int = copy(A) 
          for i,j in combinations(k,w):
            X[j] = X[i]
            X[i] = F.zero()
            A_int += X[j]*(L[m].row(j) - L[m].row(i))
            ub = min(ub, A_int.hamming_weight())
            if ub <= lb :
              return ub
        lb += 1
      w += 1
    return ub




C = codes.random_linear_code(GF(7),40,5) 
C = codes.random_linear_code(GF(17),15,4)

C = codes.random_linear_code(GF(13),30,9) 
C = codes.random_linear_code(GF(5),50,11) # le meilleur en tp est C.min < nv < ancien err car nv != C. et opti = C.

C = codes.random_linear_code(GF(7),50,7)
C = codes.random_linear_code(GF(11),50,5)
C = codes.random_linear_code(GF(17),35,6) #met du temps
C = codes.random_linear_code(GF(5),55,10)
C = codes.random_linear_code(GF(5),55,9) 

C = codes.random_linear_code(GF(5),44,5) 

C = codes.random_linear_code(GF(2),100,11)  
C = codes.random_linear_code(GF(2),100,25)  

C = codes.random_linear_code(GF(3),100,11) 
C = codes.random_linear_code(GF(23),35,6) #met du temps 


# http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/permutation.html?highlight=permutation#module-sage.combinat.permutation
# https://www.diveinto.org/python3/advanced-iterators.html
# https://git.sagemath.org/sage.git/tree/src/sage/combinat/gray_codes.py?id=3b92a85394e8fc6f7572c58d0a6ce8b1bebacce2
# Volume 4 Fascicle 2, Generating All Tuples and Permutations (2005), v+128pp. ISBN 0-201-85393-0

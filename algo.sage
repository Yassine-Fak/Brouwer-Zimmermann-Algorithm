from sage.combinat.gray_codes import combinations
import time 


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
    n = G.ncols()
    num_info_set = 0
    L = M.matrix_from_columns(range(num_info_set*k,n)).echelon_form()

    while L.rank() == k:

      a = L.pivots()
      b = range(0,num_info_set*k) + [num_info_set*k + i for i in a]
      c = copy(b)
      for i in range(n):
        if i in b:
          pass
        else:
          c = c + [i]
      d = [i + 1 for i in c]
      d = Permutation(d)
      M.permute_columns(d)
      num_info_set += 1
      L = M.matrix_from_columns(range(num_info_set*k,n)).echelon_form()

    return (M,num_info_set)


def infomation_set_brouwer(G, maxiter):

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


def infomation_set_brouwer_zimmer(G, maxiter):

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


def minimum_distance_brouwer(C,maxiter=3):

    G1 = C.generator_matrix()
    n,k = C.length(), C.dimension()
    F = C.base_field()
    G2, num_info_set = infomation_set_brouwer(G1,maxiter)
    L = list_of_system_gen_mat(G2,num_info_set,k)
    ub = n - k + 1
    lb = num_info_set
    w = 1
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
    
    q = F.cardinality()
    g = F.multiplicative_generator()
    M = [g^i for i in xrange(q-1)]
    Z = IntegerRing()

    while w <= k and lb < ub :
      X = [F.zero()]*k
      for m in xrange(num_info_set) : # pour calculer G22 = L[m]
        A = L[m].row(0)
        for i in xrange(1,w):
          A += L[m].row(i)
        a = [0]*w
        for v in xrange((q-1)^w):
          a_anc = copy(a)
          a = Z(v).digits(q-1,padto=w)
          for i in xrange(w):
            X[i] = M[a[w-1-i]]  
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


def test_rapide(): # moins de deux min.
  # (GF(),long,dim)
  L = [(11,44,6),(2,77,15),(2,100,11),(11,33,5),(3,100,11),(3,45,7),(9,50,8),(5,44,5),(25,28,5),(2,100,25),(7,40,5),(17,15,4),(7,50,7),(11,50,5),(5,55,10),(5,55,9)]
  print ("-------------------")
  for x in L :
    C = codes.random_linear_code(GF(x[0]),x[1],x[2])
    print("For {} we have : ".format(C))
    print (" ")
    print("C.minimum_distance() : ")
    a = %time C.minimum_distance()
    print a
    print (" ")
    print("minimum_distance_brouwer(C) :")
    b = %time minimum_distance_brouwer(C)
    print b
    print (" ")
    print("minimum_distance_zimmermann(C) :")
    c = %time minimum_distance_zimmermann(C)
    print c
    print ("-------------------")

def test_lent():
  # (GF(),long,dim)
  L = [(7^2,35,6)]
  print ("-------------------")
  for x in L :
    C = codes.random_linear_code(GF(x[0]),x[1],x[2])
    print("For {} we have : ".format(C))
    print (" ")
    print("C.minimum_distance() : ")
    a = %time C.minimum_distance()
    print a
    print (" ")
    print("minimum_distance_brouwer(C) :")
    b = %time minimum_distance_brouwer(C)
    print b
    print ("-------------------")


def minimum_distance_zimmermann(C,maxiter=3):

    G1 = C.generator_matrix()
    n,k = C.length(), C.dimension()
    F = C.base_field()
    G2, num_info_set = infomation_set_brouwer_zimmer(G1,maxiter)
    R = G2.matrix_from_columns(range(num_info_set*k,n))
    q = F.cardinality()
    g = F.multiplicative_generator()
    M = [g^i for i in xrange(q-1)]
    Z = IntegerRing()

    if R.rank() == 0 :
      print "Brouwer1 ! "
      return minimum_distance_brouwer(C,maxiter)

    if n - k*num_info_set >= k :
      G2, num_info_set = infomation_set_brouwer_zimmer(G1,10)
      R = G2.matrix_from_columns(range(num_info_set*k,n))
      if R.rank() == 0 :
        print "Brouwer2 ! "
        return minimum_distance_brouwer(C,10)
      print("The number of disjoint information set is : {} ".format(num_info_set))
    else :
      print("The number of disjoint information set is : {} ".format(num_info_set))

    R = G2.matrix_from_columns(range(num_info_set*k,n))
    R_pivot = R.pivots()
    r = len(R_pivot)
    b = range(0,num_info_set*k) + [num_info_set*k + i for i in R_pivot]
    c = copy(b)
    for i in range(n):
      if i in b:
        pass
      else:
        c = c + [i]
    d = [i + 1 for i in c]
    d = Permutation(d)
    G2.permute_columns(d)

    A = G2.matrix_from_columns(range(num_info_set*k,num_info_set*k+r))
    c = 0
    I = range(num_info_set*k,num_info_set*k+r)
    while A.rank() != k:
      if A.rank() < G2.matrix_from_columns(I + [c]).rank():
        A = G2.matrix_from_columns(I + [c])
        I += [c]
      c += 1
      if c == k :
        break
    
    L = list_of_system_gen_mat(G2,num_info_set,k) + [A.inverse()*G2]
    num_info_set += 1
    lb = 1
    ub = n - k + 1
    w = 1
    rel_ran = [k]*(num_info_set - 1) + [r] # The sequence of relative ranks

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
        for x in rel_ran:
          lb += max(0,(w+1)-(k-x))
        w += 1
      return ub

    while w <= k and lb < ub :
      X = [F.zero()]*k
      for m in xrange(num_info_set) : # pour calculer G22 = L[m]
        A = L[m].row(0)
        for i in xrange(1,w):
          A += L[m].row(i)
        a = [0]*w
        for v in xrange((q-1)^w):
          a_anc = copy(a)
          a = Z(v).digits(q-1,padto=w)
          for i in xrange(w):
            X[i] = M[a[w-1-i]]  
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
      for x in rel_ran:
        lb += max(0,(w+1)-(k-x))
      w += 1
    return ub


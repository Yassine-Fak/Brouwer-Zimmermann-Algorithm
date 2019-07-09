from sage.combinat.gray_codes import combinations
import time 


def systematique_form(G):
  M = copy(G)
  a = M.pivots()
  b = list(a)
  c = copy(b)
  for i in xrange(M.ncols()):
    if i in b :
      pass
    else :
      c += [i]
  d = [i+1 for i in c]
  d = Permutation(d)
  M.permute_columns(d)
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


def infomation_set_zimmermann(G, maxiter):

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


def minimum_distance_one_information_set(C):
    n,k = C.length(), C.dimension()
    F = C.base_field()
    q = F.cardinality()
    ub = n - k + 1
    w = 1
    G = systematique_form(C.generator_matrix())
    lb = 1
    if F == GF(2) :
      while w <= k and lb < ub : 
        A = G.row(0)
        for i in xrange(1,w):
          A += G.row(i)
        ub = min(ub, A.hamming_weight())
        if ub <= lb :
          return ub
        for i,j in combinations(k,w):
          A += G.row(i) + G.row(j)
          ub = min(ub, A.hamming_weight())
          if ub <= lb :
            return ub
        lb += 1 
        w += 1
      return ub

    g = F.multiplicative_generator()
    M = [g^i for i in xrange(q-1)]
    Z = IntegerRing()

    while w <= k and lb < ub :
      X = [F.zero()]*k
      A = G.row(0)
      for i in xrange(1,w):
        A += G.row(i)
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
          A += (M[a[i]] - M[a_anc[i]])*G.row(w-1-i)

        ub = min(ub, A.hamming_weight())
        if ub <= lb :
          return ub
        
        A_int = copy(A)
        for i,j in combinations(k,w):
          X[j] = X[i]
          X[i] = F.zero()
          A_int += X[j]*(G.row(j) - G.row(i))
          ub = min(ub, A_int.hamming_weight())
          if ub <= lb :
            return ub
      lb += 1
      w += 1
    return ub


def minimum_distance_brouwer(C, nb_words=900000, maxiter=5):

    n,k = C.length(), C.dimension()
    F = C.base_field()
    q = F.cardinality()
    w = 1

    if q^k <= nb_words :
      print("Here we use one information set ! ")
      return minimum_distance_one_information_set(C)

    G2, num_info_set = infomation_set_distance(C.generator_matrix(),maxiter, method = "brouwer")
    L = []
    for i in xrange(num_info_set):
      A = G2.matrix_from_columns(xrange(i*k , i*k + k))
      L += [A.inverse()*G2]
    lb = num_info_set
    ub = n - k + 1
    print("The number of disjoint information set is : {} ".format(num_info_set))
    
    if F == GF(2) :
      while w <= k and lb < ub :
        for m in xrange(0,num_info_set) : 
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
      for m in xrange(num_info_set) : 
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


def minimum_distance_zimmermann(C,maxiter=5):

    G1 = C.generator_matrix()
    n,k = C.length(), C.dimension()
    F = C.base_field()
    G2, num_info_set = infomation_set_zimmermann(G1,maxiter)
    R = G2.matrix_from_columns(range(num_info_set*k,n))
    R_pivot = R.pivots()
    r = len(R_pivot)

    if r != 0:
      b = range(num_info_set*k) + [num_info_set*k + i for i in R_pivot]
      c = copy(b)
      for i in range(n):
        if i in b:
          pass
        else:
          c += [i]
      d = [i + 1 for i in c]
      d = Permutation(d)
      G2.permute_columns(d)
      I = range(num_info_set*k,num_info_set*k+r)
      A = G2.matrix_from_columns(I)
      c = 0
      while A.rank() != k:
        if A.rank() < G2.matrix_from_columns(I + [c]).rank():
          A = G2.matrix_from_columns(I + [c])
          I += [c]
        c += 1
        if c == k :
          break    
      L = list_of_system_gen_mat(G2,num_info_set,k) + [A.inverse()*G2]
      R = [k]*num_info_set + [r]
      num_info_set += 1

    else :
      L = list_of_system_gen_mat(G2,num_info_set,k)
      R = [k]*num_info_set

    print("The number of disjoint information set is : {} ".format(num_info_set))
    lb = 1
    ub = n - k + 1
    w = 1

    if F == GF(2) :
      while w <= k and lb < ub :
        lb_int = copy(lb)
        lb = 0
        for m in xrange(num_info_set) : # pour calculer G22 = L[m]
          A = L[m].row(0)
          for i in xrange(1,w):
            A += L[m].row(i)
          ub = min(ub, A.hamming_weight())
          if ub <= lb_int :
            return ub

          for i,j in combinations(k,w):
            A += L[m].row(i) + L[m].row(j)
            ub = min(ub, A.hamming_weight())
            if ub <= lb_int :
              return ub
          lb += max(0,w+1-k+R[m])
        w += 1
      return ub
    
    q = F.cardinality()
    g = F.multiplicative_generator()
    M = [g^i for i in xrange(q-1)]
    Z = IntegerRing()
       
    while w <= k and lb < ub :
      X = [F.zero()]*k
      lb_int = copy(lb)
      lb = 0
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
          if ub <= lb_int :
            return ub
          
          A_int = copy(A) 
          for i,j in combinations(k,w):
            X[j] = X[i]
            X[i] = F.zero()
            A_int += X[j]*(L[m].row(j) - L[m].row(i))
            ub = min(ub, A_int.hamming_weight())
            if ub <= lb_int :
              return ub
        lb += max(0,w+1-k+R[m])
      w += 1
    return ub


def test_rapide(): # moins de deux min.
  # (GF(),long,dim)
  L = [(3,60,15),(3,64,16),(3,68,17),(5,44,11),(5,48,12),(7,36,9),(7,40,10),(8,36,9),(8,40,10),(9,32,8),(9,36,9),(3,45,7),(11,33,5),(3,100,11),(11,44,6),(2,77,15),(2,100,11),(9,50,8),(5,44,5),(25,28,5),(2,100,25),(7,40,5),(17,15,4),(7,50,7),(11,50,5),(5,55,10),(5,55,9)]
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
    e = %time minimum_distance_zimmermann(C)
    print e
    print ("-------------------")


def test_rapide_gf2(): # moins de deux min.
  # (GF(),long,dim)
  L = [(2,44,6),(2,77,15),(2,100,11),(2,33,5),(2,100,11),(2,45,7),(2,50,8),(2,44,5),(2,28,5),(2,100,25),(2,40,5),(2,15,4),(2,50,7),(2,50,5),(2,55,10),(2,55,9)]
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
    d = %time minimum_distance_zimmermann(C)
    print d
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
    print (" ")
    print("minimum_distance_zimmermann(C) :")
    d = %time minimum_distance_zimmermann(C)
    print d
    print ("-------------------")


def test_lent_gf2():
  # (GF(),long,dim)
  L = [(2,128,32),(2,96,32)]
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
    d = %time minimum_distance_zimmermann(C)
    print d
    print ("-------------------")


#  l'idee est que j'tulise la fonction que j'ai deja pour zimmerman 
def infomation_set_distance(G,maxiter, method = "zimmermann"):
  
    k = G.nrows()
    n = G.ncols()
    M, num_info_set = infomation_set(G)
    R = M.matrix_from_columns(range(num_info_set*k,n))
    q = n//k
    r = n%k
    i = 1
    c = 1

    if method == "brouwer" :
      while i <= maxiter and num_info_set != q :
        M_inter = copy(G)
        p = Permutations(n).random_element()
        M_inter.permute_columns(p)
        M_inter, num_info_set_inter = infomation_set(M_inter)
        if num_info_set_inter > num_info_set :
          M = copy(M_inter)
          num_info_set = num_info_set_inter
        i += 1
      return (M,num_info_set)

    # method = "zimmermann"  
    while i <= maxiter and num_info_set != q and R.rank() != r:

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
        return (M,[k]*num_info_set_inter + [r])


    # le cas ideal
    if num_info_set == q and R.rank() == r:
      list_of_ranks = [k]*num_info_set + [r]
      return (M,list_of_ranks)
    
    if R_inter.rank() == 0:
      list_of_ranks = [k]*num_info_set
      return (M,list_of_ranks)

    while i <= maxiter:
      M_inter = copy(G)
      p = Permutations(n).random_element()
      M_inter.permute_columns(p)
      M_inter, num_info_set_inter = infomation_set(M_inter)
      # je cherche a maximiser le nbre d'IS disjoint  
      if num_info_set_inter < num_info_set: # si cette condition est tjr satifaite qui ce qu'je fait 
        i += 1
        pass
      R_inter = M_inter.matrix_from_columns(range(num_info_set_inter*k,n))
      list_of_ranks = [k]*num_info_set_inter
      
      # le cas ideal
      if num_info_set_inter == q and R_inter.rank() == r:
        M = copy(M_inter)
        list_of_ranks += [r]
        return (M,list_of_ranks)

      if R_inter.rank() == 0:
        M = copy(M_inter)
        return (M,list_of_ranks)
      
      while R_inter.rank() != 0:
        R_inter_pivots = R_inter.pivots()
        len_of_R_inter_pivots = len(R_inter_pivots)
        b = range(n - R_inter.ncols()) + [n - R_inter.ncols() + i for i in R_inter_pivots]
        c = copy(b)
        for i in range(n):
          if i in b:
            pass
          else:
            c += [i]
        d = [i + 1 for i in c]
        d = Permutation(d) 
        M_inter.permute_columns(d)
        list_of_ranks += [len_of_R_inter_pivots]
        R_inter = M_inter.matrix_from_columns(range(sum(list_of_ranks),n))
        num_info_set_inter += 1
      M = copy(M_inter)
      i += 1
    return (M,list_of_ranks)

# je vais terminer la fonction des information set pour zimmerman et crer une autre pour completer les colonnes des is qui ne sont pas dijoint et apres faire les prediction  
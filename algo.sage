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
      for i in xrange(n):
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


def percentage_disjoint_infor_set(G,maxiter=1000000):
    k = G.nrows()
    n = G.ncols()
    M, num_info_set = infomation_set(G)
    q = n//k
    i = 1
    nb_permut = 1    

    while i <= maxiter and num_info_set != q :
      nb_permut += 1
      M_inter = copy(G)
      p = Permutations(n).random_element()
      M_inter.permute_columns(p)
      M_inter, num_info_set_inter = infomation_set(M_inter)
      if num_info_set_inter > num_info_set :
        M = copy(M_inter)
        num_info_set = num_info_set_inter
      i += 1
    return (num_info_set,nb_permut)


def best_permutation(C,nb_iter=1000):

    G = C.generator_matrix()
    list_of_nb_of_permu = []
    L = []
    f = open('result_permutation', 'w')
    for i in xrange(nb_iter):
      num_info_set,nb_permut = percentage_disjoint_infor_set(G,maxiter=1000000)
      list_of_nb_of_permu += [nb_permut]
      f.write("iteration : {} , num_info_set : {} , nb_permut : {} \n ".format(i,num_info_set,nb_permut))
    maximum_nb_of_permu = max(list_of_nb_of_permu)
    X = range(maximum_nb_of_permu+1)

    f.write("  \n")
    f.write("(Le nb de permut pour atteindre le nb maxi d'is, le nombre de fois que ce nombre a été affiché) \n")
    f.write("  \n")

    for i in X:
      counter = 0
      for x in list_of_nb_of_permu:
        if x == i:
          counter += 1
      L += [(i,counter)]
    f.write("{}".format(L))
    Y = []
    for i in X:
      if i == 0:
        Y += [0]
      else:
        Y += [ Y[i-1] + L[i][1] ]
        
    for i in X:
      Y[i] = float(Y[i]/nb_iter) 
    f.write("  \n")
    f.write("X = {}".format(X))
    f.write("  \n")
    f.write("Y = {}".format(Y))
    f.close()
    # list_plot(Y,color='red')
    G = list_plot(Y,plotjoined=True,color='red')
    G.save("sage.png")
    # G.show()
    x = 1
    y = 1
    for i in xrange(len(Y)):
      if Y[i] >= 0.8 and Y[i] < 0.9:
        x = X[i]
        break 

    for j in xrange(len(Y)):
      if Y[j] >= 0.9 and Y[j] < 1:
        y = X[j]
        break 
    return (maximum_nb_of_permu,y,x)

def quotient_in_function_of_the_length():
    f = open('quotient_in_function_of_the_length', 'w')
    # X = range(70,130,9)
    X = [2,3,5,7,9,11,13]
    Y = []
    f.write("  \n")
    f.write("X = {}".format(X))
    f.write("  \n")
    for n in X:
      # C = codes.random_linear_code(GF(2), n, 32)
      C = codes.random_linear_code(GF(n), 40, 10)
      Start_Time_min_dist = time.time()
      C.minimum_distance()
      Execution_Time_min_dist = time.time() - Start_Time_min_dist
      Start_Time_Zimmer = time.time()
      minimum_distance_zimmermann(C)
      Execution_Time_zimm = time.time() - Start_Time_Zimmer
      Y += [float(Execution_Time_min_dist/Execution_Time_zimm)]
    f.write("Y = {}".format(Y))
    f.close()
    G = list_plot(Y,plotjoined=True,color='red')
    G.save("quotient_in_function_of_the_length.png")


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
    # Dans le cas de zimmermann la fonction renvoie la liste des rangs

    # Le cas ideal 
    if num_info_set == q and R.rank() == r:
      if r != 0:
        return (M,[k]*num_info_set + [r])
      else : 
        return (M,[k]*num_info_set)

    # Le cas patho.
    if R.rank() == 0:
      return (M,[k]*num_info_set)
    
    # Je chercher a maximiser le nombre d'IS et le range de R
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
        if r != 0:
          return (M,[k]*num_info_set_inter + [r])
        else :
          return (M,[k]*num_info_set_inter)

    list_of_ranks = [k]*num_info_set

    if R.rank() == 0:
      return (M,list_of_ranks)
    
    # Ici je suppose qu'on a fait toutes les iterations sans atteindre la cas ideal 
    # Mnt on va essayer d'augmenter le nombre d'information set
    while R.rank() != 0 :
      R_pivots = R.pivots()
      len_of_R_pivots = len(R_pivots)
      b = range(n - R.ncols()) + [n - R.ncols() + i for i in R_pivots]
      c = copy(b)
      for i in range(n):
        if i in b:
          pass
        else:
          c += [i]
      d = [i + 1 for i in c]
      d = Permutation(d) 
      M.permute_columns(d)
      list_of_ranks += [len_of_R_pivots]
      R = M.matrix_from_columns(range(sum(list_of_ranks),n))
    return (M,list_of_ranks)


def maximum_projected_function_brouwer(ub,num_info_set):
    maximum_projected = 0
    lb = 1
    while lb < ub:
      maximum_projected += 1
      lb = num_info_set*(maximum_projected+1)
    return maximum_projected


def maximum_projected_function_zimmerman(ub,list_of_ranks,num_info_set,num_disj_info_set):
    maximum_projected = 0
    k = list_of_ranks[0]
    lb = 0 
    while lb < ub:
      maximum_projected += 1
      lb = (maximum_projected + 1)*num_disj_info_set 
      for m in xrange(num_disj_info_set,num_info_set):
        lb += max(0,maximum_projected+1-k+list_of_ranks[m])
    return maximum_projected


def delet_non_contibuting_matrix(list_of_ranks,num_info_set,num_disj_info_set,maximum_projected,verbose=True):
    k = list_of_ranks[0]
    L = copy(list_of_ranks)
    length = copy(num_info_set)
    for i in xrange(num_disj_info_set,length):
      count = 0
      for w in xrange(1,maximum_projected):
        if max(0,w+1-k+L[i]) != 0:
          break
        else: 
          count += 1
      if count == (maximum_projected-1):
        list_of_ranks.remove(L[i])
        num_info_set -= 1
        if verbose == True:
          print("Deleting non-contributing rank {} matrix".format(L[i]))
    return (list_of_ranks, num_info_set)


def list_of_system_gen_mat_zimm(M,list_of_ranks):
    n = M.ncols()
    k = M.nrows()
    l = len(list_of_ranks)
    L = []
    num_dis_is = 0
    num_over_is = 0

    for i in xrange(l) :
      if list_of_ranks[i] == k :
        num_dis_is += 1
      else :
        num_over_is += 1
    for i in xrange(num_dis_is) :
      A = M.matrix_from_columns(range(i*k, i*k + k))
      L += [A.inverse()*M]
    
    aux_ind = k*num_dis_is
    for i in xrange(-num_over_is,0) :
      I = range(aux_ind, aux_ind + list_of_ranks[i])
      A = M.matrix_from_columns(I)
      c = 0
      while A.rank() != k:
        if A.rank() < M.matrix_from_columns(I + [c]).rank():
          A = M.matrix_from_columns(I + [c])
          I += [c]
        c += 1
        if c == k :
          break
      L += [A.inverse()*M]
      aux_ind += list_of_ranks[i]
    return L


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


def minimum_distance_brouwer(C, nb_words=1, maxiter=20,verbose=True):

    n,k = C.length(), C.dimension()
    F = C.base_field()
    q = F.cardinality()
    w = 1
    nb_words = 0

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
    print("Code Minimum Weight Brouwer: length {}, dimension {}".format(n,k))
    print("The number of disjoint information set is : {} ".format(num_info_set))
    print("Lower Bound: {} , Upper Bound: {}".format(lb,ub))
    maximum_projected = copy(maximum_projected_function_brouwer(ub,num_info_set))
    print("Maximum Projected w : {}".format(maximum_projected))

    if F == GF(2) :
      while lb < ub :
        for m in xrange(0,num_info_set) : 
          A = L[m].row(0)
          for i in xrange(1,w):
            A += L[m].row(i)
          wt = A.hamming_weight()
          if ub > wt:
            ub = copy(wt)
            if verbose == True:
              print("New upper bound : {}".format(ub))
            if ub <= lb :
              return ub
            maximum_projected_inter = maximum_projected_function_brouwer(ub,num_info_set)
            if maximum_projected_inter != maximum_projected :
              maximum_projected = copy(maximum_projected_inter)
              if verbose == True :
                print("New Maximum Projected w : {}".format(maximum_projected))
          for i,j in combinations(k,w):
            A += L[m].row(i) + L[m].row(j)
            wt = A.hamming_weight()
            if ub > wt:
              ub = copy(wt)
              if verbose == True:
                print("New upper bound : {}".format(ub))
              if ub <= lb : 
                return ub
              maximum_projected_inter = maximum_projected_function_brouwer(ub,num_info_set)
              if maximum_projected_inter != maximum_projected :
                maximum_projected = copy(maximum_projected_inter)
                if verbose == True :
                  print("New Maximum Projected w : {}".format(maximum_projected))
          nb_words += binomial(k,w)  
          lb += 1
          if verbose == True:
            print("w : {}, lower: {}, upper: {}".format(w,lb,ub))
          if ub <= lb:
            if verbose == True:
              print("We have seen {} codewords ! ".format(nb_words))
            return ub
        if verbose == True:
          print("We have seen {} codewords ! ".format(nb_words))
        w += 1
      return ub
    
    q = F.cardinality()
    g = F.multiplicative_generator()
    M = [g^i for i in xrange(q-1)]
    Z = IntegerRing()

    while lb < ub :
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

          wt = A.hamming_weight()
          if ub > wt:
            ub = copy(wt)
            if verbose == True:
              print("New upper bound : {}".format(ub))
            if ub <= lb :
              return ub
            maximum_projected_inter = maximum_projected_function_brouwer(ub,num_info_set)
            if maximum_projected_inter != maximum_projected :
              maximum_projected = copy(maximum_projected_inter)
              if verbose == True :
                print("New Maximum Projected w : {}".format(maximum_projected))

          A_int = copy(A) 
          for i,j in combinations(k,w):
            X[j] = X[i]
            X[i] = F.zero()
            A_int += X[j]*(L[m].row(j) - L[m].row(i))
            wt = A_int.hamming_weight()
            if ub > wt:
              ub = copy(wt)
              if verbose == True:
                print("New upper bound : {}".format(ub))
              if ub <= lb :
                return ub
              maximum_projected_inter = maximum_projected_function_brouwer(ub,num_info_set)
              if maximum_projected_inter != maximum_projected :
                maximum_projected = copy(maximum_projected_inter)
                if verbose == True :
                  print("New Maximum Projected w : {}".format(maximum_projected))
        nb_words += binomial(k,w)*((q-1)^w)
        lb += 1
        if verbose == True:
          print("w : {}, lower: {}, upper: {}".format(w,lb,ub))
        if ub <= lb:
          if verbose == True:
            print("We have seen {} codewords ! ".format(nb_words))
          return ub
      if verbose == True:
        print("We have seen {} codewords ! ".format(nb_words))
      w += 1
    return ub


def minimum_distance_zimmermann(C,maxiter=20,verbose=True):

    G1 = C.generator_matrix()
    n,k = C.length(), C.dimension()
    F = C.base_field()
    G2, list_of_ranks = infomation_set_distance(G1,maxiter)
    num_info_set = len(list_of_ranks)
    L = list_of_system_gen_mat_zimm(G2,list_of_ranks)
    ub = n - k + 1
    w = 1
    nb_words = 0
    print("Code Minimum Weight Zimmermann: length {}, dimension {}".format(n,k))
    print("relative ranks used: {}".format(list_of_ranks))
    num_disj_info_set = 0
    for i in xrange(len(list_of_ranks)):
      if list_of_ranks[i] == k:
        num_disj_info_set += 1

    lb = num_disj_info_set         
    print("Lower Bound: {} , Upper Bound: {}".format(lb,ub))
    maximum_projected = copy(maximum_projected_function_zimmerman(ub,list_of_ranks,num_info_set,num_disj_info_set))
    print("Maximum Projected w : {}".format(maximum_projected))

    if F == GF(2) :
      while lb < ub :
        m = 0
        while m < num_info_set:
          A = L[m].row(0)
          for i in xrange(1,w):
            A += L[m].row(i)
          wt = A.hamming_weight()
          if ub > wt:
            ub = copy(wt)
            if verbose == True:
              print("New upper bound : {}".format(ub))
            if ub <= lb :
              return ub
            maximum_projected_inter = maximum_projected_function_zimmerman(ub,list_of_ranks,num_info_set,num_disj_info_set)
            if maximum_projected_inter != maximum_projected and maximum_projected_inter > 1:
              maximum_projected = copy(maximum_projected_inter)
              if verbose == True :
                print("New Maximum Projected w : {}".format(maximum_projected))
              if num_info_set != num_disj_info_set:
                list_of_ranks, num_info_set = delet_non_contibuting_matrix(list_of_ranks,num_info_set,num_disj_info_set,maximum_projected,verbose)
          for i,j in combinations(k,w):
            A += L[m].row(i) + L[m].row(j)
            wt = A.hamming_weight()
            if ub > wt:
              ub = copy(wt)
              if verbose == True:
                print("New upper bound : {}".format(ub))
              if ub <= lb :
                return ub
              maximum_projected_inter = maximum_projected_function_zimmerman(ub,list_of_ranks,num_info_set,num_disj_info_set)
              if maximum_projected_inter != maximum_projected and maximum_projected_inter > 1:
                maximum_projected = copy(maximum_projected_inter)
                if verbose == True:
                  print("New Maximum Projected w : {}".format(maximum_projected))
                if num_info_set != num_disj_info_set:
                  list_of_ranks, num_info_set = delet_non_contibuting_matrix(list_of_ranks,num_info_set,num_disj_info_set,maximum_projected,verbose)
          nb_words += binomial(k,w) 
          if max(0,w+1-k+list_of_ranks[m]) != 0:
            lb += 1
            if verbose == True:
              print("w : {}, lower: {}, upper: {}".format(w,lb,ub))
          if ub <= lb :
            if verbose == True:
              print("We have seen {} codewords ! ".format(nb_words))
            return ub 
          m += 1
        if verbose == True:
          print("We have seen {} codewords ! ".format(nb_words))
        w += 1
      return ub
    
    q = F.cardinality()
    g = F.multiplicative_generator()
    M = [g^i for i in xrange(q-1)]
    Z = IntegerRing()
       
    while lb < ub :
      X = [F.zero()]*k
      m = 0
      while m < num_info_set:
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

          wt = A.hamming_weight()
          if ub > wt:
            ub = copy(wt)
            if verbose == True:
              print("New upper bound : {}".format(ub))
            if ub <= lb :
              return ub
            maximum_projected_inter = maximum_projected_function_zimmerman(ub,list_of_ranks,num_info_set,num_disj_info_set)
            if maximum_projected_inter != maximum_projected and maximum_projected_inter > 1:
              maximum_projected = copy(maximum_projected_inter)
              if verbose == True :
                print("New Maximum Projected w : {}".format(maximum_projected))
              if num_info_set != num_disj_info_set:
                list_of_ranks, num_info_set = delet_non_contibuting_matrix(list_of_ranks,num_info_set,num_disj_info_set,maximum_projected,verbose)

          A_int = copy(A) 
          for i,j in combinations(k,w):
            X[j] = X[i]
            X[i] = F.zero()
            A_int += X[j]*(L[m].row(j) - L[m].row(i))
            wt = A_int.hamming_weight()
            if ub > wt:
              ub = copy(wt)
              if verbose == True:
                print("New upper bound : {}".format(ub))
              if ub <= lb :
                return ub
              maximum_projected_inter = maximum_projected_function_zimmerman(ub,list_of_ranks,num_info_set,num_disj_info_set)
              if maximum_projected_inter != maximum_projected and maximum_projected_inter > 1:
                maximum_projected = copy(maximum_projected_inter)
                if verbose == True:
                  print("New Maximum Projected w : {}".format(maximum_projected))
                if num_info_set != num_disj_info_set:
                  list_of_ranks, num_info_set = delet_non_contibuting_matrix(list_of_ranks,num_info_set,num_disj_info_set,maximum_projected,verbose)
        nb_words += binomial(k,w)*((q-1)^w)
        if max(0,w+1-k+list_of_ranks[m]) != 0:
          lb += 1
          if verbose == True:
            print("w : {}, lower: {}, upper: {}".format(w,lb,ub))
        if ub <= lb :
          if verbose == True:
              print("We have seen {} codewords ! ".format(nb_words))
          return ub
        m += 1
      if verbose == True:
        print("We have seen {} codewords ! ".format(nb_words))
      w += 1
    return ub


def test_rapide(): 
  # (GF(),long,dim)
  L = [(7,30,9),(3,60,15),(3,64,16),(3,68,17),(5,44,11),(5,48,12),(7,36,9),(7,40,10),(8,36,9),(8,40,10),(9,32,8),(9,36,9),(3,45,7),(11,33,5),(3,100,11),(11,44,6),(2,77,15),(2,100,11),(9,50,8),(5,44,5),(25,28,5),(2,100,25),(7,40,5),(17,15,4),(7,50,7),(11,50,5),(5,55,10),(5,55,9)]
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


def test_rapide_gf2():
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
    b = %time minimum_distance_brouwer(C,nb_words=1)
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
    # print("minimum_distance_brouwer(C) :")
    # b = %time minimum_distance_brouwer(C)
    # print b
    # print (" ")
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



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


def minimum_distance_brouwn(C):

    G1 = C.generator_matrix()
    k = G1.nrows()
    n = G1.ncols()
    G2 = infomation_set(G1)[0]
    num_info_set = infomation_set(G1)[1]
    ub = n - k + 1
    lb = num_info_set 
    L = []
    print("The number of disjoint information set is : {} ".format(num_info_set))
    
    for i in range(num_info_set):

      A = G2.matrix_from_columns(range(i*k , i*k + k))
      L = L + [A.inverse()*G2]
      print L[i]
      print " " 
    

    print L



C = codes.random_linear_code(GF(2),15,3)
G = C.generator_matrix()

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
    n = G.ncols()
    num_info_set = 0

    while M.matrix_from_columns(range(num_info_set*k,G.ncols())).rank() == k:

        a = M.matrix_from_columns(range(num_info_set*k,G.ncols())).pivots()
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

    print("the number of disjoint information set is : {} ".format(num_info_set))
    return M


def infomation_set_ancien(G):

    M = copy(G)
    k = M.nrows()
    n = M.ncols()
    res = copy(G)
    res = res*0
    num_info_set = 0

    for i in range(0,n//k):

        n = M.ncols()
        if M.rank() == k :

            num_info_set = num_info_set + 1
            a = M.pivots()
            a = list(a)
            print a
            L = M.matrix_from_columns(a)

            for l in range(k):
                for c in range(k):
                    res[l,i*k+c] = L[l,c]

            print L
            print " "
            print(L.rank())
            print " "

            b = range(n)
            c = [ i for i in b if i not in a ]
            print c
            M = M.matrix_from_columns(c)

    print("the number of disjoint information set is : {} ".format(num_info_set))
    return res.matrix_from_columns(range(num_info_set*k))


def minimum_distance_brouwn(C):
    return 0



M = random_matrix(GF(2),3,10)
C = codes.LinearCode(random_matrix(GF(5),3,20))
C = codes.LinearCode(random_matrix(GF(7**33),100,550))

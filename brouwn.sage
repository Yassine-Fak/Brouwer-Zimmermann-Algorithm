import sage.coding.linear_code
import sage.combinat.permutation


G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
C2 = LinearCode(G)
# print C2

C = codes.LinearCode(random_matrix(GF(2),3,10))
print(C)

print(C.generator_matrix())
print " "
print(C.systematic_generator_matrix())
print(C.minimum_distance())
print ""
print(C.basis())


# Comment on genere un code lineaire en utilisant uniquement le corps, long et la dim ?
# Pour permuter les colonnes on utilise ca : http://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix0.html
# permutation : http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/combinat.html#sage.combinat.combinat.CombinatorialElement
# Random : http://doc.sagemath.org/html/en/reference/misc/sage/misc/randstate.html

rnd = current_randstate().python_random()
m = rnd.randrange(100)
rnd.randrange(45,50)
print m

p = Permutations_CC(3)
p.list()
M.permute_columns(p.list()[3])

#Hello

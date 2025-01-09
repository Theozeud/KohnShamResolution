using KohnShamResolution

original_stdout = stdout
output_file = open("tests/Performance/basis.txt", "w")
redirect_stdout(output_file)

# Création des éléments puis de la base
T = Float64

# Tests integrals on elements
println("Elements")
@time intLeg = IntLegendreElements(T; ordermin = 2, ordermax = 2, normalize = false, binf = -T(1), bsup = T(1))



redirect_stdout(original_stdout)
close(output_file)
println("Performance finished")
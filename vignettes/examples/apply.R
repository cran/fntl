Rcpp::sourceCpp("apply.cpp")
X = matrix(1:12, 4, 3)
out = apply_ex(X)
print(out)


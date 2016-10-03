n = dim(dummyData)[2]

# DISTANCES = getXiCiDistances(X, obj$C)
DISTANCES = getXiDistancesRefactor(dummyData, dummyobj$C) 
# print("distances...")
# print(DISTANCES)
MEMBERSHIP = getContributionsRefactor(DISTANCES)#esta seria la funcion membership, la de contributions no esta
# CONTRIBUTIONS = getContributionsRefactor(DISTANCES)
# CONTRIBUTIONS <- contrib(MEMBERSHIP)
# print("contributions...")
# print(CONTRIBUTIONS)
contributions_phi_0 = MEMBERSHIP %*% as.matrix(dummyobj$phi_0)
# print("contributions_phi_0")
# print(contributions_phi_0)
X_PHI = dummyData %*% t(dummyobj$PHI)

phi_0_matrix <- matrix(rep(dummyobj$phi_0, times = nrow(X_PHI)), ncol = length(dummyobj$phi_0), byrow = TRUE)

PHI_0_X_PHI <- phi_0_matrix + X_PHI

X_PHI_CONTRIBUTIONS <- MEMBERSHIP * PHI_0_X_PHI
# print("X_PHI")
# print(X_PHI)
# y = apply(X_PHI, 1, sum) + contributions_phi_0
y = rowSums(X_PHI_CONTRIBUTIONS)

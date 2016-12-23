ob <- NULL
ob$C <- loadClusters(2,2)
ob$phi_0 = loadRandomVector(2)
ob$PHI = loadRandomPhi(2,2)
ob
X <- loadDataAmnfis(2)
X
amnfis.simulate(ob, X)

df <- data.frame(X)
df$y <- c(0,0,0,0,0,1,1,1,1,1)
df

mdl <- glm( y ~ . , data = df , family=binomial)

slope <- coef(mdl)[2]/(-coef(mdl)[3])
intercept <- coef(mdl)[1]/(-coef(mdl)[3])

library(lattice)
xyplot( X2 ~ X1 , data = df, groups = y,
        panel=function(...){
          panel.xyplot(...)
          panel.abline(intercept , slope)
          panel.grid(...)
        })



mdl <- amnfis(X, df$y, ob$C)
forecast <- amnfis.simulate(mdl, X, ob$C)

#' Hello world package
#'
#' AMNFIS program
#' @param k number of clusters
#' @return data
#' @export
#'
amnfis <- function(X, d, k){
  
  ################FUNCIONES AUXILIARES###################
  
  object2Params <- function(object){
    return(c(c(object$C), c(object$phi_0), c(object$PHI)))
  }
  
  params2Object <- function(object, params, k, n, m){
    i = 0
    
    object$C = matrix(data = params[i:(i + k * n)],
                      nrow = k, ncol = n)
    i = i + k * n + 1
    
    object$phi_0 = c(params[i: (i + k -1)])
    i = i + k
    
    object$PHI = matrix(data = params[i:(i + k * n - 1)],
                        nrow = k, ncol = n)
    i = i + k * n
    
    return(object)
  }
  
  fn.optim = function(v){
    obj = params2Object(obj, v, k, n, m)
    y = amnfis.simulate(obj, X)
    
    y[abs(y)<0.00001] = 0.00001
    y[abs(y)>0.99999] = 0.99999
    
    error = sum((d - y)^2)
    
    # error = -sum(d * log(y) + (1 - d) * log(1 - y))
    
    # cat(paste("error", error, "\n"))
    return(error)
  }
  
  ###################CUERPO DE LA FUNCIÓN###############
  m = dim(X)[1]
  n = dim(X)[2]
  obj = NULL
  
  #################CARGA DE PARÁMETROS ALEATORIOS(CLUSTERS, PHI_O Y PHI)###########
  obj$C = loadClusters(n, k)
  obj$phi_0 = loadRandomVector(k)
  obj$PHI = loadRandomPhi(k,n)
  ################## FIN CREACIÓN OBJECTOS ALEATORIOS#####################
  
  v = object2Params(obj)
  convergencia = FALSE
  while(convergencia == FALSE){
    a = optim(v, fn.optim, method = 'BFGS')
    if(a$convergence == 0){
      convergencia = TRUE
    }else{
      v = a$par
    }
    cat(paste("convergencia ", a$convergence))
  }
  # a = optim(v, fn.optim, control = list(maxit = 20000))
  # cat(paste("convergencia ", a$convergence))
  
  obj = params2Object(obj, a$par, k, n, m)
  return(obj)
  
}



amnfis.simulate <- function(obj, X) {
  
  n = dim(X)[2]
  
  # DISTANCES = getXiCiDistances(X, obj$C)
  DISTANCES = getXiDistancesRefactor(X, obj$C)
  # print("distances...")
  # print(DISTANCES)
  MEMBERSHIP = getContributionsRefactor(DISTANCES)#esta seria la funcion membership, la de contributions no esta
  # CONTRIBUTIONS = getContributionsRefactor(DISTANCES)
  CONTRIBUTIONS <- contrib(MEMBERSHIP)
  # print("contributions...")
  # print(CONTRIBUTIONS)
  contributions_phi_0 = CONTRIBUTIONS %*% as.matrix(obj$phi_0)
  # print("contributions_phi_0")
  # print(contributions_phi_0)
  X_PHI = X %*% t(obj$PHI)
  # print("X_PHI")
  # print(X_PHI)
  # y = apply(X_PHI, 1, sum) + contributions_phi_0
  y = rowSums(X_PHI) + contributions_phi_0
  y=1/(1+exp(-y))
  return(y)
}


loadData <- function(n){
  #print('creando la matriz de datos....')
  datos = c(0.1851488,0.3455233,-1.3720042,0.1961969,-0.2101043,-0.2215962,3.322623,-1.408931,-2.134558,-0.5122682,-0.1884232,0.5746452)
  # datos = c(-1,-1,-1,0,0,0,1,1,1,1,0,-1,1,0,-1,1,0,-1)
  # X = matrix(rnorm(12), ncol = n)#random values
  X = matrix(data = datos, ncol = n)
  return(X)
}

loadRandomVector <- function(size){
  return(rnorm(size))#random values
  # return(c(0.3773669,1.8119634))
}

loadRandomPhi <- function(k ,n){
  # d = c(-0.004729057,0.285950071,1.172539,0.1799518,0.1998691,-0.113449,-0.4013696,0.2711753)
  # phi_params = matrix(d, nrow = k, ncol = n)
  phi_params = matrix(rnorm(k * n), nrow = k, ncol = n)#random values
  return(phi_params)
}

loadClusters <- function(n, k){
  # d = c(-1.6260988,-0.5021358,1.2773631,-0.6000955,0.302947,0.7021099,1.5672107,0.757111)
  # CLUSTER = matrix(d, nrow = k, ncol = n, byrow = TRUE)
  CLUSTER = matrix(rnorm(n * k), nrow = k, ncol = n)#raandom values
  return(CLUSTER)
}

getXiCiDistances <- function(X, C){#prueba
  #   CT = t(C)#traspuesta de la matrix de clusters
  #   m = dim(X)[1]
  #   k = dim(CT)[1]
  #   mat = createMatrix(m,k)
  #   for(i in 1:m){
  #     for(j in 1:k){
  #       mat[i,j] = sum((X[i,] - CT[j,])^2)# calcula la distancia de cada dato de entrada a cada cluster
  #     }
  #   }
  D = rdist(X,C)^2
  return(D)
}

# getXiDistancesRefactor <- function(X, C){
#   m = dim(X)[1]
#   n = dim(X)[2]
#   k = dim(C)[1]
#   replByCol = rep(k, n)
#   replByRow = rep(m, n)
#   transformedX = X[,rep(1:n, replByCol)]
#   transformedC = matrix(rep(C, each=nrow(X)), nrow=m)#C[,rep(1:k, replByRow)]
#   dista = (transformedX - transformedC)^2
#   list_of_distances = split.along.dim(array(dista, c(m,k,n)), 3)
#   XtoCDistances = Reduce('+', list_of_distances)
#   return(XtoCDistances)
# }

getXiDistancesRefactor <- function(X, C){
  m = dim(X)[1]
  n = dim(X)[2]
  k = dim(C)[1]
  replByCol = rep(k, n)
  replByRow = rep(m, n)
  # print(replByCol)
  transformedX = X[,rep(1:n, replByCol)]
  # print(dim(transformedX))
  transformedC = matrix(rep(C, each=nrow(X)), nrow=m)#C[,rep(1:k, replByRow)]
  # print(dim(transformedC))
  dista = (transformedX - transformedC)^2
  distaces3D = array(dista, c(m,k,n))
  distancesList = lapply(seq(dim(distaces3D)[3]), function(x) distaces3D[ , , x])
  rsp = Reduce('+', distancesList)
  return(rsp)
}

#Recibe la matriz de distancias y calcula la contribucion de cada dato a cada cluster
getContributions <- function(distances){
  contr = apply(distances, 1, calculateContributions)#1 significa que la funcion aplica  por filas
  return(t(contr))
}

getContributionsRefactor <- function(distances){
  return(exp(-distances/rowSums(distances)))
}

contrib <- function(membership){
  return(membership/rowSums(membership))
}

calculateContributions <- function(X){
  return(exp(-X/sum(X)))
}

createMatrix <- function(m, k){
  x <- matrix()
  length(x) <- m*k
  dim(x) <- c(m,k)
  return(x)
}


#FUNCIONES DE R SOBRECARGADAS PARA AMNFIS

plot = function(x, y){
  UseMethod("plot")
}

predict = function(x){
  UseMethod("predict")
}

summary = function(model){
  UseMethod("summary")
}

plot.amnfis <- function(x,y) {
  print('graficando los datos...........')
  #lims=x["mean"]+c(-3*x["sd"],3*x["sd"])
  #h=seq(lims[1],lims[2],length=1000)
  #plot(x, y)
}

predict.amnfis <- function(x){
  print("Ejecutando la prediccion para el nuevo dato...")
}

summary.amnfis <- function(x){
  print("obteniendo el resumen del modelo...")
}

loadArtificialData = function(){
  data1 = matrix(rnorm(60, mean = 2, sd = 1), nrow = 30)
  data2 = matrix(rnorm(60, mean = 10, sd = 1), nrow = 30)
  data3 = cbind(rnorm(30, mean = 20, sd = 1), rnorm(30, mean = 1, sd = 1))
  mydata = rbind(data1, data2, data3)
}

loadArtificialClasses = function(size1, size2){
  class1 = rep(1, size1)
  class2 = rep(0, size2)
  return(c(class1, class2))
}

loadFixedData = function(){
  vec1 = c(9.6988438,10.8130845,9.960313,9.177409,9.027697,10.1697067,11.1747064,8.8430911,9.6841149,10.299404,9.9416978,11.1079718,8.5221696,9.1535845,11.5671524,10.0640491,8.4345442,10.5505362,10.7465217,10.0916271,9.5827665,9.4756067,10.8772088,10.9978583,9.48498,10.8880014,11.64234,10.3212844,9.1458854,8.6715372,19.4580347,19.0518549,19.5469181,20.0699309,19.643371,18.926044,19.8861368,18.9738247,19.5671417,21.5671072,20.4963176,22.1814699,20.1566811,20.4750505,19.5906645,19.4081632,20.202404,20.2465333,22.5983501,19.37832,22.3654685,20.2884596,19.7513829,21.4054635,19.361546,19.2219921,19.3747091,22.3756983,20.5225083,20.738964)
  vec2 = c(10.04080547,10.90807896,11.53085361,10.09048219,10.9640893,10.25135567,9.3892929,9.10601038,10.43391051,12.05139659,10.5230589,10.47608986,9.98663132,10.33134773,9.80227397,10.84523407,9.28804128,9.10730992,11.17051018,9.05305146,8.64459723,10.41679116,10.95222898,9.90234401,9.00651082,11.71960297,8.60982494,9.01929174,11.15920278,10.21887356,0.83363395,1.80395437,1.39166315,0.984936,1.72303215,1.33548702,2.63357547,3.56152338,1.98274043,0.87498201,1.14437694,-0.65187753,0.38379808,-1.53358693,0.95222896,0.67539677,-0.15998695,2.32643136,1.52044043,0.68562911,2.29842009,0.95052233,1.95219523,0.71825588,0.28510547,0.15616758,1.30461723,3.06204045,2.0675249,0.50490544)
  return(cbind(vec1, vec2))
}

loadDataCl1 <- function(){
  vec1 = c(2.403725,1.3052997,2.7338131,3.0022611,2.7473139,1.251366,1.2532291,1.7864044,3.9762732,0.7178647,2.426671,2.2872671,0.4712174,2.4314146,1.2573021,1.4596345,2.2750345,2.1867702,1.2202835,-0.6392887,3.0972038,1.6295169,1.6267859,0.3825643,3.5741614,2.4688016,3.1513287,3.1351092,1.9075098,4.9351927,9.398439,10.9528082,10.5734055,8.5020511,9.1213917,9.028043,9.1799912,8.6294745,11.2137843,9.2074718,9.828419,8.6859078,9.9444845,10.4913392,9.4520827,8.9286116,9.8799033,11.1521401,11.180486,8.9210758,10.0528612,9.8997657,11.9677087,8.9639066,9.6660183,9.8360467,9.6075269,11.4688543,8.4732205,8.461425,20.269285,21.4261585,19.9393386,18.7681032,20.0973043,18.8299044,19.8298618,19.2246872,19.2869621,18.7931123,21.3587827,21.5778189,20.3389076,20.0671534,18.4843695,20.5278038,21.3529175,19.6376293,18.6931801,20.1750783,18.2028874,18.8688394,18.5279261,20.5950625,20.6789979,19.2329008,21.224626,19.9253769,21.330975,20.9652148)
  vec2 = c(2.83312035,1.98705516,1.14535486,2.89364582,2.82606916,1.97789791,1.48671422,0.54866055,1.09799551,2.82238551,-0.07050465,1.85106387,2.31536926,1.16871121,2.51469621,2.30420803,3.10601673,1.40525553,2.41045206,2.05224656,0.68443797,3.65713039,1.87919346,2.371447,1.98832866,2.391676,3.26087993,2.27012832,1.74418101,1.59412108,11.53245831,11.32887022,8.68404138,10.35421053,9.19984142,9.73681168,10.08435248,11.50240544,9.8304405,10.47273373,10.61924175,10.10919631,8.97273746,10.72892394,10.31361009,7.95199839,10.02681751,9.65728098,10.84036915,9.84714582,10.1410311,11.1557943,9.52275695,9.73095296,10.60085567,8.22853417,9.26364452,10.28802238,9.55457077,11.23637841,0.84973674,1.71499042,0.42729809,-1.0347538,1.37927908,2.89945766,2.25129674,0.91244421,1.4043183,0.74051278,1.58766866,0.96641317,2.48679601,2.94983522,-0.0280081,0.55777286,1.71400481,-0.0105492,1.55403776,1.92140898,2.66255747,0.56418209,2.6138638,-0.65089549,0.20799128,0.45497272,1.25818523,-0.25817157,0.82367951,-0.12530707)
  return(cbind(vec1, vec2))
}

loadDataCl2 <- function(){
  vec1 = c(9.398439,10.9528082,10.5734055,8.5020511,9.1213917,9.028043,9.1799912,8.6294745,11.2137843,9.2074718,9.828419,8.6859078,9.9444845,10.4913392,9.4520827,8.9286116,9.8799033,11.1521401,11.180486,8.9210758,10.0528612,9.8997657,11.9677087,8.9639066,9.6660183,9.8360467,9.6075269,11.4688543,8.4732205,8.461425,20.269285,21.4261585,19.9393386,18.7681032,20.0973043,18.8299044,19.8298618,19.2246872,19.2869621,18.7931123,21.3587827,21.5778189,20.3389076,20.0671534,18.4843695,20.5278038,21.3529175,19.6376293,18.6931801,20.1750783,18.2028874,18.8688394,18.5279261,20.5950625,20.6789979,19.2329008,21.224626,19.9253769,21.330975,20.9652148)
  vec2 = c(11.53245831,11.32887022,8.68404138,10.35421053,9.19984142,9.73681168,10.08435248,11.50240544,9.8304405,10.47273373,10.61924175,10.10919631,8.97273746,10.72892394,10.31361009,7.95199839,10.02681751,9.65728098,10.84036915,9.84714582,10.1410311,11.1557943,9.52275695,9.73095296,10.60085567,8.22853417,9.26364452,10.28802238,9.55457077,11.23637841,0.84973674,1.71499042,0.42729809,-1.0347538,1.37927908,2.89945766,2.25129674,0.91244421,1.4043183,0.74051278,1.58766866,0.96641317,2.48679601,2.94983522,-0.0280081,0.55777286,1.71400481,-0.0105492,1.55403776,1.92140898,2.66255747,0.56418209,2.6138638,-0.65089549,0.20799128,0.45497272,1.25818523,-0.25817157,0.82367951,-0.12530707)
  return(cbind(vec1, vec2))
}

#para crear un array de matrices (diferentes dimensiones): matrices = array(dst, c(3,2,4))
#funcion para convertir un arrray de matrices en una lista matrices
split.along.dim <- function(a, n)
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])

#Para sumar las matrices dentro de una lista de matrices 
#Reduce('+', matrices)

test.apply <- function(data){
  print('aaa')
  print(data)
}
mysuma <- function(mat1, mat2){
  print("entrando a mysuma")
  print(mat1)
  print(mat2)
}

cuadrado <- function(datos){
  #   print("datos")
  #   print(datos^2)
  return(datos^2)
}

append.Rda <- function(..., file) {
  old.objects <- load(file, new.env())
  save(list = c(old.objects, ...), file = file)
}

getData <- function(df, bound) {
  # bound <- floor((nrow(df)/4)*3)         #define % of training and test set
  
  df <- df[sample(nrow(df)), ]           #sample rows 
  df.train <- df[1:bound, ]              #get training set
  df.test <- df[(bound+1):nrow(df), ]    #get test set
  return(list(df.train, df.test))
}

calculateCEC <- function(output_test_set, real_output){
  CEC <- sqrt(sum((real_output - output_test_set)^2) / length(real_output))
  return(CEC);
}

calculateUC <- function(
  output_A_from_B, output_A_from_A, output_B_from_A, output_B_from_B
){
  sum_error_A <- sum((output_A_from_B - output_A_from_A)^2)
  sum_error_B <- sum((output_B_from_A - output_B_from_B)^2)
  return(sqrt(sum_error_A + sum_error_B))
}

calculateRC <- function(output_A, output_B_from_A, output_B, output_A_from_B){
  # sum_errors_B <- (output_A - )
  # verificar si se puede realizar esta formula ya que los vectores son de diferentes tamanios
}

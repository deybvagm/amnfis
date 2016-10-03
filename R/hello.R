#' Hello world package
#'
#' AMNFIS program
#' @param k number of clusters
#' @return data
#' @export
#'
amnfis <- function(X, d, clusters){
  
  k = nrow(clusters)
  
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
  # obj$C = loadClusters(n, k) # TODO aca se deben calcular los clusters pero no de forma aleatoria
  obj$C <- clusters
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
    # cat(paste("convergencia ", a$convergence))
  }
  # a = optim(v, fn.optim, control = list(maxit = 20000))
  # cat(paste("convergencia ", a$convergence))
  
  obj = params2Object(obj, a$par, k, n, m)
  return(obj)
  
}



amnfis.simulate <- function(obj, X) {
  
  # n = dim(X)[2]
  # 
  # # DISTANCES = getXiCiDistances(X, obj$C)
  # DISTANCES = getXiDistancesRefactor(X, obj$C)
  # # print("distances...")
  # # print(DISTANCES)
  # MEMBERSHIP = getContributionsRefactor(DISTANCES)#esta seria la funcion membership, la de contributions no esta
  # # CONTRIBUTIONS = getContributionsRefactor(DISTANCES)
  # CONTRIBUTIONS <- contrib(MEMBERSHIP)
  # # print("contributions...")
  # # print(CONTRIBUTIONS)
  # contributions_phi_0 = CONTRIBUTIONS %*% as.matrix(obj$phi_0)
  # # print("contributions_phi_0")
  # # print(contributions_phi_0)
  # X_PHI = X %*% t(obj$PHI)
  # X_PHI_CONTRIBUTIONS <- CONTRIBUTIONS * X_PHI
  # # print("X_PHI")
  # # print(X_PHI)
  # # y = apply(X_PHI, 1, sum) + contributions_phi_0
  # y = rowSums(X_PHI_CONTRIBUTIONS) + contributions_phi_0
  
  n = dim(X)[2]

  # DISTANCES = getXiCiDistances(X, obj$C)
  DISTANCES = getXiDistancesRefactor(X, obj$C)
  # print("distances...")
  # print(DISTANCES)
  MEMBERSHIP = getContributionsRefactor(DISTANCES)#esta seria la funcion membership, la de contributions no esta
  # CONTRIBUTIONS = getContributionsRefactor(DISTANCES)
  # CONTRIBUTIONS <- contrib(MEMBERSHIP)
  # print("contributions...")
  # print(CONTRIBUTIONS)
  contributions_phi_0 = MEMBERSHIP %*% as.matrix(obj$phi_0)
  # print("contributions_phi_0")
  # print(contributions_phi_0)
  X_PHI = X %*% t(obj$PHI)

  phi_0_matrix <- matrix(rep(obj$phi_0, times = nrow(X_PHI)), ncol = length(obj$phi_0), byrow = TRUE)

  PHI_0_X_PHI <- phi_0_matrix + X_PHI

  X_PHI_CONTRIBUTIONS <- MEMBERSHIP * PHI_0_X_PHI
  # print("X_PHI")
  # print(X_PHI)
  # y = apply(X_PHI, 1, sum) + contributions_phi_0
  y = rowSums(X_PHI_CONTRIBUTIONS)
  
  y=1/(1+exp(-y))
  y = transform_output(y)
  return(y)
}


loadDataAmnfis <- function(n){
  #print('creando la matriz de datos....')
  datos = c(0.1851,0.3455,-1.3720,0.1961,-0.2101,-0.2215,3.3226,-1.4089,-2.1345,-0.5122,-0.1884,0.5746)
  # datos = c(-1,-1,-1,0,0,0,1,1,1,1,0,-1,1,0,-1,1,0,-1)
  # X = matrix(rnorm(12), ncol = n)#random values
  X = matrix(data = datos, ncol = n)
  return(X)
}

loadRandomVector <- function(size){
  return(rnorm(size))#random values
  # return(runif(size, min = 0, max = 1))
  # return(c(0.3773,1.8119))
}

loadRandomPhi <- function(k ,n){
  # d = c(-0.0047,0.2859,1.1725,0.1799,0.1998,-0.1134,-0.4013,0.2711)
  # phi_params = matrix(d, nrow = k, ncol = n)
  phi_params = matrix(rnorm(k * n), nrow = k, ncol = n)#random values
  # phi_params = matrix(runif(k * n), nrow = k, ncol = n)#random values
  return(phi_params)
}

loadClusters <- function(n, k){
  d = c(-1.6260,-0.5021,1.2773,-0.6000,0.3029,0.7021,1.5672,0.7571)
  CLUSTER = matrix(d, nrow = k, ncol = n, byrow = TRUE)
  # CLUSTER = matrix(rnorm(n * k), nrow = k, ncol = n)#raandom values
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
  if(is.vector(distances)){
    distances <- matrix(distances)
  }
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
  sum_error_A <- sum((output_B_from_A - output_A_from_A)^2)
  sum_error_B <- sum((output_A_from_B - output_B_from_B)^2)
  return(sqrt(sum_error_A + sum_error_B))
}

# TODO la formula no concuerda con la del paper. Validar!!
calculateRC <- function(output_A, output_B_from_A, output_B, output_A_from_B){
  sum_errors_A <- sum((output_A - output_B_from_A)^2)
  sum_errors_B <- sum((output_B - output_A_from_B)^2)
  return(((sum_errors_A / length(output_A)) + (sum_errors_B / length(output_B))) / 2)
  # verificar si se puede realizar esta formula ya que los vectores son de diferentes tamanios
}

calculateSR <- function(output_A, output_B_from_A, output_B, output_A_from_B){
  return(sqrt(calculateRC(output_A, output_B_from_A, output_B, output_A_from_B)))
}

transform_output <- function(output){
  rsp <- ifelse(as.vector(output) > 0.5, 1, 0)
  return(rsp)
}

# ejecuta un clasificador lineal para cada particion y suma sus errores
binary_classification_error <- function(formula, first_partition, second_partition) {
  return(fit_linear_classifier(formula, first_partition) + fit_linear_classifier(formula, second_partition))
}



# TODO eliminar
# fn.findclusters <- function(formula, data, n){
#   model1 <- unname(fn.splitdata(formula, data))
#   part1 <- data[model1[[1]]$part1,]
#   part2 <- data[model1[[1]]$part2,]
#   for(i in 1:n-1){
#     model_p1 <- fn.splitdata(formula, part1)
#     model_p2 <- fn.splitdata(formula, part2)
#     selected_model <- model_p1 # Supongase que es el que tiene el menor error
#     
#   }
# }







# Lamma la funcion que agrupa los datos, toma el promedio por columnas de cada grupo para definir los centroides
fn.getcentroids <- function(all_data, n, formula){
  
  fn.groupdata <- function(all_data, n, formula){
    
    if(n ==1){
      data <- list(all_data)
    }else{
      models <- list()
      
      initial_model <- fn.splitdata(formula, all_data)
      initial_model <- initial_model[[1]]
      part1 <- initial_model$data[initial_model$part1,]
      part2 <- initial_model$data[initial_model$part2,]
      
      idx_loop <- 3
      
      while(idx_loop <= n){
        # print("entro a ciclo")
        print(n)
        model_a <- fn.splitdata(formula, part1)
        model_b <- fn.splitdata(formula, part2)
        model_a <- model_a[[1]] # TODO Hacer refactor para mejorar esto. En modelo hay una lista con un solo objeto
        model_b <- model_b[[1]] # TODO Hacer refactor para mejorar esto. En modelo hay una lista con un solo objeto
        
        if(length(models) == 0){
          models <- list(model_a, model_b)
        }else{
          models[[length(models) + 1]] <- model_a
          models[[length(models) + 1]] <- model_b
        }
        errors <- lapply(models, fn.fetcherrors)
        idx_min_error <- which.min(errors)
        obj <- models[[idx_min_error]]
        models[[idx_min_error]] <- NULL
        part1 <- obj$data[obj$part1,] # TODO enviar data en la funcion que encuentra las mejores particiones
        part2 <- obj$data[obj$part2,]
        idx_loop <- idx_loop + 1
        # fn.findclusters2(part1, part2, models, n+1, formula)
        # part1 <-
        # TODO Se deberían permitir clusters de un solo dato??
        
        # models <- ifelse(length(models) == 0, list(subpart_a, subpart_b), models[[length(models) + 1]])
      }
      data <- lapply(models, fn.datafrommodels)
      data[[length(data) + 1 ]] <- part1
      data[[length(data) + 1 ]] <- part2
    }
    return(data)
  }
  
  fn.splitdata <- function(formula, data){
    original <- data
    mat <- apply(as.matrix(data[,1:ncol(data) -1]), 2, fn.getpartitions, formula = formula, data = data)
    errors <- unlist(lapply(mat, fn.fetcherrors))
    # print("errores a comparar en la capa superior")
    # print(errors)
    lowest_error_idx <- which.min(errors)
    return(unname(mat[lowest_error_idx]))
  }
  
  fn.fetcherrors <- function(object){
    # print("dato")
    # print(o$err)
    return(object$error)
  }
  
  fn.datafrommodels <- function(model){
    return(model$data)
  }
  
  fn.getpartitions <- function(v, formula, data){
    # print("formula")
    # print(formula)
    v_replicated <- fn.replicatedata(v)
    # print("v_replicated")
    # print(v_replicated)
    partitions <- t(ifelse(t(v_replicated) >= v, 1, 0))
    # print('partition')
    # print(partitions)
    fit_errors <- apply(partitions, 2, fn.fitdata, formula = formula, data = data)
    errors <- lapply(fit_errors, fn.fetcherrors)
    errors <- unlist(errors)
    lowes_error_idx <- which.min(errors)
    best_partition <- fit_errors[[lowes_error_idx]]
    # print("fit errors")
    # print(fit_errors)
    # print("min error")
    # best_partition <- min(fit_errors$error)
    # print(fit_errors)
    return(best_partition)
  }
  
  # Parte en dos conjuntos de datos, para entrenar y obtener los errores. Esto se hace por cada columna
  fn.fitdata <- function(v_logic, formula, data){
    # print("v_logic")
    # print(v_logic)
    part1 <- which(v_logic == 1, arr.ind = TRUE) # Obtiene los indices de las posiciones que tienen valor 1
    part2 <- which(v_logic == 0, arr.ind = FALSE) # Obtiene los indices de las posiciones que tienen valor 0
    # print("part1")
    # print(part1)
    # print("part2")
    # print(part2)
    data_part1 <- data[part1,] # Saca la primera particion del conjunto nde datos
    data_part2 <- data[part2,] # Saca la segunda particion del conjunto de datos
    # print("data_part1")
    # print(data_part1)
    # print("data_part2")
    # print(data_part2)
    error1 <- ifelse((nrow(data_part1) > 0 && nrow(data_part2) > 0), fit_linear_classifier(formula, data_part1), 100)#TODO Se deberia dar un error alto tambien para los conjuntos que tengan un solo dato
    error2 <- ifelse((nrow(data_part1) > 0 && nrow(data_part2) > 0), fit_linear_classifier(formula, data_part2), 100)# TODO eliminar valores quemados
    # obj <- NULL //TODO mejorar. Intentar retornanado una lista ref: http://stackoverflow.com/questions/9497114/how-do-i-make-an-array-of-classes-in-r
    # obj$error <- error1 + error2
    # obj$part1 <- part1
    # obj$part2 <- part2
    obj <- list(error = error1 + error2, part1 = part1, part2 = part2, data = data)
    return(obj)
  }
  
  # Replicalos datos de una columna del dataframe en una matriz cuadrada
  fn.replicatedata <- function(v){
    n <- length(v)
    return(matrix(rep(v, each = n), ncol = n, byrow = TRUE))
  }
  
  # Remueve la columna que representa la clase en el vector(es un vector nombrado no un dataframe lo que entra)
  fn.removeclass <- function(v){
    v <- v[1:length(v) - 1]
    return(v)
    # dataframe[, ncol(dataframe)] <- NULL
  }
  
  # Ajusta un clasificador lineal para una muestra de datos y retorna el error
  fit_linear_classifier <- function(formula, data){
    y <- data[,ncol(data)]
    response <- glm(formula, data = data)
    mse <- fn.error(y, response$fitted.values)
    return(mse)
  }
  
  # Funcion de error que calcula el MSE
  fn.error <- function(y, predictions){
    mse <- sum((y - predictions)^2)/length(y)
    return(mse)
  }
  
  groups <- fn.groupdata(all_data, n, formula)
  centroids <- lapply(groups, colMeans)
  columns <- ncol(all_data) - 1
  centroids <- lapply(centroids, fn.removeclass)
  centroid_matrix <- matrix(unlist(centroids), ncol = columns, byrow = TRUE)
  return(centroid_matrix)
}

fn.calculate_acc <- function(formula, n, df, X, d, X_test, y){
  library(doParallel)
  registerDoParallel(cores = 4)
  centroids <- fn.getcentroids(all_data = df, n = 2, formula = formula) #Obtiene los dos primeros centroides
  model <- amnfis(X = X, d = d, clusters = centroids)
  forecast <- amnfis.simulate(obj = model, X = X_test)
  acc <- length(y[y == forecast]) / length(y)
  result <- NULL
  result$acc <- c(acc)
  for (j in 3:n) {
    acc <- 0
    vec_best_data_point <- c()
    foreach (i = 1:nrow(X)) %dopar% {
      vec_data_point <- X[i,]
      centroids <- rbind(centroids, vec_data_point)
      model <- amnfis(X = X, d = d, clusters = centroids)
      forecast <- amnfis.simulate(obj = model, X = X_test)
      accuracy <- length(y[y == forecast]) / length(y)
      if(accuracy > acc){
        acc <- accuracy
        vec_best_data_point <- vec_data_point
      }
      centroids <- centroids[-nrow(centroids),] #Retira el ultimo dato para probar con otro centroide
    }
    centroids <- rbind(centroids, vec_best_data_point)
    result$acc <- c(result$acc, acc)
  }
  result$clusters <- centroids
  return(result)
}

normalize <- function(x){
  return((x - min(x)) / (max(x) - min(x)))
}











# TODO evaluar la posibilidad de hacer escalamiento de variables


# Esta seccions solo es de prueba
mytest <- function(v){
  err <- rnorm(1)
  c1 <- c(1,2,4)
  c2 <- c(3,5)
  obj <- list(err = err, c1 = c1, c2 = c2)
  
  # obj$err <- err
  # obj$c1 <- c1
  # obj$c2 <- c2
  # obj <- 2
  return(obj)
}


#para agrega objeto a lista:
# lst[[length(lst) + 1]] <- obj
# para eliminar list[[1]] <- NULL
# ref: http://stackoverflow.com/questions/17046336/here-we-go-again-append-an-element-to-a-list-in-r

# Aplica una funcion por columnas para retornar las particiones que seran probadas
getMatrix <- function(column){
  mtx <- numeric()
  for (x in column){
    mtx <- cbind(mtx, column >= x)
  }
  return(mtx)
}

# Aplica una funcion por cada columna de la matriz de particiones, realizando las particiones y ejecutando el clasificador lineal
find_partitions <- function(column, datos, y){
  datos <- cbind(datos, y)
  part1 <- matrix(datos[column==1], ncol = 4)
  part1 <- as.data.frame(part1)
  part2 <- matrix(datos[column==0], ncol = 4)
  part2 <- as.data.frame(part2)
  # response <- glm(formula = V4~., data = part1)
  # response2 <- glm(formula = V4~., data = part2)
  # print(length(part1$V4))
  # print(length(response$fitted.values))
  mse <- get_prediction_error(part1)
    # fn.error(part1$V4, response$fitted.values)
  # mse2 <- fn.error(part2$V4, response2$fitted.values)
  mse2 <- get_prediction_error(part2)
  print(mse2)
  print(mse)
  return(mse + mse2)
}

get_prediction_error <- function(partition){
  error <- ifelse(nrow(partition) == 0, 0, fit_data(partition))
  return(error)
}

fit_data <- function(partition){
  ft <- glm(formula = V4~., data = partition)
  mse <- fn.error(partition$V4, ft$fitted.values)
  return(mse)
}

result <- c()

find.betterpartition <- function(formula, dataframe, X, d, X_test, y_true){
  print('ooo')
  # res <- c()
  for(i in 1:30){
    centroids <- fn.getcentroids(all_data = dataframe, n = i, formula = formula)
    model <- amnfis(X = X, d = d, clusters = centroids)
    forecast <- amnfis.simulate(obj = model, X = X_test)
    accuracy <- length(y_true[y_true == forecast]) / length(y_true)
    result <- c(result, accuracy)
    # res <- c(res, accuracy)
  }
  return(result)
}

# ionosphere - 19 clusters, tr: 200; validation: 151
# CEC = 0.6407 
# UC = 12.12
# RC = 0.392798 TODO ojo que la formula no concuerda con la del paper
# SR = 0.62
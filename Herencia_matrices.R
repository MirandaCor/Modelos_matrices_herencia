#Relevancia: Se emplean la diagonalizacion de matrices, valores y vectores propios


#Planteamiento del problema: Calcular una matriz con padres con gentotipo A xAa.
#Las diferentes posibiliddes para los machos son A o a y para las hembras AA o aa
# con probabilidades de 1/4 equitativamente para genotipo A x AA,A x Aa, a x AA, 
#a x Aa.Obteniendo el vector planteado en la matriz A, considerar dominancia de A.

#       A × AA - A × Aa - A × aa- a × AA- a × Aa -a × aa
#A × AA  1        0.25       0      0       0       0
#A × Aa  0        0.25       0      1       0.25    0
#A × aa  0        0          0      0       0.25    0
#a × AA  0        0.25       0      0       0       0
#a × Aa  0        0.25       1      0       0.25    0
#a × aa  0        0          0      0       0.25    1

library(matlib) #libreria para funcion inv

her_crom<- function(n){ 
#matriz de probabilidades
A<- as.matrix(data.frame(c(1,0,0,0,0,0),c(0.25,0.25,0,0.25,0.25,0),
                         c(0,0,0,0,1,0),c(0,1,0,0,0,0) ,c(0,0.25,0.25,0,0.25,0.25),
                         c(0,0,0,0,0,1)))
e <- eigen(A) # valores y vectores propios
e_vec<- as.matrix(e$vectors) 
e_vec_inv<- inv(e_vec) # matriz inversa de vectores propios
d<-( e_vec_inv %*%A) %*% e_vec # P-1 * A * P
d_n<-d^n # matriz diagonal elevada a la n
x_n<- (e_vec %*% d_n )%*% e_vec_inv # P* d^n * P-1

#si suponemos que la pareja original es de tipo A x Aa entonces
#multiplicamos por (0,1,0,0,0,0)

#  fenotipo de padres =c(A×AA,A×Aa,A×aa,axAA,a×Aa,a×aa)
gen<- x_n%*% c(0,1,0,0,0,0) # fenotipo de padres AxAa

row.names(gen)<-c("A×AA","A×Aa","A×aa","axAA","a×Aa","a×aa")

colnames(gen)<-c("Genotipo_descendencia")
round (gen,2)
}

her_crom (100) # n= generaciones, 
                      

# Tenemos que en 100 generaciones la posibilidad
#de tener que pareja hermano y hermana sean del tipo A x Aa en 1/3 
#y tipo a x aa solo de 1/3


#conclusión:se comprobo que es replicable con funciones de Rstudio el modelo matematico
#de distribucion de genotipos en la descendencia.


# referencias:

#http://matema.ujaen.es/jnavas/web_modelos/pdf_mmb08_09/modelos%20matriciales.pdf
#https://cran.r-project.org/web/packages/matlib/vignettes/inv-ex1.html

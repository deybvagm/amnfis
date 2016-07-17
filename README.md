# amnfis
2016-06-25
  Los clusters se toman teniendo en cuenta los datos. ¿Los valores phi y phi_0 tambien o estos si debe ser aleatorios?. En caso que deban ser aleatorios, cuál es el rango, deben ir de acuerdo al rango de los datos? ¿Se deben escalar los valores?

  Se toman solo los dtaos de entrenamiento para detectar los centros iniciales



------------------------
-----MEJORAS------------
------------------------
1- Escalar los valores de los conjuntos de datos
2- Tener en cuenta este rango de valores de los datos escalados para generar los valores phi en este rango

-----------------------------
----OPERTUNIDAD DE MEJORA---
-----------------------------
1-Usar conjunto de validación cruzada (no entrenar siempre con los mismos datos, puede haber sesgo)
  Esto involucra hacer pruebas no con los datos ordenados sino con toma de muestras?
2-Probar con un ejemplo artificial de clasificación lineal

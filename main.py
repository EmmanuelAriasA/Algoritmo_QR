from math import sqrt  # para raiz cuadrada

# imprimir la matriz


def imprMatriz(matriz, m):
    for i in range(0, m):
        for j in range(0, m):
            # imprime el numero redondead0 por 5 cifras
            print("(", round(matriz[i][j], 5), end=")")
        print()

# --------------------------------------------------------------------------------------------


def gramschmidt(matriz, m, matrizOrtogonal):
    ki = 0
    vi = []
    qi = []
    qj = []
    gram = []
    for i in range(0, m):
        gram.append(0)
        vi.append(0)
        qi.append(0)
        qj.append(0)
    contqi = 0
    contqj = 1

    # Paso 1: q1=v1
    for i in range(0, m):
        matrizOrtogonal[i][0] = matriz[i][0]
    contqi = 1

    # paso 2:
    while(contqi < m):

        # obtener las variables que necesitamos apara aplicar la formula del PDF 29
        # vi
        for i in range(0, m):
            vi[i] = matriz[i][contqi]

        # qj
        while(contqj <= contqi):
            for i in range(0, m):
                qj[i] = matrizOrtogonal[i][contqj-1]
            contqj = contqj+1

            # ki = <vi,qj>/<qj,qj>
            ki = KI(vi, qj, m)

            for i in range(0, m):
                gram[i] = gram[i] - (ki*qj[i])

        # se suma vi+gram
        for i in range(0, m):
            qi[i] = vi[i]+gram[i]

        # agregar qi a la matriz ortogonal
        for i in range(0, m):
            matrizOrtogonal[i][contqi] = qi[i]

        # reestablacer valores para una nueva iteración
        contqi = contqi+1
        contqj = 1
        for i in range(0, m):
            gram[i] = 0

# --------------------------------------------------------------------------------------------
# Proceso de ortonormalizacion de la matriz ortogonal


def ortonormalizar(matrizOrtogonal, m, matrizOrtonormal):
    # declaracion de variables
    cont = 0
    ui = []
    wi = []
    wiN = 0
    for i in range(0, columnas):
        ui.append(0)
        wi.append(0)

    for i in range(0, m):
        # obtencion de variables para hacer ui=wi/wiNorma
        # wi
        for b in range(0, m):
            wi[b] = matrizOrtogonal[b][cont]

        # wiNorma
        for b in range(0, m):
            wiN = wiN+(wi[b]*wi[b])
        wiN = sqrt(wiN)

        # ui
        for b in range(0, m):
            ui[b] = (wi[b])/wiN

        # agregar ui a la matriz ortonormal
        for b in range(0, m):
            matrizOrtonormal[b][cont] = ui[b]

        cont = cont+1
        wiN = 0

# --------------------------------------------------------------------------------------------
# obtiene ki para gramschmidt()


def KI(vi, qj, m):
    # declaracion de variables
    ki = 0
    vector1 = vector2 = 0

    # obtiene el numerador
    for i in range(0, m):
        vector1 = vector1 + (vi[i]*qj[i])

    # obtiene el denominador
    for i in range(0, m):
        vector2 = vector2 + (qj[i]*qj[i])

    ki = vector1/vector2
    return ki


# --------------------------------------------------------------------------------------------
# sacar transpuesta
def MatrizTraspuesta(matrizTranspuesta, m, matrizBase):
    for i in range(0, m):
        for j in range(0, m):
            matrizTranspuesta[i][j] = matrizBase[j][i]


# --------------------------------------------------------------------------------------------
# multiplica dos matrices
def multiMatrices(matrizIzquierda, matrizDerecha, matrizResultado, m):
    for i in range(0, m):
        for j in range(0, m):
            for k in range(0, m):
                matrizResultado[i][j] += matrizIzquierda[i][k] * \
                    matrizDerecha[k][j]


# --------------------------------------------------------------------------------------------
# comprobar si la matriz que le damos es triangular superior
def TriangularSuperior(matriz, m):
    cont = 0
    respuesta = True
    i = 1
    while(i < m):
        cont = cont+1
        for j in range(0, cont):
            if(matriz[i][j] > 0.0000000001 or matriz[i][j] < -0.0000000001):
                respuesta = False
                break
        i = i+1
    return respuesta

# --------------------------------------------------------------------------------------------
# Algoritmo QR


def algoritmoQR(matriz, matrizOrtogonal, matrizQ, matrizQt, m):
    # definir variables
    triangular = False
    contador = 0
    contador2 = 1
    matrizR = []
    IteraciónMatrizAi = []
    matrizAi = []
    for i in range(0, m):
        matrizAi.append([0]*m)
        IteraciónMatrizAi.append([0]*m)
        matrizR.append([0]*m)

    # paso 1: Ai=A
    for j in range(0, m):
        for i in range(0, m):
            matrizAi[i][j] = matriz[i][j]

    # paso 2: Repetir hasta que Ai+1 sea triangular
    while(triangular == False):
        # a)
        # sacar la matrizQ de Ai
        gramschmidt(matrizAi, m, matrizOrtogonal)
        ortonormalizar(matrizOrtogonal, m, matrizQ)

        # sacar la matrizR de Ai
        MatrizTraspuesta(matrizQt, m, matrizQ)
        multiMatrices(matrizQt, matrizAi, matrizR, m)

        # Imprimir las matrices por cada iteración
        print("\nIteración "+str(contador2))
        print("A"+str(contador)+" = Q"+str(contador)+"*R"+str(contador))

        for i in range(0, m):
            for j in range(0, m):
                print("(", round(matrizAi[i][j], 5), end=")")
            if(i == 0):
                print(" = ", end='')
            else:
                print("   ", end='')
            for j in range(0, m):
                print("(", round(matrizQ[i][j], 5), end=")")
            if(i == 0):
                print("*", end='')
            else:
                print(" ", end='')
            for j in range(0, m):
                print("(", round(matrizR[i][j], 5), end=")")
            print()

        # b)
        # sacar Ai_mas_1 = (R)(Q)
        multiMatrices(matrizR, matrizQ, IteraciónMatrizAi, m)

        # comprobamos si la IteraciónMatrizAi ya es diagonal superior
        if(TriangularSuperior(IteraciónMatrizAi, m) == True):
            # Si es triangular superior, sale del programa
            triangular = True  # si es así, salimos de las iteraciones
        else:  # si no lo es, repetimos el proceso, IteraciónMatrizAi es la nueva Ai
            contador = contador+1
            contador2 = contador2+1
            for j in range(0, m):
                for i in range(0, m):
                    matrizAi[i][j] = IteraciónMatrizAi[i][j]
                    # limpiamos la matriz resultado y la IteraciónMatrizAi para una nueva entrada
                    # Y los datos no se empalmen
                    matrizR[i][j] = 0
                    IteraciónMatrizAi[i][j] = 0

    # los eigenvalores son los elementos de la diagonal
    for i in range(0, m):
        print(str(IteraciónMatrizAi[i][i]), ", ")


# ------------------------------------------------------main------------------------------------------------------
# declaracion de variables
filas = columnas = 0

print("Se creará una matriz de mxm")
m = int(input("ingrese el valor de m: "))
filas = columnas = m

matriz = []
matrizOrtonormal = []
matrizOrtogonal = []
matrizTranspuesta = []
matrizResultado = []
for i in range(filas):
    matriz.append([0]*columnas)
    matrizOrtogonal.append([0]*columnas)
    matrizOrtonormal.append([0]*columnas)
    matrizTranspuesta.append([0]*columnas)
    matrizResultado.append([0]*columnas)

# rellenar la matriz con vectores acomodados por columnas, no por filas
for j in range(columnas):
    for i in range(filas):
        matriz[i][j] = float(
            input("Ingresa la componente del vector "+str(j)+" posicion "+str(i)+": "))

# Imprimir matriz original
print("\nMatriz: ")
imprMatriz(matriz, m)

print("\nBase ortogonal por GramSchmidt: ")
gramschmidt(matriz, m, matrizOrtogonal)
imprMatriz(matrizOrtogonal, m)

print("\nBase ortonormal: ")
ortonormalizar(matrizOrtogonal, m, matrizOrtonormal)
imprMatriz(matrizOrtonormal, m)

print("\n eigenvalores: ")
algoritmoQR(matriz, matrizOrtogonal, matrizOrtonormal, matrizTranspuesta, m)

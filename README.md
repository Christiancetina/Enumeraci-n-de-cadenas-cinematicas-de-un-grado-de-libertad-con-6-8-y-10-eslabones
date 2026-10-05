# Documentación del Código: Síntesis Estructural y Enumeración de Inversiones

Este documento proporciona una explicación detallada del programa en Mathematica diseñado para sintetizar cadenas cinemáticas y enumerar sus inversiones (mecanismos) a partir de un número dado de eslabones ($N$) y juntas ($J$). 

El código integra cinemática de mecanismos, álgebra lineal, combinatoria y teoría de grafos. A continuación, se explica cada bloque funcional.

---

## 1. Inicialización y Conceptos Básicos de Movilidad

```mathematica
ClearAll["Global`*"];

KinematicSynthesis[Nlinks_Integer, Jjoints_Integer] := Module[
  {F, qMax, eq1, eq2, vars, families, combSpace, ...},
  
  F = 3*(Nlinks - 1) - 2*Jjoints;
  combSpace = Binomial[Nlinks*(Nlinks - 1)/2, Jjoints];
  
  If[F < 0, Print["[!] Error: Grados de libertad negativos."]; Return[]];
```

### Conceptos Involucrados:
*   **Alcance y Variables (`Module`):** En Mathematica, `Module` crea un entorno local para las variables, evitando que interfieran con otras partes del sistema.
*   **Criterio de Grübler ($F$):** Es la ecuación fundamental de la cinemática plana: $F = 3(N-1) - 2J$. Calcula los grados de libertad de un mecanismo. Si $F < 0$, la estructura está "sobrerrestringida" (es una armadura precargada), por lo que el programa se detiene usando `Return[]`.
*   **Espacio Combinatorio (`Binomial`):** Un grafo completo de $N$ nodos tiene un máximo de $\frac{N(N-1)}{2}$ aristas posibles. La función `Binomial[n, k]` calcula el coeficiente binomial (combinaciones), indicando cuántas formas matemáticas existen de ubicar $J$ juntas en esos espacios posibles. Para $N=8, J=10$, esto supera los 13 millones.

---

## 2. Cálculo de Familias Estructurales (Ecuaciones Diofánticas)

```mathematica
  qMax = If[F <= 1, Floor[(Nlinks - F + 1)/2], Min[Nlinks - F - 1, Floor[(Nlinks + F - 1)/2]]];
  vars = Table[Symbol["n" <> ToString[i]], {i, 2, qMax}];
  eq1 = Total[vars] == Nlinks;
  eq2 = Sum[i * vars[[i - 1]], {i, 2, qMax}] == 2*Jjoints;
  
  families = Solve[{eq1, eq2, And @@ Thread[vars >= 0]}, vars, Integers];
```

### Conceptos Involucrados:
*   **Orden Máximo Poligonal ($qMax$):** Un eslabón puede ser binario (2 nodos), ternario (3 nodos), cuaternario, etc. Existe un límite físico de cuántas juntas puede tener un solo eslabón dependiendo de $N$ y $F$.
*   **Ecuaciones Diofánticas:** Son ecuaciones polinómicas donde las soluciones deben ser estrictamente números enteros (no existe "medio eslabón").
    *   **Ecuación de Nodos (`eq1`):** $n_2 + n_3 + \dots + n_q = N$. La suma de la cantidad de cada tipo de eslabón debe dar el total de eslabones.
    *   **Ecuación de Aristas (`eq2`):** $2n_2 + 3n_3 + \dots + qn_q = 2J$. Como cada junta une a dos eslabones, el total de terminales disponibles en los eslabones debe ser el doble de las juntas.
*   **`Solve` e `Integers`:** Mathematica resuelve el sistema restringiendo la búsqueda a números enteros positivos (`vars >= 0`), arrojando las **Familias Estructurales** (ej. 4 binarios y 4 ternarios).

---

## 3. Pre-filtro Algebraico: Polinomio de Bôcher

```mathematica
  BocherPolynomial[adj_?MatrixQ] := Module[{n, s, a, j, r},
    n = Length[adj];
    s = Table[Tr[MatrixPower[adj, r]], {r, 1, n}];
    a = ConstantArray[0, n + 1];
    a[[1]] = 1;
    For[j = 1, j <= n, j++,
      a[[j + 1]] = -1/j * Sum[a[[j - r + 1]] * s[[r]], {r, 1, j}];
    ];
    a
  ];
```

### Conceptos Involucrados:
*   **Matriz de Adyacencia:** Representación matricial del grafo. Si el eslabón $i$ se une al $j$, la posición $(i,j)$ es $1$; si no, es $0$.
*   **Trazas y Potencias (`Tr`, `MatrixPower`):** La traza de una matriz es la suma de su diagonal principal. Elevar la matriz de adyacencia a la potencia $r$ y extraer su traza $s_r = Tr(A^r)$ revela el número de trayectorias cerradas de longitud $r$ en el grafo.
*   **Fórmula de Bôcher:** Tradicionalmente, calcular el polinomio característico requiere el determinante $|\lambda I - A|$. Bôcher proporciona una fórmula iterativa para hallar los coeficientes ($a_j$) usando únicamente las trazas.
*   **Huella Topológica:** Dos grafos isomórficos (idénticos en forma pero numerados distinto) *siempre* tienen el mismo polinomio. Comparar polinomios es computacionalmente baratísimo frente a la comparación directa de grafos ($O(V!)$).

---

## 4. Filtro Estricto de Subcadenas Rígidas (Densidad Topológica)

```mathematica
  FilterRigidSubchains[adj_?MatrixQ] := Catch[
    Module[{ePrime},
      If[VertexConnectivity[AdjacencyGraph[adj]] < 2, Throw[False]];
      Do[
        Do[
          ePrime = Total[adj[[sv, sv]], 2] / 2;
          If[2 * ePrime > 3 * k - 4, Throw[False]];
        , {sv, Subsets[Range[Nlinks], {k}]}];
      , {k, 3, Nlinks - 1}];
      True
    ]
  ];
```

### Conceptos Involucrados:
*   **Conectividad de Vértices (`VertexConnectivity`):** Un mecanismo no puede tener un eslabón que, al ser retirado, divida la máquina en dos partes sueltas. La conectividad debe ser $\ge 2$.
*   **Subgrafos y Subconjuntos (`Subsets`):** Se evalúan todos los agrupamientos posibles de $k$ eslabones dentro de la matriz (`adj[[sv, sv]]`).
*   **Desigualdad de Densidad ($2e' \le 3v' - 4$):** Concepto avanzado. Para que ninguna parte del mecanismo se bloquee formando una armadura rígida, cualquier subgrupo de vértices ($v' = k$) y sus aristas internas ($e'$) deben tener movilidad $> 0$. Al despejar la ecuación de Grübler local, se obtiene esta inecuación. Si un subgrafo la viola, se descarta todo el mecanismo usando `Throw[False]`.

---

## 5. Muestreo Estocástico y Secuencias de Grados

```mathematica
  GetDegreeSequences[famVals_] := Flatten[
    Table[ConstantArray[i + 1, famVals[[i]]], {i, 1, Length[famVals]}]
  ];
...
      Quiet@Do[
        g = RandomGraph[DegreeGraphDistribution[degSeq]];
        If[ConnectedGraphQ[g],
          adj = Normal[AdjacencyMatrix[g]];
```

### Conceptos Involucrados:
*   **Secuencia de Grados (`degSeq`):** Transforma una familia (ej. 2 binarios, 1 ternario) en una lista plana de conexiones requeridas por nodo: `{2, 2, 3}`.
*   **Distribución de Grados (`DegreeGraphDistribution`):** En lugar de generar millones de matrices al azar, se fuerza a Mathematica a generar grafos que obedezcan estrictamente la secuencia de grados calculada.
*   **Muestreo de Saturación (`RandomGraph` en un bucle):** Método de Monte Carlo o fuerza bruta inteligente. Se generan miles de grafos (5000 iteraciones). Debido a que el espacio de familias ya redujo las posibilidades, saturar el sistema con 5000 intentos asegura encontrar todas las soluciones únicas posibles sin desbordar la memoria.

---

## 6. Depuración Final de Isomorfismos

```mathematica
            If[!MemberQ[polyHashes, pHash],
               AppendTo[polyHashes, pHash];
               AppendTo[validGraphsFam, g];
            , 
               If[Count[validGraphsFam, x_ /; IsomorphicGraphQ[x, g]] == 0,
                 AppendTo[polyHashes, pHash];
                 AppendTo[validGraphsFam, g];
               ];
            ];
```

### Conceptos Involucrados:
*   **Tablas Hash (`polyHashes`):** Guarda los polinomios encontrados. Si un polinomio no está en la lista (`!MemberQ`), el grafo es 100% único.
*   **Colisiones Polinomiales:** Existen grafos raros (grafos coespectrales) que, siendo distintos, comparten el mismo polinomio de Bôcher. 
*   **Evaluación Profunda (`IsomorphicGraphQ`):** Si y solo si dos grafos comparten polinomio, se invoca esta función (que usa el algoritmo Nauty por debajo) para verificar matemáticamente si sus nodos pueden mapearse 1 a 1. Es el filtro de seguridad definitivo.

---

## 7. Enumeración de Inversiones Cinemáticas (Teoría de Grupos)

```mathematica
    Module[{g, orbits, numInversions},
      g = validChains[[k]];
      orbits = GroupOrbits[GraphAutomorphismGroup[g], VertexList[g]];
      numInversions = Length[orbits];
      totalMechanisms += numInversions;
```

### Conceptos Involucrados:
*   **Inversión Cinemática:** Es el acto de seleccionar qué eslabón de la cadena se fijará a la tierra (bancada) para crear un mecanismo operativo.
*   **Grupo de Automorfismos (`GraphAutomorphismGroup`):** Un automorfismo es una permutación de los vértices que deja al grafo viéndose exactamente igual (simetría topológica). El conjunto de todas estas permutaciones forma un grupo algebraico.
*   **Órbitas (`GroupOrbits`):** En teoría de grupos, si una simetría puede mover el Nodo A hacia la posición del Nodo B, ambos pertenecen a la misma "órbita" (clase de equivalencia).
*   **Mecanismos Únicos:** Fijar a tierra cualquier eslabón de la misma órbita producirá exactamente el mismo mecanismo físico. Por tanto, el número de mecanismos no isomórficos derivados de una cadena es exactamente igual a la cantidad de órbitas (`Length[orbits]`)..

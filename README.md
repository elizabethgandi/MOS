# MOS

## Point de départ
Générateurs créés en utilisant la méthode de la contrainte-ε et la recherche dichotomique sur le problème relaxé (relaxation linéaire).

- $I \in \{1,\dots, |I|\}$, l'ensemble des variables.
- $J \in \{1,\dots, |J|\}$, l'ensemble des contraintes.
- $A_{ij}, ~ ~ \forall i \in I, \forall j \in J$, la matrice des contraintes.
- $x_{i} \in \{0, 1\}, ~ ~ \forall i \in I$, les variables binaires.

## Première idée
À partir de chaque générateur, fixer chaque variable qui était égale à un (resp. zéro) à un (resp. zéro) et résoudre le modèle en recherchant une solution entière. De cette manière, le nombre total de variables sera grandement réduit, ne laissant comme variables que celles avec des valeurs flottantes dans le générateur.

Les résultats ont montré que cette méthode est inutile, car fixer la plupart des variables rend le modèle infaisable dans la majorité des cas (si ce n'est dans tous les cas).

## Première idée bis
À partir du générateur de la dichotomie, fixer les variables fractionnaires en variables binaires et de résoudre le problème avec un modèle d'optimisation en nombre entier mixte (MILP). Parallèlement, les variables déjà déterminées à 0 ou 1 dans le générateur sont résolues de manière continue.

## Deuxième idée
À partir de chaque générateur, calculer une solution entière qui minimise uniquement la différence entre le générateur flottant et une nouvelle solution entière calculée.

Pour cela, nous introduisons $y$, une variable binaire pour déterminer si une variable a changé de valeur par rapport au générateur.
- $I' \subset I$, inclut chaque indice de variable où le générateur a pour valeur un ou zéro.
- $y_{i} \in \{0, 1\}, ~ ~ \forall i \in I'$, est égal à $1$ si la valeur de la solution pour la variable d'indice $i$ est différente de celle du générateur. Sinon, $0$.

Pour imposer le comportement de $y$, nous ajoutons deux types de contraintes au modèle. Pour chaque variable qui était égale à un dans le générateur, une pénalité est appliquée si elle est changée à zéro :
$x_i \geq 1 - y_i$.
Réciproquement, les variables fixées à zéro dans le générateur pénaliseront la fonction si elles sont changées à un :
$x_i \leq y_i$.

La dernière étape consiste à calculer à quel point le changement d'une variable pénalisera la solution. Dans un premier temps, nous testerons différentes fonctions objectives.

- **Pénalité d'une unité (ones)** : changer une variable est pénalisé par un :
  $$\min \sum_{i \in I} y_i$$
- **Somme pondérée (wsum)** : la pénalité sera proportionnelle à chaque fonction objective :
  $$\min \sum_{j \in J} y_j \cdot \frac{c_{1j} + c_{2j}}{\max_{k \in I}(c_{1k}) + \max_{k \in I}(c_{2k})}$$
- **Intérêt SPA (wspa)** : la pénalité dépendra de l'importance de la variable dans le contexte d'un problème SPA :
  $$\min \sum_{i \in I} y_i \cdot \frac{c_{1i} + c_{2i}}{\sum_{j \in J} A_{ij}}$$

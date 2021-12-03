# TP - Lissage de maillages - Modélisations géométriques

Charly Finette

L'essentiel du code se trouve en bas du fichier "mainwindow.cpp".

## Implémentation de l'opérateur de Laplace-Beltrami

``` c++
// Opérateur de Laplace-Beltrami
// Etant donné un sommet v, renvoit le vecteur res que l'on applique au sommet v 
// pour effectuer un lissage à l'aide de l'opérateur de Laplace Beltrami

MyMesh::Point MainWindow::delta(MyMesh::VertexHandle vh)
```
Je ne vais pas détaillé plus précisément cette implémentation étant donné que:
- Je vous l'avais déjà montré durant le TP.
- Tout est réimplémenté dans l'exercice 2 de manière plus modulaire et donc plus claire.

``` c++
// Opérateur uniforme
MyMesh::Point MainWindow::deltaUniforme(MyMesh::VertexHandle vh)
```

Pour appliquer ces deux opérateurs au maillage, j'utilise :
``` c++
// Application des opérateurs sur le maillage
void MainWindow::lissage()

{
    // Pour chaque sommet du maillage
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it )
    {
       // Ici on peut choisir quel opérateur utiliser
       mesh.point(*v_it) += deltaUniforme(*v_it);
       // mesh.point(*v_it) += delta(*v_it);
    }
}
```
Cette fonction est accesible par le bouton " lissage uniforme " de mon programme. Si vous modifier le code pour voir le fonctionnement de mon opérateur de Laplace-Beltrami, il faudra quand même utiliser le bouton "lissage uniforme".
J'ai choisi de donner un bouton à l'opérateur uniforme parce qu'il fonctionne mieux.

## Implémentation des matrices de Laplace Beltrami

Je vais utiliser la librairie Eigen pour le calcul matriciel.

La méthode .idx() :permet de retourner l'id associé à chaque vertexhandle (ils sont numérotés de 0 à nbr_vertices - 1)

- implémentation différente du calcul du coef cotan : utilisation de la structure de half edge de OpenMesh

## Flou de diffusion sur un maillage bruité

- fonctionne de manière "correcte" sur la sphère mais bug sur les plus gros maillages

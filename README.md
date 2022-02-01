# TP - Lissage de maillages - Modélisations géométriques

Charly Finette

## Outils 

## Installation 

## Problématique

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

La méthode .idx() permet de retourner l'id associé à chaque vertexhandle. (ils sont numérotés de 0 à nbr_vertices - 1)

Fonction pour calculer la valeur A(v) pour v un sommet du maillage:
``` c++
// Calcul de A en fonction des aires des faces voisines du sommet
float MainWindow::A(MyMesh::VertexHandle vh)
{
    float A = 0;

    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
    {
        VertexHandle vh_neigh = *vv_it;
        HalfedgeHandle heh = mesh.find_halfedge(vh_neigh, vh);
        HalfedgeHandle heh_suiv = mesh.next_halfedge_handle(heh);
        VertexHandle vh_neigh_suiv = mesh.to_vertex_handle(heh_suiv);

        MyMesh::Point vec_0 = mesh.point(vh_neigh) - mesh.point(vh);
        MyMesh::Point vec_1 = mesh.point(vh_neigh_suiv) - mesh.point(vh);

       A += cross(vec_0, vec_1).norm() / 6;
    }

     return A;
}
```

J'utilise la formule pour calculer l'aire à partir du cross product. Comme cette méthode faussent les résultats (A est tellement petit que 1/ 2*A explose)
j'utilise un A plus simple : On divise seulement par le nombre de voisins. ( 1 / ( 2 * A / 2) = 1 / A ) 
Cette méthode transforme la sphère en une "pomme de pin".

``` c++
//Calcul de A en fonction du nombre de voisin
float MainWindow::A_alternatif(MyMesh::VertexHandle vh)
{
    int nbr_voisins = 0;

    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
    {
        nbr_voisins++;
    }

    float A = nbr_voisins / 2;

     return A;
}
``` 

J'ai recodé la partie suivante en utilisant la structure de halfedge, central à OpenMesh :
``` c++
// Etant donné un sommet v, et d'un voisin v_neigh de v, retourne la valeur du coefficient cotangenciel associé
float MainWindow::poids_cotangenciel_voisin(MyMesh::VertexHandle vh, MyMesh::VertexHandle vh_neigh)
{
    MyMesh::HalfedgeHandle heh = mesh.find_halfedge(vh, vh_neigh);

    MyMesh::HalfedgeHandle heh_prec = mesh.next_halfedge_handle(heh);
    MyMesh::VertexHandle   vh_prec = mesh.to_vertex_handle(heh_prec);

    MyMesh::HalfedgeHandle heh_op = mesh.opposite_halfedge_handle(heh);
    MyMesh::HalfedgeHandle heh_i = mesh.next_halfedge_handle(heh_op);
    MyMesh::VertexHandle vh_suiv = mesh.to_vertex_handle(heh_i);

    // On calcule l'angle alpha
    MyMesh::Point vec_alpha_0 = mesh.point(vh) - mesh.point(vh_prec);
    MyMesh::Point vec_alpha_1 = mesh.point(vh_neigh) - mesh.point(vh_prec);

    float cos_alpha = dot(vec_alpha_0,vec_alpha_1) / (norm(vec_alpha_0) * norm(vec_alpha_1));
    float cotan_alpha = cos_alpha/(sqrt(1.0f - cos_alpha*cos_alpha));

    // On calcule l'angle beta

    MyMesh::Point vec_beta_0 = mesh.point(vh) - mesh.point(vh_suiv);
    MyMesh::Point vec_beta_1 = mesh.point(vh_neigh) - mesh.point(vh_suiv);
    //qDebug() << "vec_beta_0    " << vec_beta_0[0] << " " << vec_beta_0[1] << " " << vec_beta_0[2];
    //qDebug() << "vec_beta_1    " << vec_beta_1[0] << " " << vec_beta_1[1] << " " << vec_beta_1[2];

    float cos_beta = dot(vec_beta_0, vec_beta_1) / (norm(vec_beta_0) * norm(vec_beta_1));
    float cotan_beta = cos_beta/(sqrt(1.0f - cos_beta*cos_beta));

    //qDebug() << "cos_beta = " << cos_beta;


    return cotan_alpha + cotan_beta;
}
```
Néanmoins, les résultats sont similaires : Le problème provient probablement de A. 

Implémentation de la matrice D.

``` c++
Eigen::MatrixXd MainWindow::Matrix_D()
{
    const int nbr_vertices = mesh.n_vertices();
    Eigen::MatrixXd D(nbr_vertices, nbr_vertices);

    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        for(MyMesh::VertexIter v_it2 = mesh.vertices_begin(); v_it2 != mesh.vertices_end(); ++v_it2)
    {
        MyMesh::VertexHandle vh = *v_it;
        MyMesh::VertexHandle vh2 = *v_it2;
        int i = (vh).idx();
        int j = (vh2).idx();

        // On initialise la matrice à zero
        D(i,j) = 0;

        if(i == j)
        {
            // Ici on peut choisir le A qu'on veut utiliser
            //D(i,j) =  A(vh);
            //D(i,j) = 1 / ( 2 * A(vh) );
            D(i,j) = 1 / ( 2 * A_alternatif(vh) );
        }
   }

    return D;
}
```

Implémentation de la matrice M.
``` c++
Eigen::MatrixXd MainWindow::Matrix_M()
{
    const int nbr_vertices = mesh.n_vertices();
    Eigen::MatrixXd M(nbr_vertices, nbr_vertices);

    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        for(MyMesh::VertexIter v_it2 = mesh.vertices_begin(); v_it2 != mesh.vertices_end(); ++v_it2)
    {
        MyMesh::VertexHandle vh = *v_it;
        MyMesh::VertexHandle vh2 = *v_it2;
        int i = (vh).idx();
        int j = (vh2).idx();

        M(i,j) = 0;

        if(i == j)
        {
            M(i,j) = - poids_cotangenciel_voisins(vh);
        }

        else
        {
            for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); vv_it++) // on itère sur les voisins de i
            {
                if( (*vv_it).idx() == j) // La condition vérifie que j est dans le voisinage de i
                {
                    M(i,j) = poids_cotangenciel_voisin(vh, vh2);
                    break;
                }

            }
        }
    }

    return M;

}
```


## Flou de diffusion sur un maillage bruité

Fonctionne de manière "correcte" sur la sphère mais "crash" sur les plus gros maillages.

``` c++
/ Flou de diffusion : application des matrices sur le maillage
void MainWindow::FlouDeDiffusion(float h, float lambda)
{
    const int nbr_vertices = mesh.n_vertices();
    Eigen::VectorXd X(nbr_vertices); // vecteur contenant la coordonnée x de tous les sommets
    Eigen::VectorXd Y(nbr_vertices); // vecteur contenant la coordonnée y de tous les sommets
    Eigen::VectorXd Z(nbr_vertices); // vecteur contenant la coordonnée z de tous les sommets

    //Initialisation de X,Y et Z
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        X((*v_it).idx()) = mesh.point(*v_it)[0];
        Y((*v_it).idx()) = mesh.point(*v_it)[1];
        Z((*v_it).idx()) = mesh.point(*v_it)[2];
    }

    Eigen::MatrixXd L = Matrix_D() * Matrix_M();

    X = X + lambda * h * L * X;
    Y = Y + lambda * h * L * Y;
    Z = Z + lambda * h * L * Z;

    // Application des résultats sur le mesh
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        int i = (*v_it).idx();

        mesh.point(*v_it)[0] = X[i];
        mesh.point(*v_it)[1] = Y[i];
        mesh.point(*v_it)[2] = Z[i];
    }
}
```




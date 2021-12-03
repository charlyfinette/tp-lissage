#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <Eigen/Core>
#include <Eigen/Geometry>




/* **** début de la partie boutons et IHM **** */


// exemple pour charger un fichier .obj ou .ply

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.ply *.obj)"));
    qDebug() << fileName.toUtf8().constData();
    // chargement du fichier dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_lissage_clicked()
{
    FlouDeDiffusion(1,1);
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_lissage_uniforme_clicked()
{
    lissage();
    displayMesh(&mesh);
}


void MainWindow::on_pushButton_enregistrer_clicked()
{

    if (!OpenMesh::IO::write_mesh(mesh, "/home/charly/Documents/modeles-geom/TP_LissageDeMaillage/out.ply"))
    {
      qDebug() << "write error\n";
    }
}



/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


// ********************************************************************************************************************************
// ****************************************** LISSAGE *****************************************************************************
// ********************************************************************************************************************************

// Exercice 1 : implémentation de l'opérateur de Laplace-Beltrami et de l'opérateur uniforme

// Opérateur uniforme
MyMesh::Point MainWindow::deltaUniforme(MyMesh::VertexHandle vh)
{
    MyMesh::Point res;
    int nbr_voisins = 0;

    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
    {
        res += mesh.point(*vv_it) - mesh.point(vh);
        nbr_voisins++;
    }

    res/=nbr_voisins;

    qDebug() << res[0] << res[1] << res[2] ;
    return res;
}

// Opérateur de Laplace-Beltrami
// Etant donné un sommet v, renvoit le vecteur  res avec lequel déplacer le sommet v pour effectuer un lissage à l'aide de l'opérateur de Laplace Beltrami
MyMesh::Point MainWindow::delta(MyMesh::VertexHandle vh)
{
    // On recupere ses voisins
    std::vector<MyMesh::VertexHandle> voisins;

    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
    {
         //qDebug() << "v    " << mesh.point(*vv_it)[0] << " " << mesh.point(*vv_it)[1] << " " << mesh.point(*vv_it)[2];
        voisins.push_back(*vv_it);
    }

    //qDebug() << "vh    " << mesh.point(vh)[0] << " " << mesh.point(vh)[1] << " " << mesh.point(vh)[2];

    int nbr_voisins = voisins.size();

    // On initialise le résultat à 0
    MyMesh::Point res;
    res[0] = 0.0f;
    res[1] = 0.0f;
    res[2] = 0.0f;

    // Pour chaque voisin v_i

    for(int i = 0; i < nbr_voisins ; i++)
    {
        // On calcule le cotan de alpha_i
        MyMesh::Point vec_alpha_0;
        MyMesh::Point vec_alpha_1;

        vec_alpha_0 = mesh.point(vh) - mesh.point(voisins[(i + 1) % nbr_voisins]);
        vec_alpha_1 = mesh.point(voisins[i]) - mesh.point(voisins[(i + 1) % nbr_voisins]);

        float cos_alpha_curr = dot(vec_alpha_0,vec_alpha_1) / norm(vec_alpha_0)*norm(vec_alpha_1);
        float cotan_alpha_curr = cos_alpha_curr/(sqrt(1.0f - cos_alpha_curr*cos_alpha_curr));


        // On calcule le cotan de beta_i
        MyMesh::Point vec_beta_0;
        MyMesh::Point vec_beta_1;

        vec_beta_0 = mesh.point(vh) - mesh.point(voisins[(i - 1 + nbr_voisins) % nbr_voisins]);
        vec_beta_1 = mesh.point(voisins[i]) - mesh.point(voisins[(i - 1 + nbr_voisins) % nbr_voisins]);


        float cos_beta_curr = dot(vec_beta_0,vec_beta_1) / norm(vec_beta_0)*norm(vec_beta_1);
        float cotan_beta_curr = cos_beta_curr/(sqrt(1.0f - cos_beta_curr*cos_beta_curr));


        // On calcule cotan alpha_i + cotan beta_i
        float coef_curr = cotan_alpha_curr + cotan_beta_curr;

        // On fait des tests pour vérifier que l'on a pas d'erreurs de calculs avec les cotans.
        if(std::isnan(coef_curr))
            coef_curr = 0;
        const float eps = 1e-6f;
        const float cotan_max = cos( eps ) / sin( eps );
        if(coef_curr >= cotan_max)
            coef_curr = cotan_max;

        // On calcule un élément de la somme
        MyMesh::Point point_curr = coef_curr * (mesh.point(voisins[i]) - mesh.point(vh));

        // Puis on l'ajoute à la somme
        res += point_curr;
        // qDebug() << "point_curr     " << point_curr[0] << " " << point_curr[1]  << " " << point_curr[2];

    }

   // On calcule le coefficient A
    float A = 0;

     for(int i = 0; i < nbr_voisins; i++)
     {
         MyMesh::Point vec_0 = mesh.point(vh) - mesh.point(voisins[i]);
         MyMesh::Point vec_1 = mesh.point(vh) - mesh.point(voisins[(i + 1) % nbr_voisins]);

        A += cross(vec_0, vec_1).norm() / 2;
     }

    qDebug() << "res            " << res[0] << " " << res[1]  << " " << res[2];

// Formule suggéré par mes camarades, extrait du code source de Blender
    return res * (1 / (4 * A));
}


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

// Exercice 2 : Implémentation des matrices de Laplace Beltrami

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

// Applique la fonction précédente a tous les voisins d'un sommet
float MainWindow::poids_cotangenciel_voisins(MyMesh::VertexHandle vh)
{
    float w = 0;
    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); vv_it++)
    {
        w += poids_cotangenciel_voisin(vh, *vv_it);
    }

    return w;
}

// Matrices de Laplace-Beltrami

//Implémentation de la matrice D
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

//Implémentation de la matrice M
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

// Flou de diffusion : application des matrices sur le maillage
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


#include "mainwindow.h"
#include "ui_mainwindow.h"




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
    qDebug() << "res            " << res[0] << " " << res[1]  << " " << res[2];
    return res;
}



// Opérateur de Laplace-Beltrami
float MainWindow::A(MyMesh::VertexHandle vh)
{
    std::vector<MyMesh::VertexHandle> voisins;

    for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
    {
        voisins.push_back(*vv_it);
    }

    int nbr_voisins = voisins.size();
    float A = 0;

     for(int i = 0; i < nbr_voisins; i++)
     {
         MyMesh::Point vec_0 = mesh.point(vh) - mesh.point(voisins[i]);
         MyMesh::Point vec_1 = mesh.point(vh) - mesh.point(voisins[(i + 1) % nbr_voisins]);

        A += cross(vec_0, vec_1).norm() / 6;
     }
}

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

    MyMesh::Point res;
    res[0] = 0.0f;
    res[1] = 0.0f;
    res[2] = 0.0f;

    // Pour chaque voisin v_i

    for(int i = 0; i < nbr_voisins 0; i++)
    {
        // On calcule alpha_i
        MyMesh::Point vec_alpha_0;
        MyMesh::Point vec_alpha_1;

        vec_alpha_0 = mesh.point(vh) - mesh.point(voisins[(i + 1) % nbr_voisins]);
        vec_alpha_1 = mesh.point(voisins[i]) - mesh.point(voisins[(i + 1) % nbr_voisins]);

        float cos_alpha_curr = dot(vec_alpha_0,vec_alpha_1) / norm(vec_alpha_0)*norm(vec_alpha_1);
        float cotan_alpha_curr = cos_alpha_curr/(sqrt(1.0f - cos_alpha_curr*cos_alpha_curr));


        // On calcule beta_i
        MyMesh::Point vec_beta_0;
        MyMesh::Point vec_beta_1;

        vec_beta_0 = mesh.point(vh) - mesh.point(voisins[(i - 1 + nbr_voisins) % nbr_voisins]);
        vec_beta_1 = mesh.point(voisins[i]) - mesh.point(voisins[(i - 1 + nbr_voisins) % nbr_voisins]);


        float cos_beta_curr = dot(vec_beta_0,vec_beta_1) / norm(vec_beta_0)*norm(vec_beta_1);
        float cotan_beta_curr = cos_beta_curr/(sqrt(1.0f - cos_beta_curr*cos_beta_curr));

        // On fait des tests pour vérifier que l'on a pas d'erreurs de calculs avec les cotans.
        float coef_curr = cotan_alpha_curr + cotan_beta_curr;
        //qDebug() << "dot vec_alpha    " << dot(vec_alpha_0,vec_alpha_1) ;
        //qDebug() << "dot vec_beta     " << dot(vec_beta_0,vec_beta_1) ;
        //qDebug() << "cos_beta_curr     " << cos_beta_curr;
        //qDebug() << "cos_alpha_curr     " << cos_alpha_curr;
        //qDebug() << "coef_curr     " << coef_curr;
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
    return res * (1 / (4 * A));
}



void MainWindow::lissage()
{   
    std::vector<MyMesh::Point> mem;

    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
       // qDebug() << "coord     " << mesh.point(*v_it)[0] << " " << mesh.point(*v_it)[1]  << " " << mesh.point(*v_it)[2];
       mem.push_back(delta(*v_it));
       //qDebug() << "delta             "<< delta2(*v_it)[0] << " " << delta2(*v_it)[1] << " " << delta2(*v_it)[2];
       //qDebug() << "deltaUniforme     "<< deltaUniforme(*v_it)[0] << " " << deltaUniforme(*v_it)[1] << " " << deltaUniforme(*v_it)[2];
    }

    // Pour chaque sommet du maillage
    int c = 0 ;
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it )
    {
       mesh.point(*v_it) += mem[c];
       c++;
    }

}

// Matrice de Laplace-Beltrami
std::vector<std::vector<double>> MainWindow::Matrix_D()
{
    std::vector<std::vector<double>> D;

    std::vector<VertexHandle> vertexHandles;
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it )
    {
      vertexHandles.push_back(*v_it);
    }

    const int nbr_vertices = mesh.n_vertices();

    for(int i = 0; i < nbr_vertices; i++)
    for(int j = 0; j < nbr_vertices; j++)
    {
        if(i == j)
        {
            D[i][j] = A(vertexHandles[i]);
        }
        else
        {
            D[i][j] = 0;
        }
    }

    return D;
}



std::vector<std::vector<double>> MainWindow::Matrix_M()
{
    std::vector<std::vector<double>> M;

    std::vector<VertexHandle> vertexHandles;
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it )
    {
      vertexHandles.push_back(*v_it);
    }

    const int nbr_vertices = mesh.n_vertices();

    for(int i = 0; i < nbr_vertices; i++)
    {
        std::vector<MyMesh::VertexHandle> voisins;
        for(MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vertexHandles[i]); vv_it.is_valid(); ++vv_it)
        {
            voisins.push_back(*vv_it);
        }

        for(int j = 0; j < nbr_vertices; j++)
        {
            for( const auto & neigh : voisins)
            {

            }
        }
    }

    return M;
}


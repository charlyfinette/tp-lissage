#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/writer/PLYWriter.hh>

#include <Eigen/Core>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    float A(MyMesh::VertexHandle vh);
    float A_alternatif(MyMesh::VertexHandle vh);

    float poids_cotangenciel_voisin(MyMesh::VertexHandle vh, MyMesh::VertexHandle vh_neigh);
    float poids_cotangenciel_voisins(MyMesh::VertexHandle vh);


    MyMesh::Point delta(MyMesh::VertexHandle vh);
    MyMesh::Point deltaUniforme(MyMesh::VertexHandle vh);

    Eigen::MatrixXd Matrix_D();
    Eigen::MatrixXd Matrix_M();

    void FlouDeDiffusion(float h, float lambda);

    void lissage();



private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButton_lissage_clicked();
    void on_pushButton_lissage_uniforme_clicked();
    void on_pushButton_enregistrer_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H

#include "SeismicViewer.h"

#include <vtkSphereSource.h>
#include <vtkGlyph3DMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkScalarBarActor.h>
#include <vtkAxesActor.h>
#include <vtkPCAStatistics.h>
#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include <vtkStatisticsAlgorithm.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkNew.h>

#include <iostream>
#include <cmath>
#include <iomanip>

SeismicViewer::SeismicViewer() {
    m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_window = vtkSmartPointer<vtkRenderWindow>::New();
    m_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    
    m_window->AddRenderer(m_renderer);
    m_window->SetSize(1000, 800);
    m_interactor->SetRenderWindow(m_window);
    m_renderer->SetBackground(0.9, 0.9, 0.95);
}

void SeismicViewer::VisualizeEvents(vtkSmartPointer<vtkPolyData> data) {
    if (!data) return;

    // 1. Setup the main event glyphs (spheres colored by Magnitude)
    SetupEventGlyph(data);

    // 2. Add Bounding Box
    double bounds[6];
    data->GetBounds(bounds);
    AddBoundingBox(bounds);
}

void SeismicViewer::SetupEventGlyph(vtkSmartPointer<vtkPolyData> data) {
    auto sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetRadius(10.0);
    sphere->SetThetaResolution(16);
    sphere->SetPhiResolution(16);

    auto glyph = vtkSmartPointer<vtkGlyph3DMapper>::New();
    glyph->SetInputData(data);
    glyph->SetSourceConnection(sphere->GetOutputPort());
    glyph->SetScalarModeToUsePointFieldData(); 
    glyph->SelectColorArray("Magnitude");
    glyph->SetColorModeToMapScalars();
    glyph->SetScaleFactor(1.0);
    
    // Determine magnitude range
    double magnitudeRange[2] = {-5.0, 1.0}; // Default range
    vtkDataArray* magnitudeArray = data->GetPointData()->GetArray("Magnitude");
    if (magnitudeArray) {
        magnitudeArray->GetRange(magnitudeRange);
        std::cout << "Detected Magnitude Range: [" << magnitudeRange[0] << ", " 
                  << magnitudeRange[1] << "]" << std::endl;
    } else {
        std::cerr << "Error: Magnitude array not found! Using default range." << std::endl;
    }
    glyph->SetScalarRange(magnitudeRange[0], magnitudeRange[1]); // For visualizing later

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(glyph);
    m_renderer->AddActor(actor);
    
    // Add the Scalar Bar and Orientation Widgets
    AddWidgets(glyph);
}

void SeismicViewer::AddBoundingBox(const double bounds[6]) {
    auto cube = vtkSmartPointer<vtkCubeSource>::New();
    cube->SetBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

    auto cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cubeMapper->SetInputConnection(cube->GetOutputPort());

    auto cubeActor = vtkSmartPointer<vtkActor>::New();
    cubeActor->SetMapper(cubeMapper);
    
    cubeActor->GetProperty()->SetColor(0.7, 0.7, 0.7); 
    cubeActor->GetProperty()->SetOpacity(0.1); 
    cubeActor->GetProperty()->SetRepresentationToWireframe(); 

    m_renderer->AddActor(cubeActor);
}

void SeismicViewer::AddWidgets(vtkSmartPointer<vtkGlyph3DMapper> glyph) {
    auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetTitle("Moment Magnitude");
    scalarBar->SetNumberOfLabels(5);
    scalarBar->SetLookupTable(glyph->GetLookupTable());
    scalarBar->SetOrientationToVertical();
    m_renderer->AddViewProp(scalarBar);

    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->SetTotalLength(1.5, 1.5, 1.5); 
    axes->SetShaftTypeToCylinder();
    axes->SetAxisLabels(true);

    m_widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    m_widget->SetOrientationMarker(axes);
    m_widget->SetInteractor(m_interactor);
    m_widget->SetViewport(0.0, 0.0, 0.2, 0.2); 
    m_widget->SetEnabled(1);
}


vtkSmartPointer<vtkTable> SeismicViewer::ExtractCoordinates(vtkSmartPointer<vtkPolyData> data, double center[3]) {
    auto table = vtkSmartPointer<vtkTable>::New();
    auto arrX = vtkSmartPointer<vtkDoubleArray>::New(); arrX->SetName("Northing");
    auto arrY = vtkSmartPointer<vtkDoubleArray>::New(); arrY->SetName("Easting");
    auto arrZ = vtkSmartPointer<vtkDoubleArray>::New(); arrZ->SetName("Depth");
    
    vtkIdType nPts = data->GetNumberOfPoints();
    arrX->SetNumberOfValues(nPts);
    arrY->SetNumberOfValues(nPts);
    arrZ->SetNumberOfValues(nPts);
    
    center[0] = center[1] = center[2] = 0.0;
    
    for (vtkIdType i = 0; i < nPts; ++i) {
        double p[3];
        data->GetPoint(i, p);
        arrX->SetValue(i, p[0]);
        arrY->SetValue(i, p[1]);
        arrZ->SetValue(i, p[2]);
        
        center[0] += p[0];
        center[1] += p[1];
        center[2] += p[2];
    }
    
    if (nPts > 0) {
        center[0] /= nPts;
        center[1] /= nPts;
        center[2] /= nPts;
    }
    
    table->AddColumn(arrX);
    table->AddColumn(arrY);
    table->AddColumn(arrZ);
    return table;
}


void SeismicViewer::ComputeAndVisualizePCA(vtkSmartPointer<vtkPolyData> data) {
    if (!data) return;

    double center[3]; // Output: Mean coordinates
    vtkSmartPointer<vtkTable> table = ExtractCoordinates(data, center);
    
    // PCA Calculation
    auto pca = vtkSmartPointer<vtkPCAStatistics>::New();
    pca->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, table);
    pca->SetColumnStatus("Northing", 1);
    pca->SetColumnStatus("Easting", 1);
    pca->SetColumnStatus("Depth", 1);
    pca->RequestSelectedColumns();
    pca->SetDeriveOption(true);
    pca->Update();

    // Get Results
    vtkNew<vtkDoubleArray> eigenvalues;
    pca->GetEigenvalues(eigenvalues);

    // Normalize before getting eigenvectors
    pca->SetNormalizationScheme(vtkPCAStatistics::DIAGONAL_VARIANCE);
    pca->Update();
    vtkNew<vtkDoubleArray> eigenvectors;
    pca->GetEigenvectors(eigenvectors);

    utils::printEigenvector(eigenvectors.Get());
    
    std::cout << "\n--- PCA Results ---" << std::endl;
    for (vtkIdType i = 0; i < eigenvalues->GetNumberOfTuples(); i++) {
        std::cout << "Eigenvalue " << i << " = " << std::fixed << std::setw(9)
                  << std::setprecision(6) << eigenvalues->GetValue(i) << std::endl;
    }

    double vizScaleFactor = 3.0; // 3-sigma equivalent for visualization
    VisualizeEigenvectors(center, eigenvalues, eigenvectors, vizScaleFactor);
    VisualizeConfidenceEllipsoid(center, eigenvalues, eigenvectors, vizScaleFactor);
    std::cout << "-------------------\n" << std::endl;
}

void SeismicViewer::VisualizeEigenvectors(const double center[3], vtkSmartPointer<vtkDoubleArray> eigenvalues,
        vtkSmartPointer<vtkDoubleArray> eigenvectors, double vizScaleFactor) {

    for (vtkIdType i = 0; i < eigenvalues->GetNumberOfTuples(); i++) {
        double eigenVal = eigenvalues->GetValue(i);
        double* eigenVec = eigenvectors->GetTuple(i);

        // Length = ScaleFactor * sqrt(Eigenvalue)
        double length = vizScaleFactor * std::sqrt(eigenVal);

        double endPoint[3];
        endPoint[0] = center[0] + (eigenVec[0] * length);
        endPoint[1] = center[1] + (eigenVec[1] * length);
        endPoint[2] = center[2] + (eigenVec[2] * length);

        vtkNew<vtkLineSource> lineSource;
        lineSource->SetPoint1(center);
        lineSource->SetPoint2(endPoint);

        vtkNew<vtkTubeFilter> tubeFilter;
        tubeFilter->SetInputConnection(lineSource->GetOutputPort());
        tubeFilter->SetRadius(10);
        tubeFilter->SetNumberOfSides(12);

        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(tubeFilter->GetOutputPort());

        vtkNew<vtkActor> actorEigen;
        actorEigen->SetMapper(mapper);

        // Color code: 0: Red, 1: Green, 2: Blue
        if (i == 0) actorEigen->GetProperty()->SetColor(1.0, 0.0, 0.0);
        else if (i == 1) actorEigen->GetProperty()->SetColor(0.0, 1.0, 0.0);
        else actorEigen->GetProperty()->SetColor(0.0, 0.0, 1.0);

        m_renderer->AddActor(actorEigen);
    }
}


void SeismicViewer::VisualizeConfidenceEllipsoid(const double center[3], vtkSmartPointer<vtkDoubleArray> eigenvalues, vtkSmartPointer<vtkDoubleArray> eigenvectors, double scaleFactor) {
    // 1. Setup the Transformation Matrix (Rotation, Scaling, Translation)
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->Identity();

    for (int i = 0; i < 3; ++i) {
        double vec[3];
        eigenvectors->GetTuple(i, vec);
        double val = eigenvalues->GetValue(i);
        
        // Scale factor: standard deviation * scaleFactor
        double geometryScale = scaleFactor * std::sqrt(val);

        // Set the column of the matrix
        matrix->SetElement(0, i, vec[0] * geometryScale);
        matrix->SetElement(1, i, vec[1] * geometryScale);
        matrix->SetElement(2, i, vec[2] * geometryScale);
    }

    // Set the translation
    matrix->SetElement(0, 3, center[0]);
    matrix->SetElement(1, 3, center[1]);
    matrix->SetElement(2, 3, center[2]);

    // 2. Create a sphere
    vtkNew<vtkSphereSource> source;
    source->SetRadius(1.0);
    source->SetThetaResolution(16);
    source->SetPhiResolution(16);

    // 3. Apply the transformation
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(matrix);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputConnection(source->GetOutputPort());
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    // 4. Mapper and Actor (Wireframe)
    vtkNew<vtkPolyDataMapper> ellipsoidMapper;
    ellipsoidMapper->SetInputConnection(transformFilter->GetOutputPort());

    vtkNew<vtkActor> ellipsoidActor;
    ellipsoidActor->SetMapper(ellipsoidMapper);

    ellipsoidActor->GetProperty()->SetRepresentationToWireframe();
    ellipsoidActor->GetProperty()->SetColor(0.0, 1.0, 0.0); // Green wireframe
    ellipsoidActor->GetProperty()->SetLineWidth(1.0); 
    m_renderer->AddActor(ellipsoidActor);
}


void SeismicViewer::Start() {
    m_renderer->ResetCamera();
    m_window->Render();
    m_interactor->Start();
}

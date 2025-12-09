#ifndef SEISMIC_VIEWER_H
#define SEISMIC_VIEWER_H

#include "utils.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkOrientationMarkerWidget.h>

#include <string>

class SeismicViewer {
public:
    /**
     * @brief Constructor that initializes the main VTK components.
     */
    SeismicViewer();

    /**
     * @brief Processes seismic event data, setting up the main visualization.
     * @param data The PolyData containing the seismic events (points and attributes).
     */
    void VisualizeEvents(vtkSmartPointer<vtkPolyData> data);

    /**
     * @brief Performs PCA on the data and visualizes the principal components
     * as vectors and a confidence ellipsoid.
     * @param data The PolyData input for PCA.
     */
    void ComputeAndVisualizePCA(vtkSmartPointer<vtkPolyData> data);

    /**
     * @brief Starts the VTK rendering loop.
     */
    void Start();

private:
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkRenderWindow> m_window;
    vtkSmartPointer<vtkRenderWindowInteractor> m_interactor;
    vtkSmartPointer<vtkOrientationMarkerWidget> m_widget;

    /**
     * @brief Creates the VTK glyph actor for seismic event points.
     * @param data The input PolyData.
     */
    void SetupEventGlyph(vtkSmartPointer<vtkPolyData> data);

    /**
     * @brief Adds the Bounding Box (Cube) actor to the scene.
     * @param bounds The extent [xmin, xmax, ymin, ymax, zmin, zmax] of the data.
     */
    void AddBoundingBox(const double bounds[6]);
    
    /**
     * @brief Adds the scalar bar legend and orientation widget.
     * @param glyph The vtkGlyph3DMapper used to determine the scalar range and lookup table.
     */
    void AddWidgets(vtkSmartPointer<class vtkGlyph3DMapper> glyph);

    /**
     * @brief Extracts coordinates from vtkPolyData into a vtkTable for PCA.
     * @param data The input PolyData.
     * @param center Output parameter for the calculated mean of the points.
     * @return The vtkTable with coordinate columns.
     */
    vtkSmartPointer<vtkTable> ExtractCoordinates(vtkSmartPointer<vtkPolyData> data, double center[3]);
    
    /**
     * @brief Visualizes the principal components as colored vectors.
     * @param center The mean of the data points.
     * @param eigenvalues PCA eigenvalues (variances).
     * @param eigenvectors PCA eigenvectors (directions).
     */
    void VisualizeEigenvectors(const double center[3], vtkSmartPointer<vtkDoubleArray> eigenvalues, vtkSmartPointer<vtkDoubleArray> eigenvectors, double vizScaleFactor);
    
    /**
     * @brief Visualizes the confidence ellipsoid based on PCA results.
     * @param center The mean of the data points.
     * @param eigenvalues PCA eigenvalues.
     * @param eigenvectors PCA eigenvectors.
     * @param scaleFactor Scale for the ellipsoid (e.g., 3.0 for 3-sigma).
     */
    void VisualizeConfidenceEllipsoid(const double center[3], vtkSmartPointer<vtkDoubleArray> eigenvalues, vtkSmartPointer<vtkDoubleArray> eigenvectors, double scaleFactor);
};

#endif
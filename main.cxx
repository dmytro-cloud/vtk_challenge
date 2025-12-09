#include "SeismicViewer.h"
#include "utils.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <cstdlib> // For EXIT_FAILURE, EXIT_SUCCESS

// --- Assumed External Functions (Implement these in a separate Utility file) ---

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: vtk_challenge <event_data.csv>\n";
        return EXIT_FAILURE;
    }

    // 1. Data Loading
    auto data = utils::LoadSeismicCSV(argv[1]);
    if (!data || data->GetNumberOfPoints() == 0) {
        std::cerr << "Error: Failed to load seismic data or data is empty.\n";
        return EXIT_FAILURE;
    }

    // 2. Application Setup and Visualization
    SeismicViewer viewer;

    // Setup basic scene (points and bounding box)
    viewer.VisualizeEvents(data);

    // Compute and draw PCA results (vectors and ellipsoid)
    viewer.ComputeAndVisualizePCA(data);

    // 3. Start the Render Loop
    viewer.Start();

    return EXIT_SUCCESS;
}

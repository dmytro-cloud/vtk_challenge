#include "utils.h"

namespace utils {
    vtkSmartPointer<vtkPolyData> LoadSeismicCSV(const std::string& filename)
    {
        auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
        reader->SetFileName(filename.c_str());
        // Had to change csv header because first line with names had an extra space after each name
        reader->SetFieldDelimiterCharacters(",");
        reader->SetHaveHeaders(true);
        reader->SetDetectNumericColumns(true);
        reader->Update();
    
        vtkTable* table = reader->GetOutput();
    
        // Print some information about the table
        std::cout << "Table has " << table->GetNumberOfRows() << " rows and "
        << table->GetNumberOfColumns() << " columns." << std::endl;
    
        // Optional: Print the header names
        for (vtkIdType i = 0; i < table->GetNumberOfColumns(); ++i) {
        std::cout << "Column " << i << ": " << table->GetColumnName(i) << std::endl;
        }
    
        // --- DEBUGGING STEP 1: Check the actual array type created by the reader ---
        std::cout << "\n--- Array Type Debugging ---" << std::endl;
        vtkAbstractArray* abstractNorthing = table->GetColumnByName("Northing");
        if (abstractNorthing) {
            std::cout << "Northing column type is: " << abstractNorthing->GetClassName() << std::endl;
        } else {
            std::cerr << "Error: 'Northing' column name not found at all." << std::endl;
            return nullptr;
        }
        
        // Extract columns by their header names
        vtkIntArray* northing = vtkIntArray::SafeDownCast(
            table->GetColumnByName("Northing"));
        vtkIntArray* easting = vtkIntArray::SafeDownCast(
            table->GetColumnByName("Easting"));
        vtkIntArray* depth = vtkIntArray::SafeDownCast(
            table->GetColumnByName("Depth"));
        vtkDoubleArray* magnitude = vtkDoubleArray::SafeDownCast(
            table->GetColumnByName("Moment Magnitude"));
    
        if (!northing || !easting || !depth || !magnitude)
        {
            std::cerr << "Error: Missing required columns in CSV.\n";
            return nullptr;
        }
    
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(table->GetNumberOfRows());
    
        for (vtkIdType i = 0; i < table->GetNumberOfRows(); ++i)
        {
            int x = northing->GetValue(i);
            int y = easting->GetValue(i);
            int z = depth->GetValue(i);
            // Converts to double internally (allegedly)
            points->SetPoint(i, x, y, z);
        }
    
        auto pd = vtkSmartPointer<vtkPolyData>::New();
        pd->SetPoints(points);
    
        // Attach magnitudes as scalars
        magnitude->SetName("Magnitude");
        pd->GetPointData()->AddArray(magnitude);
        pd->GetPointData()->SetScalars(magnitude);
    
        return pd;
    }
    
    void printEigenvector(vtkDoubleArray* eigenvectors) {
        for (vtkIdType i = 0; i < eigenvectors->GetNumberOfTuples(); i++)
        {
          std::cout << "Eigenvector " << i << " = ";
          double* evec = new double[eigenvectors->GetNumberOfComponents()];
          eigenvectors->GetTuple(i, evec);
          for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
          {
            if (j == 0)
              std::cout << "(";
            if (j < eigenvectors->GetNumberOfComponents() - 1)
            {
              std::cout << std::fixed << std::setw(9) << std::setprecision(6)
                        << evec[j] << ", ";
            }
            else
            {
              std::cout << std::fixed << std::setw(9) << std::setprecision(6)
                        << evec[j] << ")";
            }
          }
          delete[] evec;
          std::cout << std::endl;
        }
    }
    
}
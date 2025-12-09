#ifndef UTILS_H
#define UTILS_H

#include <vtkDelimitedTextReader.h>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <iostream>

namespace utils {

vtkSmartPointer<vtkPolyData> LoadSeismicCSV(const std::string& filename);
void printEigenvector(vtkDoubleArray* eigenvectors);

} // namespace utils

#endif
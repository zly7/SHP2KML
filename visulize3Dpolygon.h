
#ifndef VISUALIZE3DPOLYGON_H
#define VISUALIZE3DPOLYGON_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vector>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include "generate_path.h"
#include "random"
void visualizePolygonsAndPath(const std::vector<Polygon3D>& polygonPoints, const std::vector<std::vector<Point3D>>& result);
void visualizePath(const std::vector<Point3D>& result, vtkSmartPointer<vtkRenderer> renderer);
void visualizePolygons(const std::vector<Polygon3D>& polygonPoints, vtkSmartPointer<vtkRenderer> renderer, float polygonsToVisualizeProbability);

#endif // VISUALIZE3DPOLYGON_H
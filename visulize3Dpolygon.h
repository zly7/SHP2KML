
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

void visualizePolygonsAndPath(const std::vector<Polygon3D>& polygonPoints, const std::vector<Point3D>& result);

#endif // VISUALIZE3DPOLYGON_H
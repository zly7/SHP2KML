#include "visulize3Dpolygon.h"



void visualizePolygonsAndPath(const std::vector<Polygon3D>& polygonPoints, const std::vector<Point3D>& result) {
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    int totalPolygons = polygonPoints.size();
    int polygonsToVisualize = totalPolygons * 0.3;

    int counter = 0;
    // Visualize polygons
    for (const auto& polygon : polygonPoints) {
        if (counter >= polygonsToVisualize) {
            break;
        }
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (const auto& point : polygon) {
            points->InsertNextPoint(point.x, point.y, point.z);
        }

        vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(polygon.size() + 1);  // +1 to close the polygon
        for (unsigned int i = 0; i < polygon.size(); i++) {
            polyLine->GetPointIds()->SetId(i, i);
        }
        polyLine->GetPointIds()->SetId(polygon.size(), 0); 

        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        cells->InsertNextCell(polyLine);

        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetLines(cells);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        renderer->AddActor(actor);
        actor->GetProperty()->SetColor(0.0, 1.0, 0.0);  // Set polygon color to red

        counter++;
    }

    // Visualize path
    vtkSmartPointer<vtkPoints> pathPoints = vtkSmartPointer<vtkPoints>::New();
    for (const auto& point : result) {
        pathPoints->InsertNextPoint(point.x, point.y, point.z);
    }

    vtkSmartPointer<vtkPolyLine> pathPolyLine = vtkSmartPointer<vtkPolyLine>::New();
    pathPolyLine->GetPointIds()->SetNumberOfIds(result.size());
    for (unsigned int i = 0; i < result.size(); i++) {
        pathPolyLine->GetPointIds()->SetId(i, i);
    }

    vtkSmartPointer<vtkCellArray> pathCells = vtkSmartPointer<vtkCellArray>::New();
    pathCells->InsertNextCell(pathPolyLine);

    vtkSmartPointer<vtkPolyData> pathPolyData = vtkSmartPointer<vtkPolyData>::New();
    pathPolyData->SetPoints(pathPoints);
    pathPolyData->SetLines(pathCells);

    vtkSmartPointer<vtkPolyDataMapper> pathMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pathMapper->SetInputData(pathPolyData);

    vtkSmartPointer<vtkActor> pathActor = vtkSmartPointer<vtkActor>::New();
    pathActor->SetMapper(pathMapper);
    pathActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // Set path color to red

    renderer->AddActor(pathActor);


    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(0.05);  // Adjust the radius as needed

    // Use Glyph3D to place the sphere at each point in the path
    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
    glyph3D->SetInputData(pathPolyData);
    glyph3D->Update();

    // Create a mapper for the glyphs (blue dots)
    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph3D->GetOutputPort());

    // Create an actor for the glyphs
    vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(0.0, 0.0, 1.0);  // Set color to blue

    // Add the glyph actor to the renderer
    renderer->AddActor(glyphActor);

    renderer->SetBackground(1, 1, 1);  // Background color white

    // Create an axes actor
    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();

    // Create an orientation marker widget to display the axes
    vtkSmartPointer<vtkOrientationMarkerWidget> orientationMarker = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    orientationMarker->SetOutlineColor(0.9300, 0.5700, 0.1300);
    orientationMarker->SetOrientationMarker(axes);
    orientationMarker->SetInteractor(renderWindowInteractor);
    orientationMarker->SetViewport(0.0, 0.0, 0.2, 0.2);
    orientationMarker->SetEnabled(1);
    orientationMarker->InteractiveOn();

    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();
}


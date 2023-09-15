#include "generate_path.h"
#include "Point2KML_refer.h"


double distance(const Point3D& a, const Point3D& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

std::vector<Point3D> greedyTSP(const Point3D& start, const std::vector<Point3D>& points) {
    std::vector<Point3D> result;
    std::vector<bool> visited(points.size(), false); // To track visited points

    Point3D current = start;
    result.push_back(current); // Start point

    for (size_t i = 0; i < points.size(); ++i) {
        double minDist = std::numeric_limits<double>::max();
        int nextPointIndex = -1;

        for (size_t j = 0; j < points.size(); ++j) {
            if (!visited[j]) {
                double dist = distance(current, points[j]);
                if (dist < minDist) {
                    minDist = dist;
                    nextPointIndex = j;
                }
            }
        }

        if (nextPointIndex != -1) {
            visited[nextPointIndex] = true;
            current = points[nextPointIndex];
            result.push_back(current);
        }
    }

    return result;
}

std::vector<Point3D> interpolatePoints(const std::vector<Point3D>& originalResult, double assumedZ, double interval) {
    std::vector<Point3D> interpolatedResult;

    for (size_t i = 0; i < originalResult.size() - 1; ++i) {
        const Point3D& A = originalResult[i];
        const Point3D& B = originalResult[i + 1];

        double dist = distance(A, B);
        int numInterpolatedPoints = static_cast<int>(dist / interval);

        interpolatedResult.push_back(A);  // 添加起始点

        for (int j = 1; j < numInterpolatedPoints; ++j) {
            double t = static_cast<double>(j) * interval / dist;
            Point3D interpolatedPoint = {
                A.x + t * (B.x - A.x),
                A.y + t * (B.y - A.y),
                assumedZ  // 使用给定的Z坐标
            };
            interpolatedResult.push_back(interpolatedPoint);
        }
    }

    // 添加最后一个点
    if (!originalResult.empty()) {
        interpolatedResult.push_back(originalResult.back());
    }

    return interpolatedResult;
}



std::vector<Waypoint> convertToWaypoints(const std::vector<Point3D>& points3D,double centralMeridian) {
    std::vector<Waypoint> waypoints;
    for (const auto& point : points3D) {
        Waypoint waypoint;
        waypoint.lon = point.x + centralMeridian; // 经度
        waypoint.lat = point.y; // 纬度
        waypoint.alt = point.z; // 海拔
        
        // 设置默认值
        waypoint.x = 0.0;       // 如果Waypoint的x有其他用途，可以根据实际情况进行修改
        waypoint.y = 0.0;       // 如果Waypoint的y有其他用途，可以根据实际情况进行修改
        waypoint.pitch = 0.0;   // 可以根据实际情况进行修改
        waypoint.heading = 0.0; // 可以根据实际情况进行修改
        waypoint.baseLine = 0.0; // 可以根据实际情况进行修改
        waypoint.nStripID = 0;  // 可以根据实际情况进行修改

        waypoints.push_back(waypoint);
    }
    return waypoints;
}


int main() {

    GDALAllRegister();  // 注册所有的驱动

    GDALDataset *poDS;
    poDS = (GDALDataset*) GDALOpenEx("../3Dshp/BLOCK.shp", GDAL_OF_VECTOR, NULL, NULL, NULL);
    if(poDS == NULL) {
        std::cout << "Open failed. Reason: " << CPLGetLastErrorMsg() << std::endl;
        exit(1);
    }

    // 获取数据集的空间参考系统
    const char* pszWKT = poDS->GetProjectionRef();
    OGRSpatialReference oSRS;
    oSRS.importFromWkt(&pszWKT);

    // 获取中心经度
    double centralMeridian = oSRS.GetProjParm("central_meridian", -9999.0);  // 默认值为-9999.0，以便我们知道是否正确地读取了中心经度

    if (centralMeridian != -9999.0) {
        std::cout << "Central Meridian: " << centralMeridian << "°E" << std::endl;
    } else {
        std::cout << "Failed to retrieve central meridian." << std::endl;
    }

    double minX = DBL_MAX, minY = DBL_MAX, maxX = -DBL_MAX, maxY = -DBL_MAX;
    OGRLayer *poLayer;
    poLayer = poDS->GetLayer(0);

    // Create an OpenCV image for visualization. Adjust the size as needed.
    cv::Mat img = cv::Mat::zeros(5000, 5000, CV_8UC3);

    OGRFeature *poFeature;
    poLayer->ResetReading();

    std::vector<Point3D> centerPoints; //si
    std::vector<Vector3D> normalVectors; //ni
    double d_SGD = 20; //d_SGD -meters
    std::vector<Point3D> viewPoints; //vi
    std::vector<Polygon3D> polygonPoints;

    while((poFeature = poLayer->GetNextFeature()) != NULL) {
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if(poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon) {
            OGRPolygon *poPolygon = (OGRPolygon *) poGeometry;
            OGRLinearRing *poRing = poPolygon->getExteriorRing();
            Polygon3D poly;
                for(int i = 0; i < poRing->getNumPoints(); i++) {
                OGRPoint oPoint;
                poRing->getPoint(i, &oPoint);
                std::cout << "X: " << oPoint.getX() << ", Y: " << oPoint.getY() << std::endl;
                if(oPoint.getX() < minX) minX = oPoint.getX();
                if(oPoint.getY() < minY) minY = oPoint.getY();
                if(oPoint.getX() > maxX) maxX = oPoint.getX();
                if(oPoint.getY() > maxY) maxY = oPoint.getY();
                poly.points[i] = {oPoint.getX(), oPoint.getY(), oPoint.getZ()};
            }
            polygonPoints.push_back(poly);
            normalVectors.push_back(averageNormalVector(poly));
            OGRPoint centroid;
            if (poPolygon->Centroid(&centroid) != OGRERR_NONE) {
                std::cout << "Failed to compute centroid." << std::endl;
                exit(1);
            }
            centerPoints.push_back({centroid.getX(), centroid.getY(), centroid.getZ()});
            
        } else {
            std::cout << "No geometry." << std::endl;
        }

        OGRFeature::DestroyFeature(poFeature);
    }

    // generate viewponts by project
    for (size_t i = 0; i < centerPoints.size(); ++i) {
        Point3D centerPoint = centerPoints[i];
        Vector3D normalVector = normalVectors[i];
        Point3D viewPoint = {
            centerPoint.x + d_SGD * normalVector.x,
            centerPoint.y + d_SGD * normalVector.y,
            centerPoint.z + d_SGD * normalVector.z
        };
        viewPoints.push_back(viewPoint);
    }

    double d_prj = 0.958*4.5;
    double r_overlap = 0.3;
    double D_disk = d_prj * (1-r_overlap);

    //greedyTSP
    std::vector<Point3D> result = greedyTSP(centerPoints[0], centerPoints);
    double interval = 1.0;  // 1米间隔
    double assumedZValue = 30.0;  // 假设的Z坐标
    std::vector<Point3D> newResult = interpolatePoints(result, assumedZValue, interval);
    std::vector<Waypoint> waypoints = convertToWaypoints(newResult, centralMeridian);
    generateKMLFromWaypoints(waypoints, "final_output.kml");

    double scaleX = img.cols / (maxX - minX);
    double scaleY = img.rows / (maxY - minY);
    poLayer->ResetReading();
    while((poFeature = poLayer->GetNextFeature()) != NULL) {
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();

        if(poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon) {
            OGRPolygon *poPolygon = (OGRPolygon *) poGeometry;
            OGRLinearRing *poRing = poPolygon->getExteriorRing();

            std::vector<cv::Point> cvPoints;
            for(int i = 0; i < poRing->getNumPoints(); i++) {
                OGRPoint oPoint;
                poRing->getPoint(i, &oPoint);
                // Convert OGRPoint to cv::Point. Adjust the scaling and offset as needed.
                cv::Point cvPoint((oPoint.getX() - minX) * scaleX, (maxY - oPoint.getY()) * scaleY);
                cvPoints.push_back(cvPoint);
            }
            cv::polylines(img, cvPoints, true, cv::Scalar(0, 255, 0), 1);
        } else {
            std::cout << "No geometry." << std::endl;
        }

        OGRFeature::DestroyFeature(poFeature);
    }

    std::vector<cv::Point> pathPoints;
    for (const auto& point : result) {
        cv::Point cvPoint((point.x - minX) * scaleX, (maxY - point.y) * scaleY);
        pathPoints.push_back(cvPoint);
    }

    // 在图像上绘制路径
    cv::polylines(img, pathPoints, false, cv::Scalar(255, 0, 0), 2);  // 使用红色线条绘制路径

    // 如果需要在路径上的每个点处绘制小圆点，可以使用以下代码：
    for (const auto& point : pathPoints) {
        cv::circle(img, point, 3, cv::Scalar(0, 0, 255), -1);  // 使用红色圆点标记路径上的每个点
    }


    GDALClose(poDS);

    cv::imshow("Visualization", img);
    cv::waitKey(0);
    cv::imwrite("output.png", img);


    return 0;
}

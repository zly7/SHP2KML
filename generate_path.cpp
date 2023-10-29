#include "generate_path.h"
#include "Point2KML_refer.h"
#include "PoissonDiskSampler.h"
#include "visulize3Dpolygon.h"
#include <Cluster_and_direct_planning.h>


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
                double dist = Point3D::distance(current, points[j]);
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

std::vector<Point3D> interpolatePoints(const std::vector<Point3D>& originalResult, double interval) {
    std::vector<Point3D> interpolatedResult;

    for (size_t i = 0; i < originalResult.size() - 1; ++i) {
        const Point3D& A = originalResult[i];
        const Point3D& B = originalResult[i + 1];

        double dist = Point3D::distance(A, B);
        int numInterpolatedPoints = static_cast<int>(dist / interval);

        interpolatedResult.push_back(A);  // 添加起始点

        for (int j = 1; j < numInterpolatedPoints; ++j) {
            double t = static_cast<double>(j) * interval / dist;
            Point3D interpolatedPoint = {
                A.x + t * (B.x - A.x),
                A.y + t * (B.y - A.y),
                A.z + t * (B.z - A.z)
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

std::vector<std::vector<Point3D>> interpolatePoints(const std::vector<std::vector<Point3D>>& originalResult2D, double interval) {
    std::vector<std::vector<Point3D>> interpolatedResult2D;
    for(const auto& originalResult : originalResult2D) {
        std::vector<Point3D> interpolatedResult = interpolatePoints(originalResult, interval);
        interpolatedResult2D.push_back(interpolatedResult);
    }
    return interpolatedResult2D;
}




std::vector<Waypoint> convertToWaypoints(const std::vector<Point3D>& points3D,double centralMeridian) {
        std::vector<Waypoint> waypoints(points3D.size());

    // 创建空间参考
    OGRSpatialReference wGS84, cGCS114E;
    wGS84.importFromEPSG(4326);
    cGCS114E.importFromEPSG(4547);

    // 创建坐标转换对象
    OGRCoordinateTransformation* proj_trans = OGRCreateCoordinateTransformation(&cGCS114E, &wGS84);

    // 准备坐标数组
    std::vector<double> x(points3D.size()), y(points3D.size()), z(points3D.size());
    for (size_t i = 0; i < points3D.size(); ++i) {
        x[i] = points3D[i].x;
        y[i] = points3D[i].y;
        z[i] = points3D[i].z;
    }

    // 批量转换坐标
    proj_trans->Transform(points3D.size(), x.data(), y.data(), z.data());
    for (size_t i = 0; i < points3D.size(); ++i) {
        waypoints[i].lon = x[i] + centralMeridian;
        waypoints[i].lat = y[i] ;
        waypoints[i].alt = z[i];
        waypoints[i].x = 0.0;
        waypoints[i].y = 0.0;
        waypoints[i].pitch = -90;  // 相机垂直对准地面
        waypoints[i].heading = 0.0;  // 相机垂直对准地面
        waypoints[i].baseLine = 0.0;
        waypoints[i].nStripID = 0;
    }
    return waypoints;
}

std::vector<std::vector<Waypoint>> convertToWaypoints(const std::vector<std::vector<Point3D>>& points3D2D, double centralMeridian) {
    std::vector<std::vector<Waypoint>> waypoints2D;
    for(const auto& points3D : points3D2D) {
        std::vector<Waypoint> waypoints = convertToWaypoints(points3D, centralMeridian);
        waypoints2D.push_back(waypoints);
    }
    return waypoints2D;
}



bool generateKMLsFromWaypoints(const std::vector<std::vector<Waypoint>>& waypointsGroups, const std::string& baseFilename) {
    bool allSuccessful = true; // 用于跟踪所有文件是否都成功生成

    for(size_t i = 0; i < waypointsGroups.size(); ++i) {
        const auto& waypoints = waypointsGroups[i];

        // 生成文件名
        std::ostringstream filenameStream;
        filenameStream << baseFilename << (i + 1) << ".kml"; // i + 1 以确保从 1 开始编号
        std::string filename = filenameStream.str();

        // 调用 generateKMLFromWaypoints 函数
        bool success = generateKMLFromWaypoints(waypoints, filename);
        if(!success) {
            allSuccessful = false; // 如果任何一个文件生成失败，更新 allSuccessful 标志
        }
    }

    return allSuccessful; // 如果所有文件都成功生成，返回 true，否则返回 false
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
    double d_SGD = 10; //d_SGD -meters
    std::vector<Point3D> viewPoints; //vi
    std::vector<Polygon3D> polygonPoints;

    while((poFeature = poLayer->GetNextFeature()) != NULL) {
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if(poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon) {
            OGRPolygon *poPolygon = (OGRPolygon *) poGeometry;
            OGRLinearRing *poRing = poPolygon->getExteriorRing();
            Polygon3D poly;
            double totalZ = 0.0; // 用于累加所有顶点的Z值
            for(int i = 0; i < poRing->getNumPoints(); i++) {
                OGRPoint oPoint;
                poRing->getPoint(i, &oPoint);
                std::cout << "X: " << oPoint.getX() << ", Y: " << oPoint.getY() << ", Z: " << oPoint.getZ() << std::endl;
                if(oPoint.getX() < minX) minX = oPoint.getX();
                if(oPoint.getY() < minY) minY = oPoint.getY();
                if(oPoint.getX() > maxX) maxX = oPoint.getX();
                if(oPoint.getY() > maxY) maxY = oPoint.getY();
                poly.points[i] = {oPoint.getX(), oPoint.getY(), oPoint.getZ()};
                totalZ += oPoint.getZ();
            }
            double averageZ = totalZ / poRing->getNumPoints(); // 计算Z值的平均值
            polygonPoints.push_back(poly);
            normalVectors.push_back(averageNormalVector(poly));
            OGRPoint centroid;
            if (poPolygon->Centroid(&centroid) != OGRERR_NONE) {
                std::cout << "Failed to compute centroid." << std::endl;
                exit(1);
            }
            centerPoints.push_back({centroid.getX(), centroid.getY(), averageZ});
            
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
    double totalHeight = 0.0;
    for(const auto& viewPoint : viewPoints) {
        totalHeight += viewPoint.z; // 假设 z 表示高度
    }
    double averageHeightForViewpoint = totalHeight / viewPoints.size();

    double d_prj = 0.958*4.5;
    double r_overlap = 0.3;
    double D_disk = d_prj * (1-r_overlap);
    float cell_size = 2; // for hash in poisson disk sampling

    std::vector<std::vector<Point3D>> resultFrom;
    std::cout << "algorithm name: " << GlobalConstants::ALGORITHM_NAME << std::endl;
    if(strcmp(GlobalConstants::ALGORITHM_NAME,"possion_sample_and_greedyTSP")==0){

        PoissonDiskSampler sampler(cell_size, viewPoints);
        std::vector<Point3D> sample_viewpoints = sampler.PoissonDiskSampling(D_disk, viewPoints);
        
        std::vector<Point3D> resultFromTemp = greedyTSP(sample_viewpoints[0], sample_viewpoints);
        resultFrom.push_back(resultFromTemp);
    }else if (strcmp(GlobalConstants::ALGORITHM_NAME,"cluster_and_direct_planning")==0){
        int k = GlobalConstants::num_planes; 
        std::vector<std::vector<Point3D>>  clusters = kMeansClustering(viewPoints, k);
        resultFrom = generateZigzagTrajectories(clusters, D_disk, averageHeightForViewpoint);
    }else{
        std::cout<<"wrong algorithm name"<<std::endl;
    }
    
    double interval = 1.0;  // 1米间隔
    std::vector<std::vector<Point3D>> newResult = interpolatePoints(resultFrom, interval);
    std::vector<std::vector<Waypoint>> waypoints = convertToWaypoints(newResult, centralMeridian);
    generateKMLsFromWaypoints(waypoints, "final_output");

    double scaleX = img.cols / (maxX - minX);
    double scaleY = img.rows / (maxY - minY);


    visualizePolygonsAndPath(polygonPoints, resultFrom);
    
    


    GDALClose(poDS);

    cv::imshow("Visualization", img);
    cv::waitKey(0);
    cv::imwrite("output.png", img);


    return 0;
}

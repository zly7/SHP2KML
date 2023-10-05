#include "Cluster_and_direct_planning.h"


std::vector<std::vector<Point3D>> generateZigzagTrajectories(const std::vector<std::vector<Point3D>>& clusters, double d, double h) {
    std::vector<std::vector<Point3D>> trajectories;
    
    for(const auto& cluster : clusters) {
        if(cluster.empty()) continue;
        
        // Project 3D points to 2D
        std::vector<cv::Point2f> points2D;
        for(const auto& point : cluster) {
            points2D.emplace_back(point.x, point.y);
        }
        
        // Compute the minimum bounding rectangle
        cv::RotatedRect minRect = cv::minAreaRect(points2D);
        
        cv::Point2f vertices[4];
        minRect.points(vertices);

        // 计算每个顶点到原点的距离并找到最近的顶点
        cv::Point2f start = vertices[0];
        float minDistSquared = start.x * start.x + start.y * start.y; // 初始最小距离的平方

        for(int i = 1; i < 4; ++i) {
            float distSquared = vertices[i].x * vertices[i].x + vertices[i].y * vertices[i].y;
            if(distSquared < minDistSquared) {
                start = vertices[i];
                minDistSquared = distSquared;
            }
        }
        // Determine the step size along the shorter and longer edge of the rectangle
        double stepShort = std::min(minRect.size.width, minRect.size.height);
        double stepLong = std::max(minRect.size.width, minRect.size.height);
        
        // Determine the number of steps along the short and long edges
        int numStepsShort = std::ceil(stepShort / (2*d));
        int numStepsLong = std::ceil(stepLong / (2*d));
        cv::Point2f directionLong, directionShort;

        // Find the directions of the short and long edges
        float minDist1 = std::numeric_limits<float>::max();
        float minDist2 = std::numeric_limits<float>::max();
        cv::Point2f adjacentVertex1, adjacentVertex2;

        for(const auto& vertex : vertices) {
            if(vertex == start) continue;
            float distSquared = (vertex.x - start.x) * (vertex.x - start.x) + (vertex.y - start.y) * (vertex.y - start.y);
            if(distSquared < minDist1) {
                minDist2 = minDist1;
                adjacentVertex2 = adjacentVertex1;
                
                minDist1 = distSquared;
                adjacentVertex1 = vertex;
            } else if(distSquared < minDist2) {
                minDist2 = distSquared;
                adjacentVertex2 = vertex;
            }
        }

        //当在XY坐标轴的时候，一定要善于使用单位方向向量
        // Calculate the vectors from start point to the adjacent vertices
        cv::Point2f vector1 = adjacentVertex1 - start;
        cv::Point2f vector2 = adjacentVertex2 - start;

        // Calculate the squared lengths of the vectors
        float lengthSquared1 = vector1.x * vector1.x + vector1.y * vector1.y;
        float lengthSquared2 = vector2.x * vector2.x + vector2.y * vector2.y;

        // Normalize the vectors to get the unit vectors (directions)
        cv::Point2f unitVector1 = vector1 * (1.0 / std::sqrt(lengthSquared1));
        cv::Point2f unitVector2 = vector2 * (1.0 / std::sqrt(lengthSquared2));

        // Determine which vector is the direction of the long edge and which is the direction of the short edge
        if(lengthSquared1 > lengthSquared2) {
            directionLong = unitVector1;
            directionShort = unitVector2;
        } else {
            directionLong = unitVector2;
            directionShort = unitVector1;
        }

        std::vector<Point3D> zPath3D;
        zPath3D.reserve(numStepsShort * numStepsLong);

        cv::Point2f currentPoint = start;
        bool reverse = false; // Flag to determine the direction along the short edge
        float zCoordinate = h; // Assuming z-coordinate is known or fixed

        for(int i = 0; i < numStepsLong; ++i) {
            for(int j = 0; j < numStepsShort; ++j) {
                zPath3D.emplace_back(currentPoint.x, currentPoint.y, zCoordinate);
                currentPoint += (reverse ? -1 : 1) * directionShort * (2 * d);
            }
            // After completing one row, move to the next row along the long edge
            currentPoint += directionLong * (2 * d);
            // Reverse the direction along the short edge for the next row
            reverse = !reverse;
        }
        
        // Add the computed trajectory to the result
        trajectories.push_back(zPath3D);
    }
    
    return trajectories;
}

std::vector<std::vector<Point3D>> kMeansClustering(const std::vector<Point3D>& points, int k) {
    if (points.empty() || k <= 0 || k > points.size()) return {};

    // 初始化聚类中心
    std::vector<Point3D> centroids;
    for (int i = 0; i < k; ++i) {
        centroids.push_back(points[rand() % points.size()]);
    }

    std::vector<std::vector<Point3D>> clusters(k);
    bool changed = true;
    while (changed) {
        changed = false;
        // 清空当前的聚类
        for (auto& cluster : clusters) cluster.clear();

        // 将每个点分配到最近的聚类中心
        for (const auto& point : points) {
            double minDist = std::numeric_limits<double>::max();
            int minIndex = -1;
            for (int i = 0; i < k; ++i) {
                double dist = Point3D::distance(point, centroids[i]);
                if (dist < minDist) {
                    minDist = dist;
                    minIndex = i;
                }
            }
            clusters[minIndex].push_back(point);
        }

        // 重新计算聚类中心
        for (int i = 0; i < k; ++i) {
            double x = 0, y = 0, z = 0;
            for (const auto& point : clusters[i]) {
                x += point.x;
                y += point.y;
                z += point.z;
            }
            Point3D newCentroid{x / clusters[i].size(), y / clusters[i].size(), z / clusters[i].size()};
            if (Point3D::distance(newCentroid, centroids[i]) > 1e-5) {
                centroids[i] = newCentroid;
                changed = true;
            }
        }
    }
    return clusters;
}

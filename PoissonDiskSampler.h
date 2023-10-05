#ifndef POISSON_DISK_SAMPLER_H
#define POISSON_DISK_SAMPLER_H

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include "generate_path.h"

class PoissonDiskSampler {
private:
    float cell_size, minX, minY, minZ;
    int x_multiplier, y_multiplier;
    std::vector<Point3D> samplePool;
    std::unordered_map<int64_t, std::vector<Point3D>> cells;



    Point3D extractFromSamplePool(std::vector<Point3D>& samplePool);

    std::unordered_map<int64_t, std::vector<Point3D>> fillSpatialHashTable(std::vector<Point3D>& samplePool, float cellSize);

    void removeSamples(const Point3D& p, float radius, std::unordered_map<int64_t, std::vector<Point3D>>& cells, 
        std::vector<Point3D>& samplePool, int x_multiplier, int y_multiplier,float);

    static int64_t computeHash(int cellX, int cellY, int cellZ, int x_multiplier, int y_multiplier);

public:
    PoissonDiskSampler(float cell_size, std::vector<Point3D>& viewPoints);



    std::vector<Point3D> PoissonDiskSampling(float radius, std::vector<Point3D>& viewPoints);
};

#endif // POISSON_DISK_SAMPLER_H

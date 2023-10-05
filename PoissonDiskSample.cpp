#include "PoissonDiskSampler.h"



PoissonDiskSampler::PoissonDiskSampler(float cell_size, std::vector<Point3D>& viewPoints) : cell_size(cell_size), samplePool(viewPoints) {
    // 初始化 x_multiplier 和 y_multiplier
    // 初始化 cells
    // 计算整个点集的范围
    minX = std::min_element(samplePool.begin(), samplePool.end(), [](const Point3D& a, const Point3D& b) { return a.x < b.x; })->x;
    minY = std::min_element(samplePool.begin(), samplePool.end(), [](const Point3D& a, const Point3D& b) { return a.y < b.y; })->y;
    minZ = std::min_element(samplePool.begin(), samplePool.end(), [](const Point3D& a, const Point3D& b) { return a.z < b.z; })->z;

    // 计算哈希函数的乘数
    x_multiplier = std::ceil((std::max_element(samplePool.begin(), samplePool.end(), [](const Point3D& a, const Point3D& b) { return a.x < b.x; })->x - minX) /cell_size);
    y_multiplier = std::ceil((std::max_element(samplePool.begin(), samplePool.end(), [](const Point3D& a, const Point3D& b) { return a.y < b.y; })->y - minY) / cell_size);

    assert(minX > 0 && "minX should be greater than 0!");
    assert(minY > 0 && "minY should be greater than 0!");
    assert(minZ > 0 && "minZ should be greater than 0!");
    assert(x_multiplier > 0 && "x_multiplier should be greater than 0!");
    assert(y_multiplier > 0 && "y_multiplier should be greater than 0!");

}

Point3D PoissonDiskSampler::extractFromSamplePool(std::vector<Point3D>& samplePool) {
    // 生成一个随机索引
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, samplePool.size() - 1);

    int randomIndex = distrib(gen);

    // 从samplePool中取出并移除该样本点
    Point3D selectedPoint = samplePool[randomIndex];
    samplePool.erase(samplePool.begin() + randomIndex);

    return selectedPoint;
}

int64_t PoissonDiskSampler::computeHash(int cellX, int cellY, int cellZ, int x_multiplier, int y_multiplier) {
    return static_cast<int64_t>(cellZ) + 
           static_cast<int64_t>(y_multiplier) * (static_cast<int64_t>(cellY) + 
           static_cast<int64_t>(x_multiplier) * static_cast<int64_t>(cellX));
}


std::unordered_map<int64_t, std::vector<Point3D>> PoissonDiskSampler::fillSpatialHashTable(std::vector<Point3D>& samplePool, float cellSize) {
    std::unordered_map<int64_t, std::vector<Point3D>> cells;

    for (auto& point : samplePool) {
        Point3D normalizedPoint = point;  // protect the origin value,zly
        normalizedPoint.x -= minX;
        normalizedPoint.y -= minY;
        normalizedPoint.z -= minZ;

        // 计算点的哈希值
        int cellX = std::floor(normalizedPoint.x / cellSize);
        int cellY = std::floor(normalizedPoint.y / cellSize);
        int cellZ = std::floor(normalizedPoint.z / cellSize);
        int64_t hash = computeHash(cellX, cellY, cellZ, x_multiplier, y_multiplier);

        // 将点添加到相应的单元格中
        cells[hash].push_back(point);
    }
    return cells;
}

void PoissonDiskSampler::removeSamples(const Point3D& p, float radius, std::unordered_map<int64_t, std::vector<Point3D>>& cells, 
            std::vector<Point3D>& samplePool, int x_multiplier, int y_multiplier,float cellSize) {
    // 计算需要检查的单元格范围
    int minCellX = std::floor((p.x-minX - radius) / cellSize);
    int maxCellX = std::ceil((p.x-minX + radius) / cellSize);
    int minCellY = std::floor((p.y - minY - radius) / cellSize);
    int maxCellY = std::ceil((p.y - minY + radius) / cellSize);
    int minCellZ = std::floor((p.z- minZ - radius) / cellSize);
    int maxCellZ = std::ceil((p.z- minZ + radius) / cellSize);

    // 遍历所有相关单元格
    for (int x = minCellX; x <= maxCellX; ++x) {
        for (int y = minCellY; y <= maxCellY; ++y) {
            for (int z = minCellZ; z <= maxCellZ; ++z) {
                int64_t hash = computeHash(x, y, z, x_multiplier, y_multiplier);
                if (cells.find(hash) != cells.end()) {
                    auto& cell = cells[hash];
                    for (auto it = cell.begin(); it != cell.end(); /* empty */) {
                        Point3D& point = *it;
                        float distance = std::sqrt(std::pow(point.x - p.x, 2) + std::pow(point.y - p.y, 2) + std::pow(point.z - p.z, 2));
                        if (distance < radius) {
                            // 从samplePool中移除点
                            samplePool.erase(std::remove(samplePool.begin(), samplePool.end(), point), samplePool.end());
                            // 从当前单元格中移除点
                            it = cell.erase(it);
                            std::cout << "remove one point" << std::endl;
                        } else {
                            ++it;
                        }
                    }
                }
            }
        }
    }
}


std::vector<Point3D> PoissonDiskSampler::PoissonDiskSampling(float radius, std::vector<Point3D>& viewPoints) {


    // 2. 填充空间索引以快速访问样本,whcih means that you  can easily do delete operation
    cells = fillSpatialHashTable(viewPoints,cell_size);

    

    // 3. 主循环
    std::vector<Point3D> samples;
    while (!samplePool.empty()) {
        // 从样本池中提取一个有效样本
        Point3D p = extractFromSamplePool(samplePool);
        samples.push_back(p);

        // 在磁盘中移除样本
        removeSamples(p, radius, cells, samplePool,x_multiplier,y_multiplier,cell_size);
    }
    std::cout << "The final sample vector that removes points within a radius of randomly selected points in PossonDisk: " << samples.size() << std::endl;
    return samples;
}

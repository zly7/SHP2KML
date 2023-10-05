# cmake_minimum_required(VERSION 3.10) # 指定CMake的最低版本
# project(point2kml) # 指定项目的名称

# # 添加所有的.cpp文件到一个变量中
# set(SOURCE_FILES
#     Point2KML_linux.cpp
#     tinyxml2.cpp
# )

# # 创建一个执行文件
# add_executable(MyExecutableName ${SOURCE_FILES})

cmake_minimum_required(VERSION 3.10)

project(ShapefileReader)

# Find GDAL
find_package(GDAL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(ShapefileReader generate_path.cpp Point2KML_refer.cpp tinyxml2.cpp PoissonDiskSample.cpp visulize3Dpolygon.cpp)


# Link
# Link
target_link_libraries(ShapefileReader PRIVATE 
    ${GDAL_LIBRARY}
    ${OpenCV_LIBS}
    ${VTK_LIBRARIES}
)
target_include_directories(ShapefileReader PRIVATE ${GDAL_INCLUDE_DIR})


cmake_minimum_required(VERSION 3.20.0)

project(MFS_Sloshing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

add_executable(${CMAKE_PROJECT_NAME}
    ${PROJECT_SOURCE_DIR}/src/LRBF/main.cpp
    ${PROJECT_SOURCE_DIR}/src/LRBF/Collocation2D.h
    ${PROJECT_SOURCE_DIR}/src/LRBF/GMMMFSSloshing2D.h
    ${PROJECT_SOURCE_DIR}/src/LRBF/GMMMQBasis2D.h
    ${PROJECT_SOURCE_DIR}/src/LRBF/GMMRectangular.h
    ${PROJECT_SOURCE_DIR}/src/LRBF/KDTreeTsaiAdaptor.h
)

target_include_directories(${CMAKE_PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src/gmm
        ${PROJECT_SOURCE_DIR}/src/nanoflann
)

target_compile_options(${CMAKE_PROJECT_NAME}
    PUBLIC
        -D_SCL_SECURE_NO_DEPRECATE
)
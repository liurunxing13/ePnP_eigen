cmake_minimum_required(VERSION 2.8)

project(epnp_eigen)

#set(OpenCV_DIR /usr/local/share/OpenCV)
#find_package(OpenCV REQUIRED)

INCLUDE_DIRECTORIES(/usr/include/eigen3)
INCLUDE_DIRECTORIES(include)
set(src_cpp
    src/epnp.cpp
    src/main.cpp)
add_executable(${PROJECT_NAME} ${src_cpp})
#target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})

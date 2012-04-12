cmake_minimum_required(VERSION 2.8)

project( MPI_CORRELATION )


    #set(HDF5_LIBS hdf5 hdf5_cpp hdf5_hl)
    #set(HDF5_LIBS hdf5 hdf5_hl hdf5_cpp)
    find_package(HDF5 REQUIRED)
    ##include_directories(/opt/local/hadoop-home/temp/CGAL-3.8/include)
    #include_directories(/home/sagrava/work/virtual-gis/trunk/spatial/downloads/CGAL-3.8-beta1)
    #include_directories(/home/sagrava/work/virtual-gis/trunk/spatial/local/include/spatialindex)
    #include_directories(/home/sagrava/work/virtual-gis/trunk/spatial/local/include/boost)
    SET(CMAKE_CXX_COMPILER mpicxx)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
   
    set(SRCS rho.cpp h5_read.cpp h5_write.cpp )
   
    add_executable( rho ${SRCS} )
    target_link_libraries(rho ${HDF5_LIBRARIES})
    

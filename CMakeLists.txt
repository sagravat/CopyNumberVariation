cmake_minimum_required(VERSION 2.8)

project( MPI_CORRELATION )


    #set(HDF5_LIBS hdf5 hdf5_cpp hdf5_hl)
    #set(HDF5_LIBS hdf5 hdf5_hl hdf5_cpp)
    find_package(HDF5 REQUIRED)


if(TYPE)
    message("TYPE EQUAL 1")
    SET(CMAKE_CXX_COMPILER g++)
    SET(CMAKE_BUILD_TYPE Release)
   
    set(DEBUG_SRCS h5dump.cpp h5_read.cpp h5_write.cpp )
    set(HYPER_SRCS hyperread.cpp )
   
    add_executable( h5dump ${DEBUG_SRCS} )
    add_executable( hyper ${HYPER_SRCS} )

    target_link_libraries(h5dump ${HDF5_LIBRARIES})
    target_link_libraries(hyper ${HDF5_LIBRARIES})


else()
    message("TYPE NOT EQUAL 1")
    SET(CMAKE_CXX_COMPILER mpicxx)
    SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
   
    set(OMP_SRCS rho_omp.cpp h5_read.cpp h5_write.cpp h5_hyperslab.cpp )
    set(PEARSON_SRCS pearson_small_block.cpp h5_read.cpp h5_write.cpp h5_hyperslab.cpp )
    set(PEARSON_RING_SRCS pearson_block_ring.cpp h5_read.cpp h5_write.cpp h5_hyperslab.cpp )
    set(PEARSON_HYPER_SRCS pearson_hyperslab.cpp h5_read.cpp h5_write.cpp h5_hyperslab.cpp )
    set(TEST_SRCS rho_single_threaded.cpp h5_read.cpp h5_write.cpp )
    set(GISTIC_SRCS gistic.cpp h5_read.cpp h5_write.cpp )
   
    add_executable( rho_single_threaded ${TEST_SRCS} )
    add_executable( rho_omp ${OMP_SRCS} )
    add_executable( pearson_small_block ${PEARSON_SRCS} )
    add_executable( pearson_block_ring ${PEARSON_RING_SRCS} )
    add_executable( pearson_hyperslab ${PEARSON_HYPER_SRCS} )
    add_executable( gistic ${GISTIC_SRCS} )

    target_link_libraries(rho_omp ${HDF5_LIBRARIES})
    target_link_libraries(pearson_small_block ${HDF5_LIBRARIES})
    target_link_libraries(pearson_block_ring ${HDF5_LIBRARIES})
    target_link_libraries(pearson_hyperslab ${HDF5_LIBRARIES})
    target_link_libraries(rho_single_threaded ${HDF5_LIBRARIES})
    target_link_libraries(gistic ${HDF5_LIBRARIES})

endif()
    

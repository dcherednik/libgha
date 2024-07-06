CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=address -fno-omit-frame-pointer")

project(gha)
add_library(gha src/gha.c src/sle.c)

if (NOT GHA_FFT_LIB)

    set(SOURCE_FFT_LIB
        src/3rd/kissfft/kiss_fft.c
        src/3rd/kissfft/tools/kiss_fftr.c
    )
    add_library(gha_fft_impl STATIC ${SOURCE_FFT_LIB})
    set(GHA_FFT_LIB gha_fft_impl)

    target_include_directories(
        gha_fft_impl
        PRIVATE
        src/3rd/kissfft
    )

    target_include_directories(
        gha
        PRIVATE
        src/3rd/kissfft
    )

endif()

if (GHA_USE_DOUBLE_API)
    set_source_files_properties(
        src/gha.c
        src/sle.c
        src/3rd/kissfft/kiss_fft.c
        src/3rd/kissfft/tools/kiss_fftr.c
        test/main.c
        test/dtmf.c
        test/ut.c
        PROPERTIES COMPILE_FLAGS -DGHA_USE_DOUBLE_API
    )
endif()

target_include_directories(
    gha
    PRIVATE
    .
)
target_link_libraries(gha ${GHA_FFT_LIB})

add_definitions("-Wall -O2 -g")

add_executable(main test/main.c)
target_include_directories(
    main
    PRIVATE
    .
)
target_link_libraries(main gha m)

add_executable(dtmf test/dtmf.c)
target_include_directories(
    dtmf
    PRIVATE
    .
)
target_link_libraries(dtmf gha m)

add_executable(ut test/ut.c)
target_include_directories(
    ut
    PRIVATE
    src
    .
)
target_link_libraries(ut gha m)

enable_testing()
add_test(gha_test_simple_1000_0_a main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/1000hz_0.85.pcm 0 1024 0.142476 0.0000 0.850000)
add_test(gha_test_simple_1000_0_b main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/1000hz_0.85.pcm 0 1000 0.142476 0.0000 0.850000)
add_test(gha_test_simple_1000_0_c main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/1000hz_0.85.pcm 0 800 0.142476 0.0000 0.850000)
add_test(gha_test_simple_1000_90_a main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/1000hz_0.85.pcm 11 1024 0.142476 1.5670 0.850000)
add_test(gha_test_simple_1000_90_b main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/1000hz_0.85.pcm 11 1000 0.142476 1.5670 0.850000)
add_test(gha_test_simple_1000_90_c main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/1000hz_0.85.pcm 11 800 0.142476 1.5670 0.850000)


add_test(gha_test_simple_20000_0_a main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/20000hz_0.85.pcm 0 1024 2.8495171 0.0000 0.850000)
add_test(gha_test_simple_20000_0_b main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/20000hz_0.85.pcm 0 500 2.8495171 0.0000 0.850000)
add_test(gha_test_simple_20000_0_c main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/20000hz_0.85.pcm 0 128 2.8495171 0.0000 0.850000)
add_test(gha_test_simple_20000_0_d main ${CMAKE_CURRENT_SOURCE_DIR}/test/data/20000hz_0.85.pcm 0 96 2.8495171 0.0000 0.850000)


add_test(gha_test_dtmf_1 dtmf ${CMAKE_CURRENT_SOURCE_DIR}/test/data/dtmf.pcm 32 256 0.547416 0.201057 0.949511 0.200154)


add_test(ut ut)

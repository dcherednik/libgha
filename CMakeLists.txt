CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=address -fno-omit-frame-pointer")

project(gha)

set(SOURCE_FFT_LIB
    src/3rd/kissfft/kiss_fft.c
    src/3rd/kissfft/tools/kiss_fftr.c
)

set_source_files_properties(
    src/3rd/kissfft/kiss_fft.c
    src/3rd/kissfft/tools/kiss_fftr.c
    PROPERTIES COMPILE_FLAGS -Dkiss_fft_scalar=float
)

add_library(fft_impl STATIC ${SOURCE_FFT_LIB})
target_include_directories(
    fft_impl
    PRIVATE
    src/3rd/kissfft
)

add_library(gha src/gha.c)
target_include_directories(
    gha
    PRIVATE
    src/3rd
    src/3rd/kissfft
    .
)
target_link_libraries(gha fft_impl)

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


enable_testing()
add_test(gha_test_simple_1000_0_a main ../test/data/1000hz_0.85.pcm 0 1024 0.14247585 0.0000 0.850000)
add_test(gha_test_simple_1000_0_b main ../test/data/1000hz_0.85.pcm 0 1000 0.14247585 0.0000 0.850000)
add_test(gha_test_simple_1000_0_c main ../test/data/1000hz_0.85.pcm 0 800 0.14247585 0.0000 0.850000)
add_test(gha_test_simple_1000_90_a main ../test/data/1000hz_0.85.pcm 11 1024 0.14247585 1.5658 0.850000)
add_test(gha_test_simple_1000_90_b main ../test/data/1000hz_0.85.pcm 11 1000 0.14247585 1.5658 0.850000)
add_test(gha_test_simple_1000_90_c main ../test/data/1000hz_0.85.pcm 11 800 0.14247585 1.5658 0.850000)



add_test(gha_test_simple_20000_0_a main ../test/data/20000hz_0.85.pcm 0 1024 2.8495171 0.0000 0.850000)
add_test(gha_test_simple_20000_0_b main ../test/data/20000hz_0.85.pcm 0 500 2.8495171 0.0000 0.850000)
add_test(gha_test_simple_20000_0_c main ../test/data/20000hz_0.85.pcm 0 128 2.8495171 0.0000 0.850000)
add_test(gha_test_simple_20000_0_d main ../test/data/20000hz_0.85.pcm 0 64 2.8495171 0.0000 0.850000)

add_test(gha_test_dtmf_1 dtmf ../test/data/dtmf.pcm 32 256 0.547416 0.201057 0.949511 0.200154)

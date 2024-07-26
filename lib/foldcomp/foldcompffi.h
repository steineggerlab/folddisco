#pragma once
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct {
        float x, y, z;
        char atom[4];
        uint64_t atomIdx;
        char chain;
        char aa[3];
        uint64_t resIdx;
        float bfactor;
    } atom_t;

    void* foldcomp_create();
    atom_t* foldcomp_process(void* instance, const unsigned char* input, size_t length, size_t* atom_count);
    void foldcomp_free( atom_t* output);
    void foldcomp_destroy(void* instance);
#ifdef __cplusplus
}
#endif
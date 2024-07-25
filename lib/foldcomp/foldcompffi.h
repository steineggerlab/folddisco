#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct {
        char aa[3];
        char chain;
        char atom[4];
        float x, y, z;
        float bfactor;
        int atomIdx;
        int resIdx;
    } atom_t;

    void* foldcomp_create();
    atom_t* foldcomp_process(void* instance, const unsigned char* input, size_t length, size_t* atom_count);
    void foldcomp_free( atom_t* output);
    void foldcomp_destroy(void* instance);
#ifdef __cplusplus
}
#endif
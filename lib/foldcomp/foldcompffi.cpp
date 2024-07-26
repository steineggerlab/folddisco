#include "foldcompffi.h"
#include "src/foldcomp.h"

#include <sstream>

extern "C" {
    void* foldcomp_create() {
        return new Foldcomp(); 
    }

    atom_t* foldcomp_process(void* instance, const unsigned char* input, size_t length, size_t* atom_count) {
        // Make a new istream out of const char* and length; binary
        std::istringstream iss(std::string(input, input + length));
        Foldcomp* foldcomp = static_cast<Foldcomp*>(instance);

        std::vector<AtomCoordinate> atoms;
        foldcomp->read(iss);
        foldcomp->decompress(atoms);
        atom_t* atoms_new;
        *atom_count = atoms.size();
        atoms_new = (atom_t*)malloc(sizeof(atom_t) * *atom_count);
        for (size_t i = 0; i < *atom_count; i++) {
            AtomCoordinate atom = atoms[i];
            // Convert string to char array
            atoms_new[i].x = atom.coordinate.x;
            atoms_new[i].y = atom.coordinate.y;
            atoms_new[i].z = atom.coordinate.z;
            if (atom.atom.size() == 4) {
                *atoms_new[i].atom = atom.atom.c_str()[0];
            } else {
                *atoms_new[i].atom = 32;
            }
            *(atoms_new[i].atom + 1) = atom.atom.c_str()[0];
            if (atom.atom.size() == 1) {
                *(atoms_new[i].atom + 2) = 32;
                *(atoms_new[i].atom + 3) = 32;
            } else {
                *(atoms_new[i].atom + 2) = atom.atom.c_str()[1];
                if (atom.atom.size() == 2) {
                    *(atoms_new[i].atom + 3) = 32;
                } else {
                    *(atoms_new[i].atom + 3) = atom.atom.c_str()[2];
                }
            }
            atoms_new[i].atomIdx = atom.atom_index;
            atoms_new[i].chain = atom.chain.c_str()[0];
            *atoms_new[i].aa = atom.residue.c_str()[0];
            *(atoms_new[i].aa + 1) = atom.residue.c_str()[1];
            *(atoms_new[i].aa + 2) = atom.residue.c_str()[2];
            atoms_new[i].resIdx = atom.residue_index;
            atoms_new[i].bfactor = atom.tempFactor;
        }
        return atoms_new;
    }
    
    void foldcomp_free(atom_t* output) {
        free(output);
    }

    void foldcomp_destroy(void* instance) {
        delete static_cast<Foldcomp*>(instance);
    }
}
#include <iostream>
#include "FolddiscoIndex.h"
int main() {
    std::cout << "Hello, World!" << std::endl;
    // initialize index table rust -> c++  FolddiscoIndex(hashSize);
    // for all PDBS
        // rust parses PDB
        // rust extracts the hashes (dereplicate hashes)
        // rust calls the C++ code to count entires countEntries(table, hashes, id, numberOfHashes)
        // (Here we need to make sure that each thread only access one hash range at a time)

    // rust call c++ to allocate memory setup pointers (allocateEntries(table))
    // for all PDBS
        // rust parses PDB
        // rust extracts the hashes + id
        // rust call c++ store entires addEntries(table, hashes, id, numberOfHashes)
        // (Here we need to make sure that each thread only access one hash range at a time)
    // rust call c++ finish index structure finishIndex(table)

    // H1: 0, 1, 92442, 92443
    // lastHash: INT_MAX
    // 0  lastHash[0] = 0                              (offsetarray + 1)
    // 1  lastHash[1] = 1 - 0                          (offsetarray + 1)
    // 2  lastHash[2] = 92442 - 1 > SHORT_MAX + 26906  (offsetarray + 2)
    // 3  lastHash[3] = 92443 - 92442 > SHORT_MAX + 1  (offsetarray + 1)
    return 0;
}
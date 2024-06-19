//
// Created by Martin Steinegger on 6/13/24.
//

#ifndef FOLDDISCOINDEX_FOLDDISCOINDEX_H
#define FOLDDISCOINDEX_FOLDDISCOINDEX_H


class FolddiscoIndex {

    FolddiscoIndex(size_t totalHashes) {
        // allocate memory for the hash table
        offsets = new size_t[totalHashes + 1]; // 2^30
        lastId = new unsigned int[totalHashes];
        // fill lastId with UINT_MAX
        for (size_t i = 0; i < totalHashes; i++) {
            lastId[i] = UINT_MAX;
        }
        this->totalHashes = totalHashes;
    }

public:
    // make sure that only one thread access one hash range at a time
    void countEntries(size_t* hashes, size_t id, size_t numHashes) {
        for (size_t i = 0; i < numHashes; i++) {
            size_t hash = hashes[i];
            unsigned int count = 1;
            if (lastId[hash] == UINT_MAX) {
                count += (id / UCHAR_MAX);
            }
            else {
                count += ((id - lastId[hash]) / UCHAR_MAX);
            }
            __sync_fetch_and_add(&(offsets[hash]), count);
            lastId[hash] = id;
        }
    }
    // threadsafe adding
    void addEntries(size_t* hashes, size_t id, size_t numHashes) {
        for (size_t i = 0; i < numHashes; i++) {
            size_t hash = hashes[i];
            unsigned int count = 1;
            if (lastId[hash] == UINT_MAX) {
                count += (id / UCHAR_MAX);
            }
            else {
                count += ((id - lastId[hash]) / UCHAR_MAX);
            }
            size_t offset = __sync_fetch_and_add(&(offsets[hash]), count);
            lastId[hash] = id;
            // write entry compressed
            entries[offset] = TODO;
        }
    }

    // allocate entires based on offsets and setup points
    void allocateEntries() {
        size_t totalEntries = 0;
        for (size_t i = 0; i < totalHashes; i++) {
            size_t currentEntrySize = offsets[i];
            offsets[i] = totalEntries;
            totalEntries += currentEntrySize;
        }
        entries = new char[totalEntries];
        offsets[totalHashes] = totalEntries;
        for (size_t i = 0; i < totalHashes; i++) {
            lastId[i] = UINT_MAX;
        }
    }

    // finish index structure, turn offset table into back to start positions of entires
    void finishIndex() {
        for (size_t i = totalHashes; i > 0; i--) {
            offsets[i] = offsets[i - 1];
        }
        offsets[0] = 0;
    }

    // get list of DB sequences containing this k-mer
    inline char* getEntries(size_t hash, size_t* matchedListSize) {
        const ptrdiff_t diff = offsets[hash + 1] - offsets[hash];
        *matchedListSize = static_cast<size_t>(diff);
        return (entries + offsets[hash]);
    }


private:
    size_t* offsets;
    unsigned int* lastId;
    size_t totalHashes;
    char* entries;
};


#endif //FOLDDISCOINDEX_FOLDDISCOINDEX_H
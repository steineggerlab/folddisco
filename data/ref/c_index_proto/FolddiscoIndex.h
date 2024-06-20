#ifndef FOLDDISCOINDEX_FOLDDISCOINDEX_H
#define FOLDDISCOINDEX_FOLDDISCOINDEX_H

#include <stddef.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

// Define the FolddiscoIndex struct
typedef struct {
    size_t* offsets;
    unsigned int* lastId;
    size_t totalHashes;
    char* entries;
} FolddiscoIndex;

// Function to create a new FolddiscoIndex
FolddiscoIndex* create_folddisco_index(size_t totalHashes) {
    FolddiscoIndex* index = (FolddiscoIndex*)malloc(sizeof(FolddiscoIndex));
    if (!index) return NULL;

    index->offsets = (size_t*)malloc((totalHashes + 1) * sizeof(size_t));
    index->lastId = (unsigned int*)malloc(totalHashes * sizeof(unsigned int));

    if (!index->offsets || !index->lastId) {
        free(index->offsets);
        free(index->lastId);
        free(index);
        return NULL;
    }

    for (size_t i = 0; i < totalHashes; i++) {
        index->lastId[i] = UINT_MAX;
    }

    index->totalHashes = totalHashes;
    return index;
}

// Function to count entries
void count_entries(FolddiscoIndex* index, size_t* hashes, size_t id, size_t numHashes) {
    for (size_t i = 0; i < numHashes; i++) {
        size_t hash = hashes[i];
        unsigned int count = 1;
        if (index->lastId[hash] == UINT_MAX) {
            count += (id / UCHAR_MAX);
        }
        else {
            count += ((id - index->lastId[hash]) / UCHAR_MAX);
        }
        __sync_fetch_and_add(&(index->offsets[hash]), count);
        index->lastId[hash] = id;
    }
}

// Function to add entries
void add_entries(FolddiscoIndex* index, size_t* hashes, size_t id, size_t numHashes) {
    for (size_t i = 0; i < numHashes; i++) {
        size_t hash = hashes[i];
        unsigned int count = 1;
        if (index->lastId[hash] == UINT_MAX) {
            count += (id / UCHAR_MAX);
        }
        else {
            count += ((id - index->lastId[hash]) / UCHAR_MAX);
        }
        size_t offset = __sync_fetch_and_add(&(index->offsets[hash]), count);
        index->lastId[hash] = id;
        // TODO: Implement entry compression and assignment
        // index->entries[offset] = TODO;
    }
}

// Function to allocate entries
void allocate_entries(FolddiscoIndex* index) {
    size_t totalEntries = 0;
    for (size_t i = 0; i < index->totalHashes; i++) {
        size_t currentEntrySize = index->offsets[i];
        index->offsets[i] = totalEntries;
        totalEntries += currentEntrySize;
    }
    index->entries = (char*)malloc(totalEntries * sizeof(char));
    index->offsets[index->totalHashes] = totalEntries;
    for (size_t i = 0; i < index->totalHashes; i++) {
        index->lastId[i] = UINT_MAX;
    }
}

// Function to finish index
void finish_index(FolddiscoIndex* index) {
    for (size_t i = index->totalHashes; i > 0; i--) {
        index->offsets[i] = index->offsets[i - 1];
    }
    index->offsets[0] = 0;
}

// Function to get entries
char* get_entries(FolddiscoIndex* index, size_t hash, size_t* matchedListSize) {
    ptrdiff_t diff = index->offsets[hash + 1] - index->offsets[hash];
    *matchedListSize = (size_t)diff;
    return (index->entries + index->offsets[hash]);
}

// Function to free the FolddiscoIndex
void free_folddisco_index(FolddiscoIndex* index) {
    free(index->offsets);
    free(index->lastId);
    free(index->entries);
    free(index);
}

#endif // FOLDDISCOINDEX_FOLDDISCOINDEX_H
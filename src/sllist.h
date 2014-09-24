#ifndef SLLIST_H_
#define SLLIST_H_

#include "utilities.h"

// abstraction for singly linked lists
typedef struct SinglyLinkedList_st {
    struct SinglyLinkedList_st* next;
}SinglyLinkedList;

#define SllAddHead(plist, node) ((node)->next = *(plist), *(plist) = node)

// return the number of elements in the list
int SllCount(const void* const list);

// reverse order of  a list
void SllReverse(void* plist);

// remove the item from the list
void* SllRemove(void* plist, void* const node);

// free the list and set the pointer to the list to be NULL
void SllFreeList(void* plist);

// sort the linked list with qsort and a temporary array
void SllSort(void* plist,
            int(*compare)(const void* const elem1, const void* const elem2));

#endif  // SINGLE_LINK_LIST_H_

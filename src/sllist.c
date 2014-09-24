#include "sllist.h"

// return the number of elements in the list
int SllCount(const void* const list) {
    int count = 0;
    SinglyLinkedList* iter = (SinglyLinkedList*) list;

    for (; iter; ++count, iter = iter->next) {}
    return count;
}

// reverse order of  a list
void SllReverse(void* plist) {
    SinglyLinkedList** ppt = (SinglyLinkedList**)plist;
    SinglyLinkedList* newlist = NULL;
    SinglyLinkedList* el = NULL;
    SinglyLinkedList* next = NULL;

    next = *ppt;
    while (next != NULL) {
        el = next;
        next = el->next;
        el->next = newlist;
        newlist = el;
    }
    *ppt = newlist;
}

// remove the item from the list
void* SllRemove(void* plist, void* const node) {
    SinglyLinkedList* iter = *((SinglyLinkedList**) plist);
    SinglyLinkedList* t = (SinglyLinkedList*) node;
    SinglyLinkedList* pt = (SinglyLinkedList*) node;

    if (iter == t) {
        *(SinglyLinkedList**)plist = t->next;
        return t;
    }
    for (; iter && (iter != t); pt = iter, iter = iter->next) {}
    pt->next = t->next;
    return t;
}

// free the list and set the pointer to the list to be NULL
void SllFreeList(void* plist) {
    SinglyLinkedList** ppt = (SinglyLinkedList**)plist;
    SinglyLinkedList* next = *ppt;
    SinglyLinkedList* el = NULL;

    while (next != NULL) {
        el = next;
        next = el->next;
        Ckfree((char*)el);
    }
    *ppt = NULL;
}

// sort the linked list with qsort and a temporary array
void SllSort(void* plist,
            int(*compare)(const void* const elem1, const void* const elem2)) {
    SinglyLinkedList** pl = (SinglyLinkedList**)plist;
    SinglyLinkedList* list = *pl;

    int count = SllCount(list);
    if (count > 1) {
        SinglyLinkedList* el = NULL;
        SinglyLinkedList** array = NULL;
        int i = 0;

        array = CkallocOrDie(count * sizeof(*array));
        for (el = list, i = 0; el != NULL; el = el->next, i++) {
            array[i] = el;
        }
        qsort(array, count, sizeof(array[0]), compare);
        list = NULL;
        for (i = 0; i < count; i++) {
            array[i]->next = list;
            list = array[i];
        }
        Ckfree(array);
        SllReverse(&list);
        *pl = list;
    }
}

#ifndef TSORT_H_INCLUDED
#define TSORT_H_INCLUDED

struct pthread_args {
    void*               dataPointer;
    unsigned long long  dataLength;
    int                 dataSize;
    int                 (*comparator)(const void*, const void*);
};

int bisect(void* massive[], int massive_length, const void* value, int (*comparator)(const void*, const void*));
int compair_parts(const void* pointer_a, const void* pointer_b);
void* merge(struct pthread_args* args, int (*comparator)(const void*, const void*), int threads_max);
void* tsort(void* dataPointer, long long dataLength, int dataSize, int (*comparator)(const void*, const void*), int threads_max);
void* qsort_thread(void* args_pointer);

#endif // TSORT_H_INCLUDED

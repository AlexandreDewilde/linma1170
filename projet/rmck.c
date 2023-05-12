#include "rmck.h"

typedef struct node {
    struct node *next;
    int val;
} node;

typedef struct {
    node *start;
    node *end;
    int n;
} Queue;

typedef struct {
    int *array;
    int n;
    int N;
} vector;

vector *create_vector() {
    vector *vec = malloc(sizeof(vector));    
    vec->n = 0;
    vec->N = 2;
    vec->array = malloc(sizeof(int) * 2);
    return vec;
}

void delete_vector(vector *vec) {
    free(vec->array);
    free(vec);
}

void push_vector(vector *vec, int val) {
    if (vec->N <= vec->n) {
        vec->N = vec->N * 2;
        vec->array = realloc(vec->array, vec->N * sizeof(int));
    }
    vec->array[vec->n++] = val; 
}

int pop_vector(vector *vec) {
    return vec->array[vec->n-- - 1];
}

int get_vector(vector *vec, int idx) {
    return vec->array[idx];
}

void swap_vector(vector *vec, int idx_1, int idx_2) {
    int temp = vec->array[idx_1];
    vec->array[idx_1] = vec->array[idx_2];
    vec->array[idx_2] = temp;
}


Queue *create_queue() {
    Queue *q = malloc(sizeof(Queue));
    q->n = 0;
    q->start = q->end = NULL;
    return q;
}

void delete_queue(Queue *q) {
    node *current = q->start;
    for (int i = 0; i < q->n; i++) {
        node *prev = current;
        current = current->next;
        free(prev);
    }
    free(q);
}

void push_queue(Queue *q, int val) {
    q->n++;
    if (q->start == NULL) {
        q->start = q->end = malloc(sizeof(node));
        q->start->next = NULL;
        q->start->val = val;
        return;
    }
    q->end->next = malloc(sizeof(node));
    q->end = q->end->next;
    q->end->val = val;
    q->end->next = NULL;
}

int pop_queue(Queue *q) {
    q->n--;
    int val = q->start->val;
    node *prev = q->start;
    q->start = q->start->next;
    free(prev);
    return val;
}

static vector **g;

static int compare_adjency_list(const void * g1, const void * g2) {
    int idx1 = *((int*) g1); int idx2 = *((int*) g2);
    return g[idx1]->n - g[idx2]->n;
}

static int compare_int(const void *int1, const void *int2) {
    int n1 = *((const int*) int1); int n2 = *((const int*) int2);    
    return n1 - n2;
}

void rmck(int *perm, Triplet *triplets, int n_triplets, int n_nodes) {
    int max_coords = 2 * n_nodes;
    
    // Create adjency lists
    g = malloc(max_coords * sizeof(vector));
    for (int i = 0; i < max_coords; i++) {
        g[i] = create_vector();
    }

    for (int i = 0; i < n_triplets; i++) {
        push_vector(g[triplets[i].i], triplets[i].j);
    }

    int *sortedDegrees = malloc(max_coords * sizeof(int));
    for (int i = 0; i < max_coords; i++) {
        qsort(g[i]->array, g[i]->n, sizeof(int), compare_int);
        sortedDegrees[i] = i;
    }
    qsort(sortedDegrees, max_coords, sizeof(int), compare_adjency_list);

    char *visited = malloc(max_coords);
    memset(visited, 0, max_coords);

    int RIdx = max_coords - 1;

    for (int i = 0; i < max_coords; i++) {
        int x = sortedDegrees[i];
        if (visited[x]) {
            continue;
        }
        
        Queue *q = create_queue();
        push_queue(q, x);
        while (q->n) {
            int x = pop_queue(q);
            if (visited[x]) {
                continue;
            }
            perm[RIdx--] = x;
            visited[x] = 1;
            for (int i = 0; i < g[x]->n; i++) {
                if (visited[get_vector(g[x], i)]) continue;
                push_queue(q, get_vector(g[x], i));
            }
        }
        delete_queue(q);
    }

    for (int i = 0; i < max_coords; i++) {
        delete_vector(g[i]);
    }
    free(sortedDegrees);
    free(g);
    free(visited);
}

static int *global_degrees;
static int comp_degrees(const void * g1, const void * g2) {
    int x = *((int*)g1); int y = *((int*)g2);
    return global_degrees[x] - global_degrees[y];
}

void rmck_matrix(int *perm, Matrix *A) {
    int max_coords = A->m;

    int *degrees = calloc(A->m, sizeof(int));
    global_degrees = degrees;
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (fabs(A->a[i][j]) > 1e-20) {
                ++degrees[i];
            }
        }
    }
    
    int *sortedDegrees = malloc(sizeof(int)*A->m);
    for (int i = 0; i < A->m; i++) {
        sortedDegrees[i] = i;
    }
    qsort(sortedDegrees, A->m, sizeof(int), comp_degrees);

    char *visited = malloc(A->m);
    memset(visited, 0, A->m);

    int RIdx = max_coords - 1;

    for (int i = 0; i < max_coords; i++) {
        int x = sortedDegrees[i];
        if (visited[x]) {
            continue;
        }
        
        Queue *q = create_queue();
        push_queue(q, x);
        while (q->n) {
            int x = pop_queue(q);
            if (visited[x]) {
                continue;
            }
            perm[RIdx--] = x;
            visited[x] = 1;
            vector *vec = create_vector();
            for (int i = 0; i < A->m; i++) {
                if (visited[i] || fabs(A->a[x][i]) <= 1e-20) continue;
                push_vector(vec, i);
            }
            qsort(vec->array, vec->n, sizeof(int), comp_degrees);
            for (int i = 0; i < vec->n; i++) {
                push_queue(q, vec->array[i]);
            }
            delete_vector(vec);
        }
        delete_queue(q);
    }
    free(sortedDegrees);
    free(visited);
}
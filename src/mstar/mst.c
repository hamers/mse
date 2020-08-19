#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "regularization.h"

int **indexlist;

double GBSTOL;
int MAXPART;

double *errorbuf;
double *resultbuf;
double *y;
double *yscal;
double *maxerror;
double **error, **result;
double **gbsS;

//double mst_CoM_R[3];
//double mst_CoM_V[3];

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
MPI_Comm MPI_COMM_GBS_GROUP;
#endif


int Ntask;
int ThisTask;

int NumGbsGroup;
int NumTaskPerGbsGroup;
int ThisGbsGroup;
int ThisTask_in_GbsGroup;

int NumTaskPerSortGroup;
int NumSortGroup;
int ThisSortGroup;
int ThisTask_in_SortGroup;

int ok_steps;
int failed_steps;

struct ToDoList ComputationToDoList;
struct OctreeNode *OctreeRootNode;
struct OctreeNode **ParticleInNode;


double tprim1, tprim2, tprim3, tprim4, tprim5;
double time_coord;
double time_force;
double time_comm;
double time_mst;
double time_gbs;
double time_gbscomm;
clock_t wtime1, wtime0;
clock_t msttime1, msttime0;

// The functions begin here...

void my_barrier(void) {
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void die(void) {

    my_barrier();
    if (ThisTask == 0)
        printf("The code is dead because 'die()' was called..\n");
#ifdef PARALLEL
    MPI_Finalize();
#endif
    exit(0);
}

int cmp_weight_index2(const void *a, const void *b) {
    struct WeightIndex2 *a1 = (struct WeightIndex2 *)a;
    struct WeightIndex2 *a2 = (struct WeightIndex2 *)b;
    if ((*a1).weight > (*a2).weight)
        return 1;
    else if ((*a1).weight < (*a2).weight)
        return -1;
    else
        return 0;
}
int cmp_weight_index(const void *a, const void *b) {
    struct WeightIndex *a1 = (struct WeightIndex *)a;
    struct WeightIndex *a2 = (struct WeightIndex *)b;
    if ((*a1).weight > (*a2).weight)
        return 1;
    else if ((*a1).weight < (*a2).weight)
        return -1;
    else
        return 0;
}
int cmp(const void *a, const void *b) {
    struct GraphEdge *a1 = (struct GraphEdge *)a;
    struct GraphEdge *a2 = (struct GraphEdge *)b;
    if ((*a1).weight > (*a2).weight)
        return 1;
    else if ((*a1).weight < (*a2).weight)
        return -1;
    else
        return 0;
}
int cmp_int(const void *a, const void *b) {
    int int_a = *((int *)a);
    int int_b = *((int *)b);
    if (int_a == int_b)
        return 0;
    else if (int_a < int_b)
        return -1;
    else
        return 1;
}

#ifdef PARALLEL
void initialize_node_comm(void) {

    MPI_Group WORLD_GROUP, GBS_GROUP;
    MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GROUP);
    int *myranks;

    ThisGbsGroup = (int)floor(ThisTask / NumTaskPerGbsGroup);
    ThisTask_in_GbsGroup = ThisTask % NumTaskPerGbsGroup;

    myranks = calloc(NumTaskPerGbsGroup, sizeof(int));

    for (int i = 0; i < NumTaskPerGbsGroup; i++) {
        myranks[i] = ThisTask - ThisTask_in_GbsGroup + i;
    }
    for (int i = 0; i < NumGbsGroup; i++) {
        if (i == ThisGbsGroup) {
            MPI_Group_incl(WORLD_GROUP, NumTaskPerGbsGroup, myranks,
                           &GBS_GROUP);
        }
    }

    MPI_Comm_create(MPI_COMM_WORLD, GBS_GROUP, &MPI_COMM_GBS_GROUP);
    MPI_Group_free(&GBS_GROUP);
    free(myranks);

    MPI_Group_free(&WORLD_GROUP);
    MPI_Comm_size(MPI_COMM_GBS_GROUP, &NumTaskPerGbsGroup);
    MPI_Comm_rank(MPI_COMM_GBS_GROUP, &ThisTask_in_GbsGroup);
}
#endif

void initialize_mpi_or_serial(void) {

// check for reasonable input
#if defined(PARALLEL) && defined(SERIAL)
    printf(
        "Please select EITHER parallel OR serial. Currently both selected.\n");
    exit(0);
#endif
#if !defined(PARALLEL) && !defined(SERIAL)
    printf(
        "Please select EITHER parallel OR serial. Currently none selected.\n");
    exit(0);
#endif
#if defined(PRIM) && defined(PRIM_DIVIDE_CONQUER)
    printf(
        "Please select EITHER prim OR prim_divide_conquer. Currently both "
        "selected.\n");
    exit(0);
#endif
#if !defined(PRIM) && !defined(PRIM_DIVIDE_CONQUER)
    printf(
        "Please select EITHER prim OR prim_divide_conquer. Currently none "
        "selected.\n");
    exit(0);
#endif

// Init MPI if needed
#ifdef PARALLEL
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntask);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

    NumGbsGroup = 1;
    NumTaskPerGbsGroup = (int)Ntask / NumGbsGroup;

    if (NumGbsGroup * NumTaskPerGbsGroup != Ntask) {
        if (ThisTask == 0) {
            printf(
                "NumGbsGroup*NumTaskPerGbsGroup != Ntask.        NumGbsGroup = "
                "%d, NumTaskPerGbsGroup = %d.\n",
                NumGbsGroup, NumTaskPerGbsGroup);
            exit(0);
        }
    }
    if (NumTaskPerGbsGroup > 1) {
        initialize_node_comm();
    } else {
        ThisGbsGroup = ThisTask;
        ThisTask_in_GbsGroup = 0;
    }
#else
    NumGbsGroup = 1;
    Ntask = 1;
    ThisTask = 0;
    ThisGbsGroup = 0;
    ThisTask_in_GbsGroup = 0;
    NumTaskPerGbsGroup = 1;
#endif
}

int check_relative_proximity_ND_2(int v1, int v2, struct RegularizedRegion *R,
                                  int *d, int *path, int *sign) {

    int diff = abs(R->Vertex[v1].level - R->Vertex[v2].level);

    if (diff == 1) {
        if (R->Vertex[v2].parent == v1) {
            *d = 1;
            sign[0] = -1;
            path[0] = R->Vertex[v2].edgetoparent;
            return 1;
        } else if (R->Vertex[v1].parent == v2) {
            *d = 1;
            sign[0] = +1;
            path[0] = R->Vertex[v1].edgetoparent;
            return 1;
        }
        return 0;
    } else if (diff == 2) {
        if (R->Vertex[R->Vertex[v2].parent].parent == v1) {
            *d = 2;
            sign[0] = -1;
            sign[1] = -1;
            path[0] = R->Vertex[v2].edgetoparent;
            path[1] = R->Vertex[R->Vertex[v2].parent].edgetoparent;
            return 1;
        } else if (R->Vertex[R->Vertex[v1].parent].parent == v2) {
            *d = 2;
            sign[0] = +1;
            sign[1] = +1;
            path[0] = R->Vertex[v1].edgetoparent;
            path[1] = R->Vertex[R->Vertex[v1].parent].edgetoparent;
            return 1;
        }
        return 0;
    } else if (diff == 0) {
        if (R->Vertex[v1].parent == R->Vertex[v2].parent) {
            *d = 2;
            sign[0] = -1;
            sign[1] = +1;
            path[0] = R->Vertex[v2].edgetoparent;
            path[1] = R->Vertex[v1].edgetoparent;
            return 1;
        }
        return 0;
    }
    return 0;
}

int check_relative_proximity(int v1, int v2, const int Nd,
                             struct RegularizedRegion *R, int *d, int *path,
                             int *sign) {

    for (int i = 0; i < Nd; i++) {
        path[i] = -666;
        sign[i] = 0;
    }

    int proximity = 0;
    if (R->Vertex[v1].level == R->Vertex[v2].level) {

        if (R->Vertex[v1].level == 0) {
            printf("This level stuff should not happen\n");
            printf("thistask: %d --- v1 v2 %d %d ---levels: %d %d\n", ThisTask,
                   R->Vertex[v1].id, R->Vertex[v2].id, R->Vertex[v1].level,
                   R->Vertex[v2].level);
            exit(0);
        }

        int halfNd = (int)floor(0.5 * Nd);
        int *pathend = calloc(Nd, sizeof(int));

        for (int i = 0; i < halfNd; i++) {
            sign[i] = +1;
            path[i] = R->Vertex[v1].edgetoparent;
            pathend[i] = R->Vertex[v2].edgetoparent;
            v1 = R->Vertex[v1].parent;
            v2 = R->Vertex[v2].parent;
            if (v1 == v2) {
                *d = 2 * (i + 1);
                int c = 0;
                for (int m = i + 1; m < (*d); m++) {
                    sign[m] = -1;
                    path[m] = pathend[i - c];
                    c++;
                }
                proximity = 1;
                break;
            }
        }

        free(pathend);

    } else {
        *d = abs(R->Vertex[v1].level - R->Vertex[v2].level);
        if (*d <= Nd) {
            int lo, hi, thissign;
            if (R->Vertex[v1].level > R->Vertex[v2].level) {
                lo = v2;
                hi = v1;
                thissign = +1;
            } else {
                lo = v1;
                hi = v2;
                thissign = -1;
            }
            for (int i = 0; i < *d; i++) {
                path[i] = R->Vertex[hi].edgetoparent;
                sign[i] = thissign;
                hi = R->Vertex[hi].parent;
            }
            if (lo == hi) { proximity = 1; }
        }
    }

    return proximity;
}

void loop_scheduling_N2(int Nrow, int Nt, int *istart, int *jstart,
                        int *ThisBlock, int *CumNum, int task) {
    int Np = (Nrow * (Nrow - 1)) / 2;
    int block = (int)floor(Np / Nt);
    int limit = (block + 1) * Nt - Np;
    int ThisCumNum;
    if (task < limit) {
        ThisCumNum = task * block;
        *ThisBlock = block;
    } else {
        ThisCumNum = limit * block + (task - limit) * (block + 1);
        *ThisBlock = block + 1;
    }

    double drow =
        1.0 * Nrow - 0.5 -
        0.5 * sqrt((1.0 - 2.0 * Nrow) * (1.0 - 2.0 * Nrow) - 8.0 * ThisCumNum);
    int row = (int)floor(drow);
    int rowstart =
        (int)round(-0.5 * (1.0 * row + 1.0) * (row - 2.0 * Nrow) - 1.0 * Nrow);
    *istart = row;
    *jstart = ThisCumNum - rowstart + row + 1;
    *CumNum = ThisCumNum;
}

void get_v1_v2_fast(struct RegularizedRegion *R, int i, int *v1, int *v2) {
    *v1 = R->EdgeInMST[i].vertex1;
    *v2 = R->EdgeInMST[i].vertex2;
}

void get_v1_v2(int N, int target, int *v1, int *v2) {

    int final_element_on_row = 0;
    int elements_on_this_row;

    for (int row = 0; row < N - 1; row++) {
        if (row == 0) {
            elements_on_this_row = N - 2 - row;
        } else {
            elements_on_this_row = N - 1 - row;
        }
        final_element_on_row += elements_on_this_row;
        if (target <= final_element_on_row) {
            *v1 = row;
            int diff = final_element_on_row - target;
            *v2 = N - 1 - diff;
            break;
        }
    }
}

int sign(double x) {
    if (x >= 0) return 1;
    return -1;
}

double get_rmax(struct RegularizedRegion *R) {
    double rmax2 = 0;
    for (int i = 0; i < R->NumVertex; i++) {
        double ds;
        double r2 = 0;
        for (int k = 0; k < 3; k++) {
            ds = R->Pos[3 * i + k];
            r2 += ds * ds;
        }
        if (r2 > rmax2) { rmax2 = r2; }
    }
    return sqrt(rmax2);
}

void get_triplet(int index, int *triplet) {
    if (index >= 100) {
        triplet[2] = 1;
        index -= 100;
        if (index >= 10) {
            triplet[1] = 1;
            index -= 10;
            if (index == 1) {
                triplet[0] = 1;
            } else {
                triplet[0] = -1;
            }
        } else {
            triplet[1] = 0;
            if (index == 1) {
                triplet[0] = 1;
            } else {
                triplet[0] = -1;
            }
        }
    } else {
        triplet[2] = 0;
        if (index >= 10) {
            triplet[1] = 1;
            index -= 10;
            if (index == 1) {
                triplet[0] = 1;
            } else {
                triplet[0] = -1;
            }
        } else {
            triplet[1] = 0;
            if (index == 1) {
                triplet[0] = 1;
            } else {
                triplet[0] = -1;
            }
        }
    }
}

int get_triplet_index(int num) {
    if (num == 1) return 1;
    return 0;
}

int get_octant(int triplet[3]) {
    if (triplet[0] == 1) {
        if (triplet[1] == 1) {
            if (triplet[2] == -1) {
                return 6;
            } else {
                return 7;
            }
        } else {
            if (triplet[2] == -1) {
                return 4;
            } else {
                return 5;
            }
        }
    } else {
        if (triplet[1] == 1) {
            if (triplet[2] == -1) {
                return 2;
            } else {
                return 3;
            }
        } else {
            if (triplet[2] == -1) {
                return 0;
            } else {
                return 1;
            }
        }
    }
    return -1;
}

void get_triplet_from_octant(int octant, int *triplet) {

    if (octant == 0) {
        triplet[0] = -1;
        triplet[1] = -1;
        triplet[2] = -1;
        return;
    }
    if (octant == 1) {
        triplet[0] = -1;
        triplet[1] = -1;
        triplet[2] = +1;
        return;
    }
    if (octant == 2) {
        triplet[0] = -1;
        triplet[1] = +1;
        triplet[2] = -1;
        return;
    }
    if (octant == 3) {
        triplet[0] = -1;
        triplet[1] = +1;
        triplet[2] = +1;
        return;
    }
    if (octant == 4) {
        triplet[0] = +1;
        triplet[1] = -1;
        triplet[2] = -1;
        return;
    }
    if (octant == 5) {
        triplet[0] = +1;
        triplet[1] = -1;
        triplet[2] = +1;
        return;
    }
    if (octant == 6) {
        triplet[0] = +1;
        triplet[1] = +1;
        triplet[2] = -1;
        return;
    }
    if (octant == 7) {
        triplet[0] = +1;
        triplet[1] = +1;
        triplet[2] = +1;
        return;
    }
    return;
}

void list_octree_level(struct OctreeNode *Node, int *LeafSizes,
                       int **LeafParticles, int *Nleafs) {

    if (Node->IsLeaf) {
        LeafSizes[*Nleafs] = Node->Npart;
        LeafParticles[*Nleafs] = calloc(Node->Npart, sizeof(int));
        for (int i = 0; i < Node->Npart; i++) {
            LeafParticles[*Nleafs][i] = Node->ObjectList[i];
        }
        (*Nleafs)++;
    } else {
        for (int oct = 0; oct < 8; oct++) {
            if (Node->Children[oct]->ChildExists) {
                list_octree_level(Node->Children[oct], LeafSizes, LeafParticles,
                                  Nleafs);
            }
        }
    }
}

void GetNumberOfLeafs(struct OctreeNode *Node, int *Nleaf) {
    if (Node->IsLeaf) {
        (*Nleaf)++;
    } else {
        for (int oct = 0; oct < 8; oct++) {
            if (Node->Children[oct]->ChildExists) {
                GetNumberOfLeafs(Node->Children[oct], Nleaf);
            }
        }
    }
}

void free_octree(struct OctreeNode *Node) {
    for (int oct = 0; oct < 8; oct++) {
        if (Node->Children[oct]->ChildExists) {
            free_octree(Node->Children[oct]);
        }
    }
    free(Node->ObjectList);
    free(Node);
}

void allocate_new_octree_node(struct OctreeNode *Node, int oct, int *label) {

    Node->ChildExists[oct] = 1;
    Node->Children[oct] = calloc(1, sizeof(struct OctreeNode));
    Node->Children[oct]->Parent = Node;
    Node->Children[oct]->HalfSize = Node->HalfSize / 2;
    Node->Children[oct]->Npart = 0;
    Node->Children[oct]->level = Node->level + 1;
    Node->Children[oct]->IsLeaf = 0;
    Node->Children[oct]->ID = *label;

    for (int i = 0; i < 8; i++) {
        Node->Children[oct]->ChildExists[i] = 0;
    }

    (*label)++;
}

void octree(struct RegularizedRegion *R, struct OctreeNode *Node,
            int MaxLeafSize, int *label) {

    int NpartOctant[8] = {0};
    int triplet[3];
    for (int i = 0; i < Node->Npart; i++) {
        int j = Node->ObjectList[i];
        for (int k = 0; k < 3; k++) {
            triplet[k] = sign(R->Pos[3 * j + k] - Node->Center[k]);
        }
        int oct = get_octant(triplet);

        if (NpartOctant[oct] == 0) {
            allocate_new_octree_node(Node, oct, label);
        }

        indexlist[oct][NpartOctant[oct]] = j;

        NpartOctant[oct]++;
        Node->Children[oct]->Npart++;
    }

    for (int oct = 0; oct < 8; oct++) {

        if (NpartOctant[oct] == 0) { continue; }

        Node->Children[oct]->ObjectList =
            calloc(Node->Children[oct]->Npart, sizeof(int));

        for (int i = 0; i < Node->Children[oct]->Npart; i++) {
            Node->Children[oct]->ObjectList[i] = indexlist[oct][i];
        }

        get_triplet_from_octant(oct, triplet);

        for (int k = 0; k < 3; k++) {
            Node->Children[oct]->Center[k] =
                Node->Center[k] + triplet[k] * Node->Children[oct]->HalfSize;
        }

        for (int i = 0; i < Node->Children[oct]->Npart; i++) {
            ParticleInNode[Node->Children[oct]->ObjectList[i]] =
                Node->Children[oct];
        }

        if (Node->Children[oct]->Npart == 1) {
            Node->Children[oct]->IsLeaf = 1;
        }
    }

    for (int oct = 0; oct < 8; oct++) {
        if (Node->ChildExists[oct]) {
            if (Node->Children[oct]->Npart > MaxLeafSize) {
                octree(R, Node->Children[oct], MaxLeafSize, label);
            } else {
                Node->Children[oct]->IsLeaf = 1;
            }
        }
    }
}

void spatialpartition(struct RegularizedRegion *R, int MaxLeafSize) {

    double rmax = get_rmax(R);

    OctreeRootNode = calloc(1, sizeof(struct OctreeNode));

    OctreeRootNode->ID = 0;
    OctreeRootNode->Parent = NULL;
    OctreeRootNode->HalfSize = 1.001 * rmax;
    OctreeRootNode->Npart = R->NumVertex;
    OctreeRootNode->ObjectList = calloc(R->NumVertex, sizeof(int));
    OctreeRootNode->level = 0;
    for (int i = 0; i < 3; i++) {
        OctreeRootNode->Center[i] = 0.0;
    }
    for (int i = 0; i < 8; i++) {
        OctreeRootNode->ChildExists[i] = 0;
    }

    for (int i = 0; i < R->NumVertex; i++) {
        OctreeRootNode->ObjectList[i] = i;
    }

    int label = 1;
    octree(R, OctreeRootNode, MaxLeafSize, &label);
}

int get_longest(double w[3], int id[3]) {

    double L = -1;
    int ind = -1;
    for (int m = 0; m < 3; m++) {
        if (w[m] > L) {
            L = w[m];
            ind = m;
        }
    }
    return id[ind];
}

int get_edge(int Nleafs, int va, int vb) {

    int N = Nleafs - 2;
    int lo, hi;
    lo = va, hi = vb;
    if (va > vb) {
        lo = vb;
        hi = va;
    }
    int index = hi - 1;
    for (int p = 0; p < lo; p++) {
        index += N;
        N--;
    }
    return index;
}

void SolveMetaMST(struct RegularizedRegion *R, int Nleafs, int *Npart,
                  int **Particles, int *NumEdge) {

    int Nedge = Nleafs * (Nleafs - 1) / 2;
    struct WeightIndex2 *W = calloc(Nedge, sizeof(struct WeightIndex2));
    int *InMst = calloc(Nleafs, sizeof(int));

    int CumNumEdge = 0, istart = 0;
#ifdef PARALLEL
    int ThisBlock, jstart, loop;
    loop_scheduling_N2(Nleafs, Ntask, &istart, &jstart, &ThisBlock, &CumNumEdge,
                       ThisTask);
#endif

    struct OctreeNode *Node_i, *Node_j;
    int counter = 0, ind1 = -1, ind2 = -1;

    int *ind1s = calloc(Nedge, sizeof(int));
    int *ind2s = calloc(Nedge, sizeof(int));
    double *weights = calloc(Nedge, sizeof(double));

    counter = 0;
    for (int i = istart; i < Nleafs; i++) {

        Node_i = ParticleInNode[Particles[i][0]];

#ifdef PARALLEL
        if (i > istart) { jstart = i + 1; }
        for (int j = jstart; j < Nleafs; j++) {
#else
        for (int j = i + 1; j < Nleafs; j++) {
#endif

            Node_j = ParticleInNode[Particles[j][0]];

            double r2 = 0;
            for (int k = 0; k < 3; k++) {
                double ds = Node_i->Center[k] - Node_j->Center[k];
                r2 = ds * ds;
            }
            double NodeSeparation = sqrt(r2);
            double NodeSize =
                1.01 * sqrt(2) *
                (Node_i->HalfSize + Node_j->HalfSize); // with safety buffer

            if (NodeSeparation <= NodeSize) {

                double MinWeight = 1e10;
                for (int q = 0; q < Npart[i]; q++) {
                    for (int p = 0; p < Npart[j]; p++) {
                        double r2 = 0;
                        int i1 = Particles[i][q];
                        int i2 = Particles[j][p];
                        for (int k = 0; k < 3; k++) {
                            double ds = R->Pos[3 * i1 + k] - R->Pos[3 * i2 + k];
                            r2 += ds * ds;
                        }
                        double r = sqrt(r2);
                        if (r < MinWeight) {
                            MinWeight = r;
                            ind1 = i1;
                            ind2 = i2;
                        }
                    }
                }

            } else {

                double MinWeight = 1e10;
                for (int q = 0; q < Npart[i]; q++) {
                    int i1 = Particles[i][q];
                    double r2 = 0;
                    for (int k = 0; k < 3; k++) {
                        double ds = R->Pos[3 * i1 + k] - Node_j->Center[k];
                        r2 += ds * ds;
                    }
                    double r = sqrt(r2);
                    if (r < MinWeight) {
                        MinWeight = r;
                        ind1 = i1;
                    }
                }
                MinWeight = 1e10;
                for (int p = 0; p < Npart[j]; p++) {
                    int i2 = Particles[j][p];
                    double r2 = 0;
                    for (int k = 0; k < 3; k++) {
                        double ds = R->Pos[3 * i2 + k] - Node_i->Center[k];
                        r2 += ds * ds;
                    }
                    double r = sqrt(r2);
                    if (r < MinWeight) {
                        MinWeight = r;
                        ind2 = i2;
                    }
                }
            }

            r2 = 0;
            for (int k = 0; k < 3; k++) {
                double ds = R->Pos[3 * ind1 + k] - R->Pos[3 * ind2 + k];
                r2 += ds * ds;
            }
            double r = sqrt(r2);

            weights[CumNumEdge + counter] = r;
            ind1s[CumNumEdge + counter] = ind1;
            ind2s[CumNumEdge + counter] = ind2;
            counter++;

#ifdef PARALLEL
            if (counter == ThisBlock) {
                loop = 0;
                break;
            }
        }
        if (loop == 0) { break; }
#else
        }
#endif
    }

#ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, ind1s, Nedge, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ind2s, Nedge, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, weights, Nedge, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    counter = 0;
    for (int i = 0; i < Nleafs; i++) {
        for (int j = i + 1; j < Nleafs; j++) {
            W[counter].weight = weights[counter];
            W[counter].id1 = ind1s[counter];
            W[counter].id2 = ind2s[counter];
            W[counter].index1 = i;
            W[counter].index2 = j;
            counter++;
        }
    }

    free(weights);
    free(ind1s);
    free(ind2s);

    istart = 0;
    int istop = Nedge;
#ifdef PARALLEL
    int block = (int)floor(Nedge / Ntask);
    istart = ThisTask * block;
    istop = (ThisTask + 1) * block;
    if (ThisTask == Ntask - 1) istop = Nedge;
#endif

    int *flag = calloc(Nedge, sizeof(int));

    int Noff = 0;
    double w[3];
    int ind[3];
    for (int e1 = istart; e1 < istop; e1++) {

        int v1 = W[e1].index1;
        int v2 = W[e1].index2;
        ind[0] = e1;
        w[0] = W[e1].weight;

        if (w[0] < 0) continue;

        for (int v3 = 0; v3 < Nleafs; v3++) {

            if (w[0] < 0) continue;

            if (v3 == v1 || v3 == v2) continue;

            ind[1] = get_edge(Nleafs, v1, v3);
            ind[2] = get_edge(Nleafs, v2, v3);

            w[1] = W[ind[1]].weight;
            w[2] = W[ind[2]].weight;
            if (w[1] < 0 || w[2] < 0) continue;
            int elong = get_longest(w, ind);

            flag[elong] = 1;
            W[elong].weight *= -1;
        }
    }

#ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, flag, Nedge, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    for (int e1 = 0; e1 < Nedge; e1++) {
        if (flag[e1] == 0) Noff++;
    }

    struct WeightIndex2 *Wfilter =
        calloc(Nedge - Noff, sizeof(struct WeightIndex2));

    counter = 0;
    for (int e1 = 0; e1 < Nedge; e1++) {
        if (flag[e1] == 0) {
            Wfilter[counter].weight = W[e1].weight;
            Wfilter[counter].id1 = W[e1].id1;
            Wfilter[counter].id2 = W[e1].id2;
            Wfilter[counter].index1 = W[e1].index1;
            Wfilter[counter].index2 = W[e1].index2;
            counter++;
        }
    }

    Nedge -= Noff;
    free(W);
    free(flag);

    qsort(Wfilter, Nedge, sizeof(struct WeightIndex2), cmp_weight_index2);

    InMst[0] = 1;

    int remaining = Nleafs - 1;
    istart = 0;
    istop = Nedge;
#ifdef PARALLEL
    block = (int)floor(Nedge / Ntask);
    istart = ThisTask * block;
    istop = (ThisTask + 1) * block;
    if (ThisTask == Ntask - 1) istop = Nedge;
#endif

    do {
        struct WeightIndex wi;
        wi.weight = 1e10;
        wi.index = -1;

        for (int i = istart; i < istop; i++) {
            int i1 = Wfilter[i].index1;
            int i2 = Wfilter[i].index2;
            if (InMst[i1] != InMst[i2]) {
                wi.index = i;
                wi.weight = Wfilter[i].weight;
                break;
            }
        }

#ifdef PARALLEL
        MPI_Allreduce(MPI_IN_PLACE, &wi, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                      MPI_COMM_WORLD);
#endif

        InMst[Wfilter[wi.index].index1] = 1;
        InMst[Wfilter[wi.index].index2] = 1;

        R->EdgeInMST[*NumEdge].vertex1 = Wfilter[wi.index].id1;
        R->EdgeInMST[*NumEdge].vertex2 = Wfilter[wi.index].id2;
        R->EdgeInMST[*NumEdge].weight = Wfilter[wi.index].weight;
        (*NumEdge)++;

        remaining--;
    } while (remaining > 0);

    free(Wfilter);
    free(InMst);
}

void SolveSubMST(struct RegularizedRegion *R, int Npart, int *Particles,
                 int label, int *NumEdge) {

    if (Npart == 1) return;

    int Nedge = Npart * (Npart - 1) / 2;
    struct WeightIndex2 *W = calloc(Nedge, sizeof(struct WeightIndex2));
    int *InMst = calloc(Npart, sizeof(int));

    int counter = 0;
    for (int i = 0; i < Npart; i++) {
        int ind1 = Particles[i];
        for (int j = i + 1; j < Npart; j++) {
            int ind2 = Particles[j];
            double r2 = 0;
            for (int k = 0; k < 3; k++) {
                double ds = R->Pos[3 * ind1 + k] - R->Pos[3 * ind2 + k];
                r2 += ds * ds;
            }
            double r = sqrt(r2);
            W[counter].weight = r;
            W[counter].index1 = i;
            W[counter].index2 = j;
            W[counter].id1 = ind1;
            W[counter].id2 = ind2;
            counter++;
        }
    }

    qsort(W, Nedge, sizeof(struct WeightIndex2), cmp_weight_index2);

    InMst[0] = 1;

    int remaining = Npart - 1;

    do {

        for (int i = 0; i < Nedge; i++) {
            int i1 = W[i].index1;
            int i2 = W[i].index2;
            if (InMst[i1] != InMst[i2]) {
                if (InMst[i1] == 1) {
                    InMst[i2] = 1;
                } else {
                    InMst[i1] = 1;
                }

                R->EdgeInMST[*NumEdge].vertex1 = W[i].id1;
                R->EdgeInMST[*NumEdge].vertex2 = W[i].id2;
                R->EdgeInMST[*NumEdge].weight = W[i].weight;
                (*NumEdge)++;
                break;
            }
        }

        remaining--;
    } while (remaining > 0);

    free(W);
    free(InMst);
}

void SwapEdge(struct RegularizedRegion *R, int i, int j) {

    int v1 = R->EdgeInMST[i].vertex1;
    int v2 = R->EdgeInMST[i].vertex2;
    double w = R->EdgeInMST[i].weight;

    R->EdgeInMST[i].vertex1 = R->EdgeInMST[j].vertex1;
    R->EdgeInMST[i].vertex2 = R->EdgeInMST[j].vertex2;
    R->EdgeInMST[i].weight = R->EdgeInMST[j].weight;

    R->EdgeInMST[j].vertex1 = v1;
    R->EdgeInMST[j].vertex2 = v2;
    R->EdgeInMST[j].weight = w;
}

void CombineAllMST(struct RegularizedRegion *R, int NumEdge) {

    int next = -1;
    double rmin = DBL_MAX;
    for (int i = 0; i < R->NumVertex; i++) {
        double r2 = 0;
        for (int k = 0; k < 3; k++) {
            double ds = R->Pos[3 * i + k];
            r2 += ds * ds;
        }
        if (rmin > sqrt(r2)) {
            rmin = sqrt(r2);
            next = i;
        }
    }

    int counter = 0;
    R->Vertex[next].inMST = 1;

    int istop = NumEdge;

    do {

        struct WeightIndex wi;
        wi.weight = 1e10;
        wi.index = -1;

        for (int i = counter; i < istop; i++) {
            int v1 = R->EdgeInMST[i].vertex1;
            int v2 = R->EdgeInMST[i].vertex2;
            if (R->Vertex[v1].inMST != R->Vertex[v2].inMST) {
                if (R->EdgeInMST[i].weight < wi.weight) {
                    wi.weight = R->EdgeInMST[i].weight;
                    wi.index = i;
                }
            }
        }

        int new = R->EdgeInMST[wi.index].vertex1;
        int old = R->EdgeInMST[wi.index].vertex2;
        if (R->Vertex[new].inMST == 1) {
            old = R->EdgeInMST[wi.index].vertex1;
            new = R->EdgeInMST[wi.index].vertex2;
        }

        R->Vertex[new].inMST = 1;
        R->Vertex[new].parent = old;
        R->Vertex[new].level = R->Vertex[old].level + 1;
        R->Vertex[new].edgetoparent = counter;

        for (int k = 0; k < 3; k++) {
            R->State[3 * counter + k] =
                R->Pos[3 * new + k] - R->Pos[3 * old + k];
            R->State[3 * (R->NumVertex - 1) + 3 * counter + k] =
                R->Vel[3 * new + k] - R->Vel[3 * old + k];
        }

        SwapEdge(R, counter, wi.index);
        counter++;

    } while (counter < NumEdge);
}

void DivideAndConquerPrim(struct RegularizedRegion *R) {

    my_barrier();
    clock_t t0 = clock();

    for (int i = 0; i < R->NumVertex; i++) {
        R->Vertex[i].id = i;
        R->Vertex[i].degree = 0;
        R->Vertex[i].edgetoparent = -1;
        R->Vertex[i].parent = -1;
        R->Vertex[i].level = 0;
        R->Vertex[i].inMST = 0;
    }

    int Nleafs = 0;
    int MaxLeafSize = (int)floor(3 * sqrt(R->NumVertex));

    clock_t t10 = clock();
    spatialpartition(&R[0], MaxLeafSize);

    GetNumberOfLeafs(OctreeRootNode, &Nleafs);
    int *LeafSizes = calloc(Nleafs, sizeof(int));
    int **LeafParticles = calloc(Nleafs, sizeof(int *));
    clock_t t11 = clock();
    tprim1 += (double)(t11 - t10) / CLOCKS_PER_SEC;

    clock_t t20 = clock();
    Nleafs = 0;
    list_octree_level(OctreeRootNode, LeafSizes, LeafParticles, &Nleafs);
    clock_t t21 = clock();
    tprim2 += (double)(t21 - t20) / CLOCKS_PER_SEC;

    clock_t t30 = clock();
    int NumEdge = 0;
    for (int i = 0; i < Nleafs; i++) {
        SolveSubMST(R, LeafSizes[i], LeafParticles[i], i, &NumEdge);
    }
    clock_t t31 = clock();
    tprim3 += (double)(t31 - t30) / CLOCKS_PER_SEC;

    clock_t t40 = clock();
    SolveMetaMST(R, Nleafs, LeafSizes, LeafParticles, &NumEdge);
    clock_t t41 = clock();
    tprim4 += (double)(t41 - t40) / CLOCKS_PER_SEC;

    clock_t t50 = clock();
    CombineAllMST(R, NumEdge);
    clock_t t51 = clock();
    tprim5 += (double)(t51 - t50) / CLOCKS_PER_SEC;

    free_octree(OctreeRootNode);
    for (int i = 0; i < Nleafs; i++) {
        free(LeafParticles[i]);
    }
    free(LeafParticles);
    free(LeafSizes);

    clock_t t1 = clock();
    time_mst += (double)(t1 - t0) / CLOCKS_PER_SEC;
}

int EdgeSearchForPrim(struct RegularizedRegion *R, int NumLocalEdges) {

    clock_t edgesearch_time0 = clock();
    struct WeightIndex wi;

    wi.weight = 1e10;
    wi.index = -1;
    for (int i = 0; i < NumLocalEdges; i++) {
        if (R->Vertex[R->LocalEdge[i].vertex1].inMST !=
            R->Vertex[R->LocalEdge[i].vertex2].inMST) {
            wi.weight = R->LocalEdge[i].weight;
            wi.index = R->LocalEdge[i].id;
            break;
        }
    }

    clock_t edgesearch_time1 = clock();
    tprim3 += (double)(edgesearch_time1 - edgesearch_time0) / CLOCKS_PER_SEC;

    edgesearch_time0 = clock();
#ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &wi, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  MPI_COMM_WORLD);
#endif

    edgesearch_time0 = clock();
    tprim4 += (double)(edgesearch_time1 - edgesearch_time0) / CLOCKS_PER_SEC;

    return wi.index;
}

int PresentInArray(int p, int *Array, int N) {

    int result = -1;
    for (int i = 0; i < N; i++) {
        if (Array[i] == p) {
            result = i;
            break;
        }
    }
    return result;
}

int get_heaviest(double a[3], int b1, int b2, int b3) {

    int index = -1;
    double w = -1;
    for (int i = 0; i < 3; i++) {
        if (a[i] > w) {
            w = a[i];
            index = i;
        }
    }
    if (index == 0) return b1;
    if (index == 1) return b2;
    return b3;
}

void PrimMST(struct RegularizedRegion *R) {

    msttime0 = clock();

    int CumNumEdge = 0, istart = 0;
#ifdef PARALLEL
    int loop = 1, ThisBlock, jstart;
    loop_scheduling_N2(R->NumVertex, Ntask, &istart, &jstart, &ThisBlock,
                       &CumNumEdge, ThisTask);
#endif
    int NumLocalEdges = 0;
    int c = 0, counter = 0;

    double rmin = 1e10, r2;
    int last = -1;

    for (int i = 0; i < R->NumVertex; i++) {
        r2 = 0;
        for (int k = 0; k < 3; k++) {
            double ds = R->Pos[3 * i + k];
            r2 += ds * ds;
        }
        if (rmin > sqrt(r2)) {
            rmin = sqrt(r2);
            last = i;
        }
    }

    for (int i = istart; i < R->NumVertex; i++) {

#ifdef PARALLEL
        if (i > istart) { jstart = i + 1; }
        for (int j = jstart; j < R->NumVertex; j++) {
#else
        for (int j = i + 1; j < R->NumVertex; j++) {
#endif
            counter++;
            double w = 0;
            for (int k = 0; k < 3; k++) {
                double ds = R->Pos[3 * i + k] - R->Pos[3 * j + k];
                w += ds * ds;
            }
            w = sqrt(w);

            R->LocalEdge[NumLocalEdges].weight = w;
            R->LocalEdge[NumLocalEdges].vertex1 = i;
            R->LocalEdge[NumLocalEdges].vertex2 = j;
            R->LocalEdge[NumLocalEdges].id = (int)(CumNumEdge + c);

            NumLocalEdges++;
            c++;

#ifdef PARALLEL
            if (counter == ThisBlock) {
                loop = 0;
                break;
            }
        }
        if (loop == 0) { break; }
#else
        }
#endif
    }

    qsort(R->LocalEdge, NumLocalEdges, sizeof(struct GraphEdge), cmp);

    int remaining = 0;

    for (int i = 0; i < R->NumVertex; i++) {
        R->Vertex[i].id = i;
        R->Vertex[i].degree = 0;
        R->Vertex[i].edgetoparent = -1;
        R->Vertex[i].parent = -1;
        R->Vertex[i].level = 0;
        R->Vertex[i].inMST = 0;
    }

    R->Vertex[last].inMST = 1;

    c = 0;
    do {

        int NextEdge = EdgeSearchForPrim(R, NumLocalEdges);

        int v1 = -1, v2 = -1;
        get_v1_v2(R->NumVertex, NextEdge, &v1, &v2);

        R->EdgeInMST[c].vertex1 = v1;
        R->EdgeInMST[c].vertex2 = v2;

        R->Vertex[v1].degree++;
        R->Vertex[v2].degree++;

        if (R->Vertex[v1].inMST == 0) {

            R->Vertex[v1].parent = v2;
            R->Vertex[v1].level = R->Vertex[v2].level + 1;
            R->Vertex[v1].inMST = 1;
            R->Vertex[v1].edgetoparent = c;
            for (int k = 0; k < 3; k++) {
                R->State[3 * c + k] = R->Pos[3 * v1 + k] - R->Pos[3 * v2 + k];
                R->State[3 * (R->NumVertex - 1) + 3 * c + k] =
                    R->Vel[3 * v1 + k] - R->Vel[3 * v2 + k];
            }
        } else {

            R->Vertex[v2].parent = v1;
            R->Vertex[v2].level = R->Vertex[v1].level + 1;
            R->Vertex[v2].inMST = 1;
            R->Vertex[v2].edgetoparent = c;
            for (int k = 0; k < 3; k++) {
                R->State[3 * c + k] = R->Pos[3 * v2 + k] - R->Pos[3 * v1 + k];
                R->State[3 * (R->NumVertex - 1) + 3 * c + k] =
                    R->Vel[3 * v2 + k] - R->Vel[3 * v1 + k];
            }
        }
        c++;
        remaining++;
    } while (remaining < R->NumVertex - 1);

    msttime1 = clock();
    time_mst += (double)(msttime1 - msttime0) / CLOCKS_PER_SEC;
}

// computes the force function for logh
void compute_U(struct RegularizedRegion *R) {

    clock_t thistime0 = clock();

    R->U = 0;
    int istart = 0;
#ifdef PARALLEL
    int CumNumEdge, ThisBlock, jstart;
    loop_scheduling_N2(R->NumVertex, NumTaskPerGbsGroup, &istart, &jstart,
                       &ThisBlock, &CumNumEdge, ThisTask_in_GbsGroup);
    int loop = 1;
    int c = 0;
#endif

    const int Nd = GLOBAL_ND;
    int d;
    int path[2];
    int sign[2];

    for (int i = istart; i < R->NumVertex; i++) {

        double mi = R->Mass[i];
        int Li = R->Vertex[i].level;
#ifdef PARALLEL
        if (i > istart) { jstart = i + 1; }
        for (int j = jstart; j < R->NumVertex; j++) {
#else
        for (int j = i + 1; j < R->NumVertex; j++) {
#endif
            double mj = R->Mass[j];
            int Lj = R->Vertex[j].level;
            double r2 = 0;
            int proximity = 0;

            if (abs(Li - Lj) < Nd) {
                proximity =
                    check_relative_proximity(j, i, Nd, R, &d, path, sign);
            }

            if (proximity == 1) {
                for (int k = 0; k < 3; k++) {
                    double dr = 0;
                    for (int l = 0; l < d; l++) {
                        int index = path[l];
                        dr += sign[l] * R->State[3 * index + k];
                    }
                    r2 += dr * dr;
                }
            } else {
                for (int k = 0; k < 3; k++) {
                    double dr = R->Pos[3 * j + k] - R->Pos[3 * i + k];
                    r2 += dr * dr;
                }
            }

            const double invr = 1.0 / sqrt(r2);
            R->U += mi * mj * invr;
#ifdef PARALLEL
            c++;
            if (c == ThisBlock) {
                loop = 0;
                break;
            }
        }
        if (loop == 0) break;
#else
        }
#endif
    }

    clock_t ctime0 = clock();
#ifdef PARALLEL
    if (NumTaskPerGbsGroup > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &R->U, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_GBS_GROUP);
    }
#endif
    clock_t ctime1 = clock();
    time_comm += (double)(ctime1 - ctime0) / CLOCKS_PER_SEC;

    if (!isfinite(R->U)) {
        printf("U not finite.\n");
        exit(0);
    }

    R->U *= GCONST;

    clock_t thistime1 = clock();
    time_force += (double)(thistime1 - thistime0) / CLOCKS_PER_SEC;
}

// computes the force function and newtonian acceleration
void compute_U_and_Newtonian_Acc(struct RegularizedRegion *R) {

    clock_t thistime0 = clock();

    int istart, istop = R->NumVertex, jstop = R->NumVertex;
#ifdef PARALLEL
    int ThisBlock, CumNumEdge, jstart;
    loop_scheduling_N2(R->NumVertex, NumTaskPerGbsGroup, &istart, &jstart,
                       &ThisBlock, &CumNumEdge, ThisTask_in_GbsGroup);
    int loop = 1;
    int c = 0;
#else
    istart = 0;
#endif
    const int Nd = GLOBAL_ND;
    int d;
    int path[Nd];
    int sign[Nd];

    for (int i = 0; i < 3 * R->NumVertex; i++) {
        R->Acc[i] = 0;
    }

    R->U = 0;
    double dr[3], r2;

    double *Pos = R->Pos;
    double *Mass = R->Mass;
    struct GraphVertex *Vertex = R->Vertex;

    for (int i = istart; i < istop; i++) {
        const double mi = Mass[i];
        const int Li = Vertex[i].level;
#ifdef PARALLEL
        if (i > istart) jstart = i + 1;
        for (int j = jstart; j < jstop; j++) {
#else
        for (int j = i + 1; j < jstop; j++) {
#endif
            r2 = 0;
            const double mj = Mass[j];
            const int Lj = Vertex[j].level;

            int proximity = 0;
            if (abs(Li - Lj) <= Nd) {
                proximity =
                    check_relative_proximity_ND_2(j, i, R, &d, path, sign);
            }

            if (proximity) {
                for (int k = 0; k < 3; k++) {
                    dr[k] = 0;
                    for (int l = 0; l < d; l++) {
                        dr[k] += sign[l] * R->State[3 * path[l] + k];
                    }
                    r2 += dr[k] * dr[k];
                }
            } else {
                for (int k = 0; k < 3; k++) {
                    dr[k] = Pos[3 * j + k] - Pos[3 * i + k];
                    r2 += dr[k] * dr[k];
                }
            }

            const double invr = 1.0 / sqrt(r2);
            R->U += mi * mj * invr;

            const double invr3 = invr * invr * invr;

            for (int k = 0; k < 3; k++) {
                double drinvr3 = dr[k] * invr3;
                R->Acc[3 * i + k] += mj * drinvr3;
                ;
                R->Acc[3 * j + k] += -mi * drinvr3;
                ;
            }
#ifdef PARALLEL
            c++;
            if (c == ThisBlock) {
                loop = 0;
                break;
            }
#endif
        }
#ifdef PARALLEL
        if (loop == 0) break;
#endif
    }

    clock_t ctime0 = clock();
#ifdef PARALLEL
    if (NumTaskPerGbsGroup > 1) {
        MPI_Allreduce(MPI_IN_PLACE, R->Acc, 3 * R->NumVertex, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_GBS_GROUP);
        MPI_Allreduce(MPI_IN_PLACE, &R->U, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_GBS_GROUP);
    }
#endif
    clock_t ctime1 = clock();
    time_comm += (double)(ctime1 - ctime0) / CLOCKS_PER_SEC;

    int hi, lo;
    for (int i = 0; i < R->NumVertex - 1; i++) {

        int v1, v2;
        get_v1_v2_fast(R, i, &v1, &v2);

        if (R->Vertex[v1].level > R->Vertex[v2].level) {
            hi = v1;
            lo = v2;
        } else {
            hi = v2;
            lo = v1;
        }
        for (int k = 0; k < 3; k++) {
            R->MSTedgeAcc[3 * i + k] =
                GCONST * (R->Acc[3 * hi + k] - R->Acc[3 * lo + k]);
        }
    }

    R->U *= GCONST;

    clock_t thistime1 = clock();
    time_force += (double)(thistime1 - thistime0) / CLOCKS_PER_SEC;
}

void into_Cartesian_from_MST(struct RegularizedRegion *R) {

    clock_t timing0 = clock();

    int v1, v2;
    get_v1_v2_fast(R, 0, &v1, &v2);

    int first;
    if (R->Vertex[v1].level == 0) {
        first = v1;
    } else {
        first = v2;
    }
    for (int k = 0; k < 3; k++) {
        R->Pos[3 * first + k] = 0.0;
        R->Vel[3 * first + k] = 0.0;
    }
    for (int i = 0; i < R->NumVertex - 1; i++) {

        get_v1_v2_fast(R, i, &v1, &v2);

        int lo, hi;
        if (R->Vertex[v1].level > R->Vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * hi + k] = R->Pos[3 * lo + k] + R->State[3 * i + k];
            R->Vel[3 * hi + k] = R->Vel[3 * lo + k] +
                                 R->State[3 * (R->NumVertex - 1) + 3 * i + k];
        }
    }

    double M = 0, CoM[3] = {0}, CoV[3] = {0};

    for (int i = 0; i < R->NumVertex; i++) {
        M += R->Mass[i];
        for (int k = 0; k < 3; k++) {
            CoM[k] += R->Mass[i] * R->Pos[3 * i + k];
            CoV[k] += R->Mass[i] * R->Vel[3 * i + k];
        }
    }
    for (int k = 0; k < 3; k++) {
        CoM[k] /= M;
        CoV[k] /= M;
    }

    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * i + k] -= CoM[k];
            R->Vel[3 * i + k] -= CoV[k];
        }
    }

    clock_t timing1 = clock();
    time_coord += (double)(timing1 - timing0) / CLOCKS_PER_SEC;
}

void into_CoM_frame(struct RegularizedRegion *R) {

    clock_t timing0 = clock();

    for (int k = 0; k < 3; k++) {
        R->CoM_Pos[k] = 0;
        R->CoM_Vel[k] = 0;
    }
    double M = 0;
    for (int i = 0; i < R->NumVertex; i++) {
        M += R->Mass[i];
        for (int k = 0; k < 3; k++) {
            R->CoM_Pos[k] += R->Mass[i] * R->Pos[3 * i + k];
            R->CoM_Vel[k] += R->Mass[i] * R->Vel[3 * i + k];
        }
    }
    for (int k = 0; k < 3; k++) {
        R->CoM_Pos[k] /= M;
        R->CoM_Vel[k] /= M;
    }
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * i + k] -= R->CoM_Pos[k];
            R->Vel[3 * i + k] -= R->CoM_Vel[k];
        }
    }

    //for (int k = 0; k < 3; k++)
    //{
        //mst_CoM_R[k] = R->CoM_Pos[k];
        //mst_CoM_V[k] = R->CoM_Vel[k];
        //printf("INIT k %d MST TEST %g %g \n",k,R->CoM_Pos[k],mst_CoM_R[k]);
    //}


    clock_t timing1 = clock();
    time_coord += (double)(timing1 - timing0) / CLOCKS_PER_SEC;
}

void compute_T(struct RegularizedRegion *R) {
    R->T = 0;
    for (int i = 0; i < R->NumVertex; i++) {
        double dv, v2 = 0;
        for (int k = 0; k < 3; k++) {
            dv = R->Vel[3 * i + k];
            v2 += dv * dv;
        }
        R->T += 0.5 * R->Mass[i] * v2;
    }
}

void update_pos(struct RegularizedRegion *R) {

    clock_t timing0 = clock();

    int v1, v2;
    get_v1_v2_fast(R, 0, &v1, &v2);

    int first;
    if (R->Vertex[v1].level == 0) {
        first = v1;
    } else {
        first = v2;
    }
    for (int k = 0; k < 3; k++) {
        R->Pos[3 * first + k] = 0.0;
    }
    for (int i = 0; i < R->NumVertex - 1; i++) {

        get_v1_v2_fast(R, i, &v1, &v2);

        int lo, hi;
        if (R->Vertex[v1].level > R->Vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * hi + k] = R->Pos[3 * lo + k] + R->State[3 * i + k];
        }
    }

    double M = 0, CoM[3] = {0};
    for (int i = 0; i < R->NumVertex; i++) {
        M += R->Mass[i];
        for (int k = 0; k < 3; k++) {
            CoM[k] += R->Mass[i] * R->Pos[3 * i + k];
        }
    }
    for (int k = 0; k < 3; k++) {
        CoM[k] /= M;
    }
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * i + k] -= CoM[k];
        }
    }

    clock_t timing1 = clock();
    time_coord += (double)(timing1 - timing0) / CLOCKS_PER_SEC;
}

void update_vel(struct RegularizedRegion *R, double *Vel,
                enum velocity_type propagated_velocity) {

    clock_t timing0 = clock();

    int v1, v2;
    get_v1_v2_fast(R, 0, &v1, &v2);

    int first;
    if (R->Vertex[v1].level == 0) {
        first = v1;
    } else {
        first = v2;
    }
    for (int k = 0; k < 3; k++) {
        Vel[3 * first + k] = 0.0;
    }
    for (int i = 0; i < R->NumVertex - 1; i++) {

        get_v1_v2_fast(R, i, &v1, &v2);

        int lo, hi;
        if (R->Vertex[v1].level > R->Vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            if (propagated_velocity == PHYSICAL) {
                Vel[3 * hi + k] = Vel[3 * lo + k] +
                                  R->State[3 * (R->NumVertex - 1) + 3 * i + k];
            } else {
                Vel[3 * hi + k] = Vel[3 * lo + k] + R->AuxEdgeVel[3 * i + k];
            }
        }
    }

    double M = 0, CoV[3] = {0};
    for (int i = 0; i < R->NumVertex; i++) {
        M += R->Mass[i];
        for (int k = 0; k < 3; k++) {
            CoV[k] += R->Mass[i] * Vel[3 * i + k];
        }
    }
    for (int k = 0; k < 3; k++) {
        CoV[k] /= M;
    }
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            Vel[3 * i + k] -= CoV[k];
        }
    }

    clock_t timing1 = clock();
    time_coord += (double)(timing1 - timing0) / CLOCKS_PER_SEC;
}

void drift(double ds, struct RegularizedRegion *R) {

    compute_T(R);
    double dt = ds / (R->T + R->State[R->BIndex]);
    R->State[R->TimeIndex] += dt;

    for (int i = 0; i < R->NumVertex - 1; i++) {
        for (int k = 0; k < 3; k++) {
            R->State[3 * i + k] +=
                dt * R->State[3 * (R->NumVertex - 1) + 3 * i + k];
        }
    }
    update_pos(R);
}

#ifdef USE_PN

void set_auxiliary_variables(struct RegularizedRegion *R) {

    for (int i = 0; i < 3 * R->NumVertex; i++) {
        R->AuxVel[i] = R->Vel[i];
    }
    int offset = 3 * (R->NumVertex - 1);
    for (int i = 0; i < 3 * (R->NumVertex - 1); i++) {
        R->AuxEdgeVel[i] = R->State[offset + i];
    }
}

void physical_kick(double dt, struct RegularizedRegion *R) {

    int offset = 3 * (R->NumVertex - 1);
    update_vel(R, R->AuxVel, AUXILIARY);
    compute_Post_Newtonian_Acc(R, R->AuxVel);

    double dB = 0;
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            dB += R->Mass[i] * R->AuxVel[3 * i + k] * R->AccPN[3 * i + k];
        }
    }
    R->State[R->BIndex] -= dB * dt;

    for (int i = 0; i < R->NumVertex - 1; i++) {
        for (int k = 0; k < 3; k++) {
            R->State[offset + 3 * i + k] +=
                dt * (R->MSTedgeAcc[3 * i + k] + R->MSTedgeAcc_PN[3 * i + k]);
        }
    }
}

void auxiliary_kick(double dt, struct RegularizedRegion *R) {

    update_vel(R, R->Vel, PHYSICAL);
    compute_Post_Newtonian_Acc(R, R->Vel);

    for (int i = 0; i < R->NumVertex - 1; i++) {
        for (int k = 0; k < 3; k++) {
            R->AuxEdgeVel[3 * i + k] +=
                dt * (R->MSTedgeAcc[3 * i + k] + R->MSTedgeAcc_PN[3 * i + k]);
        }
    }
}
#endif

void kick(double ds, struct RegularizedRegion *R) {

    compute_U_and_Newtonian_Acc(R);
    double dt = ds / R->U;

#ifndef USE_PN
    for (int i = 0; i < R->NumVertex - 1; i++) {
        for (int k = 0; k < 3; k++) {
            R->State[3 * (R->NumVertex - 1) + 3 * i + k] +=
                dt * R->MSTedgeAcc[3 * i + k];
        }
    }
#else
    auxiliary_kick(dt / 2, R);
    physical_kick(dt, R);
    auxiliary_kick(dt / 2, R);
#endif
    update_vel(R, R->Vel, PHYSICAL);
}

double get_n_double(int i) { return (2.0 * (i + 1.0)); }
int get_n_int(int i) { return (2.0 * (i + 1.0)); }

double get_error(double a, double b) { return fabs(b - a); }

int gbs_group_with_min_workload(int *GbsGroupWorkLoad) {
    int min_workload = INT_MAX;
    int index = 0;
    for (int i = 0; i < NumGbsGroup; i++) {
        if (GbsGroupWorkLoad[i] < min_workload) {
            min_workload = GbsGroupWorkLoad[i];
            index = i;
        }
    }
    return index;
}

void divide_computational_load(struct RegularizedRegion *R) {

    int *GbsGroupWorkLoad = calloc(NumGbsGroup, sizeof(int));
    int counter = 0;

    for (int i = KMAX - 1; i >= 0; i--) {
        int ThisSubdivisionCost = get_n_int(i);
        int GroupIndex = gbs_group_with_min_workload(GbsGroupWorkLoad);
        GbsGroupWorkLoad[GroupIndex] += ThisSubdivisionCost;
        ComputationToDoList.ComputationalTask[i].GroupToPerformTask =
            GroupIndex;
        ComputationToDoList.ComputationalTask[i].ThisR = R;
        ComputationToDoList.ComputationalTask[i].NumberOfSubsteps =
            get_n_int(i);
        counter++;
    }

    ComputationToDoList.NumberOfComputationalTasks = counter;

    free(GbsGroupWorkLoad);
}

void compute_total_energy(struct RegularizedRegion *R) {

    compute_T(R);
    compute_U(R);
    R->H = R->T - R->U;
    R->B = -R->H;
    R->State[R->BIndex] = R->B;
}

void mst_leapfrog(struct RegularizedRegion *R, double Hstep,
                  int NumberOfSubsteps) {

#ifdef USE_PN
    set_auxiliary_variables(R);
#endif

    double h = Hstep / NumberOfSubsteps;

    drift(h / 2, R);
    for (int step = 0; step < NumberOfSubsteps - 1; step++) {
        kick(h, R);
        drift(h, R);
    }
    kick(h, R);
    drift(h / 2, R);
}

void copy_regularized_region(struct RegularizedRegion *Orig,
                             struct RegularizedRegion *Copy) {

    Copy->gbs_tolerance = Orig->gbs_tolerance;
    Copy->output_time_tolerance = Orig->output_time_tolerance;
    Copy->NumVertex = Orig->NumVertex;
    Copy->NumEdge = Orig->NumEdge;
    Copy->NumDynVariables = Orig->NumDynVariables;
    Copy->TimeIndex = Orig->TimeIndex;
    Copy->BIndex = Orig->BIndex;
    Copy->Hstep = Orig->Hstep;

    for (int j = 0; j < Orig->NumVertex; j++) {
        Copy->Vertex[j].id = Orig->Vertex[j].id;
        Copy->Vertex[j].degree = Orig->Vertex[j].degree;
        Copy->Vertex[j].type = Orig->Vertex[j].type;
        Copy->Vertex[j].parent = Orig->Vertex[j].parent;
        Copy->Vertex[j].edgetoparent = Orig->Vertex[j].edgetoparent;
        Copy->Vertex[j].level = Orig->Vertex[j].level;
        Copy->Vertex[j].inMST = Orig->Vertex[j].inMST;
    }
    for (int j = 0; j < Orig->NumVertex - 1; j++) {
        Copy->EdgeInMST[j].vertex1 = Orig->EdgeInMST[j].vertex1;
        Copy->EdgeInMST[j].vertex2 = Orig->EdgeInMST[j].vertex2;
    }

    for (int j = 0; j < Orig->NumDynVariables; j++) {
        Copy->State[j] = Orig->State[j];
    }
    for (int j = 0; j < 3 * Orig->NumVertex; j++) {
        Copy->Pos[j] = Orig->Pos[j];
        Copy->Vel[j] = Orig->Vel[j];
        Copy->Acc[j] = Orig->Acc[j];
        Copy->AuxVel[j] = Orig->AuxVel[j];
        Copy->AccPN[j] = Orig->AccPN[j];
    }
    for (int j = 0; j < 3 * (Orig->NumVertex - 1); j++) {
        Copy->MSTedgeAcc[j] = Orig->MSTedgeAcc[j];
        Copy->AuxEdgeVel[j] = Orig->AuxEdgeVel[j];
        Copy->MSTedgeAcc_PN[j] = Orig->MSTedgeAcc_PN[j];
    }
    for (int j = 0; j < Orig->NumVertex; j++) {
        Copy->Mass[j] = Orig->Mass[j];
    }
}

void extrapolate_single_variable(double *y, const int N, double *result,
                                 double *error, double *yscal) {

    for (int i = 0; i < N; i++) {
        gbsS[i][0] = y[i];
    }
    for (int j = 1; j < N; j++) {
        for (int i = 0; i < N - j; i++) {
            double nterm = pow(get_n_double(i + j) / get_n_double(i), 2);
            gbsS[i][j] =
                (nterm * gbsS[i + 1][j - 1] - gbsS[i][j - 1]) / (nterm - 1.0);
        }
    }
    // for (int i = 1; i < N; i++) {
    //    result[i - 1] = gbsS[0][i];
    //    error[i - 1] = fabs(gbsS[0][i] - gbsS[0][i - 1]) / yscal[i];
    //}
    result[N - 2] = gbsS[0][N - 1];
    error[N - 2] = fabs(gbsS[0][N - 1] - gbsS[0][N - 1 - 1]) / yscal[N - 1];
}

int extrapolate_all_variables(int korder) {

    int status = -1;

    struct RegularizedRegion *ThisR = NULL;
    int *SubStepDivision = calloc(korder, sizeof(int));

    int counter = 0;
    int NumDynVariables = 0;
    for (int i = 0; i < korder; i++) {
        ThisR = ComputationToDoList.ComputationalTask[i].ThisR;
        NumDynVariables = ThisR->NumDynVariables;
        SubStepDivision[counter] = i;
        counter++;
    }

    int block = (int)floor(NumDynVariables / Ntask);
    int loop_start = block * ThisTask;
    int loop_end = block * (ThisTask + 1);
    if (ThisTask == Ntask - 1) { loop_end = NumDynVariables; }

    for (int v = loop_start; v < loop_end; v++) {
        for (int i = 0; i < korder; i++) {
            y[i] = ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]]
                       .State[v];
            yscal[i] = fabs(y[i]);
        }
        extrapolate_single_variable(y, korder, result[v], error[v], yscal);
    }

    for (int i = 0; i < NumDynVariables; i++) {
        errorbuf[i] = 0;
        resultbuf[i] = 0;
    }

    for (int i = loop_start; i < loop_end; i++) {
        errorbuf[i] = error[i][korder - 2];
        resultbuf[i] = result[i][korder - 2];
    }

#ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, errorbuf, NumDynVariables, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, resultbuf, NumDynVariables, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    double gbs_error = -1;
    double gbs_tolerance =
        ComputationToDoList.ComputationalTask[SubStepDivision[0]]
            .ThisR->gbs_tolerance;

    int max_err_i = -1;
    for (int j = 0; j < NumDynVariables; j++) {
        if (errorbuf[j] >= gbs_error) {
            gbs_error = errorbuf[j];
            max_err_i = j;
        }
    }
    double Hnextfac =
        0.94 * pow(0.65 / (gbs_error / gbs_tolerance), 1.0 / (2 * korder + 1));

    if (gbs_error <= gbs_tolerance) {

        status = 1;

        for (int v = 0; v < NumDynVariables; v++) {
            ThisR->State[v] = resultbuf[v];
        }
        if (Hnextfac > 5.0) { Hnextfac = 5.0; }

        ThisR->Hstep *= Hnextfac;
        into_Cartesian_from_MST(ThisR);

        ok_steps++;

    } else {
        status = 0;
        failed_steps++;

        if (Hnextfac > 0.7) { Hnextfac = 0.7; }
        ThisR->Hstep *= Hnextfac;
    }
    free(SubStepDivision);

    return status;
}

// Integrate the subsystems using AR-MST for the desired time interval
void run_integrator(struct RegularizedRegion *R, double time_interval, double *end_time, int *stopping_condition_occurred) {


    failed_steps = 0;
    ok_steps = 0;
    time_force = 0;
    time_comm = 0;

    time_mst = 0;
    time_gbs = 0;
    time_gbscomm = 0;
    time_coord = 0;
    tprim1 = 0, tprim2 = 0, tprim3 = 0, tprim4 = 0, tprim5 = 0;


    // Needed variables
    double dt = 1e-7;
    double Hstep;
    double time;
    // Initialize the system to integrate

    into_CoM_frame(R);
    
    
#ifdef PRIM
    PrimMST(R);
#else
    DivideAndConquerPrim(R);
#endif


    compute_total_energy(R);
    R->Hstep = dt * R->U;
    R->State[R->TimeIndex] = 0;

    double percent = time_interval / 50;
    int epoch = 0;

    // Now proceed to integration
    int not_finished = 1;
    int redo = 1;
    int status;

    int i_sc;

    int i_print=0;
    double N_print=100.0;
    double f_time;
    
    do {

        do {

            redo = 0;

            // Divide the systems and substep divisions for tasks

            divide_computational_load(R);

            // Perform the leapfrog tasks
            for (int ctask = 0;
                 ctask < ComputationToDoList.NumberOfComputationalTasks;
                 ctask++) {
                if (ThisGbsGroup == ComputationToDoList.ComputationalTask[ctask]
                                        .GroupToPerformTask) {

                    // copy the system in order not to overwrite the original
                    // one
                    struct RegularizedRegion *ThisR =
                        ComputationToDoList.ComputationalTask[ctask].ThisR;

                    copy_regularized_region(
                        ThisR, &ComputationToDoList.CopyOfSingleRegion[ctask]);

                    // leapfrog the substep division
                    Hstep = ComputationToDoList.CopyOfSingleRegion[ctask].Hstep;
                    mst_leapfrog(&ComputationToDoList.CopyOfSingleRegion[ctask],
                                 Hstep,
                                 ComputationToDoList.ComputationalTask[ctask]
                                     .NumberOfSubsteps);

                    for (int i=0; i<ThisR->NumVertex; i++)
                    {
                        ThisR->Acc[i] = ComputationToDoList.CopyOfSingleRegion[ctask].Acc[i];
                    }

                }
            }

        } while (redo);

        my_barrier();
        clock_t aika0 = clock(), aika1;

#ifdef PARALLEL
        // We need to communicate the results
        for (int ctask = 0;
             ctask < ComputationToDoList.NumberOfComputationalTasks; ctask++) {
            int root = ComputationToDoList.ComputationalTask[ctask]
                           .GroupToPerformTask *
                       NumTaskPerGbsGroup;
            double *commbuf =
                ComputationToDoList.CopyOfSingleRegion[ctask].State;

            int count =
                ComputationToDoList.CopyOfSingleRegion[ctask].NumDynVariables;
            if (ThisTask != root) {
                for (int i = 0; i < count; i++) {
                    commbuf[i] = 0;
                }
            }
            MPI_Bcast(commbuf, count, MPI_DOUBLE, root, MPI_COMM_WORLD);
        }
        aika1 = clock();

        time_gbscomm += (double)(aika1 - aika0) / CLOCKS_PER_SEC;
#endif

        // Now preform the GBS extrapolation, check convergence and obtain next
        // timestep

        my_barrier();
        aika0 = clock();
        status = extrapolate_all_variables(KMAX);
        aika1 = clock();
        time_gbs += (double)(aika1 - aika0) / CLOCKS_PER_SEC;

        if (status == 1) {

            // Check whether we are ready
            time = R[0].State[R[0].TimeIndex];
            //printf("t %g\n",time);
            R[0].time = time;
            if (fabs(time - time_interval) / time_interval <
                R->output_time_tolerance) {
                not_finished = 0;
            } else {

                compute_total_energy(R);
                if (time > time_interval) {
                    R->Hstep = -R->U * (time - time_interval);
                } else if (R->Hstep / R->U + time > time_interval) {
                    R->Hstep = R->U * (time_interval - time);
                }

                /* Some debug printing */
                if (time < 0.0)
                {
                    printf("MSTAR -- WARNING time = %g < 0!\n",time);
                }
                f_time = time/time_interval;
                if ( ((int) (f_time*N_print) ) == i_print)
                {
                    printf("MSTAR -- t %g completed %.1f %%\n",time,f_time*100.0);
                    i_print+=1;
                }
                //printf(">>>> %d %g %d\n",i_print,f_time,(int) (f_time*N_print));

                /* Stopping conditions */
                int possible_stopping_condition;
                double Delta_t_stopping_condition;
                stopping_condition_function(R, &possible_stopping_condition, stopping_condition_occurred, &Delta_t_stopping_condition);

                if (*stopping_condition_occurred == 1)
                {
                    not_finished = 0;
                    printf("Stopping condition at t=%g \n",time);
                }
                if (possible_stopping_condition == 1)
                {
                    i_sc++;
                    
                    double sc_Hstep = R->U * Delta_t_stopping_condition;
                    if ( fabs(sc_Hstep) <= fabs(R->Hstep) ) // Make sure stopping condition detection does not interfere with the timestep determined above
                    {
                        R->Hstep = sc_Hstep;
                    }
                    
                    if (i_sc > 100) /* The stopping condition was probably not real; give up */
                    {
                        possible_stopping_condition = 0;
                        R->Hstep = R->U * dt;
                        i_sc = 0;
                    }
                    //printf("possible collision; i_col %d; new Hstep %g\n",i_col,R->Hstep);
                }
                else
                {
                    if (i_sc>0) /* Reset the timestep if stopping condition iterations occurred; is necessary since otherwise the code could stall indefinitely because of persistently negative timesteps */
                    {
                        R->Hstep = R->U * dt;
                        i_sc = 0;
                    }
                }
                
                int old_epoch = epoch;
                epoch = (int)floor(time / percent);
                if (epoch != old_epoch) {
#ifdef PRIM
                    PrimMST(&R[0]);
#else
                    DivideAndConquerPrim(&R[0]);
#endif
                }
            }
        }

    } while (not_finished);
    
    
    out_of_CoM_frame(R);
    
    *end_time = time;
}

void allocate_regularized_region(struct RegularizedRegion *S, int N) {

    S->gbs_tolerance = GBSTOL;

    S->output_time_tolerance = 1e-4;
    S->stopping_condition_tolerance = 1e-6;
    S->NumVertex = N;

    S->NumEdge = (N * (N - 1)) / 2;
    S->NumDynVariables = 6 * (N - 1) + 2;

    S->TimeIndex = S->NumDynVariables - 1;
    S->BIndex = S->NumDynVariables - 2;

    S->Pos = calloc(3 * N, sizeof(double));
    S->Vel = calloc(3 * N, sizeof(double));
    S->State = calloc(S->NumDynVariables, sizeof(double));
    S->Mass = calloc(N, sizeof(double));
    S->Acc = calloc(3 * N, sizeof(double));
    S->MSTedgeAcc = calloc(3 * (N - 1), sizeof(double));

    S->AccPN = calloc(3 * N, sizeof(double));
    S->AuxVel = calloc(3 * N, sizeof(double));
    S->AuxEdgeVel = calloc(3 * (N - 1), sizeof(double));
    S->MSTedgeAcc_PN = calloc(3 * (N - 1), sizeof(double));

    S->Vertex = calloc(N, sizeof(struct GraphVertex));

    double safetybuf = 1.05;
    int NumLocalEdge = (int)floor(safetybuf * S->NumEdge / Ntask);
    S->LocalEdge = calloc(NumLocalEdge, sizeof(struct GraphEdge));
    S->LocalEdgeSubset = calloc(NumLocalEdge, sizeof(struct GraphEdge));
    S->EdgeInMST = calloc(N - 1, sizeof(struct GraphEdge));
    
    S->Radius = calloc(N, sizeof(double));
    S->Stopping_Condition_Mode = calloc(N, sizeof(int));
    S->Stopping_Condition_Partner = calloc(N, sizeof(int));
    S->Stopping_Condition_Roche_Lobe_Radius = calloc(N, sizeof(double));

    S->Index = calloc(N, sizeof(int));
}

void allocate_armst_structs(struct RegularizedRegion **R, int MaxNumPart) {

    int MaxDynamicalVariables = 6 * (MaxNumPart - 1) + 2;
    *R = calloc(1, sizeof(struct RegularizedRegion));

    allocate_regularized_region(R[0], MaxNumPart);

    gbsS = calloc(KMAX, sizeof(double *));
    for (int i = 0; i < KMAX; i++) {
        gbsS[i] = calloc(KMAX, sizeof(double));
    }

    int korder = 8;
    int NumDynVariables = MaxDynamicalVariables;

    errorbuf = calloc(NumDynVariables, sizeof(double));
    resultbuf = calloc(NumDynVariables, sizeof(double));
    y = calloc(korder, sizeof(double));
    yscal = calloc(korder, sizeof(double));
    maxerror = calloc(korder - 1, sizeof(double));
    error = calloc(NumDynVariables, sizeof(double *));
    result = calloc(NumDynVariables, sizeof(double *));

    for (int i = 0; i < NumDynVariables; i++) {
        error[i] = calloc(korder - 1, sizeof(double));
        result[i] = calloc(korder - 1, sizeof(double));
    }
    indexlist = calloc(8, sizeof(int *));

    for (int i = 0; i < 8; i++) {
        indexlist[i] = calloc(MaxNumPart, sizeof(int));
    }

    ParticleInNode = calloc(MaxNumPart, sizeof(struct OctreeNode *));
    int NumAlloc = KMAX;

    ComputationToDoList.ComputationalTask =
        calloc(NumAlloc, sizeof(struct ComTask));
    for (int i = 0; i < NumAlloc; i++) {
        ComputationToDoList.ComputationalTask[i].ResultState =
            calloc(MaxDynamicalVariables, sizeof(double));
    }
    ComputationToDoList.CopyOfSingleRegion =
        calloc(NumAlloc, sizeof(struct RegularizedRegion));
    for (int i = 0; i < NumAlloc; i++) {
        allocate_regularized_region(&ComputationToDoList.CopyOfSingleRegion[i],
                                    MaxNumPart);
    }
    
   
}

// free
void free_data(struct RegularizedRegion *R) {

    for (int i = 0; i < KMAX; i++) {
        free(ComputationToDoList.ComputationalTask[i].ResultState);
    }
    free(ComputationToDoList.ComputationalTask);

    free(R->LocalEdge);
    free(R->Vertex);
    free(R->Pos);
    free(R->Vel);
    free(R->Acc);
    free(R->Mass);
    free(R->State);
    free(R->MSTedgeAcc);
    free(R->EdgeInMST);
    free(R);

    free(gbsS);
}

// Read in the test IC
void read_ic_from_file(struct RegularizedRegion *R, char *INPUTFILE) {
    FILE *fp;
    for (int task = 0; task < Ntask; task++) {
        my_barrier();
        if (ThisTask == task) {
            fp = fopen(INPUTFILE, "r");
            for (int i = 0; i < R->NumVertex; i++) {
                int status = fscanf(fp, "%d %lf %lf %lf %lf %lf %lf %lf\n",
                                    &R->Vertex[i].type, &R->Mass[i],
                                    &R->Pos[3 * i + 0], &R->Pos[3 * i + 1],
                                    &R->Pos[3 * i + 2], &R->Vel[3 * i + 0],
                                    &R->Vel[3 * i + 1], &R->Vel[3 * i + 2]);
                if (!(status > 0)) { die(); }
            }
            fclose(fp);
        }
        my_barrier();
    }
    my_barrier();
}

// The test main function
#ifdef IGNORE
int main(int argc, char *argv[]) {

    initialize_mpi_or_serial();

    my_barrier();

    GBSTOL = atof(argv[3]);
    MAXPART = atoi(argv[2]);

    failed_steps = 0;
    ok_steps = 0;
    time_force = 0;
    time_comm = 0;

    time_mst = 0;
    time_gbs = 0;
    time_gbscomm = 0;
    time_coord = 0;
    tprim1 = 0, tprim2 = 0, tprim3 = 0, tprim4 = 0, tprim5 = 0;

    my_barrier();

    /////////////////////////////////////////

    struct RegularizedRegion *R;
    int MaxNumPart = MAXPART;
    allocate_armst_structs(&R, MaxNumPart);
    read_ic_from_file(R, argv[1]);
    double IntegrationTime = 1e-5;

    // Integrate for some time interval

    wtime0 = clock();
    run_integrator(R, IntegrationTime);
    wtime1 = clock();

    // Do some output

    double elapsed = (double)(wtime1 - wtime0) / CLOCKS_PER_SEC;

    if (ThisTask == 0) {
        printf("%d %d %d %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Ntask,
               MAXPART, ok_steps + failed_steps, time_force - time_comm,
               time_comm, elapsed, time_mst, time_gbs, time_gbscomm,
               time_coord);

        printf("Total time: %.10e\n", elapsed);
#if 1
        printf("force:     %lf\n", 100 * (time_force - time_comm) / elapsed);
        printf("forcecomm: %lf\n", 100 * (time_comm) / elapsed);
        printf("prim:      %lf\n", 100 * (time_mst) / elapsed);
        printf("gbs:       %lf\n", 100 * (time_gbs) / elapsed);
        printf("divcomm:   %lf\n", 100 * (time_gbscomm) / elapsed);
        printf("coord  :   %lf\n", 100 * (time_coord) / elapsed);

        printf("\n--- prim ---\n");
#ifdef PRIM
        printf("%lf %lf %lf %lf %lf\n", 100 * tprim1 / elapsed,
               100 * tprim2 / elapsed, 100 * (tprim3) / elapsed,
               100 * tprim4 / elapsed, 100 * tprim5 / elapsed);
#else
        printf("octree           %lf\n", 100 * tprim1 / elapsed);
        printf("octree particles %lf\n", 100 * tprim2 / elapsed);
        printf("sub-msts         %lf\n", 100 * tprim3 / elapsed);
        printf("meta-mst         %lf\n", 100 * tprim4 / elapsed);
        printf("combine msts     %lf\n", 100 * tprim5 / elapsed);
#endif
        fflush(stdout);
#endif
        printf("Physical time %e\n", R->time);
        int n = MaxNumPart > 10 ? 10 : MaxNumPart;
        for (int i = 0; i < n; ++i) {
            printf(
                "Particle %d\n"
                " pos = [%e, %e, %e]\n"
                " vel = [%e, %e, %e]\n",
                i, R->Pos[3 * i + 0], R->Pos[3 * i + 1], R->Pos[3 * i + 2],
                R->Vel[3 * i + 0], R->Vel[3 * i + 1], R->Vel[3 * i + 2]);
        }
    }

    // Free stuff and finalize
    free_data(R);
#ifdef PARALLEL
    MPI_Finalize();
#endif
    return 0;
}
#endif

void stopping_condition_function(struct RegularizedRegion *R, int *possible_stopping_condition, int *stopping_condition_occurred, double *Delta_t_min)
{

    int istart, istop = R->NumVertex, jstop = R->NumVertex;
//#ifdef PARALLEL
//    int ThisBlock, CumNumEdge, jstart;
//    loop_scheduling_N2(R->NumVertex, NumTaskPerGbsGroup, &istart, &jstart,
//                       &ThisBlock, &CumNumEdge, ThisTask_in_GbsGroup);
//    int loop = 1;
//    int c = 0;
//#else
    istart = 0;
//#endif
    const int Nd = GLOBAL_ND;
    int d;
    int path[Nd];
    int sign[Nd];

    double dr[3], r, r2;
    double dv[3], da[3];
    double r_crit,r_crit_p2;
    double r_criti, r_critj;
    double rdotv,vdotv;
    double determinant;
    double Delta_t;
    
    *Delta_t_min = 1e100;
    *possible_stopping_condition = 0;
    *stopping_condition_occurred = 0;
    
    double *Pos = R->Pos;
    double *Vel = R->Vel;
    double *Acc = R->Acc;
    double *Mass = R->Mass;
    double *Radius = R->Radius;
    int *Stopping_Condition_Mode = R->Stopping_Condition_Mode;

    struct GraphVertex *Vertex = R->Vertex;

    for (int i = istart; i < istop; i++) 
    {
        R->Stopping_Condition_Partner[i] = -1; /* default: no stopping condition partner */
        R->Stopping_Condition_Roche_Lobe_Radius[i] = -1;
    }

    for (int i = istart; i < istop; i++) 
    {
        const double mi = Mass[i];
        const double radiusi = Radius[i];
        const int Li = Vertex[i].level;
        const int modei = Stopping_Condition_Mode[i];
//#ifdef PARALLEL
//        if (i > istart)
//        {
//            jstart = i + 1;
//        }
//        for (int j = jstart; j < jstop; j++)
//        {
//#else
        for (int j = i + 1; j < jstop; j++)
        {
//#endif
            r2 = 0;
            const double mj = Mass[j];
            const double radiusj = Radius[j];
            const int Lj = Vertex[j].level;
            const int modej = Stopping_Condition_Mode[j];

            int proximity = 0;
            if (abs(Li - Lj) <= Nd) {
                proximity =
                    check_relative_proximity_ND_2(j, i, R, &d, path, sign);
            }

            if (proximity) {
                for (int k = 0; k < 3; k++) {
                    dr[k] = 0;
                    for (int l = 0; l < d; l++) {
                        dr[k] += sign[l] * R->State[3 * path[l] + k];
                    }
                    r2 += dr[k] * dr[k];
                }
            } else {
                for (int k = 0; k < 3; k++) {
                    dr[k] = Pos[3 * j + k] - Pos[3 * i + k];
                    r2 += dr[k] * dr[k];
                }
            }
            
            for (int k=0; k<3; k++)
            {
                dv[k] = Vel[3 * j + k] - Vel[3 * i + k];
                da[k] = Acc[3 * j + k] - Acc[3 * i + k];
            }
            
            if (modei == 0 && modej == 0) // classical collision detection
            {
                r_crit = radiusi + radiusj;
            }
            else if (modei == 1 && modej == 1) // RLOF
            {
                r_criti = radiusi / fq_RLOF_Eggleton(mi,mj);
                r_critj = radiusj / fq_RLOF_Eggleton(mj,mi);
                r_crit = CV_max( r_criti, r_critj );
            }
            else
            {
                printf("Mixed Stopping_Condition_Mode not (yet) allowed!\n");
                exit(-1);
            }
            
            r_crit_p2 = r_crit * r_crit;
            rdotv = 0.0;
            vdotv = 0.0;
            for (int k=0; k<3; k++)
            {
                rdotv += dr[k] * dv[k];
                vdotv += dv[k] * dv[k];
            }
            r = sqrt(r2);
            
            /* Determine if a stopping condition actually occurred within the tolerance */
            if ( fabs(r - r_crit)/r < R->stopping_condition_tolerance )
            {
                *stopping_condition_occurred = 1;
                R->Stopping_Condition_Partner[i] = j;
                R->Stopping_Condition_Partner[j] = i;
                
                if (modei == 1 && modej == 1)
                {
                    if (fabs(r - r_criti)/r < R->stopping_condition_tolerance)
                    {
                        R->Stopping_Condition_Roche_Lobe_Radius[i] = r_criti * fq_RLOF_Eggleton(mi,mj);
                    }
                    else if (fabs(r - r_critj)/r < R->stopping_condition_tolerance)
                    {
                        R->Stopping_Condition_Roche_Lobe_Radius[j] = r_critj * fq_RLOF_Eggleton(mj,mi);
                    }
                }
            }

            /* Estimate if a stopping condition could have occurred in the past or future,
             * and determine the timestep needed to bring the integrator there. */
            int interpolation_method = 1;
            if (interpolation_method == 0) // straight line trajectories
            {
                determinant = rdotv * rdotv - (r2 - r_crit_p2) * vdotv;
                if (determinant >= 0.0)
                {
                    *possible_stopping_condition = 1;
                    Delta_t = ( - rdotv - sqrt(determinant) ) / vdotv;
                    if (fabs(Delta_t) < fabs(*Delta_t_min)) // determine smallest dt in case of multiple stopping conditions (unlikely)
                    {
                         *Delta_t_min = Delta_t;
                    }
                }
            }
            else if (interpolation_method == 1) // parabolic trajectories
            {
                double rdota = 0.0;
                double vdota = 0.0;
                double adota = 0.0;
                for (int k=0; k<3; k++)
                {
                    rdota += dr[k] * da[k];
                    vdota += dv[k] * da[k];
                    adota += da[k] * da[k];
                }
                double a = 0.25 * adota;
                double b = vdota;
                double c = vdotv + rdota;
                double d = 2.0 * rdotv;
                double e = r2 - r_crit_p2;
                
                double p = (8.0*a*c - 3.0*b*b)/(8.0*a*a);
                double q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d)/(8.0*a*a*a);
                
                double Delta_0 = c*c - 3.0*b*d + 12.0*a*e;
                double Delta_1 = 2.0*c*c*c - 9.0*b*c*d + 27.0*b*b*e + 27.0*a*d*d - 72.0*a*c*e;
                double Delta = -(1.0/27.0)*( Delta_1*Delta_1 - 4.0*Delta_0*Delta_0*Delta_0 );
                double Q = pow( 0.5*( Delta_1 + sqrt(Delta_1*Delta_1 - 4.0*Delta_0*Delta_0*Delta_0) ), 1.0/3.0);
                double S = 0.5 * sqrt( (-2.0/3.0)*p + (1.0/(3.0*a))*(Q + Delta_0/Q) );
//                printf("a %g b %g c %g d %g e %g\n",a,b,c,d,e);
//                printf("p %g q %g Delta_0 %g Delta_1 %g Q %g S %g Delta %g\n",p,q,Delta_0,Delta_1,Q,S,Delta);

                /* Solutions to quartic equation for Delta t */
                double x1 = -b/(4.0*a) - S - 0.5*sqrt(-4.0*S*S - 2.0*p + q/S);
                double x2 = -b/(4.0*a) - S + 0.5*sqrt(-4.0*S*S - 2.0*p + q/S);
                double x3 = -b/(4.0*a) + S - 0.5*sqrt(-4.0*S*S - 2.0*p - q/S);
                double x4 = -b/(4.0*a) + S + 0.5*sqrt(-4.0*S*S - 2.0*p - q/S);

                 /* 0-4 solutions can be real; check for realness and adopt the smallest real value if exists */
                double xmin = 1e100;
                if (x1<xmin)
                {
                    xmin = x1;
                }
                if (x2<xmin)
                {
                    xmin = x2;
                }
                if (x3<xmin)
                {
                    xmin = x3;
                }
                if (x4<xmin)
                {
                    xmin = x4;
                }
                //printf("x1 %g x2 %g x3 %g x4 %g xmin %g\n",x1,x2,x3,x4,xmin);

                if (xmin != 1e100)
                {
                    *possible_stopping_condition = 1;
                    Delta_t = xmin;
                    if (fabs(Delta_t) < fabs(*Delta_t_min)) // determine smallest dt in case of multiple stopping conditions (unlikely)
                    {
                         *Delta_t_min = Delta_t;
                    }
                }


            }
            //printf("col test r %g r_crit %g rdotv %g vdotv %g collision_occurred %d possible_collision %d Delta_t %g Delta_t_min %g stopping_condition_tolerance %g (r-r_crit)/r %g \n",r,r_crit,rdotv,vdotv,*collision_occurred,*possible_collision,Delta_t,*Delta_t_min,R->stopping_condition_tolerance,fabs(r - r_crit)/r);
        }
    }
        
//#ifdef PARALLEL
//            c++;
//            if (c == ThisBlock) {
//                loop = 0;
//                break;
//            }
//#endif
//        }
//#ifdef PARALLEL
        //if (loop == 0) break;
//#endif
//    }

    //clock_t ctime0 = clock();
//#ifdef PARALLEL
//    if (NumTaskPerGbsGroup > 1) {
//        MPI_Allreduce(MPI_IN_PLACE, R->Acc, 3 * R->NumVertex, MPI_DOUBLE,
//                      MPI_SUM, MPI_COMM_GBS_GROUP);
//        MPI_Allreduce(MPI_IN_PLACE, &R->U, 1, MPI_DOUBLE, MPI_SUM,
//                      MPI_COMM_GBS_GROUP);
//    }
//#endif

    
    //clock_t ctime1 = clock();
    //time_comm += (double)(ctime1 - ctime0) / CLOCKS_PER_SEC;

}

double fq_RLOF_Eggleton(double m1, double m2)
{
    /* m1 is the donor; m2 the accretor */
    double q = m1/m2;
    double q_p1div3 = pow(q,1.0/3.0);
    double q_p2div3 = q_p1div3*q_p1div3;
    return 0.49*q_p2div3/( 0.6*q_p2div3 + log(1.0 + q_p1div3) );
}

void out_of_CoM_frame(struct RegularizedRegion *R)
{

//    for (int k = 0; k < 3; k++)
//    {
//        printf("FIN k %d MST TEST %g \n",k,R->CoM_Pos[k]);
//    }
     
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * i + k] += R->CoM_Pos[k];
            R->Vel[3 * i + k] += R->CoM_Vel[k];
        }
    }

}


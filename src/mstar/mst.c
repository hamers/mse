#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "regularization.h"

double GBSTOL;
int MAXPART;

double *errorbuf;
double *resultbuf;
double *y;
double *yscal;
double *maxerror;
double **extrapolation_error;
double **result;
double **gbsS;

struct ToDoList ComputationToDoList;

void die(void) {
    printf("The code is dead because 'die()' was called..\n");
    //exit(0);
    error_code = 27;
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

void initialize_mpi_or_serial(void) {
#if defined(USE_PN_SPIN) && !defined(USE_PN)
    printf("PN spin terms enabled.\nPlease enable USE_PN as well. Currently USE_PN not enabled.\n");
    //exit(0);
    error_code = 28;
#endif
}

double get_timestep_estimate( struct RegularizedRegion *R ){
        int np = R->NumVertex;
        double dt_min = DBL_MAX;
        for(int i=0;i<np;i++){
                for(int j=i+1;j<np;j++){
                        double mi = R->Mass[i];
                        double mj = R->Mass[j];
                        double dr[3], r2=0;
                        for(int k=0;k<3;k++){
                                dr[k] = R->Pos[3*j+k]-R->Pos[3*i+k];
                                r2 += dr[k]*dr[k];
                        }
                        double r3 = r2*sqrt(r2);
                        double dt_freefall = 0.1*sqrt( r3/(GCONST*(mi+mj)));
                        if( dt_freefall<dt_min ) dt_min = dt_freefall;
                        double dv[3], v2=0;
                        for(int k=0;k<3;k++){
                                dv[k] = R->Vel[3*j+k]-R->Vel[3*i+k];
                                v2 += dv[k]*dv[k];
                        }
                        double dt_flyby = 0.1*sqrt(r2/v2);
                        if( dt_flyby<dt_min ) dt_min = dt_flyby;
                }
        }
        return dt_min;
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
        path[i] = -1;
        sign[i] = 0;
    }

    int proximity = 0;
    if (R->Vertex[v1].level == R->Vertex[v2].level) {

        if (R->Vertex[v1].level == 0) {
            printf("This level stuff should not happen\n");
            printf("--- v1 v2 %d %d ---levels: %d %d\n", R->Vertex[v1].id, R->Vertex[v2].id, R->Vertex[v1].level,
                   R->Vertex[v2].level);
            //exit(0);
            error_code = 29;
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

int EdgeSearchForPrim(struct RegularizedRegion *R, int NumLocalEdges) {

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

    return wi.index;
}

void PrimMST(struct RegularizedRegion *R) {

    int NumLocalEdges = 0;
    int c = 0, counter = 0;

    double rmin = DBL_MAX, r2;
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

    for (int i = 0; i < R->NumVertex; i++) {

        for (int j = i + 1; j < R->NumVertex; j++) {

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
            R->LocalEdge[NumLocalEdges].id = c;

            NumLocalEdges++;
            c++;

        }
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

#ifdef USE_PN_SPIN
    for (int k = 0; k < 3; k++) {
	    R->SpinState[k]  = R->Spin_S[3*last+k];
    }
    R->SpinStateIndex[0] = last;
#endif

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
                R->State[3*c+k] = R->Pos[3 * v1 + k] - R->Pos[3 * v2 + k];
                R->State[3*(R->NumVertex - 1) + 3 * c + k] = R->Vel[3 * v1 + k] - R->Vel[3 * v2 + k];
            }
#ifdef USE_PN_SPIN
	    for (int k = 0; k < 3; k++) {
		R->SpinState[3*(c+1)+k] = R->Spin_S[3*v1+k];
	    }
	    R->SpinStateIndex[c+1] = v1;
#endif
        } else {
            R->Vertex[v2].parent = v1;
            R->Vertex[v2].level = R->Vertex[v1].level + 1;
            R->Vertex[v2].inMST = 1;
            R->Vertex[v2].edgetoparent = c;
            for (int k = 0; k < 3; k++) {
                R->State[3 * c + k] = R->Pos[3 * v2 + k] - R->Pos[3 * v1 + k];
                R->State[3 * (R->NumVertex - 1) + 3 * c + k] = R->Vel[3 * v2 + k] - R->Vel[3 * v1 + k];
            }
#ifdef USE_PN_SPIN
	    for (int k = 0; k < 3; k++) {
		R->SpinState[3*(c+1)+k] = R->Spin_S[3*v2+k];
	    }
	    R->SpinStateIndex[c+1] = v2;
#endif
    }

    c++;
    remaining++;
    } while (remaining < R->NumVertex-1);

}

// computes the force function for logh
void compute_U(struct RegularizedRegion *R) {

    R->U = 0;
    int istart = 0;

    const int Nd = GLOBAL_ND;
    int d;
    int path[2];
    int sign[2];

    for (int i = istart; i < R->NumVertex; i++) {

        double mi = R->Mass[i];
        int Li = R->Vertex[i].level;

        for (int j = i + 1; j < R->NumVertex; j++) {
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

	    

        }
    }

    if (!isfinite(R->U)) {
        printf("mst.c -- ERROR: U not finite.\n");
        //exit(0);
        error_code = 30;
        longjmp(jump_buf,1);
    }

    R->U *= GCONST;

}

// computes the force function and newtonian acceleration
void compute_U_and_Newtonian_Acc(struct RegularizedRegion *R) {

    int istart, istop = R->NumVertex, jstop = R->NumVertex;
    istart = 0;
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

        for (int j = i + 1; j < jstop; j++) {
            r2 = 0;
            const double mj = Mass[j];
            const int Lj = Vertex[j].level;

            int proximity = 0;
            if (abs(Li - Lj) <= Nd) {
                proximity =
                    //check_relative_proximity_ND_2(j, i, R, &d, path, sign);
                    check_relative_proximity(j, i, Nd, R, &d, path, sign);
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
        }
    }

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

}


void into_Cartesian_from_MST(struct RegularizedRegion *R) {

    int v1, v2, first, lo, hi;

    get_v1_v2_fast(R, 0, &v1, &v2);

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

        if (R->Vertex[v1].level > R->Vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * hi + k] = R->Pos[3 * lo + k] + R->State[3 * i + k];
            R->Vel[3 * hi + k] = R->Vel[3 * lo + k] + R->State[3 * (R->NumVertex - 1) + 3 * i + k];
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

}

void into_CoM_frame(struct RegularizedRegion *R) {

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

}

void update_vel(struct RegularizedRegion *R, double *Vel, enum velocity_type propagated_velocity) {

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
                Vel[3 * hi + k] = Vel[3 * lo + k] + R->State[3 * (R->NumVertex - 1) + 3 * i + k];
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

#ifdef USE_PN_SPIN
    from_SpinState_to_Spin_S(R);
#endif

    for (int i = 0; i < 3 * R->NumVertex; i++) {
        R->AuxVel[i]    = R->Vel[i];
#ifdef USE_PN_SPIN
	R->AuxSpin_S[i] = R->Spin_S[i];
#endif
    }
    int offset = 3 * (R->NumVertex - 1);
    for (int i = 0; i < 3 * (R->NumVertex - 1); i++) {
        R->AuxEdgeVel[i] = R->State[offset + i];
    }
}

void physical_kick(double dt, struct RegularizedRegion *R) {

    int offset = 3 * (R->NumVertex - 1);
    update_vel(R, R->AuxVel, AUXILIARY);

#ifdef USE_PN_SPIN
    compute_Post_Newtonian_Acc(R, R->AuxVel, R->AuxSpin_S );
#else
    compute_Post_Newtonian_Acc(R, R->AuxVel );
#endif

    double dB = 0;
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            dB += R->Mass[i] * R->AuxVel[3 * i + k] * R->AccPN[3 * i + k];
        }
    }
    R->State[R->BIndex] -= dB * dt;

    for (int i = 0; i < R->NumVertex - 1; i++) {
        for (int k = 0; k < 3; k++) {
            R->State[offset + 3 * i + k] += dt * (R->MSTedgeAcc[3 * i + k] + R->MSTedgeAcc_PN[3 * i + k]);
        }
    }

#ifdef USE_PN_SPIN
    for (int i = 0;i< R->NumVertex; i++) {
	for (int k = 0; k < 3; k++) {
		R->Spin_S[3*i+k] += dt*R->Spin_dS_PN[3*i+k];
	}
   }
#endif

}

void auxiliary_kick(double dt, struct RegularizedRegion *R) {

    update_vel(R, R->Vel, PHYSICAL);

#ifdef USE_PN_SPIN
    compute_Post_Newtonian_Acc(R, R->Vel, R->Spin_S );
#else
    compute_Post_Newtonian_Acc(R, R->Vel );
#endif

    for (int i = 0;i< R->NumVertex-1; i++) {
        for (int k = 0; k < 3; k++) {
            R->AuxEdgeVel[3 * i + k] += dt * (R->MSTedgeAcc[3 * i + k] + R->MSTedgeAcc_PN[3 * i + k]);
        }
    }

#ifdef USE_PN_SPIN
    for (int i = 0;i< R->NumVertex; i++) {
	for (int k = 0; k < 3; k++) {
		R->AuxSpin_S[3*i+k] += dt*R->Spin_dS_PN[3*i+k];
	}
    }
#endif

}
#endif

void kick(double ds, struct RegularizedRegion *R) {

    compute_U_and_Newtonian_Acc(R);
    double dt = ds / R->U;

#ifndef USE_PN
    for (int i = 0; i < R->NumVertex - 1; i++) {
        for (int k = 0; k < 3; k++) {
        	R->State[3 * (R->NumVertex - 1) + 3 * i + k] += dt * R->MSTedgeAcc[3 * i + k];
        }
    }
#else

    auxiliary_kick(dt / 2, R);
    physical_kick(dt, R);
    auxiliary_kick(dt / 2, R);

#ifdef USE_PN_SPIN
    from_Spin_S_to_SpinState(R);
#endif

#endif
    update_vel(R, R->Vel, PHYSICAL);

}

double get_n_double(int i) { return (2.0 * (i + 1.0)); }
int get_n_int(int i) { return (2.0 * (i + 1.0)); }

void divide_computational_load(struct RegularizedRegion *R) {
    int counter = 0;
    for (int i = KMAX - 1; i >= 0; i--) {
        ComputationToDoList.ComputationalTask[i].ThisR = R;
        ComputationToDoList.ComputationalTask[i].NumberOfSubsteps = get_n_int(i);
        counter++;
    }
    ComputationToDoList.NumberOfComputationalTasks = counter;
}

void compute_total_energy(struct RegularizedRegion *R) {

    compute_T(R);
    compute_U(R);
    R->H = R->T - R->U;
    R->B = -R->H;
    R->State[R->BIndex] = R->B;
}

void mst_leapfrog(struct RegularizedRegion *R, double Hstep, int NumberOfSubsteps) {

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

#ifdef USE_PN_SPIN
void from_SpinState_to_Spin_S( struct RegularizedRegion *R ){
	for (int i=0; i < R->NumVertex; i++){
		int index = R->SpinStateIndex[i];
		for(int k=0;k<3;k++){
			R->Spin_S[3*index+k] = R->SpinState[3*i+k];
		}
	}
}
void from_Spin_S_to_SpinState( struct RegularizedRegion *R ){

	for (int i=0; i < R->NumVertex; i++){
		int index = R->SpinStateIndex[i];
		for(int k=0;k<3;k++){
			R->SpinState[3*i+k] = R->Spin_S[3*index+k];
		}
	}
}
#endif

void copy_regularized_region(struct RegularizedRegion *Orig,
                             struct RegularizedRegion *Copy) {

    Copy->gbs_tolerance = Orig->gbs_tolerance;
    Copy->output_time_tolerance = Orig->output_time_tolerance;
    Copy->NumVertex = Orig->NumVertex;
    Copy->NumEdge = Orig->NumEdge;
    Copy->NumDynVariables = Orig->NumDynVariables;
#ifdef USE_PN_SPIN
    Copy->NumSpinVariables = Orig->NumSpinVariables;
#endif
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

#ifdef USE_PN_SPIN
    for (int j = 0; j < Orig->NumVertex; j++) {
	Copy->SpinStateIndex[j] = Orig->SpinStateIndex[j];
    }
#endif

    for (int j = 0; j < Orig->NumDynVariables; j++) {
        Copy->State[j] = Orig->State[j];
    }
    for (int j = 0; j < 3 * Orig->NumVertex; j++) {
        Copy->Pos[j] = Orig->Pos[j];
        Copy->Vel[j] = Orig->Vel[j];
        Copy->Acc[j]    = Orig->Acc[j];
        Copy->AuxVel[j] = Orig->AuxVel[j];
        Copy->AccPN[j]  = Orig->AccPN[j];
#ifdef USE_PN_SPIN
	Copy->Spin_S[j]    = Orig->Spin_S[j];
	Copy->SpinState[j] = Orig->SpinState[j];
#endif
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
                                 double *err, double *yscal) {

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
    result[N - 2] = gbsS[0][N - 1];
    err[N - 2] = fabs(gbsS[0][N - 1] - gbsS[0][N - 1 - 1]) / yscal[N - 1];
}

int extrapolate_all_variables(int korder ) {

    int status = -1;

    struct RegularizedRegion *ThisR = NULL;
    int *SubStepDivision = calloc(korder, sizeof(int));

    int counter = 0;
    int NumDynVariables = 0, NumSpinVariables=0;
    for (int i = 0; i < korder; i++) {
        ThisR = ComputationToDoList.ComputationalTask[i].ThisR;
        NumDynVariables  = ThisR->NumDynVariables;
#ifdef USE_PN_SPIN
	NumSpinVariables  = ThisR->NumSpinVariables;
#endif
        SubStepDivision[counter] = i;
        counter++;
    }

    int NumPos  = 3*(ThisR->NumVertex-1);
    int NumVel  = NumPos;
#ifdef USE_PN_SPIN
    int NumSpin = 3*ThisR->NumVertex;
#endif

    double PosScale[korder], VelScale[korder], BScale[korder], TScale[korder];
#ifdef USE_PN_SPIN
    double SpinScale[korder];
#endif

    for (int i = 0; i < korder; i++) {
	PosScale[i]  = 0.0;
	VelScale[i]  = 0.0;
	BScale[i]    = 0.0;
	TScale[i]    = 0.0;
#ifdef USE_PN_SPIN
	SpinScale[i] = 0.0;
#endif
    }

    for (int v = 0; v < NumPos; v++) {
	for (int i = 0; i < korder; i++) {
		if( fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v] ) > PosScale[i] ) { PosScale[i] = fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v]); }
	}
    }
    for (int v = NumPos; v < NumPos+NumVel; v++) {
	for (int i = 0; i < korder; i++) {
               if( fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v] ) > VelScale[i] ) { VelScale[i] = fabs(ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v]); }
	}
    }
    int v = NumDynVariables-2;
    for (int i = 0; i < korder; i++) {
	if( fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v] ) > BScale[i] ) { BScale[i] = fabs(ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v]); }
    }
    v = NumDynVariables-1;
    for (int i = 0; i < korder; i++) {
	if( fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v] ) > TScale[i] ) { TScale[i] = fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v]); }
    }
#ifdef USE_PN_SPIN
    for (int v = 0; v < NumSpinVariables; v++) {
	for (int i = 0; i < korder; i++) {
		if( fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].SpinState[v] ) > SpinScale[i] ) { SpinScale[i] = fabs( ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].SpinState[v]); }
        }
    }
#endif

    for (int v = 0; v < NumDynVariables; v++) {
        for (int i = 0; i < korder; i++) {
            y[i] = ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].State[v];
	    if( v == NumDynVariables-1 ) { yscal[i] = TScale[i]; }
	    if( v == NumDynVariables-2 ) { yscal[i] = BScale[i]; }
	    if( v < NumPos ) { yscal[i] = PosScale[i]; }
	    if( v >= NumPos && v < NumPos+NumVel ) { yscal[i] = VelScale[i]; }

        }
        extrapolate_single_variable(y, korder, result[v], extrapolation_error[v], yscal);
    }
#ifdef USE_PN_SPIN
    for (int v = NumDynVariables; v < NumDynVariables+NumSpinVariables; v++) {
	for (int i = 0; i < korder; i++) {
		y[i] = ComputationToDoList.CopyOfSingleRegion[SubStepDivision[i]].SpinState[v-NumDynVariables];
		yscal[i] = SpinScale[i];
	}
	extrapolate_single_variable(y, korder, result[v], extrapolation_error[v], yscal);
    }
#endif

    for (int i = 0; i < NumDynVariables; i++) {
        errorbuf[i]  = extrapolation_error[i][korder - 2];
        resultbuf[i] = result[i][korder - 2];
    }
#ifdef USE_PN_SPIN
    for (int i = NumDynVariables; i < NumDynVariables+NumSpinVariables;i++) {
	errorbuf[i]  = extrapolation_error[i][korder - 2];
	resultbuf[i] = result[i][korder - 2];
    }
#endif

    double gbs_error = -1;
    double gbs_tolerance = ComputationToDoList.ComputationalTask[SubStepDivision[0]].ThisR->gbs_tolerance;

    int max_err_i = -1;
    for (int j = 0; j < NumDynVariables+NumSpinVariables; j++) {
        if (errorbuf[j] >= gbs_error) {
            gbs_error = errorbuf[j];
            max_err_i = j;
        }
    }

    double Hnextfac = 0.94 * pow(0.65 / (gbs_error / gbs_tolerance), 1.0 / (2 * korder + 1));

    if (gbs_error <= gbs_tolerance) {

        status = 1;

        for (int v = 0; v < NumDynVariables; v++) {
            ThisR->State[v] = resultbuf[v];
        }
#ifdef USE_PN_SPIN
	for (int v = NumDynVariables; v < NumDynVariables+NumSpinVariables; v++) {
	    ThisR->SpinState[v-NumDynVariables] = resultbuf[v];
	}
#endif
        if (Hnextfac > 5.0) { Hnextfac = 5.0; }

        ThisR->Hstep *= Hnextfac;

        into_Cartesian_from_MST(ThisR);
#ifdef USE_PN_SPIN
	from_SpinState_to_Spin_S(ThisR);
#endif

    } else {
        status = 0;
        if (Hnextfac > 0.7) { Hnextfac = 0.7; }
        ThisR->Hstep *= Hnextfac;
    }
    free(SubStepDivision);

    return status;


}

// Integrate the subsystems using AR-MST for the desired time interval
void run_integrator(struct RegularizedRegion *R, double time_interval, double *end_time, int *stopping_condition_occurred)
{
    int initial_collision = check_for_initial_stopping_condition(R);
    if (initial_collision == 1)
    {
        *stopping_condition_occurred = 1;
        *end_time = 0.0;
        if (MSTAR_verbose == 1)
        {
            printf("MSTAR -- initial stopping condition \n");
        }
        return;
    }


    // Needed variables
    double dt = get_timestep_estimate( R );
    double upper_dt = time_interval * 1.0e-6;
    if (dt >= upper_dt) // Want to make sure dt is not too long
    {
        dt = upper_dt;
    }
    double lower_dt = 1.0e-15; 
    if (dt <= lower_dt) // Want to make sure dt is not too small
    {
        dt = lower_dt;
    }
    
    double Hstep;
    double mst_time;

    // Initialize the system to integrate
    into_CoM_frame(R);
    PrimMST(R);

    compute_total_energy(R);
    R->Hstep = dt * R->U;
    R->State[R->TimeIndex] = 0;

    double percent = time_interval / 50;
    int epoch = 0;

    // Now proceed to integration
    int not_finished = 1;
    int endtime_iteration_phase = 0;
    int redo = 1;
    int status;

    int i_sc;

    int i_print=0;
    double N_print=10;
    double f_time;

//    check_MSTAR_system_for_distant_bodies(R);


    do {

        do {

            redo = 0;

            // Divide the systems and substep divisions for tasks
            divide_computational_load(R);

            // Perform the leapfrog tasks
            for (int ctask = 0; ctask < ComputationToDoList.NumberOfComputationalTasks; ctask++) {

                    // copy the system in order not to overwrite the original one
                    struct RegularizedRegion *ThisR = ComputationToDoList.ComputationalTask[ctask].ThisR;
                    copy_regularized_region( ThisR, &ComputationToDoList.CopyOfSingleRegion[ctask]);

                    // leapfrog the substep division
                    Hstep = ComputationToDoList.CopyOfSingleRegion[ctask].Hstep;
                    mst_leapfrog(&ComputationToDoList.CopyOfSingleRegion[ctask], Hstep, ComputationToDoList.ComputationalTask[ctask].NumberOfSubsteps);

            }

        } while (redo);

        // Now preform the GBS extrapolation, check convergence and obtain next
        // timestep

	status = extrapolate_all_variables(KMAX);

        time_t wall_time_current;
        time(&wall_time_current);
        double wall_time_diff_s = difftime(wall_time_current,wall_time_start);
        if (wall_time_diff_s > wall_time_max_s)
        {
            printf("mst.c -- wall time of %g s has exceeded allowed maximum of %g s -- system index %d \n",wall_time_diff_s,wall_time_max_s,system_index);
            error_code = 36;
            longjmp(jump_buf,1);
            return;
        }
    
        if (status == 1) {

            // Check whether we are ready
            mst_time = R[0].State[R[0].TimeIndex];
            R[0].time = mst_time;
            if ( fabs(mst_time - time_interval) / time_interval < R->output_time_tolerance) {
               not_finished = 0;
            } else {

                compute_total_energy(R);

		if( endtime_iteration_phase ){
			R->Hstep = R->U * (time_interval - mst_time);
		} else {
                	if (mst_time > time_interval) {
				endtime_iteration_phase = 1;
                    		R->Hstep = -R->U * (mst_time - time_interval);
                	} else if (R->Hstep / R->U + mst_time > time_interval) {
                    		R->Hstep = R->U * (time_interval - mst_time);
                	}
		}

#ifdef USE_PN_SPIN	// a useful two-body system debug
#if 0
		printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", mst_time, R->Pos[0],R->Pos[1],R->Pos[2],R->Pos[3],R->Pos[4],R->Pos[5],\
		R->Vel[0],R->Vel[1],R->Vel[2],R->Vel[3],R->Vel[4],R->Vel[5],R->Spin_S[0],R->Spin_S[1],R->Spin_S[2],R->Spin_S[3],R->Spin_S[4],R->Spin_S[5] );
#endif
#else
#if 0
		printf("%e %e %e %e %e %e %e %e %e %e %e %e %e\n", mst_time, R->Pos[0],R->Pos[1],R->Pos[2],R->Pos[3],R->Pos[4],R->Pos[5],\
		R->Vel[0],R->Vel[1],R->Vel[2],R->Vel[3],R->Vel[4],R->Vel[5] );
#endif
#endif
                  f_time = mst_time/time_interval;
                  if ( ((int) (f_time*N_print) ) == i_print && MSTAR_verbose == 1)
                  {
                        printf("MSTAR -- t %g completed %.1f %%\n",mst_time,f_time*100.0 );
                        i_print+=1;
                 }

                /* Stopping conditions */
                int possible_stopping_condition, *topping_condition_occurred;
                double Delta_t_stopping_condition;
                stopping_condition_function(R, &possible_stopping_condition, stopping_condition_occurred, &Delta_t_stopping_condition);

                if (*stopping_condition_occurred == 1)
                {
                    not_finished = 0;
                    if (MSTAR_verbose == 1)
                    {
                        printf("Stopping condition at t=%g \n",mst_time);
                    }
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
                        printf("mstar --  Giving up stopping condition search t %g\n",mst_time);
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
                epoch = (int)floor(mst_time / percent);
                if (epoch != old_epoch) {
                    PrimMST(&R[0]);
                }
            }
        }

    } while (not_finished);

    out_of_CoM_frame(R);

    *end_time = mst_time;

fflush(stdout);

}

void allocate_regularized_region(struct RegularizedRegion *S, int N) {

    S->gbs_tolerance = GBSTOL;

    S->output_time_tolerance = 1e-4;
    S->stopping_condition_tolerance = 1e-6;
    S->NumVertex = N;

    S->NumEdge = (N * (N - 1)) / 2;
    S->NumDynVariables = 6 * (N - 1) + 2;
#ifdef USE_PN_SPIN
    S->NumSpinVariables = 3*N;
#endif

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

#ifdef USE_PN_SPIN
    S->SpinState     = calloc(3 * N, sizeof(double));
    S->SpinStateIndex     = calloc(N, sizeof(int));
    S->Spin_S     = calloc(3 * N, sizeof(double));
    S->AuxSpin_S     = calloc(3*N, sizeof(double));
    S->Spin_dS_PN = calloc(3 * N, sizeof(double));
#endif

    S->Vertex = calloc(N, sizeof(struct GraphVertex));

    double safetybuf = 1.05;
    int NumLocalEdge = (int)floor(safetybuf * S->NumEdge);
    S->LocalEdge = calloc(NumLocalEdge, sizeof(struct GraphEdge));
    S->LocalEdgeSubset = calloc(NumLocalEdge, sizeof(struct GraphEdge));
    S->EdgeInMST = calloc(N - 1, sizeof(struct GraphEdge));

    S->Radius = calloc(N, sizeof(double));
    S->Stopping_Condition_Mode = calloc(N, sizeof(int));
    S->Stopping_Condition_Partner = calloc(N, sizeof(int));
    S->Stopping_Condition_Roche_Lobe_Radius = calloc(N, sizeof(double));

    S->MSEIndex = calloc(N, sizeof(int));
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

#ifdef USE_PN_SPIN
    int NumSpinVariables = 3*MaxNumPart;
#endif

#ifdef USE_PN_SPIN
    errorbuf      = calloc(NumDynVariables+NumSpinVariables, sizeof(double));
    resultbuf     = calloc(NumDynVariables+NumSpinVariables, sizeof(double));
    extrapolation_error         = calloc(NumDynVariables+NumSpinVariables, sizeof(double *));
    result        = calloc(NumDynVariables+NumSpinVariables, sizeof(double *));
#else
    errorbuf 	  = calloc(NumDynVariables, sizeof(double));
    resultbuf 	  = calloc(NumDynVariables, sizeof(double));
    extrapolation_error 	  = calloc(NumDynVariables, sizeof(double *));
    result 	  = calloc(NumDynVariables, sizeof(double *));
#endif

    y     = calloc(korder, sizeof(double));
    yscal = calloc(korder, sizeof(double));

#ifdef USE_PN_SPIN
    for (int i = 0; i < NumDynVariables+NumSpinVariables; i++) {
        extrapolation_error[i] = calloc(korder - 1, sizeof(double));
        result[i] = calloc(korder - 1, sizeof(double));
    }
#else
    for (int i = 0; i < NumDynVariables; i++) {
	extrapolation_error[i] = calloc(korder - 1, sizeof(double));
	result[i] = calloc(korder - 1, sizeof(double));
    }
#endif

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

#ifdef USE_PN_SPIN
    free(   R->Spin_dS_PN);
    free(   R->Spin_S);
    free(R->AuxSpin_S);
#endif

    free(R->Radius);
    free(R->Stopping_Condition_Mode);
    free(R->Stopping_Condition_Partner);
    free(R->Stopping_Condition_Roche_Lobe_Radius);

    free(R->MSEIndex);

    free(R);

    free(gbsS);
    
}

void stopping_condition_function(struct RegularizedRegion *R, int *possible_stopping_condition, int *stopping_condition_occurred, double *Delta_t_min)
{
    int istart, istop = R->NumVertex, jstop = R->NumVertex;
    istart = 0;

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
    int *MSEIndex = R->MSEIndex;
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
        const int MSEIndexi = MSEIndex[i];
        for (int j = i + 1; j < jstop; j++)
        {
            r2 = 0;
            const double mj = Mass[j];
            const double radiusj = Radius[j];
            const int Lj = Vertex[j].level;
            const int modej = Stopping_Condition_Mode[j];
            const int MSEIndexj = MSEIndex[j];

            int proximity = 0;
            if (abs(Li - Lj) <= Nd) {
                proximity =
                    //check_relative_proximity_ND_2(j, i, R, &d, path, sign);
                    check_relative_proximity(j, i, Nd, R, &d, path, sign);
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
            else if (modei == 1 && modej == 0) // RLOF + col
            {
                r_criti = radiusi / fq_RLOF_Eggleton(mi,mj);
                r_critj = radiusj;
                r_crit = CV_max( r_criti, r_critj );
            }
            else if (modei == 0 && modej == 1) // col + RLOF
            {
                r_criti = radiusi;
                r_critj = radiusj / fq_RLOF_Eggleton(mj,mi);
                r_crit = CV_max( r_criti, r_critj );
            }
            else
            {
                printf("mst.c -- stopping_condition_function -- invalid (combination) of Stopping_Condition_Mode %d and %d for particles %d and %d, respectively\n",modei,modej,i,j);
                error_code = 31;
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
                R->Stopping_Condition_Partner[i] = MSEIndexj;
                R->Stopping_Condition_Partner[j] = MSEIndexi;
                
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
            int interpolation_method = 0;
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
    for (int i = 0; i < R->NumVertex; i++) {
        for (int k = 0; k < 3; k++) {
            R->Pos[3 * i + k] += R->CoM_Pos[k];
            R->Vel[3 * i + k] += R->CoM_Vel[k];
        }
    }
}

int check_for_initial_stopping_condition(struct RegularizedRegion *R)
{
    int istart, istop = R->NumVertex, jstop = R->NumVertex;
    istart = 0;

    double dr[3], r, r2;
    double dv[3], da[3];
    double r_crit,r_crit_p2;
    double r_criti, r_critj;
    double rdotv,vdotv;
    double determinant;
    double Delta_t;
    
    int stopping_condition_occurred = 0;
    
    double *Pos = R->Pos;
    double *Vel = R->Vel;
    double *Mass = R->Mass;
    double *Radius = R->Radius;
    int *MSEIndex = R->MSEIndex;
    int *Stopping_Condition_Mode = R->Stopping_Condition_Mode;

    for (int i = istart; i < istop; i++) 
    {
        R->Stopping_Condition_Partner[i] = -1; /* default: no stopping condition partner */
        R->Stopping_Condition_Roche_Lobe_Radius[i] = -1;
    }

    for (int i = istart; i < istop; i++) 
    {
        const double mi = Mass[i];
        const double radiusi = Radius[i];
        const int modei = Stopping_Condition_Mode[i];
        const int MSEIndexi = MSEIndex[i];
        for (int j = i + 1; j < jstop; j++)
        {
            r2 = 0;
            const double mj = Mass[j];
            const double radiusj = Radius[j];
            const int modej = Stopping_Condition_Mode[j];
            const int MSEIndexj = MSEIndex[j];

            for (int k = 0; k < 3; k++) 
            {
                dr[k] = Pos[3 * j + k] - Pos[3 * i + k];
                r2 += dr[k] * dr[k];
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
            else if (modei == 1 && modej == 0) // RLOF + col
            {
                r_criti = radiusi / fq_RLOF_Eggleton(mi,mj);
                r_critj = radiusj;
                r_crit = CV_max( r_criti, r_critj );
            }
            else if (modei == 0 && modej == 1) // col + RLOF
            {
                r_criti = radiusi;
                r_critj = radiusj / fq_RLOF_Eggleton(mj,mi);
                r_crit = CV_max( r_criti, r_critj );
            }
            else
            {
                printf("mst.c -- check_for_initial_stopping_condition -- invalid (combination) of Stopping_Condition_Mode %d and %d for particles %d and %d, respectively\n",modei,modej,i,j);
                error_code = 32;
            }
            
            r_crit_p2 = r_crit * r_crit;
            r = sqrt(r2);
            
            if (r < r_crit)
            {
                stopping_condition_occurred = 1;
                R->Stopping_Condition_Partner[i] = MSEIndexj;
                R->Stopping_Condition_Partner[j] = MSEIndexi;
                
                if (modei == 1 && modej == 1)
                {
                    if (r < r_criti)
                    {
                        R->Stopping_Condition_Roche_Lobe_Radius[i] = r_criti * fq_RLOF_Eggleton(mi,mj);
                    }
                    else if (r < r_critj)
                    {
                        R->Stopping_Condition_Roche_Lobe_Radius[j] = r_critj * fq_RLOF_Eggleton(mj,mi);
                    }
                }
            }
        }
    }
    return stopping_condition_occurred;
}

//extern "C"
//{
    
#define SERIAL
//#define PARALLEL
#define PRIM
//This appears broken in this version...
//#define PRIM_DIVIDE_CONQUER

#define GLOBAL_ND 2
//#define GCONST 43007.1

#define KMAX 8
#define MAXSORTTASK 8

//#define M_PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214
//#define SPEEDOFLIGHT 299792.0

/* The two lines below: modification from A.S.H. */
#define GCONST (double) 4.0*M_PI*M_PI
#define SPEEDOFLIGHT (double) 63239.72638679138

// structs

//#undef USE_PN
#define USE_PN


/* End modifications for MSE */

#ifndef __MSTAR
#define __MSTAR

enum velocity_type { PHYSICAL, AUXILIARY };
enum variable_type { POSITION, VELOCITY, OTHER };

#ifdef IGNORE
extern int Ntask;
extern int ThisTask;

extern int NumGbsGroup;
extern int NumTaskPerGbsGroup;
extern int ThisGbsGroup;
extern int ThisTask_in_GbsGroup;

extern int NumTaskPerSortGroup;
extern int NumSortGroup;
extern int ThisSortGroup;
extern int ThisTask_in_SortGroup;
#endif

struct OctreeNode {
    int ChildExists[8];
    struct OctreeNode *Children[8];
    struct OctreeNode *Parent;
    double      Center[3];
    double      HalfSize;
    int *ObjectList;
    int Npart;
    int level;
    int IsLeaf;
    int ID;
};


struct WeightIndex {
        double weight;
        int index;
};
struct WeightIndex2 {
        double weight;
        int index1;
	int index2;
	int id1;
	int id2;
};

struct GraphEdge {
        int id;
        int vertex1;
        int vertex2;
        float weight;
};

struct GraphVertex {
        int id;
        int degree;
        int type;
        int parent;
        int edgetoparent;
        int level;
        int inMST;
	int cell[3];
};

struct RegularizedRegion {

        int NumVertex;
        int NumEdge;
	int NumDynVariables;

        double *Pos;
        double *Vel;
        double *Mass;
        double *Acc;
        double *InitialState;
        double *State;
        double *MSTedgeAcc;
        double *AuxVel;
        double *AuxEdgeVel;
	double *AccPN;
        double *MSTedgeAcc_PN;

	struct GraphEdge *Edge;
        struct GraphEdge *LocalEdge;
	struct GraphEdge *LocalEdgeSubset;
        struct GraphVertex *Vertex;
        int *MSTedgeList;
	struct GraphEdge *EdgeInMST;

        double T;
        double U;
        double H;
        double B;
        double time;

	double Udot;

	int BIndex;
	int TimeIndex;

        double CoM_Pos[3];
        double CoM_Vel[3];

	// relative error tolerance
	double gbs_tolerance;
	double output_time_tolerance;

	int korder;

	double Hstep;

    /* Stopping conditions */
    double *Radius;
    int *Stopping_Condition_Mode;
    int *Stopping_Condition_Partner; 
    double stopping_condition_tolerance;
    double *Stopping_Condition_Roche_Lobe_Radius;
    
    /* Identification (for interface with MSE) */
    int *Index;
    
};

struct ComTask {
	int GroupToPerformTask;
	struct RegularizedRegion *ThisR;
	int NumberOfSubsteps;
	double *ResultState;
};

struct ToDoList {
	int NumberOfComputationalTasks;
	struct ComTask *ComputationalTask;
	struct RegularizedRegion *CopyOfSingleRegion;
};

//extern double E0;

// functions

void compute_Post_Newtonian_Acc(struct RegularizedRegion *R, double *Vel);
int check_relative_proximity( int v1, int v2, const int Nd, struct RegularizedRegion *R, int *d, int *path, int *sign);

void stopping_condition_function(struct RegularizedRegion *R, int *possible_stopping_condition, int *stopping_condition_occurred, double *Delta_t_min);
double fq_RLOF_Eggleton(double m1, double m2);

// misc

//extern int ok_steps;
//extern int failed_steps;

#endif
//void allocate_armst_structs(struct RegularizedRegion **R, int MaxNumPart);
//void initialize_mpi_or_serial(void);

#ifndef CV_max
    #define CV_max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
//}

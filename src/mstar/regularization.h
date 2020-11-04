#include "../parameters.h"

#define GLOBAL_ND 2
#define KMAX 8

// structs

//#undef USE_PN
#define USE_PN
//#undef USE_PN_SPIN
#define USE_PN_SPIN

#ifndef __MSTAR
#define __MSTAR

enum velocity_type { PHYSICAL, AUXILIARY };
enum variable_type { POSITION, VELOCITY, OTHER };

struct WeightIndex {
        double weight;
        int index;
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
#ifdef USE_PN_SPIN
	int NumSpinVariables;
#endif
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

#ifdef USE_PN_SPIN
	double *SpinState;	
	int    *SpinStateIndex; 
	double *Spin_S;		
	double *Spin_dS_PN;	
	double *AuxSpin_S;      
#endif
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
    int *MSEIndex;
    
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


// functions

#ifdef USE_PN_SPIN
void from_Spin_S_to_SpinState( struct RegularizedRegion *R );
void from_SpinState_to_Spin_S( struct RegularizedRegion *R );
#endif
#ifdef USE_PN_SPIN
void compute_Post_Newtonian_Acc(struct RegularizedRegion *R, double *Vel, double *Spin );
#else
void compute_Post_Newtonian_Acc(struct RegularizedRegion *R, double *Vel );
#endif

int check_relative_proximity( int v1, int v2, const int Nd, struct RegularizedRegion *R, int *d, int *path, int *sign);
void stopping_condition_function(struct RegularizedRegion *R, int *possible_stopping_condition, int *stopping_condition_occurred, double *Delta_t_min);
double fq_RLOF_Eggleton(double m1, double m2);
void out_of_CoM_frame(struct RegularizedRegion *R);
int check_for_initial_stopping_condition(struct RegularizedRegion *R);

#endif

#ifndef CV_max
    #define CV_max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

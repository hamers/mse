/* MSE */

#include "evolve.h"
#include "tools.h"

extern "C"
{

double compute_orbital_period_from_semimajor_axis(double M, double a)
{
	return TWOPI*sqrt(a*a*a/(CONST_G*M));
}

double compute_semimajor_axis_from_orbital_period(double M, double P)
{
    double temp = P/(TWOPI);
    return pow( temp*temp* CONST_G*M, 1.0/3.0);
}

double generate_random_number_between_zero_and_unity()
{
    return ((double) rand() / (RAND_MAX));
}

int sample_from_3d_maxwellian_distribution(double sigma, double v[3])
{
    double u1,u2,s,theta;
    for (int k=1; k<3; k++)
    {
        u1 = generate_random_number_between_zero_and_unity();
        u2 = generate_random_number_between_zero_and_unity();
        s = sigma*sqrt(-2.0*log(1.0 - u1));
        theta = TWOPI*u2;
        v[2*k-2] = s*cos(theta);
        v[2*k-1] = s*sin(theta);
    }    

    return 0;
}

double sample_from_y_times_maxwellian_distribution(double sigma)
{
    /* Sample random variable y from a distribution dN/dy \propto y*Exp[-y^2/(2 sigma^2)] */
    
    double x = generate_random_number_between_zero_and_unity();
    double y = sigma*sqrt( -2.0*log(1.0 - x) );
    
    return y;

}

double sample_from_normal_distribution(double mu, double sigma)
{
    //std::random_device rd{};
    //std::mt19937 gen{rd()};
    //return std::normal_distribution<> d{mu,sigma};

    double u1 = generate_random_number_between_zero_and_unity();
    double u2 = generate_random_number_between_zero_and_unity();
    double result = mu + sigma * sqrt(-2.0*log(u1)) * cos(TWOPI*u2);
    return result;
}

int sample_spherical_coordinates_unit_vectors_from_isotropic_distribution(double r_hat_vec[3], double theta_hat_vec[3], double phi_hat_vec[3])
{
    double x1 = generate_random_number_between_zero_and_unity();
    double x2 = generate_random_number_between_zero_and_unity();
    double theta = acos( 2.0*x1 - 1.0 ); /* inclination */
    double phi = TWOPI*x2; /* azimuthal angle */
    
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    r_hat_vec[0] = sin_theta*cos_phi;
    r_hat_vec[1] = sin_theta*sin_phi;
    r_hat_vec[2] = cos_theta;
    
    theta_hat_vec[0] = cos_theta*cos_phi;
    theta_hat_vec[1] = cos_theta*sin_phi;
    theta_hat_vec[2] = -sin_theta;
    
    phi_hat_vec[0] = -sin_phi;
    phi_hat_vec[1] = cos_phi;
    phi_hat_vec[2] = 0.0;
    
    return 0;
}

double sample_from_power_law_distribution(double alpha, double y_lower, double y_upper)
{
    /* Sample random variable y from dN/dy ~ y^y_\alpha */

    double x = generate_random_number_between_zero_and_unity();
    double y;
    if (alpha == -1.0)
    {
        y = y_lower*pow(y_upper/y_lower,x);
    }
    else
    {
        y = pow( x*(pow(y_upper,alpha+1.0) - pow(y_lower,alpha+1.0)) + pow(y_lower,alpha+1.0), 1.0/( alpha + 1.0) );
    }
    
    return y;

}

double sample_from_Kroupa_93_imf()
{
    double m = 0;
    double x = generate_random_number_between_zero_and_unity();
    
    if (x >= 0.0 and x < kroupa_x1)
    {
        m = pow( x*kroupa_alpha1_plus_1/kroupa_C1 + kroupa_m1_pow_alpha1_plus_one, kroupa_alpha1_plus_1_pm1);
    }
    else if (x >= kroupa_x1 and x < kroupa_x2)
    {
        m = pow( (x - kroupa_x1)*kroupa_alpha2_plus_1/kroupa_C2 + kroupa_m2_pow_alpha2_plus_one, kroupa_alpha2_plus_1_pm1);
    }
    else if (x >= kroupa_x2 and x <= kroupa_x3)
    {
        m = pow( (x - kroupa_x2)*kroupa_alpha3_plus_1/kroupa_C3 + kroupa_m3_pow_alpha3_plus_one, kroupa_alpha3_plus_1_pm1);
    }
    else
    {
        printf("tools.cpp -- ERROR in sample_from_Kroupa_93_imf\n");
        exit(-1);
    }
    
    return m;
}



/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vec[3])
{
    for (int i=0; i<3; i++)
    {
        h_tot_vec[i] = 0.0;
    }
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1)
        {
            for (int i=0; i<3; i++)
            {
                h_tot_vec[i] += p->h_vec[i];
            }
        }
    }
    
    return 0;
//    printf("compute_h_tot_vector %g %g %g\n",h_tot_vec[0],h_tot_vec[1],h_tot_vec[2]);
}
    
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    double cos_INCL = cos(inclination);
    double sin_INCL = sin(inclination);
    double cos_AP = cos(argument_of_pericenter);
    double sin_AP = sin(argument_of_pericenter);
    double cos_LAN = cos(longitude_of_ascending_node);
    double sin_LAN = sin(longitude_of_ascending_node);
           
    double h = compute_h_from_a(child1_mass, child2_mass, semimajor_axis, eccentricity);

    *e_vec_x = eccentricity*(cos_LAN*cos_AP - sin_LAN*sin_AP*cos_INCL);
    *e_vec_y = eccentricity*(sin_LAN*cos_AP + cos_LAN*sin_AP*cos_INCL);
    *e_vec_z = eccentricity*(sin_AP*sin_INCL);
    
    *h_vec_x = h*sin_LAN*sin_INCL;
    *h_vec_y = -h*cos_LAN*sin_INCL;
    *h_vec_z = h*cos_INCL;

    return 0;
}

int compute_orbital_vectors_from_orbital_elements_unit(double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_hat_vec_x, double *e_hat_vec_y, double *e_hat_vec_z, double *h_hat_vec_x, double *h_hat_vec_y, double *h_hat_vec_z)
{
    double cos_INCL = cos(inclination);
    double sin_INCL = sin(inclination);
    double cos_AP = cos(argument_of_pericenter);
    double sin_AP = sin(argument_of_pericenter);
    double cos_LAN = cos(longitude_of_ascending_node);
    double sin_LAN = sin(longitude_of_ascending_node);
           
    *e_hat_vec_x = (cos_LAN*cos_AP - sin_LAN*sin_AP*cos_INCL);
    *e_hat_vec_y = (sin_LAN*cos_AP + cos_LAN*sin_AP*cos_INCL);
    *e_hat_vec_z = (sin_AP*sin_INCL);
    
    *h_hat_vec_x = sin_LAN*sin_INCL;
    *h_hat_vec_y = -cos_LAN*sin_INCL;
    *h_hat_vec_z = cos_INCL;

    return 0;
}

int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node)
{
    #ifdef DEBUG
    printf("tools.cpp -- compute_orbital_elements_from_orbital_vectors -- e_vec %g %g %g h_vec %g %g %g\n",e_vec_x,e_vec_y,e_vec_z,h_vec_x,h_vec_y,h_vec_z);
    #endif
    
    double e_vec[3] = {e_vec_x,e_vec_y,e_vec_z};
    double h_vec[3] = {h_vec_x,h_vec_y,h_vec_z};
    double eccentricity_squared = norm3_squared(e_vec);
    *eccentricity = sqrt(eccentricity_squared);
    if (*eccentricity==0.0)
    {
        *eccentricity = epsilon;
    }
    double h = norm3(h_vec);
    *semimajor_axis = compute_a_from_h(child1_mass, child2_mass, h, *eccentricity);

//    double x_vec[3] = {1.0,0.0,0.0};
//    double y_vec[3] = {0.0,1.0,0.0};
//    double z_vec[3] = {0.0,0.0,1.0};

    double h_tot = norm3(h_tot_vec);
//    printf("h_tot %g x %g y %g z %g\n",h_tot,h_tot_vec[0],h_tot_vec[1],h_tot_vec[2]);
    double x_vec[3], y_vec[3], z_vec[3];
    for (int i=0; i<3; i++)
    {
        z_vec[i] = h_tot_vec[i]/h_tot;
    }

//    printf("test %g %g %g\n",z_vec[0],z_vec[1],z_vec[2]);
    z_vec[0] = 0.0;
    z_vec[1] = 0.0;
    z_vec[2] = 1.0;

    /* the above assumes that the total angular momentum vector does not change (i.e. no SNe effects etc.) */
    
    double f = 1.0/sqrt( z_vec[0]*z_vec[0] + z_vec[2]*z_vec[2] );
    x_vec[0] = z_vec[2]*f;
    x_vec[1] = 0.0;
    x_vec[2] = -z_vec[0]*f;
    cross3(z_vec,x_vec,y_vec);

    double cos_INCL = dot3(h_vec,z_vec)/h;
    if (cos_INCL >= 1.0)
    {
        cos_INCL = 1.0 - epsilon;
    }
    if (cos_INCL <= -1.0)
    {
        cos_INCL = -1.0 + epsilon;
    }

    double LAN_vec[3],LAN_vec_unit[3];
    cross3(z_vec,h_vec,LAN_vec);
    double LAN_vec_norm = norm3(LAN_vec);
    if (LAN_vec_norm==0.0)
    {
        LAN_vec_norm = 1.0;
    }

    double e_vec_unit[3],h_vec_unit[3];

    for (int i=0; i<3; i++)
    {
        LAN_vec_unit[i] = LAN_vec[i]/LAN_vec_norm;
        e_vec_unit[i] = e_vec[i]/(*eccentricity);
        h_vec_unit[i] = h_vec[i]/h;
    }

    double sin_LAN = dot3(LAN_vec_unit,y_vec);
    double cos_LAN = dot3(LAN_vec_unit,x_vec);

    double e_vec_unit_cross_h_vec_unit[3];
    cross3(e_vec_unit,h_vec_unit,e_vec_unit_cross_h_vec_unit);

    double sin_AP = dot3(LAN_vec_unit,e_vec_unit_cross_h_vec_unit);
    double cos_AP = dot3(LAN_vec_unit,e_vec_unit);

    *inclination = acos(cos_INCL);
    *argument_of_pericenter = atan2(sin_AP,cos_AP);
    *longitude_of_ascending_node = atan2(sin_LAN,cos_LAN);
    
    return 0;

}

int get_inclination_relative_to_parent(ParticlesMap *particlesMap, int index, double *inclination_relative_to_parent)
{

    if (index > particlesMap->size())
    {
        return -1;
    }
  
    Particle *p = (*particlesMap)[index];
    if (p->is_binary == false)
    {
        *inclination_relative_to_parent = 0.0;
        return 0;
    }
    if (p->parent == -1)
    {
        *inclination_relative_to_parent = 0.0;
        return 0;
    }

    Particle *parent = (*particlesMap)[p->parent];
    

    double *h1_vec = p->h_vec;
    double *h2_vec = parent->h_vec;

    double h1 = norm3(h1_vec);
    double h2 = norm3(h2_vec);
    
    double arg = dot3(h1_vec,h2_vec)/(h1*h2);
    if (arg <= -1.0)
    {
        arg = -1.0;
    }
    if (arg >= 1.0)
    {
        arg = 1.0;
    }
    
    *inclination_relative_to_parent = acos( arg );
    
    return 0;
}

void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly)
{
    double eccentric_anomaly;
    double eccentric_anomaly_next = mean_anomaly;
    double epsilon2 = 1e-10;
    double error = 2.0*epsilon2;
    int j = 0;
    while (error > epsilon2 || j < 15)
    {
        j += 1;
        eccentric_anomaly = eccentric_anomaly_next;
        eccentric_anomaly_next = eccentric_anomaly - (eccentric_anomaly - eccentricity*sin(eccentric_anomaly) - mean_anomaly)/(1.0 - eccentricity*cos(eccentric_anomaly));
        error = fabs(eccentric_anomaly_next - eccentric_anomaly);
    }
    *cos_eccentric_anomaly = cos(eccentric_anomaly);
    *sin_eccentric_anomaly = sin(eccentric_anomaly);
}

void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly)
{
    *cos_true_anomaly = (cos_eccentric_anomaly - eccentricity)/(1.0 - eccentricity*cos_eccentric_anomaly);
    *sin_true_anomaly = sqrt(1.0 - eccentricity*eccentricity)*sin_eccentric_anomaly/(1.0 - eccentricity*cos_eccentric_anomaly);
}

double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity)
{
    double cos_eccentric_anomaly,sin_eccentric_anomaly;
    double cos_true_anomaly,sin_true_anomaly;
    
    compute_eccentric_anomaly_from_mean_anomaly(mean_anomaly,eccentricity,&cos_eccentric_anomaly,&sin_eccentric_anomaly);
    compute_true_anomaly_from_eccentric_anomaly(cos_eccentric_anomaly,sin_eccentric_anomaly,eccentricity,&cos_true_anomaly,&sin_true_anomaly);
    double true_anomaly = atan2(sin_true_anomaly,cos_true_anomaly);

    return true_anomaly;
}

double compute_mean_anomaly_from_true_anomaly(double true_anomaly, double eccentricity)
{
    double cos_true_anomaly = cos(true_anomaly);
    double sin_true_anomaly = sin(true_anomaly);
    double den = 1.0/( 1.0 + eccentricity*cos_true_anomaly );
    
    double sin_eccentric_anomaly = sqrt(1.0 - eccentricity*eccentricity)*sin_true_anomaly*den;
    double cos_eccentric_anomaly = (eccentricity + cos_true_anomaly)*den;
    double eccentric_anomaly = atan2(sin_eccentric_anomaly,cos_eccentric_anomaly);
    double mean_anomaly = eccentric_anomaly - eccentricity*sin_eccentric_anomaly;
    
    return mean_anomaly;
}

void compute_true_anomaly_from_mean_anomaly_hyperbolic_orbit(double mean_anomaly, double eccentricity,double *cos_true_anomaly,double *sin_true_anomaly)
{
    double eccentric_anomaly;
    
    double fabs_mean_anomaly = fabs(mean_anomaly);
    double sign_mean_anomaly;

    sign_mean_anomaly = copysign( 1.0, mean_anomaly);
    
    double eccentric_anomaly_next;    
    
    if (fabs_mean_anomaly < 3.0*eccentricity)
    {
        double s1 = fabs_mean_anomaly/(eccentricity-1.0);
        double s2 = pow( 6.0*fabs_mean_anomaly, 1.0/3.0);
        eccentric_anomaly_next = sign_mean_anomaly*CV_min(s1,s2);
    }
    else
    {
        eccentric_anomaly_next = sign_mean_anomaly*log(1.0 + 2.0*fabs_mean_anomaly/eccentricity);
    }
    
    double epsilon = 1e-10;
    double error = 2.0*epsilon; /* to start: anything larger than epsilon */
    int j = 0;
    while (error > epsilon)
    {
        j += 1;
        eccentric_anomaly = eccentric_anomaly_next;
        eccentric_anomaly_next = eccentric_anomaly + (eccentric_anomaly - eccentricity*sinh(eccentric_anomaly) + mean_anomaly)/(eccentricity*cosh(eccentric_anomaly) - 1.0);
        error = fabs(eccentric_anomaly_next - eccentric_anomaly);
        
        if (j > 15)
        {
            //printf("test %d %g %g %g %g %g\n",j,mean_anomaly,eccentric_anomaly,eccentric_anomaly_next,error,epsilon);
            break;
        }
    }
    
    double tau = sqrt( (eccentricity+1.0)/(eccentricity-1.0) )*tanh(0.5*eccentric_anomaly); /* tan(true_anomaly/2) */
    double tau_sq = tau*tau;
    double temp = 1.0/(1.0 + tau_sq);
    
    *cos_true_anomaly = (1.0 - tau_sq)*temp;
    *sin_true_anomaly = 2.0*tau*temp;
    
    //printf("test %g %g\n",mean_anomaly,eccentric_anomaly);
}

double sample_random_true_anomaly(double eccentricity)//,int seed)
{
    double x = generate_random_number_between_zero_and_unity();
    double mean_anomaly = (2.0*x - 1.0)*M_PI;
    double true_anomaly = compute_true_anomaly_from_mean_anomaly(mean_anomaly,eccentricity);

    return true_anomaly;
}

void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3])
{
    double total_mass = child1_mass + child2_mass;
    
    double e = norm3(e_vec);
    double h = norm3(h_vec);

    if (e==0.0)
    {
        e = epsilon;
    }

    double e_vec_unit[3],q_vec_unit[3],q_vec[3];
    cross3(h_vec,e_vec,q_vec);
    double q = norm3(q_vec);

    int i;
    for (i=0; i<3; i++)
    {        
        e_vec_unit[i] = e_vec[i]/e;
        q_vec_unit[i] = q_vec[i]/q;        
    }
    
    double e_p2 = e*e;
    double j_p2 = 1.0 - e_p2;
   
    double a = h*h*total_mass/( CONST_G*child1_mass*child2_mass*child1_mass*child2_mass*j_p2 );

    double cos_f = cos(true_anomaly);
    double sin_f = sin(true_anomaly);
    
    double r_norm = a*j_p2/(1.0 + e*cos_f);
    double v_norm = sqrt( CONST_G*total_mass/(a*j_p2) );
    
    for (i=0; i<3; i++)
    {
        r[i] = r_norm*( cos_f*e_vec_unit[i] + sin_f*q_vec_unit[i]);
        v[i] = v_norm*( -sin_f*e_vec_unit[i] + (e + cos_f)*q_vec_unit[i] );
            
        #ifdef DEBUG
        printf("tools.cpp -- from_orbital_vectors_to_cartesian -- i %d r[i] %g v[i] %g\n",i,r[i],v[i]);
        #endif
    }
}

void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3], double *true_anomaly)
{
    double total_mass = child1_mass + child2_mass;
       
    double v_dot_v = dot3(v,v);
    double r_dot_v = dot3(r,v);
    double r_norm = norm3(r);
    for (int i=0; i<3; i++)
    {
        e_vec[i] = (r[i]*v_dot_v - v[i]*r_dot_v)/(CONST_G*total_mass) - r[i]/r_norm;
    }

    double mu = child1_mass*child2_mass/total_mass;
    cross3(r,v,h_vec);
    for (int i=0; i<3; i++)
    {
        h_vec[i] *= mu;
    }
    
    double e_vec_unit[3],h_vec_unit[3],q_vec_unit[3];
    double h = norm3(h_vec);
    double e = norm3(e_vec);
    if (e==0.0)
    {
        e = epsilon;
    }
    for (int i=0; i<3; i++)
    {
        h_vec_unit[i] = h_vec[i]/h;
        e_vec_unit[i] = e_vec[i]/e;
    }
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);
    
    double cos_TA = dot3(r,e_vec_unit)/r_norm;
    double sin_TA = dot3(r,q_vec_unit)/r_norm;
    *true_anomaly = atan2( sin_TA, cos_TA );
    
    #ifdef DEBUG
    printf("tools.cpp -- from_cartesian_to_orbital_vectors -- e_vec %g %g %g h_vec %g %g %g TA %g\n",e_vec[0],e_vec[1],e_vec[2],h_vec[0],h_vec[1],h_vec[2],*true_anomaly);
    #endif
}

void compute_semimajor_axis_and_eccentricity_from_orbital_vectors(double m1, double m2, double e_vec[3], double h_vec[3], double *semimajor_axis, double *eccentricity)
{
    double e_P2 = norm3_squared(e_vec);
    double j_P2 = 1.0 - e_P2;
    double h_P2 = norm3_squared(h_vec);
    
    *eccentricity = sqrt(e_P2);
    *semimajor_axis = h_P2*(m1 + m2)/( CONST_G* m1 * m1 * m2 * m2 * j_P2 );
}

void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3])
{
    r[0] = p->R_vec[0];
    r[1] = p->R_vec[1];
    r[2] = p->R_vec[2];
    v[0] = p->V_vec[0];
    v[1] = p->V_vec[1];
    v[2] = p->V_vec[2];

}
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3])
{
    p->R_vec[0] = r[0];
    p->R_vec[1] = r[1];
    p->R_vec[2] = r[2];
    p->V_vec[0] = v[0];
    p->V_vec[1] = v[1];
    p->V_vec[2] = v[2];
}

void get_unit_vector(double vec[3], double vec_unit[3])
{
    double v = norm3(vec);
    if (v == 0.0)
    {
        v = epsilon;
    }
    
    for (int i=0; i<3; i++)
    {
        vec_unit[i] = vec[i]/v;
    }
}

void copy_particlesMap(ParticlesMap *source, ParticlesMap *target)
{
    (*target).clear();
    int index;
    int j=0;
    bool is_binary;
    ParticlesMapIterator it_p;
    for (it_p = source->begin(); it_p != source->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        //is_binary = p->is_binary;
        index = p->index;
        //printf("copy_particles index %d j %d\n",index,j);
        (*target)[index] = p;
        j++;
        //p->parent = -1;
    }
    //printf("old\n");
    //print_system(source);
    //printf("new\n");
    //print_system(target);
    
    //printf("copy_particlesMap N_s %d N_t %d \n",source->size(),target->size());    
}

void copy_all_body_properties(Particle *source, Particle *target)
{
    target->mass = source->mass;
    target->sse_initial_mass = source->sse_initial_mass;
    target->stellar_type = source->stellar_type;
    
    target->epoch = source->epoch;
    target->age = source->age;
    target->sse_main_sequence_timescale = source->sse_main_sequence_timescale;

    target->radius = source->radius;
    target->luminosity = source->luminosity;
    
    target->core_mass = source->core_mass;
    target->core_radius = source->core_radius;
    target->convective_envelope_mass = source->convective_envelope_mass;
    target->convective_envelope_radius = source->convective_envelope_radius;

    double z = source->metallicity;

    double *zpars;
    zpars = new double[20];
    zcnsts_(&z,zpars);
    target->zpars = zpars;

    for (int i=0; i<3; i++)
    {
        target->spin_vec[i] = source->spin_vec[i];
    }
    
    target->check_for_RLOF_at_pericentre = source->check_for_RLOF_at_pericentre;
}

void create_nested_system(ParticlesMap &particlesMap, int N_bodies, double *masses, int *stellar_types, double *smas, double *es, double *TAs, double *INCLs, double *APs, double *LANs)
{
    
    int N_binaries = N_bodies-1;

    int i;
    int index=0;
    for (i=0; i<N_bodies; i++)
    {
        
        Particle *p = new Particle(index, false);
        particlesMap[index] = p;
        p->mass = masses[i];
        p->stellar_type = stellar_types[i];
        p->sse_initial_mass = masses[i];
        p->evolve_as_star = true;
        
        index++;
    }

    int previous_binary;
        
    for (i=0; i<N_binaries; i++)
    {
        Particle *p = new Particle(index, true);
        (particlesMap)[index] = p;
        
        if (i==0)
        {
            p->child1 = particlesMap[0]->index;
            p->child2 = particlesMap[1]->index;
        }
        else
        {
            p->child1 = particlesMap[previous_binary]->index;
            p->child2 = particlesMap[i+1]->index;
        }
        //printf("%d %d %d %d %d\n",i,index,previous_binary,p->child1,p->child2);
        previous_binary = index;
        index++;
    }
    
    /* determine masses in all binaries */
    int N_root_finding,N_ODE_equations;
    //printf("ok\n");
    /* determine masses in all binaries */

    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);
    
    set_binary_masses_from_body_masses(&particlesMap);
    
    index=N_bodies;
    for (i=0; i<N_binaries; i++)
    {
        Particle *p = particlesMap[index];
        p->true_anomaly = TAs[i];
        //printf("DD %g %g  %g %g %g %g %g %g\n",p->child1_mass, p->child2_mass, smas[i], es[i], \
            INCLs[i], APs[i], LANs[i]);
        compute_orbital_vectors_from_orbital_elements(p->child1_mass, p->child2_mass, smas[i], es[i], \
            INCLs[i], APs[i], LANs[i], \
            &(p->e_vec[0]), &(p->e_vec[1]), &(p->e_vec[2]), &(p->h_vec[0]), &(p->h_vec[1]), &(p->h_vec[2]) );
        
        
        index++;
    }
    //initialize_stars(&particlesMap);
}

void print_system(ParticlesMap *particlesMap, int integration_flag)
{
    printf("=============================\n");
    printf("Printing system; N=%d; integration_flag=%d size = %d\n",particlesMap->size(),integration_flag,particlesMap->size());
    ParticlesMapIterator it_p;    
    if (integration_flag==0)
    {
        update_structure(particlesMap, integration_flag);
        set_up_derived_quantities(particlesMap);
    
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *p = (*it_p).second;
            
            if (p->is_binary == false)
            {
                if (p->evolve_as_star == true)
                {
                    printf("index %d -- body -- parent %d m %.15f r %g st %d mc %g minit %g menv %g epoch %g age %g rc %g renv %g lum %g Spin_freq %g Omega_crit %g\n",p->index,p->parent,p->mass,p->radius,p->stellar_type,p->core_mass,p->sse_initial_mass,p->convective_envelope_mass,p->epoch,p->age,p->core_radius,p->convective_envelope_radius,p->luminosity,norm3(p->spin_vec),compute_breakup_angular_frequency(p->mass,p->radius));
                }
                else
                {
                    printf("index %d -- body -- parent %d m %g r %g st %d \n",p->index,p->parent,p->mass,p->radius,p->stellar_type);
                }
            }
            else
            {
                printf("index %d -- binary -- parent %d child1 %d child2 %d m %g a %g e %g rp %g \n",p->index,p->parent,p->child1,p->child2,p->mass,p->a,p->e,p->a*(1.0 - p->e));
            }
        }
    }
    else
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *p = (*it_p).second;
            
            if (p->is_binary == false)
            {
                if (p->evolve_as_star == true)
                {
                    printf("index %d -- body -- parent %d m %.15f r %g st %d mc %g minit %g menv %g epoch %g age %g rc %g renv %g lum %g R_vec %g %g %g V_vec %g %g %g Spin_freq %g\n",p->index,p->parent,p->mass,p->radius,p->stellar_type,p->core_mass,p->sse_initial_mass,p->convective_envelope_mass,p->epoch,p->age,p->core_radius,p->convective_envelope_radius,p->luminosity,p->R_vec[0],p->R_vec[1],p->R_vec[2],p->V_vec[0],p->V_vec[1],p->V_vec[2],norm3(p->spin_vec));
                }
                else
                {
                    printf("index %d -- body -- parent %d m %g r %g st %d \n",p->index,p->parent,p->mass,p->radius,p->stellar_type);
                }
            }
        }
    }
}

#ifdef IGNORE
void print_bodies(ParticlesMap *particlesMap)
{
    printf("=============================\n");
    printf("Printing bodies; N=%d\n",particlesMap->size());
    
    //update_structure(particlesMap);
    //set_up_derived_quantities(particlesMap);

    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false)
        {
            if (p->evolve_as_star == true)
            {
                printf("index %d -- body -- parent %d m %.15f r %g st %d mc %g minit %g menv %g epoch %g age %g rc %g renv %g lum %g R_vec %g %g %g V_vec %g %g %g\n",p->index,p->parent,p->mass,p->radius,p->stellar_type,p->core_mass,p->sse_initial_mass,p->convective_envelope_mass,p->epoch,p->age,p->core_radius,p->convective_envelope_radius,p->luminosity,p->R_vec[0],p->R_vec[1],p->R_vec[2],p->V_vec[0],p->V_vec[1],p->V_vec[2]);
            }
            else
            {
                printf("index %d -- body -- parent %d m %g r %g st %d \n",p->index,p->parent,p->mass,p->radius,p->stellar_type);
            }
        }
    }
}
#endif

void get_parallel_and_perpendicular_vectors_and_components(double a_vec[3], double h_vec[3], double *a_par, double *a_perp, double a_perp_vec[3])
{
    *a_par = dot3(a_vec,h_vec);
    for (int i=0; i<3; i++)
    {
        a_perp_vec[i] = a_vec[i] - (*a_par) * h_vec[i];
    }
    *a_perp = norm3(a_perp_vec);
}

double get_mutual_angle(double a_vec[3], double b_vec[3])
{
    double a_norm = norm3(a_vec);
    double b_norm = norm3(b_vec);
    
    if (a_norm <= epsilon or b_norm <= epsilon)
    {
        return 0.0;
    }
    else
    {
        return acos( dot3(a_vec,b_vec) / ( a_norm * b_norm ) );
    }
}

double compute_a_from_h(double m1, double m2, double h, double e)
{
    return h * h * (m1+m2) / ( CONST_G * m1 * m1 * m2 * m2 * (1.0 - e*e ) );
}

double compute_h_from_a(double m1, double m2, double a, double e)
{
    return m1 * m2 * sqrt(CONST_G * a * (1.0 - e*e) / (m1 + m2));
}

double compute_h_dot_div_h(double m1, double m1_dot, double m2, double m2_dot, double a, double a_dot, double e, double e_dot)
{
    return m1_dot/m1 + m2_dot/m2 - c_1div2*(m1_dot + m2_dot)/(m1+m2) + c_1div2*(a_dot/a) - e*e_dot/(1.0 - e*e);
    //M_d_dot_av/M_d + M_a_dot_av/M_a - c_1div2*(M_d_dot_av + M_a_dot_av)/M + c_1div2*(da_dt/a) - e*de_dt/(1.0 - e*e);
}



bool equal_number(double x1, double x2, double tol)
{
    bool equal = false;
    if ( fabs((x1-x2)/x1) < tol)
    {
        equal = true;
    }
    return equal;
}

void check_number(double x, char *source, char *description, bool exit_on_error)
{
    double y = ((double) x);
    bool error = false;
    if (isnan(y))
    {
        error = true;
        printf("%s -- ERROR: quantity %s is NaN\n",source,description);
    }
    if (isinf(y))
    {
        error = true;
        printf("%s -- ERROR: quantity %s is infinite\n",source,description);
    }

    if (error == true and exit_on_error == true)
    {
        printf("Exiting on fatal error\n");
        exit(-1);
    }
}

int clear_particles(ParticlesMap *particlesMap)
{
    //printf("clear_internal_particles\n");
    (*particlesMap).clear();
	//highest_particle_index = 0;
    return 0;
}

void MSTAR_compute_center_of_mass_position_and_velocity(struct RegularizedRegion *R, double R_cm[3], double V_cm[3])
{
    int i,j;
    
    for (j=0; j<3; j++)
    {
        R_cm[j] = 0.0;
        V_cm[j] = 0.0;
    }

    double m;
    double m_tot=0.0;
    for (i=0; i<R->NumVertex; i++)
    {
        m = R->Mass[i];
        m_tot += m;
        for (j=0; j<3; j++)
        {
            R_cm[j] += m * R->Pos[3 * i + j];
            V_cm[j] += m * R->Vel[3 * i + j];
        }
    }

    for (j=0; j<3; j++)
    {
        R_cm[j] = R_cm[j]/m_tot;
        V_cm[j] = V_cm[j]/m_tot;
    }
}

void compute_center_of_mass_position_and_velocity(ParticlesMap *particlesMap, double R_cm[3], double V_cm[3])
{
    int i;
    double m;
    double M_tot = 0.0;
    for (i=0; i<3; i++)
    {
        R_cm[i] = 0.0;
        V_cm[i] = 0.0;
    }
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        
        if (p->is_binary == false)
        {
            m = p->mass;
            M_tot += m;
            for (i=0; i<3; i++)
            {
                R_cm[i] += m * p->R_vec[i];
                V_cm[i] += m * p->V_vec[i];
            }
        }
    }
    for (i=0; i<3; i++)
    {
        R_cm[i] /= M_tot;
        V_cm[i] /= M_tot;
    }
    printf("tools.cpp -- compute_center_of_mass_position_and_velocity -- R_CM %g %g %g V_CM %g %g %g\n",R_cm[0],R_cm[1],R_cm[2],V_cm[0],V_cm[1],V_cm[2]);
}

double find_nearest_neighbor_separation(ParticlesMap *particlesMap, int primary_index, double primary_R_vec[3])
{
    int i;
    double sep;
    double min_sep = 1.0e100;

    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        
        if (p->is_binary == false and p->index != primary_index)
        {
            sep = separation_between_vectors(primary_R_vec,p->R_vec);
            if (sep < min_sep)
            {
                min_sep = sep;
            }
        }
    }

    return min_sep;
}

void remove_binaries_from_system(ParticlesMap *particlesMap)
{
    std::vector<int> binary_indices;
    
    int i;

    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        if (p->is_binary == true)
        {
            binary_indices.push_back(p->index);
        }
        else
        {
            p->parent = -1;
            p->is_bound = false;
        }
    }
    
    for (int i=0; i<binary_indices.size(); i++)
    {
        particlesMap->erase(binary_indices[i]);
    }
}

}

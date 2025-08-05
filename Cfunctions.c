 // R CMD SHLIB -lm -lgsl -lgslcblas Cfunctions.c
 // ligne de compilation, compiler en term et après tester sur R

//--//--//--//--// Include libraries //--//--//--//--//

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//--//--//--//--// Define constants //--//--//--//--//

#include "Cconstants.c"

//--//--//--//--// Define structures //--//--//--//--//

// contains the parameters of the model
struct parameters {
  double beta;                // per capita birth rate
  double gamma;               // strength of density dependent mortality
  double sigma;               // per capita rate of cladogenesis
  double alpha;               // magnitude of the effect of gene flow on speciation
  double tau;                 // min delay for speciation
  double mui;                 // per capita migration rate between islands
  double ci;                  // mean dispersal distance for individuals between islands
  double muc;                 // per capita migration rate from the continent
  double cc;                  // mean dispersal distance for individuals from the continent
  long unsigned int Nc;       // number of individuals on the continent
  double theta;               // fundamental biodiversity number
  unsigned int nb_islands;    // number of islands in the archipelago
  double Kmax[NB_ISLAND_MAX]; // maximum carrying capacity for each island
  double ti[NB_ISLAND_MAX];   // time when each island emerges
  double pmax[NB_ISLAND_MAX]; // fraction of island life time at which its size is max
  double tm[NB_ISLAND_MAX];   // time at which island size is max
  double tl[NB_ISLAND_MAX];   // life time of each island
  double ts[NB_ISLAND_MAX];   // time of submergence of each island
  double dc[NB_ISLAND_MAX];   // initial distance of each island from the continent
  double **di;                // initial distances between islands
  double phi_max;           // maximum angle for apex angle of the island cone 
  double d;                   // increase coefficient of the apex angle of the island cone 
  double amp;                 // amplitude of sea level variations 
  double period;              // period of sea level variations 
  double phi_min;           // minimum apex angle
  double psi;               // mainland slope 
  double h_breakpoint;        // mainland height
};

// contains the variables of the continent
struct continent {
  double cumfreq[NB_SP_CONTI_MAX];    // cumulative frequencies of species of the continent
  long unsigned int nb_species;       // nb of different species on the continent
  double mucNc;                       // product muc * Nc
};

// contains the variables of the archipelago
struct archipelago {
  double K[NB_ISLAND_MAX];     // carrying capacity of each island
  double K_adj[NB_ISLAND_MAX];  // carrying capacity of each island adjusted with sea level 
  double phi[NB_ISLAND_MAX]; // apex angle of each island
  double h[NB_ISLAND_MAX];     // height of each island 
  double hz[NB_ISLAND_MAX];    // height of each island corrected with sea level
  double r[NB_ISLAND_MAX];     // theoretical radius of each island 
  double rz[NB_ISLAND_MAX];    // radius of each island corrected with sea level 
  double sl;                   // sea level 
#if ISL_ARCH == 2              // virtual carrying capacity of each island in the island-  
  double Kangle[NB_ISLAND_MAX];// archipelago 2 model (where dispersal angle is the same
  double kangle_adj[NB_ISLAND_MAX];
#endif                         // as in an archipelago
  double dc[NB_ISLAND_MAX];    // distance to the continent of each island
  double rc;                   // distance to the edge of the continent with sea level
  double **di;                 // distance between each island
  unsigned int subext[NB_ISLAND_MAX]; // equals 1 if the island is submerged and remaining individuals have been killed; equals 0 otherwise
  double local_mui[NB_ISLAND_MAX]; // per-capita local rate of inter-isl migr from each island
  double local_muc[NB_ISLAND_MAX]; // per-capita local rate of conti migration to each island
};

// contains the variables of the populations
struct population {
  unsigned int island_id;       // id of the island
  long unsigned int species_id; // id of the species
  double speciation_clock;      // waiting time to speciation ; >tmax means not a variant
  long unsigned int size;       // size of the pop
  double delta;                 // death rate of the pop
  unsigned int origin;          // 1=mig conti; 2=mig island; 3=clado
  double date_origin;           // date (time) when the pop appears
  double divergence_duration;   // duration of divergence (=tau+all delays due to gene flow)
};                              // (interpretable only for pop having speciated, ie speclock>tmax)

// contains the variables of the metapopulation
// CAUTION: metapopulation rates are not update at each event, they are computed from scratch.
// This is less efficient for computing time (although the computation cost is negligible:
// it increases by a factor 1/500), but it is required because truncating errors
// accumulate so that the sum of local rates becomes different from the total rate
struct meta_population {
  struct population pop[NB_POP_MAX]; // table of populations
  unsigned int nb_pop;      // total number of populations
  long unsigned int local_size[NB_ISLAND_MAX]; // pop size on each island
  long unsigned int total_size;                // total number of individuals
  double beta;              // total birth rate
  double delta;             // total death rate
  double sigma;             // total cladogenesis rate
  double mui;               // total inter-islands migration rate
  double muc;               // total rate of migration from the continent
};

// contains the variables of the phylogeny
struct phylogeny {
   long unsigned int *species_id; // id of the species
   double *date_origin;           // emergence time
   long unsigned int *mother;     // mother species
   double *birth;                 // birth time
   double *death;                 // death time
   long unsigned int nsp;      // number of species currently recorded in the phylogeny
};

// contains the counts useful to analyze the model but useless to run it;
// all members are count of events since last saved state; they are recorded on disk
// every wt_ss time steps; corresponding rates are n/wt_ss
struct outputs {
  long unsigned int n_mig_conti[NB_ISLAND_MAX];       // migration event from continent
  long unsigned int **n_mig_isl;                      // migration ev from island i to island j
  long unsigned int n_ana_conti_start[NB_ISLAND_MAX]; // anagenesis start after conti mig ev
  long unsigned int n_ana_isl_start[NB_ISLAND_MAX];   // anagenesis start after isl mig ev
  long unsigned int n_clado_start[NB_ISLAND_MAX];     // cladogenesis start ev
  long unsigned int n_spe_event[NB_ISLAND_MAX];       // speciation event on each island
  long unsigned int n_spe_event_ori_mig_conti[NB_ISLAND_MAX]; // spe ev of sp of origin 1
  long unsigned int n_spe_event_ori_mig_isl[NB_ISLAND_MAX];   // spe ev of sp of origin 2
  long unsigned int n_spe_event_ori_clado[NB_ISLAND_MAX];     // spe ev of sp of origin 3
  long unsigned int n_pop_ext[NB_ISLAND_MAX];         // pop extinction event
  long unsigned int n_spe_ext_isl[NB_ISLAND_MAX];     // species extinction event on island
  long unsigned int n_spe_ext_arch;                   // species extinction event on archip
};

// duration expressed in
struct duration {
  unsigned int d; // days
  unsigned int h; // hours
  unsigned int m; // minutes
  unsigned int s; // seconds
};


//--//--//--//--// functions //--//--//--//--//

// THIS FUNCTION ALLOCATES MEMORY FOR A MATRIX OF DOUBLE
void alloc_matrix_double(double ***mat, const unsigned int n, const unsigned int m) {

  *mat = malloc(n * sizeof (double **));
  if(*mat == NULL) {
    Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed\n\n");
    exit(1);
  }

  register long unsigned int i;
  for(i = 0; i < n; i++) {
    (*mat)[i] = malloc(m * sizeof (double *));
    if((*mat)[i] == NULL) {
      Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed\n\n");
      exit(1);
    }
  }
  
}

// THIS FUNCTION ALLOCATES MEMORY FOR A MATRIX OF LONG UNSIGNED INT
void alloc_matrix_luint(long unsigned int ***mat, const unsigned int n, const unsigned int m) {

  *mat = malloc(n * sizeof (long unsigned int **));
  if(*mat == NULL) {
    Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed\n\n");
    exit(1);
  }

  register long unsigned int i;
  for(i = 0; i < n; i++) {
    (*mat)[i] = malloc(m * sizeof (long unsigned int *));
    if((*mat)[i] == NULL) {
      Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed\n\n");
      exit(1);
    }
  }
  
}

// THIS FUNCTION RETURNS THE DIFFERENCE s1 - s2 (IN SECONDS) IN
// DAYS, HOURS, MINUTES AND SECONDS
struct duration time_exec(unsigned int s1, unsigned int s2) {

  struct duration t_exec;

  double duration = (double) (s2 - s1); 
  t_exec.d = (unsigned int) floor(duration / 86400.0);    // nb of days
  duration -= t_exec.d * 86400.0;
  t_exec.h = (unsigned int) floor(duration / 3600.0);     // nb of hours
  duration -= t_exec.h * 3600.0;
  t_exec.m = (unsigned int) floor(duration / 60.0);       // nb of minutes 
  t_exec.s = (unsigned int) (duration - t_exec.m * 60.0); // nb of seconds

  return t_exec;

}

// THIS FUNCTION CLEARS FILE file_name
void clear_file(char file_name[]) {

  FILE *f = NULL;
  f = fopen(file_name, "w");
  fclose(f);

}

// THIS FUNCTION MOVES FILES THAT HAVE TO BE SAVED TO THE RESULT FOLDER
void mv_files(char sim_name[], const unsigned int i_run) {

  // move output to the results folder
  FILE *term;
  term = popen("/bin/bash", "w");
  fprintf(term, "for i in $( ls out_* )\n");
  fprintf(term, "do\n");
  fprintf(term, "mv $i ${i/out_/%s.%d_}\n", sim_name, i_run);
  fprintf(term, "done\n");
  fprintf(term, "mv %s* results/\n", sim_name);
  fflush(term);
  pclose(term);

}

// THIS FUNCTION RETURNS THE MAX OF 2 DOUBLE
double min_double(double a, double b) {
  if (a < b)
    return (a);
  else
    return (b);  
}

// THIS FUNCTION RETURNS THE INDEX OF THE MIN VALUE AMONG 3 DOUBLE
unsigned int index_min_of_three(double a, double b, double c) {

  double min_val = a;
  unsigned int min_ind = 1;
  if (b < min_val) {
    min_val = b;
    min_ind = 2;
  }
  if (c < min_val) {
    min_ind = 3;
  }
  
  return(min_ind);
  
}
		  
// THIS FUNCTION INITIALIZES THE STATE OF THE CONTINENT; ADAPTATION OF THE FUNCTION rand.neutral
// FROM THE R PACKAGE untb. 
void initialize_continent(gsl_rng *rgsl, struct continent *conti, struct parameters *par) {

  long unsigned int count[NB_SP_CONTI_MAX]; // nb of indiv of each species id (initially at 0)
  long unsigned int i;
  for (i=0; i<NB_SP_CONTI_MAX; i++)
    count[i] = 0;

  long unsigned int spid = 0; // species id of the currently generated indiv
  (count[spid])++;            // 1st indiv is of species id 0
  for (i=1; i < (*par).Nc; i++) {
    // if the next indiv is of a new species id
    if (gsl_ran_flat(rgsl, 0.0, 1.0) < (*par).theta/((*par).theta + (double) i - 1.0) ) {
      spid++;
      if (spid > NB_SP_CONTI_MAX) {
	Rprintf("\n-------\n!ERROR!\n-------\nUnable to build the continent, NB_SP_CONTI_MAX too low\n\n");
	exit(1);
      }
      else
	(count[spid])++;
    }
    // if the next indiv is of a species id already existing
    else  
      if (i==1)
	(count[0])++;
      else {
	// k = index of the indiv similar to the current indiv
	long unsigned int k = (long unsigned int) floor(gsl_ran_flat(rgsl, 0.0, (double) i - 1.0));
	// ks = species id of this indiv k
	long unsigned int ks = 0;
	long unsigned int sum = count[ks];
	while (k >= sum)
	  sum += count[++ks];
	(count[ks])++;
      }
  }
  
  // convert the count as cumulative frequency
  (*conti).cumfreq[0] = count[0] / (double) (*par).Nc;
  for (i=1; i<spid; i++) 
    (*conti).cumfreq[i] = (*conti).cumfreq[i-1] + count[i] / (double) (*par).Nc;

  // record the total nb of species on the continent
  (*conti).nb_species = spid;

  // record the product muc * Nc
  (*conti).mucNc = (*par).muc * (*par).Nc;
  
}

// THIS FUNCTION CALCULATES THE CURRENT VALUE OF THE SEA LEVEL 
void sea_level_variations(struct parameters *par, struct archipelago *arch, const double t){

  (*arch).sl = (*par).amp * sin((2 * PI * t) / (*par).period);
  // Rprintf("sea level  = t: %g, amplitude: %lf, period: %lf, sl:%lf \n", t, par->amp, par->period, arch->sl);
  
  
}

// THIS FUNCTION RETURNS THE VALUE OF K FOR AN EMERGED ISLAND
double carrying_capacity(struct parameters *par, const double t, const unsigned int k_isl) {

  // the function must be called ONLY with a time corresponding to
  // a time where the island is emerged

  double f = (*par).pmax[k_isl] / (1.0-(*par).pmax[k_isl]);
  double d = ((*par).tm[k_isl] - (*par).ti[k_isl])/ ((1.0 + f)*0.1*(*par).tl[k_isl]);
  double c = f * d;
    
  return( (*par).Kmax[k_isl] * pow((t-(*par).ti[k_isl])/(*par).tl[k_isl], c) * pow(1.0 - (t-(*par).ti[k_isl])/(*par).tl[k_isl], d) / ( pow(c/(c+d), c) * pow(d/(c+d), d) ) );
  
}

// THIS FUNCTION RETURNS THE VALULE OF K FOR AN EMERGED ISLAND ADJUSTED WITH SEA LEVEL VARIATIONS
double adjusted_carrying_capacity(struct parameters *par, struct archipelago *arch, const double t, const unsigned k_isl) {

    double r = 0.0;
    double h = 0.0;    
    
    // cone parameters according to theoretical K
    (*arch).phi[k_isl] = (*par).phi_min + (*par).phi_max * (1 - exp(-(*par).d * (t - (*par).ti[k_isl]))); // apex angle based on island's age
    r = sqrt(((*arch).K[k_isl] * sin((*arch).phi[k_isl])) / PI);                                          // island's current radius based on island's apex angle and theoretical carrying capacity
    
    // normalize r and deduce h
    (*arch).r[k_isl] = rad_coef * r;
    (*arch).h[k_isl] = (*arch).r[k_isl] / tan((*arch).phi[k_isl]);
    
    // correct h with current sea level
    if ((*arch).sl >= (*arch).h[k_isl]) { // if sea level is above island's height
        (*arch).hz[k_isl] = 0;
    } else {  // if sea level is under island's height              
        (*arch).hz[k_isl] = (*arch).h[k_isl] - (*arch).sl;                             
    }
    
    // correct radius with sea level
    (*arch).rz[k_isl] = (*arch).hz[k_isl] * tan((*arch).phi[k_isl]);                 
  
    // denormalize r corrected with sea level and deduce h
    r =  (*arch).rz[k_isl] / rad_coef ;
    h =  r / tan((*arch).phi[k_isl]);
    
    return ( sqrt(r*r + h*h) * PI * r ); // Corrected K with sea level
}

// THIS FUNCTION GIVES THE DISPESAL KERNEL
double disp_kernel(const double c, const double d) {
  
  return( c / (TPI * (c*c+d*d) * sqrt(c*c+d*d)) );

}

// THIS FUNCTION INITIALIZES THE STATE OF THE ARCHIPELAGO
void initialize_archipelago(struct parameters *par, struct archipelago *arch, const double t) {
    unsigned int i, j;
  
    // No island has been submerged and their individuals killed yet
    for (i = 0; i < (*par).nb_islands; i++)
        (*arch).subext[i] = 0;

    // Set sea level
    sea_level_variations(par, arch, t); // Update sea level
    
    // set mainland edge
    (*arch).rc = tan((*par).psi) * ((*par).h_breakpoint - (*arch).sl) ;

    // Compute the initial value of K for each island
#if ISL_ARCH // Null model: a single island with size that of the equivalent archipelago
    for (i = 0; i < (*par).nb_islands; i++) { // All virtual islands at 0 but the first one which is filled 
        (*arch).K[i] = 0.0; 
        (*arch).K_adj[i] = 0.0;
        (*arch).phi[i] = 0.0;
        (*arch).h[i] = 0.0;
        (*arch).hz[i] = 0.0;
        (*arch).r[i] = 0.0;
        (*arch).rz[i] = 0.0;
    } // End for each island set to 0

    for (i = 0; i < (*par).nb_islands; i++) { // For each island
        if (t > (*par).ti[i] && t < (*par).ts[i]) { // Virtual island is emerged
            (*arch).K[0] += carrying_capacity(par, t, i); // Sum all emerged island's carrying capacities
        } // End virtual island is emerged
    } // End for each island

    (*arch).K_adj[0] = adjusted_carrying_capacity(par, arch, t, 0); // Calculate adjusted carrying capacity with sea level variations
     
#if ISL_ARCH == 2 // For the island-archipelago 2 model
    for (i = 0; i < (*par).nb_islands; i++) { // For each island
        if (t < (*par).ti[i] || t > (*par).ts[i]) { // Island is not emerged yet
            (*arch).Kangle[i] = 0.0; // Fill Kangle as if it was K in the archipelago model
            (*arch).Kangle_adj[i] = 0.0;
            (*arch).phi[i] = 0.0;
            (*arch).h[i] = 0.0;
            (*arch).hz[i] = 0.0;       
            (*arch).r[i] = 0.0;        
            (*arch).rz[i] = 0.0; 
        } else { // Island is emerged
            (*arch).Kangle[i] = carrying_capacity(par, t, i); // Theoretical island's carrying capacity
            (*arch).Kangle_adj[i] = adjusted_carrying_capacity(par, t, i); // Calculate adjusted carrying capacity with sea level variations
        } // End island is emerged
    } // End for each island
#endif

#else // Normal archipelago model    
    for (i = 0; i < (*par).nb_islands; i++) { // For each island
        if (t < (*par).ti[i] || t > (*par).ts[i]) { // If the island has not emerged yet
            (*arch).K[i] = 0.0;      
            (*arch).K_adj[i] = 0.0;
            (*arch).phi[i] = 0.0;
            (*arch).h[i] = 0.0;
            (*arch).hz[i] = 0.0;       
            (*arch).r[i] = 0.0;        
            (*arch).rz[i] = 0.0;      
        } else { // If the island has emerged
            (*arch).K[i] = carrying_capacity(par, t, i); 
            (*arch).K_adj[i] = adjusted_carrying_capacity(par, arch, t, i); // Calculate adjusted carrying capacity with sea level variations
        } // End if the island is emerged
    } // End for each island
#endif
  
    // Corrected distances from the continent and between islands 
    for (i = 0; i < (*par).nb_islands; i++) { // For each island
        (*arch).dc[i] = (*par).dc[i] - (*arch).rz[i] - (*arch).rc; // Distance from mainland + difference between theoretical radius and corrected radius 
        //Rprintf("init arch = t: %g, dc[i]: %lf, isl rad with sl: %lf \n", t, arch->dc[i], arch->rz[i]);
        
        for (j = 0; j < (*par).nb_islands; j++) { // For other neighboring islands from island i
            (*arch).di[i][j] = (*par).di[i][j] - (*arch).rz[i] - (*arch).rz[j]; // Distance between each island + difference between radii of each island
        } // End for each neighboring island from island i
    } // End for each island
      
    // Per capita rate of continental migration to each island
    for (i = 0; i < (*par).nb_islands; i++) { // For each island
#if ISL_ARCH == 2 // Use Kangle_adj to compute dispersal angles
        (*arch).local_muc[i] = sqrt((*arch).Kangle_adj[i]) / (PI_TH * (*arch).dc[i]) * disp_kernel((*par).cc, (*arch).dc[i]);
#else // Use K_adj to compute dispersal angles
        (*arch).local_muc[i] = sqrt((*arch).K_adj[i]) / (PI_TH * (*arch).dc[i]) * disp_kernel((*par).cc, (*arch).dc[i]);
#endif
    } // End for each island

    // Per capita rate of inter-island migration
    for (i = 0; i < (*par).nb_islands; i++) {   // Island where come migrants from
        double sum = 0.0;
        for (j = 0; j < (*par).nb_islands; j++) { // Island where they go
            if (i != j) 
                sum += (sqrt((*arch).K_adj[j]) / (PI_TH * (*arch).di[i][j]) * disp_kernel((*par).ci, (*arch).di[i][j]));
        }
        (*arch).local_mui[i] = (*par).mui * sum;    
    }    
}

// THIS FUNCTION INITIALIZES THE STATE OF THE METAPOPULATION
void initialize_metapop(struct parameters *par, struct archipelago *arch, struct meta_population *metapop, struct continent *conti, const double t_max) {
  
  // initialize population properties.
  // it is just to initialize the memory, because there is initially nobody 
  // ((*metapop).nb_pop=0) which means that there is no line in this table
  unsigned int k;
  for (k = 0; k < NB_POP_MAX; k++) {
    (*metapop).pop[k].island_id = 0;
    (*metapop).pop[k].species_id = 0;
    (*metapop).pop[k].speciation_clock = 0.0;
    (*metapop).pop[k].size = 0;
    (*metapop).pop[k].delta = 0.0;
    (*metapop).pop[k].origin = 0;
    (*metapop).pop[k].date_origin = 0.0;
    (*metapop).pop[k].divergence_duration = 0.0;
  }

  // initialize metapopulation properties
  (*metapop).nb_pop = 0;
  (*metapop).total_size = 0;
  (*metapop).beta = 0.0;
  (*metapop).delta = 0.0;
  (*metapop).sigma = 0.0;
  (*metapop).mui = 0.0;
  for (k = 0; k < (*par).nb_islands; k++)
    (*metapop).local_size[k] = 0;

  // initialize migration from continent: the only rate which is not 0
  double sum = 0.0;   // probability to reach at least one of the islands
  for (k=0; k<(*par).nb_islands; k++)
    sum += (*arch).local_muc[k];  
  (*metapop).muc = (*conti).mucNc * sum; // total rate of migration from the continent

}

// THIS FUNCTION INITALIZES PHYLOGENY 
void initialize_phylo(struct phylogeny *phylo) {
	// initialize species properties
	unsigned int k;
	for (k = 0; k < NB_SP_PHYLO_MAX; k++) {
		(*phylo).species_id[k] = 0;
		(*phylo).date_origin[k] = 0;
		(*phylo).mother[k] = 0;
		(*phylo).birth[k] = 0;
		(*phylo).death[k] = 0;
	}
        (*phylo).nsp = 0;
}

// THIS FUNCTION INITIALIZES THE OUTPUTS (RATES) AT ZERO
void initialize_rates(struct parameters *par, struct outputs *rates) {

  unsigned int i, j;

  (*rates).n_spe_ext_arch = 0;
  for (i=0; i<(*par).nb_islands; i++) { 
    (*rates).n_mig_conti[i] = 0;
    (*rates).n_ana_isl_start[i] = 0;
    (*rates).n_ana_conti_start[i] = 0;
    (*rates).n_clado_start[i] = 0;
    (*rates).n_spe_event[i] = 0;
    (*rates).n_spe_event_ori_mig_conti[i] = 0;
    (*rates).n_spe_event_ori_mig_isl[i] = 0;
    (*rates).n_spe_event_ori_clado[i] = 0;	
    (*rates).n_pop_ext[i] = 0;
    (*rates).n_spe_ext_isl[i] = 0;
    for (j=0; j<(*par).nb_islands; j++)
      (*rates).n_mig_isl[i][j] = 0;
  }

}

// THIS FUNCTION PRINTS THE STATE OF THE CONTINENT
void print_conti(struct continent *conti) {

  FILE *conti_file = NULL;
  conti_file = fopen("out_conti", "w");

  long unsigned int i;
  for (i=0; i<(*conti).nb_species; i++)
    fprintf(conti_file, "%lu\t%g\n", i, (*conti).cumfreq[i]);
  fclose(conti_file);

}

// THIS FUNCTION PRINTS THE STATE OF THE ARCHIPELAGO
void print_arch(struct parameters *par, struct archipelago *arch, const double t) {

  FILE *arch_file = NULL;
  arch_file = fopen("out_arch", "a");

  fprintf(arch_file, "%g\t", t);

  unsigned int i;
  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).K[i]);

  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).dc[i]);

  unsigned int j; // if nb_islands=1, no di written
  
  for (i=0; i<(*par).nb_islands; i++){
    for (j=0; j<(*par).nb_islands; j++){
      fprintf(arch_file, "%g\t", (*arch).di[i][j]);
      //Rprintf("PRINT: arch.di[%d][%d]: %g\n", i, j, (*arch).di[i][j]); 
    }
  }
    
  //for (i=0; i<((*par).nb_islands-1); i++) {
    //for (j=(i+1); j<(*par).nb_islands; j++) {
      //fprintf(arch_file, "%g\t", (*arch).di[i][j]);
      //Rprintf("PRINT: t: %g di_obs: %f for isl i: %d and isl j: %d\n", t, (*arch).di[i][j], i, j);
      //Rprintf("PRINT: arch.di[%d][%d]: %g\n", i, j, (*arch).di[i][j]); 
    //}
  //}
      
  fprintf(arch_file, "%g\t", (*arch).sl); // save sea level
  
  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).K_adj[i]);
    
  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).phi[i]);

  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).h[i]);
  
  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).hz[i]); // save each island's height 
    
  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).r[i]); // save each island's theoretical radius 
    
  for (i=0; i<(*par).nb_islands; i++) 
    fprintf(arch_file, "%g\t", (*arch).rz[i]); // save each island's adjusted radius 
    
  fprintf(arch_file, "%g\t", (*arch).rc); // save mainland egde

  fprintf(arch_file, "\n");
  
  fclose(arch_file);

}

// THIS FUNCTION PRINTS THE STATE OF THE METAPOPULATION
void print_metapop(struct meta_population *metapop, const double t) {

  FILE *metapop_file = NULL;
  metapop_file = fopen("out_metapop", "a");

  if ((*metapop).nb_pop != 0) {

    long unsigned int i;
    for (i=0; i<(*metapop).nb_pop; i++) 
      fprintf(metapop_file, "%g\t%u\t%lu\t%g\t%lu\t%u\t%g\t%g\n", t, (*metapop).pop[i].island_id, (*metapop).pop[i].species_id, (*metapop).pop[i].speciation_clock, (*metapop).pop[i].size, (*metapop).pop[i].origin, (*metapop).pop[i].date_origin, (*metapop).pop[i].divergence_duration);
  }
  else {
    fprintf(metapop_file, "%g\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", t);
  }
  fclose(metapop_file);

}

// THIS FUNCTION PRINTS THE PHYLOGENY
void print_phylo(struct phylogeny *phylo) {

  FILE *phylo_file = NULL ;
  phylo_file = fopen("out_phylo","a");

  long unsigned int i;
  for (i=0; i<(*phylo).nsp; i++)
  {
    fprintf(phylo_file, "%lu\t%g\t%lu\t%g\t%g\n", (*phylo).species_id[i],(*phylo).date_origin[i],(*phylo).mother[i],(*phylo).birth[i],(*phylo).death[i]);
  }
  //printf("Print species_id dans le fichier");
  fclose(phylo_file);
}

// THIS FUNCTION PRINTS THE RATES USEFUL TO ANALYSE THE MODEL
void print_rates(struct parameters *par, struct outputs *rates, const double t) {

  unsigned int i, j;

  // save the rates
  FILE *rates_file = NULL;
  rates_file = fopen("out_rates", "a");

  fprintf(rates_file, "%g\t", t);
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_mig_conti[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    for (j=0; j<(*par).nb_islands; j++)
      fprintf(rates_file, "%lu\t", (*rates).n_mig_isl[i][j]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_ana_isl_start[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_ana_conti_start[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_clado_start[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_spe_event[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_spe_event_ori_mig_conti[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_spe_event_ori_mig_isl[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_spe_event_ori_clado[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_pop_ext[i]);  
  for (i=0; i<(*par).nb_islands; i++)
    fprintf(rates_file, "%lu\t", (*rates).n_spe_ext_isl[i]);  
  fprintf(rates_file, "%lu\n", (*rates).n_spe_ext_arch);  

  fclose(rates_file);

  // re-initialize the rates
  initialize_rates(par, rates);

}

// THIS FUNCTION PRINTS THE CURRENT TIME
void print_time(const double t) {
  
  Rprintf("Current time: %g...\n", t);
  FILE *log_file = NULL;
  log_file = fopen("out_log", "a");
  fprintf(log_file, "Current time: %g...\n", t);
  fclose(log_file);

}

// THIS FUNCTION ADJUST THE TABLE OF POPULATIONS AFTER AN EXTINCTION
void adjust_table(struct meta_population *metapop, const unsigned int k_ext, struct outputs *rates, struct phylogeny *phylo, const double t) {

  if ((*metapop).pop[k_ext].size == 0) {    

    // keep the island id, species id and origin code of the extinct pop
   unsigned int k_isl = (*metapop).pop[k_ext].island_id;
   long unsigned int k_spe = (*metapop).pop[k_ext].species_id;
  
/*   printf("k_ext %lu \n",(*metapop).pop[k_ext].species_id);*/
    // adjust the table
   ((*metapop).nb_pop)--; // index of the last pop in the table
    if (k_ext != (*metapop).nb_pop) {
      (*metapop).pop[k_ext].island_id = (*metapop).pop[(*metapop).nb_pop].island_id;
      (*metapop).pop[k_ext].species_id = (*metapop).pop[(*metapop).nb_pop].species_id;
      (*metapop).pop[k_ext].speciation_clock = (*metapop).pop[(*metapop).nb_pop].speciation_clock;
      (*metapop).pop[k_ext].size = (*metapop).pop[(*metapop).nb_pop].size;
      (*metapop).pop[k_ext].delta = (*metapop).pop[(*metapop).nb_pop].delta;
     (*metapop).pop[k_ext].origin = (*metapop).pop[(*metapop).nb_pop].origin;
      (*metapop).pop[k_ext].date_origin = (*metapop).pop[(*metapop).nb_pop].date_origin;
      (*metapop).pop[k_ext].divergence_duration = (*metapop).pop[(*metapop).nb_pop].divergence_duration;
    }

    // record that an event of population extinction occurred
    ((*rates).n_pop_ext[k_isl])++;

    // check (and record) if this population extinction corresponds to a species extinction
    // on the island, and maybe on the archipelago
    unsigned int still_on_island = 0; // =1 if species of extinct pop still exists on island
    unsigned int still_on_archi = 0;  // =1 if species of extinct pop still exists on archip
    long unsigned int i; 
    long int j;
    for (i=0; i<(*metapop).nb_pop; i++) {
      if ((*metapop).pop[i].species_id == k_spe) {
	still_on_archi = 1; 
	if ((*metapop).pop[i].island_id == k_isl) {
	  still_on_island = 1;
        }
      }
    }
    if (still_on_island == 0) {
      ((*rates).n_spe_ext_isl[k_isl])++;
      if (still_on_archi == 0){
	((*rates).n_spe_ext_arch)++; 
	// the species is no longer on the archipelago, 
	// Search the species_id corresponding in struct phylogeny
   
	for (j=(*phylo).nsp-1;j>=0;j--){
	  if ((*phylo).species_id[j]==k_spe){
	    (*phylo).death[j]=t;
	    /*                if (t>110000) {*/
	    /*                       printf("espèce num %lu\t death time %g\n", (*phylo).species_id[j],(*phylo).death[j]);*/
	    /*                }*/
	    break;
	  }
	}
      }
     
    }
    
  }

}

// THIS FUNCTION ACHIEVES A SPECIATION EVENT FOR POP K_POP
void speciation(const double t, struct meta_population *metapop, long unsigned int *new_sp_id, struct outputs *rates, const double t_max, const long unsigned int k_pop, struct phylogeny *phylo) {
  
  // keep the island id and species id of the speciating pop
  unsigned int k_isl = (*metapop).pop[k_pop].island_id;
  long unsigned int k_spe = (*metapop).pop[k_pop].species_id;

  // Save mother species id in struct phylogeny
  (*phylo).mother[(*phylo).nsp] = (*metapop).pop[k_pop].species_id;
  //if (t>110000) {
    //printf("Speciation event, mother species : %lu\n",(*phylo).mother[(*phylo).nsp]);
  //}

  // Define the properties of the new species
  (*metapop).pop[k_pop].species_id = *new_sp_id;
  (*metapop).pop[k_pop].speciation_clock = 2.0*t_max; // no speciation next, waiting time as infinite

  // Save the new species caracteristics in struct phylogeny
  if ((*phylo).nsp == NB_SP_PHYLO_MAX) {
  Rprintf("\n-------\n!ERROR!\n-------\nNumber of species too high, NB_SP_PHYLO_MAX too low\n\n");
  exit(1);
  }
   else { 
  (*phylo).species_id[(*phylo).nsp] = *new_sp_id; //species_id
  (*phylo).date_origin[(*phylo).nsp] = (*metapop).pop[k_pop].date_origin; //date_origin
  (*phylo).birth[(*phylo).nsp] = t; //birth

   //if (t>110000) {
     //printf("Daughter species num : %lu\t mother %lu\t date origin %g\t birth %g\n",(*phylo).species_id[(*phylo).nsp],(*phylo).mother[(*phylo).nsp],(*phylo).date_origin[(*phylo).nsp],(*phylo).birth[(*phylo).nsp]);
   //}

   }

  // extinction from speciation
  // check (and record) if this speciation corresponds to a species extinction
  // on the island, and maybe on the archipelago (it occurs when the speciating population
  // is the only population of the species)
  unsigned int still_on_island = 0; // =1 if species of extinct pop still exists on island
  unsigned int still_on_archi = 0;  // =1 if species of extinct pop still exists on archip
  long unsigned int i;
  long int j;  
  for (i=0; i<(*metapop).nb_pop; i++) {
    if ((*metapop).pop[i].species_id == k_spe) {
      still_on_archi = 1; //the species still exist, there is no extinction
      //if (t>114000) {
	       //printf("still on archi = 1 \n");
      //}
      if ((*metapop).pop[i].island_id == k_isl)
	still_on_island = 1;
    }
  }
  if (still_on_island == 0) {
    ((*rates).n_spe_ext_isl[k_isl])++;
    if (still_on_archi == 0) {
      //printf("still_on_archi == 0");
      ((*rates).n_spe_ext_arch)++;
      for (j=(*phylo).nsp-1;j>=0;j--) { // Search the species_id corresponding in struct phylogeny 
	if ((*phylo).species_id[j]==k_spe) {
	  //printf("(*phylo).species_id[j] = %lu , t = %f",(*phylo).species_id[j],t);
	  (*phylo).death[j]=t;
	  break;
	}
      }
    }
  }

  // update the metapopulation properties
  (*new_sp_id)++;

  // update the number of speciation events
  ((*phylo).nsp)++;

  // record that a speciation event occurred on island (*metapop).pop[k_pop].island_id];  
  ((*rates).n_spe_event[(*metapop).pop[k_pop].island_id])++;
  if ((*metapop).pop[k_pop].origin == 1)
    ((*rates).n_spe_event_ori_mig_conti[(*metapop).pop[k_pop].island_id])++;
  if ((*metapop).pop[k_pop].origin == 2)
    ((*rates).n_spe_event_ori_mig_isl[(*metapop).pop[k_pop].island_id])++;
  if ((*metapop).pop[k_pop].origin == 3)
    ((*rates).n_spe_event_ori_clado[(*metapop).pop[k_pop].island_id])++;

  }


// THIS FUNCTION UPDATES THE LANDSCAPES, THE RATES AND THE SPECIATION TIMES
void update_landscape(struct parameters *par, struct continent *conti, struct archipelago *arch, struct meta_population *metapop, struct outputs *rates, const double t, const double wt_up, const double t_max, long unsigned int *new_sp_id, struct phylogeny *phylo) {

  unsigned int i;

#if !STATIC // the landscape and the rates are updated only for dynamic archipelagos
  
  unsigned int j;
  double sum;
  
  sea_level_variations(par, arch, t); // update sea level
  
  (*arch).rc = tan((*par).psi)*((*par).h_breakpoint-(*arch).sl); // update mainland edge
  
#if ISL_ARCH // null model: a single island with size that of the equivalent archipelago
  // no need to check here that the island-archipel-equivalent (IAE) is initially emerged,
  // it's checked in the R code calling this C program.
  // no need either to kill all remaining individuals when the IAE is submerged because
  // the simulation stops at that time (ensured in the R code)
    
  (*arch).K[0] = 0.0;        // the IAE is the first island of the table; the K of other virtual islands is initially set at 0 and is never changed 
  (*arch).K_adj[0] = 0.0; 
  (*arch).phi[0] = 0.0; 
  (*arch).h[0] = 0.0; 
  (*arch).hz[0] = 0.0; 
  (*arch).r[0] = 0.0; 
  (*arch).rz[0] = 0.0; 
  
  for (i=0; i<(*par).nb_islands; i++){ // for each island
    if (t > (*par).ti[i] && t < (*par).ts[i]){ // virtual island is emerged
      (*arch).K[0] += carrying_capacity(par, t, i); // sum carrying capacity in first island
    }
  }  
  (*arch).K_adj[0] = adjusted_carrying_capacity(par, arch,  t, 0);
  if ((*arch).K_adj[0] <= 0) { // if the first island is emerged but under sea level 
    if ((*arch).subext[i] == 0) { // if remaining individuals have yet been killed, kill them
      // look for populations in this island
      unsigned int keeplooking = 1; // =1 if there is at least 1 pop on the submerged island
      while (keeplooking == 1) {
	keeplooking = 0;
	for (j=0; j<(*metapop).nb_pop; j++) {
	  if ((*metapop).pop[j].island_id == i) {
	    // kill individuals of this population, ie update sizes
	    // local and total rates are updated once all remaining pops have been killed
	    (*metapop).total_size -= (*metapop).pop[j].size;
	    (*metapop).local_size[i] -= (*metapop).pop[j].size;
            (*metapop).pop[j].size = 0;	      	     
	    // adjust the table
	    adjust_table(metapop, j, rates, phylo, t);
	    j = (*metapop).nb_pop + 1; // +1 because adjust_table() decreases nb_pop by 1    
	    // the table has changed, start again the search for pops to remove
	    keeplooking = 1;
	  }
	} 
      } 
      // update total rates which depend on pop size (beta, sigma)
      // the other rates depending on K will be update anyway after (delta, mui, muc)
      // (*metapop).muc does not depend on local size, it does not need to be updated
      (*metapop).beta = (*par).beta * (*metapop).total_size;
      (*metapop).sigma = (*par).sigma * (*metapop).total_size;
      // island submerged and all remaining indiv have been killed
      (*arch).subext[i] = 1;	
    } 
  } // end if island is emerged but under sea level
  (*arch).dc[i] = (*par).dc[i] - (*arch).rz[i] - (*arch).rc; // adjust dc = distance from mainland - corrected radius
      
#if ISL_ARCH == 2 // for the island-archipelago 2 model, update Kangle
  for (i=0; i<(*par).nb_islands; i++) { // for each island
    if (t < (*par).ti[i] || t > (*par).ts[i]) { // if the island is not emerged yet
      (*arch).Kangle[i] = 0.0;
      (*arch).Kangle_adj[i] = 0.0; 
      (*arch).phi[i] = 0.0; 
      (*arch).h[i] = 0.0; 
      (*arch).hz[i] = 0.0; 
      (*arch).r[i] = 0.0; 
      (*arch).rz[i] = 0.0; 
    } // end if the island is not emerged yet
    else { // if the island is emerge
      (*arch).Kangle[i] = carrying_capacity(par, t, i); 
      (*arch).Kangle_adj[i] = adjusted_carrying_capacity(par, arch, t, i);
    } // end if the island is emerged
  } // end for each island
#endif
  
#else // normal archipelago model
    
  for (i=0; i<(*par).nb_islands; i++) { // for each island
  
    if (t < (*par).ti[i]) {    // if the island is not emerged 
      (*arch).K[i] = 0.0; 
      (*arch).K_adj[i] = 0.0; 
      (*arch).phi[i] = 0.0; 
      (*arch).h[i] = 0.0; 
      (*arch).hz[i] = 0.0; 
      (*arch).r[i] = 0.0; 
      (*arch).rz[i] = 0.0; 
    } // end if island not emerged yet 
    
    else if (t >= (*par).ts[i]) { // if the island is submerged
     (*arch).K[i] = 0.0; 
     (*arch).K_adj[i] = 0.0; 
     (*arch).phi[i] = 0.0; 
     (*arch).h[i] = 0.0; 
     (*arch).hz[i] = 0.0; 
     (*arch).r[i] = 0.0; 
     (*arch).rz[i] = 0.0;
     
     // if remaining individuals have yet been killed, kill them
     if ((*arch).subext[i] == 0) {
       // look for populations in this island
       unsigned int keeplooking = 1; // =1 if there is at least 1 pop on the submerged island
       while (keeplooking == 1) {
         keeplooking = 0;
	 for (j=0; j<(*metapop).nb_pop; j++) {
	   if ((*metapop).pop[j].island_id == i) {
	      
	     // kill individuals of this population, ie update sizes
	     // local and total rates are updated once all remaining pops have been killed
	     (*metapop).total_size -= (*metapop).pop[j].size;
	     (*metapop).local_size[i] -= (*metapop).pop[j].size;
	     (*metapop).pop[j].size = 0;	      	     
	      
	      // adjust the table
	      adjust_table(metapop, j, rates, phylo, t);
	      j = (*metapop).nb_pop + 1; // +1 because adjust_table() decreases nb_pop by 1
	      
	      // the table has changed, start again the search for pops to remove
	      keeplooking = 1;
	   }
	 }
       }

       // update total rates which depend on pop size (beta, sigma)
       // the other rates depending on K will be update anyway after (delta, mui, muc)
       // (*metapop).muc does not depend on local size, it does not need to be updated
       (*metapop).beta = (*par).beta * (*metapop).total_size;
       (*metapop).sigma = (*par).sigma * (*metapop).total_size;

       // island submerged and all remaining indiv have been killed
       (*arch).subext[i] = 1;	
     }
   } // end if island is submerged
   
   else { // if the island is emerged
   
     (*arch).K[i] = carrying_capacity(par, t, i); // update carrying capacity
     (*arch).K_adj[i] = adjusted_carrying_capacity(par, arch, t, i); // update carrying capacity
      
     if ((*arch).K_adj[i] <= 0.0) { // if the island is emerged but under sea level 
       if ((*arch).subext[i] == 0) { // if remaining individuals have yet been killed, kill them
	 // look for populations in this island
	 unsigned int keeplooking = 1; // =1 if there is at least 1 pop on the submerged island
	 while (keeplooking == 1) {
	   keeplooking = 0;
	   for (j=0; j<(*metapop).nb_pop; j++) {
	     if ((*metapop).pop[j].island_id == i) {
	      
	       // kill individuals of this population, ie update sizes
	       // local and total rates are updated once all remaining pops have been killed
	       (*metapop).total_size -= (*metapop).pop[j].size;
	       (*metapop).local_size[i] -= (*metapop).pop[j].size;
	       (*metapop).pop[j].size = 0;	      	     
	      
	       // adjust the table
	       adjust_table(metapop, j, rates, phylo, t);
	       j = (*metapop).nb_pop + 1; // +1 because adjust_table() decreases nb_pop by 1
	      
	       // the table has changed, start again the search for pops to remove
	       keeplooking = 1;
	     }
	   } 
	 } 
	 // update total rates which depend on pop size (beta, sigma)
	 // the other rates depending on K will be update anyway after (delta, mui, muc)
	 // (*metapop).muc does not depend on local size, it does not need to be updated
	 (*metapop).beta = (*par).beta * (*metapop).total_size;
	 (*metapop).sigma = (*par).sigma * (*metapop).total_size;
	 // island submerged and all remaining indiv have been killed
	 (*arch).subext[i] = 1;	
      } 
    } // end if island is emerged but under sea level
   } // end if the island is emerged 
  } // end for each island
  
  for (i=0; i<(*par).nb_islands; i++){ // for each island new distances
  (*arch).dc[i] = (*par).dc[i] - (*arch).rz[i] - (*arch).rc; // dc = distance from mainaland + difference between theoretical radius and corrected radius 
    for (j=0; j<(*par).nb_islands; j++) { // for neighbouring islands
      (*arch).di[i][j] = (*par).di[i][j] - (*arch).rz[i] - (*arch).rz[j]; // di = distance between each island + difference betwene th_rad and sl_rad of each islands 
          //Rprintf("TIME: %f\n d[i][j]: %f with i: %u and j: %u, should be: di_init: %f - r1: %f - r2: %f\n", t, (*arch).di[i][j], i, j, (*par).di[i][j], (*arch).rz[i], (*arch).rz[j]);
          //Rprintf("TIME: %f\n di_init: %f, di_observed: %f, rzi: %f, rzj: %f\n with i: %d and j: %d\n", t, (*par).di[i][j], (*arch).di[i][j], (*arch).rz[i], (*arch).rz[j], i, j);
          //Rprintf("TIME: %f  arch.di[%d][%d]: %g\n", t, i, j, (*arch).di[i][j]);
      if(t >=90 && t<= 110) {
        Rprintf("TIME: %f  i: %d j: %d arch.di: %g, rzi: %g, rzj: %g\n", t, i, j, (*arch).di[i][j], (*arch).rz[i], (*arch).rz[j]);
      }
          
    }
   }
  
#endif
 
  // Update the rates depending on K (delta, mui, muc)
  // beta and sigma depend on pop size, they have been update above if necessary
  // I don't update delta and mui  with the function update_rates() because all local rates
  // must be update (and because update_rates() update beta and sigma which is not needed here)

  // death rate
  (*metapop).delta = 0.0;
  for (i=0; i<(*metapop).nb_pop; i++) {
    (*metapop).pop[i].delta = (*metapop).pop[i].size * (*par).beta*exp(-(*par).gamma*(1.0 - (*metapop).local_size[(*metapop).pop[i].island_id]/(*arch).K_adj[(*metapop).pop[i].island_id]));
    (*metapop).delta += (*metapop).pop[i].delta;
  }  
 
  // rate of migration between islands
  (*metapop).mui = 0.0;
  for (i=0; i<(*par).nb_islands; i++) {   // island where come migrants from
    sum = 0.0;
    for (j=0; j<(*par).nb_islands; j++)   // island where they go
      if (i != j) 
	sum += (sqrt((*arch).K_adj[j])/(PI_TH*(*arch).di[i][j]) * disp_kernel((*par).ci,(*arch).di[i][j]));
    (*arch).local_mui[i] = (*par).mui * sum;    
    (*metapop).mui += ((*arch).local_mui[i] * (*metapop).local_size[i]);
  }

  // rate of migration from the continent
  sum = 0.0;
  for (i=0; i<(*par).nb_islands; i++) {
#if ISL_ARCH == 2 // use Kangle to compute dispersal angles    
    (*arch).local_muc[i] = sqrt((*arch).Kangle_adj[i])/(PI_TH*(*arch).dc[i]) * disp_kernel((*par).cc,(*arch).dc[i]);
#else
    (*arch).local_muc[i] = sqrt((*arch).K_adj[i])/(PI_TH*(*arch).dc[i]) * disp_kernel((*par).cc,(*arch).dc[i]);
#endif
    sum += (*arch).local_muc[i];
  }
  (*metapop).muc = (*conti).mucNc * sum;

#endif // end of if !STATIC; speciation events must be achieved even for static archipelagos
  
  // update speciation clocks and achieve speciation events
  for (i=0; i<(*metapop).nb_pop; i++) {
    (*metapop).pop[i].speciation_clock -= wt_up;
    if ((*metapop).pop[i].speciation_clock <= 0.0)
      speciation(t, metapop, new_sp_id, rates, t_max, i, phylo);
  }

}

// THIS FUNCTION UPDATES THE TOTAL RATES AFTER A LOCAL POPULATION SIZE CHANGE
// (UPDATED RATES: BETA, SIGMA, MUI AND DELTA)
void update_rates(struct parameters *par, struct archipelago *arch, struct meta_population *metapop, const unsigned int k_isl) {

  // (*metapop).muc does indep of local size, not affected by population size change
  
  // birth rate
  (*metapop).beta = (*par).beta * (*metapop).total_size;

  // cladogenesis rate
  (*metapop).sigma = (*par).sigma * (*metapop).total_size;

  // death rate (changed for all indiv on the island where pop size change occurred, k_isl)
  long unsigned int i;
  (*metapop).delta = 0.0;
  for (i=0; i<(*metapop).nb_pop; i++) {
    if ((*metapop).pop[i].island_id == k_isl)
      (*metapop).pop[i].delta = (*metapop).pop[i].size * (*par).beta*exp(-(*par).gamma*(1.0 - (*metapop).local_size[(*metapop).pop[i].island_id]/(*arch).K_adj[(*metapop).pop[i].island_id]));    
    (*metapop).delta += (*metapop).pop[i].delta;
  }

  // inter-island migration rate
  (*metapop).mui = 0.0;
  for (i=0; i<(*par).nb_islands; i++)
    (*metapop).mui += ((*arch).local_mui[i] * (*metapop).local_size[i]);


}

// THIS FUNCTION ACHIEVES A BIRTH
void birth(gsl_rng *rgsl, struct meta_population *metapop, struct parameters *par, struct archipelago *arch) {

  // choose who is giving birth
  double u = gsl_ran_flat(rgsl, 0.0, (*metapop).total_size);
  unsigned int k_pop = 0;
  long unsigned int sum = (*metapop).pop[k_pop].size;
  while (u >= sum)
    sum += (*metapop).pop[++(k_pop)].size;

  // update population and metapopulation properties
  ((*metapop).pop[k_pop].size)++;
  ((*metapop).local_size[(*metapop).pop[k_pop].island_id])++;
  ((*metapop).total_size)++;

  // update total rates
  update_rates(par, arch, metapop, (*metapop).pop[k_pop].island_id);

}

// THIS FUNCTION ACHIEVES A DEATH
void death(gsl_rng *rgsl, struct meta_population *metapop, struct outputs *rates, struct parameters *par, struct archipelago *arch, struct phylogeny *phylo, const double t) {

  // choose who is dying
  double u = gsl_ran_flat(rgsl, 0.0, (*metapop).delta);
  long unsigned int k_pop = 0;
  double sum = (*metapop).pop[k_pop].delta;
  while (u >= sum)
    sum += (*metapop).pop[++(k_pop)].delta;
    
  // update population and metapopulation properties
  ((*metapop).pop[k_pop].size)--;
  ((*metapop).local_size[(*metapop).pop[k_pop].island_id])--;
  ((*metapop).total_size)--;

  // update total rates
  update_rates(par, arch, metapop, (*metapop).pop[k_pop].island_id);

  // if there was only 1 indiv in the pop, the pop is extinct, replace it in the table
  // with the last pop of the table
  adjust_table(metapop, k_pop, rates, phylo, t);
  
}

// THIS FUNCTION ACHIEVES A CLADOGENESIS EVENT INITIATION
void cladogenesis_start(gsl_rng *rgsl, const double t, struct parameters *par, struct meta_population *metapop, struct outputs *rates, struct archipelago *arch, struct phylogeny *phylo) {

  // check that memory limits are not exceeded
  if ((*metapop).nb_pop == NB_POP_MAX) {
    Rprintf("\n-------\n!ERROR!\n-------\nNumber of populations too high, NB_POP_MAX too low\n\n");
      exit(1);
  }
  else { // a cladogenesis event is possible
 
    // choose who initating a new lineage is giving birth
    double u = gsl_ran_flat(rgsl, 0.0, (*metapop).total_size);
    unsigned int k_pop = 0;
    long unsigned int sum = (*metapop).pop[k_pop].size;
    while (u >= sum)
      sum += (*metapop).pop[++(k_pop)].size;

    unsigned int k_new = (*metapop).nb_pop; // pop nb of the new variant
  
    (*metapop).pop[k_new].island_id = (*metapop).pop[k_pop].island_id;
    (*metapop).pop[k_new].species_id = (*metapop).pop[k_pop].species_id;
    (*metapop).pop[k_new].speciation_clock = (*par).tau;
    (*metapop).pop[k_new].size = 1;
    (*metapop).pop[k_new].origin = 3;   // 1=mig conti; 2=mig island; 3=clado
    (*metapop).pop[k_new].date_origin = t;
    (*metapop).pop[k_new].divergence_duration = (*par).tau; // no gene flow yet
   
    // update population and metapopulation properties
    ((*metapop).pop[k_pop].size)--;
    ((*metapop).nb_pop)++;

    // population rates does not depend on species identity and do not need to be updated here
    // only the death rate of pops needs to be updated (different partition but same total)
    (*metapop).pop[k_pop].delta = (*metapop).pop[k_pop].size * (*par).beta*exp(-(*par).gamma*(1.0 - (*metapop).local_size[(*metapop).pop[k_pop].island_id]/(*arch).K_adj[(*metapop).pop[k_pop].island_id]));
    (*metapop).pop[k_new].delta = (*metapop).pop[k_new].size * (*par).beta*exp(-(*par).gamma*(1.0 - (*metapop).local_size[(*metapop).pop[k_new].island_id]/(*arch).K_adj[(*metapop).pop[k_new].island_id]));    
    
    // check if the indiv initiating a new lineage was the only indiv of its pop
    // if the pop is extinct, replace it in the table with the last pop of the table
    adjust_table(metapop, k_pop, rates, phylo, t);

    // record that a cladogenesis start on island (*metapop).pop[k_new].island_id occurred
    ((*rates).n_clado_start[(*metapop).pop[k_new].island_id])++;

  }
}

// THIS FUNCTION ACHIEVES A MIGRATION EVENT BETWEEN ISLANDS
void migration_islands(gsl_rng *rgsl, const double t, const double t_max, struct parameters *par, struct archipelago *arch, struct meta_population *metapop, struct outputs *rates, struct phylogeny *phylo) {
  // This function will bug if called with par.nb_islands=1. But when par.nb_islands=1,
  // metapop.mui always equals 0 (carefully checked), so that this function is is never called
  // Idem if only par.nb_islands>1 but only 1 island is emerged

  // choose the individual migrating, proportionally to
  // the frequency of all pops
  double u = gsl_ran_flat(rgsl, 0.0, (*metapop).total_size);
  unsigned int k_pop = 0;
  long unsigned int sumi = (*metapop).pop[k_pop].size;
  while (u >= sumi)
    sumi += (*metapop).pop[++(k_pop)].size;
  
  // shortcuts ; k_pop is the pop number of the migrant
  unsigned int k_ori  = (*metapop).pop[k_pop].island_id;  // island number of the origin  
  unsigned int k_spe  = (*metapop).pop[k_pop].species_id; // species number of the migrant
  
  // choose on which island the migrant arrives
  unsigned int k_isl;     // island nb where the migrant goes  
  u = gsl_ran_flat(rgsl, 0.0, (*arch).local_mui[k_ori] / (*par).mui);  // choose the destination
  if (k_ori != 0) 
    k_isl = 0;
  else
    k_isl = 1;
  double sum = sqrt((*arch).K_adj[k_isl])/(PI_TH*(*arch).di[k_ori][k_isl]) * disp_kernel((*par).ci,(*arch).di[k_ori][k_isl]); 
  while (u >= sum) {
    k_isl++;
    if (k_isl != k_ori)
      sum += (sqrt((*arch).K_adj[k_isl])/(PI_TH*(*arch).di[k_ori][k_isl]) * disp_kernel((*par).ci,(*arch).di[k_ori][k_isl])); 
  }  

  // check if this species (indep of the variant) is already present in the chosen island
  // (we retrieve the pop id of all variants of the species on this island)
  unsigned int list_k_pop[NB_POP_MAX]; // list of pops of the species on this island
  long unsigned int n_k_pop = 0;  // number of pops of the species on this island
  long unsigned int i;
  for (i=0; i<(*metapop).nb_pop; i++) {
    if ((*metapop).pop[i].island_id == k_isl && (*metapop).pop[i].species_id == k_spe) {
      list_k_pop[n_k_pop] = i;
      n_k_pop++;
    }
  }  

  // update the properties of the population on the island of destination
  
  // if the species is already present, update its properties
  if (n_k_pop != 0) {

    // choose with which variant the migrant hybridize (with prob
    // proportional to the relative pop size of each variant)
    long unsigned int all_k_pop_size = 0;
    unsigned int k_hyb; // the pop nb (in the list list_k_pop) with whom migrant hybridizes
    for (k_hyb=0; k_hyb<n_k_pop; k_hyb++)
      all_k_pop_size +=  (*metapop).pop[list_k_pop[k_hyb]].size;
    u = gsl_ran_flat(rgsl, 0.0, all_k_pop_size);
    k_hyb = 0;
    sum = (*metapop).pop[list_k_pop[k_hyb]].size;
    while (u >= sum)
      sum += (*metapop).pop[list_k_pop[++(k_hyb)]].size;
   
    // update the properties of the pop (gene flow has a delaying effect only if the pop has not yet
    // speciated, ie if its clock is not >tmax)
    if ((*metapop).pop[list_k_pop[k_hyb]].speciation_clock < t_max) {
      double spe_clock = (*metapop).pop[list_k_pop[k_hyb]].speciation_clock; // clock before mig 
      (*metapop).pop[list_k_pop[k_hyb]].speciation_clock = min_double((*metapop).pop[list_k_pop[k_hyb]].speciation_clock + ((*par).alpha / (*metapop).pop[list_k_pop[k_hyb]].size), (*par).tau); // new clock
      ((*metapop).pop[list_k_pop[k_hyb]].divergence_duration) += ((*metapop).pop[list_k_pop[k_hyb]].speciation_clock - spe_clock); // effective delay (may be < alpha/N because clock truncated at tau)
    }
    ((*metapop).pop[list_k_pop[k_hyb]].size)++;    
    ((*metapop).pop[k_pop].size)--;
    
    // update the properties of the metapop
    ((*metapop).local_size[k_isl])++;
    ((*metapop).local_size[k_ori])--;

    // update the rates
    // beta and sigma depend only on pop size, which did not change, so no need to update them
    // muc depends on the landscape only, no need to update it
    // mui depends on where are individuals and must therefore be update
    // delta must be updated for all pops on the island of origin and destination
    // I don't call update_rates() here because it does more than necessary and consider a
    // change in pop size on one island only. 
    (*metapop).mui = 0.0;
    for (i=0; i<(*par).nb_islands; i++)
      (*metapop).mui += ((*arch).local_mui[i] * (*metapop).local_size[i]);
    (*metapop).delta = 0.0;  
    for (i=0; i<(*metapop).nb_pop; i++) {
      if ((*metapop).pop[i].island_id == k_ori || (*metapop).pop[i].island_id == k_isl) 
	(*metapop).pop[i].delta = (*metapop).pop[i].size * (*par).beta*exp(-(*par).gamma*(1.0 - (*metapop).local_size[(*metapop).pop[i].island_id]/(*arch).K_adj[(*metapop).pop[i].island_id]));      
      (*metapop).delta += (*metapop).pop[i].delta;
    }

  }
  else { // if the species is new on the island, add it

    // check that memory limits are not exceeded
    if ((*metapop).nb_pop == NB_POP_MAX) {
      Rprintf("\n-------\n!ERROR!\n-------\nNumber of populations too high, NB_POP_MAX too low\n\n");
      exit(1);
    }
    else { // this migration event is possible
    
      unsigned int k_new = (*metapop).nb_pop; // k_new is defined as the next number for a new pop

      // properties of this new pop
      (*metapop).pop[k_new].island_id = k_isl;
      (*metapop).pop[k_new].species_id = k_spe;
      (*metapop).pop[k_new].speciation_clock = (*par).tau;
      (*metapop).pop[k_new].size = 1;     
      (*metapop).pop[k_new].origin = 2;                    // 1=mig conti; 2=mig island; 3=clado
      (*metapop).pop[k_new].date_origin = t;
      (*metapop).pop[k_new].divergence_duration = (*par).tau;  // no gene flow yet

      // update the properties of the metapop
      ((*metapop).pop[k_pop].size)--;
      ((*metapop).local_size[k_ori])--;
      ((*metapop).local_size[k_isl])++;
      ((*metapop).nb_pop)++;     

      // update the rates
      // beta and sigma depend only on pop size, which did not change, so no need to update them
      // muc depends on the landscape only, no need to update it
      // mui depends on where are individuals and must therefore be update
      // delta must be updated for all pops on the island of origin and destination (except the new one which is 
      // I don't call update_rates() here because it does more than necessary and consider a
      // change in pop size on one island only. 
      (*metapop).mui = 0.0;
      for (i=0; i<(*par).nb_islands; i++)
	(*metapop).mui += ((*arch).local_mui[i] * (*metapop).local_size[i]);
      (*metapop).delta = 0.0;  
      for (i=0; i<(*metapop).nb_pop; i++) {
	if ((*metapop).pop[i].island_id == k_ori || (*metapop).pop[i].island_id == k_isl) 
	  (*metapop).pop[i].delta = (*metapop).pop[i].size * (*par).beta*exp(-(*par).gamma*(1.0 - (*metapop).local_size[(*metapop).pop[i].island_id]/(*arch).K_adj[(*metapop).pop[i].island_id]));      
	(*metapop).delta += (*metapop).pop[i].delta;
      }
     
      // record that an anagenesis start due to inter-island mig to island k_isl occurred
      ((*rates).n_ana_isl_start[k_isl])++;

    }      
  }

  // check if the migrant was the only indiv of its pop (ie extinction)
  // if the pop is extinct, replace it in the table with the last pop of the table
  adjust_table(metapop, k_pop, rates, phylo, t);

  // record that a migration event from island k_ori to island k_isl occurred
  ((*rates).n_mig_isl[k_ori][k_isl])++;

}

// THIS FUNCTION ACHIEVES A MIGRATION EVENT FROM THE CONTINENT
void migration_continent(gsl_rng *rgsl, const double t, const double t_max, struct parameters *par, struct continent *conti, struct archipelago *arch, struct meta_population *metapop, struct outputs *rates) {
     
  // choose the species migrating from the continent, proportionally to
  // the frequency of each species on the continent
  double u = gsl_ran_flat(rgsl, 0.0, 1.0);
  long unsigned int k_spe = 0;       // species ID of the migrant
  while (u >= (*conti).cumfreq[k_spe])  
    k_spe++;
  
  // choose on which island the migrant arrives
  unsigned int k_isl = 0;
  double sum;
#if ISL_ARCH != 2 
  u = gsl_ran_flat(rgsl, 0.0, (*metapop).muc / (*conti).mucNc ); 
  sum = sqrt((*arch).K_adj[k_isl])/(PI_TH*(*arch).dc[k_isl]) * disp_kernel((*par).cc,(*arch).dc[k_isl]); 
  while (u >= sum) {
    k_isl++;
    sum += (*arch).local_muc[k_isl]; 
  }
#endif // for the island archipelago 2 model, k_isl = 0 is forced. Because (*arch).local_muc
       // is not consistent with the way it is computed here, and it must not be computed
       // based on Kangle since the several island are virtual in this case
  
  // check if this species (indep of the variant) is already present in the chosen island
  // (we retrieve the pop id of all variants of the species on this island)
  unsigned int k_pop;
  unsigned int list_k_pop[NB_POP_MAX]; // list of pops of the species on this island
  long unsigned int n_k_pop = 0; // number of pops of the species on this island
  for (k_pop=0; k_pop<(*metapop).nb_pop; k_pop++) {
    if ((*metapop).pop[k_pop].island_id == k_isl && (*metapop).pop[k_pop].species_id == k_spe) {
      list_k_pop[n_k_pop] = k_pop;
      n_k_pop++;
    }
  }  
	
  // if the species is already present, update its properties
  if (n_k_pop != 0) {

    // choose with which variant the migrant hybridize (with prob
    // proportional to the relative pop size of each variant)
    long unsigned int all_k_pop_size = 0;
    for (k_pop=0; k_pop<n_k_pop; k_pop++)
      all_k_pop_size +=  (*metapop).pop[list_k_pop[k_pop]].size;
    u = gsl_ran_flat(rgsl, 0.0, all_k_pop_size);
    k_pop = 0;
    sum = (*metapop).pop[list_k_pop[k_pop]].size;
    while (u >= sum)
      sum += (*metapop).pop[list_k_pop[++(k_pop)]].size;
   
    // update the properties of the pop (gene flow has a delaying effect only if the pop has not yet
    // speciated, ie if its clock is not >tmax)
    if ((*metapop).pop[list_k_pop[k_pop]].speciation_clock < t_max) {
      double spe_clock = (*metapop).pop[list_k_pop[k_pop]].speciation_clock; // clock before mig
      (*metapop).pop[list_k_pop[k_pop]].speciation_clock = min_double((*metapop).pop[list_k_pop[k_pop]].speciation_clock + ((*par).alpha / (*metapop).pop[list_k_pop[k_pop]].size), (*par).tau); // new clock
      ((*metapop).pop[list_k_pop[k_pop]].divergence_duration) += ((*metapop).pop[list_k_pop[k_pop]].speciation_clock - spe_clock); // effective delay (may be < alpha/N because clock truncated at tau)
    }
    ((*metapop).pop[list_k_pop[k_pop]].size)++;

    // update the properties of the metapop
    ((*metapop).local_size[k_isl])++;
    ((*metapop).total_size)++;

    // update total rates
    update_rates(par, arch, metapop, k_isl);

  }
  else { // if the species is new on the island, add it
    // check that memory limits are not exceeded
    if ((*metapop).nb_pop == NB_POP_MAX) {
      Rprintf("\n-------\n!ERROR!\n-------\nNumber of populations too high, NB_POP_MAX too low\n\n");
      exit(1);
    }
    else { // this migration event is possible

      // here k_pop = (*metapop).nb_pop, ie k_pop is the next number for a new pop
  
      // properties of this new pop
      (*metapop).pop[k_pop].island_id = k_isl;
      (*metapop).pop[k_pop].species_id = k_spe;
      (*metapop).pop[k_pop].speciation_clock = (*par).tau;
      (*metapop).pop[k_pop].size = 1;
      (*metapop).pop[k_pop].origin = 1;      // 1=mig conti; 2=mig island; 3=clado
      (*metapop).pop[k_pop].date_origin = t;
      (*metapop).pop[k_pop].divergence_duration = (*par).tau;  // no gene flow yet
	  
      // update the properties of the metapop
      ((*metapop).local_size[k_isl])++;
      ((*metapop).total_size)++;
      ((*metapop).nb_pop)++;

      // update total rates
      update_rates(par, arch, metapop, k_isl);
      
      // record that an anagenesis start due to mig from the conti to island k_isl occurred
      ((*rates).n_ana_conti_start[k_isl])++;

    }
  }

  // record that a migration event from the continent to island k_isl occurred
  ((*rates).n_mig_conti[k_isl])++;
}

// MAIN SIMULATION FUNCTION
void simul(int *p1, double *p2, double *p3, double *p4, double *p5, double *p6, double *p7, double *p8, double *p9, double *p10, int *p11, double *p12, int *p13, double *p14, double *p15, double *p16, double *p17, double *p18, double *p19, double *p20, double *p21, double *p22, int *p23, int *p24, double *p25, double *p26, double *p27, double *p28, double *p29, double *p30, double *p31, double *p32, double *p33)
{ 
 
  register unsigned int i, j; // index used many times
  
  // give understandable names to parameters (and cast them if necessary)
  char sim_name[100];        // name of the simulation
  char tmp[100];
  sprintf(tmp, "%d", *p1);
  strcpy(sim_name, "sim");
  strcat(sim_name, tmp);
  struct parameters par;
  par.beta = *p2;
  par.gamma = *p3;
  par.sigma = *p4;
  par.alpha = *p5;
  par.tau = *p6;
  par.mui = *p7;
  par.ci = *p8;
  par.muc = *p9;
  par.cc = *p10;
  par.Nc = (long unsigned int) *p11;
  par.theta = *p12;
  par.nb_islands = (unsigned int) *p13; 
  alloc_matrix_double(&(par.di), par.nb_islands, par.nb_islands);
  for (i=0; i<par.nb_islands; i++) {
    par.Kmax[i] = p14[i];
    par.ti[i] = p15[i];
    par.pmax[i] = p16[i];
    par.tl[i] = p17[i];
    par.tm[i] = par.ti[i] + par.pmax[i]*par.tl[i];
    par.ts[i] = par.ti[i] + par.tl[i];
    par.dc[i] = p18[i];
    for (j=0; j<par.nb_islands; j++)
      par.di[j][i] = p19[j+par.nb_islands*i];
  }
  par.phi_max = *p27;
  par.d = *p28;
  par.amp = *p29;
  par.period = *p30;
  par.phi_min = *p31;
  par.psi = *p32;
  par.h_breakpoint = *p33;
  double t_max = *p20;
  const double wt_ss = *p21;
  const double wt_up = *p22;
  const unsigned int n_run = (unsigned int) *p23;
  long unsigned int rng_seed = (long unsigned int) *p24;
#if STATIC
  const double t_eq = *p25;
#endif
  const double t_ini = *p26;

  // time when simulation starts
  struct duration t_exec;
  unsigned int t_start;
  t_start = (unsigned) time (NULL);

  // Loop on runs
  unsigned int i_run;
  for (i_run = 0; i_run < n_run; i_run++) {

    // Initialize simulation

    // time when this replicate starts
    unsigned int t_start_run;
    t_start_run = (unsigned) time (NULL);    

    // initialize the random number generator
    const gsl_rng_type *Tgsl;
    gsl_rng *rgsl;
    gsl_rng_env_setup();            // If the seed is set to 0, then
    Tgsl = gsl_rng_default;         // it is randomly chosen (equal to
    rgsl = gsl_rng_alloc (Tgsl);    // current time)
    long unsigned int rng_seed_run;
    if (rng_seed == 0)
      rng_seed_run = (unsigned) time (NULL);
    else
      rng_seed_run = rng_seed;
    gsl_rng_set (rgsl, rng_seed_run);

    // clear output files
    clear_file("out_log");
    clear_file("out_conti");
    clear_file("out_arch");
    clear_file("out_metapop");
    clear_file("out_rates");
    clear_file("out_phylo");


    // print (on screen and log file) that run begins and the seed value
    Rprintf("Begining run %d on %d of simulation AAA [%s]\n",
	   i_run + 1, n_run, sim_name);
    Rprintf("rng_seed: %lu\n", rng_seed_run);
    FILE *log_file = NULL;
    log_file = fopen("out_log", "a");
    fprintf(log_file, "Begining run %d on %d of simulation [%s]\n",
	   i_run + 1, n_run, sim_name);
    fprintf(log_file, "rng_seed:\t%lu\n", rng_seed_run);
    fclose(log_file);
    
    // initializations

    // initialize times
    double t = t_ini; // current time
    double t_ss = t_ini + wt_ss;           // next time for saving outputs
    double t_up = t_ini + wt_up;           // next time for updating the landscape and rates
#if STATIC
    t_max = t_ini + t_eq;
#endif

    // initialize the continent
    struct continent conti;
    initialize_continent(rgsl, &conti, &par);

    // initialize the archipelago
    struct archipelago arch;
    alloc_matrix_double(&(arch.di), par.nb_islands, par.nb_islands);
    initialize_archipelago(&par, &arch, t);

    // initialize the metapop
    struct meta_population metapop;      
    initialize_metapop(&par, &arch, &metapop, &conti, t_max);
    long unsigned int new_sp_id = conti.nb_species; // next id for new species

    // initialize the outputs (rates)
    struct outputs rates;
    alloc_matrix_luint(&(rates.n_mig_isl), par.nb_islands, par.nb_islands);
    initialize_rates(&par, &rates);
    
    // allocate memory for struct phylogeny
    struct phylogeny phylo;
    phylo.species_id=(long unsigned int *) malloc(NB_SP_PHYLO_MAX * sizeof (long unsigned int));
    if (phylo.species_id == NULL) // if allocation doesn't work
    { 
         Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed for phylo.species\n\n");
         exit(0);
    }
    phylo.date_origin = (double *) malloc(NB_SP_PHYLO_MAX * sizeof (double));
    if (phylo.date_origin == NULL) // if allocation doesn't work
    { 
         Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed for phylo.date_origin\n\n");
         exit(0);
    }
    phylo.mother = (long unsigned int *) malloc(NB_SP_PHYLO_MAX * sizeof (long unsigned int));
    if (phylo.mother == NULL) // if allocation doesn't work
    { 
         Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed for phylo.mother\n\n");
         exit(0);
    }
    phylo.birth = (double *) malloc(NB_SP_PHYLO_MAX * sizeof (double));
    if (phylo.birth == NULL) // if allocation doesn't work
    { 
         Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed for phylo.birth\n\n");
         exit(0);
    }
    phylo.death = (double *) malloc(NB_SP_PHYLO_MAX * sizeof (double));  
    if (phylo.death == NULL) // if allocation doesn't work
    { 
         Rprintf("\n-------\n!ERROR!\n-------\nMatrix allocation failed for phylo.death\n\n");
         exit(0);
    }
    initialize_phylo(&phylo);

    // save initial state
    print_conti(&conti);
    print_arch(&par, &arch, t);
    print_metapop(&metapop, t);
    print_rates(&par, &rates, t);

    //print_time(t);
     //Rprintf("COUCOUABCD (%g, %g, %lf)\n", t, t_max, *arch.dc, *arch.di);
             
    // Loop on time for this run
    while (t < t_max) {    
      
      // time of next event and waiting time to it
      double total_rate = metapop.beta + metapop.delta + metapop.sigma + metapop.mui + metapop.muc;
	
      double wt_ev = gsl_ran_exponential(rgsl, 1.0 / total_rate);

      // first event to occur
      unsigned int it = index_min_of_three(t_up, t_ss, t+wt_ev);
      
      // time to update the landscape, the rates, the speciation clocks and the achieve speciation events
      if (it == 1) {
	t = t_up;            // update time
	t_up += wt_up;

//	Rprintf("t=%.6f\tm.b=%.1f\tm.d=%.6f\tm.s=%.6f\tm.mi=%.10f\tm.mc=%.6f\t", t, metapop.beta, metapop.delta, metapop.sigma, metapop.mui, metapop.muc);	
//	Rprintf("\nm.np=%u\tm.ts=%lu\n",metapop.nb_pop, metapop.total_size);
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("m.ls[%u]=%lu\t", i, metapop.local_size[i]);	 
//	Rprintf("\n");
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("arch.K[%u]=%g\t", i, arch.K[i]);	 
//	Rprintf("\n");
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("locmuc[%u]=%g\t", i, arch.local_muc[i]);	 
//	Rprintf("\n");
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("locmui[%u]=%g\t", i, arch.local_mui[i]);	 
//	Rprintf("\n");
//	for (i=0; i<metapop.nb_pop; i++) 
//	  Rprintf("pop=%u\tisl=%u\tsp=%lu\tcloc=%.6f\tsiz=%lu\tdel=%.6f\n", i, metapop.pop[i].island_id, metapop.pop[i].species_id, metapop.pop[i].speciation_clock, metapop.pop[i].size, metapop.pop[i].delta);	  

//	Rprintf("update K and co\n");
	
	update_landscape(&par, &conti, &arch, &metapop, &rates, t, wt_up, t_max, &new_sp_id, &phylo);

//	Rprintf("t=%.6f\tm.b=%.1f\tm.d=%.6f\tm.s=%.6f\tm.mi=%.10f\tm.mc=%.6f\t", t, metapop.beta, metapop.delta, metapop.sigma, metapop.mui, metapop.muc);	
//	Rprintf("\nm.np=%u\tm.ts=%lu\n",metapop.nb_pop, metapop.total_size);
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("m.ls[%u]=%lu\t", i, metapop.local_size[i]);	 
//	Rprintf("\n");
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("arch.K[%u]=%g\t", i, arch.K[i]);	 
//	Rprintf("\n");
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("locmuc[%u]=%g\t", i, arch.local_muc[i]);	 
//	Rprintf("\n");
//	for (i=0; i<par.nb_islands; i++)
//	  Rprintf("locmui[%u]=%g\t", i, arch.local_mui[i]);	 
//	Rprintf("\n");
//	for (i=0; i<metapop.nb_pop; i++) 
//	  Rprintf("pop=%u\tisl=%u\tsp=%lu\tcloc=%.6f\tsiz=%lu\tdel=%.6f\n", i, metapop.pop[i].island_id, metapop.pop[i].species_id, metapop.pop[i].speciation_clock, metapop.pop[i].size, metapop.pop[i].delta);
//	Rprintf("--------------------------------------------\n");
	
      }
      
      // time to save the state of the metapopulation
      else if (it == 2) {
	t = t_ss;            // update times
	t_ss += wt_ss;
	
//	Rprintf("save output\n");
	print_arch(&par, &arch, t);
	print_metapop(&metapop, t);
	print_rates(&par, &rates, t);
//	print_time(t);
      }
      
      // an event among the 5 below is the next event
      else {

	t += wt_ev;          // update times
	
	// choose the event which occurs and do it
	double u = gsl_ran_flat(rgsl, 0.0, total_rate);
	if (u < metapop.beta) {
//	  Rprintf("birth\n");
	  birth(rgsl, &metapop, &par, &arch);	  
	}
	else if (u < metapop.beta + metapop.delta) {
//	  Rprintf("death\n");
	  death(rgsl, &metapop, &rates, &par, &arch, &phylo, t);	  	  
	}
	else if (u < metapop.beta + metapop.delta + metapop.sigma) {
//	  Rprintf("cladogenesis\n");
	  cladogenesis_start(rgsl, t, &par, &metapop, &rates, &arch, &phylo);
	}
	else if (u < metapop.beta + metapop.delta + metapop.sigma + metapop.mui) {
//	  Rprintf("migration between islands\n");
	  migration_islands(rgsl, t, t_max, &par, &arch, &metapop, &rates, &phylo);
	  
	}
	else {
//	  Rprintf("migration from continent\n");
	  migration_continent(rgsl, t, t_max, &par, &conti, &arch, &metapop, &rates);	  
	}
	
      }

    } // end of simulation loop

    // Save phylogeny
    print_phylo(&phylo);
    // printf("print phylo \n");
    // Terminate the run

    // Free allocated memory
    for (i=0; i<par.nb_islands; i++) {
      free(arch.di[i]);
      free(rates.n_mig_isl[i]);
    }
    free(arch.di);
    free(rates.n_mig_isl);
    gsl_rng_free(rgsl);

    free(phylo.species_id);
    free(phylo.date_origin);
    free(phylo.mother);
    free(phylo.birth);
    free(phylo.death);

    // Compute and print (on screen and disk) the duration of the run
    t_exec = time_exec(t_start_run, (unsigned) time (NULL));
    Rprintf("End of run %d on %d of simulation [%s] ",
	   i_run + 1, n_run, sim_name);
    Rprintf("after %dd %dh %dm %ds\n",
	     t_exec.d, t_exec.h, t_exec.m, t_exec.s);
    log_file = fopen("out_log", "a");
    fprintf(log_file, "End of run %d on %d of simulation [%s] ",
    	   i_run + 1, n_run, sim_name);
    fprintf(log_file, "after %dd %dh %dm %ds\n",
    	     t_exec.d, t_exec.h, t_exec.m, t_exec.s);
    fclose(log_file);

    // save output files
    mv_files(sim_name, i_run);
    
  } // end of the loop on runs
  
  // Free allocated memory
  for (i=0; i<par.nb_islands; i++)
    free(par.di[i]);
  free(par.di);
  

  //Print total duration on screen
  t_exec = time_exec(t_start, (unsigned) time (NULL));
  Rprintf("End of all replicates of simulation [%s] ", sim_name);
  Rprintf("after %dd %dh %dm %ds\n",
	 t_exec.d, t_exec.h, t_exec.m, t_exec.s);
  
}


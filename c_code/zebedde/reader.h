/****************************************************************************************/
/**** Define Key words for reader Last altered Oct 2006      Dave Willock ***************/
/**** Redundant keyword comments removed and extra for bio version added  ***************/
/**** Last altered Alan Lobo Sept 2008. Added new keywords, EXXON,		  ***************/
/**** DOCK_ENERGY, WANT_CONDOR FINAL_CAR								  ***************/
/****************************************************************************************/

enum { PRIME_DIRECTIVE= 1, SECONDARY_DIRECTIVE, TERTIARY_DIRECTIVE, REMAINING_DIRECTIVE };

enum {TITLE = 1, 
	PORE_FILE, STRATEGY_FILE, MAX_TEMPLATES,  VDW_SCALE,
	LIBRARY, POTS_FILE, GULP_CMD_LINE, GULP_EXTRA, 
	SEED_TEMPLATE,ANIMATION,NB_CUTOFF,CH_CUTOFF,VERBOSE,WEIGHTS,
	MODIFY,  SHAKE,	ROCK, STOP_CUTOFF,BOX, SYMMETRY,
	MINIMIZER, MOPAC, GULP, DISCOVER, C2DISCOVER, PROB_TEST, FORBIDDEN_BOND, 
	CONCENTRATION_LIMIT, CENTRE_TEMPLATE, INITIAL_MINIMIZE,
        RING_CUTOFF, CHECK_INPUT, ANALYSE, PEEK_FREQ, RANDOM,
        DEFAULTS_FILE, COST_FUNCTION, JUST_POSITION, FORCE_DIHEDRALS,
        LINE_RESTRAINT, TETHER_RESTRAINT, LOGFILE, CONVENTION, SLAB, POLY,
        TEMPERATURE, ALLOWED_TORSIONS, DOCK, MC_OPT, MC_TEMP_STEP, MC_NO_ACTIONS,
	FINAL_GULP, MONITOR, DOCK_ENERGY, WANT_CONDOR, EXXON, FINAL_CAR, NUM_DOCK, ANNEAL, UCL_};

enum {ATTEMPTS = 101, 
	STEP, TEMPLATE, INPORE, FRAGMENT, ACTION, LIMITS, /*ADDED AJWL*/ 
	NAME, FRACTION,  PATH, FORCEFIELD, MOLECULE, ARCHIVE,
	HYDROGEN, BIOSYM, PXYZ, TXYZ, STERIC2, NON_BOND2, ENERGY2, ON2, OFF2, YES2, NO2,
        ELEMENTS, CONSTRAIN, FIX, SEED, DEBUG_IT, XYZ, ZYX };

enum {MAYBE = 201, 
        ALL, ON3, OFF3, YES3, NO3, HOLD, ATOM};

enum {DEFAULT = 666};    
enum {UNIT = 998};    
enum {PARSE = 1001};

enum {BLANK_DIRECT = 999};

/*** note two defintions of non_bonded as nonb and non_ 21.10.96 DWL ***/
#define FIRST_DIRECTIVE_LIST \
	{"titl", TITLE},           {"pore", PORE_FILE},      \
	{"stra", STRATEGY_FILE},  {"maxt", MAX_TEMPLATES}, \
	{"vdw_", VDW_SCALE}, \
	{"libr", LIBRARY}, {"g_po", POTS_FILE}, \
         {"gcom", GULP_CMD_LINE}, {"gext", GULP_EXTRA}, \
        {"seed", SEED_TEMPLATE},\
	{"anim", ANIMATION}, {"nb_c", NB_CUTOFF}, {"ch_c", CH_CUTOFF}, \
	{"verb", VERBOSE}, {"weig", WEIGHTS}, \
    {"modi", MODIFY}, {"shak", SHAKE}, {"rock", ROCK}, \
    {"stop", STOP_CUTOFF}, {"gbox", BOX},\
    {"symm",SYMMETRY}, {"mini", MINIMIZER},{"mopa", MOPAC}, {"gulp", GULP}, \
    {"disc", DISCOVER}, {"c2di", C2DISCOVER}, \
    {"prob",PROB_TEST}, {"forb", FORBIDDEN_BOND}, \
	{"conc", CONCENTRATION_LIMIT}, \
    {"cent", CENTRE_TEMPLATE}, \
    {"init", INITIAL_MINIMIZE}, \
        {"ring", RING_CUTOFF}, \
        {"chec", CHECK_INPUT},\
        {"anal", ANALYSE}, \
        {"peek", PEEK_FREQ}, \
        {"rand", RANDOM}, \
        {"defa", DEFAULTS_FILE}, {"cost", COST_FUNCTION}, \
        {"just", JUST_POSITION},{"dihe", FORCE_DIHEDRALS},{"line", LINE_RESTRAINT},\
         {"teth", TETHER_RESTRAINT}, {"logf", LOGFILE}, {"conv", CONVENTION}, {"slab", SLAB},\
         {"poly", POLY}, {"temp",TEMPERATURE},{"tors", ALLOWED_TORSIONS}, \
         {"mcop", MC_OPT}, {"tste", MC_TEMP_STEP}, {"amod", MC_NO_ACTIONS}, {"fgul", FINAL_GULP},\
         {"moni", MONITOR},{"doke", DOCK_ENERGY}, {"cond", WANT_CONDOR}, {"exxo", EXXON}, {"fcar",FINAL_CAR}, \
         {"numd", NUM_DOCK},{"anne", ANNEAL}, {"ucl_",UCL_}, {"",BLANK_DIRECT}

#define SECOND_DIRECTIVE_LIST \
	    {"atte", ATTEMPTS}, {"step", STEP}, {"temp", TEMPLATE}, {"in_p", INPORE},\
	    {"frag", FRAGMENT}, {"acti", ACTION}, \
	    {"limi", LIMITS},   {"frac", FRACTION}, {"name", NAME}, \
	    {"path", PATH},     {"forc", FORCEFIELD}, \
	    {"mole", MOLECULE}, {"arch", ARCHIVE}, \
	    {"hydr", HYDROGEN},\
            {"bios", BIOSYM}, {"pxyz", PXYZ}, {"txyz", TXYZ} , {"ster", STERIC2}, \
            {"non_", NON_BOND2}, {"nonb", NON_BOND2},     {"ener", ENERGY2},\
            {"on", ON2}, {"off", OFF2}, {"yes", YES2}, {"no", NO2}, \
            {"elem", ELEMENTS}, {"cons", CONSTRAIN},     {"seed", SEED}, \
	    {"debu", DEBUG_IT}, {"xyz", XYZ}, {"zyx", ZYX}, {"dock", DOCK}, {"",BLANK_DIRECT}

/*** NOTE replacing of yes/no with on/off which solves probs  21.10.96 DWL ***/
#define THIRD_DIRECTIVE_LIST \
           {"ma", MAYBE},\
               {"al", ALL},\
               {"on", ON3}, \
               {"of", OFF3}, \
               {"ye", YES3}, \
               {"no", NO3}, \
	       {"ho", HOLD}, {"at", ATOM}, {"", BLANK_DIRECT}

#define NULL_DIRECTIVE_LIST \
       {"an", UNIT}, {"(a", UNIT}, \
       {"(c", UNIT}, {"k ", UNIT}, {"(k",  UNIT}, \
       {"#", PARSE}, {"", BLANK_DIRECT}

typedef struct 
{
  char *directive;
  int token_index;
}list;


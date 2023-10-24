#ifndef CONFIG_H
#define CONFIG_H

// DM : debug message -- disable for now
// #define DM(x) std::cerr << x
#define DM(x)

// CDM: custom debug message -- this is used only when certain DM messages
// are to be printed for help in debugging -- default it is disabled
// #define CDM(x) std::cerr << x 
#define CDM(x)

// #define CONFIG_COLLECT_SCHEDULE_GRAPH
/*
#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
#define CONFIG_PARALLEL
#endif
*/
#ifndef NDEBUG
#define TBB_USE_DEBUG 1
#endif

#endif

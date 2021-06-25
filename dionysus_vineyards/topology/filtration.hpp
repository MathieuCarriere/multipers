#include <utilities/log.h>

#ifdef LOGGING
static rlog::RLogChannel* rlFiltration =                    DEF_CHANNEL("topology/filtration/info", rlog::Log_Debug);
static rlog::RLogChannel* rlFiltrationDebug =               DEF_CHANNEL("topology/filtration/debug", rlog::Log_Debug);
#endif // LOGGING

#ifndef __MEMORY_H__
#define __MEMORY_H__

#if LOGGING     // TODO: add check for Linux (and preferably the right version of Linux)

#include "log.h"
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

static rlog::RLogChannel* rlMemory =                   DEF_CHANNEL("memory",    rlog::Log_Debug);

unsigned    report_memory()
{
    pid_t pid = getpid();
    std::stringstream smaps_name;
    smaps_name << "/proc/" << pid << "/smaps";
    std::ifstream in(smaps_name.str().c_str());

    std::string str;
    unsigned memory = 0, value;
    while (in)
    {
        in >> str;
        if (std::string(str, 0, 7) == "Private")
        {
            in >> value;
            memory += value;
        }
    }
    rLog(rlMemory, "Private memory usage: %d kB", memory);
    
    return memory;
}

#else

unsigned report_memory() { return 0; }

#endif // LOGGING

#endif // __MEMORY_H__

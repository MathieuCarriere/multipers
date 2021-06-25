#ifndef __TIMER_H__
#define __TIMER_H__

// Adapted with minor changes from http://oldmill.uchicago.edu/~wilder/Code/timer/

#include <ctime>
#include <iostream>
#include <iomanip>

class Timer
{
    public:
                    Timer(): start_time(0), acc_time(0)                     {}

        void        start();
        void        stop();
        void        check(const char* msg = 0) const;

    private:
        clock_t     start_clock;
        time_t      start_time;
        double      acc_time;

        double      elapsed_time() const;
};

// Return the total time that the timer has been in the "running" state since
// it was last "started".  For "short" time periods (less than an hour), the
// actual cpu time used is reported instead of the elapsed time.
inline double 
Timer::
elapsed_time() const
{
    time_t acc_sec = time(0) - start_time;
    if (acc_sec < 3600)
        return (clock() - start_clock) / (1.0 * CLOCKS_PER_SEC);
    else
        return (1.0 * acc_sec);
}

inline void 
Timer::
start()
{
    start_clock = clock();
    start_time = time(0);
}

inline void 
Timer::
stop()
{
    acc_time += elapsed_time();
}

// Print out an optional message followed by the current timer timing.
inline void 
Timer::
check(const char* msg) const
{
    // Print an optional message, something like "Checking timer t";
    if (msg) std::cout << msg << " : ";

    std::cout << "Elapsed time [" << std::setiosflags(std::ios::fixed)
              << std::setprecision(2) << acc_time << "] seconds\n";
}

#endif // __TIMER_H__


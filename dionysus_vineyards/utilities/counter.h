/*
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2005 -- 2007
 */

#ifndef __COUNTER_H__
#define __COUNTER_H__


#ifndef COUNTERS
    #define     GetCounter(path)        0
    #define     Count(x)
    #define     CountNum(x,y)
    #define     CountBy(x,y)
    #define     CountNumBy(x,y,z)
    #define     SetFrequency(x, freq)
    #define     SetTrigger(x, y)
    #define     Print(x)
#else // COUNTERS
    #define     GetCounter(path)            get_counter(path)
    #define     Count(x)                    do { x->count++; if ((x->count % x->frequency == 0)) x->trigger->print(); } while (0)
    #define     CountNum(x,y)               do { x->subcount[y]++; } while (0)
    #define     CountBy(x,y)                do { x->count += y; } while (0)
    #define     CountNumBy(x,y,z)           do { x->subcount[y] += z; } while (0)
    #define     SetFrequency(x, freq)       do { x->frequency = freq; } while (0)
    #define     SetTrigger(x, y)            do { x->trigger = y; } while(0)
    #define     Print(x)                    do { x->trigger->print(); } while(0)


#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <unistd.h>

class Counter
{
    public:
        typedef                 unsigned long                           CounterType;
        typedef                 std::map<std::string, Counter*>         SubCounterMap;
        typedef                 std::map<int, CounterType>              SubCountMap;

    public:
        CounterType             count;
        CounterType             frequency;
        SubCountMap             subcount;
        Counter*                trigger;

    public:
                                Counter(const std::string& full_name = "",
                                        CounterType freq = std::numeric_limits<CounterType>::max());
                                ~Counter();

        Counter*                get_child(const std::string& path, std::string::size_type pos);
        void                    print();

    private:    
        SubCounterMap           subcounters_;
        std::string             full_name_;
        
        static const char*      start_color;
        static const char*      finish_color;
        static const char       green_color[];
        static const char       normal_color[];
        static const char       empty_string[];
};
const char Counter::green_color[]       = "\033[32m";
const char Counter::normal_color[]      = "\033[0m";
const char Counter::empty_string[]      = "";
const char* Counter::start_color        = 0;
const char* Counter::finish_color       = 0;

static      Counter             rootCounter;

Counter*    get_counter(const char* path)
{
    return rootCounter.get_child(path, 0);
}

#include "counter.hpp"

#endif // COUNTERS


#endif // __COUNTER_H__

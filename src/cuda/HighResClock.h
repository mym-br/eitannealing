#include <chrono>

struct HighResClock
    {
        typedef long long                               rep;
        typedef std::nano                               period;
        typedef std::chrono::duration<rep, period>      duration;
        typedef std::chrono::time_point<high_resolution_clock>   time_point;
        static const bool is_steady = true;

        static time_point now() { return time_point::now(); }
    };

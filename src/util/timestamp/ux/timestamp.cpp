#include "../timestamp.h"
#include <sys/time.h>
#include <sys/resource.h>

unsigned long get_usec_timestamp()
{
	struct timeval t;
	struct timezone tzp;
	gettimeofday(&t, &tzp);
	return t.tv_sec*1e6 + t.tv_usec;
}

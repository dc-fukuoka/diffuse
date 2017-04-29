#include <time.h>

double dclock_(void)
{
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);
	return (double)(t.tv_sec+1e-9*t.tv_nsec);
}

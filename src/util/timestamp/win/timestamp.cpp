#include "../timestamp.h"
#include <Windows.h>
#include <sysinfoapi.h>

union qword_dword {
	unsigned long qword;
	struct {
		unsigned int low;
		unsigned int high;
	} dword;
};

unsigned long get_usec_timestamp()
{
	FILETIME time;
	union qword_dword timestamp;
	GetSystemTimeAsFileTime(&time);
	timestamp.dword.low = time.dwLowDateTime;
	timestamp.dword.high = time.dwHighDateTime;
		
	return timestamp.qword/10;
}
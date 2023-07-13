#pragma once

#include <sys/time.h>
#include <sys/resource.h>

long mem_info() {
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  return usage.ru_maxrss;
}


#pragma once

#include <iostream>
#include <map>
#include <ctime>

#ifdef __MACH__
#include <sys/time.h>
// clock_gettime is not implemented on OSX
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec* t) {
  struct timeval now;
  int rv = gettimeofday(&now, NULL);
  if (rv) return rv;
  t->tv_sec = now.tv_sec;
  t->tv_nsec = now.tv_usec * 1000;
  return 0;
}
#endif

class tm_entry {
 public:
  double total_cpu_time;
  double total_wall_time;
  unsigned long iterations;

  tm_entry() : total_cpu_time(0.), total_wall_time(0.), iterations(0) {}
  tm_entry(double cpu, double wall, unsigned long its = 1)
      : total_cpu_time(cpu), total_wall_time(wall), iterations(its) {}

  void update(const tm_entry& other) {
    total_cpu_time += other.total_cpu_time;
    total_wall_time += other.total_wall_time;
    iterations += other.iterations;
  }
};

class timing {
  bool recording;
  clock_t last_clock;
  timespec last_wall;
  std::map<std::string, tm_entry> results;

 public:
  timing() : recording(false) {}

  void start(void) {
    if (recording) return;

    recording = true;
    last_clock = clock();
    clock_gettime(CLOCK_MONOTONIC, &last_wall);
  }

  void stop(const std::string& s, bool verbose = true) {
    if (!recording) return;

    double cpu_secs, wall_secs;
    timespec now;
    now.tv_sec = 0;
    now.tv_nsec = 0;

    cpu_secs =
        static_cast<double>(clock() - last_clock) / ((double)CLOCKS_PER_SEC);

    clock_gettime(CLOCK_MONOTONIC, &now);
    wall_secs = (now.tv_sec - last_wall.tv_sec) / 1.0 +
                (now.tv_nsec - last_wall.tv_nsec) / 1000000000.0;

    if (verbose)
      std::cerr << s << ": " << cpu_secs << " s (CPU) / " << wall_secs
           << " s (wall)." << std::endl;

    tm_entry newentry(cpu_secs, wall_secs);

    std::map<std::string, tm_entry>::iterator it = results.find(s);
    if (it != results.end()) {
      newentry.update(it->second);
      results.erase(it);
    }

    results.insert(std::pair<std::string, tm_entry>(s, newentry));
    recording = false;
  }

  void show(void) {
    double total_cpu = 0., total_wall = 0.;
    std::map<std::string, tm_entry>::iterator it;

    std::cerr << "Function\t#\tCPU avg\tWall avg\n";
    for (it = results.begin(); it != results.end(); it++) {
      tm_entry e = it->second;
      std::cerr << it->first << '\t' << e.iterations << '\t'
           << e.total_cpu_time / e.iterations << '\t'
           << e.total_wall_time / e.iterations << std::endl;

      total_cpu += e.total_cpu_time;
      total_wall += e.total_wall_time;
    }
    std::cerr << "\nOverall time: " << total_cpu << " s (CPU) / " << total_wall
         << " s (wall)." << std::endl;
  }
};

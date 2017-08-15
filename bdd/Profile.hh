#ifndef PROFILE_HH_
#define PROFILE_HH_

#include <iostream>
#include <chrono>

/* class: Profile
 * helper class to measure elapsed time based on std::chrono library */
class Profile {
  public:
    std::chrono::duration<double> total;
    std::chrono::duration<double> max;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point stop;

    Profile(){
        total = std::chrono::duration<double>::zero();
        max = std::chrono::duration<double>::zero();
    }
    ~Profile(){}

    /* function: tic
     * set start time
     */
    inline void tic(void) {
      start=std::chrono::high_resolution_clock::now();
    }
    /* function: toc
     * set stop time and update total, max
     */
    inline void toc(void) {
        stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt;
        dt = std::chrono::duration_cast<std::chrono::duration<double > >(stop-start);
        total += dt;
        if (dt > max) {
            max = dt;
        }
    }

    inline void print(void) {
        std::cout << "Total elapsed time is " << total.count() << " seconds." << std::endl;
        std::cout << "Maximum elapsed time is " << max.count() << " seconds." << std::endl;
    }

};

#endif /* TICTOC_HH_ */


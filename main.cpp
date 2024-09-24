#include "sweep_decoding.h" 

int main(int argc, char *argv[])
{
    // Parse command-line arguments
    int level = std::atoi(argv[1]); // Number of iterations of punching holes in the 3D lattice
    const std::string sweepSchedule(argv[2]); // Determines sweep direction changes
    int rounds = std::atoi(argv[3]); // Number of sweep rounds
    int L = std::atoi(argv[4]); // Linear system size
    const double p = std::atof(argv[5]); // Physical qubit error rate
    const double q = std::atof(argv[6]); // Measurement error rate
    int seed = std::atoi(argv[7]); // Seed for random number generator

    // Initialize random number generator with the provided seed
    pcg32 rnEngine(seed);

    // Execute sweep decoder and print the result
    std::cout << sweep_decoder_run(level, sweepSchedule, rounds, L, p, q, rnEngine) << std::endl;

    return 0;
}

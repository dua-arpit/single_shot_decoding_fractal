#include "sweep_decoding.h"

int main(int argc, char *argv[])
{
    //argv indexing starts from 1
    int level=std::atoi(argv[1]);
    const std::string sweepSchedule(argv[2]);
    int rounds = std::atoi(argv[3]);
    int L=std::atoi(argv[4]); 
    const double p=std::atof(argv[5]);
    const double q=std::atof(argv[6]);
    // int maxruns_id=std::atoi(argv[7]);
    int seed=std::atoi(argv[7]);

    //pcg_extras::seed_seq_from<std::random_device> seed;
    pcg32 rnEngine(seed);

    std::cout<<sweep_decoder_run(level,sweepSchedule,rounds,L,p,q,rnEngine)<<std::endl;

    return 0;
}


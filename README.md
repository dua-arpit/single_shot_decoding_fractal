# single_shot_decoding_fractal
We adapt a cellular automaton decoder, sweep decoder, to the 3D surface code on a fractal lattice. The decoder is described in this paper https://arxiv.org/abs/2201.03568.

There are two important files. 

1) The sweep_decoding.h file contains all the functions for the decoding algorithm. 

2) The main.cpp executes the decoder via the function sweep_decoder_run. The function sweep_decoder_run takes in the following arguments: 

    a) level= level of the fractal lattice. level 0 corresponds to the 3D surface code. 

    b) sweepSchedule= a schedule of sweep directions. 

    c) rounds= number of rounds of error correction. 

    d) L= system size.

    e) p=physical error rate.

    f) q=measurement error rate.

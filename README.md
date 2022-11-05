# single_shot_decoding_fractal
We adapt a cellular automaton decoder, sweep decoder, to the 3D surface code on a fractal lattice. 
The sweep_decoding.h file contains all the functions for the decoding algorithm. 
The main.cpp executes the decoder. 
level= level of the fractal lattice. level 0 corresponds to the 3D surface code. 
sweepSchedule= a schedule of sweep directions. 
rounds= number of rounds of error correction. 
L= system size.
p=physical error rate.
q=measurement error rate.

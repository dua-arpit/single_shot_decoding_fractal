# single_shot_decoding_fractal
We adapt a cellular automaton decoder, sweep decoder, to the 3D surface code on a fractal lattice. 
The sweep_decoding.h file contains all the functions for the decoding algorithm. 

The main.cpp executes the decoder via the function sweep_decoder_run. The function sweep_decoder_run takes in the following arguments: 

1) level= level of the fractal lattice. level 0 corresponds to the 3D surface code. 

2) sweepSchedule= a schedule of sweep directions. 

3) rounds= number of rounds of error correction. 

4) L= system size.

5) p=physical error rate.

6) q=measurement error rate.

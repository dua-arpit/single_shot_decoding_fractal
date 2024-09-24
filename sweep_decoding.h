#include "preamble/preamble.h" 

//Structure representing the lattice with vertices, edges, faces, cubes, and logical coordinates
struct lattice
{
    vint vertices, edges, faces, cubes;
    vint logicalcoords;
};

//Convert 3D coordinates to a single integer coordinate
int_fast64_t coordinate(int_fast64_t L, vint vec)
{
    int_fast64_t x = vec[0], y = vec[1], z = vec[2];
    int_fast64_t Lp = 2 * L + 3;
    return x + y * Lp + z * Lp * Lp;
}

//Convert a single integer coordinate back to 3D coordinates
vint xyz(int_fast64_t L, int_fast64_t coordinate)
{
    int_fast64_t Lp = 2 * L + 3;
    int_fast64_t x = fmod(coordinate, Lp);
    int_fast64_t y = (fmod(coordinate, Lp * Lp) - x) / Lp;
    int_fast64_t z = (coordinate - x - y * Lp) / (Lp * Lp);
    return vint{ x, y, z };
}

//Generate the lattice structure based on size L
struct lattice lattice_fun(int_fast64_t L)
{
    struct lattice lat;
    //Iterate through lattice coordinates
    for (int_fast64_t x = 1; x < 2 * L + 2; ++x)
    {
        for (int_fast64_t y = 1; y < 2 * L + 2; ++y)
        {
            for (int_fast64_t z = 0; z < 2 * L + 3; ++z) //0 and 2L+2 for the e boundaries 
            {
                vint vec = { x, y, z };
                //Add vertices where all coordinates are odd
                if (fmod(x, 2) == 1 && fmod(y, 2) == 1 && fmod(z, 2) == 1)
                    lat.vertices.push_back(coordinate(L, vec));
                //Faces are defined using parity of coordinates below
                if ((fmod(x, 2) == 0 && fmod(y, 2) == 0 && fmod(z, 2) == 1) ||
                    (fmod(x, 2) == 0 && fmod(y, 2) == 1 && fmod(z, 2) == 0) ||
                    (fmod(x, 2) == 1 && fmod(y, 2) == 0 && fmod(z, 2) == 0))
                    lat.faces.push_back(coordinate(L, vec));
                //Add edges and coordinates of the logical operator 
                if ((fmod(x, 2) == 1 && fmod(y, 2) == 1 && fmod(z, 2) == 0) ||
                    (fmod(x, 2) == 1 && fmod(y, 2) == 0 && fmod(z, 2) == 1) ||
                    (fmod(x, 2) == 0 && fmod(y, 2) == 1 && fmod(z, 2) == 1))
                {
                    lat.edges.push_back(coordinate(L, vec));
                    if (x == 1 && y == 1) //Mark logical coordinates
                        lat.logicalcoords.push_back(coordinate(L, vec));
                }
                //Cubes are defined by all coordinates being even
                if (fmod(x, 2) == 0 && fmod(y, 2) == 0 && fmod(z, 2) == 0)
                    lat.cubes.push_back(coordinate(L, vec));
            }
        }
    }
    return lat;
}

//Check if a coordinate is inside the lattice boundaries
bool coord_inside_lattice(int_fast64_t L, int_fast64_t coord)
{
    vint xyzcoord = xyz(L, coord);
    int_fast64_t x = xyzcoord[0], y = xyzcoord[1], z = xyzcoord[2];
    return (x > 0 && x < 2 * L + 2 && y > 0 && y < 2 * L + 2 && z > -1 && z < 2 * L + 3);
}

//Find neighboring face stabilizers for a given qubit coordinate
vint neighb_facestabs_fun(int_fast64_t L, int_fast64_t qubit_coord, vint bulk_hole_coords)
{
    vint neighb_facestabs;
    vint xyzqubit_coord = xyz(L, qubit_coord);
    int_fast64_t x = xyzqubit_coord[0], y = xyzqubit_coord[1], z = xyzqubit_coord[2];
    int_fast64_t Lp = 2 * L + 3;

    //Calculate adjacent coordinates
    int_fast64_t sum_arre0 = qubit_coord + 1, sum_arre0m = qubit_coord - 1;
    int_fast64_t sum_arre1 = qubit_coord + Lp, sum_arre1m = qubit_coord - Lp;
    int_fast64_t sum_arre2 = qubit_coord + Lp * Lp, sum_arre2m = qubit_coord - Lp * Lp;

    //Check neighbors based on axis parity
    if (x % 2 == 0)
    {
        if (Type_notin_vType(bulk_hole_coords, sum_arre1) && coord_inside_lattice(L, sum_arre1))
            neighb_facestabs.push_back(sum_arre1);
        if (Type_notin_vType(bulk_hole_coords, sum_arre1m) && coord_inside_lattice(L, sum_arre1m))
            neighb_facestabs.push_back(sum_arre1m);
        if (Type_notin_vType(bulk_hole_coords, sum_arre2) && coord_inside_lattice(L, sum_arre2))
            neighb_facestabs.push_back(sum_arre2);
        if (Type_notin_vType(bulk_hole_coords, sum_arre2m) && coord_inside_lattice(L, sum_arre2m))
            neighb_facestabs.push_back(sum_arre2m);
    }
    if (y % 2 == 0)
    {
        if (Type_notin_vType(bulk_hole_coords, sum_arre0) && coord_inside_lattice(L, sum_arre0))
            neighb_facestabs.push_back(sum_arre0);
        if (Type_notin_vType(bulk_hole_coords, sum_arre0m) && coord_inside_lattice(L, sum_arre0m))
            neighb_facestabs.push_back(sum_arre0m);
        if (Type_notin_vType(bulk_hole_coords, sum_arre2) && coord_inside_lattice(L, sum_arre2))
            neighb_facestabs.push_back(sum_arre2);
        if (Type_notin_vType(bulk_hole_coords, sum_arre2m) && coord_inside_lattice(L, sum_arre2m))
            neighb_facestabs.push_back(sum_arre2m);
    }
    if (z % 2 == 0)
    {
        if (Type_notin_vType(bulk_hole_coords, sum_arre0) && coord_inside_lattice(L, sum_arre0))
            neighb_facestabs.push_back(sum_arre0);
        if (Type_notin_vType(bulk_hole_coords, sum_arre0m) && coord_inside_lattice(L, sum_arre0m))
            neighb_facestabs.push_back(sum_arre0m);
        if (Type_notin_vType(bulk_hole_coords, sum_arre1) && coord_inside_lattice(L, sum_arre1))
            neighb_facestabs.push_back(sum_arre1);
        if (Type_notin_vType(bulk_hole_coords, sum_arre1m) && coord_inside_lattice(L, sum_arre1m))
            neighb_facestabs.push_back(sum_arre1m);
    }
    return neighb_facestabs;
}

//Structure to hold fractal coordinates and related data
struct fractalcoords
{
    vvint neighb_facestabs;
    vint bulk_hole_coords, rm_sweep_indices, shelf_sweep_indices, face_stabs, qubit_coords, sweep_indices;
};

//Create fractal coordinates based on lattice, size, and level
struct fractalcoords create_fractal(struct lattice lat, int_fast64_t L, int_fast64_t level)
{
    struct fractalcoords frac_lattice;
    vvint neighb_stabs((2 * L + 3) * (2 * L + 3) * (2 * L + 3));
    vvint cornersv;
    frac_lattice.neighb_facestabs = neighb_stabs;

    if (level == 0)
    {
        frac_lattice.qubit_coords = lat.edges; //Initialize qubit coordinates
        frac_lattice.sweep_indices = lat.cubes; //Initialize sweep indices
        frac_lattice.face_stabs = lat.faces; //Initialize face stabilizers
        for (const int_fast64_t& qubit_coord : frac_lattice.qubit_coords)
            frac_lattice.neighb_facestabs[qubit_coord] = neighb_facestabs_fun(L, qubit_coord, frac_lattice.bulk_hole_coords);
        return frac_lattice;
    }
    else if (level == 1)
    {
        int_fast64_t Lv = L + 1;
        //Calculate corner coordinates for level 1
        cornersv.push_back({ int_fast64_t(round((float)Lv / 3.0)) + 1, int_fast64_t(round((float)Lv / 3.0)) + int_fast64_t(floor((float)Lv / 3.0)) });
    }
    else if (level == 2)
    {
        int_fast64_t Lv = L + 1;
        int_fast64_t x = int_fast64_t(round((float)Lv / 3.0));
        int_fast64_t y = Lv - int_fast64_t(round((float)Lv / 3.0)) - int_fast64_t(floor((float)Lv / 3.0));
        int_fast64_t iv = int_fast64_t(round((float)Lv / 3.0)) + int_fast64_t(floor((float)Lv / 3.0)) + int_fast64_t(round((float)y / 3.0)) + 1;

        vint cornerleft{ int_fast64_t(round((float)x / 3.0)) + 1, int_fast64_t(round((float)x / 3.0)) + int_fast64_t(floor((float)x / 3.0)) };
        vint cornermiddle{ int_fast64_t(round((float)Lv / 3.0)) + 1, int_fast64_t(round((float)Lv / 3.0)) + int_fast64_t(floor((float)Lv / 3.0)) };
        vint cornerright{ iv, iv + int_fast64_t(floor((float)y / 3.0)) - 1 };

        cornersv.insert(cornersv.end(), { cornerleft, cornermiddle, cornerright });
    }

    if (level != 0)
    {
        //Iterate through each corner to create bulk holes and sweep indices
        for (const vint& corner : cornersv)
        {
            int c1 = 2 * (corner[0] - 1), c2 = 2 * corner[1];
            auto holerange = range(c1, c2 + 1);

            for (const int_fast64_t& x : holerange)
            {
                for (const int_fast64_t& y : holerange)
                {
                    for (const int_fast64_t& z : holerange)
                    {
                        vint vec = { x, y, z };
                        frac_lattice.bulk_hole_coords.push_back(coordinate(L, vec)); //Add to bulk holes
                        //Remove sweep indices inside the hole
                        if ((fmod(x, 2) == 0 && fmod(y, 2) == 0 && fmod(z, 2) == 0) &&
                            (x > c1 && y > c1 && z > c2 && x < c2 && y < c2 && z < c2))
                            frac_lattice.rm_sweep_indices.push_back(coordinate(L, vec));
                        //Add sweep indices on the shelf
                        if ((fmod(x, 2) == 0 && fmod(y, 2) == 0 && fmod(z, 2) == 0) &&
                            (x == c1 || y == c1 || z == c2 || x == c2 || y == c2 || z == c2))
                            frac_lattice.shelf_sweep_indices.push_back(coordinate(L, vec));
                    }
                }
            }
        }
        //Update qubit coordinates excluding bulk holes
        for (const int_fast64_t& edge : lat.edges)
            if (Type_notin_vType(frac_lattice.bulk_hole_coords, edge))
                frac_lattice.qubit_coords.push_back(edge);

        //Update sweep indices excluding removed sweep indices
        for (const int_fast64_t& cube : lat.cubes)
            if (Type_notin_vType(frac_lattice.rm_sweep_indices, cube))
                frac_lattice.sweep_indices.push_back(cube);

        //Update face stabilizers excluding bulk holes
        for (const int_fast64_t& face : lat.faces)
            if (Type_notin_vType(frac_lattice.bulk_hole_coords, face))
                frac_lattice.face_stabs.push_back(face);

        //Assign neighboring face stabilizers
        for (const int_fast64_t& qubit_coord : frac_lattice.qubit_coords)
            frac_lattice.neighb_facestabs[qubit_coord] = neighb_facestabs_fun(L, qubit_coord, frac_lattice.bulk_hole_coords);
    }
    return frac_lattice;
}

//Find faces related to a sweep index and direction
vint faces_fun(int_fast64_t L, vint bulk_hole_coords, int_fast64_t sweep_index, vint sweep_dir, int_fast64_t pastorfuture)
{
    vint faces_porf;
    vint chg = vectimesc(pastorfuture, sweep_dir); //Calculate change based on direction
    int_fast64_t Lp = 2 * L + 3;

    //Calculate face coordinates
    int_fast64_t face0 = sweep_index + chg[0];
    int_fast64_t face1 = sweep_index + chg[1] * Lp;
    int_fast64_t face2 = sweep_index + chg[2] * Lp * Lp;

    //We only want to check if the coordinate is inside the lattice since we want to allow bulk hole future faces to be in the output
    //Check if faces are inside the lattice
    if (coord_inside_lattice(L, face0))
        faces_porf.push_back(face0);
    if (coord_inside_lattice(L, face1))
        faces_porf.push_back(face1);
    if (coord_inside_lattice(L, face2))
        faces_porf.push_back(face2);

    return faces_porf;
}

//Update syndrome for a single trail hole
vint onetrailhole_update(int_fast64_t L, vint shelf_sweep_indices, vint bulk_hole_coords, vint syndrome, vint sweep_dir)
{
    for (const int_fast64_t& sweep_index : shelf_sweep_indices)
    {
        vint syndromes_cubefaces_past, syndromes_cubefaces_ftr;
        vint faces_past = faces_fun(L, bulk_hole_coords, sweep_index, sweep_dir, -1);
        vint faces_ftr = faces_fun(L, bulk_hole_coords, sweep_index, sweep_dir, 1);
        
        //Collect syndromes for past and future faces
        for (const int_fast64_t& face_past : faces_past)
            syndromes_cubefaces_past.push_back(syndrome[face_past]);
        for (const int_fast64_t& face_ftr : faces_ftr)
            syndromes_cubefaces_ftr.push_back(syndrome[face_ftr]);

        //Check conditions to update syndrome
        if (!(one_in_vector(syndromes_cubefaces_past)) && countone_in_vector(syndromes_cubefaces_ftr) == 1)
        {
            vint zero_synd_faces_ftr;
            for (const int_fast64_t& face_ftr : faces_ftr)
            {
                if (Type_in_vType(bulk_hole_coords, face_ftr))
                {
                    //Ensure syndrome is zero for bulk hole faces
                    zero_synd_faces_ftr.push_back(face_ftr);
                }
            }
            //Update syndrome based on the number of zero syndromes
            if (zero_synd_faces_ftr.size() == 2)
            {
                syndrome[zero_synd_faces_ftr[distInt0To1(rnEngine)]] = 1;
            }
            else if (zero_synd_faces_ftr.size() == 1)
            {
                syndrome[zero_synd_faces_ftr[0]] = 1;
            }
        }
    }
    return syndrome;
}

//Structure to hold syndrome and error information
struct syndromeerror
{
    vint syndrome, error;
};

//Check trailing condition for a sweep index
bool check_trailing(int_fast64_t L, vint bulk_hole_coords, int_fast64_t sweep_index, vint syndrome, vint sweep_dir)
{
    vint syndromes_cubefaces_past, syndromes_cubefaces_ftr;
    vint faces_past = faces_fun(L, bulk_hole_coords, sweep_index, sweep_dir, -1);
    vint faces_ftr = faces_fun(L, bulk_hole_coords, sweep_index, sweep_dir, 1);
    
    //Collect syndromes for past and future faces
    for (const int_fast64_t& face_past : faces_past)
        syndromes_cubefaces_past.push_back(syndrome[face_past]);
    for (const int_fast64_t& face_ftr : faces_ftr)
        syndromes_cubefaces_ftr.push_back(syndrome[face_ftr]);

    //Determine if the trailing condition is met
    if (!(one_in_vector(syndromes_cubefaces_past)) && (countone_in_vector(syndromes_cubefaces_ftr) > 1))
        return true;
    return false;
}

//Update syndrome based on edge flips
vint syndrome_updt_fun(vint syndrome, sint edgeflips, vvint neighb_facestabs, vint bulk_hole_coords)
{
    for (const int_fast64_t& qubit_coord : edgeflips)
    {
        vint neighb_stabs = neighb_facestabs[qubit_coord];
        for (const int_fast64_t& neighb_stab : neighb_stabs)
        {
            syndrome[neighb_stab] = fmod(syndrome[neighb_stab] + 1, 2);
        }
    }
    //The above loop doesn't update the syndromes on the faces inside the hole since they were made nonzero in the onetrailhole_update, 
    //We need to put them back to 0
    //Reset syndromes inside bulk holes
    for (const int_fast64_t& bulk_hole_coord : bulk_hole_coords)
    {
        syndrome[bulk_hole_coord] = 0;
    }
    return syndrome;
}

//Perform a single sweep step
struct syndromeerror sweep_step(int_fast64_t L, int_fast64_t level, struct lattice lat, struct fractalcoords frac_lattice, struct syndromeerror synderror, vint sweep_dir, vint meas_err_faces)
{
    int_fast64_t Lp = 2 * L + 3;
    //Update syndrome for level > 0
    if (level != 0)
        synderror.syndrome = onetrailhole_update(L, frac_lattice.shelf_sweep_indices, frac_lattice.bulk_hole_coords, synderror.syndrome, sweep_dir);

    sint edgeflips;
    //Iterate through sweep indices
    for (const int_fast64_t& sweep_index : frac_lattice.sweep_indices)
    {
        if (check_trailing(L, frac_lattice.bulk_hole_coords, sweep_index, synderror.syndrome, sweep_dir))
        {
            vint faces_ftr = faces_fun(L, frac_lattice.bulk_hole_coords, sweep_index, sweep_dir, 1);
            vint nonzerosyndromes_cubefaces_ftr;
            //Collect non-zero syndromes in future faces
            for (const int_fast64_t& face_ftr : faces_ftr)
            {
                if (synderror.syndrome[face_ftr] == 1)
                    nonzerosyndromes_cubefaces_ftr.push_back(face_ftr);
            }
            //Adjust syndromes if there are three non-zero syndromes
            if (nonzerosyndromes_cubefaces_ftr.size() == 3)
            {
                nonzerosyndromes_cubefaces_ftr.erase(nonzerosyndromes_cubefaces_ftr.begin() + distInt0To2(rnEngine));
            }
            //Calculate updated edge coordinate
            int_fast64_t edge_updtd = nonzerosyndromes_cubefaces_ftr[0] + nonzerosyndromes_cubefaces_ftr[1] - sweep_index;
            //Update error and record edge flip
            synderror.error[edge_updtd] = fmod(synderror.error[edge_updtd] + 1, 2);
            updateset(edgeflips, edge_updtd);

            //we can't update the syndrome here as that changes the set of trailing indices with just the syndrome on one sweep index being updated
            //but we first want to update all the syndromes in the future of all trailing sweep indices
            //because we are not doing in parallel like you'd imagine doing in an actual physical device 
        }
    }
    //Update syndrome after processing all edge flips
    synderror.syndrome = syndrome_updt_fun(synderror.syndrome, edgeflips, frac_lattice.neighb_facestabs, frac_lattice.bulk_hole_coords);
    return synderror;
}

//Generate data errors based on qubit coordinates and error probabilities
struct syndromeerror gen_data_err(int_fast64_t L, vint qubit_coords, vvint neighb_facestabs, const double derr_prob, struct syndromeerror synderror)
{
    vint veczero((2 * L + 3) * (2 * L + 3) * (2 * L + 3));
    //Initialize synderror if needed
    for (const int_fast64_t& qubit_coord : qubit_coords)
    {
        if (distDouble0To1(rnEngine) <= derr_prob)
        {
            synderror.error[qubit_coord] = fmod(synderror.error[qubit_coord] + 1, 2); //Flip error bit
        }
        if (synderror.error[qubit_coord] == 1)
        {
            vint neighb_stabs = neighb_facestabs[qubit_coord];
            for (const int_fast64_t& neighb_face_stab : neighb_stabs)
            {
                synderror.syndrome[neighb_face_stab] = fmod(syndrome[neighb_face_stab] + 1, 2); //Update syndrome
            }
        }
    }
    return synderror;
}

//Generate measurement errors on faces based on probability
vint measerr_faces_fun(vint syndrome, vint face_stabs, const double merr_prob)
{
    vint measerr_faces;
    for (const int_fast64_t& face_stab : face_stabs)
    {
        if (distDouble0To1(rnEngine) <= merr_prob)
        {
            measerr_faces.push_back(face_stab); //Record measurement error
        }
    }
    return measerr_faces;
}

//Main decoding run function
int_fast64_t sweep_decoder_run(int_fast64_t level, const std::string& sweep_schedule, int_fast64_t rounds, int_fast64_t L, const double derr_prob, const double merr_prob, pcg32 rnEngine)
{
    auto start = std::chrono::high_resolution_clock::now(); //Start timer

    int_fast64_t timeout = 32 * L; //Set timeout limit
    int_fast64_t same_sweep_dir_rdlimit = int_fast64_t(round(log(L))); //Limit for same sweep direction
    int_fast64_t sweeprate = 1; //Sweep rate

    struct lattice lat = lattice_fun(L); //Generate lattice
    struct fractalcoords frac_lattice = create_fractal(lat, L, level); //Create fractal coordinates

    struct syndromeerror synderror;
    vint veczero((2 * L + 3) * (2 * L + 3) * (2 * L + 3));
    synderror.error = veczero; //Initialize errors
    synderror.syndrome = veczero; //Initialize syndromes 

    int_fast64_t same_sweep_dir_rdcount = 0; //Counter for sweep direction changes

    vint sweep_dir;
    //Define possible sweep directions
    vvint dir_list{ {1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1}, {-1,-1,-1}, {-1,1,1}, {1,-1,1}, {1,1,-1} };
    int_fast64_t dir_count = 0;
    sweep_dir = dir_list[dir_count % 8]; //Initialize sweep direction

    //Iterate through decoding rounds
    for (const int_fast64_t& rd_index : range(rounds))
    {
        if (same_sweep_dir_rdcount == same_sweep_dir_rdlimit)
        {
            if (sweep_schedule == "alternating")
            {
                dir_count += 1;
                sweep_dir = dir_list[dir_count % 8]; //Change direction
            }
            else if (sweep_schedule == "constant")
            {
                sweep_dir = { 1, 1, 1 }; //Keep constant direction
            }
            same_sweep_dir_rdcount = 0; //Reset counter
        }
        synderror.syndrome = veczero; //Reset syndrome //we recalculate the syndrome below
        synderror = gen_data_err(L, frac_lattice.qubit_coords, frac_lattice.neighb_facestabs, derr_prob, synderror); //Generate data errors

        vint measerr_faces;
        if (merr_prob > 0)
        {
            measerr_faces = measerr_faces_fun(synderror.syndrome, frac_lattice.face_stabs, merr_prob); //Generate measurement errors
            for (const int_fast64_t& face : measerr_faces)
                synderror.syndrome[face] = fmod(syndrome[face] + 1, 2); //Update syndrome with measurement errors
        }

        for (int_fast64_t i = 0; i < sweeprate; ++i)
            synderror = sweep_step(L, level, lat, frac_lattice, synderror, sweep_dir, meas_err_faces); //Perform sweep step

        same_sweep_dir_rdcount += 1; //Increment sweep direction counter
    }

    //Final error generation and sweep
    synderror.syndrome = veczero; 
    vint measerr_faces;
    synderror = gen_data_err(L, frac_lattice.qubit_coords, frac_lattice.neighb_facestabs, derr_prob, synderror);

    int_fast64_t same_sweep_dir_limit_timeout_session = L;
    int_fast64_t same_sweep_dir_time = 0;
    sweep_dir = { 1, 1, 1 }; //Reset sweep direction

    for (const int_fast64_t& t : range(timeout))
    {
        if (same_sweep_dir_time == same_sweep_dir_limit_timeout_session)
        {
            if (sweep_schedule == "alternating")
            {
                dir_count += 1;
                sweep_dir = dir_list[dir_count % 8]; //Change direction
            }
            else if (sweep_schedule == "constant")
            {
                sweep_dir = { 1, 1, 1 }; //Keep constant direction
            }
            same_sweep_dir_time = 0; //Reset counter
        }
        synderror = sweep_step(L, level, lat, frac_lattice, synderror, sweep_dir, measerr_faces); //Perform sweep step
        same_sweep_dir_time += 1; //Increment sweep direction timer

        //Exit if syndrome is clean
        if (!(one_in_vector(synderror.syndrome)))
            break;
    }

    //Determine logical failure
    int_fast64_t logical_failure = mod2dot_fun(lat.logicalcoords, synderror.error);
    bool syndrome_unclean = one_in_vector(synderror.syndrome);

    auto finish = std::chrono::high_resolution_clock::now(); //End timer
    std::chrono::duration<double> diff = finish - start;
    //std::cout << syndrome_unclean << "," << logical_failure << "," << diff.count() << std::endl; //Debug output
    return (syndrome_unclean || logical_failure); //Return failure status
}

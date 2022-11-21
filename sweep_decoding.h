#include "preamble/preamble.h"

struct lattice
{
    vint vertices,edges,faces,cubes;
    vint logicalcoords;
};

int_fast64_t coordinate(int_fast64_t L,vint vec)
{
    int_fast64_t x=vec[0],y=vec[1],z=vec[2];
    int_fast64_t Lp=2*L+3;
    return x+y*Lp+z*Lp*Lp;
}

vint xyz(int_fast64_t L,int_fast64_t coordinate)
{
    int_fast64_t Lp=2*L+3;
    int_fast64_t x=fmod(coordinate,Lp), y=(fmod(coordinate,Lp*Lp)-x)/Lp, z=(coordinate-x-y*Lp)/(Lp*Lp);
    vint xyz_vec{x,y,z};
    return xyz_vec;
}

struct lattice lattice_fun(int_fast64_t L)
{
    //L cubes, L_vertices=L+1
    struct lattice lat;
    for (int_fast64_t x=1; x<2*L+2;++x)
    {
        for (int_fast64_t y=1; y<2*L+2;++y)
        {
            for (int_fast64_t z=0; z<2*L+3;++z) //0 and 2L+2 for the e boundaries 
            {
                vint vec{x,y,z};
                if (fmod(x,2)==1 && fmod(y,2)==1 && fmod(z,2)==1)
                    lat.vertices.push_back(coordinate(L,vec));
                if ((fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==1) || (fmod(x,2)==0 && fmod(y,2)==1 && fmod(z,2)==0) || (fmod(x,2)==1 && fmod(y,2)==0 && fmod(z,2)==0))
                    lat.faces.push_back(coordinate(L,vec));
                if ((fmod(x,2)==1 && fmod(y,2)==1 && fmod(z,2)==0) || (fmod(x,2)==1 && fmod(y,2)==0 && fmod(z,2)==1) || (fmod(x,2)==0 && fmod(y,2)==1 && fmod(z,2)==1))
                {
                    lat.edges.push_back(coordinate(L,vec));
                    if (x==1 && y==1)
                        lat.logicalcoords.push_back(coordinate(L,vec));                    
                }
                if (fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==0)
                    lat.cubes.push_back(coordinate(L,vec));
            }
        }
    }
    return lat;
}

bool coord_inside_lattice(int_fast64_t L, int_fast64_t coord)
{
    vint xyzcoord=xyz(L,coord);
    int_fast64_t x = xyzcoord[0],y=xyzcoord[1],z=xyzcoord[2];
    return (x>0 && x<2*L+2 && y>0 && y<2*L+2 && z>-1 && z<2*L+3);
}

vint neighb_facestabs_fun(int_fast64_t L,int_fast64_t qubit_coord,vint bulk_hole_coords)
{
    vint neighb_facestabs;
    vint xyzqubit_coord=xyz(L,qubit_coord);
    int_fast64_t x = xyzqubit_coord[0],y=xyzqubit_coord[1],z=xyzqubit_coord[2];

    int_fast64_t Lp=2*L+3;
    int_fast64_t sum_arre0=qubit_coord+1, sum_arre0m=qubit_coord-1;
    int_fast64_t sum_arre1=qubit_coord+Lp, sum_arre1m=qubit_coord-Lp;
    int_fast64_t sum_arre2=qubit_coord+Lp*Lp, sum_arre2m=qubit_coord-Lp*Lp;

    if (x%2==0)
    {
        if (Type_notin_vType(bulk_hole_coords,sum_arre1) && coord_inside_lattice(L,sum_arre1))
            neighb_facestabs.push_back(sum_arre1);
        if (Type_notin_vType(bulk_hole_coords,sum_arre1m) && coord_inside_lattice(L,sum_arre1m))
            neighb_facestabs.push_back(sum_arre1m);
        if (Type_notin_vType(bulk_hole_coords,sum_arre2) && coord_inside_lattice(L,sum_arre2))
            neighb_facestabs.push_back(sum_arre2);
        if (Type_notin_vType(bulk_hole_coords,sum_arre2m) && coord_inside_lattice(L,sum_arre2m))
            neighb_facestabs.push_back(sum_arre2m);
    }
    if (y%2==0)
    {
        if (Type_notin_vType(bulk_hole_coords,sum_arre0) && coord_inside_lattice(L,sum_arre0))
            neighb_facestabs.push_back(sum_arre0);
        if (Type_notin_vType(bulk_hole_coords,sum_arre0m) && coord_inside_lattice(L,sum_arre0m))
            neighb_facestabs.push_back(sum_arre0m);
        if (Type_notin_vType(bulk_hole_coords,sum_arre2) && coord_inside_lattice(L,sum_arre2))
            neighb_facestabs.push_back(sum_arre2);
        if (Type_notin_vType(bulk_hole_coords,sum_arre2m) && coord_inside_lattice(L,sum_arre2m))
            neighb_facestabs.push_back(sum_arre2m);
    }
    if (z%2==0)
    {
        if (Type_notin_vType(bulk_hole_coords,sum_arre0) && coord_inside_lattice(L,sum_arre0))
            neighb_facestabs.push_back(sum_arre0);
        if (Type_notin_vType(bulk_hole_coords,sum_arre0m) && coord_inside_lattice(L,sum_arre0m))
            neighb_facestabs.push_back(sum_arre0m);
        if (Type_notin_vType(bulk_hole_coords,sum_arre1) && coord_inside_lattice(L,sum_arre1))
            neighb_facestabs.push_back(sum_arre1);
        if (Type_notin_vType(bulk_hole_coords,sum_arre1m) && coord_inside_lattice(L,sum_arre1m))
            neighb_facestabs.push_back(sum_arre1m);
    }
    return neighb_facestabs;
}

struct fractalcoords
{
  vvint neighb_facestabs;
  vint bulk_hole_coords,rm_sweep_indices,shelf_sweep_indices,face_stabs,qubit_coords,sweep_indices;
};

struct fractalcoords create_fractal(struct lattice lat,int_fast64_t L,int_fast64_t level)
{
    struct fractalcoords frac_lattice; 
    vvint neighb_stabs((2*L+3)*(2*L+3)*(2*L+3));
    vvint cornersv;
    frac_lattice.neighb_facestabs=neighb_stabs;
    if(level==0)
    {
        frac_lattice.qubit_coords=lat.edges;
        frac_lattice.sweep_indices=lat.cubes;
        frac_lattice.face_stabs=lat.faces;
        for (const int_fast64_t& qubit_coord:frac_lattice.qubit_coords)
            frac_lattice.neighb_facestabs[qubit_coord]=neighb_facestabs_fun(L,qubit_coord,frac_lattice.bulk_hole_coords);
        return frac_lattice;
    }
    else if (level==1)
    {
        int_fast64_t Lv=L+1;
        cornersv.push_back({int_fast64_t(round((float)Lv/(float)3))+1,int_fast64_t(round((float)Lv/(float)3))+int_fast64_t(floor((float)Lv/(float)3))});
    }
    else if (level==2)
    {
        int_fast64_t Lv=L+1;
        int_fast64_t x=int_fast64_t(round((float)Lv/(float)3));
        int_fast64_t y=Lv-int_fast64_t(round((float)Lv/(float)3))-int_fast64_t(floor((float)Lv/(float)3));
        int_fast64_t iv=int_fast64_t(round((float)Lv/(float)3))+int_fast64_t(floor((float)Lv/(float)3))+int_fast64_t(round((float)y/(float)3))+1;

        vint cornerleft{int_fast64_t(round((float)x/(float)3))+1,int_fast64_t(round((float)x/(float)3))+int_fast64_t(floor((float)x/(float)3))};
        vint cornermiddle{int_fast64_t(round((float)Lv/(float)3))+1,int_fast64_t(round((float)Lv/(float)3))+int_fast64_t(floor((float)Lv/(float)3))};
        vint cornerright{iv,iv+int_fast64_t(floor((float)y/(float)3))-1};

        cornersv.insert(cornersv.end(),{cornerleft,cornermiddle,cornerright});
    }

    if (level!=0)
    {
        for(const vint& corner:cornersv)
        {   
            int c1= 2*(corner[0]-1), c2=2*corner[1];
            auto holerange=range(c1,c2+1); 

            for(const int_fast64_t& x:holerange)
            {
                for(const int_fast64_t& y:holerange) 
                {
                    for(const int_fast64_t& z:holerange)
                    {
                        vint vec{x,y,z};
                        frac_lattice.bulk_hole_coords.push_back(coordinate(L,vec));
                        if ((fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==0) && (x>c1 && y>c1 && z>c2 && x<c2 && y<c2 && z<c2))
                            frac_lattice.rm_sweep_indices.push_back(coordinate(L,vec));
                        if ((fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==0) && (x==c1|| y==c1|| z==c2 || x==c2 || y==c2 || z==c2))
                            frac_lattice.shelf_sweep_indices.push_back(coordinate(L,vec)); 
                    }
                }
            }
        }
        for (const int_fast64_t& edge:lat.edges)
            if (Type_notin_vType(frac_lattice.bulk_hole_coords,edge))
                frac_lattice.qubit_coords.push_back(edge);

        for (const int_fast64_t& cube:lat.cubes)
            if (Type_notin_vType(frac_lattice.rm_sweep_indices,cube))
                frac_lattice.sweep_indices.push_back(cube);

        for (const int_fast64_t& face:lat.faces)
            if (Type_notin_vType(frac_lattice.bulk_hole_coords,face))
                frac_lattice.face_stabs.push_back(face);         
                
        for (const int_fast64_t& qubit_coord:frac_lattice.qubit_coords)
            frac_lattice.neighb_facestabs[qubit_coord]=neighb_facestabs_fun(L,qubit_coord,frac_lattice.bulk_hole_coords);
    }
    return frac_lattice;
}

vint faces_fun(int_fast64_t L,vint bulk_hole_coords,int_fast64_t sweep_index,vint sweep_dir,int_fast64_t pastorfuture)
{
    vint faces_porf;

    vint chg=vectimesc(pastorfuture,sweep_dir);
    int_fast64_t Lp=2*L+3;

    int_fast64_t face0,face1,face2;

    face0=sweep_index+chg[0];
    face1=sweep_index+chg[1]*Lp;
    face2=sweep_index+chg[2]*Lp*Lp;

    //we only want to check if the coordinate is inside the lattice since we want to allow bulkhole future faces to be in the output
    if (coord_inside_lattice(L,face0))
        faces_porf.push_back(face0);
    if (coord_inside_lattice(L,face1))
        faces_porf.push_back(face1);
    if (coord_inside_lattice(L,face2))
        faces_porf.push_back(face2);

    return faces_porf;
}

vint onetrailhole_update(int_fast64_t L,vint shelf_sweep_indices,vint bulk_hole_coords,vint syndrome,vint sweep_dir)
{
    for(const int_fast64_t& sweep_index:shelf_sweep_indices)
    {
        vint syndromes_cubefaces_past,syndromes_cubefaces_ftr;
        vint faces_past=faces_fun(L,bulk_hole_coords,sweep_index,sweep_dir,-1);
        vint faces_ftr=faces_fun(L,bulk_hole_coords,sweep_index,sweep_dir,1);
        for (const int_fast64_t& face_past:faces_past)
            syndromes_cubefaces_past.push_back(syndrome[face_past]);
        for (const int_fast64_t& face_ftr:faces_ftr)
            syndromes_cubefaces_ftr.push_back(syndrome[face_ftr]);
        if (!(one_in_vector(syndromes_cubefaces_past)) && countone_in_vector(syndromes_cubefaces_ftr)==1)
        {
            vint zero_synd_faces_ftr;
            for(const int_fast64_t& face_ftr:faces_ftr)
            {
                if (Type_in_vType(bulk_hole_coords,face_ftr))
                {
                    // assert(syndrome[face_ftr]==0);
                    zero_synd_faces_ftr.push_back(face_ftr);
                }
            }
            if (zero_synd_faces_ftr.size()==2) 
            {
                syndrome[zero_synd_faces_ftr[distInt0To1(rnEngine)]]=1;
            }
            else if (zero_synd_faces_ftr.size()==1) 
            {            
                syndrome[zero_synd_faces_ftr[0]]=1;
            }
        }
    }
    return syndrome;
}

struct syndromeerror
{
 vint syndrome,error;   
};

bool check_trailing(int_fast64_t L,vint bulk_hole_coords,int_fast64_t sweep_index,vint syndrome,vint sweep_dir)
{
    vint syndromes_cubefaces_past,syndromes_cubefaces_ftr;
    vint faces_past=faces_fun(L,bulk_hole_coords,sweep_index,sweep_dir,-1);
    vint faces_ftr=faces_fun(L,bulk_hole_coords,sweep_index,sweep_dir,1);
    for (const int_fast64_t& face_past:faces_past)
        syndromes_cubefaces_past.push_back(syndrome[face_past]);
    for (const int_fast64_t& face_ftr:faces_ftr)
        syndromes_cubefaces_ftr.push_back(syndrome[face_ftr]);
    if(!(one_in_vector(syndromes_cubefaces_past)) && (countone_in_vector(syndromes_cubefaces_ftr)>1))
        return true;
    return false;
}

vint syndrome_updt_fun(vint syndrome,sint edgeflips,vvint neighb_facestabs,vint bulk_hole_coords)
{
    for (const int_fast64_t& qubit_coord:edgeflips)
    {
        vint neighb_stabs=neighb_facestabs[qubit_coord];
        for (const int_fast64_t& neighb_stab:neighb_stabs)
        {
            syndrome[neighb_stab]=fmod(syndrome[neighb_stab]+1,2);
        }
    }
    //the above loop doesn't update the syndromes on the faces inside the hole, since they were made nonzero in the onetrailhole_update, we need to put them back to 0
    for(const int_fast64_t& bulk_hole_coord:bulk_hole_coords)
    {
        syndrome[bulk_hole_coord]=0;
    }
    return syndrome;
}


struct syndromeerror sweep_step(int_fast64_t L, int_fast64_t level, struct lattice lat,struct fractalcoords frac_lattice,struct syndromeerror synderror,vint sweep_dir,vint meas_err_faces)
{
    int_fast64_t Lp=2*L+3;
    if (level!=0)
        synderror.syndrome=onetrailhole_update(L,frac_lattice.shelf_sweep_indices,frac_lattice.bulk_hole_coords,synderror.syndrome,sweep_dir);

    sint edgeflips;
    for(const int_fast64_t& sweep_index:frac_lattice.sweep_indices)
    {
        if (check_trailing(L,frac_lattice.bulk_hole_coords,sweep_index,synderror.syndrome,sweep_dir))
        {
            vint faces_ftr=faces_fun(L,frac_lattice.bulk_hole_coords,sweep_index,sweep_dir,1);
            vint nonzerosyndromes_cubefaces_ftr;
            for(const int_fast64_t& face_ftr:faces_ftr)
            {
                if (synderror.syndrome[face_ftr]==1)
                    nonzerosyndromes_cubefaces_ftr.push_back(face_ftr);
            }
            // assert(nonzerosyndromes_cubefaces_ftr.size()>1);
            if (nonzerosyndromes_cubefaces_ftr.size()==3)
            {
                nonzerosyndromes_cubefaces_ftr.erase(nonzerosyndromes_cubefaces_ftr.begin()+distInt0To2(rnEngine));
            }
            int_fast64_t edge_updtd=nonzerosyndromes_cubefaces_ftr[0]+nonzerosyndromes_cubefaces_ftr[1]-sweep_index;
            // assert(Type_in_vType(frac_lattice.qubit_coords,edge_updtd));
            synderror.error[edge_updtd]=fmod(synderror.error[edge_updtd]+1,2);
            updateset(edgeflips,edge_updtd);

            //we can't update the syndrome here as that changes the set of trailing indices with just the syndrome on one sweep index being updated
            // but we first want to update all the syndromes in the future of all trailing sweep indices
            //because we are not doing in parallel like you'd imagine doing in an actual physical device 
        }
    }
    synderror.syndrome=syndrome_updt_fun(synderror.syndrome,edgeflips,frac_lattice.neighb_facestabs,frac_lattice.bulk_hole_coords);
    return synderror;
}

struct syndromeerror gen_data_err(int_fast64_t L,vint qubit_coords,vvint neighb_facestabs,const double derr_prob,struct syndromeerror synderror)
{
    vint veczero((2*L+3)*(2*L+3)*(2*L+3));
    // assert(synderror.syndrome==veczero);
    for (const int_fast64_t& qubit_coord:qubit_coords)
    {
        if (distDouble0To1(rnEngine)<= derr_prob)
        {
            synderror.error[qubit_coord]=fmod(synderror.error[qubit_coord]+1,2);
        }
        if (synderror.error[qubit_coord]==1)
        {
            vint neighb_stabs=neighb_facestabs[qubit_coord];
            for (const int_fast64_t& neighb_face_stab:neighb_stabs)
            {
                synderror.syndrome[neighb_face_stab]=fmod(synderror.syndrome[neighb_face_stab]+1,2);
            }
        }
    }
    return synderror;
}

vint measerr_faces_fun(vint syndrome,vint face_stabs,const double merr_prob)
{
    vint measerr_faces;
    for (const int_fast64_t& face_stab:face_stabs)
    {
        if (distDouble0To1(rnEngine)<= merr_prob)
        {   
            measerr_faces.push_back(face_stab);
        }
    }
    return measerr_faces;
}

int_fast64_t sweep_decoder_run(int_fast64_t level,const std::string &sweep_schedule,int_fast64_t rounds,int_fast64_t L,const double derr_prob,const double merr_prob,pcg32 rnEngine)
{
    
    auto start = std::chrono::high_resolution_clock::now();

    int_fast64_t timeout=32*L;
    int_fast64_t same_sweep_dir_rdlimit=int_fast64_t(round(log(L)));
    int_fast64_t sweeprate=1;

    struct lattice lat=lattice_fun(L);
    
    struct fractalcoords frac_lattice=create_fractal(lat,L,level);

    struct syndromeerror synderror;
    vint veczero((2*L+3)*(2*L+3)*(2*L+3));
    synderror.error=veczero;
    synderror.syndrome=veczero;

    int_fast64_t same_sweep_dir_rdcount=0;

    vint sweep_dir;
    vvint dir_list{{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1},{-1,-1,-1},{-1,1,1},{1,-1,1},{1,1,-1}};
    int_fast64_t dir_count=0;
    sweep_dir=dir_list[dir_count%8];

    for (const int_fast64_t& rd_index:range(rounds))
    {
        if(same_sweep_dir_rdcount==same_sweep_dir_rdlimit)
        {
            if (sweep_schedule=="alternating")
            {                
                dir_count+=1;
                sweep_dir=dir_list[dir_count%8];
            }
            else if (sweep_schedule=="constant")
            {
                sweep_dir={1,1,1};
            }
            same_sweep_dir_rdcount=0;
        }
        synderror.syndrome=veczero; // we recalculate the syndrome below
        synderror=gen_data_err(L,frac_lattice.qubit_coords,frac_lattice.neighb_facestabs,derr_prob,synderror);

        vint measerr_faces;
        if (merr_prob>0)
        {
            measerr_faces=measerr_faces_fun(synderror.syndrome,frac_lattice.face_stabs,merr_prob);    
            for (const int_fast64_t& face: measerr_faces)
                synderror.syndrome[face]=fmod(synderror.syndrome[face]+1,2);
        }

        for(int_fast64_t i=0; i<sweeprate;++i)
            synderror=sweep_step(L,level,lat,frac_lattice,synderror,sweep_dir,measerr_faces);
        
        same_sweep_dir_rdcount+=1;        
    }
    
    synderror.syndrome=veczero;
    vint measerr_faces;
    synderror=gen_data_err(L,frac_lattice.qubit_coords,frac_lattice.neighb_facestabs,derr_prob,synderror);

    int_fast64_t same_sweep_dir_limit_timeout_session=L;

    int_fast64_t same_sweep_dir_time=0;
    sweep_dir={1,1,1};

    for (const int_fast64_t& t:range(timeout))
    {
        if(same_sweep_dir_time==same_sweep_dir_limit_timeout_session)
        {
            if (sweep_schedule=="alternating")
            {
                dir_count+=1;
                sweep_dir=dir_list[dir_count%8];
            }
            else if (sweep_schedule=="constant")
            {
                sweep_dir={1,1,1};
            }
            same_sweep_dir_time=0;
        }
        synderror=sweep_step(L,level,lat,frac_lattice,synderror,sweep_dir,measerr_faces);

        same_sweep_dir_time+=1;

        if (!(one_in_vector(synderror.syndrome)))
            break;
    }

    int_fast64_t logical_failure = mod2dot_fun(lat.logicalcoords,synderror.error);
    bool syndrome_unclean= one_in_vector(synderror.syndrome);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = finish-start;
    // std::cout<<syndrome_unclean<<","<<logical_failure<<","<<diff.count()<<std::endl;
    return (syndrome_unclean || logical_failure);
}


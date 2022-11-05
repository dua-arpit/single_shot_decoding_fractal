#include "preamble/preamble.h"

struct lattice
{
    vint vertices,edges,faces,cubes;/*  */
    vint logical;
};

int coordinate(int L,vint vec)
{
    int x=vec[0],y=vec[1],z=vec[2];
    int Lp=(2*L+3);
    return x+y*Lp+z*Lp*Lp;
}

vint xyz(int L,int coordinate)
{
    int Lp=(2*L+3);
    int x=fmod(coordinate,Lp), y=(fmod(coordinate,(Lp*Lp))-x)/Lp, z=(coordinate-x-y*Lp)/(Lp*Lp);
    vint xyz_vec{x,y,z};
    return xyz_vec;
}

struct lattice lattice_fun(int L)
{
    struct lattice lat;
    vint logical((2*L+3)*(2*L+3)*(2*L+3));
    lat.logical=logical;
    for (int x=1; x<2*L+2;++x)
    {
        for (int y=1; y<2*L+2;++y)
        {
            for (int z=0; z<2*L+3;++z) //0 and 2L+2 for the e boundaries 
            {
                vint vec{x,y,z};
                if (fmod(x,2)==1 && fmod(y,2)==1 && fmod(z,2)==1)
                    lat.vertices.push_back(coordinate(L,vec));
                if ((fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==1) || (fmod(x,2)==0 && fmod(y,2)==1 && fmod(z,2)==0) || (fmod(x,2)==1 && fmod(y,2)==0 && fmod(z,2)==0))
                    lat.faces.push_back(coordinate(L,vec));
                if ((fmod(x,2)==1 && fmod(y,2)==1 && fmod(z,2)==0) || (fmod(x,2)==1 && fmod(y,2)==0 && fmod(z,2)==1) || (fmod(x,2)==0 && fmod(y,2)==1 && fmod(z,2)==1))
                {
                    lat.edges.push_back(coordinate(L,vec));
                    if (x==1 && y==1 && fmod(z,2)==0)
                        lat.logical[coordinate(L,vec)]=1;                    
                }
                if (fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==0)
                    lat.cubes.push_back(coordinate(L,vec));
            }
        }
    }
    return lat;
}

vint neighb_facestabs_fun(int L,int qubit_coord,vint face_stabs)
{
    vint neighb_facestabs;
    vint e0{1,0,0},e1{0,1,0},e2{0,0,1},e0m{-1,0,0},e1m{0,-1,0},e2m{0,0,-1};

    vint xyzqubit_coord=xyz(L,qubit_coord);
    int Lp=2*L+3;
    int coordmaxplus=Lp*Lp*Lp;

    int sum_arre0=qubit_coord+1, sum_arre0m=qubit_coord-1;
    int sum_arre1=qubit_coord+Lp, sum_arre1m=qubit_coord-Lp;
    int sum_arre2=qubit_coord+Lp*Lp, sum_arre2m=qubit_coord-Lp*Lp;

    if (xyzqubit_coord[0]%2==0)
    {
        if (Type_in_vType(face_stabs,sum_arre1))
            neighb_facestabs.push_back(sum_arre1);
        if (Type_in_vType(face_stabs,sum_arre1m))
            neighb_facestabs.push_back(sum_arre1m);
        if (Type_in_vType(face_stabs,sum_arre2))
            neighb_facestabs.push_back(sum_arre2);
        if (Type_in_vType(face_stabs,sum_arre2m))
            neighb_facestabs.push_back(sum_arre2m);
    }
    if (xyzqubit_coord[1]%2==0)
    {
        if (Type_in_vType(face_stabs,sum_arre0))
            neighb_facestabs.push_back(sum_arre0);
        if (Type_in_vType(face_stabs,sum_arre0m))
            neighb_facestabs.push_back(sum_arre0m);
        if (Type_in_vType(face_stabs,sum_arre2))
            neighb_facestabs.push_back(sum_arre2);
        if (Type_in_vType(face_stabs,sum_arre2m))
            neighb_facestabs.push_back(sum_arre2m);
    }
    if (xyzqubit_coord[2]%2==0)
    {
        if (Type_in_vType(face_stabs,sum_arre0))
            neighb_facestabs.push_back(sum_arre0);
        if (Type_in_vType(face_stabs,sum_arre0m))
            neighb_facestabs.push_back(sum_arre0m);
        if (Type_in_vType(face_stabs,sum_arre1))
            neighb_facestabs.push_back(sum_arre1);
        if (Type_in_vType(face_stabs,sum_arre1m))
            neighb_facestabs.push_back(sum_arre1m);
    }
    return neighb_facestabs;
}

struct fractalcoords
{
  vvint corners,neighb_facestabs;
  vint bulk_hole_coords,rm_sweep_indices,shelf_sweep_indices,face_stabs,qubit_coords,sweep_indices;
};

struct fractalcoords create_fractal(struct lattice lat,int L,int level)
{
    struct fractalcoords frac_lattice;//coordinates in the bulk of the holes && removed sweep indices due to the holes  
    vvint neighb_stabs((2*L+3)*(2*L+3)*(2*L+3));
    frac_lattice.neighb_facestabs=neighb_stabs;
    if(level==0)
    {
        frac_lattice.qubit_coords=lat.edges;
        frac_lattice.sweep_indices=lat.cubes;
        frac_lattice.face_stabs=lat.faces;
        for (const int& qubit_coord:frac_lattice.qubit_coords)
            frac_lattice.neighb_facestabs[qubit_coord]=neighb_facestabs_fun(L,qubit_coord,frac_lattice.face_stabs);
        return frac_lattice;
    }
    else if (level==1)
    {
        frac_lattice.corners.push_back({int(round(L/3)),int(round(L/3))+int(floor(L/3))+1});     
    }         
    else if (level==2)
    {
        frac_lattice.corners.push_back({int(round(L/3)),int(round(L/3))+int(floor(L/3))+1});   
        if (L==14) //L_v=15
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{2,4},{12,14}}); 
        if (L==20) //L_v=20
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{2,5},{16,19}}); 
        if (L==26)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{3,7},{21,25}}); 
        if (L==32)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{4,8},{26,30}}); 
        if (L==38)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{4,9},{30,35}}); 
    }
    else if (level==3)
    {
        frac_lattice.corners.push_back({int(round(L/3)),int(round(L/3))+int(floor(L/3))+1});            
        if (L==32)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{4,8},{26,30},{1,3},{9,11},{23,25},{31,33}});  
        if (L==38)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{4,9},{30,35},{1,3},{10,12},{27,29},{36,38}});
        if (L==44)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{5,11},{35,41},{2,4},{12,14},{32,34},{42,44}});  
        if (L==50)
            frac_lattice.corners.insert(frac_lattice.corners.end(),{{6,13},{40,46},{2,5},{14,16},{36,39},{47,49}});  
    }

    for(const vint& corner:frac_lattice.corners)
    {   
        auto holerange=range(2*corner[0],2*corner[1]-1);
        for(const int& x:holerange) //hole goes from 2*corner[0]-1 to 2*corner[1]-1
        {
            for(const int& y:holerange) 
            {
                for(const int& z:holerange)
                {
                    vint vec{x,y,z};
                    frac_lattice.bulk_hole_coords.push_back(coordinate(L,vec));
                    if ((fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==0) && (x>2*corner[0] && y>2*corner[0] && z>2*corner[0] && x<2*corner[1]-2 && y<2*corner[1]-2 && z<2*corner[1]-2))
                        frac_lattice.rm_sweep_indices.push_back(coordinate(L,vec));
                    if ((fmod(x,2)==0 && fmod(y,2)==0 && fmod(z,2)==0) && (x==2*corner[0] || y==2*corner[0] || z==2*corner[0] || x==2*corner[1]-2 || y==2*corner[1]-2 || z==2*corner[1]-2))
                        frac_lattice.shelf_sweep_indices.push_back(coordinate(L,vec)); 
                }
            }
        }
    }

    for (const int& edge:lat.edges)
        if (!(Type_in_vType(frac_lattice.bulk_hole_coords,edge)))
            frac_lattice.qubit_coords.push_back(edge);

    for (const int& cube:lat.cubes)
        if (!(Type_in_vType(frac_lattice.rm_sweep_indices,cube)))
            frac_lattice.sweep_indices.push_back(cube);   

    for (const int& face:lat.faces)
        if (!(Type_in_vType(frac_lattice.bulk_hole_coords,face)))
            frac_lattice.face_stabs.push_back(face);                
            
    for (const int& qubit_coord:frac_lattice.qubit_coords)
        frac_lattice.neighb_facestabs[qubit_coord]=neighb_facestabs_fun(L,qubit_coord,frac_lattice.face_stabs);

    return frac_lattice;
}

vint faces_fun(int L,int sweep_index,vint sweep_dir,int pastorfuture)
{
    vint faces_porf;

    vint chg=vectimesc(pastorfuture,sweep_dir);
    int Lp=2*L+3;
    int coordmaxplus=Lp*Lp*Lp;

    int face0,face1,face2;

    face0=mod(sweep_index+chg[0],coordmaxplus);
    face1=mod(sweep_index+chg[1]*Lp,coordmaxplus);
    face2=mod(sweep_index+chg[2]*Lp*Lp,coordmaxplus);

    int coordmax_face=(Lp-2)+(Lp-2)*Lp+(Lp-1)*Lp*Lp;
    int coordmin_face=1+Lp;

    if(face0<coordmax_face+1 && face0>coordmin_face-1)
        faces_porf.push_back(face0);
    if(face1<coordmax_face+1 && face1>coordmin_face-1)
        faces_porf.push_back(face1);
    if(face2<coordmax_face+1 && face2>coordmin_face-1)
        faces_porf.push_back(face2);

    return faces_porf;
}

vint onetrailhole_update(int L,vint shelf_sweep_indices,vint bulk_hole_coords,vint syndrome,vint sweep_dir)
{
    for(const int& sweep_index:shelf_sweep_indices)
    {
        vint syndromes_cubefaces_past,syndromes_cubefaces_ftr;
        vint faces_past=faces_fun(L,sweep_index,sweep_dir,-1);
        vint faces_ftr=faces_fun(L,sweep_index,sweep_dir,1);
        for (const int& face_past:faces_past)
            syndromes_cubefaces_past.push_back(syndrome[face_past]);
        for (const int& face_ftr:faces_ftr)
            syndromes_cubefaces_ftr.push_back(syndrome[face_ftr]);
        if (!(one_in_vector(syndromes_cubefaces_past)) && countone_in_vector(syndromes_cubefaces_ftr)==1)
        {
            vint zero_synd_faces_ftr;
            for(const int& face_ftr:faces_ftr)
            {
                if (Type_in_vType(bulk_hole_coords,face_ftr))
                {
                    assert(syndrome[face_ftr]==0);
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

bool check_trailing(int L,int sweep_index,vint syndrome,vint sweep_dir)
{
    vint syndromes_cubefaces_past,syndromes_cubefaces_ftr;
    vint faces_past=faces_fun(L,sweep_index,sweep_dir,-1);
    vint faces_ftr=faces_fun(L,sweep_index,sweep_dir,1);
    for (const int& face_past:faces_past)
        syndromes_cubefaces_past.push_back(syndrome[face_past]);
    for (const int& face_ftr:faces_ftr)
        syndromes_cubefaces_ftr.push_back(syndrome[face_ftr]);
    if(!(one_in_vector(syndromes_cubefaces_past)) && (countone_in_vector(syndromes_cubefaces_ftr)>1))
        return true;
    return false;
}


vint syndrome_updt_fun(int L,vint faces,vint qubit_coords,vvint neighb_facestabs,vint bulk_hole_coords,vint error,vint meas_err_faces)
{
    vint syndrome((2*L+3)*(2*L+3)*(2*L+3));
    for (const int& qubit_coord:qubit_coords)
    {
        if(error[qubit_coord]==1)
        {
            vint neighb_stabs=neighb_facestabs[qubit_coord];
            for (const int& neighb_stab:neighb_stabs)
            {
                assert(!(Type_in_vType(bulk_hole_coords,neighb_stab)));
                syndrome[neighb_stab]=fmod(syndrome[neighb_stab]+1,2);
            }
        }
    }
    for (const int& face: meas_err_faces)
        syndrome[face]=fmod(syndrome[face]+1,2);

    return syndrome;
}

struct syndromeerror sweep_step(int L,struct lattice lat,struct fractalcoords frac_lattice,struct syndromeerror synderror,vint sweep_dir,vint meas_err_faces)
{
    int Lp=2*L+3;
    int coordmaxplus=Lp*Lp*Lp;
    synderror.syndrome=onetrailhole_update(L,frac_lattice.shelf_sweep_indices,frac_lattice.bulk_hole_coords,synderror.syndrome,sweep_dir);

    for(const int& sweep_index:frac_lattice.sweep_indices)
    {
        if (check_trailing(L,sweep_index,synderror.syndrome,sweep_dir))
        {
            vint faces_ftr=faces_fun(L,sweep_index,sweep_dir,1);
            vint nonzerosyndromes_cubefaces_ftr;
            for(const int& face_ftr:faces_ftr)
            {
                if (synderror.syndrome[face_ftr]==1)
                    nonzerosyndromes_cubefaces_ftr.push_back(face_ftr);
            }

            // assert(nonzerosyndromes_cubefaces_ftr.size()!=1);

            if (nonzerosyndromes_cubefaces_ftr.size()==3)
            {
                nonzerosyndromes_cubefaces_ftr.erase(nonzerosyndromes_cubefaces_ftr.begin()+distInt0To2(rnEngine));
            }

            int edge_updtd=nonzerosyndromes_cubefaces_ftr[0]+nonzerosyndromes_cubefaces_ftr[1]-sweep_index;
            // assert(Type_in_vType(frac_lattice.qubit_coords,edge_updtd));
            synderror.error[edge_updtd]=fmod(synderror.error[edge_updtd]+1,2);
            vint neighb_stabs=frac_lattice.neighb_facestabs[edge_updtd];
            for (const int& neighb_stab:neighb_stabs)
            {
                // assert(!(Type_in_vType(frac_lattice.bulk_hole_coords,neighb_stab)));
                synderror.syndrome[neighb_stab]=fmod(synderror.syndrome[neighb_stab]+1,2);
            }
        }
    }
    for(const int& bulk_hole_coord:frac_lattice.bulk_hole_coords)
    {
        synderror.syndrome[bulk_hole_coord]=0;
    }
    // assert(synderror.syndrome==syndrome_updt_fun(L,lat.faces,frac_lattice.qubit_coords,frac_lattice.neighb_facestabs,frac_lattice.bulk_hole_coords,synderror.error,meas_err_faces));
    return synderror;
}

struct syndromeerror gen_data_err(int L,vint qubit_coords,vvint neighb_facestabs,const double derr_prob,struct syndromeerror synderror)
{
    vint veczero((2*L+3)*(2*L+3)*(2*L+3));
    // assert(synderror.syndrome==veczero);
    for (const int& qubit_coord:qubit_coords)
    {
        if (distDouble0To1(rnEngine)<= derr_prob)
        {
            synderror.error[qubit_coord]=fmod(synderror.error[qubit_coord]+1,2);
        }
        if (synderror.error[qubit_coord]==1)
        {
            vint neighb_stabs=neighb_facestabs[qubit_coord];
            for (const int& neighb_face_stab:neighb_stabs)
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
    for (const int& face_stab:face_stabs)
    {
        if (distDouble0To1(rnEngine)<= merr_prob)
        {   
            measerr_faces.push_back(face_stab);
        }
    }
    return measerr_faces;
}

int sweep_decoder_run(int level,const std::string &sweep_schedule,int rounds,int L,const double derr_prob,const double merr_prob=0)
{
    auto start = std::chrono::high_resolution_clock::now();

    int timeout=32*L;
    int same_sweep_dir_rdlimit=int(round(log(L)));
    int sweeprate=2;

    struct lattice lat=lattice_fun(L);
    
    struct fractalcoords frac_lattice=create_fractal(lat,L,level);

    struct syndromeerror synderror;
    vint veczero((2*L+3)*(2*L+3)*(2*L+3));
    synderror.error=veczero;
    synderror.syndrome=veczero;

    int same_sweep_dir_rdcount=0;

    vint sweep_dir;
    vvint dir_list{{1,1,1},{-1,1,1},{1,-1,1},{1,1,-1},{-1,-1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
    int dir_count=0;
    sweep_dir=dir_list[dir_count%8];

    for (const int& rd_index:range(rounds))
    {
        if(same_sweep_dir_rdcount==same_sweep_dir_rdlimit)
        {
            if (sweep_schedule=="alternating")
            {
                sweep_dir=dir_list[dir_count%8];
                dir_count+=1;
            }
            else if (sweep_schedule=="constant")
            {
                sweep_dir={1,1,1};
            }
            same_sweep_dir_rdcount=0;
        }
        vint veczero2((2*L+3)*(2*L+3)*(2*L+3));
        // if(rd_index==0)
        //     assert(synderror.error==veczero2);
        synderror.syndrome=veczero;
        synderror=gen_data_err(L,frac_lattice.qubit_coords,frac_lattice.neighb_facestabs,derr_prob,synderror);

        vint measerr_faces;
        if (merr_prob>0)
        {
            measerr_faces=measerr_faces_fun(synderror.syndrome,frac_lattice.face_stabs,merr_prob);    
            for (const int& face: measerr_faces)
                synderror.syndrome[face]=fmod(synderror.syndrome[face]+1,2);
        }

        for(int i=0; i<sweeprate;++i)
            synderror=sweep_step(L,lat,frac_lattice,synderror,sweep_dir,measerr_faces);
        
        same_sweep_dir_rdcount+=1;        
    }
    
    synderror.syndrome=veczero;
    vint measerr_faces;
    synderror=gen_data_err(L,frac_lattice.qubit_coords,frac_lattice.neighb_facestabs,derr_prob,synderror);

    int same_sweep_dir_limit_timeout_session=L;

    int same_sweep_dir_time=0;
    sweep_dir={1,1,1};

    for (const int& t:range(timeout))
    {
        if(same_sweep_dir_time==same_sweep_dir_limit_timeout_session)
        {
            if (sweep_schedule=="alternating")
            {
                sweep_dir=dir_list[dir_count%8];
                dir_count+=1;
            }
            else if (sweep_schedule=="constant")
            {
                sweep_dir={1,1,1};
            }
            same_sweep_dir_time=0;
        }
        synderror=sweep_step(L,lat,frac_lattice,synderror,sweep_dir,measerr_faces);

        same_sweep_dir_time+=1;

        if (!(one_in_vector(synderror.syndrome)))
            break;
    }

    int logical_failure = mod2dot_fun(lat.logical,synderror.error);
    bool syndrome_unclean= one_in_vector(synderror.syndrome);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = finish-start;
    // std::cout<<syndrome_unclean<<","<<logical_failure<<","<<diff.count()<<std::endl;

    return (syndrome_unclean || logical_failure);
}


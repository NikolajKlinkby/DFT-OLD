
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <complex> 
#include <chrono>
#include <slate/slate.hh>
# include <blas.hh> 
# include <mpi.h>


/*Multipurpose functions*/

template <typename scalar_t>
void vector_set(slate::Matrix<scalar_t> vector, int index, scalar_t value){
    int mt = vector.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector.tileMb(i) ){
            if (vector.tileIsLocal( i, 0)) {
                std::cout << "hej 3" << std::endl;
                auto tile = vector( i, 0);
                std::cout << "hej 4" << std::endl;
                auto tiledata = tile.data();
                
                tiledata[index-i*(mt+1)] = value;
                return;
                }
        }
    }
}

template <typename scalar_t>
scalar_t vector_get(slate::Matrix<scalar_t> vector, int index){
    int mt = vector.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector.tileMb(i) ){
            if (vector.tileIsLocal( i, 0)) {
                auto tile = vector( i, 0 );
                auto tiledata = tile.data();
                return tiledata[index-i*(mt+1)];
            }   
        }
    }
}

template <typename scalar_t>
scalar_t vector_get_max(slate::Matrix<scalar_t> vector){
    scalar_t max;
    for (int64_t i = 0; i < vector.mt(); i++){
        if(vector.tileIsLocal(i,0)){
            auto tile = vector(i,0);
            auto tiledata = tile.data();
            for (int64_t ii=0; ii < tile.mb(); ii++){
                if(tiledata[ii]>max){
                    max = tiledata[ii];
                }
            }
        }
    }
    return max;
}

template <typename scalar_t>
scalar_t vector_get_min(slate::Matrix<scalar_t> vector){
    scalar_t min;
    for (int64_t i = 0; i < vector.mt(); i++){
        if(vector.tileIsLocal(i,0)){
            auto tile = vector(i,0);
            auto tiledata = tile.data();
            for (int64_t ii=0; ii < tile.mb(); ii++){
                if(tiledata[ii]<min){
                    min = tiledata[ii];
                }
            }
        }
    }
    return min;
}

template <typename scalar_t>
void interp_deriv(slate::Matrix<scalar_t> x, slate::Matrix<scalar_t> y, slate::Matrix<scalar_t> y_deriv, slate::Matrix<scalar_t> x_mat, slate::Matrix<scalar_t> x_mat_deriv){
    
    /*x_mat and x_mat_deriv need to be equally created*/

    //Setting matrix
    int64_t mt = x_mat.mt();
    int64_t nt = x_mat.nt();
    int64_t n = x_mat.n();
    scalar_t x_val;

    for (int64_t i=0; i < mt; i++){
        for (int64_t j=0; j < nt; j++){
            
            

            if (x_mat.tileIsLocal(i,j) && x_mat_deriv.tileIsLocal(i,j)){
                
                

                auto tile_x = x_mat(i,j);
                auto tiledata_x = tile_x.data();
                auto tile_d = x_mat_deriv(i,j);
                auto tiledata_d = tile_d.data();

                int64_t mb = tile_x.mb();
                int64_t nb = tile_x.nb();
                int64_t stride = tile_x.stride();
                
                for (int64_t ii=0; ii < mb; ii++){
                    for (int64_t jj=0; jj < nb; jj++){
                        x_val = vector_get(x,i*(mt+1)+ii);
                        tiledata_x[ii + jj*stride] = pow(x_val , (n-1-j*(nt+1)-jj));
                        if (n-2-j*(nt+1)-jj < 0){
                            tiledata_d[ii + jj*stride] = 0;
                        }
                        else{
                            tiledata_d[ii + jj*stride] = tiledata_x[ii + jj*stride]/x_val;
                        }
                    }
                }
            }                
        }
    }
    
    
    

    //Finding derivative
    slate::Pivots pivot;
    slate::gesv(x_mat,pivot,y); //Overwrites y as output 

    slate::multiply(1.,x_mat_deriv,y,0.,y_deriv);  

}

//Defining the potentialals
template <typename scalar_t>
void ExtPot(slate::Matrix<scalar_t> lattice, slate::Matrix<scalar_t> v, double time){
    int index;
    int mt = lattice.mt();
    if (time < 1.){
        for (int64_t i = 0; i < lattice.mt(); i++){
            if (lattice.tileIsLocal( i, 0)) {
                auto tile = lattice( i, 0 );
                auto tiledata = tile.data();
                for (int64_t ii = 0; ii < tile.mb(); ++ii){
                    index = i*(mt+1)+ii;
                    if (tiledata[ii] > -3. && tiledata[ii] < 3.){
                        std::cout << "hej 1" << index << std::endl;
                        vector_set(v, index, 6. * tiledata[ii] * sin(1.*time));
                        std::cout << "nope 1" << std::endl;
                    }
                    else{
                        std::cout << "hej 2" << index << std::endl;
                        vector_set(v, index, 5. + 6. * tiledata[ii] * sin(1.*time));
                        std::cout << "nope 2" << std::endl;
                    }
                }
            } 
        }
    }
    else{
        for (int64_t i = 0; i < lattice.mt(); i++){
            if (lattice.tileIsLocal( i, 0)) {
                auto tile = lattice( i, 0 );
                auto tiledata = tile.data();
                for (int64_t ii = 0; ii < tile.mb(); ++ii){
                    index = i*(mt+1)+ii;
                    if (tiledata[ii] > -3. && tiledata[ii] < 3.){
                        vector_set(v, index, 0.);
                    }
                    else{
                        vector_set(v, index, 5.);
                    }
                }
            } 
        }
    }
}

template <typename scalar_t>
void HartreePot(slate::Matrix<scalar_t> lattice, double Delta, slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vh){
    double pot;
    int index;
    scalar_t val;
    for (int64_t j = 0; j < lattice.m(); j++){
        pot = 0;
        val = vector_get(lattice,j);
        for (int i = 0; i < lattice.mt(); i++){
            if (lattice.tileIsLocal(i,0)){
                auto tile = lattice(i,0);
                auto tiledata = tile.data();
                for (int64_t ii =0; ii < tile.mb(); ii++){
                    index = i*(lattice.mt()+1)+ii;
                    if (j != index){
                        if (val > tiledata[ii]) pot += vector_get(density,index)/(val-tiledata[ii]);
                        if (val < tiledata[ii]) pot += vector_get(density,index)/(tiledata[ii]-val);
                    }
                }
            }
        }
        vector_set(vh, j, pot*Delta);
    }
}

template <typename scalar_t>
void ExchangePot(slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vx){
    for (int64_t i = 0; i < vx.mt(); i++){
        if (vx.tileIsLocal(i,0)){
            auto tile = vx(i, 0);
            auto tiledata = tile.data();
            for (int64_t ii; ii < tile.mb(); ii++){
                //index = i*(lattice.mt()+1)+ii
                tiledata[ii] = -pow(3./M_PI*vector_get(density,i*(vx.mt()+1)+ii),1./3.);
            }
        }
    }
} 


double X(double rs, double b, double c){
    return rs + b * pow(rs, 1./2.) + c;
}

template <typename scalar_t>
void CorrelationEnergy(slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vc){
    double x0 = -0.10498; 
    double b = 3.72744;
    double c = 12.9352;
    double Q = pow(4.*c-pow(b,2.),1./2.);

    for (int i = 0; i < density.m(); i++){
        scalar_t val_i = vector_get(density, i);
        if (val_i == 0){
            vector_set(vc, i, 0.);
        }
        else{
            vector_set(vc, i, val_i * (1.-log(2.))/pow(M_PI,2.) * (
                    -log(val_i * X(1./val_i,b,c)) + 2.*b/Q * atan(Q/(2./pow(val_i,1./2.) + b)) - b*x0/X(pow(x0,2.),b,c) * (
                    log(pow(pow(1./val_i,1./2.)-x0,2.)/X(1./val_i,b,c) + 2.*(2.*x0+b)/Q * atan(Q/(2./pow(val_i,1./2.) + b))  ))  ));
        }
    }
}

template <typename scalar_t>
void CorrelationPot(slate::Matrix<scalar_t> lattice, slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vc, slate::Matrix<scalar_t> buffer, slate::Matrix<scalar_t> x_mat, slate::Matrix<scalar_t> x_mat_deriv){
    CorrelationEnergy(density, buffer);
    
    interp_deriv(lattice, buffer, vc, x_mat, x_mat_deriv);
    
}

template <typename scalar_t>
void LDAPot(slate::Matrix<scalar_t> lattice, slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vc, slate::Matrix<scalar_t> vx, slate::Matrix<scalar_t> vxc, slate::Matrix<scalar_t> buffer, slate::Matrix<scalar_t> x_mat, slate::Matrix<scalar_t> x_mat_deriv){
    ExchangePot(density, vx);
    CorrelationPot(lattice, density, vc, buffer, x_mat, x_mat_deriv);

    slate::copy(vxc, vx);
    slate::add(1.,vc,1.,vxc);
}

/* Density Functional Theory functions*/

//Generate inital guess for DFT
template <typename scalar_t>
void DFTInitDensity(int N, slate::Matrix<scalar_t> lattice, slate::Matrix<scalar_t> density){
    /*I will generate initial geusses from the infinte well*/
    scalar_t val;
    scalar_t l_0 = vector_get(lattice, 0);
    scalar_t l_1 = vector_get(lattice, lattice.m()-1 );
    scalar_t l_i;
    for (int i = 0; i < lattice.m(); i++){
        val = 0.;
        l_i = vector_get(lattice,i);
        for (int j = 0; j < N; j++){
            val += 2./(l_1-l_0)*pow(sin((j+1)*M_PI/(l_1-l_0)*(l_i-l_0)),2.);
        }
        vector_set(density, i, val);
    }
}

//Hamiltionian
template <typename scalar_t> 
void SetHamiltonian(double Delta, slate::HermitianMatrix<scalar_t> hamiltonian, slate::Matrix<scalar_t> v, slate::Matrix<scalar_t> vh, slate::Matrix<scalar_t> vxc){
    for (int64_t i=0; i < hamiltonian.mt(); i++){
        for (int64_t j=i; j < hamiltonian.nt(); j++){
            if (hamiltonian.tileIsLocal(i,j)){
                auto tile = hamiltonian(i,j);
                auto tiledata = tile.data();
                if (i == j){
                    for (int64_t ii=0; ii < tile.mb(); ii++){
                        for (int64_t jj=ii; jj < tile.nb(); jj++){ 
                            if (ii == jj){
                                tiledata[ii+jj*tile.stride()] = vector_get(v,i*(hamiltonian.mt()+1)+ii)+vector_get(vh,i*(hamiltonian.mt()+1)+ii)+vector_get(vxc,i*(hamiltonian.mt()+1)+ii)+1/pow(Delta,2);
                            }
                            else{
                                if (i*(hamiltonian.mt()+1)+ii == j*(hamiltonian.nt()+1)+jj + 1 ){
                                    tiledata[ii+jj*tile.stride()] = -1/(2 * pow(Delta,2));
                                }
                                else{
                                    tiledata[ii+jj*tile.stride()] = 0.;
                                }
                            }
                        }
                    }
                }
                else{
                    for (int64_t ii=0; ii < tile.mb(); ii++){
                        for (int64_t jj=ii; jj < tile.nb(); jj++){ 
                            if (i*(hamiltonian.mt()+1)+ii == j*(hamiltonian.nt()+1)+jj + 1 ){
                                tiledata[ii+jj*tile.stride()] = -1/(2 * pow(Delta,2));
                            }
                            else{
                                tiledata[ii+jj*tile.stride()] = 0.;
                            }
                        }
                    }
                }
            }
        }
    }
}

//Time-independent eigen solver
template <typename scalar_t>
void TIE(slate::HermitianMatrix<scalar_t> hamiltonian, std::vector<scalar_t> eval, slate::Matrix<scalar_t> evec){
    
    slate::heev(lapack::Job::AllVec,hamiltonian, eval, evec); //The solutions are normalized by default
    /*TO DO*/
    /* sort eigenvectors*/
    //Sorting the eigen values and vectors from low to high
}

//DFT generate density **and check for convergense
template <typename scalar_t>
bool DFTDensityUpdate(int N, slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> refdensity, slate::Matrix<scalar_t> evec, double threshold){
    bool succes = true;
    
    // Make refference
    slate::copy(refdensity, density);
    slate::set(1.,1.,density);
    // Set density
    slate::multiply(1.,density,evec,0.,density);
    //Check for convergence
    slate::add(1.,density,-1.,refdensity);
    if (vector_get_min(refdensity) > threshold){
        succes = false;
    }
    return succes;
}


int main(int argc , char ** argv){
    
    //Initialise MPI
    int err = 0 , mpi_provided = 0;
    err = MPI_Init_thread ( &argc , &argv , MPI_THREAD_MULTIPLE , & mpi_provided);
    assert ( err == 0 && mpi_provided == MPI_THREAD_MULTIPLE ); 

    auto start = std::chrono::high_resolution_clock::now();
    
    /* Initial parameters one can adjust **The potential is further up*/
    //Parameters
    int m = 200; //Points
    int nb = 40; //tile resolution
    double Delta = 0.05; //Spacingcout
    int N = 2; //Particles
    double threshold = 0.000001; //Threshold for convergence of SCF
    int max_it = 100; //Maximum number of iterations for SCF

    /* Solving the KS equations self consistently */
    //Initiating objects
    // Matrix(rows m, columns n, nb*nb tile block, block rows p, block columns q)
    
    //Density vector
    slate::Matrix<double> density(m,1,nb,1,1,MPI_COMM_WORLD);
    density.insertLocalTiles(); //Set tiles on CPU
    slate::Matrix<double> refdensity(m,1,nb,1,1,MPI_COMM_WORLD);
    refdensity.insertLocalTiles();

    //Lattice of system
    slate::Matrix<double> lattice(m,1,nb,1,1,MPI_COMM_WORLD);
    lattice.insertLocalTiles();

    // Potentials
    slate::Matrix<double> vxc(m,1,nb,1,1,MPI_COMM_WORLD); //Exchange correlation potential
    vxc.insertLocalTiles();
    slate::Matrix<double> vx(m,1,nb,1,1,MPI_COMM_WORLD); //Correlation potential
    vxc.insertLocalTiles();
    slate::Matrix<double> vc(m,1,nb,1,1,MPI_COMM_WORLD); //Exchange potential
    vxc.insertLocalTiles();
    slate::Matrix<double> vh(m,1,nb,1,1,MPI_COMM_WORLD); //Hartree potential
    vxc.insertLocalTiles();
    slate::Matrix<double> v(m,1,nb,1,1,MPI_COMM_WORLD); //External potential
    vxc.insertLocalTiles();

    // The Hamiltonian
    slate::HermitianMatrix<double> hamiltonian(slate::Uplo::Lower, m,nb,1,1,MPI_COMM_WORLD); 
    hamiltonian.insertLocalTiles();

    //KS eigen vektor and energies
    slate::Matrix<double> evec(m,m,nb,1,1,MPI_COMM_WORLD); //KS orbitals
    evec.insertLocalTiles();
    std::vector<double> eval;
    eval.reserve(m);

    //workspace objects
    slate::Matrix<double> mat_1(m,m,nb,1,1,MPI_COMM_WORLD); 
    mat_1.insertLocalTiles();
    slate::Matrix<double> mat_2(m,m,nb,1,1,MPI_COMM_WORLD); 
    mat_2.insertLocalTiles();
    slate::Matrix<double> buffer(m,1,nb,1,1,MPI_COMM_WORLD); 
    buffer.insertLocalTiles();

    //Defining the lattice of our system
    for(int i=0; i < m; i++){
        vector_set(lattice, i, Delta*(i-floor(1.0*m/2)));
    }

    

    int counter = 1;
    bool SelfConsistency = false; 
    

    //Calling initial guess
    DFTInitDensity(N, lattice, density);
    

    
    std::cout << "hej" << std::endl;
    ExtPot(lattice, v, 0.);
    
    /*
    std::cout << "hej" << std::endl;
    HartreePot(lattice, Delta, density, vh);
    std::cout << "hej" << std::endl;
    LDAPot(lattice, density, vc, vx, vxc,buffer,mat_1,mat_2);
    
    //Creating the hamiltonian
    SetHamiltonian(Delta, hamiltonian, v, vh, vxc);
    
    //Solving the eigenvector equation
    TIE(hamiltonian, eval, evec);
    DFTDensityUpdate(N, density, refdensity, evec, threshold);
    */

    auto stop = std::chrono::high_resolution_clock::now();

    err = MPI_Finalize ();
    assert ( err == 0 );


    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "microseconds of entire function: " << duration.count() << std::endl;

}

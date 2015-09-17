/*	main.cc		demonstration file */
//C Copyright (C) 2015 Richard Tjörnhammar
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#include "richfit.hh"
#include "iofunc.hh"

int calc_distance_matrix(gmat *A, particles coord_space) {
	int D = A->size1;
	if( A->size1 == A->size2 && A->size1 == coord_space.size()) {
		for(int i=0;i<D; i++){
			for(int j=0;j<D;j++){
				ftyp dx,dy,dz;
				dx = gsl_vector_get(coord_space[i].second,XX)-gsl_vector_get(coord_space[j].second,XX);
				dy = gsl_vector_get(coord_space[i].second,YY)-gsl_vector_get(coord_space[j].second,YY);
				dz = gsl_vector_get(coord_space[i].second,ZZ)-gsl_vector_get(coord_space[j].second,ZZ);
				ftyp len = sqrt(dx*dx+dy*dy+dz*dz);
				gsl_matrix_set( A, i, j, len );
			}
		}
		return 0;
	}else{
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: calc_distance_matrix" << std::endl;
		return 1;
	}
}


int set_coordinate_matrix( gmat *A, particles coord_space) {
	int D = A->size2;
	if( A->size2 == coord_space.size()) {
		for(int i=0;i<D; i++){
			for(int j=XX;j<=ZZ;j++){
				gsl_matrix_set( A, j, i, gsl_vector_get(coord_space[i].second,j)	);
			}
		}
		return  0;
	}else{
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: calc_distance_matrix" << std::endl;
		return -1;
	}
}

int outp_distance_matrix(gmat *A, ftyp level) {
	int D = A->size1;
	if( A->size1 == A->size2 ) {
		std::cout << " A(" << D << "," << D << ")=[" << std::endl; 
		for(int i=0; i<D; i++){
			ftyp sum=0;
			for(int j=0;j<D;j++){
				ftyp val=gsl_matrix_get(A,i,j);
				if(level!=0) {
					if(level>0) {
						val=val<level;
					}else{
						val=val>sqrt(level*level);
					}
				}
				if(val>0)
					std::cout << "\033[1;31m " << val <<"\033[0m";
				else
					std::cout << " " << val ;
				if(i!=j)
					sum += val;
			}
			std::cout << " | " << sum << std::endl;
		}
		std::cout << "];" << std::endl;
	}
	return 0;
}

std::vector<int > find_via_distance( gmat *A, ftyp level ) {
	int D = A->size1;
	std::vector<int > is_tmp;
	int sumZero = (level<0)?(1):(0);
	if (sumZero)
		level = sqrt(level*level);
	std::cout <<"INFO::LEVEL="<< level << std::endl;
	if( A->size1 == A->size2 ) {
		for(int i=0; i<D; i++){
			ftyp sum=0;
			for(int j=0;j<D;j++){
				ftyp val= gsl_matrix_get(A,i,j);
				if(level>0)
					val=val<level;
				if(i!=j)
					sum += val;
			}
			if(!sumZero)
				is_tmp.push_back(sum>1);
			else
				is_tmp.push_back(sum==0);
		}
	}
	return is_tmp;
}

int main (int argc, char **argv) {

	std::string 	filename[2];
	std::string 	ext("xyz");
	std::ifstream 	inp_molecule;
	int verbose	= 0;

	richanalysis::fileIO fIO;
	particles coord_space, carth_space;
	int N=10;

	switch( argc ) {
		case 4: {
			        std::string tmpline( argv[3] );
				filename[1] = tmpline;
                                std::size_t found = tmpline.find(ext);
                                if (found!=std::string::npos){
                                        std::cout << "INFO:: FOUND XYZ INPUT NAME = " << filename[1] << std::endl;
                                        carth_space = fIO.read_xyz(filename[1]);
                                        std::cout << "INFO:: DONE READING MODEL" << std::endl;
					if(carth_space.size()<=0) {
						std::cout << "FATAL:: FAILED TO GET COORDINATES" << std::endl;
						return(1);
					}
                                }
			}
		case 3: {
				N=atoi(argv[2]);
			}
		case 2: {
				std::string tmpline( argv[1] );
				filename[0] = tmpline;
				std::size_t found = tmpline.find(ext);
				if (found!=std::string::npos){
					std::cout << "INFO:: FOUND XYZ INPUT NAME = " << filename[0] << std::endl;
					coord_space = fIO.read_xyz(filename[0]);
					std::cout << "INFO:: DONE READING DENSITY" << std::endl;
					if(coord_space.size()<=0){
						std::cout << "FATAL:: FAILED TO GET COORDINATES" << std::endl;
						return(1);
					}
				}
			}
			break;
		default:
			std::cout << "FATAL:: FAILED TO OPEN FILE" << std::endl;
			return(1);
	}
	std::string ns=std::to_string(N);

	int D		= coord_space.size();
	int B		= carth_space.size();
	std::cout << "INFO:: GOT DIMENSION1 :: " << D << std::endl;
	std::cout << "INFO:: GOT DIMENSION2 :: " << B << std::endl;

	gmat *CRD	= gsl_matrix_calloc( DIM, D );
	gmat *CRT	= gsl_matrix_calloc( DIM, B );
	gmat *CEN	= gsl_matrix_calloc( DIM, N );
	gmat *CNT	= gsl_matrix_calloc( DIM, N );

	gmat *C0T	= gsl_matrix_calloc( DIM, N );
	gmat *C0N	= gsl_matrix_calloc( DIM, N );

	set_coordinate_matrix( CRD, coord_space );
	set_coordinate_matrix( CRT, carth_space );

	gvec *w		= gsl_vector_alloc( D );
	gvec *cw	= gsl_vector_alloc( B );
	gvec *nw	= gsl_vector_alloc( N );
	gvec *cn	= gsl_vector_alloc( N );

	gvec *z		= gsl_vector_calloc( DIM );

	gsl_vector_set_all ( w  , 1.0 );
	gsl_vector_set_all ( nw , 1.0 );
	gsl_vector_set_all ( cw , 1.0 );
	gsl_vector_set_all ( cn , 1.0 );

	richanalysis::clustering	calg;
	richanalysis::fitting		falg;
	richanalysis::tensorIO		tIO;

	calg.gsl_kmeans(CRD, w,CEN,nw);
	calg.gsl_kmeans(CRT,cw,CNT,cn);

	gsl_matrix_memcpy (C0T, CNT);
	gsl_matrix_memcpy (C0N, CEN);

	gmat *U		= gsl_matrix_calloc( DIM, DIM );
	gvec *t		= gsl_vector_calloc( DIM );
	gmat *iU	= gsl_matrix_calloc( DIM, DIM );
	gvec *it	= gsl_vector_calloc( DIM );

	std::cout << "::::INFO::::" << std::endl;
	tIO.output_matrix(CEN);

	ftyp min_rmsd=1.0E10;
	int I=0,J=0;
	for(int i=0;i<N-1;i++) {
		for(int j=i+1;j<N-1;j++) {
			gsl_matrix_memcpy (CEN,C0N);
			gsl_matrix_memcpy (CNT,C0T);
			gsl_matrix_swap_columns ( CNT, i, j );
			ftyp rmsd = falg.kabsch_fit(CEN,CNT,nw,U,t);
			if(rmsd<min_rmsd){
				I=i; J=j;
				min_rmsd=rmsd;
			}
		}
	}
	std::cout << "INFO::" << I << " " << J << std::endl;

	fIO.output_pdb("cd0n"+ns+".pdb", CEN, nw);
	fIO.output_pdb("ct0n"+ns+".pdb", CNT, cn);

	falg.invert_fit(U,t,iU,it);
	falg.apply_fit(CNT,iU,it);	// for putting it back in the original pos

	tIO.output_matrix(CNT);

	std::cout << "::::INFO::::"<< min_rmsd << std::endl;	

	if(verbose) {
		std::cout << "INFO::HAVE CENTROIDS::" <<std::endl;
		std::cout << "###########################"<<std::endl;
		tIO.output_matrix_label(CEN,nw);
		std::cout << "###########################"<<std::endl;
		tIO.output_matrix_label(CRD,w);
		std::cout << "###########################"<<std::endl;
	}

//	fIO.output_pdb("cld"+ns+".pdb", coord_space, w);

	gsl_matrix_free(CRD);
	gsl_matrix_free(CEN);
	gsl_matrix_free(C0N);
	gsl_vector_free(nw);
	gsl_vector_free(w);
	gsl_matrix_free(CRT);
	gsl_matrix_free(CNT);
 	gsl_matrix_free(C0T);
	gsl_vector_free(cw);
	gsl_vector_free(cn);
	gsl_matrix_free(U);
	gsl_vector_free(t);
	gsl_matrix_free(iU);
	gsl_vector_free(it);
	gsl_vector_free(z);
	return 0;
}

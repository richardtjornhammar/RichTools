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
#include "cluster.hh"

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

std::vector<int>
outp_distance_matrix(gmat *A, ftyp level) {
	int D = A->size1;
	std::vector<int> vi;
	if( A->size1 == A->size2 ) {
		std::cout << " A(" << D << "," << D << ")=[" << std::endl; 
		for(int i=0; i<D; i++){
			ftyp sum=0,sumb=0;
			for(int j=0;j<D;j++){
				ftyp val=gsl_matrix_get(A,i,j);
				ftyp valb=0;
				if(level!=0) {
					if(level>0) {
						valb=val<level;
					}else{
						valb=val>sqrt(level*level);
					}
				}
				if(valb>0)
					std::cout << "\033[1;31m " << valb <<"\033[0m";
				else
					std::cout << " " << valb ;
				val*=valb;
				if(i!=j){
					sum  += val;
					sumb +=valb;
				}
			}
			vi.push_back(sumb);
			std::cout << " | " << sumb << " | " << sum << std::endl;
		}
		std::cout << "];" << std::endl;
	}
	return vi;
}

std::vector<int>
compare_dist_matrices(gmat *A, gmat *B, ftyp val){

	std::vector<int> vi;
	if(A->size1==B->size1&&A->size2==B->size2&&A->size1==B->size2)
	{
		std::cout << "COMPARING..." << std::endl;
		gvec *va = gsl_vector_calloc(A->size1);
		gvec *vb = gsl_vector_calloc(B->size2);

		for(int i=0;i<A->size1;i++) 
		{
			ftyp diff=1E10,dotp;
			int I=0;
			for(int j=0;j<B->size1;j++) 
			{
				gsl_matrix_get_row ( va, A, i );
				gsl_matrix_get_row ( vb, B, j );
				gsl_vector_sub(va,vb);
				gsl_blas_ddot (va, va, &dotp);
				if(dotp<diff)
				{
					I = j;
					diff=dotp;
				}
			}
			vi.push_back(I);
		}
	}

	return vi;
}

int outp_distance_matrix(gmat *A) {
	int D = A->size1;
	if( A->size1 == A->size2 ) {
		std::cout << " A(" << D << "," << D << ")=[" << std::endl; 
		for(int i=0; i<D; i++){
			ftyp sum=0;
			for(int j=0;j<D;j++){
				ftyp val=gsl_matrix_get(A,i,j);

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
	int N		= 10;
	bool bFit	= true;

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
				bFit=(N>=0);
				if(!bFit)
					N*=-1;
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
	std::string ns	= std::to_string(N);

	int D		= coord_space.size();
	int B		= carth_space.size();

	richanalysis::cluster cl1, cl2, cl_test_d, cl_test_m;

	cl1.alloc_space(D,N); 
	cl2.alloc_space(B,N);

	cl1.set_matrix( coord_space ); 	
	cl2.set_matrix( carth_space );

	cl1.perform_clustering(); 
	cl2.perform_clustering();
	std::vector<int> den_ndx = cl1.get_labels();
	std::vector<int> crd_ndx = cl2.get_labels();
	
	richanalysis::coord_format cftool;
	std::vector<std::string > anames;
	for(int i=0;i<N;i++)
		anames.push_back("C");

	std::vector<int> idx1,idx2;

	fIO.output_pdb("mod-nofit-n"+ns+".pdb", carth_space , crd_ndx);
	fIO.output_pdb("den-nofit-n"+ns+".pdb", coord_space , den_ndx);

	if(!bFit){
// NCULSTS 1/12 with H 1/8 without
		gsl_matrix *C = gsl_matrix_alloc(DIM,N);
		gsl_vector *w = gsl_vector_alloc(N);
		cl1.copyC(C); // the centroids are unordered
		cl1.copyw(w);
		
		particles centroids = cftool.mat2par(C,w);
		gmat *dM = gsl_matrix_calloc(N,N);
		std::cout << "INFO" << centroids.size()<<std::endl;
		calc_distance_matrix(dM, centroids);
		idx1=outp_distance_matrix(dM, 1.4); // 1.39 sqrt(3)*1.39 2*1.39
		cl1.connectivity(dM,1.4);

		if(N==B){
			gsl_matrix *M = gsl_matrix_alloc(DIM,B);
			gsl_vector *v = gsl_vector_alloc(B);
			cl2.copyM(M);
			cl2.copyv(v);
			particles centrs = cftool.mat2par(M,v);
			gmat *dN = gsl_matrix_calloc(B,B);
	
			calc_distance_matrix(dN, centrs);
			idx2=outp_distance_matrix(dN, 1.4);
			cl2.connectivity(dN, 1.4);
		}
		std::cout << "IDX1 " << std::endl;
		int sum=0;
		for(int i=0;i<idx1.size();i++){
			std::cout << " " << idx1[i] ;
			sum+=idx1[i];
		}
		std::cout << "SUM " << sum << std::endl;
		sum=0;
		std::cout << "IDX2 " << std::endl;
		for(int i=0;i<idx2.size();i++){
			std::cout << " " << idx2[i] ;
			sum+=idx2[i];
		}
		std::cout << "SUM " << sum << std::endl;
		std::cout << std::endl;
		fIO.output_pdb("cen-nofit-n"+ns+".pdb", C , anames);

		return (0);
	}

	int isw = 0;

	richanalysis::cluster_node nidx;

	nidx.assign_sub(cl2,cl1);
	std::vector<int> fio_ndx, rel_ndx = nidx.find_centroid_relation();
	particles align_space = nidx.apply_rot_trans( carth_space );
	for(int i=0; i<crd_ndx.size(); i++) {
		fio_ndx.push_back( rel_ndx[crd_ndx[i]] );
	}
	fIO.output_pdb("color-rel-n"+ns+".pdb", align_space , fio_ndx);

/*
	ftyp rmsd = 0.0;
	rmsd = nidx.find_centroid_relation(cl1,cl2);

	std::vector<int> rel_ndx = nidx.get_indices();
	std::vector<int> fio_ndx;

	int icl=1;
	int nden=0, nmod=0;

	for(int i=0; i<den_ndx.size(); i++) {
		nden += den_ndx[i]==icl?1:0;
	}



	particles align_space;
	align_space = nidx.apply_rot_trans( carth_space );

	int M=N;
	std::vector<richanalysis::cluster > vmode;

	for( int ipart=0;ipart<N;ipart++){
		richanalysis::coord_format cf;
		particles px;
		px	= cf.par2par(align_space, fio_ndx, ipart);
		D	= px.size();
		richanalysis::cluster clpx;
		if( D<M )
			M=D;
		clpx.alloc_space( D, M );
		clpx.set_matrix( px );
		clpx.perform_clustering();
		clpx.find_shape();
		vmode.push_back(clpx);
	}
	richanalysis::cluster_node mrel;
	for(int i=0;i<vmode.size()-1;i++){
		for(int j=i+1;j<vmode.size();j++)
			std::pair<ftyp,ftyp > angle_distance = 
				mrel.angle_between(vmode[i],vmode[j],i,j);
	}
	fIO.output_pdb("mod"+ns+".pdb", align_space, fio_ndx);

//	HERE WE ALIGN FRAGMENTS
	particles pfrag;
	std::vector<int> ord_ndx;
	M=N;
	std::vector<richanalysis::cluster > vclus;
	for( int ipart=0;ipart<N;ipart++){
		particles px,dx,rtx;
		richanalysis::coord_format cf;

		px=cf.par2par(carth_space, fio_ndx, ipart);
		dx=cf.par2par(coord_space, den_ndx, ipart);

		D	= px.size();
		B	= dx.size();
		if( D<M || B<M )
			M = ((int)(D<B)?(D):(B));

		richanalysis::cluster clpx, cldx; 
		clpx.alloc_space( D, M );
		cldx.alloc_space( B, M );
		clpx.set_matrix( px );
		cldx.set_matrix( dx );
		cldx.perform_clustering();	
		clpx.perform_clustering();

		cldx.find_shape();
		vclus.push_back(cldx);
	
		richanalysis::cluster_node cidx;

//		if(isw) //WE ARE ALREADY IN CORRECT CENTROID SO JUST SHAPEFIT
//			rmsd = cidx.find_shape_relation(cldx,clpx);
//		else
			rmsd = cidx.find_centroid_relation(cldx,clpx);

		rtx=cidx.apply_rot_trans( px );
		pfrag=cf.app_par(rtx,pfrag);
		for( int j=0;j<rtx.size(); j++ )
			ord_ndx.push_back( ipart );
	}
	fIO.output_pdb("frag"+ns+".pdb", pfrag , ord_ndx );

	richanalysis::cluster_node crel;
	for(int i=0;i<vclus.size()-1;i++){
		for(int j=i+1;j<vclus.size();j++)
			std::pair<ftyp,ftyp > angle_distance = 
				crel.angle_between(vclus[i],vclus[j],i,j);
	}
	for(int i=0;i<vclus.size();i++){
		coord_space.push_back(vclus[i].normal());
		coord_space.push_back(vclus[i].center());
		den_ndx.push_back(i); den_ndx.push_back(i);
	}

//	THIS IS THE DENSITY WE ARE FITTING TO
	std::string alabel("H");
	fIO.output_pdb("dens"+ns+".pdb", coord_space , den_ndx , alabel);
*/
	return 0;
}

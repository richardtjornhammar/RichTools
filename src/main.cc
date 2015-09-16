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

//// C/C++ STUFF
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>

#include "richfit.hh"

int output_geometry_to_file(std::vector<std::pair<std::string,std::vector<ftyp> > > px, std::string filename){
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);

	int n=px.size();
	outp_coord << n << std::endl;
	outp_coord << "FLUSHED COORDS" << std::endl;
	for(int i=0;i<n;i++){
		outp_coord << " " << px[i].first << " " << px[i].second[XX] << " " 
			<< px[i].second[YY] << " " << px[i].second[ZZ] << std::endl; 
	}
	outp_coord.close();

	return(0);
}

int output_geometry(std::vector<std::pair<std::string,std::vector<ftyp> > > px) {
	int n=px.size();
	std::cout << n << std::endl;
	std::cout << "FLUSHED COORDS" << std::endl;
	for(int i=0;i<n;i++){
		std::cout << " " << px[i].first << " " << px[i].second[XX] << " " 
			<< px[i].second[YY] << " " << px[i].second[ZZ] << std::endl; 
	}
	return(0);
}

int output_tagged_geometry(std::vector<std::pair<std::string,std::vector<ftyp> > > px, std::vector<int > is_ring, std::string label){
	int n=px.size();
	std::cout << n << std::endl;
	std::cout << "FLUSHED RING COORDS" << std::endl;
	for(int i=0;i<n;i++){
		if(is_ring[i])
			std::cout << " " << label << " " << px[i].second[XX] << " " 
			<< px[i].second[YY] << " " << px[i].second[ZZ] << std::endl; 
		else
			std::cout << " " << px[i].first << " " << px[i].second[XX] << " " 
			<< px[i].second[YY] << " " << px[i].second[ZZ]  <<std::endl; 
	}
	return(0);	
}

int calc_distance_matrix(gmat *A, std::vector<std::pair<std::string,std::vector<ftyp> > >  coord_space) {
	int D = A->size1;
	if( A->size1 == A->size2 && A->size1 == coord_space.size()) {
		for(int i=0;i<D; i++){
			for(int j=0;j<D;j++){
				std::vector<ftyp> diff;
				diff.push_back(coord_space[i].second[XX]-coord_space[j].second[XX]);
				diff.push_back(coord_space[i].second[YY]-coord_space[j].second[YY]);
				diff.push_back(coord_space[i].second[ZZ]-coord_space[j].second[ZZ]);
				ftyp len = sqrt(diff[XX]*diff[XX]+diff[YY]*diff[YY]+diff[ZZ]*diff[ZZ]);
				gsl_matrix_set( A, i, j,	len	);
			}
		}
		return 0;	
	}else{
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: calc_distance_matrix" << std::endl;
		return 1;
	}
}


int set_coordinate_matrix( gmat *A, std::vector<std::pair<std::string,std::vector<ftyp> > >  coord_space) {
	int D = A->size2;
	if( A->size2 == coord_space.size()) {
		for(int i=0;i<D; i++){
			for(int j=XX;j<=ZZ;j++){
				gsl_matrix_set( A, j, i,	coord_space[i].second[j]	);
			}
		}
		return 0;	
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
		level=sqrt(level*level);
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


std::vector<std::pair<std::string,std::vector<ftyp> > > read_xyz(std::string filename){
	const char *c_filename = filename.c_str();
	std::ifstream inp_molecule;
	inp_molecule.open(c_filename);
	int Nparts = 0;
	std::string title;
	std::vector<ftyp> BOX;
	BOX.push_back(1.0); BOX.push_back(1.0); BOX.push_back(1.0);
	std::vector<std::pair<std::string,std::vector<ftyp> > > particles;
	if(inp_molecule.is_open()){
		std::string iline;
		getline(inp_molecule,iline);
		std::stringstream data(iline);
		data >> Nparts;
		getline(inp_molecule, iline);
		title = iline;
		std::stringstream atad(iline);
		std::string wrd; atad >> wrd;
		if( wrd == "BOX" ){
			atad >> BOX[XX]; atad >> BOX[YY]; atad >> BOX[ZZ];
		}
		while( !(inp_molecule.eof()) ){
//	DECLARATIONS
			ftyp r[3];
			std::pair<std::string, std::vector<ftyp> > svpair;
			std::string atom;
			std::vector<ftyp> rvec;
//	ACTUAL IO
			getline(inp_molecule,iline);
			std::stringstream    data(iline);
			data >> atom >> r[0] >> r[1] >> r[2];
			rvec.push_back(r[0]*BOX[XX]);
			rvec.push_back(r[1]*BOX[YY]);
			rvec.push_back(r[2]*BOX[ZZ]);
			svpair.first = atom; svpair.second = rvec;
			if(svpair.first=="")
				continue;
			particles.push_back(svpair);
		}
		inp_molecule.close();
	}else{
		std::cout << "FATAL:: COULD NOT OPEN FILE " << std::endl;
	}

	return particles;
}

int main (int argc, char **argv) {

	std::string 	filename[1];
	std::string 	ext("xyz");
	std::ifstream 	inp_molecule;
	int ioOpen[1]	=	{0};

	std::vector<std::pair<std::string, std::vector<ftyp> > > coord_space;

	switch( argc ) {
		case 2: {
				std::string tmpline( argv[1] );
				filename[0] = tmpline;
				std::size_t found = tmpline.find(ext);
				if (found!=std::string::npos){
					std::cout << "INFO:: FOUND XYZ INPUT NAME" << std::endl;
					coord_space = read_xyz(filename[0]);
				}
			}	
			break;
		default:
			std::cout << "FATAL:: FAILED TO OPEN FILE" << std::endl;
			return(1);
	}

	output_geometry(coord_space);
	int D		= coord_space.size();

	std::cout << "INFO:: GOT DIMENSION " << D << std::endl;

	gmat *CRD	= gsl_matrix_alloc( DIM, D );
	gsl_matrix_set_all ( CRD, 0.0	);	
	set_coordinate_matrix( CRD, coord_space );
	gvec *w	= gsl_vector_alloc( D );
	gsl_vector_set_all( w, 1.0 );

	int N		= 23;
	gmat *CEN	= gsl_matrix_alloc( DIM, N );
	gsl_matrix_set_all ( CEN, 0.0 );
	gvec *nw	= gsl_vector_alloc( N );
	gsl_vector_set_all(nw, 1.0);

	richanalysis::clustering	calg;
	richanalysis::tensorIO		tIO;
	calg.gsl_kmeans(CRD,w,CEN,nw);

	std::cout << "INFO::HAVE CENTROIDS::" <<std::endl;
	tIO.output_matrix_label(CEN,nw);
	
	std::cout << "###########################"<<std::endl;
	std::cout << D << std::endl;
	std::cout << "INFO::COLOURED COORDS::" <<std::endl;
	tIO.output_matrix_label(CRD,w);

	gsl_matrix_free(CRD);
	gsl_matrix_free(CEN);
	gsl_vector_free(nw);

	return 0;
}

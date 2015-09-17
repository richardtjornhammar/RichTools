/*	iofunc.hh	*/
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

#include "iofunc.hh"
#include <iomanip>

namespace richanalysis {

void 
fileIO::output_geometry(particles px) {
	int n=px.size();
	std::cout << n << std::endl;
	std::cout << "FLUSHED COORDS" << std::endl;
	for(int i=0;i<n;i++){
		std::cout << " " << px[i].first << " " << gsl_vector_get(px[i].second,XX) << " " 
			<< gsl_vector_get(px[i].second,YY) << " " << gsl_vector_get(px[i].second,ZZ) << std::endl; 
	}
}

void 
fileIO::output_geometry( particles px, std::string filename){
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);

	int n=px.size();
	outp_coord << n << std::endl;
	outp_coord << "FLUSHED COORDS" << std::endl;
	for(int i=0;i<n;i++){
		outp_coord 	<< " " << px[i].first << " " << gsl_vector_get(px[i].second,XX) << " " 
				<< gsl_vector_get(px[i].second,YY) << " " << gsl_vector_get(px[i].second,ZZ) << std::endl; 
	}
	outp_coord.close();
}

void 
fileIO::output_geometry( particles px, std::string filename, std::string label){
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);

	int n=px.size();
	outp_coord << n << std::endl;
	outp_coord << label << std::endl;
	for(int i=0;i<n;i++){
		outp_coord	<< " " << px[i].first << " " << gsl_vector_get(px[i].second,XX) << " " 
				<< gsl_vector_get(px[i].second,YY) << " " << gsl_vector_get(px[i].second,ZZ) << std::endl; 
	}
	outp_coord.close();	
}

particles 
fileIO::read_xyz(std::string filename) {
	const char *c_filename = filename.c_str();
	std::ifstream inp_molecule;
	inp_molecule.open(c_filename);
	int Nparts = 0;
	std::string title;
	std::vector<ftyp> BOX;
	BOX.push_back(1.0); BOX.push_back(1.0); BOX.push_back(1.0);
	particles ps;
	if(inp_molecule.is_open()){
		std::string iline;
		getline(inp_molecule,iline);
		std::stringstream data(iline);
		data >> Nparts;
		getline(inp_molecule, iline);
		title = iline;
		std::stringstream atad(iline);
		std::string wrd; atad >> wrd;
		if( wrd == "BOX" ){ // CAN SET A SCALE FOR COORDINATES (ORTHOGONAL SPACE)
			atad >> BOX[XX]; atad >> BOX[YY]; atad >> BOX[ZZ];
		}

		while( !(inp_molecule.eof()) ){
//	DECLARATIONS
			ftyp r[3];
			particle sp;
			std::string atom;
			gvec *rvec = gsl_vector_calloc(DIM);
//	ACTUAL IO
			getline(inp_molecule,iline);
			std::stringstream    data(iline);
			data >> atom >> r[0] >> r[1] >> r[2];
			gsl_vector_set(rvec,XX,r[XX]*BOX[XX]);
			gsl_vector_set(rvec,YY,r[YY]*BOX[YY]);
			gsl_vector_set(rvec,ZZ,r[ZZ]*BOX[ZZ]);
			sp.first = atom; sp.second = rvec;
			if(sp.first=="")
				continue;
			ps.push_back(sp);
		}
		inp_molecule.close();
	}else{
		std::cout << "FATAL:: COULD NOT OPEN FILE " << std::endl;
	}

	return ps;
}


void 
fileIO::output_pdb( std::string filename, particles px, gvec *w ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;
	fins << std::setw(5) << ((int)round(px.size())) << std::setw(5) << 0;

	if( w->size == px.size() ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0;i<px.size();i++){
			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< px[i].first	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+round(gsl_vector_get(w,i))) )
				<< std::setw(4) << 1 << "    "	<< std::setprecision(4) 
				<< std::setw(8) << gsl_vector_get(px[i].second,XX) 
				<< std::setw(8) << gsl_vector_get(px[i].second,YY)
				<< std::setw(8) << gsl_vector_get(px[i].second,ZZ) << std::setprecision(3)
				<< std::setw(6) << 1.00 << std::setw(6) << 0.00 << std::setw(12) << px[i].first << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;

	outp_coord.close();
}


void 
fileIO::output_pdb( std::string filename, gmat *M, gvec *w ) {
//NOTE:: pymol can color the clusters with CMD:: util.cbc
	const char *c_filename = filename.c_str();
	std::ofstream outp_coord;
	outp_coord.open(c_filename);
	std::string title = "COMPND    CLUSTER";
	std::string autho = "AUTHOR    GENERATED BY RICHTOOL";
	std::string fin1  = "MASTER        0    0    0    0    0    0    0    0";
	std::stringstream    fins;
	ftyp xX,yY,zZ;

	int n = M->size1, m = M->size2;
	int L = n==DIM?m:n;

	fins << std::setw(5) << L << std::setw(5) << 0;
	fins << std::setw(5) << L << std::setw(5) << 0;

	if( ((int)w->size) == L ){
		outp_coord << title << std::endl;
		outp_coord << autho << std::endl;
		for(int i=0; i<L; i++ ){

			xX=(n==DIM)?(gsl_matrix_get(M,XX,i)):(gsl_matrix_get(M,i,XX));
			yY=(n==DIM)?(gsl_matrix_get(M,YY,i)):(gsl_matrix_get(M,i,YY));
			zZ=(n==DIM)?(gsl_matrix_get(M,ZZ,i)):(gsl_matrix_get(M,i,ZZ));

			outp_coord << std::right<< std::setw(6) << "HETATM" << std::setw(5) << i+1 << std::setw(3)
				<< 'C'	<< std::setw(6) << "LIG" << std::setw(2) << ( (char) (65+round(gsl_vector_get(w,i))) )
				<< std::setw(4) << 1  << "    "	<< std::setprecision(4) 
				<< std::setw(8) << xX << std::setw(8) << yY << std::setw(8) << zZ 
				<< std::setprecision(3) << std::setw(6) << 1.00 
				<< std::setw(6) << 0.00 << std::setw(12) << 'C' << std::endl;
		}
	}
	outp_coord << fin1 << fins.str() << std::endl;
	outp_coord << "END" << std::endl;
	outp_coord.close();
}


}

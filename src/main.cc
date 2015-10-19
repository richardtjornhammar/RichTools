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

int main (int argc, char **argv) {

	std::string 	filename[2];
	std::string 	ext("xyz");
	std::ifstream 	inp_molecule;
	int verbose	= 0;

	richanalysis::fileIO fIO;
	particles coord_space, carth_space;
	int N		= 3;
	bool bFit	= true;

	switch( argc ) {
		case 3: {
			        std::string tmpline( argv[2] );
				filename[1] = tmpline;
                                std::size_t found = tmpline.find(ext);
                                if ( found != std::string::npos ){
                                        std::cout << "INFO:: FOUND XYZ INPUT NAME = " << filename[1] << std::endl;
                                        carth_space = fIO.read_xyz(filename[1]);
                                        std::cout << "INFO:: DONE READING MODEL" << std::endl;
					if( carth_space.size() <= 0 ) {
						std::cout << "FATAL:: FAILED TO GET COORDINATES" << std::endl;
						return(1);
					}
                                }
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

////DEBUG SECTION:: SINGLE PASS
/*
	richanalysis::node	n0;
	richanalysis::cluster	clu, clu2, clu1;

	n0.first .set_matrix(coord_space);
	n0.second.set_matrix(carth_space);

	gsl_matrix *gm = gsl_matrix_calloc( DIM, B );
	gsl_vector *gv = gsl_vector_calloc( B );
	n0.second.copyM(gm);
	n0.second.copyv(gv);	
	richanalysis::node_analysis	anode(n0);
	std::vector<int> clu_ndx	= anode.cluster_index_order();
	fIO.output_pdb("clu-nofit-n" + ns + ".pdb", carth_space , clu_ndx ); 
	richanalysis::fileIO fio;
	fio.output_pdb( "diagn3.pdb" , gm, gv);
	return 0;
*/
////END DEBUG

	richanalysis::layer	old_layer;
	richanalysis::layer	solved_layer;
	richanalysis::node	n0;
	if(D>=B){
		n0.first.set_matrix(coord_space);
		n0.second.set_matrix(carth_space);
	}else{
		n0.second.set_matrix(coord_space);
		n0.first.set_matrix(carth_space);
	}

	int sw = (B>=D) + 1;
	old_layer.push_back(n0);

	int C = 1;

	std::vector<int> clu_ndx,cen_ndx;
	while ( solved_layer.size() != carth_space.size() )
	{
		richanalysis::layer current_layer;
		current_layer.clear();
		while ( !old_layer.empty() )
		{
			richanalysis::node wnode	= old_layer.back();
			old_layer.pop_back();
			richanalysis::node_analysis	anode(wnode);
			std::pair<int, int> S = anode.size();

			if ( ( S.first==1 ) || ( S.second==1 ) )
			{
				solved_layer.push_back(	wnode );
			} else { 
				richanalysis::layer new_layer	= anode.get_node_layer();
				current_layer.insert( current_layer.end() , new_layer.begin() , new_layer.end() );
			}

			if ( C++ == 1 )
			{
				clu_ndx	= anode.cluster_index_order();
				cen_ndx = n0.first.get_clabels();
				std::cout << "INFO:: HAVE NDX SIZE " << clu_ndx.size() << std::endl; 
				std::vector<int> den_ndx	= n0.first.get_labels();
				fIO.output_pdb("clu-nofit-n" + ns + ".pdb", carth_space , clu_ndx ); 
				fIO.output_pdb("den-nofit-n" + ns + ".pdb", coord_space , den_ndx );
			}
		}

		while ( !current_layer.empty() )
		{
			richanalysis::node wnode	= current_layer.back(); 
			current_layer.pop_back();
			richanalysis::node_analysis 	anode(wnode);
			std::pair<int,int> S = anode.size();
			
			if ( ( S.first==1 ) || ( S.second==1 ) )
			{
				solved_layer.push_back(	wnode );
			}
			else
			{
				old_layer.push_back(	wnode );
			}
		}
	}

	richanalysis::layer_analysis	alayer(solved_layer);
	alayer.output_layer("solvedThis.pdb");

//// ALTERNATIVE
	std::cout << " ALTERNATIVE " << std::endl;
	richanalysis::particle_analysis pa;
	pa.assign_particles(coord_space, carth_space);
	std::vector< std::pair<ftyp, std::pair< int, int > > > atomlist = pa.compare_dist_matrices(2.0);
	std::cout << "INFO::"  << atomlist.size() << std::endl;

	int imax=0;
	for(int i=0 ; i< atomlist.size() ; i++ ){
		imax = ( atomlist[i].second.first )>(imax)?( atomlist[i].second.first ):(imax);
		imax = ( atomlist[i].second.second)>(imax)?( atomlist[i].second.second):(imax);
	}

	std::cout << "INFO::LARGEST::INDEX::" << imax << std::endl;
	int i1,i2;
	std::vector<int> vi1, vi2, vi_fin;

	for( int i=0 ; i<B ; i++ ) {
		vi1.push_back(1);
		vi2.push_back(1);
		vi_fin.push_back(-1);
	}

	for( int i=0 ; i<atomlist.size() ; i++ ) {
		i1 = atomlist[i].second.first ;
		i2 = atomlist[i].second.second;
		if( vi_fin[i1] < 0 && vi2[i2] ) {
			vi_fin[i1]	= i2;
			vi2[i2]		=  0;
			vi1[i1]		=  0;
			imax--;
		}
		if(imax == 0)
			break;
	}

	std::cout << "INFO::HAVE::RELATION::" << std::endl;
	for( int i=0 ; i<B ; i++ ) 
		std::cout << " " << vi_fin[i]  << std::endl;
	std::cout << std::endl;

	std::cout << "INFO" << clu_ndx.size() << std::endl;
	for( int i = 0 ; i < clu_ndx.size() ; i++ )
		std::cout << " " <<  clu_ndx[i] ;

	std::cout << std::endl;

	pa.output_result("pastuff.pdb");

	return 0;
}

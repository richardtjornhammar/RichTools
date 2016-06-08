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
#include "gradient.hh"
#include <algorithm>

bool nsort_func(std::pair<double, std::pair< int, int > > i1,std::pair<ftyp, std::pair< int, int > > i2) {
	return (i1.first<i2.first);
};

int main (int argc, char **argv) {
//
//	ARG PARSING INIT
//
	std::string			progname("richtool");
	std::string			author("\nRichard Tjörnhammar \nEmail: richard.tjornhammar@gmail.com\n");
	std::vector< std::pair< std::string, std::string > > opts;
	std::vector< std::string >	args;
	if( argc < 2 ) {
		args.push_back(argv[0]);
		args.push_back("-h");
	} else {
		for( int i=0 ; i<argc ; i++ )
			args.push_back(argv[i]);
	}
	std::pair< std::string , std::string > poop;
	poop.first  = "-iden";
	poop.second = "inpd.xyz";
	opts.push_back(poop);
	poop.first  = "-imod";
	poop.second = "inpm.xyz";
	opts.push_back(poop);
	poop.first  = "-omod";
	poop.second = "outm.pdb";
	opts.push_back(poop);
	poop.first  = "-oclu";
	poop.second = "outd.pdb";
	opts.push_back(poop);
//
	std::vector<std::string> 	iv_ext;
	iv_ext.push_back("xyz");
	std::vector<std::string> 	ov_ext;
	ov_ext.push_back("xyz");
	ov_ext.push_back("pdb");
//
//	ARG PARSE
//
	int arg = 0;
	int failed = 0;
//
	while ( ++arg < args.size()  ) {
		if	( args[arg] == "-h" ||  args[arg] =="--h" || args[arg] =="--help" || args[arg] =="-help") {
			std::cout << "Usage: " << progname << "\n\t";
			for( int opt=0 ; opt<opts.size() ; opt++ ){
				std::cout << opts[opt].first << " \t " << opts[opt].second << "\n\t";
			}
			std::cout << "-h" << "\n\t" << "--h" << "\n\t"<< "--help" << "\n\t"<< "-help" << "\n\n";
			std::cout << "Please contact" << author << "with any questions, comments or problems\n";
			failed = 2;
			break;
		} else {
			int ia0=arg, ia=++arg;
			for( int opt=0; opt<opts.size(); opt++ ) {
				if ( args[ia0] == opts[opt].first ) {
					if ( ia < args.size()  && !(args[ia][0]=='-') )
						opts[opt].second = args[ia];
					else 
						failed = 1;
					break;
				}
			}
			if( failed == 1 )
				break;
		}
	}
	if(failed) {
		if(failed==1)
			std::cout << "Unrecognized:\t" << args[arg] << "\n";
		args.clear();
		return 0;
	}
//
//	ARG INIT
//
	richanalysis::fileIO fIO;
	particles vparts[2];
	int L[2]={0,0};
	for(int i=0;i<2;i++) {
		std::string tmpline( opts[i].second );
                std::size_t found = tmpline.find(iv_ext[0]);
		if ( found != std::string::npos ) {
                	std::cout << "INFO:: FOUND VALID INPUT:: \t "	<< opts[i].first 
				<< " :: " << opts[i].second << std::endl;
               		particles p_space = fIO.read_xyz(opts[i].second);
                        std::cout << "INFO:: DONE READING COORDINATES" << std::endl ;
			L[i] = (int)p_space.size();
			if( p_space.size() <= 0 ) {
				std::cout << "FATAL:: FAILED TO GET COORDINATES" << std::endl;
				return(1);
			}
			vparts[i].swap(p_space);
		}else{
			std::cout << "FATAL:: FAILED" << std::endl;
			std::cout << "FATAL:: " <<  opts[i].first << " AND " << opts[i].second << std::endl;
			return(1);
		}
	}
//	
//	PROGRAM INIT
//
	int verbose	= 0;
	int N		= 3;
	bool bFit	= true;
	std::string ns	= std::to_string(N);

	richanalysis::node	n0;
	particles model;
	particles densi;

	if( L[0] >= L[1] ) {
		n0.first.set_matrix(  vparts[0] ); densi=vparts[0];
		n0.second.set_matrix( vparts[1] ); model=vparts[1];
	} else { 
		n0.first.set_matrix(  vparts[1] ); model=vparts[0];
		n0.second.set_matrix( vparts[0] ); densi=vparts[1];
	}
	int s=floor(model.size()/8);
	s=s<3?3:s; 
	n0.first .set_cDIM(s);	
	n0.first .realloc_centroids(s);
	n0.second.set_cDIM(s);	
	n0.second.realloc_centroids(s);

	int sw	= (L[1]>=L[0]) + 1;
	int C	= 1;

// B	ORDER DEPENDENT
	richanalysis::node_analysis nnode(n0);
	particles c_aligned = nnode.regular_fit();	// shape fit
	n0.first .find_centroids();
	n0.second.find_centroids();
//	
	std::vector< int > n0flabels = n0.first .get_labels();
	std::vector< int > n0slabels = n0.second.get_labels();
// ONLY FOR SHOW
	particles cents1 = nnode.get_centroids( 1 );
	fIO.output_pdb("cen" + std::to_string ( 1 ) + opts[2].second , cents1 ); 
 	particles cents2 = nnode.get_centroids( 2 );
	fIO.output_pdb("cen" + std::to_string ( 2 ) + opts[2].second , cents2 ); 
// E
// E
	for( int i = 0 ; i < c_aligned.size() ; i++ )
		c_aligned[i].first = model[i].first; 				// smaller returned
	std::vector<int> ndx12 = nnode.find_centroid_distance_relation(); 	// 2->1 mapping
	for(int i=0;i<ndx12.size();i++)
		std::cout << "INDEX PAIR :: \t" << i << " \t "<< ndx12[i] << std::endl;
	for(int i=0;i<n0slabels.size();i++)
		n0slabels[i] = ndx12[ n0slabels[i] ];
	n0.first .calculate_neighbors();
	n0.second.calculate_neighbors();
//	
// NEW
	nnode.nn_restraint_fit(0);
// NEW
//	
	fIO.output_pdb("nrm" + std::to_string(s)
			 + opts[2].second , c_aligned
			, n0slabels ); 
	fIO.output_pdb("nrd" + std::to_string(s)
			 + opts[3].second , densi 	
			, n0flabels ); 
//	
//	HERE USE THE NRM AS INPUT FOR THE NRD CENTROIDS
//
	richanalysis::cluster cluless( densi , c_aligned );	//	THIS CONTAINS FRAGMENTED PSEUDO SOLUTION
	int Ncl = cluless.length_C();

	cluless.calc_distance_matrix( -1 );

	cluless.print_centroids("centroids.pdb");
	// TRY:	ORDINARY ORDERED BEST FIT
	//	ALGINMENT NOT SHAPE FIT TO THE CENTROIDS

	particles pcs = cluless.get_centroids();
	richanalysis::particle_analysis pa;
	pa.assign_particles( pcs , c_aligned );
//	pa.remove_centroids();
	pa.density_model_hybrid( pcs, c_aligned );
	pa.print_model("realigned.pdb");
//
	return 0;
}
//
//// DEBUG SECTION:: SINGLE PASS
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
/*
	old_layer.push_back(n0);
	std::vector<int> clu_ndx,cen_ndx;
	while ( solved_layer.size() !=  vparts[1].size() )
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
				fIO.output_pdb("clu-nofit-n" + ns + ".pdb", vparts[1] , clu_ndx ); 
				fIO.output_pdb("den-nofit-n" + ns + ".pdb", vparts[0] , den_ndx );
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
*/

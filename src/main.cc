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
	poop.first = "-iden";
	opts.push_back(poop);
	poop.first = "-imod";
	opts.push_back(poop);
	poop.first = "-omod";
	opts.push_back(poop);
	poop.first = "-oclu";
	opts.push_back(poop);

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
				std::cout << opts[opt].first << " \t [string] " << "\n\t";
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
//
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

	if( L[0] >= L[1] ) {
		n0.first.set_matrix(  vparts[0] );
		n0.second.set_matrix( vparts[1] ); 
	}else{
		n0.first.set_matrix(  vparts[1] );
		n0.second.set_matrix( vparts[0] );
	}

	int sw = (L[1]>=L[0]) + 1;
	old_layer.push_back(n0);
	int C = 1;
/*
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
	std::cout << "::ALTERNATIVE::" << std::endl;
	n0.first.find_centroids();
	n0.second.find_centroids();
	richanalysis::node_analysis nnode(n0);
//	
//	particles c_aligned = nnode.regular_fit();		// shape fit
//	for( int i=0 ; i<c_aligned.size() ; i++ )
//		c_aligned[i].first = carth_space[i].first;
//	fIO.output_pdb("test-can" + ns + ".pdb",  c_aligned );

//	PARTICLE ANALYSIS
	int D = L[0] , B = L[1];
	std::cout << "INFO::PA:: " << D << " & " << B << std::endl;
	richanalysis::particle_analysis pa,ps;
	richanalysis::cluster cl4dg;

	int NN=D>B?D:B;
	int MM=D>B?B:D;
	cl4dg.alloc_space(NN,MM);
	if( L[0] > L[1] ) {
		cl4dg.set_matrix(vparts[0]);
	}
	else {
		cl4dg.set_matrix(vparts[1]);
	}

//	cl4dg.print_all();
	richanalysis::coord_format cfp;
	
	gsl_matrix *pos_d	= gsl_matrix_calloc( DIM, NN );
	gsl_matrix *pos_c	= gsl_matrix_calloc( DIM, MM );
	gsl_vector *sig		= gsl_vector_calloc( NN );
	gsl_vector *cnt		= gsl_vector_calloc( NN );
	gsl_vector *cid		= gsl_vector_calloc( NN );
	gsl_vector *rho2	= gsl_vector_calloc( NN );

	cl4dg.copyC(pos_c);		// centroids	
	cl4dg.copyM(pos_d);		// clustered	
 	cl4dg.copyv(cid);

  	for ( int j = 0; j < NN; j++ )
	{
		int id = gsl_vector_get(cid,j);
		gsl_vector_view column = gsl_matrix_column ( pos_d, j  );
		gsl_vector_view col_mu = gsl_matrix_column ( pos_c, id );
		double d,d0,cnt0;
		gsl_vector_sub( &column.vector, &col_mu.vector );

		d	= gsl_blas_dnrm2 ( &column.vector );
		d0	= gsl_vector_get (  sig , id	);
		cnt0	= gsl_vector_get (  cnt , id	);
		gsl_vector_set( sig , id, d0 + d	);
		gsl_vector_set( cnt , id, cnt0 + 1.0	);
	}

	double Er=0.0;
	for (int j = 0; j < MM; j++)
	{
		Er += gsl_vector_get(cnt,j);
//		std::cout	<< "INFO::STD::" << sqrt(gsl_vector_get(sig,j)/gsl_vector_get(cnt,j)) 
//				<< " " << gsl_vector_get(cnt,j) << std::endl;
	}
	std::cout << "INFO::EREST:: "<< Er/((float)MM) << std::endl;
	particles red_parts = cfp.mat2par( pos_c , sig );

	fIO.output_pdb("test-cents_" + ns + ".pdb",  red_parts );
	gvec *r1 = gsl_vector_alloc( DIM );
	gvec *r2 = gsl_vector_alloc( DIM );
	gmat *DISTM = gsl_matrix_calloc(MM,MM);

	double tbond=1.5;
	tbond*=tbond;
	double tprox=0.9;
	tprox*=tprox;

	for(int i=0;i<MM; i++){
		double b = gsl_vector_get(sig,i);
		b*=0.1;
		for(int j=0;j<MM;j++){
			gsl_matrix_get_col(r1,pos_c,i);
			gsl_matrix_get_col(r2,pos_c,j);
			gsl_vector_sub(r1,r2);
			ftyp len = gsl_blas_dnrm2(r1);
			len	*= len;
			//len	 = len<tbond?tbond:len;
			len	+= (i>j)?(len+b):(len-b);
			gsl_matrix_set( DISTM, i, j, len );
		}
	}

	particles new_parts;
	richanalysis::particle_analysis single_reduction;
	new_parts = single_reduction.assign_via_distmatrix(DISTM);

	if( L[0] > L[1] ) 
	{
		for(int i=0;i<MM;i++)
			new_parts[i].first = vparts[1][i].first ;
	}
	else
	{
		for(int i=0;i<MM;i++)
			new_parts[i].first = vparts[0][i].first ;
	}

	fIO.output_geometry( new_parts, "new_parts.xyz");

	return 0;
}

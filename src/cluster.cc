/*	cluster.cc	*/
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

#include <iostream>
#include "richtypes.h"
#include "cluster.hh"
#include <algorithm>
#include <random>

using namespace richanalysis;

void
particle_analysis::output_result( std::string filename ) {
	if( complete() ) {
		int N = model_.size();
		particles ps;
		for(int i=0 ; i<N ; i++ ) {
			gvec	*w = gsl_vector_alloc(DIM);
			gsl_matrix_get_col(w,C_,i);
			particle p;
			std::string label	=	"C" ;
			p.first			=	label;
			p.second		= 	gsl_vector_alloc(DIM);
			gsl_vector_memcpy(p.second,w);
			ps.push_back(p);
			gsl_vector_free(w);
		}
		output_pdb( filename , ps );
	}
}

bool particle_sort_func(std::pair<ftyp, std::pair< int, int > > i1,std::pair<ftyp, std::pair< int, int > > i2) {
	return (i1.first<i2.first);
};

particles
particle_analysis::assign_via_distmatrix( gmat *A ) {

	int dspace	= DIM;
	int N		= A->size1;
	gsl_matrix *X	= gsl_matrix_calloc( N , dspace );
	gsl_matrix *Y	= gsl_matrix_calloc( dspace, N );
	gsl_vector *b   = gsl_vector_calloc( N );
	gsl_vector_set_all( b , 1.0 );

	if( !complete() && !single() && (b->size==A->size1) )
	{
		if(A->size1!=A->size2) {
			std::cout << "INFO::ERROR::DIMENSIONS" << std::endl;
		}

		A_	= gsl_matrix_calloc( A->size1 , A->size2 );
		gsl_matrix_memcpy( A_ , A );
		B_	= gsl_matrix_calloc( A->size1 , A->size2 );
		gsl_matrix_memcpy( B_ , A );

		gsl_matrix *D	= gsl_matrix_calloc( N , N );
		gsl_matrix_memcpy( D  , A );

		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				double dij	= gsl_matrix_get(D,i,j);
				double dNj	= gsl_matrix_get(D,N-1,j);
				double diN	= gsl_matrix_get(D,i,N-1);
				double nDij	= (diN+dNj-dij)*0.5;
				gsl_matrix_set ( D, i, j , nDij );
			}
		}

		gsl_matrix *U	= gsl_matrix_alloc( N, N );
		gsl_matrix *V	= gsl_matrix_alloc( N, N );
		gsl_vector *S	= gsl_vector_alloc( N );
		gsl_vector *wrk = gsl_vector_alloc( N );
		gsl_linalg_SV_decomp ( D, V, S, wrk ); 		
		gsl_matrix_memcpy(  U, D );
	
		gsl_matrix *E	= gsl_matrix_calloc( N, dspace );
		for(int i=0;i<dspace;i++)
			gsl_matrix_set( E, i, i, sqrt(gsl_vector_get(S,i)) );

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
				1.0, U, E, 0.0, X );

		C_	= gsl_matrix_calloc( dspace, N );
		vc_	= gsl_vector_calloc( N );
		gsl_vector_set_all( vc_ ,  1.0 );

		M_	= gsl_matrix_calloc( dspace, N );
		wc_	= gsl_vector_calloc( N );
		gsl_vector_set_all( wc_ ,  1.0 );

		gsl_matrix_transpose_memcpy( C_, X );
		gsl_matrix_transpose_memcpy( M_, X );

		bSingleSet_ = true;
	}
	else 
	{
		std::cout << "INFO::ERROR::ALREADY::COMPLETED" << std::endl;
	}
	gsl_matrix_transpose_memcpy( Y, X );
	std::cout << "INFO::HAVE " << X->size1 << " AND " << X->size2 << std::endl;
	richanalysis::coord_format cfp;
	particles red_parts = cfp.mat2par( Y , b );
	std::cout << "INFO::HAVE " << red_parts.size() << std::endl;

	return red_parts;

}


particles
particle_analysis::assign_via_distmatrix( gmat *A , gvec *b ) {

	int dspace	= DIM;
	int N		= A->size1;
	gsl_matrix *X	= gsl_matrix_calloc( N , dspace );

	if( !complete() && !single() && (b->size==A->size1) )
	{
		int N	= A->size1;
		if(A->size1!=A->size2) {
			std::cout << "INFO::ERROR::DIMENSIONS" << std::endl;
		}

		A_	= gsl_matrix_calloc( A->size1 , A->size2 );
		gsl_matrix_memcpy( A_ , A );
		B_	= gsl_matrix_calloc( A->size1 , A->size2 );
		gsl_matrix_memcpy( B_ , A );

		gsl_matrix *D	= gsl_matrix_calloc( N , N );
		gsl_matrix_memcpy( D  , A );

		// ZERO CHECKING NOT IN ORIGINAL ALGORITHM
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				if( i==j )
					continue;
				double BF 	 = gsl_vector_get(b,i) * gsl_vector_get(b,j) ; 	// s2
				double Aval	 = gsl_matrix_get(A,i,j);
				Aval		*= Aval; 					// d2
				if( i>j )
					Aval += BF;
				if( i<j )
					Aval -= BF;
				gsl_matrix_set( D,i,j, (Aval>=0)?(Aval):(0.0) );
			}
		}

		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				double dij	= gsl_matrix_get(D,i,j);
				double dNj	= gsl_matrix_get(D,N-1,j);
				double diN	= gsl_matrix_get(D,i,N-1);
				double nDij	= (diN+dNj-dij)*0.5;
				gsl_matrix_set ( D, i, j , nDij );
			}
		}

		gsl_matrix *U	= gsl_matrix_alloc( N, N );
		gsl_matrix *V	= gsl_matrix_alloc( N, N );
		gsl_vector *S	= gsl_vector_alloc( N );
		gsl_vector *wrk = gsl_vector_alloc( N );
		gsl_linalg_SV_decomp ( D, V, S, wrk ); 		
		gsl_matrix_memcpy(  U, D );
		
		gsl_matrix *E	= gsl_matrix_calloc( N, dspace );
		for(int i=0;i<dspace;i++)
			gsl_matrix_set( E, i, i, sqrt(gsl_vector_get(S,i)) );

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
				1.0, U, E, 0.0, X );

		C_	= gsl_matrix_calloc( dspace, N );
		vc_	= gsl_vector_calloc( N );
		gsl_vector_set_all( vc_ ,  1.0 );

		M_	= gsl_matrix_calloc( dspace, N );
		wc_	= gsl_vector_calloc( N );
		gsl_vector_set_all( wc_ ,  1.0 );

		gsl_matrix_transpose_memcpy( C_, X );
		gsl_matrix_transpose_memcpy( M_, X );

		bSingleSet_ = true;
	}
	else 
	{
		std::cout << "INFO::ERROR::ALREADY::COMPLETED" << std::endl;
	}

	richanalysis::coord_format cfp;
	particles red_parts = cfp.mat2par( X , b );
	return red_parts;

}

void
particle_analysis::print_model( std::string filename ) {
		if( bAssigned_ ) {
			richanalysis::fileIO fIO;
			fIO.output_pdb( filename , model_ );
		}
}

void
particle_analysis::print_density( std::string filename ) {
		if( bAssigned_ ) {
			richanalysis::fileIO fIO;
			fIO.output_pdb( filename , density_ );
		}
}

void
particle_analysis::assign_matrices(void) {
	std::vector< std::pair<ftyp, std::pair< int, int > > > vi;
	rel_idx_.clear();
	
	if( complete() ) {
		int N		= model_.size();
		A_		= gsl_matrix_calloc(model_.size(),model_.size());
		B_ 		= gsl_matrix_calloc(model_.size(),model_.size());
		C_		= gsl_matrix_calloc(DIM,model_.size());
		gmat *P_	= gsl_matrix_calloc(DIM,model_.size());
		M_		= gsl_matrix_calloc(DIM,density_.size());
		vc_		= gsl_vector_calloc(density_.size());
		wc_		= gsl_vector_calloc(model_.size());
		gvec *r		= gsl_vector_calloc( DIM );
		gvec *rc	= gsl_vector_calloc( DIM );
		cen_m_		= gsl_vector_calloc( DIM );
		cen_d_		= gsl_vector_calloc( DIM );
		rc	= gsl_vector_calloc( DIM );
		// IF SAME SIZE THEN DONT CLUSTER
		if( !( density_.size() == model_.size() ) ){
			for(int i=0;i<density_.size();i++) { 
				gsl_vector_memcpy( r , density_[i].second ); 
				gsl_matrix_set_col(M_, i, r);
			}	
			gsl_kmeans( M_ , vc_, C_, wc_ );
		}else{
			gsl_matrix_memcpy( C_ , M_  );
			gsl_vector_memcpy( wc_, vc_ );
		}
		
		// REGULAR ALIGNMENT TO
		richanalysis::fitting rfit;
		U_ 	= gsl_matrix_calloc( DIM , DIM );
		t_ 	= gsl_vector_calloc( DIM );
		float diff	= rfit.kabsch_fit( P_ , C_ , U_ , t_ );
		bHaveUt_	= true;
		
		rfit.apply_fit	( model_ , U_ , t_ );

		for(int i=0;i<N;i++) {
			double min_d = 1e10;
			gsl_vector_memcpy( r , model_[i].second ); 
			// gsl_matrix_set_col( P_, i, r );
			for( int j=0; j<N; j++ ) {
				gsl_vector_memcpy 	( r , model_[i].second ); 
				gsl_matrix_get_col	( rc, C_ , j ); 
				gsl_vector_sub	  	( rc, r );
				double len	= 	gsl_blas_dnrm2( rc );
				len		= 	len*len;
				if( len < min_d ) {
					std::pair<ftyp, std::pair<int, int> > rel;
					rel.first = len; rel.second.first=i; rel.second.second=j;
					vi.push_back(rel);
					min_d=len;
				}
			}
		}
		calc_distance_matrices();
//		output_matrix(C_);
//		output_matrix(P_);
		bMatrices_	= true;
	} else {
		std::cout <<"INFO::IHAVENOTHING!" <<std::endl;
	}
/*
	std::sort (vi.begin(), vi.end(), particle_sort_func);
	std::vector<int> vi1, vi2, vi_fin;
	for( int i=0 ; i<model_.size() ; i++ ) {
		vi1.push_back(1);
		vi2.push_back(1);
		rel_idx_.push_back(-1);
	}
	int i1,i2;
	for( int i=0 ; i<vi.size() ; i++ ) {
		i1 = vi[i].second.first ;
		i2 = vi[i].second.second;
		if( rel_idx_[i2] < 0 && vi1[i1] ) {
			rel_idx_[i2]	= i1;
			vi1[i1]		=  0;
			vi2[i2]		=  0;
		}
	}
*/
}

int 			
particle_analysis::calc_distance_matrix( gmat *D, particles p ){
	if( D->size1==p.size() && D->size1==D->size2 ) {
		gvec *r1	= gsl_vector_calloc( DIM );
		gvec *r2	= gsl_vector_calloc( DIM );
		for( int i=0 ; i<p.size() ; i++ ) {
			gsl_vector_memcpy( r1 , p[i].second ); 
			for( int j=0 ; j<p.size() ; j++ ) {
				gsl_vector_memcpy( r2 , p[j].second );
				gsl_vector_sub(r2,r1);
				ftyp len = gsl_blas_dnrm2(r2);
				gsl_matrix_set( D, i, j, len );
			}
		}
		gsl_vector_free(r1);
		gsl_vector_free(r2);
	} else {
		std::cout << "DIMENSION ERROR" << std::endl;
	}
	return 0;
}

void
particle_analysis::remove_centroids() {
	if( model_.size()>0 && density_.size()>0 ) {
		gsl_vector *rs = gsl_vector_calloc(DIM);
		gsl_vector *rc = gsl_vector_calloc(DIM);
		gsl_vector *rm = gsl_vector_calloc(DIM);
		gsl_vector *rd = gsl_vector_calloc(DIM);

		for( int i=0 ; i<model_.size() ; i++ ) {
			gsl_vector_memcpy	(rm , model_[i].second);
			gsl_vector_scale	(rm , 1.0/model_.size());
			gsl_vector_add		(rs , rm);
		}
		for( int i=0 ; i<model_.size() ; i++ ) 
			gsl_vector_sub( model_[i].second	, rs );

		for( int i=0 ; i<density_.size() ; i++ ) {
			gsl_vector_memcpy	(rd , density_[i].second);
			gsl_vector_scale	(rd , 1.0/density_.size());
			gsl_vector_add		(rc , rd);
		}
		for( int i=0 ; i<density_.size() ; i++ ) 
			gsl_vector_sub( density_[i].second	, rc );

		gsl_vector_memcpy( cen_m_ , rs );
		gsl_vector_memcpy( cen_d_ , rc );
	}else{
		std::cout << "AINT DOING JACK" << std::endl;
	}
}

void
particle_analysis::density_model_integer_run( particles pd , particles pm ) 
{
	std::mt19937 generator(1871865);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	double dice_roll = distribution(generator);
	if( pd.size() == pm.size() ) {
		gsl_matrix *D0 = gsl_matrix_calloc(pm.size(),pm.size());
		gsl_matrix *D1 = gsl_matrix_calloc(pm.size(),pm.size());
		gsl_matrix *D2 = gsl_matrix_calloc(pm.size(),pm.size());
		copyA(D0);
		gsl_vector *rm = gsl_vector_calloc(DIM);
 		gsl_vector *rd = gsl_vector_calloc(DIM);
		gsl_vector *rk = gsl_vector_calloc(DIM);
		gsl_vector *rr = gsl_vector_calloc(DIM);

		calc_distance_matrix( D1 , pd ); 
		calc_distance_matrix( D2 , pm );

		ftyp beta	= 1.0E-4;

		int N = (int)( pd.size() );
		gsl_vector *d1i = gsl_vector_calloc(N);
		gsl_vector *d0i = gsl_vector_calloc(N);
		gsl_vector *d1j = gsl_vector_calloc(N);
		gsl_vector *d0j = gsl_vector_calloc(N);
		int nsw=1e3;
		for(int isw=0;isw<nsw;isw++)
		{
			beta += 1E-4*((float)isw/nsw);
			int I = (int)(distribution(generator)*(pd.size()-1.0) );
			int J = (int)(distribution(generator)*(pd.size()-1.0) );
			J = (I==J)?(I+1)>=N?(0):(I+1):(J);
			//std::cout << "INFO:: " << I << " " << J << " " << N << std::endl;
			particle ptmp=pd[I];
			pd[I]=pd[J];pd[J]=ptmp;
			calc_distance_matrix( 	D1 , pd 	);
			gsl_matrix_get_row(d1i , D1 , I);
			gsl_matrix_get_row(d0i , D0 , I);
			gsl_matrix_get_row(d1j , D1 , J);
			gsl_matrix_get_row(d0j , D0 , J);	

			gsl_vector_sub	 ( d1i , d0i  );
			gsl_vector_sub	 ( d1j , d0j  );
			ftyp lenj = gsl_blas_dnrm2( d1j );
			ftyp leni = gsl_blas_dnrm2( d1i );
			leni*=leni;
			lenj*=lenj;
			
			if(distribution(generator)>exp(-beta*(leni + lenj) )){
				particle ptmp=pd[I];
				pd[I]=pd[J];pd[J]=ptmp;		
			}
		}
		std::cout << "DONE" << std::endl;

		gsl_vector_free(d1i);
		gsl_vector_free(d0i);
		gsl_vector_free(d0j);
		gsl_vector_free(d1j);
		gsl_vector_free(rm);	gsl_vector_free(rd);
		gsl_vector_free(rk);	gsl_vector_free(rr);
		gsl_matrix_free(D0);	gsl_matrix_free(D1);

		for(int i=0;i<pd.size();i++) {
			pd[i].first=pm[i].first;
			density_[i]=pd[i];
		}
		
	}	
}

void
particle_analysis::density_model_hybrid( particles pd , particles pm ) 
{

	std::mt19937 generator(1871865);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	double dice_roll = distribution(generator);

//  HOW TO ROLL THE DICE
//  dice_roll = distribution(generator);

	if( pd.size() == pm.size() ) {
		gsl_matrix *D0 = gsl_matrix_calloc(pm.size(),pm.size());
		gsl_matrix *D1 = gsl_matrix_calloc(pm.size(),pm.size());
		gsl_matrix *D2 = gsl_matrix_calloc(pm.size(),pm.size());
		copyA(D0);
		gsl_vector *rm = gsl_vector_calloc(DIM);
 		gsl_vector *rd = gsl_vector_calloc(DIM);
		gsl_vector *rk = gsl_vector_calloc(DIM);
		gsl_vector *rr = gsl_vector_calloc(DIM);

		calc_distance_matrix( D1 , pd ); 
		ftyp cutoff 	= 5.0;
		ftyp dr 	= 5.0E-1;
		ftyp beta	= 1.0E-3;
		int NSWEEPS	= 1E4;
		for( int ISWEEP=0 ; ISWEEP<NSWEEPS ; ISWEEP++ ) {
			ftyp fsweep	= ISWEEP/((float)NSWEEPS);
			ftyp s2m	= cutoff*cutoff;
			ftyp invgau_m	= 1.0/sqrt(2*s2m*M_PI);
			ftyp s2d	= cutoff*cutoff*2.0*2.0; 
			ftyp invgau_d	= 1.0/sqrt(2*s2d*M_PI);

			// FINDS NEAREST CENTROIDS ( CAN BE DEGENERATE )
			for( int i=0 ; i<pm.size() ; i++ ) {
				ftyp len_match	= 1.0E6;
				int	ikeep	= 0;
				gsl_vector_memcpy( rm , pm[i].second );
				for( int j=0 ; j<pd.size() ; j++ ) {
					gsl_vector_memcpy( rd , pd[j].second );
					gsl_vector_sub	 ( rd , rm  );
					ftyp len = gsl_blas_dnrm2(rd);
					if( len < len_match ) {
						ikeep = j;
						len_match	= len;
						ftyp ds		= dr>len?dr:len;	
						gsl_vector_memcpy( rk , rd );
						gsl_vector_set	(  rr , XX , (2.0*distribution(generator)-1.0)*ds );
						gsl_vector_set	(  rr , YY , (2.0*distribution(generator)-1.0)*ds );
						gsl_vector_set	(  rr , ZZ , (2.0*distribution(generator)-1.0)*ds );
						gsl_vector_add  (  rk , rr );
					}
				}
				gsl_vector_add  (  rk , rm );
				ftyp sum0 = 0.0, sum1 = 0.0, sum2=0.0;
				for( int k=0 ; k<pm.size() ; k++ ) { // SELF ENERGY
					if( k!=i ) { 
						gsl_vector_memcpy ( rm , pm[k].second );
						gsl_vector_sub	  ( rm , rk );
						ftyp len0	= gsl_blas_dnrm2(rm); // CORRESPONDS TO ONE ROW OF D0
						ftyp ilen0	= 1.0/len0;
						ftyp d0		= gsl_matrix_get(D0,i,k);
						len0 -= d0; len0 *= len0;
						sum0 -= exp(-len0*0.5/s2m )/d0/d0/d0; 
					}
				}
				sum0*=(fsweep+1.0)*4.0;
				for( int k=0 ; k<pd.size() ; k++ ) {
					gsl_vector_memcpy ( rm , pd[k].second );
					gsl_vector_memcpy ( rr , pd[k].second );
					gsl_vector_sub	  ( rr , pm[k].second );
					gsl_vector_sub	  ( rm , rk );
					ftyp len1 = gsl_blas_dnrm2(rm);
					ftyp len2 = gsl_blas_dnrm2(rr);
					len1 *= len1; len2 *= len2;
					sum1-=exp(-len1*0.5/s2d )*invgau_d; 
					sum2-=exp(-len2*0.5/s2d )*invgau_d; 
				}
				ftyp x	 = (1.05-fsweep);
				sum1	*= x*x;
				sum2	*= x*x;
				if( distribution(generator) <= exp( -1.0*( (sum1+sum0) - sum2) )*beta )
					gsl_vector_memcpy	( pm[i].second , rk );
/*
				ftyp summ = 0.0, sumd = 0.0;
				for( int k=0 ; k<pm.size() ; k++ ) {
					ftyp c1  = gsl_matrix_get(D0,i,k);
					ftyp c2  = gsl_matrix_get(D1,ikeep,k);
					summ	+= (c1<cutoff?c1:cutoff) ;
					sumd	+= (c2<cutoff?c2:cutoff) ;
				}

				gsl_vector_memcpy( rm , pm[  i  ].second );
				gsl_vector_memcpy( rd , pd[ikeep].second );
				ftyp delta 	 = abs(summ-sumd);
				gsl_vector_scale	( rd , 5.0E-4*exp(-delta) ); //exp(-delta*1.0E-8) );// exp(-delta)
				gsl_vector_add		( rm , rd );
*/

			}
		}
		gsl_vector_free(rm);	gsl_vector_free(rd);
		gsl_vector_free(rk);	gsl_vector_free(rr);
		gsl_matrix_free(D0);	gsl_matrix_free(D1);
	} else {
		std::cout << "CANNOT CREATE HYBRID MODEL" << std::endl;
	}

	for( int i=0 ; i<pm.size() ; i++ ) {
		gsl_vector_memcpy( model_[i].second , pm[i].second); 
	}
}

int 
cluster::calc_distance_matrix(int i) {
	if( isSet() ) {
		//if( A_!=NULL  ) 
		//	gsl_matrix_free(A_);

		int D	= ( i>=0 ) ? (M_->size2) : (C_->size2);
		A_ 	= gsl_matrix_calloc( D, D );
		gsl_matrix *MAT = gsl_matrix_calloc( DIM , D );

		if(i>=0)
			gsl_matrix_memcpy( MAT , M_ );
		else
			gsl_matrix_memcpy( MAT , C_ );

		gvec *r1 = gsl_vector_calloc( DIM );
		gvec *r2 = gsl_vector_calloc( DIM );

		ftyp xmin, x2min, X, X2;
		for(int i=0;i<D; i++){
			xmin=1E10;
			for(int j=0;j<D;j++){
				gsl_matrix_get_col(r1,MAT,i);
				gsl_matrix_get_col(r2,MAT,j);
				gsl_vector_sub(r1,r2);
				ftyp len = gsl_blas_dnrm2(r1);
				gsl_matrix_set( A_, i, j, len );
				if(i!=j && len<xmin) {
					xmin	= len;
					x2min	= len*len;
				}
			}
			X += xmin; X2 += x2min;
		}
		ftyp n   = A_->size2;
		maxbond_ = X/n+sqrt( (X2-X*X/n)/n );
		return 0;
	}else{
		std::cout << "ERROR::PARTICLE SETS NOT ASSIGNED" << std::endl;
		return 1;
	}
}


int 
particle_analysis::calc_distance_matrices(void) {

	if( complete() ) {
		if(	A_	) // CLUSTERED DENSITY
			gsl_matrix_free(A_);
		A_  = gsl_matrix_alloc( C_->size2, C_->size2 );

		gvec *r1 = gsl_vector_alloc( DIM );
		gvec *r2 = gsl_vector_alloc( DIM );
		
		int D = A_->size2;

		for(int i=0;i<D; i++){
			for(int j=0;j<D;j++){
				gsl_matrix_get_col(r1,C_,i);
				gsl_matrix_get_col(r2,C_,j);
				gsl_vector_sub(r1,r2);
				ftyp len = gsl_blas_dnrm2(r1);
				gsl_matrix_set( A_, i, j, len );
			}
		}

		if(B_) // MODEL
			gsl_matrix_free(B_);
		B_  = gsl_matrix_alloc( C_->size2, C_->size2 );

		for(int i=0;i<D; i++){
			for(int j=0;j<D;j++){
				gsl_vector_memcpy(r1,model_[i].second);
				gsl_vector_memcpy(r2,model_[j].second);
				gsl_vector_sub(r1,r2);
				ftyp len = gsl_blas_dnrm2(r1);
				gsl_matrix_set( B_, i, j, len );
			}
		}
		return 0;
	}else{
		std::cout << "ERROR::PARTICLE SETS NOT ASSIGNED" << std::endl;
		return 1;
	}
}

std::vector<int>
particle_analysis::outp_distance_matrix( gmat *A ) {
	int D = A->size1;
	ftyp level = 1.5;
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
					std::cout << "\033[1;31m " << val <<"\033[0m";
				else
					std::cout << " " << val ;
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
particle_analysis::outp_distance_matrix( gmat *A, ftyp level ) {
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

std::vector< std::pair<ftyp, std::pair< int, int > > >
particle_analysis::compare_dist_matrices(gmat *A, gmat *B, ftyp val) {

	std::vector< std::pair<ftyp, std::pair< int, int > > > vi;

	if(A->size1==B->size1&&A->size2==B->size2&&A->size1==B->size2)
	{
		gvec *va = gsl_vector_calloc(A->size1);
		gvec *vb = gsl_vector_calloc(B->size2);
		for(int i=0;i<A->size1;i++) 
		{
			ftyp diff=1E10,dotn2;
			int I=0;
			for(int j=0;j<B->size1;j++) 
			{
				gsl_matrix_get_row ( va, A, i );
				gsl_matrix_get_row ( vb, B, j );
				gsl_vector_sub( va, vb );
				diff=gsl_blas_dnrm2(va);
				std::pair<ftyp, std::pair<int, int> > rel;
				rel.first = diff; rel.second.first=i; rel.second.second=j;
				vi.push_back(rel);
			}
		}
	}
	std::sort (vi.begin(), vi.end(), particle_sort_func);
	return vi;
}

/*
			void			outp_distance_matrix(gmat *A);
void outp_distance_matrix(gmat *A) {
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
}
*/

std::vector<int > 
particle_analysis::find_via_distance( gmat *A, ftyp level ) {
	int D = A->size1;
	std::vector<int > is_tmp;
	int sumZero = (level<0)?(1):(0);
	if (sumZero)
		level = sqrt(level*level);
	std::cout <<"INFO::FIND VIA DISTANCE::LEVEL "<< level << std::endl;
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

cluster::cluster(particles m, particles c){
	int D1=m.size(),D2=c.size();
	alloc_space( D1 , D2 );

	for(int i=0;i<D1; i++) {
			gsl_matrix_set_col ( M_ , i, m[i].second ) ;
	}
	for(int i=0;i<D2; i++) {
			gsl_matrix_set_col ( C_ , i, c[i].second ) ;
			centroid_names_.push_back(c[i].first);
	}
	gsl_kmeans0( M_ , vc_ , C_ , wc_ );
}

void
cluster::print_centroids( std::string filename ){
		int N = C_->size2;
		particles ps;
		for(int i=0 ; i<N ; i++ ) {
			gvec	*w = gsl_vector_alloc(DIM);
			gsl_matrix_get_col(w,C_,i);
			particle p;
			if(centroid_names_.size()==C_->size2) {
				std::string label	=	centroid_names_[i];
				p.first			=	label;
			}
			else {
				std::string label	=	"C" ;
				p.first			=	label;
			}
			p.second		= 	gsl_vector_alloc(DIM);
			gsl_vector_memcpy(p.second,w);
			ps.push_back(p);
			gsl_vector_free(w);
		}
		richanalysis::fileIO fIO;
		fIO.output_pdb( filename , ps );
}


void
cluster::alloc_space( int D1, int D2 ) {
//D1	COORDINATE	SPACE
//D2	CENTROID	SPACE
//	std::cout << "ALLOC:: " << D1 << " " << D2 << std::endl;
	if( D1 > 0 && D2 > 0 ) {
		M_	= gsl_matrix_calloc( DIM, D1 );	
		C_	= gsl_matrix_calloc( DIM, D2 );	 
		vc_	= gsl_vector_calloc( D1 );
		wc_	= gsl_vector_calloc( D2 );
		ws_	= gsl_vector_calloc( D2 ); // centroid sigma
		bSet_   = true;
	} else {
		std::cout << "ERROR IN ALLOCATION ROUTINE:: "<< D1 << ", " << D2 << std::endl;
	}
}

void
cluster::realloc_centroids( int D2 ) {
	if( D2 > 0 && bSet_ ) {
		gsl_matrix_free(C_ );
		gsl_vector_free(wc_);
		gsl_vector_free(ws_);
		C_	= gsl_matrix_calloc( DIM, D2 );	 
		wc_	= gsl_vector_calloc( D2 );
		ws_	= gsl_vector_calloc( D2 ); // centroid sigma
	}
}

int 
cluster::set_matrix( particles coord_space ) {
	int D		= coord_space.size();
	int reval	= -1;

	int rD = get_cDIM(); 

	if( bSet_ == 0 ) {
//		std::cout << "INFO::ALLOCATING "<<D<<" "<<rD<<std::endl;
		alloc_space( D, rD );
	}

	if( M_->size2 == coord_space.size() && M_->size1 == DIM && bSet_ ) {
		// std::cout << "INFO::ASSIGNING::COORDS "<< D <<" "<<rD<<std::endl;
		id inInfo;
		idlabels_.clear();
		for(int i=0;i<D; i++)	{
			inInfo.first	= i;
			inInfo.second	= coord_space[i].first;
			idlabels_.push_back ( inInfo );
			gsl_matrix_set_col  ( M_ , i, coord_space[i].second );
		}
		reval =0; 
	} else {
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: set_matrix" 
			<< D << " " << rD << " " << M_->size2 << std::endl;
	}
	if(rD>1)
		find_centroids();

	return  reval;
}

bool sort_func(std::pair<ftyp, std::pair< int, int > > i1,std::pair<ftyp, std::pair< int, int > > i2) {
	return (i1.first<i2.first);
};

std::vector< int >
node_analysis::find_centroid_distance_relation( void ) {
	int verbose=0;
	gmat *C1 = gsl_matrix_calloc(	DIM, parents_.first.length_C()	);
	gmat *C2 = gsl_matrix_calloc(	DIM, parents_.second.length_C()	);
	gvec *ci = gsl_vector_calloc(	parents_.second.length_M()	);

	parents_.second	.copyv(ci);
	parents_.first	.copyC(C1);
	parents_.second	.copyC(C2);

	std::vector< std::pair<ftyp, std::pair< int, int > > > vi;

	if(C1->size2!=C2->size2)
		std::cout << "ERROR::ERROR::WILL::MALFUNCTION::" << std::endl;

	gvec *va  = gsl_vector_calloc( DIM );
	gvec *vb  = gsl_vector_calloc( DIM );
	gvec *v1c = gsl_vector_calloc( DIM );
	gvec *v2c = gsl_vector_calloc( DIM );

	for(int i=0 ; i<C1->size2 ; i++ ) {
		ftyp diff=1E10,dotn2;
		for( int j=0 ; j<C2->size2 ; j++ ) {
			gsl_matrix_get_col ( va, C1, i );
			gsl_matrix_get_col ( vb, C2, j );
			gsl_vector_sub( va, vb );
			diff = gsl_blas_dnrm2(va);
			std::pair<ftyp, std::pair<int, int> > rel;
			rel.first = diff*diff; rel.second.first=i; rel.second.second=j;
			vi.push_back(rel);
		}
	}
	std::sort (vi.begin(), vi.end(), sort_func);

	if(verbose)
		for( int i=0 ; i<vi.size() ; i++ )
			std::cout << "INFO:: \t" << i << " \t " 
			<< vi[i].second.first << "   " 
			<< vi[i].second.second << std::endl; 

	std::vector< std::pair<ftyp, std::pair< int, int > > > vir;
	std::vector<int> vi_fin, vi0;
	while( vir.size() < C1->size2 ) {
		for( int i=0 ; i<vi.size() ; i++ ) {
			bool exists=false;
			int I = vi[i].second.first, J = vi[i].second.second;
			for( int j=0; j<vir.size(); j++)
				if( vir[j].second.first == I || vir[j].second.second == J)
					exists=true;
			if(!exists) {
				vir.push_back(vi[i]);
				vi_fin.push_back(-1);
			}
		}
	}

	for( int i=0 ; i<vir.size() ; i++ ){
		if(verbose)
			std::cout << "FOUND:: \t" << i << " \t " 
			<< vir[i].second.first << "   " 
			<< vir[i].second.second << std::endl; 
		vi_fin[vir[i].second.second]=vir[i].second.first;
	}

	MDidx_.clear();
	for(int i=0; i<vi_fin.size() ;i++){
		MDidx_.push_back(vi_fin[i]);
	}
	bMDlabels_=true;

	return vi_fin;
};

bool sort_f00(std::pair< int, int >  i1, std::pair< int, int > i2) {
	return (i1.second < i2.second);
};

particles
node_analysis::centroid_frag_fit( void ) {
	std::vector< int > cd_idx	=  find_centroid_distance_relation();
//	std::vector< int > cd_idx	=  find_centroid_relation();
	gmat *M1 = gsl_matrix_calloc( DIM, parents_.first .length_M() );
	gmat *M2 = gsl_matrix_calloc( DIM, parents_.second.length_M() );
	gvec *v1 = gsl_vector_alloc ( parents_.first .length_M() );
	gvec *v2 = gsl_vector_alloc ( parents_.second.length_M() );
	parents_.second.copyv(v2);
	parents_.first .copyv(v1);
	parents_.second.copyM(M2);
	parents_.first .copyM(M1);

	std::vector< std::pair<int,int> > i2c, c2i;
	for( int i=0 ; i<parents_.second.length_M() ; i++ ){
		std::pair<int,int> ip;
		ip.first = i; ip.second = gsl_vector_get(v2,i);
		i2c.push_back(ip);
	}

	for( int ic=0 ; ic<cd_idx.size() ; ic++ )
		for( int i=0 ; i < i2c.size() ; i++ )
			if(i2c[i].second==ic)
				c2i.push_back(i2c[i]);

	/*std::cout << " I2C :: " << i2c.size() << " ][ ";
	for( int ic=0 ; ic < i2c.size() ; ic++ )
		std::cout << " ( " << i2c[ic].first << ", " << i2c[ic].second << " ) " ;
	std::cout << std::endl;
	std::cout << " C2I :: " << c2i.size() << " ][ ";
	for( int ic=0 ; ic < c2i.size() ; ic++ )
		std::cout << " ( " << c2i[ic].first << ", " << c2i[ic].second << " ) " ;
	std::cout << std::endl;*/

	gmat *U 	= gsl_matrix_calloc( DIM, DIM );
	gmat *iU 	= gsl_matrix_calloc( DIM, DIM );
	gvec *it 	= gsl_vector_calloc( DIM );
	gvec *t 	= gsl_vector_calloc( DIM );

	particles parts;
	richanalysis::coord_format cf;
	for( int ic=0 ; ic<cd_idx.size() ; ic++ ){
		particles dp	= cf.truncmat(M1,v1,cd_idx[ic]);
 		particles mp	= cf.truncmat(M2,v2,ic);
		
		gsl_matrix *P	= gsl_matrix_alloc(DIM,mp.size());
		gsl_matrix *Q	= gsl_matrix_alloc(DIM,dp.size());
		for(int i=0;i<dp.size();i++) {
			gsl_matrix_set_col(Q,i,dp[i].second);
		}
		for(int i=0;i<mp.size();i++) {
			gsl_matrix_set_col(P,i,mp[i].second);
		}
		ftyp rmsd	= shape_fit( P, Q, U, t, 1 ); 
		std::cout << "INFO::RMSD::" << rmsd << std::endl;
		invert_fit	( U, t, iU, it );
		apply_fit	( mp, iU, it ); 
		parts.insert( parts.end() , mp.begin() , mp.end() );
		gsl_matrix_free(P);
		gsl_matrix_free(Q);		
	}

	particles o_parts;
	for(int i=0;i<parents_.second.length_M();i++){
		o_parts.push_back( parts[c2i[i].first] );		
	}

	gsl_matrix_free( U );
	gsl_vector_free( t );
	gsl_matrix_free(iU );
	gsl_vector_free(it );
	gsl_matrix_free(M1 );
	gsl_matrix_free(M2 );
	gsl_vector_free(v1 );
	gsl_vector_free(v2 );


	return o_parts;
};

particles
node_analysis::centroid_to_nearest(void) {
	gmat *C1 = gsl_matrix_calloc(DIM, parents_.first.length_C() );
	gmat *C2 = gsl_matrix_calloc(DIM, parents_.second.length_C() );
	gvec *ci = gsl_vector_calloc(parents_.second.length_M());

	parents_.second.copyv(ci);

	parents_.first.copyC(C1);
	parents_.second.copyC(C2);

	std::vector< std::pair<ftyp, std::pair< int, int > > > vi;

	if(C1->size2!=C2->size2)
		std::cout << "ERROR::ERROR::WILL::MALFUNCTION::" << std::endl;

	gvec *va  = gsl_vector_calloc( DIM );
	gvec *vb  = gsl_vector_calloc( DIM );
	gvec *v1c = gsl_vector_calloc( DIM );
	gvec *v2c = gsl_vector_calloc( DIM );

	for(int i=0 ; i<C1->size2 ; i++ ) {
		ftyp diff=1E10,dotn2;
		for( int j=0 ; j<C2->size2 ; j++ ) {
			gsl_matrix_get_col ( va, C1, i );
			gsl_matrix_get_col ( vb, C2, j );
			gsl_vector_sub( va, vb );
			diff=gsl_blas_dnrm2(va);
			std::pair<ftyp, std::pair<int, int> > rel;
			rel.first = diff*diff; rel.second.first=i; rel.second.second=j;
			vi.push_back(rel);
		}
	}
	std::sort (vi.begin(), vi.end(), sort_func);

	std::vector<int> vi1, vi2, vi_fin;
	for( int i=0 ; i<C1->size2 ; i++ ) {
		vi1.push_back(1);
		vi2.push_back(1);
		vi_fin.push_back(-1);
	}
	int i1,i2;
	for( int i=0 ; i<vi.size() ; i++ ) {
		i1 = vi[i].second.first ;
		i2 = vi[i].second.second;
		if( vi_fin[i2] < 0 && vi1[i1] ) {
			vi_fin[i2]	= i1;
			vi1[i1]		=  0;
			vi2[i2]		=  0;
		}
	}

	gsl_matrix *P	= gsl_matrix_calloc(DIM,parents_.second.length_M());
	gsl_vector *v	= gsl_vector_calloc(parents_.second.length_M());

	parents_.second.copyM(P);
	parents_.second.copyv(v);

	richanalysis::coord_format cf;
//	for( int i=0 ; i<vi_fin.size() ; i++ )
//		std::cout << "INFO:: " << vi_fin[i] << std::endl;
	particles pf;
	pf 		= cf.mat2par ( P, v );
	for( int i=0 ; i<pf.size() ; i++ ){
		pf[i].first = "C";
		int ci	= gsl_vector_get(v,i);
		int ti	= vi_fin[ci];
	}
//// yeah... I know...
	for(int i=0 ; i<C1->size2 ; i++ ) {
		ftyp diff=1E10;
		gsl_matrix_get_col	( va, C1, vi_fin[i] );
		gsl_matrix_get_col	( vb, C2, i );
		gsl_vector_sub		( va, vb );
		// HERE GO THROUGH THE ATOMS
		for( int j=0; j<pf.size() ; j++ ) {
			if(gsl_vector_get(ci,j)==i)
				gsl_vector_add(pf[j].second,va);
		}
	}

	gsl_matrix_free(C2);
	gsl_matrix_free(C1);
	gsl_vector_free(va);
	gsl_vector_free(vb);

	return pf;
}

particles
node_analysis::seeded_centroids() { 

	particles pf;

	gmat *M = gsl_matrix_calloc(DIM, parents_.first.length_M() );
	gmat *C = gsl_matrix_calloc(DIM, parents_.second.length_M() );
	gvec *v = gsl_vector_calloc(parents_.first.length_M()  );
	gvec *w = gsl_vector_calloc(parents_.second.length_M() );

	parents_.first.copyM(M);
	parents_.second.copyM(C);
	parents_.first.gsl_seeded_kmeans( M , v , C , w );

	richanalysis::coord_format cf;
	gsl_vector_set_all(w,1);
	pf 	= cf.mat2par (C, w);

	gsl_matrix_free(M);
	gsl_matrix_free(C);
	gsl_vector_free(v);
	gsl_vector_free(w);

	return pf;
}


void
cluster::print_neighbors(void) {
	if(bHaveNN_) {
		for( int i=0 ; i<nnList_.size() ; i++ ) {
			std::cout << "CLUSTER " << i << " HAS " 
				<< nnList_[i][0] << " NEIGHBORS:  "  
				<< std::endl;
			for( int j=1 ; j<nnList_[i].size() ;j++ ) {
				std::cout << " " << nnList_[i][j];
			}
			std::cout <<  std::endl;
		}
	} else {
		std::cout << "ERROR::NO LIST" << std::endl;
	}
}

void
cluster::calculate_nn_and_contacts(void)
{
	int verbose 	= 0 ;

	calc_distance_matrix(0);
	gmat *D		= gsl_matrix_calloc( length_M(), length_M() );
	gvec *v		= gsl_vector_calloc( length_M() );
//
	//HERE NOW
//
	copyv(v);
	copyA(D);

	gsl_matrix_free(D);
	gsl_vector_free(v);
}

void
cluster::calculate_neighbors(void)
{
	int verbose	= 0;
	int   NC 	= length_C();
	gmat *Ci 	= gsl_matrix_calloc(DIM,length_C() );
	gmat *Mi 	= gsl_matrix_calloc(DIM,length_M() );
	gvec *sigmas 	= gsl_vector_calloc( length_C()  );
	gvec *m2c 	= gsl_vector_calloc( length_M()  );

	gvec *rc	= gsl_vector_calloc(DIM);
	gvec *dr	= gsl_vector_calloc(DIM);
	gvec *rhat	= gsl_vector_calloc(DIM);
	gvec *rI	= gsl_vector_calloc(DIM);
	gvec *ri	= gsl_vector_calloc(DIM);
	gvec *rj	= gsl_vector_calloc(DIM);
	gvec *rw	= gsl_vector_calloc(DIM);

	copyws(sigmas); 
	copyC(Ci); copyM(Mi); copyv(m2c);
	nnList_.clear(); 
	bHaveNN_=false;

	if(verbose){
		std::cout << "HAVE SIGMAS::";
		for( int i=0 ; i<sigmas->size;i++)
			std::cout << gsl_vector_get(sigmas,i) << " ";
		std::cout << std::endl;
		std::cout << "HAVE LABELS::";
		for( int i=0 ; i<m2c->size;i++)
			std::cout << gsl_vector_get(m2c,i) << " ";
		std::cout << std::endl;
	}

//	i <-> j <==> j <-> i
// INDEX -> NR_OF_NEIGHBORS+LIST

	for(int i=0;i<Ci->size2;i++) {
		std::vector< int > vid;
		vid.push_back(0);
		nnList_.push_back(vid);
		for(int j=0;j<Ci->size2;j++) {
			if( i==j )
				continue;
			gsl_matrix_get_col	( ri, Ci, i 	);
			gsl_matrix_get_col	( rj, Ci, j 	);
			gsl_vector_memcpy	( rc, ri 	);
			gsl_vector_add		( rc, rj 	);
			gsl_vector_scale	( rc, 0.5 	); 	// MIDPOINT
			gsl_vector_memcpy	( rhat, ri 	);	
			gsl_vector_sub		( rhat, rj 	);
			double nrd 	= 	gsl_blas_dnrm2( rhat );
			gsl_vector_scale	( rhat, 1.0/nrd );	// DIRECTOR

			double tol = 0.0;
			double sig = gsl_vector_get(sigmas,i); 

			for(int I=0;I<Mi->size2;I++){
				if( gsl_vector_get(m2c,I) != i )
					continue;
				gsl_matrix_get_col	( rI, Mi, I );
				gsl_vector_sub	  	( rI, rc );
				gsl_blas_ddot( rI, rhat, &tol ); 
				if(verbose)
					std::cout << "TOL::" << tol << std::endl;
				if( tol<sig*2.0 && tol >= 0.0 ) {			// HIT!		
					nnList_[i][0]++;
					nnList_[i].push_back(j);
					break;				// STOP LOOKING 
				}	
			}	
		}
	}

//	CHECK FOR INCOMPLETE NEIGHBOR PAIRS HERE (NOT IMPLEMENTED)

	if(verbose) { 
		for(int i=0;i<nnList_.size();i++) {
			std::cout << "CLUSTER " << i << " HAS SIGMA " << gsl_vector_get(sigmas,i) 
				<<" AND "<< nnList_[i][0] << " NUMBER OF NEIGHBORS:\t "  << std::endl;
			for(int j=1;j<nnList_[i].size();j++){
				std::cout << " " << nnList_[i][j];
			}
			std::cout <<  std::endl;
		}
	}
	bHaveNN_=true;

	gsl_matrix_free(	Ci	);
	gsl_matrix_free(	Mi	);
	gsl_vector_free(	rc	);
	gsl_vector_free(	dr	);
	gsl_vector_free(	rI	);
	gsl_vector_free(	ri	);
	gsl_vector_free(	rj	);
	gsl_vector_free(	rhat	);
	gsl_vector_free(	rw	);
	gsl_vector_free(	sigmas	);
	gsl_vector_free(	m2c	);

}

particles
node_analysis::get_centroids( int sw ) 
{
	particles pfc;
	richanalysis::coord_format cf;
	switch(sw){
		case 1:	{
				gvec *w = gsl_vector_calloc(	 parents_.first .length_C() );
				gsl_vector_set_all(w,1);
				gmat *C = gsl_matrix_calloc(DIM, parents_.first .length_C() );
				parents_.first .copyC(C);		
				pfc 	= cf.mat2par (C, w);
				gsl_matrix_free(C);
				gsl_vector_free(w);
			}
			break;
		case 2: {
				gvec *w = gsl_vector_calloc(	 parents_.second.length_C() );
				gsl_vector_set_all(w,1);
				gmat *C = gsl_matrix_calloc(DIM, parents_.second.length_C() );
				parents_.second.copyC(C);
				pfc 	= cf.mat2par (C, w);
				gsl_matrix_free(C);
				gsl_vector_free(w);
			}
			break;
		default:
			break;
	}

	for(int i=0;i<pfc.size();i++)
		if(sw==1)
			pfc[i].first="Ar";
		else
			pfc[i].first="U";
	return pfc;
}


int
cluster::find_centroids( ) {
	if( bSet_ ) {
		int M=vc_->size; int N=wc_->size;
		gsl_kmeans( M_ , vc_ , C_ , wc_ );
		NperC_.clear();
		gsl_vector_set_all( ws_, 0.0 );
		for(int i=0;i<N;i++) {
			int numi=0;
			for( int j=0 ; j<vc_->size ; j++ ) {
				numi+=(gsl_vector_get(vc_,j)==i)?1:0;
			}
			NperC_.push_back(numi);
		}
		double X[N], X2[N], Z[N], S[N];
//		std::cout << "INFO CENTROID SIGMA B : "<< N << " & " << M << std::endl;
		for( int j=0; j<M ; j++ ) 
		{
			gvec *r_dump 	= gsl_vector_calloc( DIM );
			gvec *r_cent 	= gsl_vector_calloc( DIM );
			gvec *r_temp 	= gsl_vector_calloc( DIM );
			int ic 		= gsl_vector_get( vc_, j );
			X [ic]  	= 0.0; 
			X2[ic]		= 0.0;

			gsl_matrix_get_col (r_dump, M_,  j );
			gsl_matrix_get_col (r_cent, C_, ic );
			gsl_vector_sub     (r_dump, r_cent );
			ftyp dlen = gsl_blas_dnrm2( r_dump );
			X[ic] += dlen; X2[ic] += dlen*dlen; 
			Z[ic] += 1.0 ;
		}
		for( int i=0 ; i<N ; i++ ) {
			if( Z[i] > 1.0 )	{
				S[i] = sqrt( X2[i]/(Z[i]-1.0) );
			} else {
				S[i] = 0.1;
				//std::cout << "INFO::UNHEALTHY::CLUSTER" << std::endl;
			}
			gsl_vector_set(ws_,i,S[i]);
		}
		//std::cout << "INFO CENTROID SIGMA E" << std::endl;
	}
	return NperC_.size();
}

std::vector<int>	
cluster::get_labels( void ){
	std::vector<int> ndx;
	if( bSet_ )
		for(int i=0;i<vc_->size;i++)
			ndx.push_back(gsl_vector_get(vc_,i));
	return ndx;
}

std::vector<int>	
cluster::get_clabels( void ){
	std::vector<int> ndx;
	if( bSet_ )
		for(int i=0;i<wc_->size;i++)
			ndx.push_back(gsl_vector_get(wc_,i));
	return ndx;
}

particles
cluster::get_model( void ){
	richanalysis::coord_format cfp;
	if( bSet_ ) {
		particles red_parts = cfp.mat2par( M_ , vc_ );
		return red_parts;
	} else {
		particles red_parts;
		return red_parts;
	}
}

particles
cluster::get_centroids( void ) {
	richanalysis::coord_format cfp;
	if( bSet_ ) {
		particles red_parts = cfp.mat2par( C_ , wc_ );
		return red_parts;
	} else {
		particles red_parts;
		return red_parts;
	}
}


bool logic_desc(std::pair<int, double> id1 ,std::pair<int, double> id2) {
	return id1.second>id2.second;
}

bool logic_asce(std::pair<int, double> id1 ,std::pair<int, double> id2) {
	return id1.second<id2.second;
}


void
cluster::order_centroids(void) {
	// SANITY CHECK 
	if( isSet() ) {
		int N 		= length_C();
		particles cents	= get_centroids( );
		gsl_vector *r0	= gsl_vector_calloc(DIM);
		int n		= ((float)N);
		for(int i=0;i<N;i++) {
			gsl_vector_add(r0, cents[i].second);
		}
		if(n==0){
			std::cout << "CANNOT PROCEED WITH CENTROID ORDERING" << std::endl;
			exit(-1);
		}
		gsl_vector_scale( r0, 1.0/n );
		std::vector< std::pair<int, double> > vpid;
		for(int i=0;i<N;i++){
			std::pair<int,double> pid;
			gsl_vector_sub(cents[i].second,r0);
			pid.first  = i; 
			pid.second = gsl_blas_dnrm2(cents[i].second);
			vpid.push_back(pid);
		}
		std::sort( vpid.begin(), vpid.end(), logic_asce );
		if(o_idx_.size()>0)
			o_idx_.clear();
		for(int i=0;i<N;i++) {
			std::pair< int, int > spii;
			spii.first	= i;
			spii.second	= vpid[i].first;
			o_idx_.push_back(spii);			
		}
	}
}

bool
node_analysis::assign_node( node n ) {
	std::cout << "INFO:: ASSIGNING NODE" << std::endl;
	n.first.find_centroids();
	n.second.find_centroids();
//	std::cout << "INFOINFO::HAVE::length_C:: " << n.second.length_C() << std::endl;
	cluster c1 = n.first;
	cluster c2 = n.second;
	std::vector<int> cindx;

	int N, M, L, K;
	int rDIM=c1.get_cDIM();
//	std::cout << "INFOINFO::HAVE::CDIM:: " << rDIM << std::endl;
	if( c1.get_cDIM() != c2.get_cDIM() ) {
		std::cout << "ERROR WITH DIMENSION OF CLUSTER PAIRS" << std::endl;
		bNode_= false; bLayer_= false;
		return false;
	}

	if( c1.isSet() && c2.isSet() ) {
		parents_.first  = c1;
		parents_.second = c2;
		L = c1.length_M();
		N = c2.length_M();
		M = c2.length_C();
		K = c1.length_C();
//		std::cout << "INFOINFO::HAVE::M:: " << M << std::endl;
		if( K == M && M == rDIM && rDIM>1 ) {

			std::vector<int> vi;
			vi = find_centroid_relation(); 

			gmat *M1 = gsl_matrix_calloc(DIM,L);
			gvec *v1 = gsl_vector_calloc(L);
			gmat *M2 = gsl_matrix_calloc(DIM,N);
			gvec *v2 = gsl_vector_calloc(N);
			c2.copyM(M2); c2.copyv(v2);
			c1.copyM(M1); c1.copyv(v1);
			ids IDX,IDY;
			IDX = c1.getIDs();
			IDY = c2.getIDs();
			for(int i=0;i<v2->size;i++){
				double val = gsl_vector_get(v2,i);
				int ival = (int)round(val);
				int nval = vi[ival];
				gsl_vector_set(v2,i,nval);
				cindx.push_back(nval);
			}
			parents_.second.writev(v2);
			cidx_.swap(cindx);
			bLayer_ = false;
			for( int ipart=0 ; ipart<M ; ipart++ ) { // M CLUSTERS
				richanalysis::coord_format cf;
				particles	px, py;
				ids		idx, idy;
				px	= cf.truncmat( M1, v1, ipart);
				idx	= cf.truncIDs(IDX, v1, ipart);
				py	= cf.truncmat( M2, v2, ipart);
				idy	= cf.truncIDs(IDY, v2, ipart);
				int D	= px.size();
				int B	= py.size();
				int rD	= (B<DIM)?((B>0)?(B):(-1)):((D<DIM)?((D>0)?(D):(-1)):(DIM));
				node child;
				if( rD > 0 ) {
					richanalysis::cluster clpx,clpy;
					clpx.set_cDIM(rD); clpy.set_cDIM(rD);	
					clpx.set_matrix( px );
					clpy.set_matrix( py );
					child.first	= clpx;
					child.first.setIDs(idx);
					child.second	= clpy;
					child.second.setIDs(idy);
					bLayer_ = true;	
				}
				children_.push_back(child);
			}
		} else {
			bLayer_= false;
		}
		bNode_ = true;
	}else{
		bNode_ = false;
	}

	return bNode_;
}

particles
node_analysis::nn_restraint_fit( int verbose ) {

	particles pf;
	cluster c1 	= parents_.first;
	cluster c2 	= parents_.second;

	gmat *P 	= gsl_matrix_calloc( DIM, c2.length_M());
	gmat *Q 	= gsl_matrix_calloc( DIM, c1.length_M());

	gmat *C 	= gsl_matrix_calloc( DIM, c2.length_C());
	gmat *K 	= gsl_matrix_calloc( DIM, c1.length_C());

	gmat *U 	= gsl_matrix_calloc( DIM, DIM );
	gmat *iU 	= gsl_matrix_calloc( DIM, DIM );
	gvec *it 	= gsl_vector_calloc( DIM );
	gvec *t 	= gsl_vector_calloc( DIM );
	gvec *vl 	= gsl_vector_calloc( c2.length_M() );

	richanalysis::coord_format cf;

	gsl_vector_set_all(vl,1);
	c2.copyM(P); c2.copyC(C);
	c1.copyM(Q); c1.copyC(K);

//	HERE WE HAVE NEEDED DATA NOW FIND RESTRAINT POINTS FOR THE MODEL FIT		
//	GOAL:	CHECK MUTUAL NEIGHBOR CORRESPONDANCE					
//		CALCULATE A VECTOR POSITION (INTEGER OF DENSITY CLOSE TO BORDER) 		
//		FIND THE (INTEGER) MODEL PART THAT BEST CORRESPONDS TO SUCH A MIDPOINT		NEW ROUTINE
//		FIND THE BEST FIT THAT MINIMIZES THIS DISTANCE				OVERLOADED SHAPE FIT

	if( !c1.haveNN() || !c2.haveNN() ) {
		c1.calculate_neighbors();
		c2.calculate_neighbors();
	}

	//std::cout << "INFO::MODEL::MAX::NN::DIST:: " << c2.max_dist() << std::endl; 

	std::vector< std::vector<  int  > > nc1 = c1.get_neighbors();
	std::vector< std::vector<  int  > > nc2 = c2.get_neighbors();
	std::vector<int> ndx21 = find_centroid_distance_relation();
	std::vector<int> ndx12;
	for(int i=0;i<ndx21.size();i++) // init
		ndx12.push_back(0.0);
	for(int i=0;i<ndx21.size();i++) // set
		ndx12[ndx21[i]]=i;
	if( verbose ) { 
		std::cout << "INFO::" << std::endl;
		for(int i=0;i<ndx21.size();i++){
			std::cout << " " << i << " | " << ndx21[i] << " " << ndx12[i] << std::endl;
		}
	}
//	DO TWO CLUSTERS HAVE A COMMON BORDER?
	if( c1.length_C() != c2.length_C() )
		std::cout << "ERROR::ILL FORMATED NODE" << std::endl;
	if( c1.length_C() != nc2.size() )
		std::cout << "ERROR::BAD NEIGHBOR LIST" << std::endl;

	if( verbose ) {
//	FOR PRINTING NEIGHBOR LIST...
		c1.print_neighbors(); // LOCAL INDEXING
		c2.print_neighbors(); // LOCAL INDEXING
	}

	std::vector< std::pair< int, int > > common;
	for( int i=0 ; i<c1.length_C() ; i++ ) {
		if( nc1[i][0]>0 ) {
			for( int j=1 ; j<=nc2[ndx21[i]][0] ; j++ ) {
				for(int k=0;k<nc1[i][0];k++) {
					if(nc2[ndx21[i]][j]==nc1[i][1+k]) {
						std::pair< int, int > ip;
						ip.first=i;ip.second=nc2[ndx21[i]][j];
						common.push_back(ip);
					}
				}
			}
		}
	}

//	CREATE THE COMMON BORDER LIST
	std::vector< std::vector< int > > common_nn;
	for( int i=0 ; i<c1.length_C() ; i++ ) {
		std::vector< int > nns;
		nns.push_back(0);
		common_nn.push_back(nns);
	}
	for( int i=0 ; i<common.size() ; i++ ) {
		bool bSkip=false;
		for( int j=1 ; j<=common_nn[ common[i].first ][0] ; j++ ) {
		//	MIGHT HAVE SYMMETRY ORDER ISSUE HERE, NOPE
			if( common_nn[ common[i].first ][j] == common[i].second ){
				bSkip=true; break;
			}
		}
		if( !bSkip ) {
			common_nn[common[i].first ][0]++;
			common_nn[common[i].first ].push_back(common[i].second);
			common_nn[common[i].second][0]++;
			common_nn[common[i].second].push_back(common[i].first);
		}

	}

	if( verbose ) { 
		for( int i=0 ; i<common_nn.size() ; i++ ) {
			std::cout << "INFO:: \t"<< i << "\t < " << common_nn[i][0] << " > " ;
			for( int j=1 ; j<=common_nn[i][0] ; j++ ) {
				std::cout << " " << common_nn[i][j];
			}
			std::cout << std::endl;
		}
	}

//	WHICH INDEX IN THE MODEL IS CLOSEST TO THE SPECIFIC COMMON BORDER?
//	get midpoint from MBIG calc mindist in set MSMALL to midpoint
//	NOTE common_nn[i][0] IS (N)UMBER OF RESTRAINTS FOR i
//		common_nn[i][j] IS THE MODEL NEIGHBOR CLUSTER INDEX
//		N TIMES: CALCULATE THE MIDPOINT AND LOOP THIS MODEL 
//			CLUSTER RETURN CLOSEST MODEL PART
// 
	gvec *r0 	= gsl_vector_calloc( DIM );
	gvec *r1 	= gsl_vector_calloc( DIM );
	gvec *r2 	= gsl_vector_calloc( DIM );
	gvec *rm 	= gsl_vector_calloc( DIM );
	std::vector< int > m_labels = c2.get_labels();
//	HERE WE MIGHT HAVE THE WRONG CORRESPONDANCE 1->2 NOT 2->1, TRY BOTH
	std::cout << "INFO::" << C->size2 << " AND " << K->size2 << std::endl;
	std::cout << "INFO::" << P->size2 << " AND " << Q->size2 << std::endl;
	std::cout << "INFO::" << m_labels.size() << std::endl;
	for( int ic=0 ; ic<common_nn.size() ; ic++ ) {
		gsl_matrix_get_col(r0,C,ic);
		for( int ir=1 ; ir<=common_nn[ic][0] ; ir++ ) {
			int ri = common_nn[ic][ir];
			gsl_matrix_get_col	( r1, C, ri );
			gsl_vector_memcpy	( r2, r0  );
			gsl_vector_add		( r2, r1  );
			gsl_vector_scale	( r2, 0.5 );
			gsl_vector_memcpy	( rm, r2  );
			double  mlen = 1E10;
			int	imin = 0;
			for(int im=0;im<P->size2;im++) {
				if(ndx21[m_labels[im]] != ic )
					continue;
				gsl_matrix_get_col	( r1, P, im );
				gsl_vector_memcpy	( r2, rm  );
				gsl_vector_sub		( r2, r1  );
				double len	= gsl_blas_dnrm2( r2 );
				if( len < mlen ) {
					imin = im; 
					mlen = len;
				}
			}
			std::cout << ic << " " << imin << std::endl;
			std::cout << "::INFO::FOUND::"<< std::endl; 
			//std::cout << "INFO::FOUND::IDX " << imin << " WITH LENGTH " << mlen << std::endl; 
			gsl_matrix_get_col	( r1, P, imin );
			richanalysis::tensorIO tIO;
			tIO.output_vector(r1);
		}
	}
	//output all of these restraint positions

	gsl_vector_free(r0);
	gsl_vector_free(r1);
	gsl_vector_free(r2);
	gsl_vector_free(rm);
//	BELOW IS A SCRATCH FIELD
//	pf  = cf.mat2par (P, vl);
//	if(verbose)
//		std::cout << "<< INFO >>" << c1.length_M() << " " << c2.length_M() << std::endl;
/*
//	DO ALL EIGHT
			ftyp	shape_fit	(	gmat *P , gmat *Q ,	// IN
					  		gvec *w1, gvec *w2, 	// IN
					  		gmat *U , gvec *t ,	// OUT 
					  		int II );
*/
/*
	ftyp rmsd	= shape_fit( P, Q, U, t ); 
	output_matrix(U);
	invert_fit	( U, t, iU, it );
	output_matrix(iU);
	apply_fit( pf, iU, it ); 

	bool overwrite_model=true;
	for( int i=0 ; i<pf.size() ; i++ ){
		pf[i].first = "C";
	}

	if( overwrite_model ) {
		parents_.second.set_matrix(pf);
		parents_.second.find_centroids();
	}
*/
	return pf;
}

particles
node_analysis::regular_fit( int choice ) {

	particles pf;
	cluster c1 	= parents_.first;
	cluster c2 	= parents_.second;

	gmat *P 	= gsl_matrix_calloc( DIM, c2.length_M());
	gmat *Q 	= gsl_matrix_calloc( DIM, c1.length_M());
	gmat *U 	= gsl_matrix_calloc( DIM, DIM );
	gmat *iU 	= gsl_matrix_calloc( DIM, DIM );
	gvec *it 	= gsl_vector_calloc( DIM );
	gvec *t 	= gsl_vector_calloc( DIM );
	gvec *vl 	= gsl_vector_calloc( c2.length_M() );

	richanalysis::coord_format cf;

	gsl_vector_set_all(vl,1);
	c2.copyM(P);
	c1.copyM(Q);

	pf 		= cf.mat2par (P, vl);

	std::cout << "INFO > " << c1.length_M() << " " << c2.length_M() << std::endl;

	ftyp rmsd	= shape_fit( P, Q, U, t, choice ); 
	output_matrix(U);
	invert_fit	( U, t, iU, it );
	output_matrix(iU);
	apply_fit( pf, iU, it ); 

	bool overwrite_model=true;
	for( int i=0 ; i<pf.size() ; i++ ){
		pf[i].first = "C";
	}

	if( overwrite_model ) {
		parents_.second.set_matrix(pf);
		parents_.second.find_centroids();
	}

	return pf;
}


particles
node_analysis::regular_fit(void) {

	particles pf;
	cluster c1 	= parents_.first;
	cluster c2 	= parents_.second;

	gmat *P 	= gsl_matrix_calloc( DIM, c2.length_M());
	gmat *Q 	= gsl_matrix_calloc( DIM, c1.length_M());
	gmat *U 	= gsl_matrix_calloc( DIM, DIM );
	gmat *iU 	= gsl_matrix_calloc( DIM, DIM );
	gvec *it 	= gsl_vector_calloc( DIM );
	gvec *t 	= gsl_vector_calloc( DIM );
	gvec *vl 	= gsl_vector_calloc( c2.length_M() );

	richanalysis::coord_format cf;

	gsl_vector_set_all(vl,1);
	c2.copyM(P);
	c1.copyM(Q);

	pf 		= cf.mat2par (P, vl);

	std::cout << "INFO_INFO_INFO>" << c1.length_M() << " " << c2.length_M() << std::endl;

	ftyp rmsd	= shape_fit( P, Q, U, t, -1); // THIS NOW WORKS:: NEEDED N4 NORM
	output_matrix(U);
	invert_fit	( U, t, iU, it );
	output_matrix(iU);
	apply_fit( pf, iU, it ); 

	bool overwrite_model=true;
	for( int i=0 ; i<pf.size() ; i++ ){
		pf[i].first = "C";
	}

	if( overwrite_model ) {
		parents_.second.set_matrix(pf);
		parents_.second.find_centroids();
	}

	return pf;
}

std::vector<int> 
node_analysis::find_centroid_relation( void ) {
	ftyp min_rmsd	= 1.0E10;

	cluster c1 = parents_.first;
	cluster c2 = parents_.second;

	int N  = c1.length_C();

	if( N != c2.length_C() )
		std::cout << "ERROR IN IDX REL" << std::endl;

	gmat *U		= gsl_matrix_calloc( DIM, DIM );
	gvec *t		= gsl_vector_calloc( DIM );

	int J		= 0;

	gmat *C0T	= gsl_matrix_calloc( DIM, N );
	gmat *C0N	= gsl_matrix_calloc( DIM, N );
	gmat *CNT	= gsl_matrix_calloc( DIM, N );
	gmat *CEN	= gsl_matrix_calloc( DIM, N );
	gvec *gv 	= gsl_vector_calloc( DIM );

	c1.copyC(C0N);
	c2.copyC(C0T);

	gmat *Ut = gsl_matrix_calloc( DIM, DIM );
	gvec *tt = gsl_vector_calloc( DIM );

	std::vector<int> iv;
	for(int i=0;i<N;i++)
		iv.push_back(i);

	std::vector<std::vector<int> > imv = all_permutations(iv);

	for(int j=0;j<imv.size();j++) {
		gsl_matrix_memcpy (CEN, C0N);
		for(int i=0;i<iv.size();i++) {
			gsl_matrix_get_col (gv, C0T, i);
			gsl_matrix_set_col (CNT, imv[j][i], gv );
		}
		ftyp rmsd = kabsch_fit(CNT,CEN,U,t);
		if(rmsd<min_rmsd) {
			J=j; 
			min_rmsd=rmsd;	
			gsl_matrix_memcpy (Ut, U);
			gsl_vector_memcpy (tt, t);
		}
	}

	if( idx_.size()!=0 ) {
		std::cout << "INFO::ERROR WITH IDX_" << std::endl;
		idx_.clear();
		iidx_.clear();
	}

	for(int i=0;i<iv.size();i++)
		iidx_.push_back(0);

	for(int i=0;i<iv.size();i++) {
		idx_.push_back(imv[J][i]);
		iv[i]=imv[J][i];
		iidx_[imv[J][i]]=i;
	}
	
	gsl_matrix_free(CNT);
	gsl_matrix_free(C0N);
	gsl_matrix_free(C0T);
	gsl_matrix_free(U);
	gsl_vector_free(t);
	gsl_vector_free(gv);

	return iv;
}

void
layer_analysis::output_layer ( std::string filename ) {
	if( bSet_ ) {
		int N = own_.size();
		bool bFirst=true;
		particles ps;

		ids glob_ids;
		for(int i=0;i<N;i++){
			ids	ids1, ids2;
			ids1		= own_[i].first.getIDs();
			ids2		= own_[i].second.getIDs();
			bFirst 		= own_[i].first.length_M()>=own_[i].second.length_M();
			cluster c_d 	= bFirst?own_[i].first:own_[i].second;
			int	N_d 	= bFirst?own_[i].first.length_M():own_[i].second.length_M();
			gsl_matrix *m	= gsl_matrix_alloc(DIM,N_d);
			gsl_vector *v	= gsl_vector_calloc(DIM);
			gsl_vector *w	= gsl_vector_calloc(DIM);
			c_d.copyM(m);
			for(int j=0;j<N_d;j++){
				gsl_matrix_get_col ( v, m, j );
				gsl_vector_add( w, v );
			}
			gsl_vector_scale(w,1.0/((float)(N_d)));
			particle p;
			id identity	=	bFirst?(ids2[0]):(ids1[0]);
			p.first  	=	identity.second;
			glob_ids.push_back(	identity	);
			p.second 	= gsl_vector_alloc(DIM);
			gsl_vector_memcpy(p.second,w);
			ps.push_back(p);
			gsl_matrix_free(m);
			gsl_vector_free(v);
			gsl_vector_free(w);
		}
		output_pdb( filename , ps );
//		for(int i=0;i<glob_ids.size();i++)
//			std::cout << " SCRAMBLED " << glob_ids[i].first << " " << glob_ids[i].second << std::endl;
	}
}

int 
clustering::connectivity(gsl_matrix *B, ftyp val) {
	int i,j,k,q,N,C=0,min;

	int nr_sq=B->size1;
	if(B->size1!=B->size2)
		return(-1);

	std::vector<int>	res;
	std::vector<int>	nvisi;

	std::vector<int>	s;
	std::vector<int>	NN;
	std::vector<int>	ndx;

	N=nr_sq;
	res.push_back(0);
	for(i=0; i<N; i++ ){
   		nvisi.push_back(i+1);
		res.push_back(0);res.push_back(0);
		ndx.push_back(i);
	}

	while(!ndx.empty()){
		i=ndx.back(); ndx.pop_back(); 
		NN.clear();

		if(nvisi[i]>0){
			C--;
			for(j=0;j<N;j++){
				if( gsl_matrix_get(B,j,i)<val==1 ) {
					NN.push_back(j);
				}
			}
			while(!NN.empty()){
				k=NN.back(); NN.pop_back();
				nvisi[k]=C;
				for(j=0;j<N;j++){
	   				if( gsl_matrix_get(B,j,k)<val==1 ){
	     					for(q=0;q<N;q++){
							if(nvisi[q]==j+1){
								NN.push_back(q);
							}
						}
					}
				}
			}
		}
	}
   
	std::cout <<"INFO " << C << "\ncluster data:\n";
	std::vector<int> Nc; // NUMBER OF POINTS IN EACH CLUSTER
	for(i=0;i<-1*C;i++)
		Nc.push_back(0);

	for(q=0;q<N;q++){
		res[q*2+1]=q;
		res[q*2]=nvisi[q]-C;
		Nc[res[q*2]]++;
		std::cout << " " << res[q*2] << " " << res[2*q+1] << std::endl ;
   	}
	// HOW MANY IN EACH?
	for(i=0;i<-1*C;i++)
		std::cout << "CLUSTER " << i << " HAS " << Nc[i] << " ELEMENTS" << std::endl ;

	return(0);   
}

int 
clustering::gsl_seeded_kmeans(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw ){
	int NN=dat->size2, MM=dat->size1, KK=cent->size2;

	gvec *labels = gsl_vector_calloc(NN);
	gvec *counts = gsl_vector_calloc(KK);
	gmat *tmp_ce = gsl_matrix_calloc(MM,KK);

	if( ((int)cent->size1)!=MM || !gsl_vector_isnonneg(w) || ((int)w->size)!=NN || ((int)nw->size)!=KK )	{
		std::cout << "ERROR::BAD::KMEANS" << std::endl;
		return -1;
	}
	gsl_vector_set_zero(nw);
	gsl_vector_set_zero(w);

	int h, i, j;
	ftyp old_error, error, TOL=1E-10; 

	do {
		old_error = error, error = 0; 
		for (i = 0; i < KK; i++ ) {
			gsl_vector_set(counts,i,0);
			for (j = XX; j <= ZZ; j++){
				gsl_matrix_set(tmp_ce,j,i,0.0);
			}
 		}
		for (h = 0; h < NN; h++) {
			ftyp min_distance = 1E30;
			for (i = 0; i < KK; i++) {
				ftyp distance = 0;
				for ( j = XX; j<=ZZ ; j++ ) {
					distance += square( gsl_matrix_get( dat,j,h ) - gsl_matrix_get( cent,j,i ) );
				}
				if (distance < min_distance) {
	 				gsl_vector_set(labels,h,i); min_distance = distance; 
				} 
			} 
			for ( j=XX ; j<=ZZ ; j++ ){
 				gsl_matrix_set( tmp_ce, j, gsl_vector_get(labels,h),
				 gsl_matrix_get( dat, j, h ) + gsl_matrix_get(tmp_ce, j, gsl_vector_get(labels,h)) );
			}
			gsl_vector_set(counts,gsl_vector_get(labels,h),1.0+gsl_vector_get(counts,gsl_vector_get(labels,h)));
	 		error += min_distance; 
		}
	 	for (i = 0; i < KK; i++) {
	 		for ( j=XX ; j<=ZZ ; j++ ) {
				gsl_matrix_set(cent,j,i,
				gsl_vector_get(counts,i)?(gsl_matrix_get(tmp_ce,j,i)/gsl_vector_get(counts,i)):(gsl_matrix_get(tmp_ce,j,i)));
	 		}
	 	}
	} while ( fabs(error - old_error) > TOL );	// WHILE THEY ARE MOVING

	ftyp wi=0.0, nwh=0.0;
	for( i=0 ; i<NN ; i++) {
		h 	= gsl_vector_get(labels,i);
		wi	= gsl_vector_get(w,i);
		nwh	= gsl_vector_get(nw,h);
		gsl_vector_set(w,i,h);		// MIGHT NOT WANT TO OVERWRITE THIS
		gsl_vector_set(nw,h,nwh+wi); 		// NOT NORMALIZED HERE
	}

	gsl_vector_free(labels); 
	gsl_vector_free(counts);
	gsl_matrix_free(tmp_ce);

	return 0;
}


int 
clustering::gsl_kmeans0(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw ){
	int NN=dat->size2, MM=dat->size1, KK=cent->size2;

	gvec *labels = gsl_vector_calloc(NN);
	gvec *counts = gsl_vector_calloc(KK);
	gmat *tmp_ce = gsl_matrix_calloc(MM,KK);

	if( ((int)cent->size1)!=MM || !gsl_vector_isnonneg(w) || ((int)w->size)!=NN || ((int)nw->size)!=KK )	{
		std::cout << "ERROR::BAD::KMEANS" << std::endl;
		return -1;
	}
	gsl_vector_set_zero(nw);
	gsl_vector_set_zero(w);

	int h, i, j;
	ftyp old_error, error = 1E30, TOL=1E-10; 
/*
 *	std::vector<int> myvector;
 *	for ( i=0; i<NN; ++i ) myvector.push_back(i); 
 *	int steps = 1+NN%2;
 *	for ( i=0; i < steps ; i++ )
 *		random_shuffle ( myvector.begin(), myvector.end() );
 *	i=0;
 *	for ( h=0 ; h<KK ; h++ ){
 *		for ( j=XX ; j<=ZZ ; j++ ) { 
 *			gsl_matrix_set( cent, j, h , gsl_matrix_get( dat, j, myvector[h] ) );
 *		} 
 *	}
 */

	do {
		old_error = error, error = 0; 
		for (i = 0; i < KK; i++ ) {
			gsl_vector_set( counts, i, 0 );
			for (j = XX; j <= ZZ; j++){
				gsl_matrix_set( tmp_ce, j, i, 0.0 );
			}
 		}
		for (h = 0; h < NN; h++) {
			ftyp min_distance = 1E30;
			for (i = 0; i < KK; i++) {
				ftyp distance = 0;
				for ( j = XX; j<=ZZ ; j++ ) {
					distance += square( gsl_matrix_get( dat,j,h ) - gsl_matrix_get( cent,j,i ) );
				}
				if (distance < min_distance) {
	 				gsl_vector_set(labels,h,i); min_distance = distance; 
				} 
			} 
			for ( j=XX ; j<=ZZ ; j++ ){
 				gsl_matrix_set( tmp_ce, j, gsl_vector_get(labels,h),
				 gsl_matrix_get( dat, j, h ) + gsl_matrix_get(tmp_ce, j, gsl_vector_get(labels,h)) );
			}
			gsl_vector_set(counts,gsl_vector_get(labels,h),1.0+gsl_vector_get(counts,gsl_vector_get(labels,h)));
	 		error += min_distance; 
		}
	 	for (i = 0; i < KK; i++) {
	 		for ( j=XX ; j<=ZZ ; j++ ) {
				gsl_matrix_set(cent,j,i,
				gsl_vector_get(counts,i)?(gsl_matrix_get(tmp_ce,j,i)/gsl_vector_get(counts,i)):(gsl_matrix_get(tmp_ce,j,i)));
	 		}
	 	}
	} while ( fabs(error - old_error) > TOL );	// WHILE THEY ARE MOVING

	ftyp wi=0.0, nwh=0.0;
	for( i=0 ; i<NN ; i++) {
		h 	= gsl_vector_get(labels,i);
		wi	= gsl_vector_get(w,i);
		nwh	= gsl_vector_get(nw,h);
		gsl_vector_set(w,i,h);		// MIGHT NOT WANT TO OVERWRITE THIS
		gsl_vector_set(nw,h,nwh+wi); 		// NOT NORMALIZED HERE
	}

	gsl_vector_free(labels); 
	gsl_vector_free(counts);
	gsl_matrix_free(tmp_ce);

	return 0;
}


int 
clustering::gsl_kmeans(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw ){
	int NN=dat->size2, MM=dat->size1, KK=cent->size2;

	gvec *labels = gsl_vector_calloc(NN);
	gvec *counts = gsl_vector_calloc(KK);
	gmat *tmp_ce = gsl_matrix_calloc(MM,KK);

	if( ((int)cent->size1)!=MM || !gsl_vector_isnonneg(w) || ((int)w->size)!=NN || ((int)nw->size)!=KK )	{
		std::cout << "ERROR::BAD::KMEANS" << std::endl;
		return -1;
	}
	gsl_vector_set_zero(nw);
	gsl_vector_set_zero(w);

	int h, i, j;
	ftyp old_error, error = 1E30, TOL=1E-10; 

	std::vector<int> myvector;
	for ( i=0; i<NN; ++i ) myvector.push_back(i); 
	int steps = 1+NN%2;
	for ( i=0; i < steps ; i++ )
		random_shuffle ( myvector.begin(), myvector.end() );
	
	i=0;
	for ( h=0 ; h<KK ; h++ ){
		for ( j=XX ; j<=ZZ ; j++ ) { 
			gsl_matrix_set( cent, j, h , gsl_matrix_get( dat, j, myvector[h] ) );
		} 
	}

	do {
		old_error = error, error = 0; 
		for (i = 0; i < KK; i++ ) {
			gsl_vector_set( counts, i, 0 );
			for (j = XX; j <= ZZ; j++){
				gsl_matrix_set( tmp_ce, j, i, 0.0 );
			}
 		}
		for (h = 0; h < NN; h++) {
			ftyp min_distance = 1E30;
			for (i = 0; i < KK; i++) {
				ftyp distance = 0;
				for ( j = XX; j<=ZZ ; j++ ) {
					distance += square( gsl_matrix_get( dat,j,h ) - gsl_matrix_get( cent,j,i ) );
				}
				if (distance < min_distance) {
	 				gsl_vector_set(labels,h,i); min_distance = distance; 
				} 
			} 
			for ( j=XX ; j<=ZZ ; j++ ){
 				gsl_matrix_set( tmp_ce, j, gsl_vector_get(labels,h),
				 gsl_matrix_get( dat, j, h ) + gsl_matrix_get(tmp_ce, j, gsl_vector_get(labels,h)) );
			}
			gsl_vector_set(counts,gsl_vector_get(labels,h),1.0+gsl_vector_get(counts,gsl_vector_get(labels,h)));
	 		error += min_distance; 
		}
	 	for (i = 0; i < KK; i++) {
	 		for ( j=XX ; j<=ZZ ; j++ ) {
				gsl_matrix_set(cent,j,i,
				gsl_vector_get(counts,i)?(gsl_matrix_get(tmp_ce,j,i)/gsl_vector_get(counts,i)):(gsl_matrix_get(tmp_ce,j,i)));
	 		}
	 	}
	} while ( fabs(error - old_error) > TOL );	// WHILE THEY ARE MOVING

	ftyp wi=0.0, nwh=0.0;
	for( i=0 ; i<NN ; i++) {
		h 	= gsl_vector_get(labels,i);
		wi	= gsl_vector_get(w,i);
		nwh	= gsl_vector_get(nw,h);
		gsl_vector_set(w,i,h);		// MIGHT NOT WANT TO OVERWRITE THIS
		gsl_vector_set(nw,h,nwh+wi); 		// NOT NORMALIZED HERE
	}

	gsl_vector_free(labels); 
	gsl_vector_free(counts);
	gsl_matrix_free(tmp_ce);

	return 0;
}

/*
int 
clustering::gslkmeans(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw ) {
	int NN=dat->size2, MM=dat->size1, KK=cent->size2;
//
//	std::mt19937 generator;
//	std::uniform_int_distribution<int> distribution(1,KK-1);
//	int roll = distribution(generator);
//
	gvec *labels = gsl_vector_calloc(NN);
	gvec *counts = gsl_vector_calloc(KK);
	gmat *tmp_ce = gsl_matrix_calloc(MM,KK);

	gvec *rt0 = gsl_vector_calloc(DIM);
	gvec *rt1 = gsl_vector_calloc(DIM);

	if( ((int)cent->size1)!=MM || !gsl_vector_isnonneg(w) || ((int)w->size)!=NN || ((int)nw->size)!=KK )	{
		std::cout << "ERROR::BAD::KMEANS" << std::endl;
		return -1;
	}
	gsl_vector_set_zero(nw);
	gsl_vector_set_zero(w);

	int h, i, j;
	ftyp old_error, error = 1E30, TOL=1E-20; 

	std::vector<int> myvector;
	for ( i=0; i<NN; ++i ) myvector.push_back(i); 
//	int step = floor(NN/KK);
//	for ( i=0; i<KK; i++ ) myvector.push_back(distribution(generator));
	for ( i=0; i < KK ; i++ )
		random_shuffle ( myvector.begin(), myvector.end() );
	
	i=0;
	for ( h=0 ; h<KK ; h++ ) {
		gsl_matrix_get_col( rt0	, dat, myvector[h] );
		gsl_matrix_set_col(cent	, h, rt0 );
	}
	do {
		old_error = error, error = 0; 
		for (i = 0; i < KK; i++ ) {
			gsl_vector_set( counts, i, 0 );
			gsl_matrix_set_col(tmp_ce, i, rt0 );
 		}
		for (h = 0; h < NN; h++) {
			ftyp min_distance = 1E30;
			for (i = 0; i < KK; i++) {
				gsl_matrix_get_col( rt0	,  dat, h );
				gsl_matrix_get_col( rt1	, cent, i );
				gsl_vector_sub( rt0, rt1 );
				ftyp distance = gsl_blas_dnrm2( rt0 ); 
				if ( distance < min_distance ) {
	 				gsl_vector_set(labels,h,i); 
					min_distance = distance; 
				}
			} 
			int lh = gsl_vector_get(labels,h);
			gsl_matrix_get_col( rt0	, dat, h );
			gsl_matrix_get_col( rt1	, tmp_ce, lh );
			gsl_vector_add( rt0, rt1 );
			gsl_matrix_set_col( tmp_ce, lh, rt0 );
			gsl_vector_set(counts,lh,1.0+gsl_vector_get(counts,lh));
	 		error += min_distance; 
		}
	 	for (i = 0; i < KK; i++) {
			ftyp cnts = gsl_vector_get(counts,i);
			gsl_matrix_get_col( rt0, tmp_ce , i);
			if(cnts)
				gsl_vector_scale(rt0,1.0/cnts);
			gsl_matrix_set_col(cent,i,rt0);
	 	}
	} while ( fabs(error - old_error) > TOL );	// WHILE THEY ARE MOVING

	ftyp wi=0.0, nwh=0.0;
	for( i=0 ; i<NN ; i++) {
		h 	= gsl_vector_get(labels,i);
		wi	= gsl_vector_get(w,i);
		nwh	= gsl_vector_get(nw,h);
		gsl_vector_set(w,i,h);		// MIGHT NOT WANT TO OVERWRITE THIS
		gsl_vector_set(nw,h,nwh+wi); 	// NOT NORMALIZED HERE
	}

	gsl_vector_free(labels); 
	gsl_vector_free(counts);
	gsl_matrix_free(tmp_ce);

	return 0;
}
*/

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

using namespace richanalysis;

void
cluster::alloc_space( int D1, int D2 ) {
//D1	COORDINATE	SPACE
//D2	CENTROID	SPACE
	if(D1>=D2){
		M_= gsl_matrix_calloc( DIM, D1 );
		bSet_[0] = 1;
		C_= gsl_matrix_calloc( DIM, D2 );
		bSet_[1] = 1;
		vc_	= gsl_vector_calloc( D1 );
		bSet_[2] = 1;
		wc_	= gsl_vector_calloc( D2 );
		bSet_[3] = 1;
	} else {
		std::cout << "ERROR IN ALLOCATION ROUTINE:: "<< D1 << ", " << D2 << std::endl;
	}
}

int
cluster::perform_clustering ( void ){
	if( bSet_[0] && bSet_[1] && bSet_[2] && bSet_[3] )
		return gsl_kmeans(M_, vc_, C_, wc_);
}

void
cluster::seM ( int i , int j , ftyp val){
	if(i<M_->size1&&j<M_->size2&&i>=0&&j>=0)
		gsl_matrix_set(M_,i,j,val);
}

void
cluster::seC ( int i , int j , ftyp val){
	if(i<C_->size1&&j<C_->size2&&i>=0&&j>=0)
		gsl_matrix_set(C_,i,j,val);
}

void
cluster::sew ( int i, ftyp val ) {
	if(i<wc_->size&&i>=0)
		gsl_vector_set(wc_,i,val);
}

void
cluster::sev ( int i, ftyp val ) {
	if(i<vc_->size&&i>=0)
		gsl_vector_set(vc_,i,val);
}

int 
cluster::set_matrix( particles coord_space ) {
	int D = M_->size2;
	if( M_->size2 == coord_space.size() && M_->size1 == DIM ) {
		for(int i=0;i<D; i++){
			for(int j=XX;j<=ZZ;j++){
				seM( j, i, gsl_vector_get(coord_space[i].second,j) );
			}
		}
		gsl_vector_set_all ( vc_ , 1.0 ); // potentially set with label
		return  0;
	}else{
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: calc_distance_matrix" << std::endl;
		return -1;
	}
}


int 
clustering::gsl_kmeans(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw ){
	int NN=dat->size2, MM=dat->size1, KK=cent->size2;

	gvec *labels = gsl_vector_alloc(NN);
	gvec *counts = gsl_vector_alloc(KK);
	gmat *tmp_ce = gsl_matrix_alloc(MM,KK);

	if( ((int)cent->size1)!=MM || !gsl_vector_isnonneg(w) || ((int)w->size)!=NN || ((int)nw->size)!=KK )	
		return -1;
	gsl_vector_set_zero(nw);

	int h, i, j;
	ftyp old_error, error = 1E30, TOL=1E-4; 

	std::vector<int> myvector;
	for (i=0; i<NN; ++i) myvector.push_back(i); 
	random_shuffle ( myvector.begin(), myvector.end() );

	i=0;
	for ( h=0 ; h<KK ; h++ ){
		for ( j=XX ; j<=ZZ ; j++ ){ 
			gsl_matrix_set( cent, j, h , gsl_matrix_get( dat, j, myvector[h] ) );
		} 
	}
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
		gsl_vector_set(w,i,h);			// MIGHT NOT WANT TO OVERWRITE THIS
		gsl_vector_set(nw,h,nwh+wi); 		// NOT NORMALIZED HERE
	}

	gsl_vector_free(labels); 
	gsl_vector_free(counts);
	gsl_matrix_free(tmp_ce);

	return 0;
}

ftyp
node_indices::find_index_relation(cluster c1, cluster c2){
	ftyp min_rmsd	= 1.0E10;

	int N = c1.length_C();

	if(N != c2.length_C())
		std::cout <<"ERROR IN IDX REL"<<std::endl;

	gmat *U		= gsl_matrix_calloc( DIM, DIM );
	gvec *t		= gsl_vector_calloc( DIM );

	int J		= 0;

	gmat *C0T	= gsl_matrix_calloc( DIM, N );
	gmat *C0N	= gsl_matrix_calloc( DIM, N );
	gmat *CNT	= gsl_matrix_calloc( DIM, N );
	gmat *CEN	= gsl_matrix_calloc( DIM, N );
	gvec *gv 	= gsl_vector_calloc( DIM );

	c1.copyC(C0T);
	c2.copyC(C0N);

	std::vector<int> iv;
	for(int i=0;i<N;i++)
		iv.push_back(i);
	std::vector<std::vector<int> > imv = all_permutations(iv);

	for(int j=0;j<imv.size();j++){
		gsl_matrix_memcpy (CEN, C0N);
		for(int i=0;i<iv.size();i++) {
			gsl_matrix_get_col (gv, C0T, i);
			gsl_matrix_set_col (CNT, imv[j][i], gv );
		}
		ftyp rmsd = kabsch_fit(CEN,CNT,U,t);
		if(rmsd<min_rmsd){
			J=j; 
			min_rmsd=rmsd;
		}
	}

	for(int i=0;i<iv.size();i++)
		idx_.push_back(imv[J][i]);

	gsl_matrix_free(CNT);
	gsl_matrix_free(C0N);
	gsl_matrix_free(C0T);
	gsl_matrix_free(U);
	gsl_vector_free(t);
	gsl_vector_free(gv);

	return min_rmsd;
}

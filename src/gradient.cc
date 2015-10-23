/*	gradient.cc	*/
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

#include "richtypes.h"
#include "gradient.hh"
#include <iostream>
#include <algorithm>

using namespace richanalysis;

void
gradient::assign_matrices( particles p1, particles p2 ) {

	if(p1.size()<p2.size()){
		p1.swap(p2);
	}
	int N=p1.size(), M=p2.size();
	// MODEL IS IN p2
	if(M_) 
		gsl_matrix_free(M_);
	M_	= gsl_matrix_alloc(DIM,M);
	if(D_) 
		gsl_matrix_free(D_);
	D_	= gsl_matrix_alloc(DIM,N);
	if(F_) 
		gsl_matrix_free(F_);
	F_	= gsl_matrix_calloc(DIM,M);
	if(Mmu_) 
		gsl_matrix_free(Mmu_);
	Mmu_	= gsl_matrix_calloc(DIM,M); 
	if(Dmu_) 
		gsl_matrix_free(Dmu_);
	Dmu_	= gsl_matrix_calloc(DIM,M);
	if(Msig_) 
		gsl_vector_free(Msig_);
	Msig_	= gsl_vector_calloc(M);
	if(Dsig_) 
		gsl_vector_free(Dsig_);
	Dsig_	= gsl_vector_calloc(M);

	gvec *vt	= gsl_vector_calloc(M);
	gvec *wt	= gsl_vector_calloc(N); 	// HAS A MU INDEX

	std::vector<double> cnts;

	gvec *r1 = gsl_vector_alloc( DIM );
	gvec *r2 = gsl_vector_alloc( DIM );
	for(int i=0;i<M;i++){
		double min_val=1e10;
		for(int j=0;j<M;j++){
			if( i==j )
				continue;
			gsl_vector_memcpy(r1,p2[i].second);
			gsl_vector_memcpy(r2,p2[j].second);
			gsl_vector_sub(r1,r2);
			double len = gsl_blas_dnrm2(r1);
			if( len < min_val ) {
				min_val=sqrt(len);
			}
		}
		gsl_matrix_set_col(M_,i,p2[i].second);
		gsl_vector_set(Msig_,i,min_val/2.0);
		cnts.push_back(0.0);
	}

	for(int i=0;i<N;i++)
		gsl_matrix_set_col(D_,i,p1[i].second);
	gsl_kmeans( D_, wt, Dmu_, vt);

	gsl_vector *r = gsl_vector_calloc(DIM);
	gsl_vector *p = gsl_vector_calloc(DIM);
	for(int i=0;i<N;i++) {
		int j = gsl_vector_get(wt,i);
		gsl_matrix_get_col(r, Dmu_, j);
		gsl_matrix_get_col(p, D_, i);
		gsl_vector_sub	( p, r );
		double len = gsl_blas_dnrm2( p );
		double sig = gsl_vector_get( Dsig_ , j );
		gsl_vector_set(Dsig_,j,sig+len); 
		cnts[j]+=1.0;
	}
	for(int i=0;i<M;i++) {	
		double sig = gsl_vector_get( Dsig_, i );
		gsl_vector_set(Dsig_,i,sqrt( sig/cnts[i]) ); 
	}
	bSet_ = true;

//	OUTPUT
	bool bOut=true;
	if(bOut){
		output_vector(Dsig_);
		output_vector(Msig_);
	}
	gsl_vector_free(vt); gsl_vector_free( r); gsl_vector_free( p);
	gsl_vector_free(wt); gsl_vector_free(r1); gsl_vector_free(r2);
}

/*	richfit.cc	*/
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
#include <sstream>

#include "richfit.hh"
#include <math.h>
#include <algorithm>

using namespace richanalysis;

ftyp 
linalg::v_sum(gvec *v) {
	ftyp s=0.0;
	if(v->size>0)
		for(unsigned int i=0;i<v->size;i++)
			s+=gsl_vector_get(v,i);
	return s;
}

int 
linalg::calc_centroid( gmat *P, gvec *w, gvec *p0  ){
////	WEIGHTED CENTROID WHERE w IS NORMALIZED
	ftyp sum=0.0;
	for( unsigned int i=0 ; i<w->size ; i++)
		sum+=gsl_vector_get(w,i);
	if( p0->size == P->size1 && w->size == P->size2 ) {
		for( unsigned int i=0 ; i<p0->size ; i++ ) {
			ftyp val=0.0;
			for( unsigned int j=0 ; j<w->size; j++ ){
				val+=gsl_matrix_get(P,i,j)*gsl_vector_get(w,j)/sum;
			}
			gsl_vector_set(p0,i,val);
		}
	}else{
		return 1;	
	}	

	return 0;
}

int 
linalg::calc_vmprod( gvec *w, gmat *P, gmat *wP){
	if( w->size == P->size2 && wP->size1 == P->size1 && wP->size2 == P->size2 ) {
		for( unsigned int i=0 ; i<P->size1; i++ ){
			for( unsigned int j=0 ; j<P->size2; j++ )
				gsl_matrix_set(wP, i, j, gsl_matrix_get(P,i,j)*gsl_vector_get(w,j));
		}
		return 0;
	}else{
		return 1;
	}
}

int 
linalg::center_matrix(gmat *P, gvec *p0){
////	CENTER COORDINATES P USING CENTROID p0
	if( P->size1 == p0->size && P->size2>0 ){
		for(unsigned int i=0;i<P->size1;i++)
			for(unsigned int j=0;j<P->size2;j++)
				gsl_matrix_set(P, i, j,
		 (gsl_matrix_get(P, i, j)-gsl_vector_get(p0, i)) );
		return 0;
	}else{
		return 1;
	}
}

int 
linalg::invert_matrix(gmat *A, gmat *invA) {
	int sign;
	gsl_permutation *p 	= gsl_permutation_alloc(A->size1);
	gsl_matrix *tmp_invA 	= gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix *tmpA 	= gsl_matrix_alloc(A->size1, A->size2);	
	gsl_matrix_memcpy(tmpA , A);

	gsl_linalg_LU_decomp (tmpA,p,&sign);
	gsl_linalg_LU_invert (tmpA,p,tmp_invA);
	gsl_matrix_memcpy(invA , tmp_invA);

	gsl_permutation_free(p);
	gsl_matrix_free(tmpA);
	gsl_matrix_free(tmp_invA);

	return sign;
}

ftyp 
linalg::get_det(gmat *A) {
	ftyp det=0.0;
	int signum;

	gsl_permutation *p = gsl_permutation_alloc(A->size1);

	gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix_memcpy(tmpA , A);

	gsl_linalg_LU_decomp(tmpA , p , &signum);
	det = gsl_linalg_LU_det(tmpA , signum);
	gsl_permutation_free(p);
	gsl_matrix_free(tmpA);

	return det;
}

int 
linalg::svd_dec(gmat *r, gmat *u, gvec *s, gmat *v) {
	int M=r->size1,N=r->size2;
	if(!(		r->size1==u->size1 && r->size2==u->size2 
		&&	r->size2==v->size1 && r->size2==v->size2
		&&	s->size ==v->size1 )) {
		std::cout << "ERROR::CANNOT DO SVD" << std::endl;
		std::cout << "ERROR:: DIMS M,N = " << M << ", "<< N << std::endl;
		std::cout << "ERROR:: U(M,N) \t = " <<u->size1 << ", " << u->size2<< std::endl;
		std::cout << "ERROR:: S(N)   \t = " <<s->size << std::endl;
		std::cout << "ERROR:: V(N,N) \t = " << v->size1<< ", " << v->size2<< std::endl;
		exit(1);
	}
	gsl_vector *w = gsl_vector_calloc(N);	
	gsl_matrix_memcpy ( u, r );
	gsl_linalg_SV_decomp( u, v, s, w );
	gsl_vector_free(w);
	return 0;
}

int
linalg::svd_rec(gmat *r, gmat *u, gvec *s, gmat *v) {
	int M=r->size1,N=r->size2;
	if(!(		r->size1==u->size1 && r->size2==u->size2 
		&&	r->size2==v->size1 && r->size2==v->size2
		&&	s->size ==v->size1 )) {
		std::cout << "ERROR::CANNOT DO SVD" << std::endl;
		std::cout << "ERROR:: DIMS M,N = " << M << ", "<< N << std::endl;
		std::cout << "ERROR:: U(M,N) \t = " <<u->size1 << ", " << u->size2<< std::endl;
		std::cout << "ERROR:: S(N)   \t = " <<s->size << std::endl;
		std::cout << "ERROR:: V(N,N) \t = " << v->size1<< ", " << v->size2<< std::endl;
		exit(1);
	}
	gsl_matrix *t = gsl_matrix_calloc(M,N);
	gsl_matrix *z = gsl_matrix_calloc(N,N);

	gsl_matrix_memcpy ( r, u );
	for(int i=0;i<N;i++){
		gsl_matrix_set( z,i,i, gsl_vector_get(s,i));
	}
//	dgemm cannot handle &u=&t
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, u, z, 0.0, t);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, t, v, 0.0, u);
	gsl_matrix_memcpy ( t, u );
	gsl_matrix_memcpy ( u, r );
	gsl_matrix_memcpy ( r, t );
	gsl_matrix_free(t);
	gsl_matrix_free(z);
	return 0;
}


ftyp 
fitting::shape_fit(	gmat *P , gmat *Q ,	// IN
			gvec *w1, gvec *w2, 	// IN
			gmat *U , gvec *t ,	// OUT 
			int II ) {	// CASE VARIABLE (DETERMINANT)

	if( !(	   P->size1 == Q->size1 && P->size2 == w1->size
		&& U->size1 == U->size2 && U->size1 == P->size1 
		&& Q->size1 == U->size1 && Q->size2 == w2->size
		&& t->size  == P->size1 && P->size1 <= P->size2 ) )
		return(-1.0); 				// PROBLEMS WITH THE DIMENSIONS

	int 	L = ((int)P->size2),
		D = ((int)P->size1),
		LL= ((int)Q->size2);			// MATRIX DIMS

	if( !gsl_vector_ispos(w1) || !gsl_vector_ispos(w2) )
		return(-2.0); 				// PROBLEMS WITH VALUE

	ftyp wsum1 = v_sum(w1), wsum2 = v_sum(w2);
	if( wsum1 <= 0.0 || wsum2 <= 0.0 ) {
		return(-2.0); 				// PROBLEMS WITH VALUE
	} else if( !(wsum1 == 1.0 && wsum2 == 1.0) ) {	// THIS SHOULD BE DONE
		gsl_vector_scale(w1,1.0/wsum1);
		gsl_vector_scale(w2,1.0/wsum2);
	}
	gsl_vector *p0 = gsl_vector_alloc(D);
	gsl_vector *q0 = gsl_vector_alloc(D);

	calc_centroid(P,w1,p0); center_matrix(P,p0);
	calc_centroid(Q,w2,q0); center_matrix(Q,q0);

	gsl_matrix *w1P = gsl_matrix_alloc(P->size1,P->size2);
	calc_vmprod( w1, P, w1P );
	gsl_matrix *w2Q = gsl_matrix_alloc(Q->size1,Q->size2);
	calc_vmprod( w2, Q, w2Q );

	gsl_matrix *TMP = gsl_matrix_alloc( D, D );
	gsl_matrix *C   = gsl_matrix_alloc( D, D );
	gsl_matrix *V   = gsl_matrix_alloc( D, D );

	gsl_matrix *EYE = gsl_matrix_alloc( D, D );
	gsl_matrix_set_identity( EYE );
	ftyp det;

//	THE P AND Q INPUTS ARE STORED IN DxN MATRIX NOT NxD
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, w1P, w1P, 0.0, C);

////	THIS IS WHERE THE DECOMPOSITION OCCURS
	gsl_matrix *U1	= gsl_matrix_alloc( D, D );
	gsl_matrix *V1	= gsl_matrix_alloc( D, D );
	gsl_vector *S1	= gsl_vector_alloc( D );
	gsl_vector *wrk1= gsl_vector_alloc( D );

	gsl_linalg_SV_decomp ( C, V1, S1, wrk1 );
	gsl_matrix_memcpy( U1, C );	// HAVE ROTATION

	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, w2Q, w2Q, 0.0, C);
//	AND THE SECOND SET
	gsl_matrix *U2	= gsl_matrix_alloc( D, D );
	gsl_matrix *V2	= gsl_matrix_alloc( D, D );
	gsl_vector *S2	= gsl_vector_alloc( D );
	gsl_vector *wrk2= gsl_vector_alloc( D );

	gsl_linalg_SV_decomp ( C, V2, S2, wrk2 );
	gsl_matrix_memcpy( U2, C );	// HAVE ROTATION

////	ROTATION	FOR THE MODEL (Q)
//	TAKE NOTE OF THE RESULTING SIGN ON U
	gsl_vector_memcpy(t, p0);
	switch(II){
		case  0: break;
		case  1: gsl_matrix_set  ( EYE,  XX, XX, -1.0 ); gsl_matrix_set( EYE, YY, YY, -1.0 ); break;
		case  2: gsl_matrix_set  ( EYE,  XX, XX, -1.0 ); gsl_matrix_set( EYE, ZZ, ZZ, -1.0 ); break;
		case  3: gsl_matrix_set  ( EYE,  YY, YY, -1.0 ); gsl_matrix_set( EYE, ZZ, ZZ, -1.0 ); break;
		case  4: gsl_matrix_set  ( EYE,  XX, XX, -1.0 ); break;
		case  5: gsl_matrix_set  ( EYE,  YY, YY, -1.0 ); break;
		case  6: gsl_matrix_set  ( EYE,  ZZ, ZZ, -1.0 ); break;
		case  7: gsl_matrix_scale( EYE, -1.0); break;
		default: break;
	}
	gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, EYE, V2, 0.0, TMP ); // CblasNoTrans, CblasTrans
	gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, V1, TMP, 0.0, C ); // CblasNoTrans, CblasNoTrans	
	gsl_matrix_memcpy( U , C );
	gsl_blas_dgemv( CblasNoTrans, -1.0, U, q0, 1.0, t );

	gsl_matrix_free(EYE);
	gsl_matrix_free(TMP);
	gsl_matrix_free( C );
	gsl_matrix_free( V );
	gsl_matrix_free(w1P);	gsl_matrix_free(w2Q);
	gsl_matrix_free(U1);	gsl_matrix_free(U2);
	gsl_matrix_free(V1);	gsl_matrix_free(V2);
	gsl_vector_free(S1);	gsl_vector_free(S2);
	gsl_vector_free(wrk1);	gsl_vector_free(wrk2);
	gsl_vector_free(p0);
	gsl_vector_free(q0);

	return (0.0);
}


ftyp 
fitting::shape_fit(	gmat *P , gmat *Q ,	// IN
			gmat *U , gvec *t , int type	) {	
	if( !(	   P->size1 == Q->size1 
		&& U->size1 == U->size2 && U->size1 == P->size1 
		&& Q->size1 == U->size1 
		&& t->size  == P->size1 && P->size1 <= P->size2 ) )
		return(-1.0); 				// PROBLEMS WITH THE DIMENSIONS

	int 	L = ((int)P->size2),
		D = ((int)P->size1),
		LL= ((int)Q->size2);			// MATRIX DIMS

	ftyp wsum1 = P->size2;
	ftyp wsum2 = Q->size2;
	if( wsum1 <= 0.0 || wsum2 <= 0.0 ) {
		return(-2.0); 				// PROBLEMS WITH VALUE
	}

	gsl_vector *w1 = gsl_vector_calloc((int)P->size2);
	gsl_vector *w2 = gsl_vector_calloc((int)Q->size2);
	gsl_vector_set_all( w1 , 1.0/wsum1 );
	gsl_vector_set_all( w2 , 1.0/wsum2 );

	gsl_vector *p0 = gsl_vector_alloc(D);
	gsl_vector *q0 = gsl_vector_alloc(D);

	calc_centroid(P,w1,p0); center_matrix(P,p0);
	calc_centroid(Q,w2,q0); center_matrix(Q,q0);

	gsl_matrix *P_0 = gsl_matrix_alloc(P->size1,P->size2);
	gsl_matrix *P_1 = gsl_matrix_alloc(P->size1,P->size2);
	gsl_matrix *Q_0 = gsl_matrix_alloc(Q->size1,Q->size2);

	gsl_matrix *w1P = gsl_matrix_alloc(P->size1,P->size2);
	calc_vmprod( w1, P, w1P );
	gsl_matrix *w2Q = gsl_matrix_alloc(Q->size1,Q->size2);
	calc_vmprod( w2, Q, w2Q );
	gsl_matrix_memcpy( P_0, w1P );
	gsl_matrix_memcpy( Q_0, w2Q );

	gsl_matrix *DIF	= gsl_matrix_alloc( D, D );
	gsl_matrix *sqP = gsl_matrix_calloc( D, D );
	gsl_matrix *sqQ = gsl_matrix_calloc( D, D );
	gsl_matrix *sqP0 = gsl_matrix_calloc( D, D );
	gsl_matrix *sqQ0 = gsl_matrix_calloc( D, D );

	gsl_matrix *TMP = gsl_matrix_alloc( D, D );
	gsl_matrix *C   = gsl_matrix_alloc( D, D );
	gsl_matrix *V   = gsl_matrix_alloc( D, D );

	gsl_matrix *EYE = gsl_matrix_alloc( D, D );

//	THE P AND Q INPUTS ARE STORED IN DxN MATRIX NOT NxD
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, w1P, w1P, 0.0, C);
	gsl_matrix_memcpy( sqP, C );

////	THIS IS WHERE THE DECOMPOSITION OCCURS
	gsl_matrix *U1	= gsl_matrix_alloc( D, D );
	gsl_vector *S1	= gsl_vector_alloc( D );
	gsl_matrix *V1	= gsl_matrix_alloc( D, D );

	gsl_vector *wrk1= gsl_vector_alloc( D );

	gsl_linalg_SV_decomp ( C, V1, S1, wrk1 ); // HAVE SHAPE
	gsl_matrix_memcpy( U1, C );	

	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, w2Q, w2Q, 0.0, C);
	gsl_matrix_memcpy( sqQ, C );

//	AND THE SECOND SET
	gsl_matrix *U2	= gsl_matrix_alloc( D, D );
	gsl_vector *S2	= gsl_vector_alloc( D );
	gsl_matrix *V2	= gsl_matrix_alloc( D, D );
	gsl_vector *wrk2= gsl_vector_alloc( D );

	gsl_linalg_SV_decomp ( C, V2, S2, wrk2 );// HAVE SHAPE
	gsl_matrix_memcpy( U2, C ); 

////	ROTATION	FOR THE MODEL (Q)
//	TAKE NOTE OF THE RESULTING SIGN ON U
	gsl_matrix *Ut = gsl_matrix_calloc( D, D );
	gsl_vector *tt = gsl_vector_calloc( D );

	ftyp cnt=(ftyp)DIM;

	ftyp rmsd	= 0.0 , min_rmsd	= type==2?1e10:0.0;
	ftyp total_rmsd	= 0.0 , total_min_rmsd	= 1e10;
	ftyp xi 	= 1.0 , gau_norm 	= 1.0/sqrt(2.0*M_PI)/xi;

	double dval,dv1,dv2;
	dv1	= get_det(V1);
	dv2	= get_det(V2);
	dval	= dv1*dv2;
		
	for(int II=0; II<8; II++)
	{
		gsl_vector_memcpy( tt, p0 );
		gsl_matrix_set_identity( EYE );
		gsl_matrix_memcpy( P_0, w1P );
		gsl_matrix_memcpy( Q_0, w2Q );
		gsl_matrix_memcpy( sqP0, sqP );
		gsl_matrix_memcpy( sqQ0, sqQ );

		switch(II){
			case  0: break;
			case  1: gsl_matrix_set  ( EYE,  XX, XX, -1.0 ); gsl_matrix_set( EYE, YY, YY, -1.0 ); break;
			case  2: gsl_matrix_set  ( EYE,  XX, XX, -1.0 ); gsl_matrix_set( EYE, ZZ, ZZ, -1.0 ); break;
			case  3: gsl_matrix_set  ( EYE,  YY, YY, -1.0 ); gsl_matrix_set( EYE, ZZ, ZZ, -1.0 ); break;
			case  4: gsl_matrix_set  ( EYE,  XX, XX, -1.0 ); break;
			case  5: gsl_matrix_set  ( EYE,  YY, YY, -1.0 ); break;
			case  6: gsl_matrix_set  ( EYE,  ZZ, ZZ, -1.0 ); break;
			case  7: gsl_matrix_scale( EYE, -1.0); break;
			case  8: // for playing around
				gsl_matrix_set  ( EYE,  XX, XX, dv2 );
				gsl_matrix_set  ( EYE,  YY, YY, dv1 );
				gsl_matrix_set  ( EYE,  ZZ, ZZ, dval ); 
				break;
			default: break;
		}
//	BELOW CORRESPONDS TO WAHBA	M
		if( type<0 ) {
			gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, EYE, V2, 0.0, TMP ); 
			gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, V1, TMP, 0.0, C ); 
		} else {	
//	BELOW CORRESPONDS TO KABSCH	M^T
			gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, EYE, V1, 0.0, TMP ); 
			gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, V2, TMP, 0.0, C ); 	
		}
		int dC = get_det(C);

		gsl_matrix_memcpy( Ut , C );
		gsl_blas_dgemv( CblasNoTrans, -1.0, Ut, q0, 1.0, tt );

		switch(type) {
			case 2:	{
			gsl_matrix_memcpy( DIF, sqQ0 );
			gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, Ut, sqP0, -1.0, DIF );
			rmsd=0.0;
			for( int i=0 ; i<DIM ; i++ ){
				rmsd += (	square(gsl_matrix_get(DIF,XX,i))
					+	square(gsl_matrix_get(DIF,YY,i))
					+	square(gsl_matrix_get(DIF,ZZ,i))  );
			}
			rmsd /= cnt;
			if( rmsd < min_rmsd ) {
				gsl_matrix_memcpy( U , Ut );
				gsl_vector_memcpy( t , tt );
				min_rmsd = rmsd;
			}
				} break;
			case 1: {
			gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, Ut, P_0, 0.0, P_1 );
			ftyp tot_cnt	= 0.0;
			min_rmsd 	= 0.0;
			rmsd		= 0.0;
			for(int i=0;i<P_1->size2;i++) {
				for(int j=0;j<Q_0->size2;j++) {
					double x2	 = square( gsl_matrix_get(P_1,XX,i)-gsl_matrix_get(Q_0,XX,j) );
					double y2	 = square( gsl_matrix_get(P_1,YY,i)-gsl_matrix_get(Q_0,YY,j) );
					double z2	 = square( gsl_matrix_get(P_1,ZZ,i)-gsl_matrix_get(Q_0,ZZ,j) );
					double r2	 = x2+y2+z2;
					min_rmsd	+= r2*r2*r2 ;
					tot_cnt		+= 1.0;
				}
			}
			min_rmsd/=tot_cnt;
			if(  min_rmsd < total_min_rmsd  ) {
				gsl_matrix_memcpy( U , Ut );
				gsl_vector_memcpy( t , tt );
				total_min_rmsd = min_rmsd;
			}
				}; break;
			default:{
			gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, Ut, P_0, 0.0, P_1 );
			ftyp tot_cnt	= 0.0;
			min_rmsd	= 0.0;
			rmsd		= 0.0;
			for(int i=0; i<P_1->size2;i++ ) {	
				for(int j=0;j<Q_0->size2;j++) {
					double x2	= square( gsl_matrix_get(P_1,XX,i)-gsl_matrix_get(Q_0,XX,j) );
					double y2	= square( gsl_matrix_get(P_1,YY,i)-gsl_matrix_get(Q_0,YY,j) );
					double z2	= square( gsl_matrix_get(P_1,ZZ,i)-gsl_matrix_get(Q_0,ZZ,j) );
					double r2	= (x2+y2+z2);
					min_rmsd	+= r2*r2*r2 ;
					rmsd		+= r2;
					tot_cnt		+= 1.0;
				}
			}
			min_rmsd/=tot_cnt;
			if(  min_rmsd < total_min_rmsd ) {
				gsl_matrix_memcpy( U , Ut );
				gsl_vector_memcpy( t , tt );
				total_min_rmsd = min_rmsd;
			}			
				} break;
		}	
	}

	gsl_matrix_free(EYE);	gsl_matrix_free(sqP);
	gsl_matrix_free(TMP);	gsl_matrix_free(sqQ);
	gsl_matrix_free( C );	gsl_matrix_free(DIF);
	gsl_matrix_free( V );
	gsl_matrix_free(Ut);	gsl_vector_free(tt);
	gsl_matrix_free(w1P);	gsl_matrix_free(w2Q);
	gsl_matrix_free(U1);	gsl_matrix_free(U2);
	gsl_matrix_free(V1);	gsl_matrix_free(V2);
	gsl_vector_free(S1);	gsl_vector_free(S2);
	gsl_vector_free(wrk1);	gsl_vector_free(wrk2);
	gsl_vector_free(p0);	gsl_vector_free(w1);
	gsl_vector_free(q0);	gsl_vector_free(w2);

	return min_rmsd;
}

void
fitting::invert_fit( gmat *U, gvec *t, gmat *invU, gvec *invt ) {
	if( invU->size1!=U->size1 || invU->size2!=U->size2 || t->size!=invt->size)
		std::cout << "ERROR IN DIMENSIONS"<< std::endl;
	invert_matrix(U,invU);
	gsl_blas_dgemv(CblasNoTrans, -1.0, invU, t, 0.0, invt);	
}


ftyp 
fitting::kabsch_fit(	gsl_matrix *P, gsl_matrix *Q,			// IN
			gsl_matrix *U, gsl_vector *t ) {		// OUT
//// NOTE ON IO
//	IN 	P(D,M), Q(D,N), w(M)	
//	OUT	U(D,N), t(D), returns rmsd>0 or error<0
	if( !(	   P->size1 == Q->size1 && P->size2 == Q->size2
		&& U->size1 == U->size2 && U->size1 == P->size1 
		&& t->size  == P->size1 ) )
		return(-1.0); 				// PROBLEMS WITH THE DIMENSIONS
	int L = ((int)P->size2), D = ((int)P->size1);	// STORED AS COLUMN MATRIX
	int LL= ((int)Q->size2);			// L CANNOT BE LARGER THAN LL

	gsl_vector *p0 = gsl_vector_alloc(D);
	gsl_vector *q0 = gsl_vector_alloc(D);
	gsl_vector *w  = gsl_vector_alloc(((int)P->size2));

	double wsum=P->size2;
	gsl_vector_set_all(w,1.0/wsum);

	calc_centroid(P,w,p0); center_matrix(P,p0);
	calc_centroid(Q,w,q0); center_matrix(Q,q0);

	gsl_matrix *C	= gsl_matrix_alloc(D,D);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, P, Q, 0.0, C);

	gsl_matrix *W	= gsl_matrix_alloc( D, D );
	gsl_matrix *V	= gsl_matrix_alloc( D, D );
	gsl_vector *S	= gsl_vector_alloc( D );
	gsl_vector *work= gsl_vector_alloc( D );

	gsl_linalg_SV_decomp( C, V, S, work );
	gsl_matrix_memcpy( W, C );
	
	gsl_matrix *EYE = gsl_matrix_alloc( D, D );
	gsl_matrix *TMP = gsl_matrix_alloc( D, D );
	
	gsl_matrix_set_identity( EYE );
	
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, W, V, 0.0, C);
	double det = get_det(C);
	if (det < 0){	// FLIP IT!
		gsl_matrix_set(EYE,ZZ,ZZ,-1); 
		gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, EYE, V, 0.0, TMP);
		gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, W, TMP, 0.0, C);
	}
	//gsl_matrix_transpose_memcpy(U,C);
	gsl_matrix_memcpy(U,C);

	gsl_blas_dgemv(CblasNoTrans, -1.0, U, p0, 1.0, q0);
	gsl_vector_memcpy(t,q0); 

	gsl_matrix *DIFF = gsl_matrix_alloc( D, L );

	double rmsd = 0.0;
	gsl_matrix_memcpy( DIFF, Q );
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, U, P, -1.0, DIFF);
	for( int i=0 ; i<L ; i++ ){
		rmsd += (  square(gsl_matrix_get(DIFF,XX,i))
			 + square(gsl_matrix_get(DIFF,YY,i))
			 + square(gsl_matrix_get(DIFF,ZZ,i)) );
	}
	rmsd/=wsum;

	gsl_matrix_free(V); gsl_matrix_free(W); gsl_matrix_free(TMP); gsl_matrix_free(EYE); gsl_matrix_free(C);
	gsl_vector_free(p0); gsl_vector_free(q0); gsl_vector_free(w); gsl_vector_free(S);gsl_vector_free(work);

	return sqrt(rmsd);
}

ftyp 
fitting::kabsch_fit(	gsl_matrix *P, gsl_matrix *Q, gsl_vector *w,	// IN
			gsl_matrix *U, gsl_vector *t ) {		// OUT
//// NOTE ON IO
//	IN 	P(D,M), Q(D,N), w(M)	
//	OUT	U(D,N), t(D), returns rmsd>0 or error<0
	if( !(	   P->size1 == Q->size1 && P->size2 == Q->size2
		&& U->size1 == U->size2 && U->size1 == P->size1 
		&& P->size2 == w->size  && t->size  == P->size1 ) )
		return(-1.0); 				// PROBLEMS WITH THE DIMENSIONS
	int L = ((int)P->size2), D = ((int)P->size1);	// STORED AS COLUMN MATRIX
	int LL= ((int)Q->size2);			// L CANNOT BE LARGER THAN LL
	if( !gsl_vector_ispos(w) )
		return(-2.0); 				// PROBLEMS WITH VALUE
	double wsum=v_sum(w);
	if( wsum <= 0.0 )
		return(-2.0); 				// PROBLEMS WITH VALUE
	gsl_vector_scale(w,1.0/wsum);
	gsl_vector *p0 = gsl_vector_alloc(D);
	gsl_vector *q0 = gsl_vector_alloc(D);

	calc_centroid(P,w,p0); center_matrix(P,p0);
	calc_centroid(Q,w,q0); center_matrix(Q,q0);

	gsl_matrix *C	= gsl_matrix_alloc(D,D);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, P, Q, 0.0, C);

	gsl_matrix *W	= gsl_matrix_alloc( D, D );
	gsl_matrix *V	= gsl_matrix_alloc( D, D );
	gsl_vector *S	= gsl_vector_alloc( D );
	gsl_vector *work= gsl_vector_alloc( D );

	gsl_linalg_SV_decomp( C, W, S, work );
	gsl_matrix_memcpy( V, C );
	
	gsl_matrix *EYE = gsl_matrix_alloc( D, D );
	gsl_matrix *TMP = gsl_matrix_alloc( D, D );
	
	gsl_matrix_set_identity( EYE );
	
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, V, W, 0.0, C);
	double det = get_det(C);

	if (det < 0){	// FLIP IT!
		gsl_matrix_set(EYE,D-1,D-1,-1);
		gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, EYE, W, 0.0, TMP);
		gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, V, TMP, 0.0, C);
	}
	gsl_matrix_transpose_memcpy(U,C);

	gsl_blas_dgemv(CblasNoTrans, -1.0, U, p0, 1.0, q0);
	gsl_vector_memcpy(t,q0); 

	gsl_matrix *DIFF = gsl_matrix_alloc( D, L );

	double rmsd = 0.0;

	gsl_matrix_memcpy( DIFF, Q );
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, U, P, -1.0, DIFF);
	for( int i=0 ; i<L ; i++ ){
		rmsd += gsl_vector_get(w,i)*( 
				   square(gsl_matrix_get(DIFF,XX,i))
				 + square(gsl_matrix_get(DIFF,YY,i))
				 + square(gsl_matrix_get(DIFF,ZZ,i)) );
	}
	gsl_matrix_free(V); gsl_matrix_free(W); gsl_matrix_free(TMP); gsl_matrix_free(EYE); gsl_matrix_free(C);
	gsl_vector_free(p0);gsl_vector_free(q0); gsl_vector_free(S);gsl_vector_free(work);

	return sqrt(rmsd);
}	

void 
fitting::apply_fit( gmat *M, gmat *U, gvec *t) {
	gsl_vector *pos = gsl_vector_alloc( 3 );
	gsl_vector *res = gsl_vector_alloc( 3 );
	int n = M->size1, m=M->size2;
	int L = n>m?n:m;

	for( int i=0 ; i<L ; i++ ){
		if(n>m){
			gsl_vector_set( pos, XX , gsl_matrix_get( M, i, XX ) );
			gsl_vector_set( pos, YY , gsl_matrix_get( M, i, YY ) );
			gsl_vector_set( pos, ZZ , gsl_matrix_get( M, i, ZZ ) );
		}else{
			gsl_vector_set( pos, XX , gsl_matrix_get( M, XX, i ) );
			gsl_vector_set( pos, YY , gsl_matrix_get( M, YY, i ) );
			gsl_vector_set( pos, ZZ , gsl_matrix_get( M, ZZ, i ) );
		}

		gsl_blas_dgemv( CblasNoTrans, 1.0, U, pos, 0.0, res );

		if(n>m){
			gsl_matrix_set(M,i,XX,gsl_vector_get(res,XX)+gsl_vector_get(t,XX));
			gsl_matrix_set(M,i,YY,gsl_vector_get(res,YY)+gsl_vector_get(t,YY));
			gsl_matrix_set(M,i,ZZ,gsl_vector_get(res,ZZ)+gsl_vector_get(t,ZZ));
		}else{
			gsl_matrix_set(M,XX,i,gsl_vector_get(res,XX)+gsl_vector_get(t,XX));
			gsl_matrix_set(M,YY,i,gsl_vector_get(res,YY)+gsl_vector_get(t,YY));
			gsl_matrix_set(M,ZZ,i,gsl_vector_get(res,ZZ)+gsl_vector_get(t,ZZ));
		}
	}
	gsl_vector_free(pos);
	gsl_vector_free(res);
}

void 
fitting::apply_fit( particles px, gmat *U, gvec *t) {
	gsl_vector *pos = gsl_vector_calloc( DIM );
	gsl_vector *res = gsl_vector_calloc( DIM );
	gsl_vector *com = gsl_vector_calloc( DIM );
	int L = ((int)px.size());

	for( int i=0 ; i<L ; i++ ) {
		gsl_vector_set( pos, XX, gsl_vector_get( px[i].second,XX ) );
		gsl_vector_set( pos, YY, gsl_vector_get( px[i].second,YY ) );
		gsl_vector_set( pos, ZZ, gsl_vector_get( px[i].second,ZZ ) );
		
		gsl_blas_dgemv(CblasNoTrans, 1.0, U, pos, 0.0, res);

		gsl_vector_set( px[i].second , XX, gsl_vector_get(res,XX)+gsl_vector_get(t,XX) );
		gsl_vector_set( px[i].second , YY, gsl_vector_get(res,YY)+gsl_vector_get(t,YY) );
		gsl_vector_set( px[i].second , ZZ, gsl_vector_get(res,ZZ)+gsl_vector_get(t,ZZ) );
	}
	gsl_vector_free(pos);
	gsl_vector_free(res);
	gsl_vector_free(com);
}



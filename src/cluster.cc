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
	if( D1 > 0 && D2 > 0 ) {
		M_	= gsl_matrix_calloc( DIM, D1 );	// MODEL
		C_	= gsl_matrix_calloc( DIM, D2 );	// DATA
		vc_	= gsl_vector_calloc( D1 );
		wc_	= gsl_vector_calloc( D2 );
		bSet_   = true;
	} else {
		std::cout << "ERROR IN ALLOCATION ROUTINE:: "<< D1 << ", " << D2 << std::endl;
	}
}

int 
cluster::set_matrix( particles coord_space ) {
	int D		= coord_space.size();
	int reval	=	-1;

	int rD = get_cDIM();
	if( bSet_ == 0 ) {
		alloc_space(D,rD);
	}
	if( M_->size2 == coord_space.size() && M_->size1 == DIM && bSet_ ) {
		id inInfo;
		idlabels_.clear();
		for(int i=0;i<D; i++) {
			inInfo.first	= i;
			inInfo.second	= coord_space[i].first;
			idlabels_.push_back ( inInfo );
			gsl_matrix_set_col  ( M_ , i, coord_space[i].second );
		}
		gsl_vector_set_all ( vc_ , 1.0 );
		reval =0; 
	} else {
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: set_matrix" 
			<< D << " " << rD << " " << M_->size2 << std::endl;
	}
	if(rD>1)
		find_centroids();

	return  reval;
}

int
cluster::find_centroids ( void ){
	if( bSet_ ) {
		gsl_kmeans( M_, vc_, C_, wc_ );
		int N=wc_->size;
		int M=vc_->size;
		NperC_.clear();
		for(int i=0;i<N;i++){
			int numi=0;
			for( int j=0 ; j<vc_->size ; j++ ) {
				numi+=(gsl_vector_get(vc_,j)==i)?1:0;
			}
			NperC_.push_back(numi);
		}
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

bool
node_analysis::assign_node( node n ) {

	n.first.find_centroids();
	n.second.find_centroids();
	cluster c1 = n.first;
	cluster c2 = n.second;
	std::vector<int> cindx;
	
	int N, M, L, K;
	int rDIM=c1.get_cDIM();

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

	for(int j=0;j<imv.size();j++){
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

	if(idx_.size()!=0) {
		std::cout << "INFO::ERROR WITH IDX_" << std::endl;
		idx_.clear();
		iidx_.clear();
	}
	for(int i=0;i<iv.size();i++)
		iidx_.push_back(0);

	for(int i=0;i<iv.size();i++){
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


bool
node_analysis::find_index_orders( void )
{
	std::vector<int> at_order, tmp_order, cl_order;

	glob_idx_order_.clear();

	if ( idx_.size() == children_.size() ) 
	{
		// std::cout<< "DOING THIS1" << std::endl;
		cidx_.clear();
		for ( int i=0 ; i<idx_.size() ; i++ ) 
		{
			for ( int j=0 ; j<children_[i].second.length_M() ; j++ ) 
			{
				cidx_.push_back(	iidx_[i] 	);
				cl_order.push_back(	iidx_[i] 	);
			}
		}
	}

	if( cidx_.size() == parents_.second.length_M() )
	{
/*
		// std::cout << "DOING THIS2" << std::endl;
		gsl_vector *gv = gsl_vector_calloc( parents_.second.length_M() );
		parents_.second.copyv(gv);

		int N=gv->size;
		for(int i=0;i<N;i++) 
		{
			tmp_order.push_back( (int)gsl_vector_get(gv,i) );
		}

		for(int i=0;i<N;i++)
		{
			int c = tmp_order[i];
			for(int j=0;j<N;j++) 
			{
				int cj = ( cidx_[j]>=0 )?(iidx_[cidx_[j]]):(-1);
				if( c == cj ) 
				{
					at_order.push_back(j);
					glob_idx_order_.push_back(j);
					cidx_[j]=(cidx_[j]+1)*(-1);
					break;
				}
			}
		}

		//cl_order.swap(		cidx_	);
		//at_order.swap( glob_idx_order_ 	);
*/
		return true;
	}

	return false;
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

void
cluster0::alloc_space( int D1, int D2 ) {
//D1	COORDINATE	SPACE
//D2	CENTROID	SPACE
//	std::cout << "INFO::ALLOC::INFO " << D1 << " " << D2 << " " << DIM << std::endl;
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
cluster0::perform_clustering ( void ){
	if( bSet_[0] && bSet_[1] && bSet_[2] && bSet_[3] ) {
		gsl_kmeans(M_, vc_, C_, wc_);
		int N=wc_->size;
		int M=vc_->size;
		NperC_.clear();
		for(int i=0;i<N;i++){
			int numi=0;
			for( int j=0 ; j<vc_->size ; j++ ) {
				numi+=(gsl_vector_get(vc_,j)==i)?1:0;
			}
			NperC_.push_back(numi);
		}
	}
}

int
cluster0::perform_clustering ( ftyp cutoff ){
	if( bSet_[0] && bSet_[1] && bSet_[2] && bSet_[3] ) {
		gsl_kmeans(M_, vc_, C_, wc_, cutoff);
		int N=wc_->size;
		int M=vc_->size;
		NperC_.clear();
		for(int i=0;i<N;i++){
			int numi=0;
			for( int j=0 ; j<vc_->size ; j++ ) {
				numi+=(gsl_vector_get(vc_,j)==i)?1:0;
			}
			NperC_.push_back(numi);
		}
	}
}

void
cluster0::seM ( int i , int j , ftyp val){
	if(i<M_->size1&&j<M_->size2&&i>=0&&j>=0)
		gsl_matrix_set(M_,i,j,val);
}

void
cluster0::seC ( int i , int j , ftyp val){
	if(i<C_->size1&&j<C_->size2&&i>=0&&j>=0)
		gsl_matrix_set(C_,i,j,val);
}

void
cluster0::sew ( int i, ftyp val ) {
	if(i<wc_->size&&i>=0)
		gsl_vector_set(wc_,i,val);
}

void
cluster0::sev ( int i, ftyp val ) {
	if(i<vc_->size&&i>=0)
		gsl_vector_set(vc_,i,val);
}

std::vector<int>	
cluster0::get_labels( void ){
	std::vector<int> ndx;
	if( bSet_[0] && bSet_[1] && bSet_[2] && bSet_[3] )
		for(int i=0;i<vc_->size;i++)
			ndx.push_back(gsl_vector_get(vc_,i));
	return ndx;
}

int 
cluster0::set_matrix( particles coord_space ) {
	int D = M_->size2;
	if( M_->size2 == coord_space.size() && M_->size1 == DIM ) {
		for(int i=0;i<D; i++){
			for(int j=XX;j<=ZZ;j++){
				seM( j, i, gsl_vector_get(coord_space[i].second,j) );
			}
		}
		gsl_vector_set_all ( vc_ , 1.0 ); // potentially set with label
		bSet_[0] = 1; bSet_[1] = 1;
		return  0;
	}else{
		std::cout << "ERROR:: BAD CONTAINER LENGTHS:: calc_distance_matrix" << std::endl;
		return -1;
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
clustering::gsl_kmeans(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw ){
	int NN=dat->size2, MM=dat->size1, KK=cent->size2;

	gvec *labels = gsl_vector_alloc(NN);
	gvec *counts = gsl_vector_alloc(KK);
	gmat *tmp_ce = gsl_matrix_alloc(MM,KK);

	if( ((int)cent->size1)!=MM || !gsl_vector_isnonneg(w) || ((int)w->size)!=NN || ((int)nw->size)!=KK )	
		return -1;
	gsl_vector_set_zero(nw);

	int h, i, j;
	ftyp old_error, error = 1E30, TOL=1E-8; 

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


int 
clustering::gsl_kmeans(gmat *dat, gsl_vector *w, gmat *cent, gsl_vector *nw, ftyp cutoff){
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
				if ( distance < min_distance   ) { // && distance > square(cutoff)
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

std::vector<int>			
cluster_node0::find_simple_relation( void ){
	std::vector<std::pair<int,int > > residx;
	std::vector<int> idx, sidx;
	std::vector<ftyp > rmsds;

	if(!subSet()) {
		std::cout << "INFO::SUBNOTSET::" << std::endl;
		return idx;
	}
	int N = vc1_.size();
	int M = vc2_.size();
	std::cout << "INFO:: HAVE SUBCLUSTERS OF SIZE::" << N << " AND " << M << std::endl;
	if(N!=M) {
		return idx;
	}

	for(int i=0;i<N;i++) {
	double rmsd0=100.0;
	for(int j=0;j<M;j++) {
		gmat *P = gsl_matrix_calloc(DIM,vc2_[i].length_M());
		vc2_[i].copyM(P);
		gmat *Q = gsl_matrix_calloc(DIM,vc1_[j].length_M());
		vc1_[j].copyM(Q);
		gmat *U = gsl_matrix_calloc(DIM,DIM);
		gvec *t = gsl_vector_calloc(DIM);

		ftyp rmsd=shape_fit( P, Q, U, t, 2);

		std::pair<int,int> pi;
		pi.first=i; pi.second=j;
		rmsds.push_back(rmsd);
		residx.push_back(pi);
	}
	}

	idx.push_back(-1); idx.push_back(-1); idx.push_back(-1);
	sidx.push_back(-1); sidx.push_back(-1); sidx.push_back(-1);
	while( N>0 ) {
		ftyp rms	= 1E10; int J=0;
		for(int i=0;i<residx.size();i++){
			if( rms>rmsds[i] ) {
				rms=rmsds[i];
				J	= i;
			}
		}
		rmsds[J]=1E10;
		if( idx[residx[J].first]==-1 && sidx[residx[J].second]==-1 ){ // 
			idx[residx[J].first]=residx[J].second;
			sidx[residx[J].second]=residx[J].first;
			N--;
		}
	}
	std::cout << "INFO::HAVE::N_IDX::" << idx.size() << std::endl;

	idx_.clear();
	for(int i=0;i<idx.size();i++){
		idx_.push_back(idx[i]);
		std::cout << "INFO:: " << idx[i] << std::endl;
	}
	return idx;
}

void
cluster_node0::assign_sub( cluster0 c1, cluster0 c2 ) {

	int N,M;
	if(c1.isSet()&&c2.isSet()){
		parents_.first  = c1;
		parents_.second = c2;

		N=c1.length_M();
		M=c1.length_C();
		gmat *M1 = gsl_matrix_calloc(DIM,N);
		gvec *v1 = gsl_vector_calloc(N);
		c1.copyM(M1); c1.copyv(v1);
		for( int ipart=0; ipart<M; ipart++ ){
			richanalysis::coord_format cf;
			particles px;
			px	= cf.truncmat(M1,v1,ipart);
			int D	= px.size();
			richanalysis::cluster0 clpx;
			clpx.alloc_space( D, 3 );	// THESE CENTOIDS SHOULD NOT BE USED
			clpx.set_matrix( px );
			clpx.perform_clustering();
			vc1_.push_back(clpx);
		}

		N=c2.length_M();
		M=c1.length_C();
		gmat *M2 = gsl_matrix_calloc(DIM,N);
		gvec *v2 = gsl_vector_calloc(N);
		c2.copyM(M2); c2.copyv(v2);
		for( int ipart=0; ipart<M; ipart++ ){
			richanalysis::coord_format cf;
			particles px;
			px	= cf.truncmat(M2,v2,ipart);
			int D	= px.size();
			richanalysis::cluster0 clpx;
			clpx.alloc_space( D, 3 );	// THESE CENTOIDS SHOULD NOT BE USED
			clpx.set_matrix( px );
			clpx.perform_clustering();
			vc2_.push_back(clpx);
		}
		subSet_=1;

	}
}

particles 
cluster_node0::apply_rot_trans( particles px , int I){

	particles rpx;
	for(int i=0;i<px.size();i++){
		particle p0;
		p0.first	= px[i].first;
		p0.second	= gsl_vector_calloc(DIM);
		gsl_vector_memcpy( p0.second, px[i].second);
		rpx.push_back( p0 );
	}

	if(bUtSet_){
		if(I>=0)
			apply_fit( rpx,  U_, t_	 );
		if(I<0)
			apply_fit( rpx, iU_, it_ );
	}

	return rpx;
}

particles 
cluster_node0::apply_rot_trans( particles px ){
	particles rpx;
	for(int i=0;i<px.size();i++){
		particle p0;
		p0.first	= px[i].first;
		p0.second	= gsl_vector_calloc(DIM);
		gsl_vector_memcpy( p0.second, px[i].second);
		rpx.push_back( p0 );
	}

	if(bUtSet_){
		if(sgn_>=0)
			apply_fit( rpx,  U_, t_	 );
		if(sgn_<0)
			apply_fit( rpx, iU_, it_ );
	}
	return rpx;
}

void
cluster_node0::invert_transform( void ){
	if( iU_->size1==U_->size1 && iU_->size2==U_->size2 && t_->size==it_->size){
		invert_matrix(U_,iU_);
		gsl_blas_dgemv(CblasNoTrans, -1.0, iU_, t_, 0.0, it_);
	}
}

std::vector<int> 
cluster_node0::find_centroid_relation( void ){
	ftyp min_rmsd	= 1.0E10;

	cluster0 c1 = parents_.first;
	cluster0 c2 = parents_.second;

	int N = c1.length_C();

	if( N != c2.length_C())
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

	for(int j=0;j<imv.size();j++){
		gsl_matrix_memcpy (CEN, C0N);
		for(int i=0;i<iv.size();i++) {
			gsl_matrix_get_col (gv, C0T, i);
			gsl_matrix_set_col (CNT, imv[j][i], gv );
		}
		ftyp rmsd = kabsch_fit(CNT,CEN,U,t);
		if(rmsd<min_rmsd){
			J=j; 
			min_rmsd=rmsd;	
			gsl_matrix_memcpy (Ut, U);
			gsl_vector_memcpy (tt, t);
		}
	}

	if(idx_.size()!=0) {
		std::cout << "INFO::ERROR WITH IDX_" << std::endl;
		idx_.clear();
		iidx_.clear();
	}
	for(int i=0;i<iv.size();i++)
		iidx_.push_back(0);

	for(int i=0;i<iv.size();i++){
		idx_.push_back(imv[J][i]);
		iv[i]=imv[J][i];
		iidx_[imv[J][i]]=i;
	}

	bDirRel_=2*(c1.length_M()>c2.length_M())-1;

	if(!have_transform())
	{
		gsl_matrix_memcpy (U_, Ut);
		gsl_vector_memcpy (t_, tt);
		invert_transform();
		bUtSet_=1;
		sgn_=1;
	}

	gsl_matrix_free(CNT);
	gsl_matrix_free(C0N);
	gsl_matrix_free(C0T);
	gsl_matrix_free(U);
	gsl_vector_free(t);
	gsl_vector_free(gv);

	return iv;
}

ftyp
cluster0::find_shape(){
//	Do a single shape fit

	if( isSet() ) { 
		int N 		= length_M();
		gmat *MEN	= gsl_matrix_calloc( DIM, N );
		copyM(MEN);
	
		gmat *Utmp	= gsl_matrix_calloc( DIM, DIM );
		gvec *ttmp	= gsl_vector_calloc( DIM );
		
		ftyp min_rmsd 	= shape(MEN,Utmp,ttmp);

		gsl_matrix_memcpy (Uc_, Utmp);
		gsl_vector_memcpy (tc_, ttmp);
		bPCset_		=	1;

		gsl_matrix_free(MEN);
		gsl_matrix_free(Utmp); 
		gsl_vector_free(ttmp);
	}

	return 0.0;
}

ftyp 
cluster0::shape(		gmat *P, gmat *U , gvec *t 	) 
{	
	if( !(	   P->size1 
		&& U->size1 == U->size2 && U->size1 == P->size1 
		&& t->size  == P->size1 && P->size1 <= P->size2 ) )
		return(-1.0); 		// PROBLEMS WITH THE DIMENSIONS

	int 	L = ((int)P->size2),
		D = ((int)P->size1);	// MATRIX DIMS

	ftyp wsum = P->size2;
	if( wsum <= 0.0 ) {
		return(-2.0); 		// PROBLEMS WITH VALUE
	}

	gsl_vector *w = gsl_vector_alloc((int)P->size2);
	gsl_vector_set_all( w , 1.0/wsum );
	gsl_vector *p0 = gsl_vector_alloc(D);

	calc_centroid(P,w,p0); center_matrix(P,p0);
	gsl_matrix *wP = gsl_matrix_alloc(P->size1,P->size2);
	calc_vmprod( w, P, wP );

	gsl_matrix *W	= gsl_matrix_alloc( L, D );
	gsl_matrix *V   = gsl_matrix_alloc( D, D );
	gsl_matrix *C   = gsl_matrix_alloc( L, D );
	gsl_matrix_transpose_memcpy (C, wP);
	gsl_vector *S	= gsl_vector_alloc( D );
	gsl_vector *wrk = gsl_vector_alloc( D );

	gsl_matrix *TMP = gsl_matrix_alloc( D, D );
	gsl_matrix *EYE = gsl_matrix_alloc( D, D );

	gsl_linalg_SV_decomp ( C, V, S, wrk );
	gsl_matrix_memcpy( W, C );

	output_vector(S);

	gsl_matrix_transpose_memcpy(U,V);
	gsl_vector_memcpy(t, p0);

	gsl_vector_free(w);
	gsl_vector_free(p0);	
	gsl_matrix_free( wP);	
	gsl_matrix_free( W );	
	gsl_matrix_free( V );
	gsl_matrix_free( C );
	gsl_vector_free( S );	
	gsl_vector_free(wrk);	
	gsl_matrix_free(TMP);		
	gsl_matrix_free(EYE);	

	return 0.0;
}

particle
cluster0::normal( void ){
	particle par;
	if( has_shape() ){
		gmat *Uc1	= gsl_matrix_calloc(DIM,DIM);
		gvec *tc	= gsl_vector_calloc(DIM);
		copyUc(Uc1);
		copytc(tc);
		gvec *a = gsl_vector_calloc(DIM);
 		gsl_matrix_get_row (a, Uc1, ZZ);
		par.first ="Ag";
		gsl_vector_add(a,tc);
		par.second=a;
		return par;
	}else{
		return par;
	}
}
particle
cluster0::center( void ){
	if( has_shape() ){
		gvec *tc	= gsl_vector_calloc(DIM);
		copytc(tc);
		particle par;
		par.first ="Ag";
		par.second=tc;
		return par;
	}
}

ftyp			
cluster_node0::angle_between(cluster0 c1, cluster0 c2){

	if( c1.has_shape() && c2.has_shape() ){
		ftyp a,b,c,x,y,z;
		gmat *Uc1 = gsl_matrix_calloc(DIM,DIM);
		c1.copyUc(Uc1);
		a = gsl_matrix_get(Uc1,XX,ZZ); 
		b = gsl_matrix_get(Uc1,YY,ZZ); 
		c = gsl_matrix_get(Uc1,ZZ,ZZ);
		gmat *Uc2 = gsl_matrix_calloc(DIM,DIM);
		c2.copyUc(Uc2);
		x = gsl_matrix_get(Uc2,XX,ZZ); 
		y = gsl_matrix_get(Uc2,YY,ZZ); 
		z = gsl_matrix_get(Uc2,ZZ,ZZ);

		gvec *tc1	= gsl_vector_calloc(DIM);
		c1.copytc(tc1);
		gvec *tc2	= gsl_vector_calloc(DIM);
		c2.copytc(tc2);
		gsl_vector_sub(tc2,tc1);

		ftyp dx = gsl_vector_get(tc2,XX);
		ftyp dy = gsl_vector_get(tc2,YY);
		ftyp dz = gsl_vector_get(tc2,ZZ);

		ftyp dlen0=1.0/sqrt(dx*dx+dy*dy+dz*dz);
		ftyp ilen1=1.0/sqrt(a*a+b*b+c*c);
		ftyp ilen2=1.0/sqrt(x*x+y*y+z*z);

		ftyp angle=atan2( (dx*(b*z-c*y)+dy*(c*x-a*z)+dz*(a*y-b*x))*dlen0*ilen1*ilen2 , (a*x+b*y+c*z)*ilen1*ilen2 );

		gsl_matrix_free(Uc1); gsl_matrix_free(Uc2);
		std::cout << "INFO::TORSION  " << angle*180/M_PI << std::endl;
		return angle;
	}else{
		return -1000.0;
	}
}


std::pair<ftyp,ftyp >			
cluster_node0::angle_between(cluster0 c1, cluster0 c2, int I, int J){

	std::pair<ftyp,ftyp > ang_dist;

	if( c1.has_shape() && c2.has_shape() ){
		ftyp a,b,c,x,y,z;
		gmat *Uc1 = gsl_matrix_calloc(DIM,DIM);
		gmat *Uc2 = gsl_matrix_calloc(DIM,DIM);
		gvec *tc1 = gsl_vector_calloc(DIM);
		gvec *tc2 = gsl_vector_calloc(DIM);
		c1.copyUc(Uc1);
		c2.copyUc(Uc2);
		c1.copytc(tc1);
		c2.copytc(tc2);
		gsl_vector_sub(tc1,tc2);
		{ 
			a = gsl_matrix_get(Uc1,XX,ZZ); 
			b = gsl_matrix_get(Uc1,YY,ZZ); 
			c = gsl_matrix_get(Uc1,ZZ,ZZ);

			x = gsl_matrix_get(Uc2,XX,ZZ); 
			y = gsl_matrix_get(Uc2,YY,ZZ); 
			z = gsl_matrix_get(Uc2,ZZ,ZZ);

			ftyp dx = gsl_vector_get(tc1,XX);
			ftyp dy = gsl_vector_get(tc1,YY);
			ftyp dz = gsl_vector_get(tc1,ZZ);

			ftyp dlen0=1.0/sqrt(dx*dx+dy*dy+dz*dz);
			ftyp ilen1=1.0/sqrt(a*a+b*b+c*c);
			ftyp ilen2=1.0/sqrt(x*x+y*y+z*z);

			ftyp angle	= atan2( (dx*(b*z-c*y)+dy*(c*x-a*z)+dz*(a*y-b*x))*dlen0*ilen1*ilen2 , (a*x+b*y+c*z)*ilen1*ilen2 );
			ang_dist.first	= angle*180.0/M_PI;
			ang_dist.second	= 1.0/dlen0;
		}
		gsl_vector_free(tc1); gsl_vector_free(tc2);
		gsl_matrix_free(Uc1); gsl_matrix_free(Uc2);
	}else{
		ang_dist.first=0.0; ang_dist.second=-1.0;
	}
	return ang_dist;
}

std::vector<int>
cluster_node0::global_fragment_order( void )
{
	std::vector<int> at_order, tmp_order;

	if( cidx_.size() == parents_.second.length_M() )
	{
		gsl_vector *gv = gsl_vector_calloc( parents_.second.length_M() );
		parents_.second.copyv(gv);

		int N=gv->size;
		for(int i=0;i<N;i++) {
			tmp_order.push_back( (int)gsl_vector_get(gv,i) );
		}

		for(int i=0;i<N;i++){
			int c = tmp_order[i];
			for(int j=0;j<N;j++) {
				int cj = ( cidx_[j]>=0 )?(iidx_[cidx_[j]]):(-1);
				if( c == cj ) {
					at_order.push_back(j);
					cidx_[j]=(cidx_[j]+1)*(-1);
					break;
				}
			}
		}
		for(int i=0;i<N;i++) {
			cidx_[i] = (-1)*(cidx_[i])-1;
		}
	}

	return at_order;
}

std::pair<particles, std::vector<int> >
cluster_node0::apply_fragment_trans( particles coords , std::vector<int> p_ndx)
{
	std::pair<particles, std::vector<int> > p_id;
	if( idx_.size() == vc1_.size() && vc2_.size() == vc1_.size() && p_ndx.size()==coords.size() ) {
		particles swap_trans;
		std::vector<int> cidx;
		for(int i=0;i<idx_.size();i++) {

			richanalysis::coord_format cf;
			particles p_i = cf.par2par( coords, p_ndx, idx_[i] );

			if( p_i.size() != vc2_[i].length_M() )
				std::cout << "INFO::FORMAT::ERROR " <<  p_i.size() << " AND " << vc2_[i].length_M() << std::endl;

			gmat *P = gsl_matrix_calloc( DIM,	vc1_[idx_[i]].length_M()	);
			gmat *Q = gsl_matrix_calloc( DIM,	vc2_[i].length_M()		);

			vc1_[idx_[i]].copyM(P); // DATA  IN P 
			vc2_[i].copyM(Q);	// MODEL IN Q 

			gmat *U = gsl_matrix_calloc(DIM,DIM);
			gvec *t = gsl_vector_calloc(DIM);
			ftyp rmsd = shape_fit( P, Q, U, t, 1 );

			gsl_matrix_memcpy (U_, U);
			gsl_vector_memcpy (t_, t);
			invert_transform();
			particles p_tmp = apply_rot_trans( p_i , 1 );

			p_tmp.swap(p_i);
			p_tmp.clear();

			bUtSet_=1;
			sgn_=1;

			gsl_matrix_free(P);
			gsl_matrix_free(Q);
			gsl_matrix_free(U); 
			gsl_vector_free(t);
			for(int j=0;j<p_i.size();j++){
				swap_trans.push_back(p_i[j]);
				cidx.push_back( idx_[i] );
				cidx_.push_back(idx_[i]);
			}
		}
		cidx.swap(p_ndx);
		cidx.clear();
		p_id.first=swap_trans;
	}else{
		p_id.first=coords;
	}
	p_id.second = p_ndx;

	return p_id;
}

particles 
cluster_node0::assign_particles(void){
	particles ps;
	if(idx_.size() == vc2_.size()  ) {
		for(int i=0;i<vc2_.size();i++){
			int N = vc2_[i].length_M();
			gsl_matrix *M=gsl_matrix_calloc(DIM,N);
			gsl_vector *gv =gsl_vector_calloc(DIM);
			vc2_[i].copyM(M);
			for(int j=0;j<N;j++){
				particle p;
				p.second = gsl_vector_calloc(DIM);
				gsl_matrix_get_col (gv, M, j);
				gsl_vector_memcpy(p.second, gv);
				ps.push_back(p);
			}
			gsl_matrix_free(M);
			gsl_vector_free(gv);
		}
	}
	return ps;
}

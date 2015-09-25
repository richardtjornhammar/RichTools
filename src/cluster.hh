/*	cluster.hh	*/
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
#include "simple.hh"
#include "iofunc.hh"
#include "richfit.hh"

namespace richanalysis {

	class clustering : public simple_ops {
		public:
		//! null constructor
			inline clustering() {} //!< Null constructor

		//! performs k-means clustering on the specified data
			int gsl_kmeans( gmat *, gvec *, gmat *, gvec * );
			int gsl_kmeans( gmat *, gvec *, gmat *, gvec *, ftyp );
	};

	class cluster : public clustering, public tensorIO {
		public:
			cluster() { bSet_[0]=0; bSet_[1]=0; bSet_[2]=0; bSet_[3]=0; };
			std::vector<int>	get_labels  ( void );
			void			alloc_space ( int, int );
			void			seM ( int, int, ftyp );
			void			sev ( int, ftyp );
			void			seC ( int, int, ftyp );
			void			sew ( int, ftyp );
			int 			set_matrix( particles );  
			void			print_all( void ) {
							output_matrix( M_ );
							output_vector(vc_ );
							output_matrix( C_ );
							output_vector(wc_ );
						} ;
			int			perform_clustering( void );
			int			perform_clustering( ftyp );
			int			length_C(void){ return C_->size2; };
			int			length_M(void){ return M_->size2; };
			void	copyC(gmat *C0){ if(C0->size1==C_->size1&&C0->size2==C_->size2) { gsl_matrix_memcpy(C0, C_); } };
			void	copyM(gmat *M0){ if(M0->size1==M_->size1&&M0->size2==M_->size2) { gsl_matrix_memcpy(M0, M_); } };
			void	copyv(gvec *v0){ if(v0->size ==vc_->size) { gsl_vector_memcpy(v0, vc_); } };
			void	copyw(gvec *w0){ if(w0->size ==wc_->size) { gsl_vector_memcpy(w0, wc_); } };
			int	NpC(int i){ if(NperC_.size()>0&&i<NperC_.size()){ return NperC_[i]; } };
			~cluster() {}
	private:
		gmat	*M_; 		// 0 THE ORIGINAL COORDINATES 
		gmat	*C_; 		// 1 THE CENTROIDS
		gvec	*vc_;		// 0 UNSORTED LABELS
		gvec	*wc_;		// 1 UNSORTED LABELS
		int	bSet_[4];
		std::vector<int>	NperC_;
	};

	class node_indices : public fitting {
		public:
			node_indices() { bDirRel_=0; bUtSet_=0; 
					U_  = gsl_matrix_calloc( DIM, DIM );  t_ = gsl_vector_calloc( DIM ); 
					iU_ = gsl_matrix_calloc( DIM, DIM ); it_ = gsl_vector_calloc( DIM ); sgn_=1; }

			ftyp			find_centroid_relation( cluster c1, cluster c2 );
			ftyp			find_shape_trans( cluster c1, cluster c2 );

			std::vector<int> 	get_indices(void){ return idx_; };
			int			direction_relation(void){ return bDirRel_;};
			int			have_transform(void){ return bUtSet_; };
			particles 		apply_rot_trans( particles , int );
			particles 		apply_rot_trans( particles );
			void			invert_transform(void);
			void			printUt(void){  output_matrix(U_); output_vector(t_); };
			void			printiUt(void){ output_matrix(iU_);output_vector(it_);};
		private:
			int bDirRel_;
			std::vector<int> idx_;
			int N_,M_,I_,J_;
			int bUtSet_;
			int sgn_;
			gmat *U_,*iU_;
			gvec *t_,*it_;
	};

}

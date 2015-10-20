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

	class clustering : public simple_ops , public tensorIO  {
		public:
		//! null constructor
			inline clustering() {} //!< Null constructor

		//! performs k-means clustering on the specified data
			int gsl_kmeans( gmat *, gvec *, gmat *, gvec * );

		//! performs my seeded k-means clustering on the specified data
			int gsl_seeded_kmeans( gmat *, gvec *, gmat *, gvec * );

		//! performs my own connectivity clustering algorithm
			int connectivity( gmat * , ftyp val );
	};

	class cluster :  public clustering {
		public:
			cluster() { 	bSet_ = false; rDIM_=DIM; };
			void		alloc_space ( int, int );
			int		set_matrix( particles ); 
			int		find_centroids( void );
			std::vector<int>	get_labels( void );
			std::vector<int>	get_clabels( void );

			void		copyC(gmat *C0) { 
						if(C0->size1==C_->size1&&C0->size2==C_->size2) { gsl_matrix_memcpy(C0, C_); } 
					};
			void		copyM(gmat *M0) { 
						if(M0->size1==M_->size1&&M0->size2==M_->size2) { gsl_matrix_memcpy(M0, M_); } 
					};
			void		copyv(gvec *v0) { 
						if(v0->size ==vc_->size) { gsl_vector_memcpy(v0, vc_); } 
					};
			void		copyw(gvec *w0) { 
						if(w0->size ==wc_->size) { gsl_vector_memcpy(w0, wc_); } 
					};

			void		writev(gvec *v0) { 
						if(v0->size ==vc_->size) { gsl_vector_memcpy(vc_, v0); } 
					};

			int		get_cDIM(void)		{ return rDIM_; 	};	//	centroid dimension
			void		set_cDIM(int rDIM)	{ rDIM_=rDIM; 
						if(bSet_) { 
							gsl_matrix_free(C_); 
							C_=gsl_matrix_calloc(DIM,rDIM); 
							gsl_vector_free(wc_); 
							gsl_vector_calloc(rDIM); }		};	//	centroid dimension
			int		length_C(void)		{ return C_->size2;	};
			int		length_M(void)		{ return M_->size2;	};
			bool		isSet(void)		{ return ( bSet_ );	};	
			ids		getIDs(void)		{ return (idlabels_);	};
			void		setIDs(ids vid)		{ idlabels_ = vid ;	};
			void		print_all( void ) {
							output_matrix( M_ );
							output_vector(vc_ );
							output_matrix( C_ );
							output_vector(wc_ );
						} ;
		private:
			gmat	*M_; 		// 0 THE ORIGINAL COORDINATES 
			gmat	*C_; 		// 1 THE CENTROIDS
			gvec	*vc_;		// 0 UNSORTED LABELS
			gvec	*wc_;		// 1 UNSORTED LABELS
			int	bSet_;
			ids	idlabels_;
			std::vector< int >	NperC_;
			int	rDIM_;
	};
	
	typedef std::pair< cluster, cluster > node;
	typedef std::vector< node > layer;

	class node_analysis : public fitting, public fileIO  {
		public:
			node_analysis() { 		bNode_ = false; bLayer_ = false; 	};
			node_analysis( node n ) { 	bNode_ = assign_node( n ); 		};
			std::vector< int >	find_centroid_relation( void );
			particles		regular_fit(void);
			particles		seeded_centroids( void );
			std::vector< int >	global_index_order ( void ) { return glob_idx_order_; };
			std::vector< int >	cluster_index_order( void ) { return cidx_; };
			std::pair<int,int>	size(void){ std::pair<int,int> pi; 
							pi.first=parents_.first.length_M();
							pi.second=parents_.second.length_M();
							return pi; };
			void			copyv(gvec *v0, int sw) { 
							if(sw==1) { parents_.first .copyv(v0); }
							if(sw==2) { parents_.second.copyv(v0); } 
						};
			void			copyw(gvec *w0, int sw) { 
							if(sw==1) { parents_.first .copyw(w0); }
							if(sw==2) { parents_.second.copyw(w0); } 
						};
			void			centroid_to_nearest(void);
			bool			assign_node( node n );
			bool			allSet()	{ return ( bNode_ && bLayer_ ); };
			bool			haveNode()	{ return ( bNode_); };
			bool			haveLayer()	{ return ( bLayer_ ); };
			layer			get_node_layer( void ){ return children_; };	
			void			print_all(void){ parents_.first.print_all();parents_.second.print_all();};		
			~node_analysis(){};
		private:
			node	parents_;
			layer	children_;
			bool bNode_; bool bLayer_;
			std::vector<int>  idx_;
			std::vector<int> iidx_;
			std::vector<int> cidx_;
			std::vector<int> glob_idx_order_;
	};

	class layer_analysis : public fileIO {
		public:
			layer_analysis(){ bSet_ = false; };
			layer_analysis(layer l){ bSet_ = true; own_=l; };
			void output_layer( std::string );
		private:
			layer	own_;
			bool	bSet_;
	};

	class particle_analysis : public clustering , public fileIO {
		public:
			particle_analysis(){ bAssigned_ = false; bMatrices_ = false; };
			int 			calc_distance_matrices(void);
			std::vector< int >	find_via_distance	( gmat *A, ftyp level );
			std::vector<int> 	outp_distance_matrix	( gmat *A, ftyp level );
			std::vector<int> 	outp_distance_matrix	( ftyp level ) { 
						return outp_distance_matrix( A_, level ); };
			std::vector< std::pair<ftyp, std::pair< int, int > > >	compare_dist_matrices(gmat *A, gmat *B, ftyp val);
			std::vector< std::pair<ftyp, std::pair< int, int > > >	compare_dist_matrices(ftyp val) { 
						 return compare_dist_matrices(A_,B_,val); };
			void			assign_particles( particles pd, particles pm )	{
				if( pd.size() > pm.size() ) {
					density_.clear(); density_.insert( density_.end() , pd.begin() , pd.end() );
					model_.clear(); model_.insert( model_.end() , pm.begin() , pm.end() );
				} else {
					density_.clear(); density_.insert( density_.end() , pm.begin() , pm.end() );
					model_.clear(); model_.insert( model_.end() , pd.begin() , pd.end() );					
				} 
				bAssigned_ = true; 
				assign_matrices(); };
			void			assign_matrices( void );
			int	check(){ return M_->size2; };
			void	output_result( std::string );
			bool			complete( void ) { return bAssigned_; };
			bool			matrices( void ) { return bMatrices_; };
			particles		output_reduced_density(void);
		private:
			bool bAssigned_, bMatrices_;
			particles density_;
			particles model_;
			ftyp	rcut_;
			gmat	*A_, *B_;
			gmat	*C_, *M_;			
			gvec	*vc_;		// 0 UNSORTED LABELS
			gvec	*wc_;		// 1 UNSORTED LABELS
			std::vector< std::pair<ftyp, std::pair< int, int > > > relation_;
	};

	class cluster0 : public clustering, public linalg {
		public:
			cluster0() { 	bSet_[0]=0; bSet_[1]=0; bSet_[2]=0; bSet_[3]=0; bPCset_=0;
					Uc_  = gsl_matrix_calloc( DIM, DIM );	tc_ = gsl_vector_calloc( DIM ); 
					Axis_ = gsl_vector_calloc( DIM );	Euler_ = gsl_vector_calloc( DIM ); };

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
			int			length_C(void){ return C_->size2; };
			int			length_M(void){ return M_->size2; };
			void	copyC(gmat *C0)		{ if(C0->size1==C_->size1&&C0->size2==C_->size2) { gsl_matrix_memcpy(C0, C_); } };
			void	copyM(gmat *M0)		{ if(M0->size1==M_->size1&&M0->size2==M_->size2) { gsl_matrix_memcpy(M0, M_); } };
			void	setM( gmat *M0)		{ if(M0->size1==M_->size1&&M0->size2==M_->size2) { gsl_matrix_memcpy(M_, M0); } };
			void	copyUc(gmat *U0)	{ if(U0->size1==Uc_->size1&&U0->size2==Uc_->size2) { gsl_matrix_memcpy(U0, Uc_); } };
			void	copytc(gvec *tc0)	{ if(tc0->size ==tc_->size) { gsl_vector_memcpy(tc0, tc_); } };
			void	copyv(gvec *v0)		{ if(v0->size ==vc_->size) { gsl_vector_memcpy(v0, vc_); } };
			void	copyw(gvec *w0)		{ if(w0->size ==wc_->size) { gsl_vector_memcpy(w0, wc_); } };
			int	NpC(int i)		{ if(NperC_.size()>0&&i<NperC_.size()){ return NperC_[i]; } };
			bool			isSet()	{return bSet_[0]>0&&bSet_[1]>0&&bSet_[2]>0&&bSet_[3]>0;}
			ftyp			find_shape();
		//! Calculates the shape
			ftyp			shape(	gmat *P , gmat *U , gvec *t 	);
			void			print_shape( void ) { output_matrix( Uc_ ); };
			bool			has_shape(void){return (bPCset_==1); };
			particle		normal(void);
			particle		center(void);
			~cluster0() {}
	private:
		gmat	*M_; 		// 0 THE ORIGINAL COORDINATES 
		gmat	*C_; 		// 1 THE CENTROIDS
		gvec	*vc_;		// 0 UNSORTED LABELS
		gvec	*wc_;		// 1 UNSORTED LABELS
		int	bSet_[4];
		std::vector<int>	NperC_;
//	shape stuff
		gmat *Uc_;
		gvec *tc_;
		gvec *Euler_;
		gvec *Axis_;	
		ftyp angle_;
		int bPCset_;
	};

	class cluster_node0 : public fitting {
		public:
			cluster_node0() { bDirRel_=0; bUtSet_=0; subSet_=0;
					U_  = gsl_matrix_calloc( DIM, DIM );  t_ = gsl_vector_calloc( DIM ); 
					iU_ = gsl_matrix_calloc( DIM, DIM ); it_ = gsl_vector_calloc( DIM ); sgn_=1; }

			void			assign_sub( cluster0 c1, cluster0 c2 );
			bool			subSet() { return (subSet_==1); };

			std::vector<int>	find_centroid_relation( void );

			std::vector<int>	find_simple_relation( void );

			std::pair< particles,
			std::vector<int> >	apply_fragment_trans( particles, std::vector<int> );

			std::vector<int> 	get_indices(void){ return idx_; };
			std::vector<int> 	global_fragment_order(void);
			int			direction_relation(void){ return bDirRel_;};
			int			have_transform(void){ return bUtSet_; };
			particles 		apply_rot_trans( particles , int );
			particles 		apply_rot_trans( particles );
			particles		assign_particles( void ); // packs subclusters into a particles struct
			void			invert_transform(void);
			void			printUt(void){  output_matrix(U_); output_vector(t_); };
			void			printiUt(void){ output_matrix(iU_);output_vector(it_);};
			ftyp			angle_between(cluster0 c1, cluster0 c2);
			std::pair<ftyp, ftyp >	angle_between(cluster0 c1, cluster0 c2, int i, int j);
		private:
			int subSet_;
			std::pair<richanalysis::cluster0, richanalysis::cluster0 > parents_;
			std::vector<richanalysis::cluster0 > vc1_, vc2_; // subclusters
			int bDirRel_;
			std::vector<int> idx_;		// clusters
			std::vector<int> iidx_;		// clusters
			std::vector<int> cidx_;
			int N_,M_,I_,J_;
			int bUtSet_;
			int sgn_;
			gmat *U_,*iU_;
			gvec *t_,*it_;
	};

}

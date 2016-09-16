function[U, r, lrms] = unordered_align( P, Q )

	sz1	= size(P) ;
	sz2 	= size(Q) ;

	if (length(sz1) ~= 2 || length(sz2) ~= 2)
		error 'P and Q must be matrices' ;
	end

	D	= sz1(1)           % dimension of space
	N	= sz1(2)           % number of points
	M	= sz2(2)

	n 	= ones(1,M)/M ;
	m 	= ones(1,N)/N ;

	p0	= P*m' ;		% the centroid of P
	q0	= Q*n' ;		% the centroid of Q

	v1	= ones(1,N);		% row vector of N ones
	v2	= ones(1,M);

	P 	= P - p0*v1 ;		% translating P to center the origin
	Q 	= Q - q0*v2 ;		% translating Q to center the origin

	Pdm	= bsxfun(@times,m,P)
	Qdn 	= bsxfun(@times,n,Q)
 
	sP	= (Pdm*Pdm')		;
	sQ	= (Qdn*Qdn')		;
	C 	= sP*sQ' 		; 
	size_C	= size(C)		;
	[V,S,W] =  svd(C)		;	
	size_V	= size(V) 		;
	size_S	= size(S) 		;
	size_W	= size(W) 		;
	I 	=  eye(D)		;

	if (det(V*W') < 0)  	% more numerically stable than using (det(C) < 0)
		I(D,D) = -1 ;
	end
	rvmin	= 1e8;
	imin 	= 1;
	niter	= 1;
	if( D==3 )
		niter = 8
	end
	U 	= size( C );
	r 	= size(1,D);
	i_deadp = [ 1 ] ;
	i_deadq = [ 40 , 41 , 42 ] ;

	for i=1:niter
		if(i==1)
	                I(1,1)= 1; I(2,2)= 1; I(3,3)= 1;
		end
                if(i==2)
        	        I(1,1)= 1; I(2,2)= 1; I(3,3)=-1;
		end
                if(i==3)
                	I(1,1)= 1; I(2,2)=-1; I(3,3)= 1;
		end
                if(i==4)
			I(1,1)=-1; I(2,2)= 1; I(3,3)= 1;
		end
                if(i==5)
  	              	I(1,1)= 1; I(2,2)=-1; I(3,3)=-1;
		end
                if(i==6)
                	I(1,1)=-1; I(2,2)= 1; I(3,3)=-1;
		end
                if(i==7)
                	I(1,1)=-1; I(2,2)=-1; I(3,3)= 1;
		end
                if( i==8 )
                	I(1,1)=-1; I(2,2)=-1; I(3,3)=-1;
		end
		Ut  	= W*I*V'	;
		rt  	= q0 - Ut*p0 	;
		Diff	= Ut*sP - sQ	;     		% P, Q already centered
		i
		lrms	= sqrt(sum(sum(Diff.*Diff))/N);  	% (for the non-weighted case)
		up	= Ut*P; nq = Q;
		tot_sum = 0;
		if( true )
			ndead	= length(i_deadp) ;
			mdead	= length(i_deadq) ;
			for iu  = 1:ndead
				dmin = 1e8;
				for iq = 1:mdead
					dr  = up (:,i_deadp(iu) ) - nq (:,i_deadq(iq) );
					len = sqrt( sum( dr.*dr ) );
					if len < dmin
						dmin=len;
					end
				end
				tot_sum += dmin;
			end
		end

		if( true )
			for iu = 1:N
				dmin   = 1e8;
				for iq = 1:M
					dr  = up(:,iu)-nq(:,iq);
					len = sqrt( sum(dr.*dr) );
					if len < dmin
						dmin = len;
					end
				end
				tot_sum += dmin;
			end
		end

		lrms = tot_sum/N	;
		rv   = Ut*P+rt		;

		if( lrms < rvmin )
			i
			U	= Ut		
			r 	= rt	
			rvmin	= lrms	;
		end

	end
end

clear all
load distg.dat
H=distg;
N=length(H)
BFACTOR=20
bf=sqrt(BFACTOR/10.0/2.0)
bf=0.0
d2ij=zeros(N,N);
for i=1:N
	for j=1:N
		d = sum((H(i,:)-H(j,:)).^2);
		if (j>i)
			d+=bf;
		end
		if (j<i)
			d-=bf;
		end
		d2ij(i,j)=d;
	end
end

D=zeros(N,N);

for i=1:N
	for j=1:N
		D(i,j) = (d2ij(i,N)+d2ij(j,N)-d2ij(i,j))*0.5;
	end
end

[ U , S , V ]	= svd(D);
E		= S^0.5;
e		= E(:,1:3);
X		= U*e;

labels		= ones(N,1);
labels		= [ 'C' 'H' 'H' 'H' 'C' 'H' 'H' 'C' 'H' 'H' 'C' 'H' 'H' 'C' 'H' 'H' 'C' 'H' 'H' 'H' ]';
test		= [ labels ones(N,1)*"  " num2str(X(:,1:3)) ]
save test.dat test


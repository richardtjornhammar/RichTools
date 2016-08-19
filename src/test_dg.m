clear all
H=[
  -23.7100   24.1000   85.4400
  -22.5600   23.7600   85.6500
  -21.5500   24.6200   85.3800
  -22.2600   22.4200   86.1900
  -23.2900   21.5300   86.4800
  -20.9300   22.0300   86.4300
  -20.7100   20.7600   86.9400
  -21.7900   19.9300   87.1900
  -23.0300   20.3300   86.9600
  -24.1300   19.4200   87.2500
  -23.7400   18.0500   87.0000
  -24.4900   19.4600   88.7500
  -23.3700   19.8900   89.5200
  -24.8500   18.0000   89.0900
  -23.9600   17.4800   90.0800
  -24.6600   17.2400   87.7500
  -24.0800   15.8500   88.0100
  -23.9600   15.1600   86.7600
  -23.3400   13.7100   87.1000
  -21.9600   13.8700   87.6300
  -24.1800   13.0300   88.1100
  -23.2900   12.8200   85.7600
  -23.1900   11.2800   86.2200
  -21.8100   11.0000   86.7000
  -24.1500   11.0300   87.3200
  -23.5300   10.3200   84.9800
  -23.5400    8.9800   85.4800
  -23.8600    8.0100   84.3400
  -23.9800    6.5760   84.8900
  -23.2800    6.4460   86.1300
  -23.3000    5.7330   83.7800
  -22.7300    4.5360   84.3100
  -22.2000    6.7130   83.3000
  -22.7900    8.0170   83.3800
  -21.8100    6.4120   81.9200
  -20.8500    5.5220   81.5200
  -20.8300    5.5670   80.1200
  -21.7700    6.4720   79.7400
  -22.3400    6.9680   80.8000
  -20.0100    4.6970   82.1500
  -19.1800    3.9390   81.4700
];

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

dim=3;
[ U , S , V ]	= svd(D);
E		= S^0.5;
e		= E(:,1:exp_dim);
X		= U*e;

[ U, S ,V ]	= svd(U,1);
Z=U*S(:,1:3);
Y 		= U'*U*S*V'

labels		= ones(N,1)*10;
%test		= [ labels ones(N,1)*"  " num2str(X(:,1:exp_dim)) ]
%save test.dat test


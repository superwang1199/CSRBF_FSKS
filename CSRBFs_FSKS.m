clear all
tic
%RBF
rbf = @(e,r) r.^8.*(66*spones(r)-154*r+121*r.^2-32*r.^3);
dxrbf =  @(e,r,dx) -22*e^2.*dx.*r.^7.*(16*r.^2-39*r+24*spones(r));
dyrbf =  @(e,r,dy) -22*e^2.*dy.*r.^7.*(16*r.^2-39*r+24*spones(r));
dzrbf =  @(e,r,dz) -22*e^2.*dz.*r.^7.*(16*r.^2-39*r+24*spones(r));
dxxrbf = @(e,r,dx) 528*e^4*dx.^2.*r.^6.*(7*spones(r)-6*r)...
                   -22*e^2*r.^7.*(24*spones(r)-39*r+16*r.^2);
dyyrbf = @(e,r,dy) 528*e^4*dy.^2.*r.^6.*(7*spones(r)-6*r)...
                   -22*e^2*r.^7.*(24*spones(r)-39*r+16*r.^2);
dzzrbf=@(e,r,dz) 528*e^4*dz.^2.*r.^6.*(7*spones(r)-6*r)...
                   -22*e^2*r.^7.*(24*spones(r)-39*r+16*r.^2);
Te = @(x,y,z,t) (0.5+exp(-t)).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z);
Ti = @(x,y,z,t)  exp(-t).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z);
Tr = @(x,y,z,t)  exp(-t).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z);
%% 
F1 = @(x,y,z,t) -exp(-t).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z)-pi.*pi.*(0.5+exp(-t)).^2.*(cos(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z).^2 ...
                +sin(pi.*x).^2.*cos(pi.*y).^2.*sin(pi.*z).^2+sin(pi.*x).^2.*sin(pi.*y).^2.*cos(pi.*z).^2)+3.*pi.*pi.*(0.5+exp(-t)).^2 ...
                .*sin(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z).^2+0.5.*sin(pi.*x).*sin(pi.*y).*sin(pi.*z)+sin(pi.*x).^4.*sin(pi.*y).^4.* ...
                sin(pi.*z).^4.*((0.5+exp(-t)).^4-exp(-4.*t));
            
F2 = @(x,y,z,t) -(0.5+exp(-t)).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z)-pi.*pi.*exp(-2.*t).*(cos(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z).^2 ...
                +sin(pi.*x).^2.*cos(pi.*y).^2.*sin(pi.*z).^2+sin(pi.*x).^2.*sin(pi.*y).^2.*cos(pi.*z).^2)+3.*pi.*pi.*exp(-2.*t) ...
                .*sin(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z).^2;

F3 = @(x,y,z,t) (-3.*exp(-4.*t)-(0.5+exp(-t)).^4).*sin(pi.*x).^4.*sin(pi.*y).^4.*sin(pi.*z).^4-16.*pi.*pi.*exp(-5.*t).*sin(pi.*x) ...
                .^3.*sin(pi.*y).^3.*sin(pi.*z).^3.*(cos(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z).^2 +sin(pi.*x).^2.*cos(pi.*y).^2.*sin(pi.*z).^2 ...
                +sin(pi.*x).^2.*sin(pi.*y).^2.*cos(pi.*z).^2)+12.*pi.*pi.*exp(-5.*t).*sin(pi.*x).^5.*sin(pi.*y).^5.*sin(pi.*z).^5;
%% 
rho=1; a_a=1; w_ei=1; ck=1; c_ve=1; c_vi=1;
K=4;
n=10;
dt = 1e-4;
sn=2^K+1;
ep = [0.6 0.6 0.6];
N1=[0 27 125 729 4913 ];  
N5=[27 125 729 4913 8000];
%% 
name = sprintf('yuanzhu_%d',N5(K));load(name);
xdata=dsites1(:,1)/2; ydata=dsites1(:,2)/2; zdata=dsites1(:,3)/2; 
intdata = [xdata,ydata,zdata];  %ÄÚµã
if sn==3
    sn=2;
end
bdydata=[];
t1=linspace(-1,1,sn+2)';
t2=linspace(0,2*pi,4*(sn+1)+1)';
for j=2:length(t1)-1
    for i=1:length(t2)-1  
        x1=cos(t2(i));
        y1=sin(t2(i));
        z1=t1(j);
        bdydata=[bdydata;x1 y1 z1];
    end
end
t3=linspace(0,1,(sn+1)/2+1)';
bdydata=[bdydata;0 0 1;0 0 -1];
s=0;
for j=2:length(t3)
    s=s+2;
    t4=linspace(0,2*pi,(4*s)+1)';
    for i=1:4*s
        x2=t3(j)*cos(t4(i));
        y2=t3(j)*sin(t4(i));
        z2=1; z3=-1;
        bdydata=[bdydata;x2 y2 z2;x2 y2 z3];
    end
end
bdydata=bdydata/2;
xbdydata=bdydata(:,1);  ybdydata=bdydata(:,2);  zbdydata=bdydata(:,3);
alldata=[intdata;bdydata];
epoints = alldata;
name = sprintf('eg1_epoints');load(name);
%%
for j=1:K
    xctrs{j}=[xdata(N1(j)+1:N1(j+1)) ydata(N1(j)+1:N1(j+1)) zdata(N1(j)+1:N1(j+1))];
    DM = DistanceMatrixCSRBF(epoints,xctrs{j},ep(j));
    EM{j} = rbf(ep(j),DM);  
    a = DistanceMatrixCSRBF(intdata,xctrs{j},ep(j));
    aa = DistanceMatrixCSRBF(alldata,xctrs{j},ep(j));
    A{j}=rbf(ep(j),a);    %A
    AA{j}=rbf(ep(j),aa);
    m = cell2mat(xctrs(1,j));
    ax = DifferenceMatrix(intdata(:,1),m(:,1));
    AX{j}=dxrbf(ep(j),a,ax);     %AX
    AXX{j}=dxxrbf(ep(j),a,ax);     %AXX
    ay=DifferenceMatrix(intdata(:,2),m(:,2));
    AY{j}=dyrbf(ep(j),a,ay);     %AY
    AYY{j}=dyyrbf(ep(j),a,ay);     %AYY
    az=DifferenceMatrix(intdata(:,3),m(:,3));
    AZ{j}=dzrbf(ep(j),a,az);     %AZ
    AZZ{j}=dzzrbf(ep(j),a,az);     %AZZ
    bdy = DistanceMatrixCSRBF(bdydata,xctrs{j},ep(j));
    BDY{j} = rbf(ep(j),bdy);    
    if(j==1)
        nEM=EM{1};
        nA=A{1};
        nAA=AA{1};
        nAX=AX{1};
        nAY=AY{1};
        nAZ=AZ{1};
        nAXX=AXX{1};
        nAYY=AYY{1};
        nAZZ=AZZ{1};
        nBDY=BDY{1};
    end
    if(j>1)
        nEM=[nEM EM{j}];
        nA=[nA A{j}];
        nAA=[nAA AA{j}];
        nAX=[nAX AX{j}];
        nAY=[nAY AY{j}];
        nAZ=[nAZ AZ{j}];
        nAXX=[nAXX AXX{j}];
        nAYY=[nAYY AYY{j}];
        nAZZ=[nAZZ AZZ{j}];
        nBDY=[nBDY BDY{j}];
    end
end
Te0=Te(alldata(:,1),alldata(:,2),alldata(:,3),0);
Ti0=Ti(alldata(:,1),alldata(:,2),alldata(:,3),0);
Tr0=Tr(alldata(:,1),alldata(:,2),alldata(:,3),0);
u=nAA\Te0;  
v=nAA\Ti0;  
w=nAA\Tr0;  
clear nAA;  clear Te0;  clear Ti0;  clear Tr0;
for k=1:n
    Lu = dt.*((nAX*u).*nAX+(nAY*u).*nAY+(nAZ*u).*nAZ);
    Lv = dt.*((nAX*v).*nAX+(nAY*v).*nAY+(nAZ*v).*nAZ);
    Lw = 16*dt.*(nA*w).^3.*((nAX*w).*nAX+(nAY*w).*nAY+(nAZ*w).*nAZ);
    Hu = dt.*(nA*u).*(nAXX+nAYY+nAZZ);
    Hv = dt.*(nA*v).*(nAXX+nAYY+nAZZ);
    Hw = 4*dt.*(nA*w).^4.*(nAXX+nAYY+nAZZ);
    g1 = Te(bdydata(:,1),bdydata(:,2),bdydata(:,3),dt*k);
    g2 = Ti(bdydata(:,1),bdydata(:,2),bdydata(:,3),dt*k);
    g3 = Tr(bdydata(:,1),bdydata(:,2),bdydata(:,3),dt*k);
    b1 = dt*F1(intdata(:,1),intdata(:,2),intdata(:,3),dt*k);
    b2 = dt*F2(intdata(:,1),intdata(:,2),intdata(:,3),dt*k);
    b3 = dt*F3(intdata(:,1),intdata(:,2),intdata(:,3),dt*k);
    u = [rho*c_ve*nA-Lu-Hu+dt*rho*w_ei*nA;nBDY]\([rho*c_ve*nA-dt*a_a*ck*(nA*u).^3.*nA;zeros(size(nBDY))]*u ...
        +[dt*rho*w_ei*nA*v+dt*a_a*ck*(nA*w).^4+b1;g1]);
    v = [rho*c_vi*nA-Lv-Hv+dt*rho*w_ei*nA;nBDY]\([rho*c_vi*nA;zeros(size(nBDY))]*v+[dt*rho*w_ei*nA*u+b2;g2]);
    w =[4*(nA*w).^3.*nA-Lw-Hw;nBDY]\([4*(nA*w).^3.*nA-dt*ck*(nA*w).^3.*nA;zeros(size(nBDY))]*w+[dt*ck*(nA*u).^4+b3 ;g3]);   %((Tr)^n)^4
    uu(:,k)= u;  vv(:,k)= v;  ww(:,k)= w;
end
toc
pfTe = nEM*u;  
pfTi = nEM*v;
pfTr = nEM*w;
  function DM = DistanceMatrixCSRBF(dsites,ctrs,ep)
  DM=DistanceMatrix(dsites,ctrs);
  DM=1-ep*DM;
  DM(DM<0)=0;
  end
  function DM = DifferenceMatrix(datacoord,centercoord)
  [dr,cc] = ndgrid(datacoord(:),centercoord(:));
  DM = dr-cc;
  end
  

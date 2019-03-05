% biome kmeans
close all;
clear all;

% define scale factor
% SST, d(MLD), NO3, SiO3, PO4, Chl-a
scl = [1 1 1 1 1 1];

% coordinates
lon=-179.5:179.5;
lat=-89.5:89.5;
[LON,LAT]=meshgrid(lon,lat);

% loading basin masks
load globalmask.mat;

% load mixed layer
load mxl.mat;

% load log10-chla
load chla.mat;

% load NO3
load no3.mat;

%load SiO3
load sil.mat;

%load PO4
load po4.mat;

%load sst
load sst.mat;

% 2D field to 1D vector
cnt=0;
for i=1:360
    for j=1:180
        if ~isnan(sst(i,j)+mxl(i,j)+chla(i,j)+no3(i,j)+sil(i,j)+po4(i,j))
            cnt=cnt+1;
            xind(cnt)=i;
            yind(cnt)=j;
            sst0(cnt)=sst(i,j);
            basin(cnt)=mask(i,j);
            mlx0(cnt)=mxl(i,j);
            chl0(cnt)=chla(i,j);
            nit0(cnt)=no3(i,j);
            sil0(cnt)=sil(i,j);
            po40(cnt)=po4(i,j);
        end
    end
end

% X is the dimensional data matrix. Each column is different variable
K=find(basin==1);
X(:,1)=sst0(K);
X(:,2)=mlx0(K);
X(:,3)=nit0(K);
X(:,4)=sil0(K);
X(:,5)=po40(K);
X(:,6)=chl0(K);
X(isnan(X)|isinf(abs(X)))=0; % for safety

% Standardize the data matrix. All variables are normalized and scaled. 
for l=1:length(scl)
      mu(l)=mean(X(:,l));
      SD(l)=std(X(:,l));
      z(:,l)=(X(:,l)-mu(l))/SD(l)*scl(l);
end

% start the calculation. loop over for the number of K-clusters
for kval = 2:20;

    % cluster!
    disp(['working on kval = ',num2str(kval)])
    [IDX,C]=kmeans(z,kval,'MaxIter',200);

    % re-build the cluster matrix by placing the centeroid values
    for m=1:length(IDX)
       zc(m,:)=C(IDX(m),:);
    end

    % variance explained by the centeroid reconstruction
    n=kval-1;
    Varexp(n)=var(zc(:));
    Nc(n)=kval;

    % record the cluster as a map
    tmp=NaN*ones(360,180);
    for m=1:length(K)
        tmp(xind(K(m)),yind(K(m)))=IDX(m);
        Index{n}.Map=tmp;
    end

    % calculate the area covered by each cluster
    dA1=(6.37e6*pi/180)^2;
    dA=dA1*cos(LAT'*pi/180).*mask;
    for l=1:kval
       tmp=dA;
       tmp(Index{n}.Map~=l)=NaN;
       area(l,n)=nansum(nansum(tmp));
    end

    tmpc=NaN*zeros(length(scl),kval);
    % centroids
    tmpc = diag(SD./scl)*C'+repmat(mu',[1 kval]);
    %tmpc(1,:)=asin(tmpc(1,:))*180/pi;
    Index{n}.CenteroidDim=tmpc;
    Index{n}.CenteroidND=C';

end

save kmeans_stats.mat Varexp Nc Index area;


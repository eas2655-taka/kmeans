% biome kmeans
close all;
clear all;

scl = 1;

% coordinates
lon=-179.5:179.5;
lat=-89.5:89.5;
[LON,LAT]=meshgrid(lon,lat);

% loading basin masks
load globalmask.mat;

% load log10-chla
load chla.mat;

% 2D field to 1D vector
cnt=0;
for i=1:360
    for j=1:180
        if ~isnan(chla(i,j))
            cnt=cnt+1;
            xind(cnt)=i;
            yind(cnt)=j;
            basin(cnt)=mask(i,j);
            chl0(cnt)=chla(i,j);
        end
    end
end

% X is the dimensional
K=find(basin==1);
X(:,1)=chl0(K);
X(isnan(X)|isinf(abs(X)))=0; % for safety

% Standardize the data 
for l=1:length(scl)
      mu(l)=mean(X(:,l));
      SD(l)=std(X(:,l));
      z(:,l)=(X(:,l)-mu(l))/SD(l)*scl(l);
end

% start the calculation. loop over for the number of K-clusters
for kval = 2:10;

    % cluster!
    disp(['working on kval = ',num2str(kval)])
    %seedval=rand(kval,1)*2 - 1;
    %[IDX,C]=kmeans(z,kval,'MaxIter',200,'Start',seedval);
    [IDX,C]=kmeans(z,kval,'MaxIter',200,'Replicates',5,'Display','final');

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
    Index{n}.CentroidDim=tmpc;
    Index{n}.CentroidND=C';

end

save chla_stats.mat Varexp Nc Index area;


%% Name:
% "LatticeHopping_NaYbF4"
%% Author:
% Yunxiang Zhang, Fudan University, 200438, PR China
% zyx@fudan.edu.cn
% Updated versions can be found on https://github.com/YunxiangZhangLab/LatticeHopping
%% Purpose:
% MC simulation of energy migration over Yb3+ sublattice in hexagonal phase
% NaYbF4 nanocrystals.
%% Keywords:
% NaYbF4 upconversion; 
% energy migration; 
% Dexter energy transfer;
% exchange mechanism;
% hexagonal phase;
% core-shell nanoparticle;
% Monte Carlo simulation...
%% Description:
% This piece of code simulates the energy migration trajectories for a
% quantum of photon energy absorbed by each individual sensitizer Yb3+ ions
% followed by its hopping over the Yb3+ sublattice according to the
% exchange mechanism until it is terminated at the boundary of core-shell
% interface. Core/shell diameters may be changed to different values and
% number of grid points need to be adjusted accordingly to cover the
% core-shell structure. Core:Shell volume ratio can be tuned this way and
% verifed by checking the number of RE3+ sites. To simulate inward energy
% migration from shell to core please set "inward=1" or otherwise set
% "inward=0" for outward migration. Simulation results were recorded as
% HDF5 files in the same directory. 
%% Copyright (C) 2021~2022, Yunxiang Zhang
% This software is provided as is without any warranty whatsoever.
% Permission to use, copy, modify, and distribute modified or
% unmodified copies is granted, provided this copyright and disclaimer
% are included unchanged. 
%% work dir
% set up work directory
path='/Users/yxzhang/Documents/MATLAB/latticehopping/';
%% const
% set up lattice constants and core-shell structure size configuration.
% core and shell compartment were designed to have almost the same volume.
roc=5.93;% unit cell length
roa=3.47;% unit cell height
dcore=114.61; % core diameter
dshell=144.40; % shell diameter

%% coordinates init
% set up lattice grid index (x, y, h) where x and y were cube coordinates
% for indexing tiled hexagonals and h were coordinates along Z-direction.

%number of grid points
nx=29;
ny=29;
nh=89;
%
tc=[];
% set up (x,y,h) indexing range in integers
% x -14:1:14; y -14:1:14; h -44:1:44
x=round(linspace(-14,14,29));

y=round(linspace(-14,14,29));

h=round(linspace(-44,44,89));

% 3D matrices of dimension 29*29*89
% xx,yy,hh were all grid point coordinates in (x,y,h) convention
% id were grid point identity matrix
% css were grid point core-shell attributes matrix

% allocating/defining xx,yy,hh,id,css  
xx=round(ones(29,29,89));
yy=xx;hh=xx;
id=xx;css=xx;

% initializing xx,yy,hh
for i=1:1:29
    xx(i,:,:)=i-15;
    for j=1:1:29
        yy(i,j,:)=j-15;
        for k=1:1:89
            hh(i,j,k)=k-45;
        end
    end
end
%% 
% initializing id
% each grid point has a unique identity of integer value
id=round(xx*10000+yy*100+hh);
% diameters of the sphere where each grid point is sitting at
dd=(xx.^2+yy.^2+xx.*yy-(-1).^hh.*(xx+yy)+1/3)*(2.0*roc/roa)^2+(hh).^2;
% assigning core-shell attibutes to each grid point according to its
% position: inside core -- 0; in the shell -- 1; outside the shell --2.
css(find(dd<=(dcore/roa)^2))=0;
css(find(dd>(dcore/roa)^2 & dd<=(dshell/roa)^2))=1;
css(find(dd>(dshell/roa)^2))=2;
% calculating Euclidean coordinate for each grid point.
X=(xx-yy)/2.0*roc;
Y=sqrt(3)*roc*((xx+yy)/2.0-(-1).^hh/3.0);
Z=hh*roa/2.0;

% initialize core, shell and inert(outside the shell) indexing matrices
core=find(css==0);
shell=find(css==1);
inert=find(css==2);

% count number of lattice grid in core and shell compartment
n_core=numel(find(css==0));
n_shell=numel(find(css==1));
n_inert=numel(find(css==2));

%shell Yb; --> core (Lu Er);
%Shell Yb; 75%grid Yb 25% Na
%shell Na->-1; half of odd grids goes to Na

%random assignment of Na/Yb
%coreNaYb=core(find(mod(id(core),2)==1));
%N_NaYb=numel(coreNaYb);
%css(coreNaYb(find(rand(N_NaYb,1)<=0.5)))=-1;

%fixed assignment
% css(core(find(mod(id(core),4)==1)))=-1;%option 1 for core
% css(shell(find(mod(id(shell),4)==1)))=-1;%option 1 for shell

% css(core(find(mod(id(core),4)==3)))=-1;%option 2 for core
% css(shell(find(mod(id(shell),4)==3)))=-1;%option 2 for shell

%intercalating grid occupancy assignment such that 1f sites were 50%
%occupied, the other half should be occupied by Na+ ions (css-->-1).
css(core(find(mod(xx(core)+yy(core),2)==0 & mod(id(core),4)==1)))=-1;
css(core(find(mod(xx(core)+yy(core),2)==1 & mod(id(core),4)==3)))=-1;

css(shell(find(mod(xx(shell)+yy(shell),2)==0 & mod(id(shell),4)==1)))=-1;
css(shell(find(mod(xx(shell)+yy(shell),2)==1 & mod(id(shell),4)==3)))=-1;

css(inert(find(mod(xx(inert)+yy(inert),2)==0 & mod(id(inert),4)==1)))=-1;
css(inert(find(mod(xx(inert)+yy(inert),2)==1 & mod(id(inert),4)==3)))=-1;

% find indices of RE3+ lattice sites in core and shell compartments.
coyb=find(css==0);
shyb=find(css==1);
inyb=find(css==2);
% number of RE3+ lattice sites in core and shell compartments.
n_coyb=numel(find(css==0));
n_shyb=numel(find(css==1));
n_inyb=numel(find(css==2));

% initializing neighbor identity matrix for all lattice grid points. each
% grid point (x, y, h) could potentially have at most 8 neighbors with site
% identities nb(1:8, x, y, h).
nb=round(ones(8, 29 , 29 , 89));

nb(1,:,:,:)=xx*10000+yy*100+hh+2;
nb(2,:,:,:)=xx*10000+yy*100+hh-2;
nb(3,:,:,:)=(xx+(-1).^(hh+1))*10000+(yy+(-1).^(hh+1))*100+hh+1;
nb(4,:,:,:)=(xx+(-1).^(hh-1))*10000+(yy+(-1).^(hh-1))*100+hh-1;
nb(5,:,:,:)=xx*10000+(yy+(-1).^(hh+1))*100+hh+1;
nb(6,:,:,:)=xx*10000+(yy+(-1).^(hh-1))*100+hh-1;
nb(7,:,:,:)=(xx+(-1).^(hh+1))*10000+yy*100+hh+1;
nb(8,:,:,:)=(xx+(-1).^(hh-1))*10000+yy*100+hh-1;


%% grid points, RE3+ sites counting
% initialize core-shell attibutes of neighbor matrix for compartment
% identification
csnb=nb*0+2; % inert attributes by default
csnb(find(ismembertol(nb,id(core),1e-12)))=0; % neighbors are inside the core
csnb(find(ismembertol(nb,id(shell),1e-12)))=1; % neighbors are in the shell

% break down analysis of RE site numbers 
%coremember who has all neighbours in core
core_core=coyb(find(sum(csnb(:,coyb).^2,1)==0));
%core member who has a neighbour in shell
core_edge=coyb(find(sum(csnb(:,coyb).^2,1)~=0));
%active shell member who has all neighbors in active shell
shell_shell=shyb(find(sum((csnb(:,shyb)-1).^2,1)==0));
%active shell member who has a neighbor not in active shell
shell_edge=shyb(find(sum((csnb(:,shyb)-1).^2,1)~=0));
%active shell member who has a neighbour in core
shell_inner=shell_edge(find(prod(csnb(:,shell_edge),1)==0));
%active shell member who has a neighbour in inert outer shell
shell_outer=shell_edge(find(prod(csnb(:,shell_edge),1)~=0));
%inert outer shell member who has a neighbour in active shell
inert_edge=inyb(find(sum((csnb(:,inyb)-2).^2,1)~=0));

disp(['number of grid points in core, shell, inert compartment']);
disp([n_core n_shell n_inert]);
disp(['number of RE3+ sites in core, shell, inert compartment']);
disp([n_coyb n_shyb n_inyb]);
disp(['number of RE3+ sites of different parts of the nanostructure']);
disp([numel(core_core) numel(core_edge) numel(shell_inner) numel(shell_shell) numel(shell_outer) numel(inert_edge)])

%% core-shell vis
% visualization of RE lattice sites in 3D for core-shell nanostructure
figure(16)
scatter3(X(core_core),Y(core_core),Z(core_core),Z(core_core)*0+1,'.','m','fill');
hold on
scatter3(X(core_edge),Y(core_edge),Z(core_edge),Z(core_edge)*0+9,'o','m','fill');
%scatter3(X(shell_edge),Y(shell_edge),Z(shell_edge),Z(shell_edge)*0+5,'o','m','fill');
%scatter3(X(shell_inner),Y(shell_inner),Z(shell_inner),Z(shell_inner)*0+5,'o','r','fill');
%scatter3(X(shell_outer),Y(shell_outer),Z(shell_outer),Z(shell_outer)*0+5,'o','b','fill');
%scatter3(X(inert_edge),Y(inert_edge),Z(inert_edge),Z(inert_edge)*0+5,'o','m','fill');
%scatter3(X(shell_shell),Y(shell_shell),Z(shell_shell),Z(shell_shell)*0+5,'o','g','fill');
xlabel('x');
ylabel('y');
zlabel('z');

axis equal;

hold on

vis=find(X(shell)<0 | Y(shell)<0 | Z(shell)<0);
scatter3(X(shell(vis)),Y(shell(vis)),Z(shell(vis)),Z(shell(vis))*0+16,'.','b');
hold off;

%%radii analysis of various interfaces
r_core_edge=sqrt(X(core_edge).^2+Y(core_edge).^2+Z(core_edge).^2);
rce=mean(r_core_edge)
std(r_core_edge)

r_shell_inner=sqrt(X(shell_inner).^2+Y(shell_inner).^2+Z(shell_inner).^2);
rsi=mean(r_shell_inner)
std(r_shell_inner)


r_shell_outer=sqrt(X(shell_outer).^2+Y(shell_outer).^2+Z(shell_outer).^2);
rso=mean(r_shell_outer)
std(r_shell_outer)


r_inert_edge=sqrt(X(inert_edge).^2+Y(inert_edge).^2+Z(inert_edge).^2);
rie=mean(r_inert_edge)
std(r_inert_edge)

disp(['core diameter, shell diameter:'])
disp([dcore,dshell])
disp(['core shell diameter estimation from interface radii averaging:'])
disp([rce+rsi, rso+rie])



% display number of RE3+ sites in core or shell compartments
% n_core
% disp(numel(find(css==0)));
% n_shell
% disp(numel(find(css==1)));

%% nbcs init
% initialize core-shell attibutes of neighbor matrix for lattice hopping
nbcs=nb*0-1; % inert attributes by default occupied by Na+
nbcs(find(ismembertol(nb,id(coyb),1e-12)))=0; % neighbors are inside the core
nbcs(find(ismembertol(nb,id(shyb),1e-12)))=1; % neighbors are in the active shell
nbcs(find(ismembertol(nb,id(inyb),1e-12)))=2; % neighbors are in the outer inert compartment

% %%
% %alternative way to initialize core-shell attibutes of neighbor matrix for
% %lattice hopping
% tic
% nbcs=nb*0;
% 
% for i=1:29
%     for j=1:29
%         for k=1:89
%             
% nbid=nb(:,i,j,k);
% xn=round(nbid/10000);
% yn=round((nbid-xn*10000)/100);
% hn=nbid-xn*10000-yn*100+45;
% yn=yn+15;xn=xn+15;
% 
% xn(find(xn<1))=1;xn(find(xn>29))=29;
% yn(find(yn<1))=1;yn(find(yn>29))=29;
% hn(find(hn<1))=1;hn(find(hn>89))=89;
% 
% nbcs(:,i,j,k)=css(sub2ind([29 29 89],xn,yn,hn));
% 
%         end
%     end
% end
% toc
%% monte carlo simulation
%

tic

%hopping probability calculation based on Dexter energy transfer
p38=1.d0/(1.d0/exp(-2.d0*(3.84d0-3.47d0)/0.3d0)*2.d0+6.d0);

p12=(1.d0-p38*6.d0)/2.d0;

pp=[p12,p12,p38,p38,p38,p38,p38,p38];

% inward migration flag: 
% 0--simulation starts from Yb3+ sites in core and stops when reaching
% core-shell interface;
% 1--simulation starts from Yb3+ sites in active shell and stops when
% reaching core-shell interface.

%inward=0;
inward=1;
%sensitizing compartment
sen=0;

%looping index limits & sensitizing sites
nsim=1000%k <-number of simulations (recommend nsim >= 1000 times)
if inward
    nss=n_shyb;%j <-number of Yb3+ sites for absorbing photons(shell)
    sen=shyb;
else
    nss=n_coyb;%j <-number of Yb3+ sites for absorbing photons(core)
    sen=coyb;
end
nstep=5000%i <-maximum number of hopping steps

%trajectory steps init
steps=ones(nss,nsim)*0;

% main loop starts here
for k=1:nsim
    
    % trajectory matrices init
    tr=ones(nss,nstep)*0;

    for j=1:nss

        % trajectory init
        start=sen(j);
        cur=start;

        for i=1:nstep

            tr(j,i)=cur;
            neighbour=nb(:,cur);

            vld=find(nbcs(:,cur)==0 | nbcs(:,cur)==1);

            nvld=numel(vld);

            pv=cumsum([0 pp(vld)/sum(pp(vld))]);

            p=rand;

            for l=1:nvld

                if (p>pv(l) & p <= pv(l+1))

                    icur=neighbour(vld(l));

                end

            end

            cur=find(id==icur);

            if css(cur)==1-inward

                steps(j,k)=i;
                break; 

            end


        end%i



    end%j

    [ms, ims]=max(steps(:,k));
    trx=tr(:,1:ms);

    c=clock;
    txtTimeStamp = sprintf('_%4d_%02d_%02d_%02d%02d%02d', ...
            c(1),c(2),c(3),c(4),c(5),round(c(6)));

    %record trajectories
    if inward
        trf=[path 'trsP6_' num2str(k,'%05i') txtTimeStamp '.h5']

        h5create(trf, '/trs', [nss ms], 'ChunkSize', [nss 1], 'Deflate', 9)
        h5write(trf,'/trs', uint32(trx));

    else
        trf=[path 'trcP6_' num2str(k,'%05i') txtTimeStamp '.h5']

        h5create(trf, '/trc', [nss ms], 'ChunkSize', [nss 1], 'Deflate', 9)
        h5write(trf,'/trc', uint32(trx));
    end
    % trmx=tr(:,ims);
    % 
    % %trc=tr(1:ms,:);
    % %trj=trc(:,100);trj=trj(find(trj~=0));
    % 
    % figure(2)
    % 
    % trj=trmx(find(trmx~=0));
    % 
    % plot3(X(trj),Y(trj),Z(trj),'ro-')
    % xlabel('x');
    % ylabel('y');
    % zlabel('z');
    % axis equal;
    % 
    % hold on 
    % scatter3(X(core_edge),Y(core_edge),Z(core_edge),Z(core_edge)*0+9,'o','m','fill');
    % scatter3(X(shell_outer),Y(shell_outer),Z(shell_outer),Z(shell_outer)*0+5,'o','b','fill');
    % hold off
    % 
    % 
    % 

    figure(3)
    histogram(steps(:,k))
    median(steps(:,k))
    round(mean(steps(:,k)))



end%k


%% simulation output
%save simulation data to hdf5 format
saveh5=1;

c=clock;

if saveh5
    sprintf('Saving Migration Steps(%dx%d, uint16) to HDF5 ...', ...
        nss, nsim);
    txtTimeStamp = sprintf('_%4d_%02d_%02d_%02d%02d%02d', ...
        c(1),c(2),c(3),c(4),c(5),round(c(6)));
    h5path=path;
    if inward
        h5name=['shellstepsP6',txtTimeStamp];
    else
        h5name=['corestepsP6',txtTimeStamp];
    end
    savfile=fullfile(h5path ,[h5name,'.h5']);

    if inward
        h5create(savfile,'/shellsteps', [nss, nsim], 'ChunkSize', [nss 1], 'Deflate', 9);
        h5write(savfile,'/shellsteps', uint16(steps));
    else
        h5create(savfile,'/coresteps', [nss, nsim], 'ChunkSize', [nss 1], 'Deflate', 9);
        h5write(savfile,'/coresteps', uint16(steps));
    end
    
    h5create(savfile,'/id', [29 29 89], 'ChunkSize', [29 29 1], 'Deflate', 9);
    h5write(savfile,'/id', int32(id))

    h5create(savfile,'/css', [29 29 89], 'ChunkSize', [29 29 1], 'Deflate', 9);
    h5write(savfile,'/css', int8(css))
    
    h5create(savfile,'/nb', [8 29 29 89], 'ChunkSize', [8 29 29 1], 'Deflate', 9);
    h5write(savfile,'/nb', int32(nb))

    h5create(savfile,'/nbcs', [8 29 29 89], 'ChunkSize', [8 29 29 1], 'Deflate', 9);
    h5write(savfile,'/nbcs', int8(nbcs))
   

    h5create(savfile, '/ctime', 6); 
    h5write(savfile,'/ctime', c);     
    
    sprintf('Done writing simulation output to HDF5!');
end




toc

load handel
sound(y,Fs)






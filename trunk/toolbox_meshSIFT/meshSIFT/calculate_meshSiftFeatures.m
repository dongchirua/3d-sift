function featurelist = calculate_meshSiftFeatures(vertex,faces,poi,params)
% function to calculate the feature vectors for the points of intrest (poi)
% of a mesh
%
% -----
% INPUT
% -----
% vertex: Mx3 matrix containing all vertices
% faces:  Nx3 matrix containing all faces of the mesh
% poi: Px2 matrix containing the index of the point of intrest and the
%      assigned scale: poi = [index scale]
% params (optional): structure of optional parameters. if none given,
%                    defaults will be used
% params.filename: will be assigned to every structure 
% params.vertices: Mx3xQ matrix containing the vertices of the smoothed
%                  meshes
% params.scales: vector containing the scales of intrest
%
% params.patchmode = 1;  manner of distributing patches around the point of intrest: 1 = square, 2 = logpolar
% params.orientationmode = 1;  mode to calculate the normal (based on normals(1) or on Umax(2)) 
% params.fscale = 3;  factor to scale overall size of the feature
% params.forient = 3;  factor for the size of the neighborhood concerned in defining the orientation. 
% params.fsigorient = 1.5;  factor for the sigma of the weighting factors for defining the orientation. 
% params.fsighist = 1.5;  factor for the sigma of the weighting factors for the histogram
% params.fsighistangle = 0.5;  factor for the sigma for the angle of the weighting factors for the histogram
%                                  (only if patchmode = logpolar)
% params.nbinsorient = 360;  # bins for the histogram to select the orientation. 
% params.nbinsfeature = 8;  # bins for the histogram in the feature vector. 
% params.fdist = 1.5;  factor for the distance from the keypoint, where to center the local patch for the histogram
% params.fpatch = 1.25;  factor for the radius for the local patch of the histogram 
% params.normmode = 1;  mode to normalise the histograms
%                           mode 0 = no clipping
%                           mode 1 = only clipping
%                           mode 2 = clipping and normalising again
%                           mode 3 = only clipping and replacing NaN's by zeros.
% params.curvature_smoothing = 3;  size of smoothing ring for curvature calculation 
% params.progbar = 1 % show progressbar
%
% ------
% OUTPUT
% ------
% featurelist = cell containing for every feature a structure with
%     filename = name of the file
%     kpcoordinates = coordinates of the keypoint
%     kpindex = index of the keypoint
%     kpscale = scale of the keypoint
%     orientation = canonical orientation assigned to the keypoint
%     feature = the assigned feature vector
%     numberOfFeatures = number of features assigned to this keypoint (in
%         case of multiple orientations)
%
% Author: Chris Maes
% 2010/04

%% init
if nargin < 3
    error('not enough parameters given')
end
filename = getoptions(params,'filename','unknown');
scales = getoptions(params,'scales',[]);
vertices = getoptions(params,'vertices',[]);

if size(poi,2)<2
    error('poi must be a Nx2 matrix')
elseif size(poi,2)>2
    warning('MATLAB:paramAmbiguous','poi must be a Nx2 matrix, only the first two columns will be used')
    poi = poi(:,1:2);
end
if isempty(scales)
    scales = unique(poi(:,2));
end
if isempty(vertices)
    vertices = calculate_mesh_smoothing(vertex,faces,scales);
end

%% parameters
% info about the parameters can be found at the top of the m-file or by
% typing 'help calculate_meshSiftFeatures'
progbar = getoptions(params,'progbar',1); % show progressbar
patchmode = getoptions(params,'patchmode',1);
orientationmode = getoptions(params,'orientationmode',1);
fscale = getoptions(params,'fscale',3);
forient = getoptions(params,'forient',3);
fsigorient = getoptions(params,'fsigorient',1.5);
fsighist = getoptions(params,'fsighist',1.5);
fsighistangle = getoptions(params,'fsighistangle',0.5);
nbinsorient = getoptions(params,'nbinsorient',360);
nbinsfeature = getoptions(params,'nbinsfeature',8);
fdist = getoptions(params,'fdist',1.5);
fpatch = getoptions(params,'fpatch',1.25);
normmode = getoptions(params,'normmode',1);
options.curvature_smoothing = getoptions(params,'curvature_smoothing',3);
sfigs = getoptions(params,'sfigs',0); % show figures?

%% initialisation
options.verb = 0; % don't show progressbar for compute_curvature
featurelist = cell(0);
curvinfo = cell(numel(scales),3);
normals = cell(numel(scales),1);
nvert = size(vertex,1);
orientedges = -pi:2*pi/nbinsorient:pi; % edges for the histogram
gfilt = exp(-(-8:8).^2./(2*17^2)); % params sigma = 17, 17bins. these params define following lines as well
sedges = -1:2/nbinsfeature:1; % edges for the shape index histogram
oedges = -pi:2*pi/nbinsfeature:pi; % edges for the orientation histogram

%% preparative calculations
scalesoi = unique(poi(:,2));
for sc = scalesoi'
    ns = find(scales==sc);
    vert = vertices(:,:,ns);
    normals{ns} = compute_normal(vert,faces); % compute normals to vertices and faces
    [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vert,faces,options); % compute k1, k2, H, K,...
    curvinfo{ns,1} = 2/pi*atan((Cmin+Cmax)./(Cmin-Cmax)); % shape index
    curvinfo{ns,2} = sqrt(2*Cmean.^2-Cgauss); % curvedness
    curvinfo{ns,3} = Umax;
end
    
clear scalesoi ns sc Umin Umax Cmin Cmax Cmean Cgauss Normal vert
%% calculations
npoi = size(poi,1);
for ipoi = 1:npoi
    %% init
    if progbar
        progressbar(ipoi,npoi);
    end
    ioi = poi(ipoi,1);
    scale = poi(ipoi,2);
    ns = find(scales==scale);
    vert = vertices(:,:,ns);
    coi = vert(ioi,:);
    scaledscale = fscale*scale; %factor for the scale (influances overall scale)
    sigorient = fsigorient*scaledscale; % sigma for the gaussian used to define the orientation.
    radius = forient*scaledscale; % radius of the regio to concern
    distr = fdist*scaledscale; % distance to go from the keypoint for the centers of the histograms.
    patchr = fpatch*scaledscale; % radius for every patch to concern
    sighist = fsighist*scaledscale; % sigma for the gaussian weighted histograms (radius in both square and logpolar)
    sighistangle = fsighistangle*pi/2; % sigma for the gaussian weighted histograms in logpolar mode (angles)
    if sfigs
        figure()
        [x,y,z] = xyz(vertex);
        trisurf(faces,x,y,z);view(0,90);colormap(gray(256));axis off
        shading interp
        hold on
        [x,y,z] = xyz(vertex(ioi,:));
        scatter3(x,y,z,200,'ro','filled')
        hold off
    end
    
%% calculate orientation
    L = inf*ones(nvert,1);
    L(sum((vert-repmat(coi,nvert,1)).^2,2)>radius^2) = -inf;
    options.constraint_map = L; % don't calculate distances if the euclidic distance is greater than radius
    Dref = perform_fast_marching_mesh(vert, faces, ioi,options); % geodetic distance to the point of intrest
    Dref(Dref>radius)=inf;
    nbhi = find(Dref<radius & Dref~=0); % indices of points in neighborhood, poi excluded
    if numel(nbhi)>1
        gwf = exp(-Dref.^2./(2*sigorient^2)) ; % gaussian weighting factors, defined on the geodetical distance

        nref = dot(normals{ns}(:,nbhi)',repmat(gwf(nbhi),1,3),1); % guassian weighted version of the normal
        nref = nref'/sqrt(sum(nref.^2));%normalize reference normal

        tangles = real(acos(sum(normals{ns}.*repmat(nref,1,nvert),1))); % angles between normals and reference normal

        mproj = eye(3) - nref*nref'; % projection matrix onto plane of reference
        if orientationmode ~= 2 % orientationmode == 1
            projs = mproj*normals{ns}; % projections of the normals onto the plane of reference
        else
            projs = mproj*curvinfo{ns,3}; % projections of direction of maximal curvature onto plane of reference
        end
        if sfigs
            figure();
            tempnbhi = [nbhi; ioi];
            [x,y,z] = xyz(vert(tempnbhi,:));
            dx = max(x)-min(x);
            dy = max(y)-min(y);
            dfactor = .15;
            X=[min(x)-dfactor*dx max(x)+dfactor*dx;min(x)-dfactor*dx max(x)+dfactor*dx];
            Y=[min(y)-dfactor*dy min(y)-dfactor*dy;max(y)+dfactor*dy max(y)+dfactor*dy];
%             X=[min(x) max(x);min(x) max(x)];
%             Y=[min(y) min(y);max(y) max(y)];
            [x,y,z] = xyz(vert(ioi,:)); 
            [u,v,w] = xyz(12*nref');
            Z=(-u*(X-x)-v*(Y-y))/w+z;
            surf(X,Y,Z)
            shading flat
            alpha(0.2)
            
            shw(vert(tempnbhi,:),1);hold on
            quiver3(x,y,z,u,v,w,0,'r','lineWidth',2)
            [x,y,z] = xyz(vert(tempnbhi,:));
            [u,v,w] = xyz(normals{ns}(:,tempnbhi)');
            quiver3(x,y,z,u,v,w,'k','lineWidth',1.5);
            
            cosa = nref(3)/sqrt(nref(3)^2+nref(2)^2); % angle of rotation in XZ-plane
            sina = cosa*nref(2)/nref(3);
            m1 = [1 0 0 0; 0 cosa -sina 0; 0 sina cosa 1;0 0 0 1];
            cosa = sqrt(nref(3)^2+nref(2)^2); % angle of rotation in YZ-plane
            sina = nref(1);
            m2 = [cosa 0 -sina 0; 0 1 0 0; sina 0 cosa 1;0 0 0 1];
            t1 = m2*m1;
            tempvert = t1*[vert(tempnbhi,:) ones(numel(tempnbhi),1)]';
            mtvert = tempvert(3,end);
            m3 = [1 0 0 0; 0 1 0 0; 0 0 0 mtvert;0 0 0 1];
            posproj = inv(t1)*m3*tempvert;
            
            [x,y,z] = xyz(posproj');
            [u,v,w] = xyz(projs(:,tempnbhi)');
            quiver3(x,y,z,u,v,w,'g','lineWidth',1.5);
            axis off
            hold off
            clear posproj x y z u v w X Y Z cosa sina m1 m2 m3 mtransf mtvert t1 tempnbhi
        end
        projs = projs./repmat(sqrt(sum(projs.^2,1)),3,1);%normalize projections
        n0 = projs(:,find(projs(1,:)~=0 & ~isnan(projs(1,:)),1,'first')); % direction in the plane with degree zero - arbitrary chosen.
        try
            sangles = acos(sum(projs.*repmat(n0,1,nvert),1)); % calculation of slant angles resp to n0
            signs = cross(projs,repmat(n0,1,nvert));
            signs = double(sign(signs(1,:))==sign(nref(1))); % signs = 1 if the angle is positive clockwise
            signs(signs==0)=-1; % signs = -1 means a negative angle clockwise. (so counterclockwise)
            sangles = sangles.*signs;
            nout = whist(sangles(nbhi),gwf(nbhi)'.*tangles(nbhi),orientedges); % calculate weighted histogram as close as possible to LO2008
        catch
            nout =  zeros(1,length(orientedges));
        end
        clear signs mproj sangles
        
        ohist = nout;
        for i = 1:3 % 3 times convolution
            ohist = conv(ohist,gfilt);
            ohist(9:16) = ohist(9:16)+ohist(nbinsorient+9:nbinsorient+16);
            ohist(nbinsorient+1:nbinsorient+8) = ohist(nbinsorient+1:nbinsorient+8) + ohist(1:8);
            ohist = ohist(9:nbinsorient+8);
        end
        clear i nout
        imaxs = maxima(ohist,0.8*max(ohist)); % look for all maxima with orientation > 80% of maximum
        thetas = zeros(size(imaxs));
        for mi = 1:numel(imaxs)
            imax = imaxs(mi);
            if imax==1
                imax = 2;
                ths = interp1(1:2:5,orientedges(imax:imax+2),2:2:4); %angles of maximum and neighbours
                ths = [2*ths(1)-ths(2) ths]; %#ok<AGROW>
                histv = ohist([nbinsorient 1 2]);
            elseif imax==nbinsorient
                imax = imax-1;
                ths = interp1(1:2:5,orientedges(imax-1:imax+1),2:2:4); %angles of maximum and neighbours
                ths = [ths 2*ths(2)-ths(1)]; %#ok<AGROW>
                histv = ohist([nbinsorient-1 nbinsorient 1]);
            else
                ths = interp1(1:2:7,orientedges(imax-1:imax+2),2:2:6); %angles of maximum and neighbours
                histv = ohist(imax-1:imax+1);
            end
            abc = [ths'.^2 ths' ones(3,1)]\histv';
            thetas(mi) = -abc(2)/2/abc(1);
        end
        clear mi imax ths histv abc imaxs
        orientations = zeros(length(thetas),3);
        sangles = zeros(length(thetas),nvert);
        cosa = nref(3)/sqrt(nref(3)^2+nref(2)^2); % angle of rotation in XZ-plane
        sina = cosa*nref(2)/nref(3);
        m1 = [1 0 0; 0 cosa -sina; 0 sina cosa];
        cosa = sqrt(nref(3)^2+nref(2)^2); % angle of rotation in YZ-plane
        sina = nref(1);
        m2 = [cosa 0 -sina; 0 1 0; sina 0 cosa];
        t1 = m2*m1;
        t1i = inv(t1);
        for thi = 1:length(thetas)
            theta = thetas(thi);
            cosa = cos(theta); % angle of rotation in XY-plane
            sina =sin(theta);
            m3 = [cosa sina 0; -sina cosa 0; 0 0 1];
            mtransf = t1i*m3*t1;
            orientation = mtransf*n0;
            orientations(thi,:) = orientation;
            if sfigs
                figure(1);
                hold on;
                orvector = orientation*15;
                tempcoi = vertex(ioi,:)
                quiver3(tempcoi(1),tempcoi(2),tempcoi(3),orvector(1),orvector(2),orvector(3),0,'linewidth',1.5)
                axis off
                hold off
                clear orvector tempcoi
            end
        %% Calculate slant angles with reference to the orientations
            mproj = eye(3) - nref*nref'; % projection matrix onto plane of reference
            projs = mproj*normals{ns}; % projections of the normals onto the plane of reference
            projs = projs./repmat(sqrt(sum(projs.^2,1)),3,1);%normalize projections
            sangle = acos(sum(projs.*repmat(orientation,1,nvert),1)); % calculation of slant angles resp to orientation
            signs = cross(projs,repmat(orientation,1,nvert));
            signs = double(sign(signs(1,:))==sign(nref(1))); % signs = 1 if the angle is positive clockwise
            signs(signs==0)=-1; % signs = -1 means a negative angle clockwise. (so counterclockwise)
            sangles(thi,:) = sangle.*signs;
        end
        clear sina cosa m1 m2 m3 mtransf thi theta orientation mproj projs sangle signs

%% Feature description
        dists = Dref(nbhi); % geodetic distances for the neigborhood
        nbh = vert(nbhi,:)'; % the neigborhood around the point of intrest
        nbh = nbh-repmat(coi',1,numel(nbhi)); % vector from point of intrest to every vertex
        nors = size(orientations,1); % number of orientations
        for ior = 1:nors
            orientation = orientations(ior,:);
            feature = zeros(9*2*nbinsfeature,1); % initialize feature vector
%% square patchmode
            if patchmode ~= 2
                icenters = zeros(9,1);
                icenters(1) = ioi;
                if sfigs
                    [x,y,z] = xyz(vert(nbhi,:));
                    shw(vert(nbhi,:));hold on
                end
                for i = 0:7
                    theta = i*pi/4;
                    if rem(i,2)==0
                        dr = repmat(distr,numel(nbhi),1);
                    else
                        dr = repmat(sqrt(2)*distr,numel(nbhi),1);
                    end
                    cosa = cos(theta); % angle of rotation in XY-plane
                    sina = sin(theta);
                    m3 = [cosa sina 0; -sina cosa 0; 0 0 1];
                    mtransf = t1i*m3*t1;
                    dir = mtransf*orientation'; % direction where to search the center
                    nproj = cross(nref,dir); % normal of the plane formed by the reference normal and the search direction.
                    mproj = eye(3) - nproj*nproj'; % matrix to project all vectors onto the plane
                    projs = mproj*nbh; % projections of all points in the neighborhood onto the plane 
                    cosas = sum(projs.*nbh,1)./sqrt(sum(projs.^2,1).*sum(nbh.^2,1)); % cosines of the angles between the points and the canonical orientation
                    signs = cross(nbh,repmat(nref,1,numel(nbhi)));
                    signs = cross(signs,repmat(dir,1,numel(nbhi)));
                    signs = double(sign(signs(1,:))==sign(nref(1))); % signs = 1 if the nbh-vector points in the search direction.
                    cosas = cosas.*signs;
                    l = dr.^2 + dists.^2 - 2*dr.*dists.*cosas'; % squares of the distances between a vertex and the center we would like to find.
                    [m,il] = min(l); % il = index of the vertex closest to the center we would like to find.
                    icenters(2+i) = nbhi(il);
                    if sfigs
                        orvector = dir*dr(1);
                        quiver3(coi(1),coi(2),coi(3),orvector(1),orvector(2),orvector(3),0,'linewidth',2)
                        clear orvector
                    end
                end
                if sfigs
                    shw(vert(icenters,:),1,100,'ro','filled');
                    view(nref)
                    axis off
                end
                clear i theta dr cosa sina m3 mtransf dir nproj mproj projs cosas signs l m il
                Dc = inf*ones(nvert,9);
                L = inf*ones(nvert,1);
                L(sum((vert-repmat(coi,nvert,1)).^2,2)>distr^2) = -inf;
                options.constraint_map = L; % don't calculate distances if the euclidic distance is greater than distr
                Dc(:,1) = perform_fast_marching_mesh(vert,faces,icenters(1),options);
                L = inf*ones(nvert,1);
                L(sum((vert-repmat(coi,nvert,1)).^2,2)>9*distr^2) = -inf;
                options.constraint_map = L; % don't calculate distances if the euclidic distance is greater than 3*distr
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(2),options);
                Dc(Q==icenters(2),2) = D(Q==icenters(2));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(3),options);
                Dc(Q==icenters(3),3) = D(Q==icenters(3));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(4),options);
                Dc(Q==icenters(4),4) = D(Q==icenters(4));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(5),options);
                Dc(Q==icenters(5),5) = D(Q==icenters(5));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(6),options);
                Dc(Q==icenters(6),6) = D(Q==icenters(6));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(7),options);
                Dc(Q==icenters(7),7) = D(Q==icenters(7));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(8),options);
                Dc(Q==icenters(8),8) = D(Q==icenters(8));
                [D SS Q] = perform_fast_marching_mesh(vert,faces,icenters(9),options);
                Dc(Q==icenters(9),9) = D(Q==icenters(9));   
                for ic = 1:numel(icenters)
                    icenter = icenters(ic);
                    gwfc = exp(-Dc(:,ic).^2./(2*sighist^2)) ; % gaussian weighting factors, defined on the geodetical distance
                    nbhic = (Dc(:,ic)<patchr & Dref<radius);
                    shist = whist(curvinfo{ns,1}(nbhic),gwf(nbhic).*gwfc(nbhic).*curvinfo{ns,2}(nbhic),sedges); % histogram of shape index values
                    shist = normclip(shist,normmode); % normalised and clipped histogram
                    ohist = whist(sangles(ior,nbhic),gwf(nbhic)'.*gwfc(nbhic)'.*tangles(nbhic),oedges); % calculate weighted histogram. TODO: same form as in beginning!
                    ohist = normclip(ohist,normmode); % normalised and clipped histogram
                    hist = normclip([shist ohist],0); % normalised and clipped histogram % changed 3/3/2010
                    feature((ic-1)*16+1:(ic-1)*16+16) = hist;
                end
                clear L Dc D SS Q ic icenter gwfc nbhic shist ohist hist
                featurelist{end+1} = struct('filename',filename,...
                    'kpindex',ioi,'kpcoordinates',vertex(ioi,:),'orientation',orientation,...
                    'feature',feature,'numberOfFeatures',nors,'scale',scale); %#ok<AGROW>
%% logpolar patchmode
            else
                % center
                nbhic = Dref<patchr;
                shist = whist(curvinfo{ns,1}(nbhic),gwf(nbhic).*gwf(nbhic).*curvinfo{ns,2}(nbhic),sedges); % histogram of shape index values
                shist = normclip(shist,normmode); % normalised and clipped histogram
                ohist = whist(sangles(ior,nbhic),gwf(nbhic)'.*gwf(nbhic)'.*tangles(nbhic),oedges); % calculate weighted histogram.
                ohist = normclip(ohist,normmode); % normalised and clipped histogram
                hist = normclip([shist ohist],0); % normalised and clipped histogram
                feature(1:16) = hist;
                for ia = 0:3
                    histangle = pi/2*ia;
                    sanglesoi = sangles(ior,:)-histangle;
                    sanglesoi(sanglesoi<-pi) = sanglesoi(sanglesoi<-pi)+2*pi;
                    nbhica = sanglesoi>(-pi/2) & sanglesoi<(pi/2);
                    gwfa = exp(-sanglesoi.^2./(2*sighistangle^2)) ; % gaussian weighting factors, defined on the geodetical distance
                    for ir = 1:2
                        dr = ir*distr;
                        gwfr = exp(-(Dref-dr).^2./(2*sighist^2)) ; % gaussian weighting factors, defined on the geodetical distance
                        nbhic = nbhica' & Dref>(ir-1)*distr & Dref<(ir+1)*distr;
                        shist = whist(curvinfo{ns,1}(nbhic),gwf(nbhic).*gwfr(nbhic).*gwfa(nbhic)'.*curvinfo{ns,2}(nbhic),sedges); % histogram of shape index values
                        shist = normclip(shist,normmode); % normalised and clipped histogram
                        ohist = whist(sangles(ior,nbhic),gwf(nbhic)'.*gwfr(nbhic)'.*gwfa(nbhic).*tangles(nbhic),oedges); % calculate weighted histogram.
                        ohist = normclip(ohist,normmode); % normalised and clipped histogram
                        hist = normclip([shist ohist],0); % normalised and clipped histogram % changed 3/3/2010
                        feature(ia*32+ir*16+1:ia*32+(ir+1)*16) = hist;
                    end % radius 1 and 2
                end % angles 0 to 3pi/2
                clear ia histangle nbhica sanglesoi gwfa ir dr gwfr nbhic shist ohist hist
                featurelist{end+1} = struct('filename',filename,...
                    'kpindex',ioi,'kpcoordinates',vertex(ioi,:),'orientation',orientation,...
                    'feature',feature,'numberOfFeatures',nors,'scale',scale); %#ok<AGROW>
            end %patchmode
        end
        clear ior nors dists nbh
        if sfigs
            figure();
            bar(feature)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normhist = normclip(hist,mode)
% mode 0 = no clipping
% mode 1 = only clipping
% mode 2 = clipping and normalising again
% mode 3 = only clipping and replacing NaN's by zeros.
if nargin < 2
    mode = 1;
end
lim = 1/sqrt(numel(hist));
normhist = hist/sqrt(sum(hist.^2));
if mode == 0
    return
end
ind = normhist>lim;
normhist(ind) = lim;
if mode == 1
    return
elseif mode == 2
    normhist = normhist/sqrt(sum(normhist.^2));
    return
elseif mode == 3
    normhist(isnan(normhist))=0;
    return
end
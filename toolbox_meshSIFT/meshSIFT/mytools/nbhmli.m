function indices = nbhmli(faces,index,r)

% generates the neigborhood of a point at a given index in the vertex
% array.
% vertex = nx3 matrix with the vertex points from a mesh
% faces = mx3 matrix with the corner indices of the faces of a mesh
% index = index of point of intrest in vertex
% r = the radius of the neighborhood
% nbh = list of indices in the neighborhood of vertex(index)
% Author: Chris Maes
% 2009/10

%% calculation indices
if r == 1
    r1 = find(faces(:,1)==index)';
    r2 = find(faces(:,2)==index)';
    r3 = find(faces(:,3)==index)';
    nel = 2*(length(r1)+length(r2)+length(r3));
    indices = zeros(1,nel); %#ok<AGROW>
    j = 1;
    for r = r1
        if isempty(find(indices==faces(r,2),1))
            indices(j) = faces(r,2);
            j = j+1;
        end
        if isempty(find(indices==faces(r,3),1))
            indices(j) = faces(r,3);
            j = j+1;
        end
    end
    for r = r2
        if isempty(find(indices==faces(r,1),1))
            indices(j) = faces(r,1);
            j = j+1;
        end
        if isempty(find(indices==faces(r,3),1))
            indices(j) = faces(r,3);
            j = j+1;
        end
    end
    for r = r3
        if isempty(find(indices==faces(r,1),1))
            indices(j) = faces(r,1);
            j = j+1;
        end
        if isempty(find(indices==faces(r,2),1))
            indices(j) = faces(r,2);
            j = j+1;
        end
    end
    indices = indices(1:j-1);
else
    indices = index;
    for n = 1:r
        for index = indices
            ind = nbhmli(faces,index,1);
            cleanind = zeros(size(ind));
            j = 1;
            for i = ind
                if isempty(find(indices==i,1))
                    cleanind(j) = i;
                    j = j+1;
                end
            end
            cleanind = cleanind(1:j-1);
            if j~=1
                indices = [indices cleanind]; %#ok<AGROW>
            end
        end
    end
end

end

    
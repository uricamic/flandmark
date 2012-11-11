function [ cnt, S, bad_idx ] = paramStats( annotation_struct, options, verbose )
%PARAMSTATS Collects statistics about possible positions of each landmark
%point (component)
%   Detailed explanation goes here
%  
% INPUT:
%   image      ... array of cells, each cell containes bbox found by detector and filename
%   bw         ... vector [w; h] contains size of base window. 
%   components ... matrix of size 2xM, M is count of components [w0 ... wM-1; h0 ... hM-1]. Order of components is important!
%   bw_margin  ... vector [w; h] base window margin in percents of original bbox.
%   image_path ... path to directory containing image database. Must contain subdir /mgt with xml annotations. 
%   verbose    ... true for image output with landmarks and components depicted
% 
% OUTPUT:
%   cnt     ... count of images that passed test
%   S       ... matrix 4xM [x0_min; y0_min; x0_max; y0_max; ... xM-1_min; yM-1_min; xM-1_max; yM-1_max] bbox for each component
%   bad_idx ... indices to image array of images that does not passed this
%               test
% 
% 16-07-10 Michal Uricar
% 10-04-11 Michal Uricar, estimation of S0 fixed (face center intead of nose)
% 12-07-11 Michal Uricar, corners dataset
% 21-03-12 Michal Uricar, LFW annotation in one file

    if (nargin < 3)
        verbose = false;
    end;       
    
    bad_idx = [];
    
    % count of components
    M = size(options.components, 2);
    
    cnt = 0;
    S = [inf(2, M); -inf(2, M)];
    
    % for each image in db
    N = annotation_struct.N;
    for i = 1 : N        
        [Iframe, Annotation, I, itmp.bbox, bbox, OrigPoints] = getImageFrame(options, i, annotation_struct);
        
        % !!! Annotation is empty at i == 578 !!!
        if (isempty(Annotation))
            bad_idx = [bad_idx i];
            continue;
        end;
        
        Points = prepareS0gt(Annotation.P, options);    % transform nose to the center of face
        Points = Points(:, options.comselect);          % extract relevant points only (name list in options.compnames)
        Points(:, M) = Annotation.P(:, 10);             % copy back original nose position

        %% Check for each component if it fits in normalized frame
        flag = true;        
        for j = 1 : M
            if ( ((Points(1, j) - options.components(1, j)/2) < 0) || ((Points(1, j) + options.components(1, j)/2) > options.bw(1)) || ...
                 ((Points(2, j) - options.components(2, j)/2) < 0) || ((Points(2, j) + options.components(2, j)/2) > options.bw(2)) )
                flag = false;
            end;
        end;
        
        % increase number of valid images if flag is true
        if (flag)
            cnt = cnt + 1;
            fprintf('%.2f%% - Passed: tested image no.%d file %s...\n', i*100/N, i, Annotation.image.filename);
            
            %% recompute AABB for S
            for j = 1 : M
                bb = [Points(1, j) - options.components(1, j)/2 Points(2, j) - options.components(2, j)/2 ...
                      Points(1, j) + options.components(1, j)/2 Points(2, j) + options.components(2, j)/2 ];
                S(1, j) = min(S(1, j), bb(1));
                S(2, j) = min(S(2, j), bb(2));
                S(3, j) = max(S(3, j), bb(3));
                S(4, j) = max(S(4, j), bb(4));
            end;
            
        else 
            fprintf('%.2f%% - Not passed: tested image no.%d file %s...\n', i*100/N, i, Annotation.image.filename);
            bad_idx = [bad_idx i];
        end;        
        
        %% Visualization
        if (verbose)
            aabb = makeAABB(itmp.bbox);
            aabb2 = makeAABB(bbox);

            % show original image with annotation
            figure(1);
            subplot(2, 2, 1);
            imshow(I, []); hold on;
            % show ground truth points
%             for k = 1 : M
            plot(OrigPoints(1, :), OrigPoints(2, :), 'r.', 'LineWidth', 2, 'MarkerSize', 10);
%             end;
            % show aabb found by detector
            line([aabb(1,:) aabb(1,1)], [aabb(2,:) aabb(2,1)], 'color', 'b');
            line([aabb2(1,:) aabb2(1,1)], [aabb2(2,:) aabb2(2,1)], 'color', 'y');
            hold off;

            % show normalized frame
            subplot(2, 2, 2)
            imshow(Iframe, []); hold on;
            for k = 1 : M
                plot(Points(1, k), Points(2, k), 'r.', 'LineWidth', 2, 'MarkerSize', 10);
                bb = makeAABB([Points(1, k) - options.components(1, k)/2 Points(2, k) - options.components(2, k)/2 ...
                               Points(1, k) + options.components(1, k)/2 Points(2, k) + options.components(2, k)/2 ]);
                line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', 'y');
            end;
            hold off;
            
            % show recomputed S also
            if (max(S(1, :)) < inf && min(S(3, :)) > -inf)                
                subplot(2, 2, 3)
                imshow(Iframe, []); hold on;
                for k = 1 : M
                    plot(Points(1, k), Points(2, k), 'r.', 'LineWidth', 2, 'MarkerSize', 10);
                    bb = makeAABB(S(:, k));
                    line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', 'b');
                    bb = makeAABB([Points(1, k) - options.components(1, k)/2 Points(2, k) - options.components(2, k)/2 ...
                                   Points(1, k) + options.components(1, k)/2 Points(2, k) + options.components(2, k)/2 ]);
                    line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', 'y');
                end;
                hold off;
            end;
            
%             saveas(gcf, ['./img/' 'bw_' Annotation.image.filename]);
%             close gcf;
        end;
    end;
end


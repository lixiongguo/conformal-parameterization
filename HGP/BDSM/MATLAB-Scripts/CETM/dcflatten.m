function [vt, numbrokentriangs, flaterr, els, output] = dcflatten(f, nv, edgelens, maxNumIter, fid, varargin)

if numel(varargin)>0
    Theta = varargin{1};
end

if numel(varargin)>1
    ufixed = varargin{2};
end

if numel(varargin)>2
    bLSCM = varargin{3};
else
    bLSCM = 1;
end

output = {};

% count and display number of faces and vertices (and number of angles)
numfaces = size(f, 1);
numangles = 3 * numfaces;
numvertices = nv;
fprintf(fid, '%u faces, %u vertices.\n', numfaces, numvertices);


% calculate array loglength of same dimensions as f
% loglength(m, n) is the length of the edge of triangle m opposite
% its nth vertex (n = 1,2,3). 
loglength = reallog(edgelens);

% Calculate the (directed) adjacency matrix. adjacencymatrix(m,n) = 1 if the oriented
% boundary of a triangle contains the directed edge from vertex m to vertex
% n, and 0 otherwise. This matrix is not quite symmetric due to boundary edges.
adjacencymatrix = sparse([f(:,1); f(:,2); f(:,3)], ...
                         [f(:,2); f(:,3); f(:,1)], ...
    	                 ones(3 * numfaces, 1), ...
                         numvertices, numvertices, 3 * numfaces);
if any(any(adjacencymatrix > 1))
    error('Triangles must be oriented consistently.')
end
% Use adjacencymatrix to find the boundaryvertices.
boundaryvertices = full( any(adjacencymatrix ~= adjacencymatrix', 2) );

% adjacencymatrix is not needed anymore.
clear adjacencymatrix;
interiorvertices = ~boundaryvertices;
fprintf(fid, '%u boundary vertices, %u interior vertices.\n', nnz(boundaryvertices), nnz(interiorvertices));

% u_per_triangle is a (3 * numfaces) by numvertices matrix which is used to
% distribute the u-values. u_per_triangles(m, n) is 1 if f(m) == n, otherwise 0. 
% (Here, f is indexed as linear vector.)
u_per_triangle = sparse(1 : numangles, f(:), ones(numangles,1), numangles, numvertices, numangles);
% Allocate numfaces by 3 matrices which are used in dcfunctional.
upt = zeros(numfaces, 3);
angles = zeros(numfaces, 3);
ct = zeros(numfaces, 3);

% The following variables are used for the statistics that outfun displays.
numbrokentriangs = uint32(0);
thisfunctionvalue = 0;
lastfunctionvalue = 0;

% The discrete conformal functional with gradient and hessian.
    function [val, grad, hess] = dcfunctional(u, Theta)
        % Cout broken triangs for the statistics.
        numbrokentriangs = 0;
        % upt(m,n) is the u-value of the nth vertex in triangle n.
        upt(:) = u_per_triangle * u;
        % newloglenth(m,n) is the new logarithmic length of the edge
        % opposite the nth vertex in the mth triangle.
        newloglength = loglength + upt(:, [2, 3, 1]) + upt(:, [3, 1, 2]);
        % angles(m,n) is the angle at the nth vertex of the mth triangle.
        % ct(m,n) is the corresponding cotan (or zero if triangle is broken).
%         [angles(:, 1), angles(:, 2), angles(:, 3), ct(:, 1), ct(:, 2), ct(:, 3)] = ...
%             arrayfun(@triangle_angles, newloglength(:,1), newloglength(:,2), newloglength(:,3));
        
        [angles, numbrokentriangs]= meshAnglesFromFaceEdgeLen2( exp(newloglength*2) );
        ct = cot( angles );
        ct(any(ct>2e8 | ct<-2e8, 2),:) = 0;  % avoid NAN
        
        % Calculate the value of the functional.
%         val = 2 * pi * sum(u) + -pi * sum(upt(:)) + sum(angles(:) .* newloglength(:)) + 0.5 * sum(clausen(2 * angles(:)));
        val = dot(Theta, u) + -pi * sum(upt(:)) + sum(angles(:) .* newloglength(:)) + 0.5 * sum(clausen(2 * angles(:)));
        % Bookkeeping for the statistics
        lastfunctionvalue = thisfunctionvalue;
        thisfunctionvalue = val;
        % Calculate the gradient.
        grad = Theta - (u_per_triangle' * angles(:));
        % Build the Hessian.
%         ii = [ f(:, 1);  f(:, 2); f(:, 3); f(:, 1); f(:, 2); f(:, 2); f(:, 3); f(:, 3); f(:, 1)];
%         jj = [ f(:, 1);  f(:, 2); f(:, 3); f(:, 2); f(:, 1); f(:, 3); f(:, 2); f(:, 1); f(:, 3)];
%         hh = [ct(:, 2) + ct(:, 3); ...
%             ct(:, 3) + ct(:, 1); ...
%             ct(:, 1) + ct(:, 2); ...
%             -ct(:, 3); ...
%             -ct(:, 3); ...
%             -ct(:, 1); ...
%             -ct(:, 1); ...
%             -ct(:, 2); ...
%             -ct(:, 2)];
%         hess = sparse(ii, jj, hh, numvertices, numvertices);
        
        hess = sparse(f(:, [1 2 3 1 2 2 3 3 1]), f(:, [1 2 3 2 1 3 2 1 3]), [ct(:,[2 3 1])+ct(:,[3 1 2]) -ct(:, [3 3 1 1 2 2])], numvertices, numvertices);

%         hess0 = - sparse( f(:,[2 3 1 3 1 2]), f(:,[3 1 2 2 3 1]), [ct ct], nv, nv );
%         hess0 = spdiags(-sum(hess0,2), 0, hess0);
    end

    function [alpha, beta, gamma, cota, cotb, cotc] = triangle_angles(loga, logb, logc)
        a = exp(loga);
        b = exp(logb);
        c = exp(logc);
        s0 =  a + b + c;
        s1 = -a + b + c;
        s2 =  a - b + c;
        s3 =  a + b - c;
        if s1 <= 0 || s2 <= 0 || s3 <= 0
            numbrokentriangs = numbrokentriangs + uint32(1);
            alpha = pi * (s1 <= 0);
            beta  = pi * (s2 <= 0);
            gamma = pi * (s3 <= 0);
            cota = 0;
            cotb = 0;
            cotc = 0;
            return;
        end
        alpha = 2 * atan(realsqrt(s2 * s3 / (s1 * s0)));
        beta  = 2 * atan(realsqrt(s3 * s1 / (s2 * s0)));
        gamma = 2 * atan(realsqrt(s1 * s2 / (s3 * s0)));
        p = 0.5 / realsqrt(s1 * s2 * s3 * s0);
        cota = p * (s1 * s0 - s2 * s3);
        cotb = p * (s2 * s0 - s3 * s1);
        cotc = p * (s3 * s0 - s1 * s2);
    end

% allocate vector u used in targetfunction.
u = zeros(numvertices, 1);

if ~exist('Theta', 'var') || isempty(Theta)  % default, free boundary condition
    Theta = ones(numvertices, 1)*2*pi;
%     Theta(boundaryvertices) = 0;
%     Theta = sparse(Theta);

%     ufix = sparse(numvertices, 1);
    ufixed = [find(boundaryvertices) zeros( nnz(boundaryvertices), 1)];
end

if ~exist('ufixed', 'var')
    ufixed = zeros(0, 2);
end

% Zohar: I'm disabling this to handle closed surfaces.
%assert( all( abs(Theta(interiorvertices)-2*pi)<eps*10 ) );

fixid = ufixed(:,1);
% varid = find(Theta);
varid = setdiff(1:numvertices, fixid);

u(fixid) = ufixed(:,2);

    % Clip boundaryvertices out of dcfunctional.
    function [y, g, h] = targetfunction(x)
        u(varid) = x;
        [y, g, h] = dcfunctional(u, Theta);
       
        g = g(varid);
        h = h(varid, varid);
    end

% Prepare for the minimization.
xstart = zeros(numel(varid), 1);
tolgrad = 1e-7;

    % Output function which displays statistics and provides stopping criterion.
    function stop = outfun(x, optimValues, state)
        stop = false;
        %return;
        switch state
            case 'init'
                fprintf(fid, '\n          func value    inf norm of                       broken       cg\n');  
                fprintf(fid,   'iter       increase      gradient      max x    min x    triangles    iter\n\n');
            case 'iter'
                fprintf(fid, '%4u    %12g    %11g    %5.1g    %5.1g    %5u       %4u\n', ...
                        optimValues.iteration, ...
                        thisfunctionvalue - lastfunctionvalue, ...
                        optimValues.firstorderopt, ...
                        max(x), ...
                        min(x), ...
                        numbrokentriangs, ...
                        optimValues.cgiterations);
                if (norm(optimValues.gradient, Inf) <= tolgrad)
                    stop = true;
                    fprintf(fid, 'Max norm of gradient <= %g.\n\n', tolgrad);
                end
        end
    end

    function h = hess_target(x, lambda)
        [~, ~, h] = targetfunction(x);
    end

    function g = get_obj_grad(f, x)
        [~, g] = f(x);
    end

    function f = get_obj_fun(f, x)
        f = f(x);
    end
    
% Set optimization options. TolFun and TolX are set to 0 because stopping
% criterion is provided by outfun.
options = optimset(...
    'MaxIter', maxNumIter, ...
    ...% 'Algorithm', 'interior-point',...
    'Algorithm', 'trust-region-reflective',...
    'HessFcn', @hess_target,...
    'GradObj', 'on', ...
    ...%'Hessian', 'on', ...
    'Hessian', 'user-supplied', ...
    'LargeScale', 'on', ...
    ... %'DerivativeCheck', 'on', ...
    ... %'FunValCheck', 'off', ...
    'TolFun', 1e-10, ...
    ...%'StepTolerance', 1e-6, ...
    'TolX', 0, ... 
    'TolPCG', 1.0e-3, ...
    'PrecondBandWidth', Inf, ...
    'OutputFcn', @outfun, ...
    ...%'Display', 'iter-detailed', ...
    'Display', 'off', ...
    ...%'Diagnostics', 'on', ...
    'UseParallel','always');

if numel(varid) == 0
    disp('All u is fixed');
else
    % Minimize!
    lb = -1e3*ones(size(xstart));
    ub = 1e3*ones(size(xstart));
    if 0
        [xsol, fval, exitflag, output] = fminunc(@targetfunction, xstart, options);
        %[xsol] = minFunc(@targetfunction, xstart, options); 
    elseif 1
        [xsol, fval, exitflag, output] = fmincon(@targetfunction, xstart, [], [], [], [], lb, ub, [], options);
    else % OptiToolbox
        options2 = optiset(... % run optiset to view all options
            'solver', 'ipopt', ...
            ... % 'solverOpts', ipoptset('linear_solver', 'pardiso', 'max_soc', 10), ... % run ipoptset for all options
            ... % 'solver', 'lbfgsb', ...
            'maxfeval', 1e9, ...
            'maxtime', 1e9, ...
            'maxiter', maxNumIter, ...
            ... % 'iterfun', @iterfun, ...
            ... % 'derivCheck', 'on', ...
            'warnings', 'all', ... 
            'display', 'iter' ...
            );
        g_obj = @(x) get_obj_grad(@targetfunction, x);
        f_obj = @(x) get_obj_fun(@targetfunction, x);    
        h_obj = @(x, sigma, lambda) sigma*hess_target(x, lambda);
        %Opt = opti('fun', f_obj, 'grad', g_obj, 'hess', h_obj, 'ndec', length(xstart), 'options', options2, 'bounds', lb, ub);
        Opt = opti('fun', f_obj, 'grad', g_obj, 'ndec', length(xstart), 'options', options2, 'bounds', lb, ub);
        [xsol, fval, exitflag, output] = solve(Opt, xstart);
    end
    output.exitflag = exitflag;
    fprintf(fid, [output.message '\n']);
    u(varid) = xsol;
end

% Don't need these anymore.
clear upt u_per_triangle angles ct

els = exp(loglength + u(f(:, [2 3 1])) + u(f(:, [3 1 2])));

if bLSCM
    % Lay out the flat mesh.
    fprintf(fid, 'laying out flattened mesh ... \n');
    vt = lscm(els, f);
else
    vt = zeros(nv,2);
end

flaterr = norm( 1 - sqrt( meshFaceEdgeLen2s(vt, f) )./els, 'fro' );




% % triangforedge(m, n) is the triang containing directed edge from vertex m
% % to vertex n, or 0 if no such edge exists.
% triangforedge = sparse([f(:,1), f(:,2), f(:,3)], ...
%                        [f(:,2), f(:,3), f(:,1)], ...
%                        [1:numfaces, 1:numfaces, 1:numfaces], numvertices, numvertices, 3 * numfaces);
% % edgelength(m, n) is the length of directed edge from vertex m to n.
% edgelength = sparse([f(:,1); f(:,2); f(:,3)], ...
%                     [f(:,2); f(:,3); f(:,1)], ...
%                     exp([loglength(:, 3) + u(f(:, 1)) + u(f(:, 2)); ...
%                          loglength(:, 1) + u(f(:, 2)) + u(f(:, 3)); ...
%                          loglength(:, 2) + u(f(:, 3)) + u(f(:, 1))]));
% % Allocate vt for vertex coordinates of flat mesh. Third coordinate is
% % zero. It is there because this facilitates displaying the flat mesh.
% vt = zeros(numvertices, 3);
% % edgeslopte(m, n) is to hold the slope angle of directed edge from vertex
% % m to n.
% edgeslope = double(triangforedge | triangforedge'); % sparse matrix with given sparsity pattern.                     
% 
% traversedualspanningtree(@travroot, @travleft, @travright);
% 
%     function traversedualspanningtree(traverserootedge, traverseleftedge, traverserightedge)
%         % init edge queue
%         edgequeue.size = numfaces;
%         edgequeue.data = zeros([2, numfaces], 'uint32');
%         edgequeue.i1 = uint32(0);
%         edgequeue.i2 = uint32(0);
% 
%         function pushedge(edge)
%             if edgequeue.i2 - edgequeue.i1 >= edgequeue.size
%                 error('Edge queue is full.');
%             end
%             edgequeue.data(:, mod(edgequeue.i2, edgequeue.size) + 1) = edge;
%             edgequeue.i2 = edgequeue.i2 + 1;
%         end
% 
%         function edge = popedge()
%             if edgequeue.i1 == edgequeue.i2
%                 error('Edge queue is empty.');
%             end
%             edge = edgequeue.data(:, edgequeue.i1 + 1);
%             edgequeue.i1 = edgequeue.i1 + 1;
%             if edgequeue.i1 >= edgequeue.size
%                 edgequeue.i1 = edgequeue.i1 - edgequeue.size;
%                 edgequeue.i2 = edgequeue.i2 - edgequeue.size;
%             end
%         end
% 
%         facetag = false(numfaces, 1);
%         roottriang = uint32(1);
%         rootedge = f(roottriang, [1,2]);
%         facetag(roottriang) = true;
%         pushedge(rootedge);
%         traverserootedge(rootedge);
%         oppedge = rootedge([2,1]);
%         oppface = triangforedge(oppedge(1), oppedge(2));
%         if (oppface > 0)
%             facetag(oppface) = true;
%             pushedge(oppedge);
%         end
% 
%         while edgequeue.i1 ~= edgequeue.i2 % edge queue not empty
%             edge = popedge();
%             face = triangforedge(edge(1), edge(2));
%             switch f(face, 1)
%                 case edge(1)
%                     v3 = f(face, 3);
%                 case edge(2)
%                     v3 = f(face, 2);
%                 otherwise
%                     v3 = f(face, 1);
%             end
%             leftedge = [edge(1); v3];
%             leftface = triangforedge(leftedge(1), leftedge(2));
%             rightedge = [v3; edge(2)];
%             rightface = triangforedge(rightedge(1), rightedge(2));
%             if (leftface > 0 && ~facetag(leftface))
%                 facetag(leftface) = true;
%                 pushedge(leftedge);
%             end
%             traverseleftedge(leftedge, edge);
%             if (rightface > 0 && ~facetag(rightface))
%                 facetag(rightface) = true;
%                 pushedge(rightedge);
%             end
%             traverserightedge(rightedge, edge);
%         end
%     end
% 
%     function travroot(edge)
%         i1 = edge(1);
%         i2 = edge(2);
%         edgeslope(i1, i2) = 0;
%         edgeslope(i2, i1) = pi;
%         x = edgelength(i1, i2);
%         vt(edge, :) = [0, 0, 0; 
%                        x, 0, 0];
%     end
% 
%     function travleft(edge2, edge1)
%         i1 = edge1(1);
%         i2 = edge1(2);
%         i3 = edge2(2);
% 
%         c = full(edgelength(i1, i2)); % without the full, realsqrt below complains.
%         a = full(edgelength(i2, i3));
%         b = full(edgelength(i3, i1));
%         alpha = 2 * atan(realsqrt(max((a - b + c) * (a + b - c) / ((-a + b + c) * (a + b + c)), 0)));
%         slope = edgeslope(i1, i2) + alpha;
%         edgeslope(i1, i3) = slope;
%         edgeslope(i3, i1) = slope - pi;
%         vt(i3, :) = vt(i1, :) + b * [cos(slope), sin(slope), 0];
%     end
% 
%     function travright(edge2, edge1)
%         i1 = edge1(1);
%         i2 = edge1(2);
%         i3 = edge2(1);
% 
%         c = full(edgelength(i1, i2)); % without the full, realsqrt below complains.
%         a = full(edgelength(i2, i3));
%         b = full(edgelength(i3, i1));
%         beta = 2 * atan(realsqrt(max((-a + b + c) * (a + b - c) / ((a - b + c) * (a + b + c)), 0)));
%         slope = edgeslope(i1, i2) - beta;
%         edgeslope(i3, i2) = slope;
%         edgeslope(i2, i3) = slope + pi;
%         vt(i3, :) = vt(i2, :) - a * [cos(slope), sin(slope), 0];
%     end
% 
% fprintf(fid, 'done.\n');
% 
% 
% % Don't need the 3rd vt coordinate any longer.
% vt(:,3) = [];

end % of function dcflatten


% SUBFUNCTIONS ----------------------------------------

% ------------------------------------------------------------------

function y = clausen(x)
%CLAUSEN Clausen's integral

% take equivalent x-value between -pi and pi
x = mod(x + pi, 2 * pi) - pi;

zerox = (x == 0);
smallx = (~zerox & abs(x) <= 2.0944);
bigx = ~(zerox | smallx);

x(bigx) = x(bigx) - pi * sign(x(bigx));
xx = x .* x;

y = zeros(size(x));

y(smallx) = ((((((((((((                    ...
    2.3257441143020875e-22 * xx(smallx)     ...
    + 1.0887357368300848e-20) .* xx(smallx)  ...
    + 5.178258806090624e-19) .* xx(smallx)   ...
    + 2.5105444608999545e-17) .* xx(smallx)  ...
    + 1.2462059912950672e-15) .* xx(smallx)  ...
    + 6.372636443183181e-14) .* xx(smallx)   ...
    + 3.387301370953521e-12) .* xx(smallx)   ...
    + 1.8978869988971e-10) .* xx(smallx)     ...
    + 1.1482216343327455e-8) .* xx(smallx)   ...
    + 7.873519778281683e-7) .* xx(smallx)    ...
    + 0.00006944444444444444) .* xx(smallx)  ...
    + 0.013888888888888888) .* xx(smallx)    ...
    - reallog(abs(x(smallx))) + 1.0) .* x(smallx);

y(bigx) = ((((((((((((                                ...
				  3.901950904063069e-15 * xx(bigx)    ...
				+ 4.566487567193635e-14) .* xx(bigx)   ...
				+ 5.429792727596476e-13) .* xx(bigx)   ...
				+ 6.5812165661369675e-12) .* xx(bigx)  ...
				+ 8.167010963952222e-11) .* xx(bigx)   ...
				+ 1.0440290284867003e-9) .* xx(bigx)   ...
				+ 1.3870999114054669e-8) .* xx(bigx)   ...
				+ 1.941538399871733e-7) .* xx(bigx)    ...
				+ 2.927965167548501e-6) .* xx(bigx)    ...
				+ 0.0000496031746031746) .* xx(bigx)   ...
				+ 0.0010416666666666667) .* xx(bigx)   ...
				+ 0.041666666666666664) .* xx(bigx)    ...
				+ -0.693147180559945) .* x(bigx);
            
end % of function clausen

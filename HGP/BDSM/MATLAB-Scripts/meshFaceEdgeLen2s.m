function l_sq = meshFaceEdgeLen2s(x, t)
%for each triangle it gives the three edge lengths *squared*

if ~isreal(x)
    x = [real(x) imag(x)];
end

frownorm2 = @(M) sum(M.^2, 2);

l_sq = frownorm2( x(t(:, [2 3 1]), :) - x(t(:, [3 1 2]), :) );
l_sq = reshape( l_sq, [], 3 );

% l_sq = reshape( sum((x(t(:, [2 3 1]), :) - x(t(:, [3 1 2]), :)).^2, 2), [], 3 );
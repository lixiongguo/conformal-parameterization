% angles from edgeLen^2
% given three edge lens a, b, c of a triangle, the angle opposite to a is:
% alpha = acos((b^2 + c^2 - a^2)/(2*b*c))
function [angles, nbroktris] = meshAnglesFromFaceEdgeLen2(el2, bFixRange)

if nargin < 2
    bFixRange = 1;
end
    
%coss = el2*[-1 1 1; 1 -1 1; 1 1 -1]*0.5 ./ sqrt(el2(:, [2 3 1]).*el2(:, [3 1 2]));
coss = el2*[-1 1 1; 1 -1 1; 1 1 -1]*0.5 ./ (el2(:, [2 3 1]) .* el2(:, [3 1 2])).^0.5; % yalmip

% assert( all(abs(coss(:))<1 + eps) );

if bFixRange
    bti = any(abs(coss)>1 + eps, 2);
    nbroktris = sum(bti);
    coss( bti, : ) = 1;
else
    nbroktris = 0;
end

% bti1 =  coss(:)>1 + eps;
% bti2 = -coss(:)>1 + eps;
% nbroktris = sum(bti1) + sum(bti2);
% coss( bti1 ) = 1;
% coss( bti2 ) = -1;

if bFixRange
    coss = max( min( coss, 1 ), -1 );
end

angles = acos(coss);



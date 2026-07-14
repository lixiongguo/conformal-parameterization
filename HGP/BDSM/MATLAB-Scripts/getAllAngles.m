function [ e_FE ] = getAllAngles( ELS, mesh_FE )

% We populate e_FE by the indices of _FE and the values of _ELS 
e_FE = arrayfun(@(y) ELS(y), mesh_FE);
e_FE = sqrt(e_FE);
% Turn the e_FE to a table so we can use 'rowfun' to calculate the angles
e_FE_table = array2table(e_FE);
e_FE = rowfun(@(a,b,c) extractAngles(a,b,c),e_FE_table);
e_FE = table2array(e_FE);

end


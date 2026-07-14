%% --------------------------------------------------------------------%
%                                                                      %
%                   Create graphs and save data                        %
%                                                                      %
%----------------------------------------------------------------------%

% calc the genus by Euler = 1 -(V-E+F+B)/2
[~,y,~] = networkComponents(mesh_BoundaryAdjacency);
y(y==1)=[];
numOfBoundaries = ~isempty(y) * size(y,1);
genus = 1 - (size(mesh_V,1) - size(mesh_EV,1) + size(mesh_FE,1) + numOfBoundaries)/2;

close all;
figure(1)
ax=gca;
titleString = strcat(meshName, '-   V=', int2str(size(mesh_V,1)), ',  F=', int2str(size(mesh_FE,1)), ',  Genus=', int2str(genus));

set(gcf, 'Position', get(0, 'Screensize'));
hold on


% Plot the difference sum of curvature, from prescribed solved by CVX.
s1 = subplot(2,3,1:2);       
plot(1:numel(ELS_norm), ELS_norm)
grid on;
title('Norm $L_{\infty}$ of curvature on interior vertices from prescribed','Interpreter','latex');
xlabel('Iteration number','Interpreter','latex');
ylabel('Norm','Interpreter','latex');

s2 = subplot(2,3,4:5);       
plot(1:numel(Energy), Energy)
grid on;

title('Energy = $\max \left ( \sigma(a)^2,\frac{1}{\sigma(b)^2}  \right ) + \lambda_l \cdot \left \| J_\kappa ^l - (\overrightarrow{\kappa}^0 - \overrightarrow{\kappa}(e_l)) \right \| _\infty $','Interpreter','latex');
xlabel('Iteration number','Interpreter','latex');
ylabel('Energy','Interpreter','latex');


% 
subplot(2,3,[3 6]);    
% actual_curvature_interior = sum(TEST_Actual_Curvature(final_interiorVertices,:),2) - 2 * pi;
% 
% 
% finalcurv = pi * ones(length(mesh_V),1);
% finalcurv(final_interiorVertices) = 2*pi;
% actual_curvature_interior = sum(TEST_Actual_Curvature(:,:),2) - finalcurv;

plot(1:numel(interiorVertices), K(interiorVertices,end) - target_curvature(interiorVertices));
%plot(1:length(mesh_V),target_curvature(:) - actual_curvature_interior);


grid on;
title('Curvature differece(prescribed/UV)','Interpreter','latex');
xlabel('Vertex number','Interpreter','latex');
ylabel('Curvature difference','Interpreter','latex');

myFig = gcf;

% %%
% figure(2)
% pathToUV = strcat('C:\code\SpaceDeformation2D\results\UVs\', meshName, '.jpg')
% im=imread(pathToUV);
% imshow(im)

%%

pathToSaveVars = strcat('C:\code\SpaceDeformation2D\results\mat-files\final-results\', meshName, '.mat')
pathToSaveFig = strcat('C:\code\SpaceDeformation2D\results\mat-files\final-results\', meshName, '.jpg')
saveas(gcf, pathToSaveFig);
save(pathToSaveVars);


%%

if length(unique(TEST_orientationResults)) == 1
    orientationText = 'Successful';
else
    orientationText = 'Failed';
end

filename = strcat('C:\code\SpaceDeformation2D\results\statistics\', meshName, '.xlsx')


A = {'Mesh Name', meshName;'','';'',''; 'Number of Faces', size(mesh_FE,1); 'Number of Vertices', size(mesh_V,1); 'Genus', genus;'','';'','';
    'Run Time', timer; 'Final Energy', Energy(end); 'Max Distortion', max(tau_sqrt); 'Average Distortion', mean(tau_sqrt);
    'Orientation Test over all triangles', orientationText};
sheet = 1;
xlswrite(filename,A,sheet)


% Change the columns width so its readable
ExcelApp = actxserver('excel.application');
ExcelApp.Visible=1;
NewWorkbook = ExcelApp.Workbooks.Open(filename);   
NewSheet = NewWorkbook.Sheets.Item(1);
NewRange = NewSheet.Range('A1');
NewRange.ColumnWidth = 40;
NewRange = NewSheet.Range('B1');
NewRange.ColumnWidth = 40;

pastePlotToExcel('D24', ExcelApp, myFig);
close (1);
clear myFig;

%%
figure(10);
bins = 10;
a = hist(tau_sqrt, [0.5:0.5:5]);
barh(a,0.5)
histSum = sum(hist(tau_sqrt, bins));
for j = 1:size(a,2)
    text(a(1,j)-0.1,j, cellstr(strcat({'  '} ,num2str(((100*a(1,j))/histSum),'%0.2f%%'))), 'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10);
end
ax = gca;
ax.XLim = ([0 sum(hist(tau_sqrt, bins))]);
ax.XTick = [0 histSum/10 histSum*2/10 histSum*3/10 histSum*4/10 histSum*5/10 histSum*6/10 histSum*7/10 histSum*8/10 histSum*9/10 histSum];
ax.XTickLabel = {'0%','10%','20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'};
ax.YTickLabel = {0.5:0.5:5};
pastePlotToExcel('D1', ExcelApp, gcf);
close(gcf);

figure(9)
caxis manual;
caxis([0 5]);
colorbar;
colormap('jet');
pastePlotToExcel('M2', ExcelApp, gcf);
close(gcf);

%% Close the file and kill the process

NewWorkbook.Save;
NewWorkbook.Close;
ExcelApp.Quit;
ExcelApp.delete;
system('taskkill /F /IM EXCEL.EXE');






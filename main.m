clear; clc; close all
addpath('./functions/') % Load all the function definitions 

%Load the data
load('./Data/growingSphere.mat')

%Plotting parameters
canView1 = [218,13];  fntSz = 30; axLim = 1.1;
% Inputs for OECS computation
t0 = 1; regulFac = 0;

% select data at t0
[~,IdXMesh] = min(abs(timeArr-t0));
v1ct = v{1,IdXMesh}';v2ct = v{2,IdXMesh}';v3ct = v{3,IdXMesh}';
xct = x{IdXMesh}'; yct = y{IdXMesh}'; zct = z{IdXMesh}';
Trianct = TrianT{IdXMesh};

% Compute the strain rate eigenvalues and eigenvectors
[s2,s1,eig2vert,eig1vert] = OECScompute(xct,yct,zct,Trianct,v1ct,v2ct,v3ct,regulFac);

% Plot the quantities
f= figure('color','w','Units','normalized','OuterPosition',[0.0206 0.1481 0.9823 0.7558]);
ax0 = subplot(1,3,1); % The velocity field 
vMag = sqrt(v1ct.^2+v2ct.^2+v3ct.^2);
p = trisurf(Trianct,xct,yct,zct,vMag,'FaceAlpha',1,'Edgecolor','none'); hold on
p.SpecularStrength = 0.2; p.AmbientStrength = 0.3;c = colorbar; colormap(ax0,'turbo');shading interp; view(canView1);camlight
axis equal off; camva(10);c.FontSize = fntSz;c.Location = 'southoutside';
title(sprintf('$$ \\mathbf{v}(\\mathbf{x},%d$$)',t0),'Interpreter','latex','FontSize',fntSz); hold off
xlim([-axLim axLim]); ylim([-axLim axLim]); zlim([-axLim axLim]);

subplot(1,3,2) % The s2 field - repellers
p = trisurf(Trianct,xct,yct,zct,s2,'FaceAlpha',1,'Edgecolor','none'); hold on 
p.SpecularStrength = 0.2; p.AmbientStrength = 0.3;
quiver3(xct,yct,zct,eig2vert(1,:),eig2vert(2,:),eig2vert(3,:),...
    'k','ShowArrowHead','off','LineWidth',1); hold on; shading interp; 
c = colorbar; axis equal off; camva(10);c.FontSize = fntSz;view(canView1);c.Location = 'southoutside';camlight
title(sprintf('$$ s_2(\\mathbf{x},%d),\\mathbf{e_2}(\\mathbf{x},%d)$$',t0,t0),'Interpreter','latex','FontSize',fntSz); hold off
xlim([-axLim axLim]); ylim([-axLim axLim]); zlim([-axLim axLim]);

subplot(1,3,3) % The s1 field - attractors 
p = trisurf(Trianct,xct,yct,zct,s1,'FaceAlpha',1,'Edgecolor','none'); hold on 
p.SpecularStrength = 0.2; p.AmbientStrength = 0.3;
quiver3(xct,yct,zct,eig1vert(1,:),eig1vert(2,:),eig1vert(3,:),...
    'k','ShowArrowHead','off','LineWidth',1); hold on; shading interp; 
c = colorbar; axis equal off; camva(10);c.FontSize = fntSz;view(canView1);c.Location = 'southoutside';camlight
title(sprintf('$$ s_1(\\mathbf{x},%d),\\mathbf{e_1}(\\mathbf{x},%d)$$',t0,t0),'Interpreter','latex','FontSize',fntSz); hold off
xlim([-axLim axLim]); ylim([-axLim axLim]); zlim([-axLim axLim]);

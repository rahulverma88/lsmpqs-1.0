function [k_avg, k_min, k_max] = findCurvatureNearIntfc(data,g,mask)

% function [k_avg k_min k_max] = displayCurvature(data,g,mask,curv,fid,plot)
% data - level set function, will plot curvature values in 'curv' on its level = 0
% g    - grid
% mask - level set function defining pore and grain space
% curv - curvature data
% fid  - (optional) file id of an outfile
%       Outfile, if supplied, will be appended with curvature average, min, max.
% plot - (optional) 1 if you want a figure with curvature plot, defaults to 1
 dist = -2*g.dx(1);
 curv = curvatureSecond(g,data);
 S = isosurface(data, 0);   %s is structure containing vertices and faces of the isosurface
 fld_intfc = S.vertices;
 x_fld = fld_intfc(:,1);
 y_fld = fld_intfc(:,2);  %X and Y might be SWITCHED
 z_fld = fld_intfc(:,3);
 curv(mask > dist) = 0;
 curv_interp = interp3(curv,x_fld,y_fld,z_fld);
 mask_interp = interp3(mask,x_fld,y_fld,z_fld);
 %K = diag(curv_interp);
 %mask_interp = diag(mask_interp);
 % fluid interface are points on boundary of data, that are not on
 % grain(mask) boundary
 
 I =  (mask_interp < dist) & (curv_interp > 0) ;
 K1 = curv_interp(I);
 
 k_avg = mean(K1) ;
 k_min = min(K1) ;
 k_max = max(K1);

 %}
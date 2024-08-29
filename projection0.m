function proj2d = projection0(data3d, param, iview)

angle_rad = param.deg(iview)/360*2*pi;
proj2d = zeros(param.nu, param.nv, 'single');
lineLength = 4; % source length in mm
nLinePoints = 10; % number of points along the line source
linePositions = linspace(-lineLength/2, lineLength/2, nLinePoints);

for ip = 1:nLinePoints
    sourceZ = param.sourcez + linePositions(ip);
    sourceY = param.sourcey;
    sourceX = param.sourcex;
    [vv, uu] = meshgrid(-param.vs, -param.us);
    [xx, yy] = meshgrid(param.xs, param.ys);
    rx = (((xx.*cos(angle_rad) - yy.*sin(angle_rad)) - xx(1,1))/param.dx + 1);
    ry = (((xx.*sin(angle_rad) + yy.*cos(angle_rad)) - yy(1,1))/param.dy + 1);
    
    temp3D = zeros(size(data3d), 'like', data3d);  
        for iz = 1:param.nz   
            tempSlice = interp2(data3d(:,:,iz), rx, ry, param.interptype);
            tempSlice(isnan(tempSlice)) = 0;
            temp3D(:,:,iz) = tempSlice;  
        end   
    temp3D = permute(temp3D,[1 3 2]);
    [yy, zz] = meshgrid(param.ys, param.zs);

    for ix = round(param.nx/2 + sourceX/param.dx):round(param.DSD/param.dx)
        Ratio = (param.xs(ix) - sourceX) / (param.DSD - sourceX);
     
    
        geomratio =1;
        shifted_yy = (yy - sourceY / geomratio);
        shifted_zz = (zz - sourceZ / geomratio);
        pu = uu * Ratio;
        pv = vv * Ratio;    
        pu = (pu - shifted_yy(1,1)) / param.dy + 1; 
        pv = (pv - shifted_zz(1,1)) / param.dz + 1; 
        tmp = interp2(single(temp3D(:,:,ix)), single(pv), single(pu), param.interptype);      
        tmp(isnan(tmp)) = 0;   
        proj2d = proj2d + tmp';
    end
end

proj2d = proj2d / nLinePoints; % Normalize the projection
end
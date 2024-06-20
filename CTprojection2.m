function [ proj ] = CTprojection2( img, param )
%CTPROJECTION Summary of this function goes here
%   Detailed explanation goes here

proj = zeros(param.nu, param.nv, param.nProj,'single');

for i = 1:param.nProj    
    proj(:,:,i) = projection2(img,param,i);
end
    
end

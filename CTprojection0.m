function [ proj ] = CTprojection0( img, param )
%CTPROJECTION Summary of this function goes here
%   Detailed explanation goes here

proj = zeros(param.nu, param.nv, param.nProj,'single');

for i = 1:param.nProj    
    proj(:,:,i) = projection0(img,param,i);
end
    
end

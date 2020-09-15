function [FM]=Flow_Routing_Step2(E,T,Qrdt) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   T = flow_matrix(E, R) computes a sparse linear system representing flow from
%   pixel to pixel in the DEM represented by the matrix of height values, E.  R
%   is the matrix of pixel flow directions as computed by dem_flow.  
%   T is  numel(E)-by-numel(E).  The value T(k,l) is the negative of the fractional
%   flow from pixel l to pixel k, where pixels are numbered columnwise. For
%   example, if E is 15-by-10, then T is 150-by-150, and T(17,18) is the
%   negative of the fractional flow from pixel 18 (row 3, column 2) to pixel 17
%   (row 2, column 2).
%   [i,j] = ind2sub([m,n],l); flow into --->  [i,j] = ind2sub([m,n],k); 
Qrdt = reshape(Qrdt,numel(E),1); 
[m,n]=size(E); 
%T = (-T + speye(m*n,n*m)); %%% Flow Contribution to that pixel  
%FM = sum(T,2); %%%%% Sum of contribution from upslope pixels  
%FM=full(FM);
FM=(-T + speye(m*n,n*m))*Qrdt; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FM = reshape(FM, size(E));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 

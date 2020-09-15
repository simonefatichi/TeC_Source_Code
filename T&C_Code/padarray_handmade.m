function [Enl_Row]=padarray_handmade(M,srow,scol)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(M);
Enl_Col=[M(:,1)*ones(1,scol), M, M(:,n)*ones(1,scol)]; 
Enl_Row = [ones(srow,1)*Enl_Col(1,:) ; Enl_Col ; ones(srow,1)*Enl_Col(end,:) ];
end 
%%%%%%%%%%%%%%%
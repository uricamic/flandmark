function [ res ] = crossprod( a, b )
%CROSSPROD Summary of this function goes here
%   Detailed explanation goes here

    Ax = [  0,   -a(3),   a(2); 
           a(3),     0,  -a(1);
          -a(2),  a(1),      0;];
      
   res = Ax * b(:);

end


function [ outside_option ] = outside(z)
% Worker's outside option

global b0 b1

outside_option = b0 * z^b1;

end


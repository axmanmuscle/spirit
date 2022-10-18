function out = spirit_conv(array, weights)
%fill_points helper function to fill in missing k-space points
% This function fills in the missing k-space data for a given kernel
% We will always be starting from the first point such that the kernel fits
% fully
%
% i.e. in a 5x5 grid with a 3x3 kernel, we start at [2 2]:
% x x x o o
% x x x o o
% x x x o o
% o o o o o
% o o o o o
%
% Similarly, we only go to points such that the kernel fits in fully
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%     array: the full 3D (ny x nx x ncoils) kspace data
%
%     kernel size: [ky kx] size of the SPIRiT kernel
%
%     coil: integer, which coil number we're working on
%
%     weights: the previously solved-for interpolation weights for the
%     kernel
%
% *********************
%   Output Variables:
% ********************* 
%
%    out: a 2D (ny x nx) array of reconstructed k-space data, only for the
%    specified coil

nCoils = size(array, 3);
out = 0*array;

for coil = 1:nCoils
  wi = weights(:, :, :, coil);
  for getCoil = 1:nCoils
    wij = wi(:, :, getCoil);
    res = filter2(wij, array(:, :, getCoil), 'same');
    out(:, :, coil) = out(:, :, coil) + res;
  end
end

end
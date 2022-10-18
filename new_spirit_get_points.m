function out = new_spirit_get_points(pt_idx, array, kernel_sz)
%get_points_spirit helper function to retrieve correct points
% this function is used both to set up the weights and get the appropriate
% points to interpolate with when filling in k-space
% slightly different than the GRAPPA version because of the points we need
% to grab
%
% What about if we know that the kernel is square and just need the size
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%     pt_idx: A 3 element vector containing the [row column coilidx] of a specific
%     point. Coil index matters here since for solving for the weights, we
%     use the point at [row column] from the other coils.
%         NOTE: this is different than grappa
%
%     array: the full 3D (ny x nx x ncoils) kspace data
%
%     kernel_sz: a 2 element array containing the [row column] size of the
%     kernel. for spirit we'll always use this kernel completely filled in
%     except for the center point
%
% *********************
%   Output Variables:
% ********************* 
%
%    out: the set of points for the input pt_idx over all the coils

ncoils = size(array, 3);
py = pt_idx(1);
px = pt_idx(2);
pc = pt_idx(3);

k_y = kernel_sz(1);
k_x = kernel_sz(2);

kdy = (k_y - 1) / 2;
kdx = (k_x - 1) / 2;

kcy = (k_y + 1) / 2;
kcx = (k_x + 1) / 2;

ypts = py-kdy:py+kdy;
xpts = px-kdx:px+kdx;

pts = array(ypts, xpts, :);

out = pts(:);


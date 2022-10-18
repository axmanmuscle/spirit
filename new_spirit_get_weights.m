function [W, A, B] = new_spirit_get_weights(acr, kernel_sz)
%new_spirit_get_weights lets try this again

% auto-calibration region dimensions
szAcrY = size(acr, 1);
szAcrX = size(acr, 2);
ncoils = size(acr, 3);

% kernel dimensions
szKerY = kernel_sz(1);
szKerX = kernel_sz(2);

n = szKerY*szKerX;

% centering
center_y = (szKerY - 1)/2;
center_x = (szKerX - 1)/2;

% number of sliding kernel fits in the ACR
numfits_x = szAcrX - szKerX + 1;
numfits_y = szAcrY - szKerY + 1;

npoints = szKerY*szKerX*ncoils;
W = zeros([kernel_sz ncoils ncoils]);
A = zeros([numfits_x*numfits_y npoints]);

for coilIdx = 1:ncoils
  B = zeros( [numfits_x*numfits_y 1] );
  b_idx = 1;
  for x = 1:numfits_x
    for y = 1:numfits_y
      B(b_idx) = acr( y + center_y, x + center_x, coilIdx);
      pts = new_spirit_get_points([y+center_y x+center_x coilIdx], acr, kernel_sz);
      A(b_idx, :) = pts;
      b_idx = b_idx + 1;
    end
  end

  dk = zeros([szKerY szKerX ncoils]); 
  dk((end+1)/2, (end+1)/2, coilIdx) = 1;
  idxY = find(dk);
  smp = ones(size(dk));
  smp(idxY) = 0;
  idxA = find(smp);

  A2 = A(:, idxA);
  rk = A2 \ B;
  
  ker = zeros(size(dk));
  ker(idxA) = rk;

  W(:, :, :, coilIdx) = ker;


end
end
function run_spirit()
addpath('../grappa/');
load("data.mat", 'd1');

acr_sz = [31 31];
kernel_sz = [3 3];

acr = get_acr(d1, acr_sz);
weights = new_spirit_get_weights(acr, kernel_sz);

sm = grappa_samplingmask(size(d1), acr_sz, 2, 'both');
d2 = bsxfun( @times, d1, sm);

x0 = grappa(d2, kernel_sz, acr_sz);

idx_acq = d2~=0;
idx_nacq = d2 == 0;

check_adjoints = 0;
if check_adjoints
  [o, e] = checkAdjoint(d2, @applyG);
  if o ~= 1
    error('G adjoint wrong');
  end

  apA = @( in, op ) applyA(in, size(d2), idx_nacq, op);

  [o2, e2] = checkAdjoint(d2(idx_nacq), apA);
  if o2 ~= 1
    error('A adjoint wrong');
  end
end

% set up b for Ax = b
b = applyG(d2, 'notransp') - d2;
b = -b(:);

% initial guess of GRAPPA?
x_nacq = lsqr(@applyA, b, 1e-6, 20, [], [], x0(idx_nacq), size(d2), idx_nacq);
x_out = d2;
x_out(idx_nacq) = x_nacq;

figure; imshowscale(mri_reconSSQ(x_out));

check_lustig = 1;

if check_lustig
  addpath(genpath('../ESPIRiT'));
  kernel_test = calibSPIRiT(acr, kernel_sz, size(d1, 3), 0);

  err_ker = norm(kernel_test(:) - weights(:));
  if abs(err_ker) > 1e-6
    error('Difference in kernel weights');
  end

  GOP = SPIRiT(kernel_test, 'conv', [size(d1, 1), size(d1, 2)]);

  %%% GOP * data = (G - I)data
  %%% GOP' * data = (G - I)* data
  oi = spirit_conv(d2, weights);
  oadj = spirit_conv_adj(d2, weights);
  test1 = GOP * d2;
  test2 = oi - d2;
  err_recon = norm(test1(:) - test2(:)) / norm(test2(:));
  if abs(err_recon) > 1e-6
    error('Difference in reconstruction (G-I)x');
  end

  test3 = GOP' * d2;
  test4 = oadj - d2;
  err_recon_adj = norm(test3(:) - test4(:)) / norm(test4(:));
  if abs(err_recon_adj) > 1e-6
    error('Difference in reconstruction (G-I)*x');
  end

  [res, ~] = cgSPIRiT(d2, GOP, 20, 0, x0);
  err_total = norm(res(:) - x_out(:))/norm(x_out(:));

  figure; imshowscale(mri_reconSSQ(res));
  if abs(err_total) > 1e-6
    error('Difference in total reconstruction after LSQR.');
  end
end

  function out = applyG(in, op)
    if strcmp(op, 'transp')
      out = spirit_conv_adj(in, weights);
    else
      out = spirit_conv(in, weights);
    end
  end

  function out = applyA(in, array_sz, idx_nacq, op)
    if strcmp(op, 'transp')
      tmp = reshape(in, array_sz);
      out = applyG(tmp, op) - tmp;
      out = out(idx_nacq);
    else
      tmp = zeros(array_sz);
      tmp(idx_nacq) = in;
      out = applyG(tmp, op) - tmp;
      out = out(:);
    end
  end

end
function run_spirit_fista()
addpath('../grappa/');
load("data.mat", 'd1');

d1 = d1./norm(d1(:));
acr_sz = [31 31];
kernel_sz = [3 3];

acr = get_acr(d1, acr_sz);
weights = new_spirit_get_weights(acr, kernel_sz);

sm = grappa_samplingmask(size(d1), acr_sz, 2, 'horiz');
d2 = bsxfun( @times, d1, sm);

x0 = grappa(d2, kernel_sz, acr_sz);

idx_acq = d2~=0;
idx_nacq = d2 == 0;

y_data = d2(idx_acq);
eps = 1e-2;

checkProx = 1;
if checkProx
  proxTest = @(in, sc) projB0(in, sqrt(eps));
  test1 = proxh(x0, 0);
  test2 = proxAffine(proxTest, x0, @applyD, -1*y_data, 1);
  
  err = norm(test1(:) - test2(:))/norm(test2(:));
  if abs(err) > 1e-6
    error('Something wrong with proximal operator');
  end
end

%[xStar, objVal, relDiffs] = fista(0*x0(:), @grad_g, @proxh, 'g', @normG, 'h', @h, 'N', 100, 'verbose', true, 't', 0.05 );
[xStar, objVal, relDiffs] = fista_wLS(x0(:), @normG, @grad_g, @proxh, 'h', @h, 'N', 100, 'verbose', true, 't0', 0.01);

figure; imshowscale(mri_reconSSQ(d1));
figure; imshowscale(mri_reconSSQ( reshape(xStar, size(d1))));

x3 = norm(xStar(idx_acq) - y_data);
disp(x3)

  function out = proxh(in, sc)
    din = in(idx_acq) - y_data;
    tmp_proj = projB0(din, sqrt(eps));
    tmp2 = din - tmp_proj;
    out = zeros(size(d2));
    out(idx_acq) = tmp2;
    out = in(:) - out(:);
  end

  function out = applyG(in, op)
    in = reshape(in, size(d2));
  if strcmp(op, 'transp')
    out = spirit_conv_adj(in, weights);
  else
    out = spirit_conv(in, weights);
  end
  out = out(:);
  end

  function out = h(x)
    out = ind(x, y_data, sqrt(eps));
  end

  function out = ind(x, y, r)
    dx = x(idx_acq);
    if norm(dx - y) < r
      out = 0;
    else
      out = Inf;
    end
  end

  function out = projB(x, y, r)
    dxy = x - y;
    if norm(dxy) < r
      out = x;
    else
      out = dxy * (r / norm(dxy)) + y;
    end
  end

function out = projB0(x, r)
    if norm(x) < r
      out = x;
    else
      out = x * (r / norm(x));
    end
  end

  function out = normG(in)
    tmpg = applyG(in, 'notransp');
    out = 0.5 * norm(tmpg(:) - in(:))^2;
  end

  function out = grad_g(in)
    tmpgrad = applyG(in, 'notransp') - in;
    out = applyG(tmpgrad, 'transp') - tmpgrad;
    out = out(:);
  end

  function out = applyD(in, op)
    if strcmp(op, 'transp')
      out = 0*d2;
      out(idx_acq) = in;
    else
      out = in(idx_acq);
    end
  end

end
function out = proxAffine( prox, x, applyA, b, alpha )
  out = x - alpha * applyA( applyA( x, 'notransp' ) + b - ...
    prox( applyA( x, 'notransp' ) + b, 1/alpha ), 'transp' );
end
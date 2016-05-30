%isabel_data must be transposed
isabel_x = hdf5read('isabel_2d.h5','/Velocity/X-comp')';
isabel_y = hdf5read('isabel_2d.h5','/Velocity/Y-comp')';
isabel_fld = zeros(2, size(isabel_x, 1), size(isabel_x, 2));
isabel_fld(1,:,:) = isabel_x;
isabel_fld(2,:,:) = isabel_y;

metsim1_x = hdf5read('metsim1_2d.h5','/Velocity/X-comp');
metsim1_y = hdf5read('metsim1_2d.h5','/Velocity/Y-comp');
metsim1_fld = zeros(2, size(metsim1_x, 1), size(metsim1_x, 2));
metsim1_fld(1,:,:) = metsim1_x;
metsim1_fld(2,:,:) = metsim1_y;
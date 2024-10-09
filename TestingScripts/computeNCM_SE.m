%Compute strain energy and NCM from 3D traction field

testfield = force_vector{5};

%% Strain Energy

%SE is computed as 1/2 * integral(T(r)*u(r) dxdydz)




%% Net Contractile Moment

% Computed as a function of centroid, positions, and forces
% Mij = integral(ri * Fj(r) dr) where ri is cartestion component of
% position vector with respect to center of mass and Fj are cartesian
% coordinates of the force (traction * facet area) at a facet. 
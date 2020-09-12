function [ KE, PE ] = getEnergy(  pos, vel, mass, G )
%GETENERGY Get kinetic energy (KE) and potential energy (PE) of simulation
%   pos is N x 3 matrix of positions
%   vel is N x 3 matrix of velocities
%   mass is an N x 1 vector of masses
%   G is Newton's Gravitational constant
%   KE is the kinetic energy of the system
%   PE is the potential energy of the system

% Kinetic Energy:
KE = 0.5 * sum(sum( mass .* vel.^2 ));


% Potential Energy:

% positions r = [x,y,z] for all particles
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);

% matrix that stores all pairwise particle separations: r_j - r_i
dx = x' - x;
dy = y' - y;
dz = z' - z;

% matrix that stores r for all particle pairwise particle separations 
r = sqrt(dx.^2 + dy.^2 + dz.^2);

% sum over upper triangle, to count each interaction only once
PE = G *  sum(sum(triu(-(mass*mass')./r,1)));

end


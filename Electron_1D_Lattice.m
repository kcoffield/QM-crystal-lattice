% The quantum mechanical problem of the motion of an electron through a lattice
% of identical atoms. This program finds the matrix H (the Hamiltonian operator)
% so that the eigen values (=allowed energies) may be found from the time-
% independent Schroedinger equation HPhi=EPhi. Now, H = -((hbar)^2)/(2mass)*(del^2)
% + V, where hbar is the quantum constant, mass is the mass of the electron, del is
% the delta operator, and V is the potential energy function that determines the
% motion of the partical. This program approximates the second derivative (del^2)Phi
% as the difference (1/d^2)[Phi(x+d)-2Phi(x)+Phi(x-d)] and this formula is used
% recursively in the first matrix calculated in this code, 'A'. The 2nd matrix
% calculated in this code, 'B' represents the 2nd term in the Hamiltonian, or V.
% if we let c = (hbar^2)/(2*mass), then the working equation for this program is 
% (1/c)*H = (1/d^2)*[-Phi(n-1)+2*Phi(n)-Phi(n+1)]+(1/c)V(n)Phi(n) = [(1/d^2)*A + B]Phi.

k = 10;			% number of atoms
hbar = 1.05457e-34;   mass = 9.10938e-31;   c = (hbar^2)/(2*mass);
m = 3e-10;		% distance between atoms
s = 20;			% subintervals between atoms
dx = m/s;		% discretize the distance m
M = k*m/dx;		% size of MxM matrix
L = k*m;		% total length of crystal
x = 0:dx:L;		% Parse the interval [0,L]
xx = mod(x,m);	% due to periodicity of V(x)...
V0 = 100*hbar;	% The depth of the potential wells
V = @(x) V0*cos(x*2*pi/m); % The potential energy function
A = zeros(M,M);	% Define the matrix that will store kinetic part of Hamiltonian
	


% This first matrix 'A' will represent the first term in the Hamiltonian, or the 
% second spatial derivative term.
	% Define the size of the matrix 'A'
		A = zeros(M,M);
	% Boundary Conditions + something else???
		A(1,1) = 2;  A(M,M)   = 2;
		A(1,2) = -1; A(M,M-1) = -1;
		A(1,M) = -1; A(M,1)   = -1;
	% Loop through A for difference rule for 2nd derivative
		for j = 2:M-1
		  for n = 1:M
		    A(j,j-1) = -1;
		    A(j,j)   = 2;
		    A(j,j+1) = -1;
		  end
		end



% This second matrix 'B' will represent the second term in the Hamiltonian, or
% the potential, V.
	% Assign the size of matrix 'B'
		B = zeros(M,M);
	% Loop through B for potential values
		for j = 1:M
			for n = 1:M
				xxj = xx(j);
		 		B(j,j) = V(xxj);
		 	end
		end



% Scaled by c, the Hamiltonian for the time-independent Schroedinger equation,
% H*phi = E*phi is
H = c*(1/(dx^2))*A + B;
% where E will be the eigenvalues. Plot these values and compare to V(max) and
% to theory. For the current configuration (23 Nov 2017), V(max) = 0.

E = eig(H);
xx = xx(1:length(xx)-1);
stairs(xx,E)
hold on
V_max = max(V(x));   V_max(1:length(xx)) = V_max;
plot(xx,V_max)
% FEM1D_abar()
% We consider the FE approximation of 
%  -(a(x) u'(x))' = f(x) a.s., for x in [0,1],
% subject to u(0) = g a.s. and a(1)*u'(1) = h a.s. Here a(x) =a(w,x) is a 
% random field where the average (integral) over each element is supplied. 
% The force term f(x) is constant for x in [0,1]. 
%
% We use piecewise linear elements on numFemNodes (given) nodes and 
% numFemNodes-1 elements to approximate the solution. Quadrature of a(x) is
% carried out over a finer mesh using numQuadPts (given) points per element 
% (inclusive of the element's end points).
%
% INPUT:
%  numFemNodes   - [number] number of grid points for the coarse grid 
%                  (inclusive of interval end points, i.e. [0,1] by n grid 
%                  points gives coarse mesh size of dx = 1/(n-1))
%  abar          - [numFemNodes-1 x 1 vector] spatial average on each elem
%  force         - [numFemNodes x 1 vector] of forces
%  condg         - [number] Dirichlet (essential) boundary condition at x=0
%  condh         - [number] Neumann (natural) boundary condition at x=1
%
% OUTPUT:
%  U             - [numFemNodes x 1 vector] of displacements

function [U] = FEM1D_abar(numFemNodes, abar, force, condg, condh)                         
% *** INPUT ***
%
% FE MESH PROPERTIES
%  nn         - number of nodes (= number of points of the coarse mesh)
%  ne         - number of elements (= nn - 1)
nn = numFemNodes;
ne = numFemNodes-1;
%
% *** PRE-PROCESSING ***
%
%
% ALLOCATE DATE STRUCTURES
%  U         - nn x 1 vector of displacements (the approximate solution)
U = zeros(nn,1);
%
% GENERATE FINITE ELEMENT MESH
% The following structures are used to describe the FE discretization. The
% finite element approximation is carried out on a coarse mesh [node,
% element] and a finer mesh [nodef, connectf] is used for quadrature.
%  node      - matrix of node coordinates (local), i.e. node(i) = x_i 
%  element   - is the matrix of node connectivity, i.e. element(i,:) = e_i,
%              a vector of the nodes (in global coord:s) that comprise 
%              element i
%  dir       - global coord:s relative to the coarse mesh [node, element] 
%              of Dirichlet nodes (here, chosen to be the left end point)
%  free      - active degrees of freedom (n in node not in dir)
[~, elem] = generateMesh(0,1,nn,ne);
dir = elem(1,1);                          
free = setdiff((1:nn)',dir);            
%
% *** PROCESSING ***
%
%  K          - global stiffness matrix
%  D          - unknown displacement vector
%  F          - global force vector
%  Fe,Fn,Fd   - The global force vector F = Fe + Fn - Fd where Fe is the
%               element force vector, Fn is the force vector associated to
%               the natural boundary condition, Fd is the force vector
%               associated to the essential boundary condition.
%
% COMPUTE DISPLACEMENTS U
%  wn      - natural boundary condition at the right hand end point
wn = zeros(nn, 1); wn(end) = 1;
% FN can be computed independently of the choice of shape of basis function
Fh = wn*condh;
% Assemble K, Fe, and Fd
[K, Fe, Fg] = assembleGlobal(elem, nn, abar, force, condg);
% Assemble F
F = Fe + Fh - Fg;
% Compute displacements (suffices to solve for independent nodes)
D = K(free,free)\F(free);
% Save displacements for each sample
U(free,1) = D;
% Apply Dirichlet/essential boundary condition
%U(dir,1) = g; % this is taken care of during assembleGlobal() by Fd
end
%
% GENERATE MESH
% Generate vector of node coordinates and a matrix of node connectivity for
% an interval [l,L] having numNodes nodes and numElements elements.
% Supply:
%  l             - left hand end point 
%  L             - right hand end point
%  numNodes      - number of nodes on [l,L]
%  numElements   - number of elements on [l,L]
function [node, connectivity] = generateMesh(l,L,numNodes,numElements)
nodesPerElement = ceil(numNodes/numElements);
node = linspace(l, L, numNodes)';
step = (1:floor(numNodes/numElements):floor(numNodes/numElements)*numElements)';
connectivity = zeros(numElements,nodesPerElement);
connectivity(:,1) = step;
for i = 1:nodesPerElement-1
   connectivity(:,i+1) = step + i.*ones(numElements,1);
end
end
%
% ASSEMBLE GLOBAL
% Assemble global stiffness matrix, K, and distributed forces vector, Fe
function [K, Fe, Fd] = assembleGlobal(elem, nn, a, f, g)
% COMPUTE ELEMENTWISE: K, Fe and Fd
%  B1     - change of coordinates matrix 
%  B2     - assembling Fe (assume f(x) is constant on each element)
%  A      - length of element (same for each uniform element)
ne = size(elem,1);
K = sparse(nn,nn);
Fe = zeros(nn,1);
B1 = [1,-1;-1,1];
B2 = [2,1;1,2];
A = 1/(nn-1); % H
wd = zeros(nn,1); wd(1) = -1;
Fd = wd*g*a(1)/A; % ei = 1
for ei = 1:ne % loop over all elements; we only solve for free nodes
    %  dof    - global degrees of freedom map
    dof = elem(ei,:);
    % Assemble local stiffness matrix sparse(r,s,t,nn,nn)
    C1 = B1*a(ei)/A;    C2 = [dof;dof]; 
    r = C2(:);    s = [dof,dof]';    t = C1(:);
    % Assemble K, globabl stiffness matrix
    K = K + sparse(r,s,t,nn,nn); % OR: K(dof,dof) = K(dof,dof) + B1*(1/A)*a(ei);
    % Assemble Fe; if 1st element then must apply Dirichlet condition
    Fe(dof) = Fe(dof) + B2*[f(elem(ei,1));f(elem(ei,2))].*(A/6); 
end
end
%%% END OF FILE
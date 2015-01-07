% This script demonstrates the manual computation of integrals of shape
% functions and their derivatives over mesh elements and compares the
% results with toast's built-in integration functions.


%% Load a mesh consisting of linear triangles
mesh = toastMesh('../../test/2D/meshes/circle25_32.msh');

% Pick an element to work on and extract its geometry
el = mesh.Element(10);
[vtx,idx,eltp] = el.Data();

% Compute natural coordinates
a1 = vtx(2,1)*vtx(3,2) - vtx(3,1)*vtx(2,2);
a2 = vtx(3,1)*vtx(1,2) - vtx(1,1)*vtx(3,2);
a3 = vtx(1,1)*vtx(2,2) - vtx(2,1)*vtx(1,2);
b1 = vtx(2,2) - vtx(3,2);
b2 = vtx(3,2) - vtx(1,2);
b3 = vtx(1,2) - vtx(2,2);
c1 = vtx(3,1) - vtx(2,1);
c2 = vtx(1,1) - vtx(3,1);
c3 = vtx(2,1) - vtx(1,1);

% Compute element size
elsize = 0.5 * (a1+a2+a3);

% This can also be done with a built-in function:
elsize_toast = el.Size();
err = abs(elsize-elsize_toast);


%% Test 1:
% Let's compute the integral of the shape functions associated
% with each of the 3 element nodes over the element, utilising an
% analytic integration rule.
intf_analytic = elsize * [1/3; 1/3; 1/3];

% Now do the same with a numerical quadrature rule
qwght = 1/6 * [1 1 1];
qabsc = [1/2,0; 0,1/2; 1/2,1/2]';
fun = el.ShapeFun(qabsc);
der = el.ShapeDer(qabsc);
intf_numeric = zeros(3,1);
for q=1:3
    jac = der(:,:,q)*vtx;
    v = qwght(q) * det(jac);
    intf_numeric = intf_numeric + v*fun(:,q);
end

% Now compare the results with toast's built-in function
intf_toast = el.Mat('F');
err_analytic = norm(intf_analytic-intf_toast);
err_numeric  = norm(intf_numeric -intf_toast);
fprintf ('Test of IntF (integral of shape function over element):\n');
fprintf ('Discrepancy between toast and analytic integral: %e\n', err_analytic);
fprintf ('Discrepancy between toast and numerical integral: %e\n', err_numeric);


%% Test 2:
%Next, the integral of the product of two shape functions
intff_analytic = elsize/12 * [2,1,1; 1,2,1; 1,1,2];

% And the numerical quadrature version
intff_numeric = zeros(3,3);
for q=1:3
    jac = der(:,:,q)*vtx;
    v = qwght(q) * det(jac);
    intff_numeric = intff_numeric + v * (fun(:,q) * fun(:,q)');
end

% And compare with the toast function
intff_toast = el.Mat('FF');
err_analytic = norm(intff_analytic-intff_toast);
err_numeric  = norm(intff_numeric -intff_toast);
fprintf ('\nTest of IntFF (integral of product of shape functions):\n');
fprintf ('Discrepancy between toast and analytic integral: %e\n', err_analytic);
fprintf ('Discrepancy between toast and numerical integral: %e\n', err_numeric);


%% Test 3:
% Next, the integral of the product of two shape function derivatives
intdd_analytic = 1/(4*elsize) * [b1^2  + c1^2,  b1*b2 + c1*c2, b1*b3 + c1*c3; ...
                        b1*b2 + c1*c2, b2^2  + c2^2,  b2*b3 + c2*c3; ...
                        b1*b3 + c1*c3, b2*b3 + c2*c3, b3^2  + c3^2 ];
                    
% And the numerical quadrature version
intdd_numeric = zeros(3,3);
for q=1:3
    jac = der(:,:,q)*vtx;
    ijac = inv(jac);
    v = qwght(q) * det(jac);
    ider = ijac*der(:,:,q);
    intdd_numeric = intdd_numeric + v * (ider' * ider);
end

% Now compare the results with toast's built-in function
intdd_toast = el.Mat('DD');
err_analytic = norm(intdd_analytic-intdd_toast);
err_numeric  = norm(intdd_numeric -intdd_toast);
fprintf ('\nTest of IntDD (integral of product of shape function derivatives):\n');
fprintf ('Discrepancy between toast and analytic integral: %e\n', err_analytic);
fprintf ('Discrepancy between toast and numerical integral: %e\n', err_numeric);

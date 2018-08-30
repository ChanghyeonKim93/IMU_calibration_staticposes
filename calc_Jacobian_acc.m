function J = calc_Jacobian_acc(acc, theta)
global grav

M = size(acc,2);
J = zeros(M,9);

a_yz = theta(1);
a_zy = theta(2);
a_zx = theta(3);
sx = theta(4);
sy = theta(5);
sz = theta(6);
bx = theta(7);
by = theta(8);
bz = theta(9);

for i = 1:M
    axs = acc(1,i);
    ays = acc(2,i);
    azs = acc(3,i);
    
    ax = sx*(axs+bx) - sy*a_yz*(ays+by) + sz*a_zy*(azs+bz);
    ay = sy*(ays+by) - sz*a_zx*(azs+bz);
    az = sz*(azs+bz);
    
    J(i,1) = -2*ax*(-sy*(ays+by));
    J(i,2) = -2*ax*( sz*(azs+bz));
    J(i,3) = -2*ay*(-sz*(azs+bz));
    
    J(i,4) = -2*ax*( axs+bx);
    J(i,5) = -2*ax*(-a_yz*(ays+by))-2*ay*(ays+by);
    J(i,6) = -2*ax*( a_zy*(azs+bz))-2*ay*(a_zx*(azs+bz))-2*az*(azs+bz);
    
    J(i,7) = -2*ax*sx;
    J(i,8) = -2*ax*(-sy*a_yz)-2*ay*sy;
    J(i,9) = -2*ax*( sz*a_zy)-2*ay*(-sz*a_zx)-2*az*(sz);
end

end
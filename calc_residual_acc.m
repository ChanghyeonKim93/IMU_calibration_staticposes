function r = calc_residual_acc(acc, theta)
global grav

M = size(acc,2);
r = zeros(M,1);

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
    r(i) = grav*grav - ax*ax - ay*ay - az*az;
end

end
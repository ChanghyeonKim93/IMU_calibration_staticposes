function R = quat_int_from_zero_RK4(w, t)

len = length(t);
if(length(w) ~= length(t))
   assert('Difference length of inputs !\n'); 
end

q = zeros(4,len);
q(:,1) = [1,0,0,0].';

for i = 1:len-1
    dt = t(i+1) - t(i);
    w_tmp = w(:,i);
    k1 = quat_derivative_kch(q(:,i),         w_tmp);
    k2 = quat_derivative_kch(q(:,i)+k1*dt/2, w_tmp);
    k3 = quat_derivative_kch(q(:,i)+k2*dt/2, w_tmp);
    k4 = quat_derivative_kch(q(:,i)+k3*dt,   w_tmp);
    
    q(:,i+1) = q(:,i) + 1/6*dt*(k1 + 2*k2 + 2*k3 + k4);    
    q(:,i+1) = q(:,i+1)/norm(q(:,i+1));

end

R = quat_to_dcm_kch(q(:,end));

end
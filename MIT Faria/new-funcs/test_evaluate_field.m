function test_evaluate_field(p)


temp = fft2(1./(1+p.xx.*p.xx + p.yy.*p.yy));
temp = temp + fft2(1./(1+(p.xx+p.Lx).*(p.xx+p.Lx) + p.yy.*p.yy));
temp = temp + fft2(1./(1+(p.xx-p.Lx).*(p.xx-p.Lx) + p.yy.*p.yy));
temp = temp + fft2(1./(1+p.xx.*p.xx + (p.yy+p.Ly).*(p.yy+p.Ly)));
temp = temp + fft2(1./(1+p.xx.*p.xx + (p.yy-p.Ly).*(p.yy-p.Ly)));

function outthing = evalthing(x,y,p)
    outthing = 1/(1+x^2+y^2) + 1/(1+(x+p.Lx)^2+y^2) + 1/(1+(x-p.Lx)^2+y^2) + + 1/(1+x^2+(y + p.Ly)^2) + 1/(1+x^2+(y - p.Ly)^2);
end


error = zeros(size(-20:0.1:20));
error1 = zeros(size(-20:0.1:20));
i = 1;
yval = 0;
for xval = -20:0.1:20
    error(i) = abs(evaluate_field(temp, yval, xval, p) - evalthing(xval,yval,p));
    error1(i) = abs(evaluate_field_old(temp, yval, xval, p) - evalthing(xval,yval,p));
    i = i + 1;
end

plot(-20:0.1:20, error, 'g')
hold on;
plot(-20:0.1:20, error1,'r')
end


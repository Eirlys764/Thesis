function Xn = side(gr, gFov1, gFov2, n)

ss = zeros(1,n);
for i=1:n
    sig(i) = gr(:,i)'*(gFov2-gFov1);
    if sig(i) == 0 % bisection
        ss(i)=0;
    elseif sig(i) < 0 %Fov1
        ss(i)=1;
    elseif sig(i) > 0 %Fov2
        ss(i)=2;
    end
end

cs0 = countj(0, ss);
cs1 = countj(1, ss);
cs2 = countj(2, ss);
Xn =-1;
if 2 <= max(cs1, cs2) && max(cs1, cs2)==n-cs0 && n-cs0 <=n
    Xn = 1;
elseif n-1 <= cs0 <= n
    Xn = 0;
end

end

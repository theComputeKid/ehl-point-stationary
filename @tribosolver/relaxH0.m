function relaxH0(obj,l)

Lk = obj.Levels(l);
p = Lk.Results.p;
dx = Lk.Domain.dx;
dy = Lk.Domain.dy;
fb = Lk.fb;
relH0 = obj.Relaxation.h0;
h0 = Lk.Results.h0;

dF = gather(sum(p .* dx .* dy, "all") + fb);
h0 = h0 + relH0*dF;

for i = 1:obj.Domain.mgl
    obj.Levels(i).Results.h0 = h0;
    obj.Levels(i).calcDeformation;
end

end
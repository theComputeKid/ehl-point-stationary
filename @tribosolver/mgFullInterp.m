function mgFullInterp(obj,l)

Lk = obj.Levels(l);
Lkc = obj.Levels(l-1);

x = Lk.Domain.x;
xc = Lkc.Domain.x;

y = Lk.Domain.y;
yc = Lkc.Domain.y;

pc = Lkc.Results.p;

obj.Levels(l).Results.p = (interp2(xc',yc',pc',x',y',"makima"))';
obj.Levels(l).Results.p(obj.Levels(l).Results.p < 0) = 0;

end
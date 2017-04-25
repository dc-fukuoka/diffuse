if (exist("i")==0 || i<0) i=0

print i

set zrange [0:0.16]
splot "fort.2222" index i ti ""
pause 0.1
i = i+1
if (i<10) reread

i = -1

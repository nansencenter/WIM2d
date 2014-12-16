
x1 = rand(150,50);
x2 = rand(150,50);
x3 = rand(150,50);
x4 = rand(150,50);
x5 = rand(150,50);
x6 = rand(150,50);
ip = [1,2];

%% call the mex
[y1,y2,y3,y4,y5]  = WIM2d_run_io_mex(x1,x2,x3,x4,x5,x6,ip);
if 0
   %%outputs=0
   norm(y1)
   norm(y2)
   norm(y3)
   norm(y4)
elseif 0
   %%outputs=inputs^2
   norm(y1-x1.^2)
   norm(y2-x2.^2)
   norm(y3-x3.^2)
   norm(y4-x4.^2-x5.^2-x6.^2)
   norm(y5-ip.^2)
elseif 1
   %%outputs=inputs^2
   norm(y1-x1.^2)
   norm(y2-x2.^2)
   norm(y3-x3.^2)
   norm(y4-x4.^2-x5.^2-x6.^2)
   norm(y5-x5.^2)
end

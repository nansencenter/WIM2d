function masks = make_masks_puv(pmask,nbdy)
%% trying to copy bigrid.f in HYCOM code, so we can use the advection code in there
%% INPUTS (eg):
%% pmask: [70x70 double] - 0 on land, 1 on water
%% nbdy  = 4 - no of ghost cells
%% OUTPUTS (eg):
%% masks = 
%%     pmask: [78x78 double] - padded with nbdy ghost cells, so size is (nx+2*nbdy)*(ny+2*nbdy)
%%                           - these are set to land (pmask=0) at the moment
%%     umask: [78x78 double]
%%     vmask: [78x78 double]
%%       isp: [78x1 double]  - no of sections in j-th row (using pmask)
%%       ifp: [78x2 double]  - size is (ny+2*nbdy)*ms, where ms is the max no of sections;
%%                             ifp(j,k) gives the first index of the k-th section
%%       ilp: [78x2 double]  - ilp(j,k) gives the last  index of the k-th section
%%       isu: [78x1 double]  - no of sections in j-th row (using umask)
%%       ifu: [78x2 double]
%%       ilu: [78x2 double]
%%       isv: [78x1 double]  - no of sections in j-th row (using vmask)
%%       ifv: [78x2 double]
%%       ilv: [78x2 double]

%% extend pmask with zeros
%% TODO: do something different with open boundaries?
[ii,jj]     = size(pmask);
masks.pmask = zeros(ii+2*nbdy,jj+2*nbdy);
masks.pmask(1+nbdy:ii+nbdy,1+nbdy:jj+nbdy)   = pmask;
masks.umask = 0*masks.pmask;
masks.vmask = 0*masks.pmask;
clear pmask;

for j=1:jj+2*nbdy
   for i=2:ii+2*nbdy
      %% leave as zero for i=1
      if (masks.pmask(i-1,j)>0&masks.pmask(i,j)>0)
         masks.umask(i,j)  = 1;
      end
   end
end
for j=2:jj+2*nbdy
   %% leave as zero for j=1
   for i=1:ii+2*nbdy
      if (masks.pmask(i,j-1)>0&masks.pmask(i,j)>0)
         masks.vmask(i,j)  = 1;
      end
   end
end

% %may not need masks.qmask
% for j=2:jj+2*nbdy
%    for i=2:ii+2*nbdy
%       %% leave as zero for i,j=1 - check this is OK?
%       %% - open boundary conditions?
%       if (min([masks.pmask(i,j),masks.pmask(i-1,j),masks.pmask(i,j-1),masks.pmask(i-1,j-1)])>0)
%          masks.qmask(i,j)  = 1;
%       elseif ((masks.pmask(i,j)>0&masks.pmask(i-1,j-1)>0) | (masks.pmask(i-1,j)>0&masks.pmask(i,j-1)>0))
%          masks.qmask(i,j)  = 1;
%       end
%    end
% end

%%get sections for each type of point
[masks.ifp,masks.ilp,masks.isp]  = indxi(masks.pmask,nbdy);
[masks.ifu,masks.ilu,masks.isu]  = indxi(masks.umask,nbdy);
[masks.ifv,masks.ilv,masks.isv]  = indxi(masks.vmask,nbdy);
%[masks.ifq,masks.ilq,masks.isq]  = indxi(masks.qmask,nbdy);

return;


function [IF,IL,IS] = indxi(ipt,nbdy)

% --- input array ipt contains 1 at grid point locations, 0 elsewhere
% --- output is arrays IF, IL, IS  where
% --- IF(j,k) gives row index of first point in column j for k-th section
% --- IL(j,k) gives row index of last point
% --- IS(j) gives number of sections in column j (maximum: ms)

ms       = 20;%max no of sections
[ii,jj]  = size(ipt);
ii       = ii-2*nbdy;
jj       = jj-2*nbdy;
IS       = zeros(jj+2*nbdy,1);
IF       = zeros(jj+2*nbdy,ms);
IL       = zeros(jj+2*nbdy,ms);

for j_=1-nbdy:jj+nbdy%loop over all j
   j     = j_+nbdy;
   k     = 1;
   last  = ipt(1,j);
   if (last == 1)
      IF(j,k) = 1;
   end
   for i_=2-nbdy:ii+nbdy
      i  = i_+nbdy;
      if (last == 1 & ipt(i,j) == 0)
         IL(j,k)  = i_-1;
         k        = k+1;
      elseif (last == 0 & ipt(i,j) == 1)
         if (k > ms)
            error('ms too small');
         end
         IF(j,k) = i_;
      end
      last = ipt(i,j);
   end
   if (last == 1)
      IL(j,k)  = ii+nbdy;
      IS(j)    = k;
   else
      IS(j) = k-1;
   end

end

ms = max(IS);
IF = IF(:,1:ms);
IL = IL(:,1:ms);

return

import sys,os
sys.path.append('../../py_funs')
from matplotlib import pyplot as plt
import numpy as np
import fns_get_data as Fdat
import fns_plot_data as Fplt

odir   = 'out'
files  = os.listdir(odir)
figdir = odir+'/figs'
if not os.path.exists(figdir):
    os.mkdir(figdir)

spacing = 20 # save a figure every 20 files
count   = spacing

if 1:# other parameters (should they be in .b files?)
    CFL = .4
    alp = 0 # atten coeff

for bfil in files:
    if ('ADVweno'==bfil[:7]) and ('.b'==bfil[-2:]):
        # check if .b file
        if count==spacing:
            count = 1 # reset count

            # open .b file to get info:
            bf  = open(odir+'/'+bfil,'r')
            lin = bf.readline() # no of  records
            lin = bf.readline() # nx
            nx  = int(lin.split()[0])
            lin = bf.readline() # ny
            ny  = int(lin.split()[0])
            lin = bf.readline() # OPT
            OPT = int(lin.split()[0])
            bf.close()

            # open .a file to get data:
            basename = bfil[:-2]
            afil     = odir+'/'+basename+'.a' 
            aid      = open(afil,'rb')
            X        = Fdat.get_array(aid,nx,ny)
            Y        = Fdat.get_array(aid,nx,ny)
            h        = Fdat.get_array(aid,nx,ny) # advected quantity
            aid.close()

            # plot and save figure
            fig = Fplt.cmap_3d(X/1.e3,Y/1.e3,h,['$x$, km','$y$, km','h'])

            #############################################################
            # plot test curves on top
            # - get some parameters from X,Y
            dx = X[1,0]-X[0,0]
            dy = Y[0,1]-Y[0,0]
            x0 = X.min()
            y0 = Y.min()
            xm = .5*(X.max()-x0)
            ym = .5*(Y.max()-y0)
            xx = X[:,0]
            yy = Y[0,:]
            n  = int(basename[-3:])

            #############################################################
            if OPT==1:
                uc = 30. # const speed m/s
                xc = x0+5*xm/3.
                # theta = 180. # deg straight across
                theta = 135. # deg
                u     = 0*X+uc*np.cos(np.pi/180.*theta)
                v     = 0*X+uc*np.sin(np.pi/180.*theta)
                #
                dt = CFL*dx/uc
                nt = 2*xm/(uc*dt)
                #
                x1 = xc+uc*np.cos(np.pi*theta/180)*n*dt
                y1 = y0+uc*np.sin(np.pi*theta/180)*n*dt
                x2 = X.max()+uc*np.cos(np.pi*theta/180)*n*dt
                #
                col = 'g'
                plt.plot(x1/1.e3+0*yy[yy>y1],yy[yy>y1]/1e3,col,linewidth=2)
                plt.plot(x2/1.e3+0*yy[yy>y1],yy[yy>y1]/1e3,col,linewidth=2)
                plt.plot(np.array([x1,x2])/1.e3,np.array([y1,y1])/1.e3,col,linewidth=2)
                plt.ylim([y0/1.e3,Y.max()/1.e3])
                plt.xlim([x0/1.e3,X.max()/1.e3])

                if 1:
                    # plot section
                    y3     = 0*xx
                    jc     = np.logical_and(xx>=x1,xx<=x2)
                    y3[jc] = ym/2.
                    fac    = ym/2/1.e3
                    plt.plot(xx/1.e3,y3/1.e3,'k')
                    fac = fac*np.exp(n*alp*dt)
                    plt.plot(xx/1.e3,h[:,-1]*fac,'--k')

            #############################################################
            # elif OPT==2:
            #     uc = 30. # const speed m/s
            #     xc = 2*xm/3.
            #     u  = -uc*X/xm
            #     v  = 0*X
            #     #
            #     dt = CFL*dx/uc
            #     nt = 2*xm/(uc*dt)

            #############################################################
            elif OPT==3:
                Rm = xm/3.
                Ym = ym/12.
                #
                angrot = (1./20.)*np.pi/180 # radian/s
                u      = 0*X
                v      = 0*X
                #
                u[R<Rm*1.45] = -Y[R<Rm*1.45]*angrot
                v[R<Rm*1.45] =  X[R<Rm*1.45]*angrot
                #
                max_speed = Rm*angrot
                dt        = CFL*dx/max_speed
                nt        = round(2*np.pi/(angrot*dt))
                dtheta    = dt*angrot
                #
                x3 = x1*np.cos(n*dtheta)-y1*np.sin(n*dtheta)
                y3 = y1*np.cos(n*dtheta)+x1*np.sin(n*dtheta)
                plt.plot(x3/1.e3,y3/1.e3,'r')
                x4 = x2*np.cos(n*dtheta)-y2*np.sin(n*dtheta)
                y4 = y2*np.cos(n*dtheta)+x2*np.sin(n*dtheta)
                plt.plot(x4/1.e3,y4/1.e3,'m')
            #############################################################

            figname  = figdir+'/'+basename+'.png'
            plt.savefig(figname)
            plt.close()
            fig.clf()
            print('saving to: '+figname)
        else:
            count = count+1

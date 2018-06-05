#Domain [-300 300]^2 km
#Origin (41.3209371228N, 289.46309961W)
#Projection Lambert
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import h5py as hp
import scipy.ndimage.filters as filter
initday = 14
inittime = 12-4+4 #Run Initialization(UTC) - 4hrs for EDT Conversion + 4hrs for integration time
stepsize = 15 #minutes
'''
amaxits = 7
rmaxits = 7
astd = 2
rstd = 2
atol =-0.002
rtol =-0.002
'''
amaxits = 9
rmaxits = 9
astd = 1
rstd = 1
atol =-0.005
rtol =-0.005
parallelres = 7
meridianres = 7
dx = 0.3
dy = 0.3

print "Loading Data"
f=hp.File('attFTLEOutput.mat','r')
attdata = f['F'][:,:,:]
f.close()
dim = attdata.shape 
f=hp.File('repFTLEOutput.mat','r')
repdata = f['F'][:,:,:]
f.close()

mydata = np.genfromtxt('Downloads\myfile.csv', delimiter=',')
print "Data is in"
# Organize velocity data
uvar = mydata[:,0]
vvar = mydata[:,1]
del mydata
u = np.empty(dim)
v = np.empty(dim)
index = 0
for t in range(dimn[0]):
    for y in range(dim[1]):
        for x in range(dim[2]):
            u[t,x,y] = uvar[index]
            v[t,x,y] = vvar[index]
            index+=1
del uvar, vvar
#Simplectic Matrix
J = np.array([[0, 1], [-1, 0]])

#Plot set up
plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
origin = [41.3209371228, -70.53690039]
print "Begin Map"
m = Basemap(width=120000,height=96000,\
    rsphere=(6378137.00,6356752.3142),\
    resolution='f',area_thresh=0.,projection='lcc',\
    lat_1=35.,lat_0=origin[0],lon_0=origin[1])
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.linspace(m.llcrnrlat,m.urcrnrlat,parallelres),labels=[True,False,False,False])
m.drawmeridians(np.linspace(m.llcrnrlon,m.urcrnrlon,meridianres),labels=[False,False,False,True])
m.drawstates()

#initialize ridge collections
print "Initializing Contours"
x = np.linspace(0, m.urcrnrx, dim[1])
y = np.linspace(0, m.urcrnry, dim[2])
xx, yy = np.meshgrid(x, y)
aridge = m.contour(xx,yy,np.transpose(attdata[0,:,:]),levels=[-1])
rridge = m.contour(xx,yy,np.transpose(repdata[0,:,:]),levels=[-1])

A = np.empty(dim)
for t in range(dim[0]):
    duy, dux = np.gradient(u[t,:,:],dx,dy,edge_order=2)
    dvy, dvx = np.gradient(v[t,:,:],dx,dy,edge_order=2)
    for m in range(dim[1]):
        for n in range(dim[2]):
            Utemp = np.array([u[t, n, m], v[t, n, m]])
            Grad = np.array([[dux[n, m], duy[n, m]], [dvx[n, m], dvy[n, m]]])
            S = 0.5*(Grad + np.transpose(Grad))
            A[t,n, m] = np.dot(Utemp, 
                                 np.dot(np.dot(np.transpose(J), 
                                               np.dot(S, J)), Utemp))
            /np.dot(Utemp, Utemp)


repulsionquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0:,:])),shading='gouraud',cmap='viridis')
print "Begin Loop"
for t in range(dim[0]):
    for c in aridge.collections:
        c.remove()
    for c in rridge.collections:
        c.remove()
    attftle = attdata[t,:,:]
    repftle = repdata[t,:,:]
    its = 0
    while its < amaxits:
        attftle = filter.gaussian_filter(attftle,sigma=astd)
        its+=1
    its = 0
    while its < rmaxits:
        repftle = filter.gaussian_filter(repftle,sigma=rstd)
        its+=1

    adx, ady = np.gradient(attftle,dx,dy,edge_order=2)
    adxdx, adydx = np.gradient(adx,dx,dy,edge_order=2)
    adxdy, adydy = np.gradient(ady,dx,dy,edge_order=2)
    
    rdx, rdy = np.gradient(repftle,dx,dy,edge_order=2)
    rdxdx, rdydx = np.gradient(rdx,dx,dy,edge_order=2)
    rdxdy, rdydy = np.gradient(rdy,dx,dy,edge_order=2)

    adirdiv = np.empty([dim[1],dim[2]])
    amineig = np.empty([dim[1],dim[2]])

    rdirdiv = np.empty([dim[1],dim[2]])
    rmineig = np.empty([dim[1],dim[2]])

    for i in range(dim[1]):
        for j in range(dim[2]):
            aeig = np.linalg.eig([[adxdx[i,j],adxdy[i,j]],[adydx[i,j],adydy[i,j]]])
            aeigmin =  np.argmin(aeig[0])
            adirdiv[i,j] = np.dot(aeig[1][:,aeigmin],[adx[i,j],ady[i,j]])
            amineig[i,j] = aeig[0][aeigmin]
            
            reig = np.linalg.eig([[rdxdx[i,j],rdxdy[i,j]],[rdydx[i,j],rdydy[i,j]]])
            reigmin =  np.argmin(reig[0])
            rdirdiv[i,j] = np.dot(reig[1][:,reigmin],[rdx[i,j],rdy[i,j]])
            rmineig[i,j] = reig[0][reigmin]
    repulsionquad.set_array(np.ravel(np.transpose(A[:,:])))
    apotridge = np.ma.masked_where(amineig>=atol,adirdiv)
    aridge = m.contour(xx, yy, np.transpose(apotridge),levels =[0],colors='blue')
    rpotridge = np.ma.masked_where(rmineig>=rtol,rdirdiv)
    rridge = m.contour(xx, yy, np.transpose(rpotridge),levels =[0],colors='red')
    minute = stepsize * t
    h, minute = divmod(minute,60)
    x, y = m(origin[1],origin[0])
    print x,y
    m.scatter(x,y,marker='*',color='g',s=20*16)
    plt.annotate('Tower',xy=(x-0.05*x,y+0.03*y),size=15)
    plt.title("FTLE, 8-{0}-2017 {1:02d}{2:02d} EDT".format(initday+(inittime+h)//24, (inittime+h)%24, minute))
    plt.savefig('MV_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')

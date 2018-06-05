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
initday = 4
inittime = 0-4+4 #Run Initialization(UTC) - 4hrs for EDT Conversion + 4hrs for integration time
stepsize = 15 #minutes
'''
amaxits = 7
rmaxits = 7
astd = 2
rstd = 2
atol =-0.002
rtol =-0.002
'''
amaxits = 3
rmaxits = 9
astd = 1
rstd = 1
atol = 0.00012
rtol =-0.005
parallelres = 7
meridianres = 7
dx = 0.3
dy = 0.3

star = [37.19838, -80.57834]

cdict = {'red':  [(0.0, 0.0000, 0.0000),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 1.0000, 1.0000)],
        'green': [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.5450, 0.5450)],
        'blue':  [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.0000, 0.0000)]}
plt.register_cmap(name='CO', data=cdict)

print "Loading Data"
f=hp.File('FTLEOutput.mat','r')
#f=hp.File('fulltimeFTLE.mat','r')
attdata = f['F'][:,:,:]
f.close()
dim = attdata.shape 
xlen = dim[1]/3
ylen = dim[2]
tlen = dim[0]
ftle = attdata[:,0::3,:]
print np.max(ftle)
eig1 = attdata[:,1::3,:]
eig2 = attdata[:,2::3,:]
del attdata, dim
#f=hp.File('repFTLEOutput.mat','r')
#repdata = f['F'][:,:,:]
#f.close()

'''
mydata = np.genfromtxt('myfile.csv', delimiter=',')
#mydata = np.genfromtxt('fulltimevel.csv', delimiter=',')
print "Data is in"
# Organize velocity data
uvar = mydata[:,0]
vvar = mydata[:,1]
del mydata
u = np.empty([tlen,xlen,ylen])
v = np.empty([tlen,xlen,ylen])
index = 0
for t in range(tlen):
    for y in range(ylen):
        for x in range(xlen):
            u[t,x,y] = uvar[index]
            v[t,x,y] = vvar[index]
            index+=1
del uvar, vvar
#Simplectic Matrix
J = np.array([[0, 1], [-1, 0]])
'''
#Plot set up
plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
#origin = [41.3209371228, -70.53690039]
origin = [37.208, -80.5803]
print "Begin Map"
m = Basemap(width=77700,height=76800,\
    rsphere=(6378137.00,6356752.3142),\
    resolution='f',area_thresh=0.,projection='lcc',\
    lat_1=35.,lat_0=origin[0],lon_0=origin[1])
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.linspace(m.llcrnrlat,m.urcrnrlat,parallelres),labels=[True,False,False,False])
m.drawmeridians(np.linspace(m.llcrnrlon,m.urcrnrlon,meridianres),labels=[False,False,False,True])
m.drawstates()
m.drawrivers()

#initialize ridge collections
print "Initializing Contours"
x = np.linspace(0, m.urcrnrx, xlen)
y = np.linspace(0, m.urcrnry, ylen)
xx, yy = np.meshgrid(x, y)

attdata = ftle
del ftle

aridge = m.contour(xx,yy,np.transpose(attdata[0,:,:]),levels=[-1])
#ridge = m.contour(xx,yy,np.transpose(repdata[0,:,:]),levels=[-1])
print "Calculating Rhodot"
'''
A = np.empty([tlen,xlen,ylen])
for t in range(tlen):
    duy, dux = np.gradient(u[t,:,:],dx,dy,edge_order=2)
    dvy, dvx = np.gradient(v[t,:,:],dx,dy,edge_order=2)
    for j in range(ylen):
        for n in range(xlen):
            Utemp = np.array([u[t, n, j], v[t, n, j]])
            Grad = np.array([[dux[n, j], duy[n, j]], [dvx[n, j], dvy[n, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            A[t,n, j] = np.dot(Utemp, np.dot(np.dot(np.transpose(J), np.dot(S, J)), Utemp))/np.dot(Utemp, Utemp)

del duy,dux,dvy,dvx,Utemp,Grad,J,S
'''
F = np.load('SWVA_eulerian.npz')
A = F['arr_0'] #Min Eigenvalue of S
del F
colormin = 2.0/3.0*A.min()
colormax = 2.0/3.0*A.max()
#colorlevel = 2/3.0*np.min(np.fabs([colormin,colormax]))
print colormin, colormax
#colorlevel = np.min(np.fabs([colormin,colormax]))
repulsionquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0,:,:])),vmin=colormin,\
        vmax=colormax,shading='gouraud',cmap='CO')
clb = plt.colorbar(repulsionquad,fraction=0.037, pad=0.02)
print "Begin Loop"
for t in range(tlen):
    for c in aridge.collections:
        c.remove()
    #for c in rridge.collections:
    #    c.remove()
    attftle = attdata[t,:,:]
    atteig1 = eig1[t,:,:]
    atteig2 = eig2[t,:,:]
    #repftle = repdata[t,:,:]
    its = 0
    while its < amaxits:
        attftle = filter.gaussian_filter(attftle,sigma=astd)
        atteig1 = filter.gaussian_filter(atteig1,sigma=astd)
        atteig2 = filter.gaussian_filter(atteig2,sigma=astd)
        its+=1
    '''
    its = 0
    while its < rmaxits:
        repftle = filter.gaussian_filter(repftle,sigma=rstd)
        its+=1
    '''
    adx, ady = np.gradient(attftle,dx,dy,edge_order=2)
    #adxdx, adydx = np.gradient(adx,dx,dy,edge_order=2)
    #adxdy, adydy = np.gradient(ady,dx,dy,edge_order=2)
    
    #rdx, rdy = np.gradient(repftle,dx,dy,edge_order=2)
    #rdxdx, rdydx = np.gradient(rdx,dx,dy,edge_order=2)
    #rdxdy, rdydy = np.gradient(rdy,dx,dy,edge_order=2)

    adirdiv = np.empty([xlen,ylen])
    #amineig = np.empty([xlen,ylen])

    #rdirdiv = np.empty([xlen,ylen])
    #rmineig = np.empty([xlen,ylen])

    for i in range(xlen):
        for j in range(ylen):
            #aeig = np.linalg.eig([[adxdx[i,j],adxdy[i,j]],[adydx[i,j],adydy[i,j]]])
            #aeigmin =  np.argmin(aeig[0])
            adirdiv[i,j] = np.dot([adx[i,j],ady[i,j]],[atteig1[i,j],atteig2[i,j]])
            #amineig[i,j] = aeig[0][aeigmin]
            
            #reig = np.linalg.eig([[rdxdx[i,j],rdxdy[i,j]],[rdydx[i,j],rdydy[i,j]]])
            #reigmin =  np.argmin(reig[0])
            #rdirdiv[i,j] = np.dot(reig[1][:,reigmin],[rdx[i,j],rdy[i,j]])
            #rmineig[i,j] = reig[0][reigmin]
    repulsionquad.set_array(np.ravel(np.transpose(A[t,:,:])))
    apotridge = np.ma.masked_where(attftle<=atol,adirdiv)
    aridge = m.contour(xx, yy, np.transpose(apotridge),levels =[0],colors='blue')
    #rpotridge = np.ma.masked_where(rmineig>=rtol,rdirdiv)
    #rridge = m.contour(xx, yy, np.transpose(rpotridge),levels =[0],colors='red')
    minute = stepsize * t
    print t
    h, minute = divmod(minute,60)
    x, y = m(star[1],star[0])
    m.scatter(x,y,marker='*',color='g',s=20*16)
    plt.annotate('Lidar',xy=(x-0.05*x,y+0.03*y),size=15)
    plt.title("FTLE, 9-{0}-2017 {1:02d}{2:02d} EDT".format(initday+(inittime+h)//24, (inittime+h)%24, minute))
    plt.savefig('MV_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
import matplotlib
matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import Tkinter as tk
import h5py as hp
from mpl_toolkits.basemap import Basemap
import scipy.ndimage.filters as filter
timestep = 0
domain = 'MV'
if domain == 'MV':
    originlat = 41.3209371228
    originlon =-70.53690039
    xres = 112
    yres = 96
    xstep = 0.3
    ystep = 0.3
    mapres ='f'
    areathresh = 0
    meridianres = 6
    parallelres = 6
    defaultstepsize = 15
elif domain =='CA':
    originlat = 36.7964
    originlon =-120.822
    xres = 2000
    yres = 2000
    xstep = 1.2
    ystep = 1.2
    mapres ='i'
    areathresh = 100
    meridianres = 6
    parallelres = 6
    defaultstepsize = 15
elif domain == 'KF':    
    originlat = 37.208
    originlon =-80.5803
    xres = 77.700
    yres = 76.800
    xstep = 0.3
    ystep = 0.3
    mapres ='i'
    areathresh = 100
    meridianres = 6
    parallelres = 6
    defaultstepsize = 15

usetime = 0

#f=hp.File('attFTLEOutput.mat','r')
f=hp.File('FTLEOutput.mat','r')
#f=hp.File('repFTLEOutput.mat','r')
ftledata = f['F'][:,::3,:]
eig1data = f['F'][:,1::3,:]
eig2data = f['F'][:,2::3,:]
f.close()

print np.max(ftledata)

#Button Functions-----------------------------------------------------------
def _nextTime(timeincrement):
    global ridges
    for c in ridges.collections:
        c.remove()
    global timestep
    timestep = timestep + timeincrement
    ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])
    quad.set_array(np.ravel(np.transpose(np.squeeze(ftledata[timestep,:,:]))))
    timeentry.delete(0,tk.END)
    timeentry.insert(tk.END,timestep)
    if usetime == 1:
        month, day, year, hour, minute = _calculateTime()
        ax.set_title("{0:02d}-{1:02d}-{2:04d},  {3:02d}:{4:02d}".format(month, day, year, hour, minute))    
    else:
        ax.set_title("Time Step {0}".format(timestep))
    canvas.draw()

def _jumpTime():
    global ridges
    for c in ridges.collections:
        c.remove()
    newtimestep = int(timeentry.get())
    global timestep
    timestep = newtimestep
    ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])
    quad.set_array(np.ravel(np.transpose(np.squeeze(ftledata[timestep,:,:]))))
    if usetime == 1:
        month, day, year, hour, minute = _calculateTime()
        ax.set_title("{0:02d}-{1:02d}-{2:04d},  {3:02d}:{4:02d}".format(month, day, year, hour, minute))    
    else:
        ax.set_title("Time Step {0}".format(timestep))
    timeentry.delete(0,tk.END)
    timeentry.insert(tk.END,timestep)
    canvas.draw()

def _runFilter():
    global ridges
    for c in ridges.collections:
        c.remove()
    ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])
    filterdata = np.squeeze(ftledata[timestep,:,:])
    maxits = int(iterationentry.get())
    print maxits
    iterations = 0
    stddev = float(deviationentry.get())
    print stddev
    while iterations < maxits:
        print iterations
        filterdata = filter.gaussian_filter(filterdata,sigma=stddev)
        iterations += 1
    quad.set_array(np.ravel(np.transpose(filterdata[:,:])))
    canvas.draw()
    
def _extractLCS():
    global ridges
    #ind = 0
    for c in ridges.collections:
        #if ind == 0
         #   print c.get_segments()
         #   ind += 1
        c.remove()
    filterftledata = np.squeeze(ftledata[timestep,:,:])
    filtereig1data = np.squeeze(eig1data[timestep,:,:])
    filtereig2data = np.squeeze(eig2data[timestep,:,:])
    maxits = int(iterationentry.get())
    print maxits
    iterations = 0
    stddev = float(deviationentry.get())
    print stddev
    while iterations < maxits:
        print iterations
        filterftledata = filter.gaussian_filter(filterftledata,sigma=stddev)
        filtereig1data = filter.gaussian_filter(filtereig1data,sigma=stddev)
        filtereig2data = filter.gaussian_filter(filtereig2data,sigma=stddev)
        iterations += 1
    dx, dy = np.gradient(filterftledata,xstep,ystep)
    #dxdx, dydx = np.gradient(dx,xstep,ystep)
    #dxdy, dydy = np.gradient(dy,xstep,ystep) 
    dirdiv = np.empty([dim[1],dim[2]])
    #minimumeig = np.empty([dim[1],dim[2]])
    for i in range(dim[1]):
        for j in range(dim[2]):
            #eig = np.linalg.eig([[dxdx[i,j],dxdy[i,j]],[dydx[i,j],dydy[i,j]]])
            #eigmin =  np.argmin(eig[0])
            dirdiv[i,j] = np.dot([filtereig1data[i,j],filtereig2data[i,j]],[dx[i,j],dy[i,j]])
            #minimumeig[i,j] = eig[0][eigmin]    
    #tol = float(eigthresholdentry.get())
    fthresh = float(ftlethresholdentry.get())
    #potridge = np.ma.masked_where(minimumeig>=tol,dirdiv)
    potridge = np.ma.masked_where(ftledata[timestep,:,:]<fthresh,dirdiv)
    ridges = m.contour(xx, yy, np.transpose(potridge),levels =[0],colors='orange')
    canvas.draw()
    
def _calculateTime():
    timezones={'UTC':0,'EST':-5,'EDT':-4,'CST':-6,'CDT':-5,'MST':-7,'MDT':-6,'PST':-8,'PDT':-7}
    hr_adj = timezones[datatime.get()]-timezones[desiredtime.get()]
    initialtime = initialtimeentry.get()
    initmin = int(initialtime[2:])
    inithr = int(initialtime[:2])
    minute = int(timeentry.get())*float(stepsizeentry.get())+initmin
    hour, minute = divmod(minute,60)
    hour = hour + inithr - hr_adj
    day, hour = divmod(hour,24)
    initialdate = initialdateentry.get()
    initmon = int(initialdate[:2])
    initday = int(initialdate[2:4])
    inityear = int(initialdate[4:])
    day = day + initday
    month = initmon
    year = inityear
    return int(month), int(day), int(year), int(hour), int(minute)

def setTitle():
    if usetime == 1:
        month, day, year, hour, minute = _calculateTime()
        ax.set_title("{0:02d}-{1:02d}-{2:04d},  {3:02d}:{4:02d}".format(month, day, year, hour, minute))
    else:
        ax.set_title("Time Step {0}".format(timestep))
        
def _timeFunc():
    global usetime
    if usetime == 0:
        usetime = 1
    else:
        usetime = 0
        
    setTitle()
    canvas.draw()
    
#Main -----------------------------------------------------------------------
root = tk.Tk()
root.wm_title("Embedding in TK")
#root.filename = tkFileDialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
#print (root.filename)

#Interface Objects-----------------------------------------------------------
initialdatelabel = tk.Label(root, text="Initial Date (MMDDYYYY)")
initialdateentry = tk.Entry(root)
initialdateentry.insert(tk.END,'01011970')

initialtimelabel = tk.Label(root, text="Initial Time (HHMM)")
initialtimeentry = tk.Entry(root)
initialtimeentry.insert(tk.END,'0000')

datatime = tk.StringVar(root)
datatime.set('UTC') # default value
datatimezonelabel = tk.Label(root, text="Data Time Zone")
datatimezone = tk.OptionMenu(root, datatime, 'UTC','EST','EDT','CST','CDT','MST','MDT','PST','PDT')

desiredtime = tk.StringVar(root)
desiredtime.set('UTC') # default value
desiredtimezonelabel = tk.Label(root, text="Desired Time Zone")
desiredtimezone = tk.OptionMenu(root, desiredtime, 'UTC','EST','EDT','CST','CDT','MST','MDT','PST','PDT')

stepsizelabel = tk.Label(root, text="Step Size (minutes)")
stepsizeentry = tk.Entry(root)
stepsizeentry.insert(tk.END,defaultstepsize)

timelabel = tk.Label(root, text="Time Step")
timeentry = tk.Entry(root)
timeentry.insert(tk.END,timestep)

deviationlabel = tk.Label(root, text="Filter Standard Deviation")
deviationentry = tk.Entry(root)
deviationentry.insert(tk.END,0)

iterationlabel = tk.Label(root, text="Filter Iteration")
iterationentry = tk.Entry(root)
iterationentry.insert(tk.END,0)

eigthresholdlabel = tk.Label(root, text="Eigenvalue Threshold")
eigthresholdentry = tk.Entry(root)
eigthresholdentry.insert(tk.END,'-0.0')

ftlethresholdlabel = tk.Label(root, text="FTLE Threshold, max = {0:0.03f}".format(np.max(ftledata)))
ftlethresholdentry = tk.Entry(root)
ftlethresholdentry.insert(tk.END,'0.0')

nextbutton = tk.Button(root, text='Next Step', command= lambda: _nextTime(1))
previousbutton = tk.Button(root, text='Previous Step', command= lambda: _nextTime(-1))
jumpbutton = tk.Button(root, text='Jump to Time', command= lambda: _jumpTime())
filterbutton = tk.Button(root, text='Run Filter', command= lambda: _runFilter())
extractbutton = tk.Button(root, text='Extract LCS', command= lambda: _extractLCS())#attdata,aridges))
timebutton = tk.Button(root, text='Toggle Time', command= lambda: _timeFunc())

#Object Placement-----------------------------------------------------------
initialdatelabel.grid(row=0,column=0,sticky="e")
initialdateentry.grid(row=0,column=1,padx=10)

initialtimelabel.grid(row=1,column=0,sticky="e")
initialtimeentry.grid(row=1,column=1,padx=10)

datatimezonelabel.grid(row=2,column=0,sticky='e')
datatimezone.grid(row=2,column=1,padx=10,sticky='w')

desiredtimezonelabel.grid(row=3,column=0,sticky='e')
desiredtimezone.grid(row=3,column=1,padx=10,sticky='w')

timelabel.grid(row=4,column=0,sticky="e")
timeentry.grid(row=4,column=1,padx=10)

stepsizelabel.grid(row=5,column=0,sticky="e")
stepsizeentry.grid(row=5,column=1,padx=10)

deviationlabel.grid(row=6,column=0,sticky="e")
deviationentry.grid(row=6,column=1,padx=10)

iterationlabel.grid(row=7,column=0,sticky="e")
iterationentry.grid(row=7,column=1,padx=10)

eigthresholdlabel.grid(row=8,column=0,sticky="e")
eigthresholdentry.grid(row=8,column=1,padx=10)

previousbutton.grid(row=9,column=0,sticky="e")
nextbutton.grid(row=9,column=1,sticky="w",padx=10)

jumpbutton.grid(row=10,column=0,sticky="e")
filterbutton.grid(row=10,column=1,sticky="w",padx=10)

extractbutton.grid(row=11,column=1,sticky="w",padx=10)
timebutton.grid(row=12,column=1,sticky="w",padx=10)

ftlethresholdlabel.grid(row=13,column=0,sticky="e")
ftlethresholdentry.grid(row=13,column=1,padx=10)
# Adding graph to Canvas--------------------------------------------


dim = ftledata.shape
fig = Figure(figsize=(10, 7), dpi=100)
ax = fig.add_subplot(111)

m = Basemap(width=xres*1000,height=yres*1000,\
    rsphere=(6378137.00,6356752.3142),\
    resolution=mapres,area_thresh=areathresh,projection='lcc',\
    lat_1=35.,lat_0=originlat,lon_0=originlon,ax=ax)

def format_coord(i, j):
    #return m(i,j,inverse=True)
    return i/1000-1000,j/1000-1000

ax.format_coord = format_coord

m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.linspace(m.llcrnrlat,m.urcrnrlat,parallelres),labels=[True,False,False,False])
m.drawmeridians(np.linspace(m.llcrnrlon,m.urcrnrlon,meridianres),labels=[False,False,False,True])
m.drawstates()
m.drawrivers()
x = np.linspace(0, m.urcrnrx, dim[1])
y = np.linspace(0, m.urcrnry, dim[2])
xx, yy = np.meshgrid(x, y)

quad = m.pcolormesh(xx,yy, np.transpose(np.squeeze(ftledata[timestep,:,:])),shading='gouraud',cmap='viridis')
#quad = m.contourf(xx,yy, np.transpose(np.squeeze(ftledata[timestep,:,:])),shading='gouraud',cmap='viridis')
cbar = fig.colorbar(quad,ax=ax,pad=0.01)
ax.set_title('Time Step {0}'.format(timestep))
ridges = m.contour(xx,yy,np.transpose(np.squeeze(ftledata[timestep,:,:])),levels =[-1])

canvas = FigureCanvasTkAgg(fig,root)
canvas.show() 
canvas.get_tk_widget().grid(row=1,column=3,rowspan=40)


toolbar_frame = tk.Frame(root)
toolbar_frame.grid(row=0,column=3,sticky="w")
toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame) 
toolbar.update() 


#root.mainloop()
tk.mainloop()
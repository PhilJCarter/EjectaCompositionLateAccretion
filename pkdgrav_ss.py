### pkdgrav_ss.py

"""Classes and functions for accessing PKDGRAV data."""

import struct
import numpy as npy
import os

class ssHeader:
    """
    PKDGRAV ss file header class
    
    t - time (year/(2*pi))
    N - number of particles
    MagicNumber - ss magic number, indicates normal or reduced format [-1]
    """
    def __init__(self, t=0, N=0, magic=0):
        self.t = t
        self.N = N
        self.MagicNumber = magic
        
class ss:
    """
    PKDGRVAV ss file class
    
    header - ss file header
    m - particle mass (solar masses)
    r - particle radius (AU)
    x,y,z - particle coordinates (AU)
    vx,vy,vz - particle velocity componenets (AU/year/(2*pi))
    sx,sy,sz - particle spin components [spin is not used]
    org_idx - original particle index
    color - particle color, indicates type
                [3  = planetesimal
                 2  = Jupiter
                 11 = Saturn]    
    (id - particle ID)
    (cf - core fraction)
    (origin - mass origin hsitogram)
    (corigin - core origin histogram)
        
    
    Methods:
    
    read(<file>,extras=False)
        - read ss file <file>
        - if extras is True, will try to read origin histograms
    
    calcOE(mu=1.)
        - calculates orbital elements from cartesian coords
          and velocities
        - assumes a central object of 1 solar mass by default (mu)
        """
    
    def __init__(self):
        self.header = ssHeader()
        self.m = npy.zeros(self.header.N)
        self.r = npy.zeros(self.header.N)
        self.x = npy.zeros(self.header.N)
        self.y = npy.zeros(self.header.N)
        self.z = npy.zeros(self.header.N)
        self.vx = npy.zeros(self.header.N)
        self.vy = npy.zeros(self.header.N)
        self.vz = npy.zeros(self.header.N)
        self.sx = npy.zeros(self.header.N)
        self.sy = npy.zeros(self.header.N)
        self.sz = npy.zeros(self.header.N)
        self.color = npy.zeros(self.header.N).astype(int)
        self.org_idx = npy.zeros(self.header.N).astype(int)

    def read(self, fname, extras=False):
        dbytes = 8
        ibytes = 4
        with open(fname, 'rb') as fo:
            self.header.t = struct.unpack('>d', fo.read(dbytes))[0]
            self.header.N = struct.unpack('>i', fo.read(ibytes))[0]
            self.header.MagicNumber = struct.unpack('>i',fo.read(ibytes))[0]

            # initialise arrays
            self.m = npy.zeros(self.header.N)
            self.r = npy.zeros(self.header.N)
            self.x = npy.zeros(self.header.N)
            self.y = npy.zeros(self.header.N)
            self.z = npy.zeros(self.header.N)
            self.vx = npy.zeros(self.header.N)
            self.vy = npy.zeros(self.header.N)
            self.vz = npy.zeros(self.header.N)
            self.sx = npy.zeros(self.header.N)
            self.sy = npy.zeros(self.header.N)
            self.sz = npy.zeros(self.header.N)
            self.color = npy.zeros(self.header.N).astype(int)
            self.org_idx = npy.zeros(self.header.N).astype(int)

            # read data
            for j in range(0, self.header.N):
                self.m[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.r[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.x[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.y[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.z[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.vx[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.vy[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.vz[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.sx[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.sy[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.sz[j] = struct.unpack('>d', fo.read(dbytes))[0]
                self.color[j] = struct.unpack('>i', fo.read(ibytes))[0]
                self.org_idx[j] = struct.unpack('>i', fo.read(ibytes))[0]
        
        
        if extras:  #try to load origin histograms, id and core fraction
            if os.path.isfile(fname+'.iord'):
                self.id = npy.loadtxt(fname+'.iord',skiprows=1,unpack=True).astype(int)
            if os.path.isfile(fname+'.cfrac'):
                self.cf = npy.loadtxt(fname+'.cfrac',skiprows=1,unpack=True)
            if os.path.isfile(fname+'.origin_bins'):
                self.origin = npy.transpose(npy.loadtxt(fname+'.origin_bins',skiprows=1,unpack=True))
            if os.path.isfile(fname+'.core_origin'):
                self.corigin = npy.transpose(npy.loadtxt(fname+'.core_origin',skiprows=1,unpack=True))


 
    def calcOE(self,mu=1.):
        d=npy.sqrt(self.x**2+self.y**2+self.z**2)
        v2=self.vx**2+self.vy**2+self.vz**2
        rv=self.x*self.vx+self.y*self.vy+self.z*self.vz
        hx=self.y*self.vz - self.z*self.vy
        hy=self.z*self.vx - self.x*self.vz
        hz=self.x*self.vy - self.y*self.vx
        hxy=hx**2+hy**2
        h2=hxy+hz**2

        self.i=npy.arctan2(npy.sqrt(hxy),hz)

        E=0.5*v2 - mu/d
        self.a=npy.where(E<0,-0.5*mu/E,npy.where(E>0,0.5*mu/E,2*mu/v2))
        self.e=npy.where(E<0,npy.sqrt(1.-h2/(mu*self.a),where=E<0),npy.where(E>0,npy.sqrt(1.+h2/(mu*self.a)),1.))



class collisions:
    """
    pkdgrav collision data class
    
    t - time (yr)
    idt - target id
    idp - projectile id
    mt - target mass (solar masses)
    mp - projectile mass (solar masses)
    coret - target core mass (solar masses)
    corep - projectile core mass (solar masses)
    v - impact velocity (cm/s)
    vesc - escape velocity (cm/s)
    b - impact parameter
    type - numerical collision outcome classification
        [1.0 = perfect merge
         2.0 = partial accretion
         3.2 = Hit-and-run - projectile intact
         3.6 = Hit-and-run - projectile disrupted
         3.8 = Hit-and-run - projectile supercat disrupted
         4.0 - Erosion
         5.0 - Supercatastrophic disruption]
    nfrag - number of post-collision remnants
    mlr - largest remnant mass (solar masses)
    mslr - second largest remnant mass (solar masses)
    corelr - largest remnant core mass (solar masses)
    coreslr - second largest remnant core mass (solar masses)
    idlr - largest remnant id
    idslr - second largest remnant id
    m2dust - mass placed in unresolved debris
    alr - semi-major axis largest remnant orbit (au)
    elr - eccentricity largest remnant orbit
    ilr - inclination largest remnant orbit
    aslr - semi-major axis second largest remnant orbit (au)
    eslr - eccentricity second largest remnant orbit
    islr - inclination second largest remnant orbit
    IDrem1 - id of 1st remnant
    IDrem2 - id of 2nd remnant
    iprevt - index of previous collision of target
    iprevp - index of previous collision of projectile
    
    
    Methods:
    
    load(path='./',parfile='ss.par')
        - load collision data
    get_history(<pID>)
        - extract collisions experienced by particle <pID>
        - returns array of collision indices in history of <pID>
    """
    def __init__(self,path='./'):
        self.t = None
        self.idt = None
        self.idp = None
        self.mt = None
        self.mp = None
        self.coret = None
        self.corep = None
        self.v = None
        self.vesc = None
        self.b = None
        self.type = None
        self.nfrag = None
        self.mlr = None
        self.mslr = None
        self.corelr = None
        self.coreslr = None
        self.idlr = None
        self.idslr = None
        self.m2dust = None
        self.alr = None
        self.elr = None
        self.ilr = None
        self.aslr = None
        self.eslr = None
        self.islr = None
        self.IDrem1 = None
        self.IDrem2 = None
        self.iprevt = None
        self.iprevp = None
        
        self.load(path=path)

    def load(self,path='./',parfile='ss.par'):
        fname1 = 'collisions_makeearth.txt'
        fname2 = 'collisions_vQ.txt'
        (t,id1,id2,m1,m2,self.type,nfrag,self.mlr,self.mslr,mc1,mc2,
            self.corelr,self.coreslr,rt_exp,rp_exp,
                self.b) = npy.loadtxt(path+fname1, unpack=True)
        (self.v,self.vesc,self.alr,self.elr,self.ilr,self.aslr,
            self.eslr,self.islr,self.m2dust) = npy.loadtxt(path+fname2, usecols=(0,1,9,10,11,12,13,14,15),unpack=True)

        self.idt=npy.where(m1>=m2,id1,id2).astype(int)
        self.idp=npy.where(m1>=m2,id2,id1).astype(int)
        self.mt=npy.where(m1>=m2,m1,m2)
        self.mp=npy.where(m1>=m2,m2,m1)
        self.coret=npy.where(m1>=m2,mc1,mc2)
        self.corep=npy.where(m1>=m2,mc2,mc1)
        
        self.t = t/(2*npy.pi)
        self.nfrag = nfrag.astype(int)
        
        (iprev1,iprev2,targID,projID,self.IDrem1,
            self.IDrem2) = npy.loadtxt(path+'collarray.txt',unpack=True,dtype='int')
        self.iprevt = npy.where(m1>=m2,iprev1,iprev2)
        self.iprevp = npy.where(m1>=m2,iprev2,iprev1)

        
        # if H&R - proj intact, or only 1 fragment, lr has target id, otherwise has ID of collider 1. 
        self.idlr = npy.where( (npy.abs(self.type-3.2)<0.00001) + (self.nfrag < 2), self.idt, id1 )
        self.idslr = npy.where( self.idlr == self.idt, self.idp, self.idt )
        
        self.idlr = npy.where( self.mlr > 0, self.idlr, -1000 )
        self.idslr = npy.where( self.mslr > 0, self.idslr, -1000 )

      
    
        

    def get_history(self,pID):    
        objcol = npy.array([]).astype(int)
        
        j = npy.where((self.IDrem1==pID)+(self.IDrem2==pID))
        j = j[0]
        if len(j)>0:
            objcol = npy.append(objcol, j[-1])
    
        oldl = 0
        newl = 1
        n = 0
        prevn = 0
        while newl>oldl:
            oldl = len(objcol)
            for l in range (prevn,len(objcol)):
                ic = self.iprevt[objcol[l]]
                if ic >= 0 and ic not in objcol:
                    objcol = npy.append(objcol,ic)
                    n += 1
                ic = self.iprevp[objcol[l]]
                if ic >= 0 and ic not in objcol:
                    objcol = npy.append(objcol,ic)
                    n += 1
            prevn = oldl
            newl = len(objcol)
    
        print('Found', newl, 'collisions for object:', pID)
        return objcol





def readparam(fname):
    """
    readparam(<file>)
        - read pkdgrav parameter file
        - returns a dictionary of parameter names and values
    """
    data = {}
    with open(fname) as fo:
        for line in fo.readlines():
            if line[0] == '#':
                continue
            else:
                key,val = (line.replace(' ','').replace('\t','').split('#')[0]).split('=')
                try:
                    val = float(val)
                except:
                    val = str(val)
                data[key] = val
                #print key,val
    return data



if __name__ == "__main__":
    import sys
    import numpy as npy
    import matplotlib.pyplot as plt

    file = 'ssic.ss'
    if len(sys.argv) > 1:
        file = sys.argv[1]

    pkdss = ss()
    pkdss.read(file,extras=True)

    print(pkdss.header.N)   # print number particles
    print(pkdss.m[pkdss.color==2])   # print mass of Jupiter
    

    # calculate orbital elements and plot a vs e
    d=npy.sqrt(pkdss.x**2+pkdss.y**2+pkdss.z**2)
    v2=pkdss.vx**2+pkdss.vy**2+pkdss.vz**2
    rv=pkdss.x*pkdss.vx+pkdss.y*pkdss.vy+pkdss.z*pkdss.vz
    hx=pkdss.y*pkdss.vz - pkdss.z*pkdss.vy
    hy=pkdss.z*pkdss.vx - pkdss.x*pkdss.vz
    hz=pkdss.x*pkdss.vy - pkdss.y*pkdss.vx
    hxy=hx**2+hy**2
    h2=hxy+hz**2

    i=npy.arctan2(npy.sqrt(hxy),hz)

    mu = 1.
    E=0.5*v2 - mu/d
    a=npy.where(E<0,-0.5*mu/E,npy.where(E>0,0.5*mu/E,2*mu/v2))
    e=npy.where(E<0,npy.sqrt(1-h2/(mu*a)),npy.where(E>0,npy.sqrt(1+h2/(mu*a)),1.))

    
    plt.scatter(a,e)
    plt.axis([1.,15.,0.,1.])
    plt.xlabel('a (au)')
    plt.ylabel('e')
    plt.show()
           

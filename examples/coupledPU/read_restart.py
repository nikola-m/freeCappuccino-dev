#!/usr/bin/env python3
"""Read freeCappuccino 'restart' (gfortran unformatted sequential). Usage: read_restart.py restart nx ny
Prints Ghia Re=100 centreline comparison (uniform square cavity, cell order i fastest)."""
import sys, struct, numpy as np
def records(fn):
    b=open(fn,'rb').read(); o=0; out=[]
    while o<len(b):
        n=struct.unpack('<i',b[o:o+4])[0]; out.append(b[o+4:o+4+n]); o+=8+n
    return out
def load(fn):
    r=records(fn); k=0
    it=struct.unpack('<id',r[0][:12]) if len(r[0])>=12 else None
    # records: (itime,time) [gradpcmf] flmass u v w p ...
    d=[np.frombuffer(x,dtype='<f8') for x in r[1:]]
    return r,d
if __name__=='__main__':
    fn,nx,ny=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
    r,d=load(fn)
    # find u: first record after flmass whose length > ncells
    nc=nx*ny
    big=[x for x in d if len(x)>=nc]
    u,v,w,p=big[1],big[2],big[3],big[4]
    U=u[:nc].reshape(ny,nx); V=v[:nc].reshape(ny,nx)
    xc=(np.arange(nx)+.5)/nx; yc=(np.arange(ny)+.5)/ny
    gy=np.array([0,.0547,.0625,.0703,.1016,.1719,.2813,.4531,.5,.6172,.7344,.8516,.9531,.9609,.9688,.9766,1])
    gu=np.array([0,-.03717,-.04192,-.04775,-.06434,-.10150,-.15662,-.21090,-.20581,-.13641,.00332,.23151,.68717,.73722,.78871,.84123,1])
    um=np.interp(gy,yc,U[:,nx//2-1]*.5+U[:,nx//2]*.5)
    print('max |u_cl - Ghia| =',np.abs(um[1:-1]-gu[1:-1]).max())

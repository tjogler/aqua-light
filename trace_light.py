#!/usr/bin/env python

import os,sys,glob
import numpy as np
from math import degrees
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pylab as P

class led(object):
    '''
    class that makes leds and reads light distribution from a cvs file
    File mus have 7 header lines and then tabulated X,Y
    In case anyone prefers other formats please add corresponding function
    '''
    def __init__(self,filename,wavelength,power):
        self.wavelength=wavelength
        self.fileName=filename
        self.power=power
        self.read_data()
        
    def read_data(self):
        x=[]
        y=[]
        try:
            f=open(self.fileName,'r+')
        except:
            print 'ERROR file %s coud not be opened'%self.fileName
            exit()
        a=f.readlines()
        a=a[6:]
        for line in a:
            val=line.split(',')
            x.append(val[0])
            y.append(val[1])
        self.angle=np.array(x,dtype=float)
        hi=1./100*np.array(y,dtype=float)
        integral=0.5*np.trapz(hi,x=self.angle)
        self.intensity=self.power/integral*hi
        
class lamp(object):
    '''
    class that creates a lamp object, this class also calculates the light profile at any given level
    The lamp consists of led objects and the coordinates of each led must be passed as x,y lists or arrays
    '''
    def __init__(self,x=[],y=[],z=0,led=[],angle=[],bind=[]):
        self.x=x
        self.y=y
        self.z=z
        self.leds=led
        self.angle=angle #should be in rad
        self.area=-1
        self.area_z=0
        self.bind=bind
        self.intensity=-1
        self.overIntens=-1
        self.phisteps=1.
        self.alphasteps=1.
        self.nWater=1.33 #no wavelength dependence yet but that's anyway minor effect
        self.nGlass=1.52 #no wavelength dependence yet but that's anyway minor effect
        self.calc_refraction_angle(1.0,self.nWater)
        self._init_rays()
        
        #print self.angle*180./np.pi

    def _init_rays(self):
        c=0
        self.rays=[] #contains all the direction vectors for each light source position
        ray_list=[]
        for x in self.x: #loop over all lights
            y=self.y[c]
            alpha=self.angle[c]
            #print alpha
            phisteps=self.phisteps
            nRaysPerPhi=int(np.ceil(alpha/(self.alphasteps*np.pi/180.)))
            nRays=int(np.ceil(alpha/(self.alphasteps*np.pi/180.))*phisteps)
            alphastep=alpha/nRaysPerPhi #angle stepsize in rad
            print 'Number of rays per lightsource: %s'%nRays
            print '%s rays per phi and %s phi angles'%(alpha/alphastep,nRaysPerPhi)
            print self.leds[c].angle
            for a in np.linspace(0,alpha,nRaysPerPhi): #loop over all angles in single light cone
                totRef=self.calc_total_reflection(a,self.nWater,self.nGlass)
                angle=a*180/np.pi
                intens=np.interp(angle,self.leds[c].angle,self.leds[c].intensity)/nRaysPerPhi
                for p in np.linspace(0.,2.*np.pi-(2.*np.pi/phisteps),phisteps-1): # calc for each y position in light cone
                    #print 'alpha: %s   phi: %s'%(a,p)
                    #if np.sqrt(np.tan(a+alphastep)**2+yy**2) < self.radius[c]:
                    r=np.tan(a)*1
                    ray=[np.array([r*np.cos(p),r*np.sin(p),-1]),np.array([x,y,self.z]),intens,a,totRef]#[0] direction, [1] start point
                    ray_list.append(ray)
                    #if a==0:
                     #   break #a==0 means you have r==0 and hence all rays are the same, so to save computing time only one is calculated
            c+=1
        self.rays=ray_list
        
        

    def _set_phisteps(self,phistep):
        self.phisteps=phistep
            
    def _set_alphasteps(self,alphastep):
        self.alphasteps=alphastep        

    def calc_intersection(self,plane,ray,limits):
        inter=False
        for p in plane:
            k=np.dot(p[0],(ray[1]-p[1]))
            #print 'p,k: ',p,k
            try:
                t=-k/np.dot(p[0],ray[0])
                #print 't: ',t
            except:
                t=np.nan
                S=np.nan
            if (t!=np.nan) and (t>0): #all rays must propergate downwards (t>0)
                #print 't: %s'%t
                S=ray[0]*t+ray[1]
                #print 'S:',S
                #print 'S %s'%S
                #print (S<=limits).all(), (S>=-limits).all()
                if (S<=limits).all() and (S>=-limits).all() :
                    #print type(S)
                    #print 'S: %s  Limit:%s'%(S,limits)
                    inter=True
                    break
        #print S,p,limits,(S<=limits).all() and (S>=-limits).all()
        return inter,S,p

    def calc_new_ray(self,iPoint,ray,plane):
        #print ray[3]*180./np.pi
        if ray[3]*180./np.pi > 35.15: #no reflection but just refraction and ray leaves VOI
            newray=-999 # -999 flaqed for deletion 
            return newray
        new=plane[0].copy()   #determines which component of the vector has to be reflected
        #print plane
        np.place(new,new==0,1)
        #print new
        newdirect=ray[0]*new # relfects the required vector component without changing the others
        newray=[newdirect,iPoint,ray[2],ray[3],ray[4]] # the intensity is not changed when reflected that should be modified in the future 
        return newray

    def get_intensity_level(self,planes,level,limits):
        planes[4][1]=np.array([1,1,level])
        limits[2]=abs(level) # need to change the z limit to the current level z so that the intersect function only breaks if ray hits the bottom plane within x_limit,y_limit
        counter =0
        temp_ray_list=[]
        intensX=[]
        intensY=[]
        intensP=[]
        bookkeeping=[]
        
        print 'Get_intensity_level: %s Rays at level'%len(self.rays)
        ntrc=0
        for r in self.rays: # check for each ray if it intersects with any plane

            res,S,p=self.calc_intersection(planes,r,limits)
            
            if res and S[2]!=level and r[4]:
                newRay=self.calc_new_ray(S,r,p)
                if newRay!=-999:
                    #print counter
                    #print  self.rays[counter]
                    #print  newRay
                    #print '\n'
                    self.rays[counter]=newRay
                    temp_ray_list.append(newRay)
                    bookkeeping.append(counter)
                else:
                    del self.rays[counter]
            if res and S[2]!=level and not r[4]:
                ntrc+=1
                del self.rays[counter]
            if res and S[2]==level:
                #print counter
                intensX.append(S[0])
                intensY.append(S[1])
                intensP.append(r[2])
            if not res:
                del self.rays[counter]
            counter+=1
                
        while temp_ray_list!=[]:
            
            #print len(temp_ray_list)
            #print temp_ray_list
            c=0
            book=[]
            for t in temp_ray_list:
                #print t
                res,S,p=self.calc_intersection(planes,t,limits)
                if res and S[2]!=level:
                    newRay=self.calc_new_ray(S,t,p)
                    if newRay!=-999:
                        self.rays[bookkeeping[c]]=newRay
                        temp_ray_list[c]=newRay
                        
                    else:
                        del self.rays[bookkeeping[c]]
                        
                if res and S[2]==level:
                    intensX.append(S[0])
                    intensY.append(S[1])
                    intensP.append(r[2])
                    del temp_ray_list[c]
                    #print temp_ray_list
                if not res:
                    del temp_ray_list[c]
                    del self.rays[bookkeeping[c]]
                c+=1
                #print counter

        print '%s Rays left volume because of no total reflection'%ntrc
        #intensDist=np.histogram2d(intensX,intensY,weights=intensP,bins=100)
        return np.array(intensX),np.array(intensY),np.array(intensP)
                    
                
    def calc_refraction_angle(self,n1,n2):
        self.angle=np.arcsin(n1*np.sin(self.angle)/n2)

    def calc_total_reflection(self,alpha,n1,n2):
        beta=90-(np.arcsin(n1/n2*np.sin(alpha)))*180/np.pi
        thetaC=np.arcsin(n1/n2)*180/np.pi
        if beta>thetaC:
            return True
        else:
            return False
        
        
    def get_ind_area(self,z):
        self.area_z=z
        self.area=np.pi*(z*np.tan(self.angle))**2.
        
    def get_radius(self,z):
         self.radius=abs(z*np.tan(self.angle))
        
    def get_intens(self,z):
        self.get_ind_area(z)
        self.radius=abs(z*np.tan(self.angle))
        self.intensity=self.power/self.area
        
    def get_overlay_intens(self,z,levelx,levely):
        self.overIntens=np.meshgrid(levelx,levely)
        self.intenity=self.get_intens(z)
            
        x2,y2 = np.meshgrid(levelx,levely)
        self.overIntens=np.zeros_like(x2)
        counter=0
        self.intensX=levelx
        self.intensY=levely
        for x in self.x:
            mask=(np.sqrt((levelx-x)**2.+(levely-self.y[counter])**2.)<self.radius[counter])
            self.overIntens+=self.intensity[counter]*(np.sqrt((x2-x)**2.+(y2-self.y[counter])**2.)<self.radius[counter])
            counter+=1
            
    def set_ray_resolution(self,nalpha,nphi):
        '''
        Sets the number of steps in alpha and in phi
        Number of rays should not be too large because of computation time, usually 2e5 should suffice, just use a bit larger bins when plotting the light distribution
        '''
        self._set_alphasteps(nalpha)
        self._set_phisteps(nphi)
        self._init_rays()
        
        
class tank(object):

    '''
    Class that defines the aquarium tanks dimensions
    '''

    
    def __init__(self,l=38,w=38,h=43,water=40,lamp=50,sand=2):
        self.length=l
        self.width=w
        self.height=h
        self.lampZreal=lamp
        self.water=water
        self.waterZ=water-self.lampZreal
        self.sand=sand
        self.sandZ=sand-self.lampZreal
        self.set_limits()
        self.set_planes()

    def set_lamp(self,height):
        self.lampZreal=height
        self.waterZ=self.water-self.lampZreal
        self.sandZ=self.sand-self.lampZreal
        print self.lampZreal,self.waterZ,self.sandZ
        self.set_limits()
        self.set_planes()

    def set_limits(self):
        #coordinate system has origin in center of tank
        self.limitX=self.length/2.
        self.limitY=self.width/2.
        self.limitZ=self.sandZ/2.
        self.limits=np.array([self.limitX,self.limitY,self.limitZ])

    def set_planes(self):
        '''
        define the bottom and sides of the tank as planes, so five planes
        planes are defined as [normal_vector,point1,point2]
        n(p1-p2)=n1*(p11-p21)+n2(p12-p22)+n3(p13-p23)
           2
           _ 
        1 |_|3  bottom 5
           4
        '''
        self.planes=[]
        n_xlimit=np.array([1,0,0])
        n_ylimit=np.array([0,1,0])
        n_zlimit=np.array([0,0,1])

        p1_1=np.array([-1*self.limitX,5.,5.])
        plane1=[-1.*n_xlimit,p1_1]

        p2_1=np.array([3.,self.limitY,3.])
        plane2=[-1.*n_ylimit,p2_1]

        p3_1=np.array([self.limitX,5.,5.])
        plane3=[-1.*n_xlimit,p3_1]

        p4_1=np.array([3.,-self.limitY,3.])
        plane4=[-1.*n_ylimit,p4_1]

        p5_1=np.array([1.,1.,self.limitZ])
        plane5=[1.*n_zlimit,p5_1]

        self.planes=[plane1,plane2,plane3,plane4,plane5]
        
        
def calc_light(dim,fnRb,fnWhite):

    '''
    script that sets up the required objects and does the plotting of the light distribution at the specified levels
    Requires as input the files with the led properties

    Currently messed up change to config file style at some point

    Add a GUI for ease of use 
    '''
    
    rbLed=led(filename=fnRb,wavelength=460,power=3.)
    wLed=led(filename=fnWhite,wavelength=520,power=3.)
    
    nano60=tank()
    nano60.set_lamp(nano60.water+20)
    x=[-8,-6,-4,-2,0,2,4,6,8,-8,-6,-4,-2,0,2,4,6,8]
    #x=-1.*np.array([-8.,-8.,-8.,-6.,-6.,-6.,-4.,-4.,-4.,4.,4.,4.,6.,6.,6.,8.,8.,8.])
    #x=np.array([-14,-10.,-6.,-2.,2,6.,10.,14,-14.,-10.,-6.,-2.,2,6.,10,14,-14,-10.,-6.,-2.,2,6.,10,14])
    #y=np.array([-11,-11.,-11.,-8.,-8.,-8.,-4.,-4.,-4.,0.,0.,0.,4.,4.,4.,8.,8.,8.,12.,12.,12.,-14.,-14.,-14.])
    #y=[-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,0.,0.,0.,0.,0.,0.,0.,0.,10.,10.,10.,10.,10.,10.,10.,10.]
    #x=[-8.,0.,8.,-8.,0.,8.,-8.,0.,8.,-8.,0.,8.,-8.,0.,8.,-8.,0.,8.,-8.,0.,8.,-8.,0.,8.]
    #y1=[-5.,-5.,-5.,-5.,-5.,-5.,-5.,-5.,-5.]
    y=-1.*np.array([-5.,-5.,-5.,-5.,-5.,-5.,-5.,-5.,-5.,5.,5.,5.,5.,5.,5.,5.,5.,5.])
    #y=[-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,10.,10.,10.,10.,10.,10.,10.,10.,10.]
    #y=[-2.,0.,2,-2.,0.,2,-2.,0.,2,-2.,0.,2,-2.,0.,2,-2.,0.,2,]
    #angle=np.array([115,115,115,115,115,115,115,115,115])*0.5*np.pi/180.
    angle=np.ones(18)*180*0.5*np.pi/180.
    #x=[18.]
    #y=[0.]
    #angle=np.array([90])*0.5*np.pi/180.
    z=0.0
    leds=[rbLed,wLed,rbLed,wLed,rbLed,wLed,rbLed,wLed,rbLed,wLed,rbLed,wLed,rbLed,wLed,rbLed,wLed,rbLed,wLed]
    lamps=lamp(x=y,y=x,z=z,led=leds,angle=angle)
    #print lamp.x
    numBin=1.4*nano60.length*10
    sizex=np.linspace(-1.2*nano60.length/2.,1.2*nano60.length/2.,numBin)
    sizey=np.linspace(-1.2*nano60.width/2.,1.2*nano60.width/2.,numBin)
    #print sizex

    zlevel=np.linspace(nano60.waterZ,nano60.sandZ,int(abs(nano60.waterZ-nano60.sandZ)/5.))
    #zlevel=[nano60.waterZ]
    #print zlevel
    fig=P.figure()
    counter=0

    lamps.set_ray_resolution(1.,180.)
    
    for z in zlevel:
        print 'calculating level %i at waterdepth %f'%(counter,abs(z-nano60.waterZ))
        #lamp.get_overlay_intens(z,sizex,sizey)
        #ax=fig.add_axes([0.175,0.15,0.8,0.8])
        X,Y,H=lamps.get_intensity_level(nano60.planes,z,nano60.limits)
        if counter==0:
            #print X,Y,H
            colorlevel=np.linspace(np.array(H).min(),np.array(H).max(),20)

        ax=fig.add_subplot(3,3,counter+1)
        counter+=1
        ax.add_patch(patches.Rectangle((-nano60.length/2.,-nano60.width/2.), width=nano60.length, height=nano60.width,linewidth=2,fill=False))
        #ax=plt.contourf(sizex,sizey,lamp.overIntens,levels=colorlevel)
        print 'Size x: %s  min,max: %s, %s'%(np.size(X),X.min(),X.max())
        print 'Size y: %s  min,max: %s, %s'%(np.size(Y),Y.min(),Y.max())
        print 'Size w: %s  min,max: %s, %s'%(np.size(H),H.min(),H.max())
        #hist=np.histogram2d(X,Y,weights=H,bins=[np.linspace(-6.3,6.3,20),np.linspace(-6.3,6.3,20)])
        #print hist
        #ax=plt.imshow(hist, interpolation='nearest')
        tsHistMatplot,xedges,yedges=np.histogram2d(X,Y,weights=H,bins=[np.linspace(-25,25,50),np.linspace(-25,25,50)])
        tsHistMatplot.shape,xedges.shape,yedges.shape
        #print xedges
        #print yedges
        extension=[yedges[0],yedges[-1],xedges[0],xedges[-1]]
    
        tsImgMplot=ax.imshow(tsHistMatplot,origin="lower",extent=extension,interpolation='nearest',vmin=0,vmax=0.2)#lamp.intensRay*500)
        #lower,upper=ax.xlim()
        #ax=plt.pcolormesh(X, Y, H)
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(tsImgMplot,cax=cax)
    fig.show()


    var=raw_input("Press enter to exit")
    
parser=argparse.ArgumentParser()
parser.add_argument("--dimension",help="bxtxh Breite Tiefe Hoehe",type=str,default="38x38x43")
parser.add_argument("--fnWhite",help="Filename containing intensity data for white led",type=str,default="cree_xpe2_cold_white_fine.csv")
parser.add_argument("--fnRb",help="Filename containing intensity data for royal blue led",type=str,default="cree_xte_royal_blue_fine.csv")
args=parser.parse_args()

calc_light(args.dimension,args.fnWhite,args.fnRb)

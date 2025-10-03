# galcrash_app.py - Full simulation GUI with star disks and interactive controls
from GalaxGui import GalaxySimGUI 
import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication)
import pyqtgraph as pg
from pyqtgraph import mkPen, mkBrush
class Galaxy:
    def __init__(self, mass, pos, vel, ahalo, vhalo, rthalo):
        self.mass = mass
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.ahalo = ahalo
        self.vhalo = vhalo
        self.rthalo = rthalo
        self.acc = np.full((3,1),0.0)
    
    def setposvel(self, pos, vel):
        self.pos = pos
        self.vel = vel

    def scalemass(self, massFact):
        self.mass = self.mass * massFact
        self.vhalo = self.vhalo * massFact * 0.25
        self.ahalo = self.ahalo * massFact * 0.1
        a2 = -self.mass/(self.vhalo ** 2)
        a1 = -2 * self.ahalo * self.mass/(self.vhalo ** 2 )
        a0 = -self.mass * (self.ahalo**2)/(self.vhalo**2)
        q = a1/3.0 - (a2**2)/9.0
        r = (a1 * a2 - 3.0 * a0)/6.0 - (a2**3)/27.0

        s1 = (r + np.sqrt((q**3) + (r**2)))**0.3333
        s2 = (r - np.sqrt((q**3) + (r**2)))**0.3333
        self.rthalo = (s1 +s2)-a2/3.0

    def MoveGalaxy(self, dt):
        newpos = self.pos + self.vel * dt + 0.5 * self.acc * dt**2
        newvel = self.vel + self.acc * dt
        self.pos = newpos
        self.vel = newvel
    def update_acc(self, opos):
        G = 1.0
        dr = opos - self.pos
        r = np.sqrt(np.sum(dr**2, axis=0))
        Accmag = -(G * self.InMass(r)) / (r**2)
        calcacc = (Accmag * dr) / r
       
        return calcacc
    def Potential(self, pos):
        G = 1.0
        dr = pos - self.pos
        r = np.sqrt(np.sum(dr**2, axis=0))
        pot = G * self.InMass(r) / r
        
        return pot    
    def InMass(self, r):
        indices = r < self.rthalo
        intmass = np.full(r.shape,0.0)
        if intmass[indices].shape != (0,):
            intmass[indices] = (self.vhalo**2)*(r[indices]**3) / ((self.ahalo+r[indices])**2)
        if intmass[~indices].shape != (0,):
            intmass[~indices] = self.mass
    
        return intmass
    def Density(self, pos):
        r_inner = pos* 0.999
        r_outer = pos * 1.001
        m_inner = self.InMass(r_inner)
        m_outer = self.InMass(r_outer)
        dm = m_outer - m_inner
        vol = (4/3)*np.pi*((r_outer**3) - (r_inner**3))
        density = dm / vol

        return density
    
    def DynFric(self, pmass, ppos, pvel):
        G = 1.0
        ingamma = 3
        dv = pvel - self.vel
        v = np.linalg.norm(dv)
        dr = ppos - self.pos
        r = np.linalg.norm(dr)
        g_rho = self.Density(r)
        fric_mag = 4.0 * np.pi * G * ingamma * pmass * g_rho * v/((1+v)**3)
        friction = dv / v * fric_mag
        
        return friction
class Star(Galaxy):

    def __init__(self, mass, ahalo, vhalo, rthalo, pos, vel, dsize, theta, phi, n):
        super().__init__(mass, pos, vel, ahalo, vhalo, rthalo)
        self.dsize = dsize
        self.theta = theta
        self.phi = phi
        self.n = n
        
        self.spos = np.full((3,self.n),0.0)
        self.svel = np.full((3,self.n),0.0)
        self.sacc = np.full((3,self.n),0.0)

    def move(self, dt):
        nspos = self.spos + self.svel * dt + 0.5 * self.sacc * dt**2
        nsvel = self.svel + self.sacc * dt
        self.spos = nspos
        self.svel = nsvel
    def InitStar(self):
        cosphi = np.cos(self.phi)
        sinphi = np.sin(self.phi)
        costheta = np.cos(self.theta)
        sintheta = np.sin(self.theta)

        for i in range(self.n):
            bad = True
            while bad:
                xtry = self.dsize * (1. - 2. * np.random.random())
                ytry = self.dsize * (1. - 2. * np.random.random())
                rtry = np.sqrt(xtry**2 + ytry**2)
                if rtry < self.dsize:
                    bad = False
                
            ztry = 0.0
            xrot = xtry * cosphi + ytry * sinphi * costheta + ztry * sinphi * sintheta
            yrot = -xtry * sinphi + ytry * cosphi * costheta + ztry * cosphi * sintheta
            zrot = -ytry * sintheta + ztry * costheta
            rot = np.array([xrot, yrot, zrot])
            self.spos[:,i] = rot + self.pos.reshape(-1)

            vcirc = np.sqrt(self.InMass(rtry) / rtry)

            vxtry = -vcirc * ytry / rtry
            vytry = vcirc * xtry / rtry
            vztry = 0.0
            
            vxrot = vxtry * cosphi + vytry * sinphi * costheta + vztry * sinphi * sintheta
            vyrot = -vxtry * sinphi + vytry * cosphi * costheta + vztry * cosphi * sintheta
            vzrot = -vytry * sintheta + vztry * costheta

            vrot = np.array([vxrot, vyrot, vzrot])
            self.svel[:,i] = vrot + self.vel.reshape(-1)
            self.sacc = np.full((1,3),0.0)

    def scalemass(self, massFact):
        self.dsize = self.dsize * np.sqrt(massFact)
        super().scalemass(massFact)
class orbit:
    def __init__(self, energy, rp, tp, ecc, m1, m2, pos1, pos2, vel1, vel2,):
        self.energy = energy
        self.rp = rp
        self.tp = tp
        self.m1 = m1
        self.m2 = m2
        self.ecc = ecc
        self.pos1 = pos1
        self.pos2 = pos2
        self.vel1 = vel1
        self.vel2 = vel2
        self.initorbit()

    def initorbit(self):
        mu = self.m1 + self.m2
        p = 2 * self.rp
        nhat = np.sqrt(mu/p**3)
        cots = 3.0 * nhat * self.tp
        s = np.arctan(1/cots)
        cottheta = (1. /np.tan(s/2.))**(1/3)
        theta = np.arctan(1./cottheta)
        tanfon2 = 2./np.tan(theta*2.)
        r = (p/2.)*(1+tanfon2**2)

        vel = np.sqrt(2.*mu/r)
        sinsqphi = p/(2.*r)
        phi = np.arcsin(np.sqrt(sinsqphi))
        f = 2. * np.arctan(tanfon2)
        xc = -r * np.cos(f)
        yc = r * np.sin(f)
        vxc = vel * np.cos(f+phi) 
        vyc = -vel * np.sin(f+phi)
        xcom = self.m2*xc/(self.m1+self.m2)
        ycom = self.m2*yc/(self.m1+self.m2)
        vxcom = self.m2*vxc/(self.m1+self.m2)
        vycom = self.m2*vyc/(self.m1+self.m2)   

        self.pos1 = np.array([[-xcom], [-ycom], [0.0]])
        self.pos2 = np.array([[xc-xcom], [yc-ycom], [0.0]])   
        self.vel1 = np.array([[-vxcom], [-vycom], [0.0]])
        self.vel2 = np.array([[vxc-vxcom], [vyc-vycom], [0.0]])

class GUI(GalaxySimGUI):
    def __init__(self):
        super(GUI,self).__init__()
        self.show()

        self.btn_start.clicked.connect(self.MakeGalaxy)
        self.btn_start.clicked.connect(self.MakeOrbit)
        self.btn_start.clicked.connect(self.update_plot)

        self.btn_reset.clicked.connect(self.reset)
        
        self.time = 0.0

        self.Dist=[]
        self.Vel=[]
        self.Time = []

        self.graphicsView.setXRange(-20, 20)
        self.graphicsView.setYRange(-20, 20)
        self.graphicsView.plotItem.hideAxis('bottom')
        self.graphicsView.plotItem.hideAxis('left')

        self.velplot.setXRange(0, 2000)
        self.velplot.setYRange(0, 1000)  
        self.velplot.plotItem.setLabel('bottom', 'Time[Myrs]')
        self.velplot.plotItem.setLabel('left', 'Relative Velocity[km/s]')


        self.pen1 = mkPen(color=(255, 0, 0))
        self.brush1 = mkBrush('r')
        self.pen2 = mkPen(color=(155, 41, 163))
        self.pen3 = mkPen(color=(255, 255, 255))
        self.pen4 = mkPen(color=(85, 170, 0))
        self.pen5 = mkPen(color=(184, 6, 0))
        self.pen6 = mkPen(color=(0, 255, 0))
        self.brush2 = mkBrush('g')

        self.velplot.plotItem.getAxis('bottom').setPen(self.pen3)
        self.velplot.plotItem.getAxis('left').setPen(self.pen2)

        self.distplot.setXRange(0, 2000)
        self.distplot.setYRange(0, 200)
        
        self.distplot.plotItem.hideAxis('bottom')
        self.distplot.plotItem.setLabel('left', 'Galaxy Separation[kpc]')
        self.distplot.plotItem.getAxis('left').setPen(self.pen4)

        self.distlabel.setStyleSheet("color:rgb( 85, 170, 0)")
        self.distlabel.setText("GalaxySeparation: ")
        self.vellabel.setStyleSheet("color:rgb(155, 41, 163)")
        self.vellabel.setText("RelativeVelocity: ")
        self.timelabel.setStyleSheet("color:white")
        self.timelabel.setText("Time: ")

    def MakeGalaxy(self):
        galmass = 4.8
        ahalo = 0.1
        vhalo = 1.0
        rthalo = 5.0
        pos = np.full((3,1),0.0)
        vel = np.full((3,1),0.0)
        disksize = 2.5

        galtheta = int(self.green_theta.value())
        galphi = int(self.green_phi.value())
        comptheta = int(self.red_theta.value())
        compphi = int(self.red_phi.value())

        galn = int(0.5*self.star_count.value())
        compn = int(0.5*self.star_count.value())
        self.gal = Star(galmass, ahalo, vhalo, rthalo, pos, vel, disksize, galtheta, galphi, galn)
        self.comp = Star(galmass, ahalo, vhalo, rthalo, pos, vel, disksize, comptheta, compphi, compn)

        if self.cb_big_halos.isChecked():
            self.gal.rthalo = 20.0
            self.gal.mass = (self.gal.vhalo**2)*(self.gal.rthalo**3) / ((self.gal.ahalo+self.gal.rthalo)**2)
            self.comp.rthalo = self.comp.rthalo * 4.0
            self.comp.mass = (self.comp.vhalo**2)*(self.comp.rthalo**3) / ((self.comp.ahalo+self.comp.rthalo)**2)
        else:
            self.gal.rthalo = 5.0
            self.gal.mass = 4.8


        massrat = float(self.red_mass.value())
        self.comp.scalemass(massrat)

    def MakeOrbit(self):
        energy = 0.0
        rp = 3.0
        ecce=1.0
        tp = float(self.peri.value())

        self.crash = orbit(energy, rp, tp, ecce, self.gal.mass, self.comp.mass, self.gal.pos, self.comp.pos, self.gal.vel, self.comp.vel)
        self.gal.setposvel(self.crash.pos1, self.crash.vel1)
        self.comp.setposvel(self.crash.pos2, self.crash.vel2)

        if self.cb_red_centered.isChecked():
            self.graphicsView.setXRange(self.gal.pos[0,0]-20, self.gal.pos[0,0]+20)
            self.graphicsView.setYRange(self.gal.pos[1,0]-20, self.gal.pos[1,0]+20)
            
        Gal = pg.ScatterPlotItem(brush=self.brush2, size = 6)
        Gal.setData([self.gal.pos[0,0]], [self.gal.pos[1,0]])
        
        Comp = pg.ScatterPlotItem(brush=self.brush1, size = 6)
        Comp.setData([self.comp.pos[0,0]], [self.comp.pos[1,0]])
        
        self.graphicsView.addItem(Gal)
        self.graphicsView.addItem(Comp)
        self.gal.InitStar() 
        self.comp.InitStar()
        self.graphicsView.plotItem.plot(self.gal.spos[0,:], self.gal.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen6)
        self.graphicsView.plotItem.plot(self.comp.spos[0,:], self.comp.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen1, fill = True)
        dist = 3.5 * np.linalg.norm(self.gal.pos - self.comp.pos)
        self.distlabel.setText("Galaxy Separation: {:.2f} kpc".format(dist))   
        vel = 250 * np.linalg.norm(self.gal.vel - self.comp.vel) 
        self.vellabel.setText("Relative Velocity: {:.2f} km/s".format(vel))
        self.timelabel.setText("Time: 0.0 Myrs")

    def update_plot(self):
        self.btn_reset.clicked = True
        dtime = self.sim_spd.value()
        if self.btn_reset.clicked:
            while self.btn_reset.clicked:
                self.graphicsView.plotItem.clear()
                self.distplot.plotItem.clear()
                self.velplot.plotItem.clear()

                self.gal.acc = self.comp.update_acc(self.gal.pos)
                self.comp.acc = self.gal.update_acc(self.comp.pos)
                dist = 3.5 * np.linalg.norm(self.gal.pos - self.comp.pos)

                if self.cb_friction.isChecked():
                    self.gal.acc = self.gal.acc + self.comp.DynFric(self.gal.InMass(dist/3.5), self.gal.pos, self.gal.vel)
                    self.comp.acc = self.comp.acc + self.gal.DynFric(self.comp.InMass(dist/3.5), self.comp.pos, self.comp.vel)

                comacc=((self.gal.acc*self.gal.mass) + (self.comp.acc*self.comp.mass))/(self.gal.mass+self.comp.mass)
                self.gal.acc = self.gal.acc - comacc
                self.comp.acc = self.comp.acc - comacc
                self.gal.sacc = self.gal.update_acc(self.gal.spos)+self.comp.update_acc(self.gal.spos)
                self.comp.sacc = self.comp.update_acc(self.comp.spos)+self.gal.update_acc(self.comp.spos)

                self.gal.MoveGalaxy(dtime)
                self.comp.MoveGalaxy(dtime)
                self.gal.move(dtime)
                self.comp.move(dtime)
                
                dist = 3.5 * np.linalg.norm(self.gal.pos - self.comp.pos)
                vel = 250 * np.linalg.norm(self.gal.vel - self.comp.vel)

                if self.cb_red_centered.isChecked():
                    self.graphicsView.setXRange(self.gal.pos[0,0]-20, self.gal.pos[0,0]+20)
                    self.graphicsView.setYRange(self.gal.pos[1,0]-20, self.gal.pos[1,0]+20)
                if self.cb_green_centered.isChecked():
                    self.graphicsView.setXRange(self.comp.pos[0,0]-20, self.comp.pos[0,0]+20)
                    self.graphicsView.setYRange(self.comp.pos[1,0]-20, self.comp.pos[1,0]+20)
                Gal = pg.ScatterPlotItem(brush=self.brush2, size = 6, symbol ='o')
                Gal.setData([self.gal.pos[0,0]], [self.gal.pos[1,0]])
                Comp = pg.ScatterPlotItem(brush=self.brush1, size = 6, symbol = 'o')
                Comp.setData([self.comp.pos[0,0]], [self.comp.pos[1,0]])
                self.graphicsView.addItem(Gal)
                self.graphicsView.addItem(Comp)
                self.graphicsView.plotItem.plot(self.gal.spos[0,:], self.gal.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen6)
                self.graphicsView.plotItem.plot(self.comp.spos[0,:], self.comp.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen1, fill = True)
                self.Dist.append(dist)
                self.Vel.append(vel)
                self.Time.append(self.time * 12)
                self.distlabel.setText("Galaxy Separation: {:.2f} kpc".format(dist))   
                self.vellabel.setText("Relative Velocity: {:.2f} km/s".format(vel))
                self.timelabel.setText("Time: {:.2f} Myrs".format(self.time * 12))

                QApplication.processEvents()

                self.time += dtime

        else:
            None

    def reset(self):
        self.btn_pause = False
        self.btn_reset.clicked = False
        self.graphicsView.plotItem.clear()
        self.graphicsView.setXRange(-20, 20)
        self.graphicsView.setYRange(-20, 20)
        self.distplot.plotItem.clear()
        self.velplot.plotItem.clear()
        self.gal = None
        self.comp = None
        self.time = 0.0
        self.Dist = []
        self.Vel = []
        self.Time = []
        self.distlabel.setText("Galaxy Separation: 0.0 kpc")
        self.vellabel.setText("Relative Velocity: 0.0 km/s")
        self.timelabel.setText("Time: 0.0 Myrs")

    def pause(self):
        self.btn_pause = False

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = GUI() 
    window.show()
    sys.exit(app.exec_())


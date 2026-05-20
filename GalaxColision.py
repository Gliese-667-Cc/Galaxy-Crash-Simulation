# galcrash_app.py - Full simulation GUI with star disks and interactive controls
from GalaxGui import GalaxySimGUI 
import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication)
from PyQt5.QtCore import QTimer
import pyqtgraph as pg
from pyqtgraph import mkPen, mkBrush
import numpy as np

class Galaxy:
    def __init__(self, mass, a_scale, pos, vel, disksize, theta, phi, n):
        self.mass = float(mass)
        self.a = float(a_scale) # The Plummer scale length (softens the core)
        
        self.pos = np.array(pos, dtype=float).reshape(3, 1)
        self.vel = np.array(vel, dtype=float).reshape(3, 1)
        
        self.n = n
        self.spos = np.zeros((3, self.n))
        self.svel = np.zeros((3, self.n))
        self.sacc = np.zeros((3, self.n))
        
        self._init_disk(disksize, theta, phi)

    def setposvel(self, pos, vel):
        """Moves the center AND shifts all stars seamlessly to match."""
        new_pos = np.array(pos).reshape(3,1)
        new_vel = np.array(vel).reshape(3,1)
        
        self.spos += (new_pos - self.pos)
        self.svel += (new_vel - self.vel)
        self.pos = new_pos
        self.vel = new_vel

    def _init_disk(self, R_max, theta_deg, phi_deg):
        """Generates a perfectly stable disk using Plummer mechanics."""
        # 1. Uniform area distribution in 2D disk
        r = R_max * np.sqrt(np.random.random(self.n))
        theta_angle = 2.0 * np.pi * np.random.random(self.n)
        
        x = r * np.cos(theta_angle)
        y = r * np.sin(theta_angle)
        z = np.zeros(self.n)
        
        # 2. Perfect circular velocity for Plummer potential
        # v = sqrt(G * M * r^2 / (r^2 + a^2)^(3/2))
        v_circ = np.sqrt(self.mass * r**2 / (r**2 + self.a**2)**1.5)
        
        vx = -v_circ * np.sin(theta_angle)
        vy =  v_circ * np.cos(theta_angle)
        vz = np.zeros(self.n)
        
        # 3. 3D Rotation Matrices
        th = np.radians(theta_deg)
        ph = np.radians(phi_deg)
        
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(th), -np.sin(th)],
                       [0, np.sin(th),  np.cos(th)]])
        Rz = np.array([[np.cos(ph), -np.sin(ph), 0],
                       [np.sin(ph),  np.cos(ph), 0],
                       [0, 0, 1]])
        Rot = Rz @ Rx
        
        # 4. Apply rotation and translate to center
        self.spos = Rot @ np.vstack((x, y, z)) + self.pos
        self.svel = Rot @ np.vstack((vx, vy, vz)) + self.vel

    def get_gravity(self, target_pos):
        """Calculates acceleration exerted BY this galaxy on ANY target positions."""
        dr = target_pos - self.pos
        r2 = np.sum(dr**2, axis=0)
        
        # Plummer acceleration: a = -G * M * dr / (r^2 + a^2)^1.5
        denom = (r2 + self.a**2)**1.5
        return -self.mass * dr / denom
        
    def scalemass(self, massFact):
        if massFact == 1.0:
            return
        self.mass *= massFact
        self.a *= (massFact ** 0.333) # Scale length increases slightly with massclass orbit:

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
        self.btn_pause.clicked.connect(self.pause)
        self.btn_reset.clicked.connect(self.reset)
        
        self.time = 0.0

        self.Dist=[]
        self.Vel=[]
        self.Time = []

        # Control flags for simulation loop
        self.running = False
        self.paused = False

        self.graphicsView.setXRange(-20, 20)
        self.graphicsView.setYRange(-20, 20)
        self.graphicsView.plotItem.hideAxis('bottom')
        self.graphicsView.plotItem.hideAxis('left')

        self.velplot.setXRange(0, 2000)
        self.velplot.setYRange(0, 1000)  
        self.velplot.plotItem.setLabel('bottom', 'Time[Myrs]')
        self.velplot.plotItem.setLabel('left', 'Relative Velocity[km/s]')


        self.pen1 = mkPen(color=(255,0,0))
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

        # plotting curves for live update
        self.dist_curve = self.distplot.plot([], [], pen=self.pen4)
        self.vel_curve = self.velplot.plot([], [], pen=self.pen2)

        self.distlabel.setStyleSheet("color:rgb( 85, 170, 0)")
        self.distlabel.setText("GalaxySeparation: ")
        self.vellabel.setStyleSheet("color:rgb(155, 41, 163)")
        self.vellabel.setText("RelativeVelocity: ")
        self.timelabel.setStyleSheet("color:white")
        self.timelabel.setText("Time: ")

    def MakeGalaxy(self):
        galmass = 4.8
        a_scale = 1.5 
        disksize = 3.0
        pos = np.full((3,1),0.0)
        vel = np.full((3,1),0.0)

        red_theta = int(self.red_theta.value())
        red_phi = int(self.red_phi.value())
        green_theta = int(self.green_theta.value())
        green_phi = int(self.green_phi.value())

        galn = int(0.5*self.star_count.value())
        compn = int(0.5*self.star_count.value())
        
        self.gal = Galaxy(galmass, a_scale, pos, vel, disksize, red_theta, red_phi, galn)
        self.comp = Galaxy(galmass, a_scale, pos, vel, disksize, green_theta, green_phi, compn)

        massrat = float(self.red_mass.value())
        self.gal.scalemass(massrat)

    def MakeOrbit(self):
        energy = 0.0
        rp = 3.0
        ecce=0.5
        tp = float(self.peri.value())

        self.crash = orbit(energy, rp, tp, ecce, self.gal.mass, self.comp.mass, self.gal.pos, self.comp.pos, self.gal.vel, self.comp.vel)
        self.gal.setposvel(self.crash.pos1, self.crash.vel1)
        self.comp.setposvel(self.crash.pos2, self.crash.vel2)

        if self.cb_red_centered.isChecked():
            self.graphicsView.setXRange(self.gal.pos[0,0]-20, self.gal.pos[0,0]+20)
            self.graphicsView.setYRange(self.gal.pos[1,0]-20, self.gal.pos[1,0]+20)
            
        Gal = pg.ScatterPlotItem(brush=self.brush1, size = 6)
        Gal.setData([self.gal.pos[0,0]], [self.gal.pos[1,0]])
        
        Comp = pg.ScatterPlotItem(brush=self.brush2, size = 6)
        Comp.setData([self.comp.pos[0,0]], [self.comp.pos[1,0]])
        
        self.graphicsView.addItem(Gal)
        self.graphicsView.addItem(Comp)
        self.gal.InitStar() 
        self.comp.InitStar()
        self.graphicsView.plotItem.plot(self.gal.spos[0,:], self.gal.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen1)
        self.graphicsView.plotItem.plot(self.comp.spos[0,:], self.comp.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen6)
        dist = 3.5 * np.linalg.norm(self.gal.pos - self.comp.pos)
        self.distlabel.setText("Galaxy Separation: {:.2f} kpc".format(dist))   
        vel = 250 * np.linalg.norm(self.gal.vel - self.comp.vel) 
        self.vellabel.setText("Relative Velocity: {:.2f} km/s".format(vel))
        self.timelabel.setText("Time: 0.0 Myrs")

    def update_plot(self):
        # start a QTimer-driven simulation loop (single-step per timeout)
        if getattr(self, 'gal', None) is None or getattr(self, 'comp', None) is None:
            return

        # create timer if needed
        if not hasattr(self, 'timer'):
            self.timer = QTimer(self)
            self.timer.timeout.connect(self.step)

        # set running state and start timer (interval in ms)
        self.running = True
        self.paused = False
        interval = 16  # target ~60 FPS; step size uses sim_spd value
        self.timer.start(interval)

    def step(self):
        if not getattr(self, 'running', False) or getattr(self, 'paused', False):
            return
        if getattr(self, 'gal', None) is None or getattr(self, 'comp', None) is None:
            return

        # 1. Adaptive timestep calculation
        base_dt = float(self.sim_spd.value())
        dist_now = 3.5 * np.linalg.norm(self.gal.pos - self.comp.pos)
        scale = max(0.02, min(1.0, dist_now / 5.0))
        dtime = base_dt * scale

        # --- SYMPLECTIC STEP 1: Half-step position update ---
        self.gal.pos += self.gal.vel * (0.5 * dtime)
        self.comp.pos += self.comp.vel * (0.5 * dtime)
        self.gal.spos += self.gal.svel * (0.5 * dtime)
        self.comp.spos += self.comp.svel * (0.5 * dtime)

        # --- SYMPLECTIC STEP 2: Accelerations at half-step positions ---
        a_g = self.comp.get_gravity(self.gal.pos)
        a_c = self.gal.get_gravity(self.comp.pos)

        # Dynamically add friction if enabled
        if self.cb_friction.isChecked():
            v_rel = self.gal.vel - self.comp.vel
            v_rel_mag = np.linalg.norm(v_rel)
            if v_rel_mag > 0:
                friction_strength = 0.1 * (self.gal.mass * self.comp.mass) / (dist_now**2 + 0.1)
                friction_acc = -friction_strength * (v_rel / v_rel_mag)
                a_g += friction_acc / self.gal.mass
                a_c -= friction_acc / self.comp.mass
        # Per-star accelerations
        self.gal.sacc = self.gal.get_gravity(self.gal.spos) + self.comp.get_gravity(self.gal.spos)
        self.comp.sacc = self.comp.get_gravity(self.comp.spos) + self.gal.get_gravity(self.comp.spos)

        # --- SYMPLECTIC STEP 3: Full-step velocity update ---
        self.gal.vel += a_g * dtime
        self.comp.vel += a_c * dtime
        self.gal.svel += self.gal.sacc * dtime
        self.comp.svel += self.comp.sacc * dtime

        # --- SYMPLECTIC STEP 4: Final half-step position update ---
        self.gal.pos += self.gal.vel * (0.5 * dtime)
        self.comp.pos += self.comp.vel * (0.5 * dtime)
        self.gal.spos += self.gal.svel * (0.5 * dtime)
        self.comp.spos += self.comp.svel * (0.5 * dtime)

        # --- REDRAW LOGIC ---
        dist = 3.5 * np.linalg.norm(self.gal.pos - self.comp.pos)
        vel = 250 * np.linalg.norm(self.gal.vel - self.comp.vel)

        if self.cb_red_centered.isChecked():
            self.graphicsView.setXRange(self.gal.pos[0,0]-20, self.gal.pos[0,0]+20)
            self.graphicsView.setYRange(self.gal.pos[1,0]-20, self.gal.pos[1,0]+20)
        if self.cb_green_centered.isChecked():
            self.graphicsView.setXRange(self.comp.pos[0,0]-20, self.comp.pos[0,0]+20)
            self.graphicsView.setYRange(self.comp.pos[1,0]-20, self.comp.pos[1,0]+20)

        self.graphicsView.plotItem.clear()

        Gal = pg.ScatterPlotItem(brush=self.brush2, size=6, symbol='o')
        Gal.setData([self.gal.pos[0,0]], [self.gal.pos[1,0]])
        Comp = pg.ScatterPlotItem(brush=self.brush1, size=6, symbol='o')
        Comp.setData([self.comp.pos[0,0]], [self.comp.pos[1,0]])
        
        self.graphicsView.addItem(Gal)
        self.graphicsView.addItem(Comp)
        self.graphicsView.plotItem.plot(self.gal.spos[0,:], self.gal.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen6)
        self.graphicsView.plotItem.plot(self.comp.spos[0,:], self.comp.spos[1,:], pen=None, symbol='t1', symbolSize=3, symbolPen=self.pen1, fill=True)
        
        self.Dist.append(dist)
        self.Vel.append(vel)
        self.Time.append(self.time * 12)
        
        try:
            self.dist_curve.setData(self.Time, self.Dist)
            self.vel_curve.setData(self.Time, self.Vel)
        except Exception:
            pass

        self.distlabel.setText("Galaxy Separation: {:.2f} kpc".format(dist))
        self.vellabel.setText("Relative Velocity: {:.2f} km/s".format(vel))
        self.timelabel.setText("Time: {:.2f} Myrs".format(self.time * 12))

        self.time += dtime
    def reset(self):
        # stop simulation loop and clear state
        self.running = False
        self.paused = False
        if hasattr(self, 'timer') and self.timer.isActive():
            self.timer.stop()
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
        # toggle pause by stopping/starting the timer
        if hasattr(self, 'timer'):
            if self.timer.isActive():
                self.timer.stop()
                self.paused = True
                try:
                    self.btn_pause.setText('Resume')
                except Exception:
                    pass
            else:
                self.timer.start()
                self.paused = False
                try:
                    self.btn_pause.setText('Pause')
                except Exception:
                    pass
        else:
            # fallback: toggle paused flag
            self.paused = not getattr(self, 'paused', False)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = GUI() 
    window.show()
    sys.exit(app.exec_())


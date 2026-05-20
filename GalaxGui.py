from PyQt5.QtWidgets import (QApplication, QWidget, QLabel, QSpinBox, QDoubleSpinBox, QCheckBox, QPushButton, QVBoxLayout, QHBoxLayout)
import pyqtgraph as pg
import sys



class GalaxySimGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Galaxy Collision Simulator")
        self.setGeometry(100, 100, 1200, 800)
        self.setStyleSheet("background-color: black; color: white;")
        self.initUI()

    def initUI(self):
        main_layout = QHBoxLayout(self)

        # Left control panel
        control_panel = QVBoxLayout()

        # Red Galaxy Controls
        control_panel.addWidget(QLabel("Red Theta:"))
        self.red_theta = QSpinBox()
        self.red_theta.setRange(0, 360)
        control_panel.addWidget(self.red_theta)

        control_panel.addWidget(QLabel("Red Phi:"))
        self.red_phi = QSpinBox()
        self.red_phi.setRange(0, 360)
        control_panel.addWidget(self.red_phi)

        # Green Galaxy Controls
        control_panel.addWidget(QLabel("Green Theta:"))
        self.green_theta = QSpinBox()
        self.green_theta.setRange(0, 360)
        control_panel.addWidget(self.green_theta)

        control_panel.addWidget(QLabel("Green Phi:"))
        self.green_phi = QSpinBox()
        self.green_phi.setRange(0, 360)
        control_panel.addWidget(self.green_phi)

        # Parameters
        control_panel.addWidget(QLabel("Peri [kpc]:"))
        self.peri = QDoubleSpinBox()
        self.peri.setDecimals(1)
        self.peri.setValue(10.5)
        control_panel.addWidget(self.peri)

        control_panel.addWidget(QLabel("Green Galaxy Mass (w.r.t Red):"))
        self.red_mass = QDoubleSpinBox()
        self.red_mass.setDecimals(2)
        self.red_mass.setValue(1.00)
        control_panel.addWidget(self.red_mass)

        control_panel.addWidget(QLabel("Number of Stars:"))
        self.star_count = QSpinBox()
        self.star_count.setRange(1, 10000)
        self.star_count.setValue(250)
        control_panel.addWidget(self.star_count)

        control_panel.addWidget(QLabel("Simulation Speed:"))
        self.sim_spd = QDoubleSpinBox()
        self.sim_spd.setRange(0.01, 100.0)
        self.sim_spd.setValue(0.01)
        control_panel.addWidget(self.sim_spd)


        # Checkboxes
        self.cb_friction = QCheckBox("Friction")
        self.cb_big_halos = QCheckBox("Big Halos")
        self.cb_red_centered = QCheckBox("Red Centered")
        self.cb_green_centered = QCheckBox("Green Centered")
        control_panel.addWidget(self.cb_friction)
        control_panel.addWidget(self.cb_big_halos)
        control_panel.addWidget(self.cb_red_centered)
        control_panel.addWidget(self.cb_green_centered)

        # Buttons
        self.btn_start = QPushButton("Start")
        self.btn_pause = QPushButton("Pause")
        self.btn_reset = QPushButton("Reset")
        control_panel.addWidget(self.btn_start)
        control_panel.addWidget(self.btn_pause)
        control_panel.addWidget(self.btn_reset)

        main_layout.addLayout(control_panel)

        # Right graphics and labels
        right_panel = QVBoxLayout()
        self.top_label = QLabel("TextLabel")
        right_panel.addWidget(self.top_label)

        # Use PyQtGraph for galaxy view
        self.graphicsView = pg.PlotWidget()
        self.graphicsView.setMinimumHeight(500)
        self.graphicsView.setBackground('k')
        self.graphicsView.plotItem.hideAxis('left')
        self.graphicsView.plotItem.hideAxis('bottom')
        right_panel.addWidget(self.graphicsView)

        # Stack distPlot over velPlot
        stacked_plots = QVBoxLayout()

        # Distance plot
        self.distplot = pg.PlotWidget()
        self.distplot.setBackground('k')
        self.distplot.plotItem.setLabel('left', 'Distance')
        self.distplot.plotItem.setLabel('bottom', 'Time')
        stacked_plots.addWidget(self.distplot)

        # Velocity plot
        self.velplot = pg.PlotWidget()
        self.velplot.setBackground('k')
        self.velplot.plotItem.setLabel('left', 'Velocity')
        self.velplot.plotItem.setLabel('bottom', 'Time')
        stacked_plots.addWidget(self.velplot)

        right_panel.addLayout(stacked_plots)

        # Labels for distance and time
        labels_layout = QHBoxLayout()
        self.distlabel = QLabel("Distance: 0.0 kpc")
        self.timelabel = QLabel("Time: 0.0 Myr")
        self.vellabel = QLabel("Velocity: 0.0 km/s")
        labels_layout.addWidget(self.distlabel)
        labels_layout.addWidget(self.timelabel)
        labels_layout.addWidget(self.vellabel)
        right_panel.addLayout(labels_layout)

        main_layout.addLayout(right_panel)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = GalaxySimGUI()
    window.show()
    sys.exit(app.exec_())

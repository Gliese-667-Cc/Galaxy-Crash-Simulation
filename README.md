#Galactic Collision Simulator

##Abstract
This simulation models the collision of two galaxies, incorporating the complex 
effects of gravitational friction during star-star interactions alongside the 
physical influence of dark matter halos. Derived from established astrophysical 
research, the implementation utilizes a simplified model of galactic cores driven 
\by Newtonian mechanics to accurately visualize these massive cosmic events.
----------------------------------------------------------------------------------
##Current Implementation
The current iteration of the codebase is developed entirely in Python, prioritizing 
rapid prototyping and mathematical accuracy:

Physics & Computation: Leverages NumPy for highly efficient numerical operations, 
                       vector mathematics, and state updates.

Visualization: Utilizes PyQtGraph for real-time rendering, plotting, and interactive
               graphical output.
------------------------------------------------------------------------------------
##Vision & Future Roadmap
Future development is focused on significantly enhancing both the computational depth 
and rendering performance of the engine. The planned roadmap includes:

Physics Architecture Upgrade: Transition the computational framework from Newtonian 
                              mechanics to Lagrangian mechanics. By utilizing potential
                              vector fields, the engine will handle highly complex orbital
                              systems with greater accuracy and numerical stability while
                              optimizing code complexity.

Performance Migration: Port the core engine to C++ and utilize OpenGL for rendering, vastly
                        augmenting the simulation's performance, frame rate, and particle scale.

Interactive UI: Integrate imGUI to construct a real-time control panel. This slider-based 
                interface will enable users to dynamically adjust initial parameters, facilitating
                comprehensive, scenario-based studies of galactic collision dynamics.

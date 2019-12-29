# Computational Physics

This repository holds the code I've written and relevant physics to solving the problems for the Computational Physics class I took in Spring 2019 at Duke. The homework problems come from the textbook Computational Physics by Mark Newman. The final project was one of my own choosing.


### Notes on Final Project
The final project solves Navier-Stokes equations (2nd order partial differential equations) using a FTCS(forward-time centered-space) algorithm for laminar fluid flow through a channel. It does this for  unobstructed channels and for two kinds of obstructions.

For the final project, the main functions used to solve the systems of equations are in

    channel_flow.py

This is called by

    Nair Final Project.ipynb

where most of the graphs and numerical analysis are done.

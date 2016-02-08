# Explicit Runge--Kutta Methods in C++

This is a collection of a few of the more common Runge-Kutta integration schemes.

I've hard-coded a few of the simple schemes (Euler, Mid-Point, "Classical" Runge--Kutta). 
In addition, I've included code for computing a general-form Runge--Kutta method from its Butcher table.

## Code Structure:

This code should compile nicely using g++ on unix platforms by calling "make" in terminal.

I have a driver-file called main.cpp that shows how to call each method, and provides two simple examples for dynamics functions.

For now, I just have the integration methods write to a file. 
There is a Matlab script that can be used to plot the solution and compare it to ode45 (Matlab's variable-step version of RK45).

## Hard-Coded methods:
- Euler's Method 
- MidPoint Method
- Classical Runge--Kutta

## Butcher-Table methods:
- RK2 (same as the Mid-Point method)
- RK4A (same as the classical Runge--Kutta)
- RK4B (4th-order "3/8 Rule")
- RK45 (Runge--Kutta--Fehlberg, 5th-order method)
- RK5 (5th-order Runge--Kutta, from paper by Fehlberg)
- RK10 (10th-order Runge--Kutta, from paper by Feagin)





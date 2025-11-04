# Pandemic Modeling with SEEIIR / SEIR ODEs

This project implements and explores a deterministic compartmental disease-spread model based on the SEEIIR (extended SEIR) framework, solving the model using custom ODE solvers and analyzing epidemic dynamics under changing contact rates and import of infected individuals.

The project simulates:

- Baseline infectious disease spread

- Time-varying contact/infection rate β(t)

- Effects of interventions and behavioral changes

- Multiple epidemic waves via imported cases

All models are written from scratch using a modular OOP structure with inheritance for model variants.
All data is based on data from Folkehelse instituttet. 

# Project Proposal

## Summary
We will parallelize an N-bodies problem with a Barnes-Hut Tree algorithm.

## Background
One potential use for an N-bodies algorithm is to model the interactions between atoms. Assuming that the model is correct, one could then run a simulation that predicts the equilibrium structure of a given system with high accuracy, which could be useful for predicting the structures of novel materials and substances.

The force at work between the atoms in this type of problem is electronic. Atoms consist of both positively and negatively charged particles, so they experience both an attraction to each other (between the oppositely charged particles) and a repulsion (between the same charged particles). The repulsion force tends to be much larger, but also fall off more quickly over large distances, resulting in atoms settling into equilibrium distances from each other where the attractive and repulsive forces are approximately equal and opposite.

So the overall algorithm for this type of simulation is as follows, for every timestep:
- Calculate the net force acting on every particle based on its position relative to all of the other particles
- Add in a random component based on the temperature of the system (molecules are always moving randomly)
- Calculate new positions by applying these forces over some discrete amount of time

## Challenge
As discussed in class, the main challenge with a problem like this is how to distribute the work in order to minimize the load imbalance and maximize the locality.

## Resources
We’ll use the GHC clusters to develop and test the code for this project.

We’re planning to use pthreads to implement our parallelism.

## Goals and Deliverables
75%: We will implement a working, parallel solution for our N-bodies problem using a Barnes-Hut type algorithm

100%: In addition to our solution working it will also have a complex algorithm for semi-static workload assignment

125%: On top of the algorithm for balancing the workload our solution will also optimize for locality using a method similar to the “cost zones” we discussed in lecture

## Schedule
3/21/22-3/23/22: create project proposal

3/28/22-4/3/22: complete outline of code for 75% deliverable and finish brainstorming

4/4/22-4/10/22: code and debug the 75% deliverable

4/11/22 **(Checkpoint)**: Be at 75% deliverable

4/12/22-4/17/22: start brainstorming on 100% deliverable and have code outline created

4/18/22-4/24/22: code and debug 100% deliverable

4/25/22-4/28/22: Be at 100% deliverable and complete writeup

4/29/22: Project due

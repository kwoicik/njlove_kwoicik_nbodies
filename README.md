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

We’re planning to use OpenMP to parallelize the program.

## Algorithm

For every iteration, create a parallel for loop indexing over number of threads. Within each thread:
1. Decide on a partition of the tree to claim, mark it as claimed somehow.
2. Calculate the forces on each particle in its partition.
3. Using the force on this particle and the "time" the iteration takes, calculate a new position.
4. Add the updated particle position to a new tree for next iteration.

So for every iteration, we construct a new tree out of the updated particle positions, and we are only READING the old positions.

## Goals and Deliverables
75%: We will implement a working, parallel solution for our N-bodies problem using a Barnes-Hut type algorithm

100%: In addition to our solution working it will also have a complex algorithm for semi-static workload assignment

125%: On top of the algorithm for balancing the workload our solution will also optimize for locality using a method similar to the “cost zones” we discussed in lecture

## Milestone

Thus far, we have a fully operational input file parser which takes in a text file containing information about bodies and uses Barnes-Hut to create a tree that can then be further parsed and modified by our algorithm. 

Overall, we haven't quite hit our 75% deliverable, but most of the work left towards it is debugging and we don't believe that will set us back too much. In short, we still think we will be able to be at our 100% goal by the end of the semester.

During the poster session, we plan to have a live demo of our visualizer displaying our data. Provided that the code runs fast enough, we will compute it live too, but if all else fails, we will just visualize an aleady computed nbodies run throughout all of its iterations until we hit convergence. This visualization will be plotted in matplotlib via python and will display where all of the bodies are within our bounding box. 

We don't have many concerns at this time, if anything, we expect the hard part should come when implementing semi-static workload assignment or "cost zones" if we get that far. Overall, we're feeling good!

## Schedule
3/21/22-3/23/22: create project proposal

3/28/22-4/3/22: complete outline of code for 75% deliverable and finish brainstorming

4/4/22-4/10/22: code and debug the 75% deliverable

4/11/22 **(Checkpoint)**: Be at 75% deliverable

4/12/22-4/17/22: (kelly) finish writing the output file writing code as well as the visualizer, (nathan) complete all the physics code and finish piecing together the parallel portion of the code

4/18/22-4/19/22: (kelly and nathan) brainstorm how to implement semi-static assignment on top of the pre-existing parallel code

4/20/22-4/25/22: (kelly and nathan) finish implementing the semi-static assignment, (kelly) start compiling datasets to display for the demo

4/26/22-4/28/22: (kelly and nathan) finish writing the deliverable for the project, (kelly) get the demo visualizations saved

4/29/22: Project due

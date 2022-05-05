
/*
 * Header File
 *
 * The goal of our project is to use a Barnes-Hut Algorithm to perform
 * atomic simulations in parallel. We'll do it in 2 dimensions to keep
 * life slightly simpler.
 *
 * Contains definitions for data structures, library imports, etc.
 *
 * Nathan Love (njlove) and Kelly Woicik (kwoicik)
 */

#ifndef __BARNESHUT_H__
#define __BARNESHUT_H__

// The "bodies" we'll be simulating will be atoms, rather than the more
// typical planets and stars
typedef struct
{
    // Position, self-explanatory
    double x;
    double y;

    // Velocity, hopefully also self-explanatory
    // Note that velocity is only meaningful for individual bodies, there's
    // no point in tracking it for conglomerates because they won't be moved
    // as a conglomerate. So it's **undefined** for conglomerates.
    double vx;
    double vy;

    // Where "w" stands for weight, as it contributes to the center of
    // mass, charge, whatever
    double w;

    // Where "c" stands for cost concerning the number of nodes this one has
    // interacted with in the previous iteration of the simulation
    int c;

    // Each individual body will need to have a unique ID so we can be sure
    // we don't calculate its interaction with itself. Conglomerates will
    // have an id of -1
    short id;

    // Indicates whether this part of the tree has been "claimed" by a thread
    short claimed;
} body_t;

// Each node will have up to 4 children. If one of them is NULL then
// that indicates that quadrant is empty. When all are empty we are at a
// leaf node, so total is an individual particle.
typedef struct node
{
    body_t total;
    int left;
    int right;
    int top;
    int bot;
    struct node* nw;
    struct node* ne;
    struct node* se;
    struct node* sw;
} node_t;

// Partitions will just be singly-linked lists
typedef struct link
{
    node_t* subtree;
    struct link* next;
} link_t;

typedef struct partition
{
    link_t* head;
} partition_t;

#endif

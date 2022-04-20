
/*
 * Implementation File
 *
 * Contains most of the logic for the actual computation, but not the
 * I/O or how to parse command-line options. We're assuming all of those
 * will be provided to us from main.
 */

#include "barnesHut.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

// Allocate and return a new, empty Barnes-Hut tree node
node_t* newTree(double left, double right, double bot, double top)
{
    // Allocate memory
    node_t* result = (node_t*) malloc(sizeof(node_t));

    // Init center of property
    result->total.x = 0;
    result->total.y = 0;
    result->total.w = 0;
    result->total.id = -1;

    // Init bounds
    result->left = left;
    result->right = right;
    result->bot = bot;
    result->top = top;

    // Init children
    result->nw = NULL;
    result->ne = NULL;
    result->sw = NULL;
    result->se = NULL;

    return result;
}

// Check whether a tree node is a leaf or interior node
bool isLeaf(node_t* tree)
{
    return (tree->nw == NULL && tree->ne == NULL &&
            tree->sw == NULL && tree->se == NULL);
}

// Check whether a node is an empty leaf
bool isEmpty(node_t* tree)
{
    return (isLeaf(tree) && tree->total.id == -1);
}

// Returns the quadrant in the tree that this coordinate falls into
node_t* quadrantFor(node_t* tree, double x, double y)
{
    double midx = tree->left + ((tree->right - tree->left) / 2);
    double  midy = tree-> bot + ((tree->top - tree->bot) / 2);
    // Western half
    if (x <= midx)
    {
        if (y <= midy)
        {
            return tree->sw;
        }

        return tree->nw;
    }

    // Eastern half
    else
    {
        if (y <= midy)
        {
            return tree->se;
        }

        return tree->ne;
    }
}

// Splits a leaf node into 4 new leaf nodes, automatically taking care
// of any particle it may have contained
void splitTree(node_t* tree)
{
    // Here for debugging, can be removed later
    if (!isLeaf(tree))
    {
        printf("Error, splitting a node that isn't a leaf\n");
        abort();
    }

    // Do the actual splitting
    double midx = tree->left + ((tree->right - tree->left) / 2);
    double midy = tree-> bot + ((tree->top - tree->bot) / 2);
    tree->nw = newTree(tree->left, midx, midy + 1, tree->top);
    tree->ne = newTree(midx + 1, tree->right, midy + 1, tree->top);
    tree->sw = newTree(tree->left, midx, tree->bot, midy);
    tree->se = newTree(midx + 1, tree-> right, tree->bot, midy);

    // Handle the moving of a particle if we had one
    // We always should, but I'll check just in case
    if (tree->total.id != -1)
    {
        node_t* newLeaf = quadrantFor(tree, tree->total.x, tree->total.y);
        newLeaf->total.x = tree->total.x;
        newLeaf->total.y = tree->total.y;
        newLeaf->total.w = tree->total.w;
        newLeaf->total.id = tree->total.id;
        tree->total.id = -1;
    }
}

// Insert a body into a tree
void insertBody(node_t* tree, double x, double y, double vx,
        double vy, double w, int id)
{
    // Base case; empty leaf node, insert here
    if (isEmpty(tree))
    {
        tree->total.x = x;
        tree->total.y = y;
        tree->total.vx = vx;
        tree->total.vy = vy;
        tree->total.w = w;
        tree->total.id = id;
        return;
    }

    // Non-empty leaf node, split it
    else if (isLeaf(tree))
    {
        splitTree(tree);
    }

    // This should now be an interior node if it wasn't already; recurse
    node_t* next = quadrantFor(tree, x, y);
    insertBody(next, x, y, vx, vy, w, id);

    // Update the totals for this node
    tree->total.x = ((w * x) + (tree->total.x * tree->total.w)) / (w + tree->total.w);
    tree->total.y = ((w * y) + (tree->total.y * tree->total.w)) / (w + tree->total.w);
    tree->total.w += w;
}

// Searches for and returns the node containing the body with the given id
node_t* findBody(node_t* tree, int id)
{
    // Base case, return this tree if the id matches otherwise NULL
    if (isLeaf(tree))
    {
        return (tree->total.id == id) ? tree : NULL;
    }

    node_t* res;
    res = findBody(tree->nw, id);
    if (res != NULL)
    {
        return res;
    }

    res = findBody(tree->ne, id);
    if (res != NULL)
    {
        return res;
    }

    res = findBody(tree->se, id);
    {
        return res;
    }

    return findBody(tree->sw, id);
}

// Frees a tree
void freeTree(node_t* tree)
{
    //TODO: free the tree
}

// Returns the distance between two points
double distanceBetween(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0));
}

// Wrapper that does distance between nodes, to save a little typing
double nodeDistance(node_t* t1, node_t* t2)
{
    return distanceBetween(t1->total.x, t1->total.y, t2->total.x, t2->total.y);
}

// Calculate the force between two bodies based on electric interactions,
// based on the Morse Potential
//
// I know it isn't super readable, it's just a straight-up physics equation
double forceBetween(double r, double r0, double u, double alpha)
{
    return -2 * u * alpha * (exp(-2 * alpha * (r - r0)) - exp(-1 * alpha * (r - r0)));
}

// Calculate the net force acting on a body at the given position from every
// body in the given tree
void netForce(node_t* tree, int id, double x, double y, double r0,
        double u, double alpha, double cutoff, double* netx, double* nety)
{
    if (isEmpty(tree))
    {
        return;
    }

    double d = distanceBetween(x, y, tree->total.x, tree->total.y);
    double ratio = (tree->right - tree->left) / d;

    // If we're under the ratio we can approximate
    // But also if we're at a leaf we have to stop here anyways
    if (ratio < cutoff || isLeaf(tree))
    {
        if (tree->total.id != id)
        {
            double f = forceBetween(d, r0, u, alpha);
            *netx += f * ((tree->total.x - x) / d);
            *nety += f * ((tree->total.y - y) / d);
        }
        return;
    }

    // Otherwise recurse!
    netForce(tree->nw, id, x, y, r0, u, alpha, cutoff, netx, nety);
    netForce(tree->ne, id, x, y, r0, u, alpha, cutoff, netx, nety);
    netForce(tree->se, id, x, y, r0, u, alpha, cutoff, netx, nety);
    netForce(tree->sw, id, x, y, r0, u, alpha, cutoff, netx, nety);
}

// Parses the input file and creates the tree
// TODO: should probably parallelize this in the final version
void parse_input_file(FILE* input, node_t** res, double* dimx, double* dimy, int* nBodies)
{
    fscanf(input, "%lf %lf\n", dimx, dimy);
    fscanf(input, "%d\n", nBodies);

    // initialize tree
    double left = - (*dimx / 2);
    double right = *dimx / 2;
    double bot = - (*dimy / 2);
    double top = *dimy / 2;
    *res = newTree(left, right, bot, top);

    // place all bodies into tree
    for (int i = 0; i < *nBodies; i++)
    {
        double x, y, vx, vy, w;
        fscanf(input, "%lf %lf %lf %lf %lf\n", &x, &y, &vx, &vy, &w);
        insertBody(*res, x, y, vx, vy, w, i);
    }
}

// Writes the positions of particles to the given output file
void write_output(FILE* output, node_t* tree, int nBodies)
{
    for (int i = 0; i < nBodies; i++)
    {
        node_t* cur = findBody(tree, i);
        fprintf(output, "%f %f\n", cur->total.x, cur->total.y);
    }
}

/*
 * Helpers for reading input arguments, shamelessly stolen right out of
 * Assignments 3/4
 */

static int _argc;
static const char **_argv;

const char *get_option_string(const char *option_name, const char *default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return _argv[i + 1];
    return default_value;
}

int get_option_int(const char *option_name, int default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return atoi(_argv[i + 1]);
    return default_value;
}

float get_option_float(const char *option_name, float default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return (float)atof(_argv[i + 1]);
    return default_value;
}

static void show_help(const char *program_path) {
    printf("Usage: %s OPTIONS\n", program_path);
    printf("\n");
    printf("OPTIONS:\n");
    printf("\t-i <input_filename> (required)\n");
    printf("\t-o <output_filename> (default out.txt)\n");
    printf("\t-t <num_of_threads> (default 1)\n");
    printf("\t-n <num_of_iterations> (default 10)\n");
    printf("\t-h <cutoff_ratio> (default 0.0)\n");
    printf("\t-u <morse_constant> (default 1.0)\n");
    printf("\t-r <equilibrium_radius> (default 1.0)\n");
    printf("\t-a <alpha> (default 1.0)\n");
    printf("\t-s <timestep> (default 1.0)\n");
}

/*
 * Main function--
 *
 * 1. Parse inputs
 * 2. Run simulation
 * 3. Write output file
 */
int main(int argc, const char *argv[])
{
    _argc = argc - 1;
    _argv = argv + 1;

    // Program parameters
    const char* inputFilename = get_option_string("-i", NULL);
    const char* outputFilename = get_option_string("-o", "out.txt");
    int numThreads = get_option_int("-t", 1);
    int numIters = get_option_int("-n", 10);
    double cutoffRatio = get_option_float("-h", 0.0);

    // Physical constants for simulation
    double morseConst = get_option_float("-u", 1.0);
    double eqbmRadius = get_option_float("-r", 1.0);
    double alpha = get_option_float("-a", 1.0);
    double stepSize = get_option_float("-s", 1.0);

    // Open input file
    FILE* input = fopen(inputFilename, "r");
    if (!input)
    {
        printf("Error: Unable to read file: %s, exiting...\n", inputFilename);
        return 1;
    }

    // Setup initial state from input file
    node_t* tree;
    double dimx, dimy;
    int nBodies;
    parse_input_file(input, &tree, &dimx, &dimy, &nBodies);
    double xmin = - (dimx / 2);
    double xmax = dimx / 2;
    double ymin = - (dimy / 2);
    double ymax = dimy / 2;

    // Open output file
    // Only once, so we can write multiple times without erasing
    FILE* output = fopen(outputFilename, "w");
    if (!output)
    {
        printf("Error: Unable to write output file: %s, exiting...\n", outputFilename);
        return 1;
    }
    fprintf(output, "%f %f\n", dimx, dimy);
    fprintf(output, "%d\n", nBodies);

    //TODO: need some stuff to setup and track timing

    // Simulation loop-- executes numIters times
    // 1. Create new, empty tree
    // 2. Fill in with the new positions of bodies
    // 3. Free old tree and set to new one
    for (int iter = 0; iter < numIters; iter++)
    {
        node_t* nextTree = newTree(xmin, xmax, ymin, ymax);
        for (int body = 0; body < nBodies; body++)
        {
            // Find current position and stuff
            node_t* cur = findBody(tree, body);
            double xf, yf, dvx, dvy, dx, dy;
            xf = 0;
            yf = 0;

            // Calculate force based on all other particles
            netForce(tree, body, cur->total.x, cur->total.y, eqbmRadius,
                    morseConst, alpha, cutoffRatio, &xf, &yf);

            // Update velocity assuming constant force during timestep
            dvx = xf * stepSize;
            dvy = yf * stepSize;

            // Update position based on average of old and new velocity
            dx = (cur->total.vx + (dvx / 2)) * stepSize;
            dy = (cur->total.vy + (dvy / 2)) * stepSize;

            // We know the changes, insert into new tree
            insertBody(nextTree, cur->total.x + dx, cur->total.y + dy,
                    cur->total.vx + dvx, cur->total.vy + dvy, cur->total.w, body);
        }

        freeTree(tree);
        tree = nextTree;

        write_output(output, tree, nBodies);
    }

    // Write outputs and do reporting
    // TODO: add that stuff
    return 0;
}


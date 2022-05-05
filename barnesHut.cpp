
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
#include <chrono>

// Allocate and return a new, empty Barnes-Hut tree node
node_t* newTree(double left, double right, double bot, double top)
{
    // Allocate memory
    node_t* result = (node_t*) malloc(sizeof(node_t));

    // Init center of property
    result->total.x = 0;
    result->total.y = 0;
    result->total.c = 0;
    result->total.w = 0;
    result->total.id = -1;
    result->total.claimed = -1;

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
    double midx = (tree->right + tree->left) / 2;
    double midy = (tree->top + tree->bot) / 2;
    tree->nw = newTree(tree->left, midx, midy, tree->top);
    tree->ne = newTree(midx, tree->right, midy, tree->top);
    tree->sw = newTree(tree->left, midx, tree->bot, midy);
    tree->se = newTree(midx, tree-> right, tree->bot, midy);

    // Handle the moving of a particle if we had one
    // We always should, but I'll check just in case
    if (tree->total.id != -1)
    {
        node_t* newLeaf = quadrantFor(tree, tree->total.x, tree->total.y);
        newLeaf->total.x = tree->total.x;
        newLeaf->total.y = tree->total.y;
        newLeaf->total.vx = tree->total.vx;
        newLeaf->total.vy = tree->total.vy;
        newLeaf->total.w = tree->total.w;
        newLeaf->total.id = tree->total.id;
        tree->total.id = -1;
    }
}

// Insert a body into a tree
void insertBody(node_t* tree, double x, double y, double vx,
        double vy, double w, int c, short id)
{
    /*
    printf("Inserting into tree: %lf %lf %lf %lf\n",
            tree->left, tree->right, tree->bot, tree->top);
            */
    // Base case; empty leaf node, insert here
    if (isEmpty(tree))
    {
        tree->total.x = x;
        tree->total.y = y;
        tree->total.vx = vx;
        tree->total.vy = vy;
        tree->total.c = c;
        tree->total.w = w;
        tree->total.id = id;
        return;
    }

    // Non-empty leaf node, split it
    else if (isLeaf(tree))
    {
        if (tree->total.id == id)
        {
            printf("Error: double-insertion\n");
            abort();
        }
        splitTree(tree);
    }

    // This should now be an interior node if it wasn't already; recurse
    node_t* next = quadrantFor(tree, x, y);
    insertBody(next, x, y, vx, vy, w, c, id);

    // Update the totals for this node
    tree->total.x = ((w * x) + (tree->total.x * tree->total.w)) / (w + tree->total.w);
    tree->total.y = ((w * y) + (tree->total.y * tree->total.w)) / (w + tree->total.w);
    tree->total.w += w;
    tree->total.c += c;
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
    if (res != NULL)
    {
        return res;
    }

    return findBody(tree->sw, id);
}

// Frees a tree
void freeTree(node_t* tree)
{
    if (isLeaf(tree))
    {
        free(tree);
        return;
    }

    freeTree(tree->nw);
    freeTree(tree->ne);
    freeTree(tree->se);
    freeTree(tree->sw);
    free(tree);
    return;
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
short netForce(node_t* tree, int id, double x, double y, double r0, double u,
        double alpha, double cutoff, double* netx, double* nety)
{
    if (isEmpty(tree))
    {
        return 0;
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
        return 1;
    }

    // Otherwise recurse!
    short res = 1;
    res += netForce(tree->nw, id, x, y, r0, u, alpha, cutoff, netx, nety);
    res += netForce(tree->ne, id, x, y, r0, u, alpha, cutoff, netx, nety);
    res += netForce(tree->se, id, x, y, r0, u, alpha, cutoff, netx, nety);
    res += netForce(tree->sw, id, x, y, r0, u, alpha, cutoff, netx, nety);
    return res;
}

// Parses the input file and creates the tree
void parse_input_file(const char* inputFilename, node_t** res, double* dimx,
        double* dimy, int* nBodies)
{
    FILE* input = fopen(inputFilename, "r");
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
        insertBody(*res, x, y, vx, vy, w, 1, i);
    }

    fclose(input);
}

// Parses the input file and copies the bodies into some arrays
void parseCopyInput(const char* inputFilename, double* xcopy, double* ycopy,
        double* vxcopy, double* vycopy, double* wcopy)
{
    double dimx, dimy;
    int nBodies;
    FILE* input = fopen(inputFilename, "r");
    fscanf(input, "%lf %lf\n", &dimx, &dimy);
    fscanf(input, "%d\n", &nBodies);

    // place all bodies into tree
    for (int i = 0; i < nBodies; i++)
    {
        fscanf(input, "%lf %lf %lf %lf %lf\n", &xcopy[i], &ycopy[i],
                &vxcopy[i], &vycopy[i], &wcopy[i]);
    }

    fclose(input);
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
    printf("\t-v <program_version> (1 = sequential, 2 = naive parallel, 3 = cost zones parallel)\n");
    printf("\t-l <num_laps> (default 10)\n");
    printf("\t-t <num_of_threads> (default 1)\n");
    printf("\t-n <num_of_iterations> (default 10)\n");
    printf("\t-h <cutoff_ratio> (default 0.0)\n");
    printf("\t-u <morse_constant> (default 1.0)\n");
    printf("\t-r <equilibrium_radius> (default 1.0)\n");
    printf("\t-a <alpha> (default 1.0)\n");
    printf("\t-s <timestep> (default 0.001)\n");
}

// Runs the simulation sequentially, with an added array to
// avoid searching for specific nodes
void runSimSequential(const char* inputFilename, const char* outputFilename, int numIters, double cutoffRatio,
        double morseConst, double eqbmRadius, double alpha, double stepSize)
{
    // Setup initial state from input file
    node_t* tree;
    double dimx, dimy;
    int nBodies;
    parse_input_file(inputFilename, &tree, &dimx, &dimy, &nBodies);

    double xmin = - (dimx / 2);
    double xmax = dimx / 2;
    double ymin = - (dimy / 2);
    double ymax = dimy / 2;

    double* xs = (double*) calloc(sizeof(double), nBodies);
    double* ys = (double*) calloc(sizeof(double), nBodies);
    double* vxs = (double*) calloc(sizeof(double), nBodies);
    double* vys = (double*) calloc(sizeof(double), nBodies);
    double* ws = (double*) calloc(sizeof(double), nBodies);
    parseCopyInput(inputFilename, xs, ys, vxs, vys, ws);

    // Open output file
    FILE* output = fopen(outputFilename, "w");

    // Add sim info to top of output file
    fprintf(output, "%lf %lf\n", dimx, dimy);
    fprintf(output, "%d\n", nBodies);

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
            double curx = xs[body];
            double cury = ys[body];
            double curvx = vxs[body];
            double curvy = vys[body];
            double curw = ws[body];

            double xf, yf, dvx, dvy, dx, dy;
            xf = 0;
            yf = 0;

            // Calculate force based on all other particles
            short c = netForce(tree, body, curx, cury, eqbmRadius,
                    morseConst, alpha, cutoffRatio, &xf, &yf);

            // Update velocity assuming constant force during timestep
            dvx = xf * stepSize;
            dvy = yf * stepSize;

            // Update position based on average of old and new velocity
            dx = (curvx + (dvx / 2)) * stepSize;
            dy = (curvy + (dvy / 2)) * stepSize;

            double newx = curx + dx;
            newx = (newx > xmax) ? newx - dimx : newx;
            newx = (newx < xmin) ? newx + dimx : newx;
            double newy = cury + dy;
            newy = (newy > ymax) ? newy - dimy : newy;
            newy = (newy < ymin) ? newy + dimy : newy;

            // We know the changes, insert into new tree
            insertBody(nextTree, newx, newy, curvx + dvx,
                    curvy + dvy, curw, c, body);

            // Update our arrays
            xs[body] = newx;
            ys[body] = newy;
            vxs[body] = curvx + dvx;
            vys[body] = curvy + dvy;
            ws[body] = curw;
        }

        freeTree(tree);
        tree = nextTree;

        write_output(output, tree, nBodies);
    }

    freeTree(tree);
    free(xs);
    free(ys);
    free(vxs);
    free(vys);
    fclose(output);
}

// Runs the simulation in parallel where each thread computes the
// same number of bodies each iteration. This version also uses the
// extra array
void runSimParallel(const char* inputFilename, const char* outputFilename, int numThreads, int numIters,
        double cutoffRatio, double morseConst, double eqbmRadius, double alpha, double stepSize)
{
    // Stuff to time the threads
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;
    double threadTimes[numThreads];


    // Setup initial state from input file
    node_t* tree;
    double dimx, dimy;
    int nBodies;
    parse_input_file(inputFilename, &tree, &dimx, &dimy, &nBodies);
    double xmin = - (dimx / 2);
    double xmax = dimx / 2;
    double ymin = - (dimy / 2);
    double ymax = dimy / 2;

    double* xs = (double*) calloc(sizeof(double), nBodies);
    double* ys = (double*) calloc(sizeof(double), nBodies);
    double* vxs = (double*) calloc(sizeof(double), nBodies);
    double* vys = (double*) calloc(sizeof(double), nBodies);
    double* ws = (double*) calloc(sizeof(double), nBodies);
    parseCopyInput(inputFilename, xs, ys, vxs, vys, ws);

    // Open output file
    FILE* output = fopen(outputFilename, "w");

    // Add sim info to top of output file
    fprintf(output, "%lf %lf\n", dimx, dimy);
    fprintf(output, "%d\n", nBodies);

    // Simulation loop-- executes numIters times
    // 1. Create new, empty tree
    // 2. Fill in with the new positions of bodies
    // 3. Free old tree and set to new one
    for (int iter = 0; iter < numIters; iter++)
    {
        node_t* nextTree = newTree(xmin, xmax, ymin, ymax);

        int thread;
        #pragma omp parallel for default(shared) private(thread) schedule(dynamic)
        for (thread = 0; thread < numThreads; thread++)
        {
            if (iter == 0)
            {
                threadTimes[thread] = 0;
            }
            auto start = Clock::now();

            for (int body = thread; body < nBodies; body += numThreads)
            {
                // Find current position and stuff
                double curx = xs[body];
                double cury = ys[body];
                double curvx = vxs[body];
                double curvy = vys[body];
                double curw = ws[body];

                double xf, yf, dvx, dvy, dx, dy;
                xf = 0;
                yf = 0;

                // Calculate force based on all other particles
                short c = netForce(tree, body, curx, cury, eqbmRadius,
                        morseConst, alpha, cutoffRatio, &xf, &yf);

                // Update velocity assuming constant force during timestep
                dvx = xf * stepSize;
                dvy = yf * stepSize;

                // Update position based on average of old and new velocity
                dx = (curvx + (dvx / 2)) * stepSize;
                dy = (curvy + (dvy / 2)) * stepSize;

                double newx = curx + dx;
                newx = (newx > xmax) ? newx - dimx : newx;
                newx = (newx < xmin) ? newx + dimx : newx;
                double newy = cury + dy;
                newy = (newy > ymax) ? newy - dimy : newy;
                newy = (newy < ymin) ? newy + dimy : newy;

                // Update our arrays
                xs[body] = newx;
                ys[body] = newy;
                vxs[body] = curvx + dvx;
                vys[body] = curvy + dvy;
                ws[body] = curw;

                // We know the changes, insert into new tree
                #pragma omp critical
                {
                    insertBody(nextTree, newx, newy, curvx + dvx,
                            curvy + dvy, curw, c, body);
                }
            }

            threadTimes[thread] += duration_cast<dsec>(Clock::now() - start).count();
        }

        freeTree(tree);
        tree = nextTree;

        write_output(output, tree, nBodies);
    }

    // Report workload balancing information
    for (int i = 0; i < numThreads; i++)
    {
        printf("Elapsed time for thread %d: %lf\n", i, threadTimes[i]);
    }

    freeTree(tree);
    free(xs);
    free(ys);
    free(vxs);
    free(vys);
    fclose(output);
}

// Extra helpers for partitioning trees by cost, needed for the last version

// Starts a new partition
partition_t* newPartition()
{
    partition_t* res = (partition_t*) malloc(sizeof(partition_t));
    res->head = NULL;
    return res;
}

// Inserts another subtree into a partition
void partitionInsert(partition_t* p, node_t* t)
{
    link_t* newLink = (link_t*) malloc(sizeof(link_t));
    newLink->subtree = t;
    newLink->next = p->head;
    p->head = newLink;
}

// Frees a list
void listFree(link_t* l)
{
    if (l == NULL)
    {
        return;
    }

    listFree(l->next);
    free(l);
}

// Frees a partition
void partitionFree(partition_t* p)
{
    listFree(p->head);
    free(p);
}

// Claims a partition of a tree and returns its cost
int claimPartition(partition_t* part, node_t* tree, int limit, short thread)
{
    // Base case 1: we hit limit, tree is NULL or we can't claim anymore
    if (limit == 0 || tree == NULL || tree->total.claimed != -1)
    {
        return 0;
    }

    // Base case 2: tree is small enough to claim
    if (tree->total.c < limit)
    {
        partitionInsert(part, tree);
        tree->total.claimed = thread;
        return tree->total.c;
    }

    // Otherwise, try claiming from each of its children
    int claim = claimPartition(part, tree->nw, limit, thread);
    claim += claimPartition(part, tree->ne, limit - claim, thread);
    claim += claimPartition(part, tree->se, limit - claim, thread);
    claim += claimPartition(part, tree->sw, limit - claim, thread);
    return claim;
}

// Compute changes for all bodies in subtree based on oldTree, then add them
// to newTree
void computePartition(node_t* subtree, node_t* oldTree, node_t* newTree, short thread, double eqbmRadius,
        double morseConst, double alpha, double cutoffRatio, double stepSize, double dimx, double dimy)
{
    // Base case 1: empty subtree or claimed by another, do nothing
    if (isEmpty(subtree) ||
            (subtree->total.claimed != -1 && subtree->total.claimed != thread))
    {
        return;
    }

    // Base case 2: body, do computation
    if (isLeaf(subtree))
    {

        // Find current position and stuff
        double xf, yf, dvx, dvy, dx, dy;
        xf = 0;
        yf = 0;

        // Calculate force based on all other particles
        short c = netForce(oldTree, subtree->total.id, subtree->total.x, subtree->total.y,
                eqbmRadius, morseConst, alpha, cutoffRatio, &xf, &yf);

        // Update velocity assuming constant force during timestep
        dvx = xf * stepSize;
        dvy = yf * stepSize;

        // Update position based on average of old and new velocity
        dx = (subtree->total.vx + (dvx / 2)) * stepSize;
        dy = (subtree->total.vy + (dvy / 2)) * stepSize;

        // Recompute the bounds
        double xmin = - (dimx / 2);
        double xmax = dimx / 2;
        double ymin = - (dimy / 2);
        double ymax = dimy / 2;

        double newx = subtree->total.x + dx;
        newx = (newx > xmax) ? newx - dimx : newx;
        newx = (newx < xmin) ? newx + dimx : newx;
        double newy = subtree->total.y + dy;
        newy = (newy > ymax) ? newy - dimy : newy;
        newy = (newy < ymin) ? newy + dimy : newy;

        // We know the changes, insert into new tree
        #pragma omp critical ( insertion )
        {
            insertBody(newTree, newx, newy, subtree->total.vx + dvx,
                    subtree->total.vy + dvy, subtree->total.w, c, subtree->total.id);
        }
        return;
    }

    // Recursive case
    computePartition(subtree->nw, oldTree, newTree, thread, eqbmRadius,
            morseConst, alpha, cutoffRatio, stepSize, dimx, dimy);
    computePartition(subtree->ne, oldTree, newTree, thread, eqbmRadius,
            morseConst, alpha, cutoffRatio, stepSize, dimx, dimy);
    computePartition(subtree->se, oldTree, newTree, thread, eqbmRadius,
            morseConst, alpha, cutoffRatio, stepSize, dimx, dimy);
    computePartition(subtree->sw, oldTree, newTree, thread, eqbmRadius,
            morseConst, alpha, cutoffRatio, stepSize, dimx, dimy);
}

// Call computePartition for all of the subtrees in the given partition
void computeAll(partition_t* p, node_t* oldTree, node_t* newTree, short thread, double eqbmRadius,
        double morseConst, double alpha, double cutoffRatio, double stepSize, double dimx, double dimy)
{
    if (p != NULL)
    {
        for (link_t* curPar = p->head; curPar != NULL; curPar = curPar->next)
        {
            computePartition(curPar->subtree, oldTree, newTree, thread, eqbmRadius,
                    morseConst, alpha, cutoffRatio, stepSize, dimx, dimy);
        }
    }
}

// Runs the simulation in parallel where each thread computes a
// sub tree that is theoretically about equal in cost
//
// This is accomplished by keeping track of how many nodes each body
// interacted with on the previous iteration, and then having each
// thread pick a partition of the tree with approximately equal interactions
void runSimCosts(const char* inputFilename, const char* outputFilename, int numThreads, int numIters,
        double cutoffRatio, double morseConst, double eqbmRadius, double alpha, double stepSize)
{
    // Stuff to time the threads
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;
    double threadTimes[numThreads];

    // Setup initial state from input file
    node_t* tree;
    double dimx, dimy;
    int nBodies;
    parse_input_file(inputFilename, &tree, &dimx, &dimy, &nBodies);
    double xmin = - (dimx / 2);
    double xmax = dimx / 2;
    double ymin = - (dimy / 2);
    double ymax = dimy / 2;

    // Open input/output files
    FILE* output = fopen(outputFilename, "w");

    // Add sim info to top of output file
    fprintf(output, "%lf %lf\n", dimx, dimy);
    fprintf(output, "%d\n", nBodies);

    // Simulation loop-- executes numIters times
    // 1. Create new, empty tree
    // 2. Fill in with the new positions of bodies
    // 3. Free old tree and set to new one
    for (int iter = 0; iter < numIters; iter++)
    {
        node_t* nextTree = newTree(xmin, xmax, ymin, ymax);

        int thread;
        #pragma omp parallel for default(shared) private(thread) schedule(dynamic)
        for (thread = 0; thread < numThreads; thread++)
        {
            if (iter == 0)
            {
                threadTimes[thread] = 0;
            }
            auto start = Clock::now();

            // Claim partitions of (hopefully) equal cost
            partition_t* part = newPartition();
            int limit = (tree->total.c / numThreads) + 1;

            #pragma omp critical ( partition )
            {
                claimPartition(part, tree, limit, thread);
            }

            // Do the computation for our claimed partition
            computeAll(part, tree, nextTree, thread, eqbmRadius, morseConst,
                    alpha, cutoffRatio, stepSize, dimx, dimy);

            partitionFree(part);

            threadTimes[thread] += duration_cast<dsec>(Clock::now() - start).count();
        }

        // Clean up any stragglers
        partition_t* temp = newPartition();
        partitionInsert(temp, tree);
        computeAll(temp, tree, nextTree, -1, eqbmRadius, morseConst,
                alpha, cutoffRatio, stepSize, dimx, dimy);
        partitionFree(temp);

        freeTree(tree);
        tree = nextTree;

        write_output(output, tree, nBodies);
    }

    // Report workload balancing information
    for (int i = 0; i < numThreads; i++)
    {
        printf("Elapsed time for thread %d: %lf\n", i, threadTimes[i]);
    }

    freeTree(tree);
    fclose(output);
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
    // Stuff to time the program
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;
    double totalTime = 0;

    // Program parameters
    _argc = argc - 1;
    _argv = argv + 1;

    const char* inputFilename = get_option_string("-i", NULL);
    const char* outputFilename = get_option_string("-o", "out.txt");
    int numThreads = get_option_int("-t", 1);
    int numIters = get_option_int("-n", 1000);
    double cutoffRatio = get_option_float("-h", 0.0);
    int ver = get_option_int("-v", 1);
    int laps = get_option_int("-l", 10);

    // Test open input file
    FILE* input = fopen(inputFilename, "r");
    if (!input)
    {
        printf("Error: Unable to read input file: %s, exiting...\n", inputFilename);
        show_help("barnesHut");
        return 1;
    }
    fclose(input);

    // Test open output file
    FILE* output = fopen(outputFilename, "w");
    if (!output)
    {
        printf("Error: Unable to write output file: %s, exiting...\n", outputFilename);
        return 1;
    }
    fclose(output);

    // Physical constants for simulation
    double morseConst = get_option_float("-u", 1.0);
    double eqbmRadius = get_option_float("-r", 1.0);
    double alpha = get_option_float("-a", 1.0);
    double stepSize = get_option_float("-s", 0.001);

    for (int lap = 0; lap < laps + 1; lap++)
    {
        auto start = Clock::now();

        switch(ver)
        {
            case 1:
                runSimSequential(inputFilename, outputFilename, numIters,
                        cutoffRatio, morseConst, eqbmRadius, alpha, stepSize);
                break;

            case 2:
                runSimParallel(inputFilename, outputFilename, numThreads, numIters,
                        cutoffRatio, morseConst, eqbmRadius, alpha, stepSize);
                break;

            case 3:
                runSimCosts(inputFilename, outputFilename, numThreads, numIters,
                        cutoffRatio, morseConst, eqbmRadius, alpha, stepSize);
                break;
        }

        // Do an extra, untimed lap first to warm up caches
        if (lap != 0)
            {
                totalTime += duration_cast<dsec>(Clock::now() - start).count();
            }
    }

    // Report timing
    printf("Average time per run across %d runs: %lf\n", laps, totalTime / laps);
    return 0;
}


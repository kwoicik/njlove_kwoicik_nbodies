/**
 * Parses Input file and creates a Barnes-Hut implemented tree of the data.
 * Also writes all bodies to an output file to later be visualized.
 * Nathan Love (njlove), Kelly Woicik (kwoicik)
 */

#include "parse_input_file.h"

#include <cstdio>
#include <cstdlib>
#include <unistd.h>

// calculates the center of mass given a list of bodies
int calculate_cmass(body_t* bodies, int num_bodies, bool x)
{
    int sum = 0;
    for (int i = 0; i < num_bodies; i++)
    {
        if (x) {
            sum += bodies[i].x;
        }
        else {
            sum += bodies[i].y;
        }
    }
    return sum / num_bodies;
}

// checks whether a body is in the given bounding box
bool body_in_bounding_box(int x1, int y1, int x2, int y2, body_t body)
{
    return body.x >= x1 && body.y <= y1 && body.x <= x2 && body.y >= y2;
}

// finds which node the given body will be traversing next
int find_next_node_from_body(node_t *node, body_t body)
{
    bool in_zero = body_in_bounding_box(node->nodes[0]->bounds[0], 
                                        node->nodes[0]->bounds[1],
                                        node->nodes[0]->bounds[2],
                                        node->nodes[0]->bounds[3],
                                        body);
    bool in_one = body_in_bounding_box(node->nodes[1]->bounds[0], 
                                        node->nodes[1]->bounds[1],
                                        node->nodes[1]->bounds[2],
                                        node->nodes[1]->bounds[3],
                                        body);
    bool in_two = body_in_bounding_box(node->nodes[2]->bounds[0], 
                                        node->nodes[2]->bounds[1],
                                        node->nodes[2]->bounds[2],
                                        node->nodes[2]->bounds[3],
                                        body);
    bool in_three = body_in_bounding_box(node->nodes[3]->bounds[0], 
                                        node->nodes[3]->bounds[1],
                                        node->nodes[3]->bounds[2],
                                        node->nodes[3]->bounds[3],
                                        body);
    return (in_zero) ? 0 : (in_one) ? 1 : (in_two) ? 2 : 3;
}

// adds a body to a bodies list and updates the body count
void add_to_bodies_list(node_t *node, body_t body)
{
    body_t* new_bodies = (body_t *) malloc((node->num_bodies + 1) * sizeof(body_t));
    for (int i = 0; i < node->num_bodies; i++)
    {
        new_bodies[i] = node->bodies[i];
    }
    new_bodies[node->num_bodies] = body;
    free(node->bodies);
    node->bodies = new_bodies;
    node->num_bodies += 1;
}

// calculates the new bounding box of a leaf node from its parents
void update_leaf_bounding_box(node_t *node, int i)
{
    int old_x1, old_y1, old_x2, old_y2;
    int mid_x, mid_y;
    old_x1 = node->bounds[0];
    old_y1 = node->bounds[1];
    old_x2 = node->bounds[2];
    old_y2 = node->bounds[3];
    mid_x = old_x1 + (old_x2 - old_x1)/2; 
    mid_y = old_y1 - (old_y1 - old_y2)/2;

    int x1, y1, x2, y2;
    x1 = (i == 0 || i == 2) ? old_x1 : mid_x;
    y1 = (i == 0 || i == 1) ? old_y1 : mid_y;
    x2 = (i == 0 || i == 2) ? mid_x : old_x2;
    y2 = (i == 0 || i == 1) ? mid_y : old_y2;

    node->nodes[i]->bounds[0] = x1;
    node->nodes[i]->bounds[1] = y1;
    node->nodes[i]->bounds[2] = x2;
    node->nodes[i]->bounds[3] = y2;
    printf("\ni: %d\nx1, y1, x2, y2 -> %d, %d, %d, %d\n", i, x1, y1, x2, y2);
}

// places a body on the tree
void place_body(node_t *node, body_t body)
{
    printf("body %d, %d\n", body.x, body.y);
    if (node->is_leaf && node->num_bodies == 0)
    {
        printf("entered empty leaf node\n");
        // we now claim this leaf with our body
        node->num_bodies = 1;

        body_t* new_bodies = (body_t *) malloc(sizeof(body_t));
        new_bodies[0] = body;
        node->bodies = new_bodies;

        node->cmass_x = body.x;
        node->cmass_y = body.y;
        return;
    }
    else if (node->is_leaf && node->num_bodies != 0)
    {
        printf("entered non-empty leaf node\n");
        // we need to split up this current leaf into 4 children leaves
        node->is_leaf = false;

        // update bodies list
        add_to_bodies_list(node, body);

        // calculate new center of mass
        node->cmass_x = calculate_cmass(node->bodies, node->num_bodies, true);
        node->cmass_y = calculate_cmass(node->bodies, node->num_bodies, false);

        // create new nodes
        for (int i = 0; i < 4; i++) 
        {
            node->nodes[i] = (node_t *) malloc(sizeof(node_t));
            node->nodes[i]->is_leaf = true;
            node->nodes[i]->num_bodies = 0;
            node->nodes[i]->cmass_x = 0;
            node->nodes[i]->cmass_y = 0;
            
            update_leaf_bounding_box(node, i);
        }

        // place both bodies
        for (int i = 0; i < node->num_bodies; i++)
        {
            body_t cur_body = node->bodies[i];
            int node_to_traverse_next = find_next_node_from_body(node, cur_body);
            printf("node to enter: %d\n", node_to_traverse_next);
            place_body(node->nodes[node_to_traverse_next], cur_body);
        }
        return;
        
    }
    else {
        // update current interior node
        add_to_bodies_list(node, body);
        node->cmass_x = calculate_cmass(node->bodies, node->num_bodies, true);
        node->cmass_y = calculate_cmass(node->bodies, node->num_bodies, false);

        // we keep parsing until we get a leaf
        int node_to_traverse_next = find_next_node_from_body(node, body);
        place_body((node->nodes[node_to_traverse_next]), body);
        return;
    }
}

// parses the input file and creates the tree
int parse_input_file(const char *input_filename, node_t* root)  
{
    FILE *input = fopen(input_filename, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", input_filename);
        return 1;
    }

    int dim_x, dim_y;
    int num_bodies;

    fscanf(input, "%d %d\n", &dim_x, &dim_y);
    fscanf(input, "%d\n", &num_bodies);

    // initialize tree;
    root->num_bodies = 0;
    root->is_leaf = true;
    root->cmass_x = 0;
    root->cmass_y = 0;
    root->bounds[0] = 0;
    root->bounds[1] = dim_y;
    root->bounds[2] = dim_x;
    root->bounds[3] = 0;

    // place all bodies into tree
    for (int i = 0; i < num_bodies; i++)
    {
        body_t cur;
        fscanf(input, "%d %d\n", &cur.x, &cur.y);

        place_body(root, cur);
    }
    
    return 0;
}

void write_output_file(int dim_x, int dim_y, const char *output_filename, node_t *root)
{
    FILE *output = fopen(output_filename, "w");

    fprintf(output, "%d %d\n", dim_x, dim_y);
    fprintf(output, "%d\n", root->num_bodies);

    for (int i = 0; i < root->num_bodies; i++)
    {
        fprintf(output, "%d %d\n", root->bodies[i].x, root->bodies[i].y);
    }
}

int main(int argc, char *argv[])
{
    int opt = 0;
    const char *input_filename;

    opt = getopt(argc, argv, "f:");
    switch (opt) {
        case 'f':
            input_filename = optarg;
            break;
        default:
            printf("error");
    }

    node_t* root = (node_t*) malloc(sizeof(node_t));

    parse_input_file(input_filename, root);

    printf("num_bodies: %d\n", root->num_bodies);
    printf("bod1 %d, %d\n", root->bodies[0].x, root->bodies[0].y);
    printf("bod2 %d, %d\n", root->bodies[1].x, root->bodies[1].y);
    printf("cmass %d, %d\n", root->cmass_x, root->cmass_y);
    printf("node1 bods %d, %d\n", root->nodes[1]->bodies[0].x, root->nodes[1]->bodies[0].y);
    printf("node2 bods %d, %d\n", root->nodes[2]->bodies[0].x, root->nodes[2]->bodies[0].y);

    FILE *input = fopen(input_filename, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", input_filename);
        return 1;
    }

    int dim_x, dim_y;
    fscanf(input, "%d %d\n", &dim_x, &dim_y);

    write_output_file(dim_x, dim_y, "output_test1.txt", root);

    return 0;
}


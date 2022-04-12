/**
 * Parses Input file and creates a Barnes-Hut implemented tree of the data
 * Nathan Love (njlove), Kelly Woicik (kwoicik)
 */

 typedef struct {
    int x;
    int y;
 } body_t;

 typedef struct Node {
    body_t* bodies; // pointer to what will be an array of bodies (includes all children)
    int num_bodies; // includes all children
    bool is_leaf;
    int cmass_x;
    int cmass_y;
    struct Node* nodes[4]; // pointers to children
    int bounds[4]; // contains bounding box of tree node [x1 y1 x2 y2]
                   // bounding bodes denote top left corner and bottom right corner
                   // the origin is at the bottom left of the screen
 } node_t;
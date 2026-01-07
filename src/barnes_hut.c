#include "barnes_hut.h"

/** Barnes-Hut algorithm
 * Based on the clustered particles, uterlize the spacial and time locality to improve the cache performance.
 * - Build the local tree for each cluster, parallelize the instertion of the particles => no significant improvement.
 * - Or build the global tree
 *  -------------------------------------------------- */
void barnes_hut(double* pos_x, double* pos_y, double* mass, int N, int cluster[][N], double* region, int* clusters_size, int k, double* fx, double* fy, int n_threads, double THETA_MAX) {   
    /* initialize the forces */
    for (int i = 0; i < N; i++) {
        fx[i] = 0;
        fy[i] = 0;
    }

    /*
        // Each cluster has a local tree
        TNode* tTree[k];
        # pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < k; i++) {
            tTree[i] = clusters_size[i] == 0 ? NULL : create_new_TNode(-1, region[0], region[1], region[2], region[3]);
            tTree[i]->PID = cluster[i][0];
            tTree[i]->pos_x = pos_x[cluster[i][0]];
            tTree[i]->pos_y = pos_y[cluster[i][0]];
            for(int j = 1; j < clusters_size[i]; j++){
                insert(tTree[i], pos_x[cluster[i][j]], pos_y[cluster[i][j]], mass[cluster[i][j]], cluster[i][j]);
            }
        }

        // compute the forces
        for(int i = 0; i < k; i++){
            #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
            for(int j = 0; j < clusters_size[i]; j++){
                // each particle traverses all the trees
                int PID = cluster[i][j];
                for(int kk = 0; kk < k; kk++){
                    compute_force(pos_x[PID], pos_y[PID], mass[PID], PID, tTree[kk], &fx[PID], &fy[PID], N, THETA_MAX);
                }
            }
        }

        // free trees
        for (int  i = 0; i < k; i++){
            destroy(tTree[i]);
        }
    */
    

    /* Build the global tree */
    /* initialize the root */
    TNode* tTree = create_new_TNode(-1, region[0], region[1], region[2], region[3]);
    tTree->PID = 0;
    tTree->pos_x = pos_x[0];
    tTree->pos_y = pos_y[0];
    tTree->mass = mass[0];
    /* insert each particle one by one */
    for(int i=1; i < N; i++){
        insert(tTree, pos_x[i], pos_y[i], mass[i], i);
    }

    /* Compute the forces */
    for(int i = 0; i < k; i++){
        #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
        for(int j = 0; j < clusters_size[i]; j++){
            /* each particle traverses all the trees */
            int PID = cluster[i][j];
            compute_force(pos_x[PID], pos_y[PID], mass[PID], PID, tTree, &fx[PID], &fy[PID], N, THETA_MAX);
        }
    }

    /* Free the tree */
    destroy(tTree);
    
}

/** Create a new tree node
 * According the position of the subsquare, determine the region of it.
 * -------------------------------------------------- */
TNode* create_new_TNode(int index, double LB, double RB, double DB, double UB) {
    TNode* new_TNode = malloc(sizeof(TNode));
    new_TNode->child[0] = NULL;
    new_TNode->child[1] = NULL;
    new_TNode->child[2] = NULL;
    new_TNode->child[3] = NULL;

    double mid_x = 0.5 * (LB + RB);
    double mid_y = 0.5 * (DB + UB);
    switch (index) {
        case -1:
            new_TNode->LB = LB;
            new_TNode->RB = RB;
            new_TNode->DB = DB;
            new_TNode->UB = UB;
            break;
        case 0:
            new_TNode->LB = LB;
            new_TNode->RB = mid_x;
            new_TNode->DB = DB;
            new_TNode->UB = mid_y;
            break;
        case 1:
            new_TNode->LB = LB;
            new_TNode->RB = mid_x;
            new_TNode->DB = mid_y;
            new_TNode->UB = UB;
            break;
        case 2:
            new_TNode->LB = mid_x;
            new_TNode->RB = RB;
            new_TNode->DB = DB;
            new_TNode->UB = mid_y;
            break;
        case 3:
            new_TNode->LB = mid_x;
            new_TNode->RB = RB;
            new_TNode->DB = mid_y;
            new_TNode->UB = UB;
            break;
    }
    return new_TNode;
}

/** Insert the particle into the tree
 * If the node represents a particle, then split the node and insert the particle into the corresponding child node.
 * If the node is representing a big object, then update the massive and the centroid of the node.
 * - determine the region the particle belongs to and insert it into the corresponding child node.
 * - if the child node is NULL, then create a new node and insert the particle into it.
 * - else, recursively insert the particle into the child node.
 * -------------------------------------------------- */
int insert(TNode* tNode, double pos_x, double pos_y, double mass, int PID) {
    double mid_x = 0.5 * (tNode->LB + tNode->RB);
    double mid_y = 0.5 * (tNode->UB + tNode->DB);

    if (tNode->PID != -1) {
        if ((pos_x == tNode->pos_x) && (pos_y == tNode->pos_y)) {
            printf("Two particles are detected at the same location and the simulation terminates.\n");
            exit(0);
        }
        int index = (tNode->pos_y > mid_y) + 2 * (tNode->pos_x > mid_x);
        tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
        tNode->child[index]->PID = tNode->PID;
        tNode->child[index]->pos_x = tNode->pos_x;
        tNode->child[index]->pos_y = tNode->pos_y;
        tNode->child[index]->mass = tNode->mass;
        tNode->PID = -1;
    }

    /* Update the massive and the centroid */
    double new_mass = tNode->mass + mass;
    tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
    tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
    tNode->mass = new_mass;

    int index = (pos_y > mid_y) + 2 * (pos_x > mid_x);
    if (tNode->child[index] == NULL) {
        tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
        tNode->child[index]->pos_x = pos_x;
        tNode->child[index]->pos_y = pos_y;
        tNode->child[index]->mass = mass;
        tNode->child[index]->PID = PID;
    } else {
        insert(tNode->child[index], pos_x, pos_y, mass, PID);
    }
    return 0;
}



void _compute_force(double pos_x, double pos_y, double mass, TNode* tNode, double* fx, double* fy, double G) {
    double r_x = pos_x - tNode->pos_x;
    double r_y = pos_y - tNode->pos_y;
    double r_squared = r_x * r_x + r_y * r_y;
    double r_plummer = sqrt(r_squared) + EPSILON_O;
    double force_factor = -G * mass * tNode->mass / (r_plummer * r_plummer * r_plummer);
    *fx += force_factor * r_x;
    *fy += force_factor * r_y;
}


/** Compute the force between the particle and the node
 * If the node is a leaf node, directly calculate the force.
 * If the node is a big object, calculate the theta to check if the node is far enough, then calculate the force.
 * - calculate the distance between the particle and the node.
 * - calculate the force factor and update the force.
 * -------------------------------------------------- */
void compute_force(double pos_x, double pos_y, double mass, int PID, TNode* tNode, double* fx, double* fy, int N, double THETA_MAX) {
    double G = 100.0 / N;
    /* The node is empty or the node is the particle itself */
    if (!tNode || tNode->PID == PID) {
        return;
    }

    /* The node is a leaf node, directly calculate the force
       or calculate theta to check if the node is far enough, then calculate the force */
    if (tNode->PID != -1) {
        _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
        
    } else {
        double mid_x = 0.5 * (tNode->LB + tNode->RB);
        double mid_y = 0.5 * (tNode->DB + tNode->UB);
        double width = (tNode->RB - tNode->LB);
        double theta = width / sqrt((pos_x - mid_x) * (pos_x - mid_x) +
                                    (pos_y - mid_y) * (pos_y - mid_y));

        if (theta <= THETA_MAX) {
            _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
        } else {
            compute_force(pos_x, pos_y, mass, PID, tNode->child[0], fx, fy, N, THETA_MAX);
            compute_force(pos_x, pos_y, mass, PID, tNode->child[1], fx, fy, N, THETA_MAX);
            compute_force(pos_x, pos_y, mass, PID, tNode->child[2], fx, fy, N, THETA_MAX);
            compute_force(pos_x, pos_y, mass, PID, tNode->child[3], fx, fy, N, THETA_MAX);
        }
    }
}


/** Destroy the tree recursively
 *  -------------------------------------------------- */
void destroy(TNode* root) {
    if (root == NULL) {
        return;
    } else if (root->PID == -1) {
        destroy(root->child[0]);
        destroy(root->child[1]);
        destroy(root->child[2]);
        destroy(root->child[3]);
    }
    free(root);
}
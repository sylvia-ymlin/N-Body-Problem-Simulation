#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define THETA_MAX 0.5
#define TWO_ORDER 1

double G;
int N;

typedef struct PNode {
    double pos_x;
    double pos_y;
    double mass;
} PNode;

typedef struct TNode {
    double LB, RB, DB, UB;
    struct TNode* child[4];
    PNode* particle;
    bool is_leaf;
} TNode;

TNode* create_new_TNode(int index, double LB, double RB, double DB, double UB);
TNode* buildTree(PNode* particles, int N, float LB, float RB, float DB, float UB, int n_threads);
int insert(TNode* qNode, PNode* particle);                                                         
void barnesHut(PNode* particle, TNode* qNode, double* fx, double* fy);
void preOrder(TNode* tNode);
void destroy(TNode* root);

static double get_wall_seconds() {  // Returns the current time in seconds, when openmon is not used
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

int main(int argc, char* argv[]) {

    #ifdef _OPENMP
    double time_tol = omp_get_wtime();
    #else
    double time_tol = get_wall_seconds();
    #endif

    if (argc != 6) {
        printf("You should enter the following parameters in order:\n");
        printf("N filname nsteps delta_t n_threads\n");
        return 1;
    }

    int N = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int n_threads = atoi(argv[5]);

    // int N = 3;
    // char* filename = "./input_data/sun_and_planets_N_3.gal";
    // int nsteps = 1000000;
    // double delta_t = 1e-8;

    FILE* data_file = fopen(filename, "rb");
    if (data_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    PNode particles[N];
    double vx[N];
    double vy[N];
    double brightness[N];
    double fx[N];
    double fy[N];
    double mass_inver[N];
    double acc_x[N];
    double acc_y[N];

    // Read the data from the file
    for (int i = 0; i < N; i++) {
        fread(&particles[i], sizeof(PNode), 1, data_file);
        mass_inver[i] = 1.0 / particles[i].mass;
        fread(&vx[i], sizeof(double), 1, data_file);
        fread(&vy[i], sizeof(double), 1, data_file);
        fread(&brightness[i], sizeof(double), 1, data_file);
        acc_x[i] = 0.0;
        acc_y[i] = 0.0;
    }

    fclose(data_file);

    double LB = -1.0;
    double RB = 2.0;
    double DB = -1.0;
    double UB = 2.0;

    G = 100.0 / N;

    // initial acceleration
    #if TWO_ORDER
    for (int i = 0; i < N; i++) {
        fx[i] = 0;
        fy[i] = 0;
    }
    TNode* tTree = buildTree(particles, N, LB, RB, DB, UB, n_threads);
    for (int i = 0; i < N; i++) {
        barnesHut(&particles[i], tTree, &fx[i], &fy[i]);
    }
    destroy(tTree);
    #endif

    // Time integration
    for (int step = 0; step < nsteps; step++) {

        #if TWO_ORDER
        for (int i = 0; i < N; i++) {
            acc_x[i] = fx[i] * mass_inver[i];
            acc_y[i] = fy[i] * mass_inver[i];
            particles[i].pos_x += delta_t * vx[i] + 0.5 * delta_t * delta_t * acc_x[i];
            particles[i].pos_y += delta_t * vy[i] + 0.5 * delta_t * delta_t * acc_y[i];
            if(particles[i].pos_x < LB || particles[i].pos_x > RB || particles[i].pos_y < DB || particles[i].pos_y > UB) {
                printf("At least one particle is out of the region, and the simulation has been terminated.\n");
                exit(0);
            }
        }
        #endif

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
        #endif
        for (int i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
        }

        TNode* tTree = buildTree(particles, N, LB, RB, DB, UB, n_threads);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
        #endif
        for (int i = 0; i < N; i++) {
            barnesHut(&particles[i], tTree, &fx[i], &fy[i]);
        }

        // update
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
        #endif
        for (int i = 0; i < N; i++) {
            #if TWO_ORDER
            vx[i] += 0.5 * delta_t * (fx[i] * mass_inver[i] + acc_x[i]);
            vy[i] += 0.5 * delta_t * (fy[i] * mass_inver[i] + acc_y[i]);
            #else
            vx[i] += delta_t * fx[i] * mass_inver[i];
            vy[i] += delta_t * fy[i] * mass_inver[i];
            particles[i].pos_x += delta_t * vx[i];
            particles[i].pos_y += delta_t * vy[i];
            if (particles[i].pos_x < LB || particles[i].pos_x > RB || particles[i].pos_y < DB || particles[i].pos_y > UB) {
                printf("At least one particle is out of the region, and the simulation has been terminated.\n");
                exit(0);
            }
            #endif
        }

        destroy(tTree);
    }

    FILE* rfile = fopen("result.gal", "w");

    if (rfile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fwrite(&particles[i], sizeof(PNode), 1, rfile);
        fwrite(&vx[i], sizeof(double), 1, rfile);
        fwrite(&vy[i], sizeof(double), 1, rfile);
        fwrite(&brightness[i], sizeof(double), 1, rfile);
    }

    #ifdef _OPENMP
    printf("f_std tests took %7.8f wall seconds.\n", omp_get_wtime() - time_tol);
    #else
    printf("f_std tests took %7.8f wall seconds.\n", get_wall_seconds() - time_tol);
    #endif

    return 0;
}


TNode* buildTree(PNode* particles, int N, float LB, float RB, float DB, float UB, int n_threads) {
    // Divide the region into 4 squares
    float square1[4] = {LB, 0.5 * (LB + RB), DB, 0.5 * (DB + UB)};
    float square2[4] = {0.5 * (LB + RB), RB, DB, 0.5 * (DB + UB)};
    float square3[4] = {LB, 0.5 * (LB + RB), 0.5 * (DB + UB), UB};
    float square4[4] = {0.5 * (LB + RB), RB, 0.5 * (DB + UB), UB};
    // Initialize the groups of particles
    int group1_ofParticles[N];
    int N1 = 0;
    int group2_ofParticles[N];
    int N2 = 0;
    int group3_ofParticles[N];
    int N3 = 0;
    int group4_ofParticles[N];
    int N4 = 0;
    // Group the particles
    for (int i = 0; i < N; i++) {
        if (particles[i].pos_x == 0.5 * (LB + RB) && particles[i].pos_y == 0.5 * (DB + UB)) {
            group1_ofParticles[N1] = i;
            N1++;
        } else if (particles[i].pos_x >= square1[0] && particles[i].pos_x <= square1[1] &&
            particles[i].pos_y >= square1[2] && particles[i].pos_y <= square1[3]) {
            group1_ofParticles[N1] = i;
            N1++;
        }else if (particles[i].pos_x > square2[0] && particles[i].pos_x <= square2[1] &&
            particles[i].pos_y >= square2[2] && particles[i].pos_y < square2[3]) {
            group2_ofParticles[N2] = i;
            N2++;
        }else if (particles[i].pos_x >= square3[0] && particles[i].pos_x < square3[1] &&
            particles[i].pos_y > square3[2] && particles[i].pos_y <= square3[3]) {
            group3_ofParticles[N3] = i;
            N3++;
        }else{
            group4_ofParticles[N4] = i;
            N4++;
        }
    }
    // Create and initialize the subtrees
    TNode* subTrees[4];
    for (int i = 0; i < 4; i++) {
        subTrees[i] = NULL;
    }

    #ifdef _OPENMP
    #pragma omp parallel num_threads(n_threads)
    #endif
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                if (N1 > 0) {
                    subTrees[0] = create_new_TNode(-1, square1[0], square1[1], square1[2], square1[3]);
                    subTrees[0]->particle = &particles[group1_ofParticles[0]];
                }
                for (int i = 1; i < N1; i++) {
                    insert(subTrees[0], &particles[group1_ofParticles[i]]);
                }
                preOrder(subTrees[0]);
            }

            #pragma omp section
            {
                if (N2 > 0) {
                    subTrees[1] = create_new_TNode(-1, square2[0], square2[1], square2[2], square2[3]);
                    subTrees[1]->particle = &particles[group2_ofParticles[0]];
                }
                for (int i = 1; i < N2; i++) {
                    insert(subTrees[1], &particles[group2_ofParticles[i]]);
                }
                preOrder(subTrees[1]);
            }

            #pragma omp section
            {
                if (N3 > 0) {
                    subTrees[2] = create_new_TNode(-1, square3[0], square3[1], square3[2], square3[3]);
                    subTrees[2]->particle = &particles[group3_ofParticles[0]];
                }
                for (int i = 1; i < N3; i++) {
                    insert(subTrees[2], &particles[group3_ofParticles[i]]);
                }
                preOrder(subTrees[2]);
            }

            #pragma omp section
            {
                if (N4 > 0) {
                    subTrees[3] = create_new_TNode(-1, square4[0], square4[1], square4[2], square4[3]);
                    subTrees[3]->particle = &particles[group4_ofParticles[0]];
                }
                for (int i = 1; i < N4; i++) {
                    insert(subTrees[3], &particles[group4_ofParticles[i]]);
                }
                preOrder(subTrees[3]);
            }
        }
    }
    

    TNode* tTree = create_new_TNode(-1, LB, RB, DB, UB);
    tTree->is_leaf = 0;
    tTree->particle = malloc(sizeof(PNode));
    for (int i = 0; i < 4; i++) {
        if (subTrees[i] != NULL) {
            tTree->child[i] = subTrees[i];
            tTree->particle->mass += subTrees[i]->particle->mass;
            tTree->particle->pos_x += subTrees[i]->particle->mass * subTrees[i]->particle->pos_x;
            tTree->particle->pos_y += subTrees[i]->particle->mass * subTrees[i]->particle->pos_y;
        }
    }
    tTree->particle->pos_x /= tTree->particle->mass;
    tTree->particle->pos_y /= tTree->particle->mass;

    return tTree;
}

TNode* create_new_TNode(int index, double LB, double RB, double DB, double UB) {
    TNode* new_TNode = malloc(sizeof(TNode));
    if (new_TNode != NULL) {
        for (int i = 0; i < 4; i++) {
            new_TNode->child[i] = NULL;
        }
        new_TNode->is_leaf = 1;
        double mid_x = 0.5 * (LB + RB);
        double mid_y = 0.5 * (DB + UB);
        new_TNode->LB = (index == -1) * LB + (index == 0 || index == 1) * LB + (index == 2 || index == 3) * mid_x;
        new_TNode->RB = (index == -1) * RB + (index == 0 || index == 1) * mid_x + (index == 2 || index == 3) * RB;
        new_TNode->DB = (index == -1) * DB + (index == 0 || index == 2) * DB + (index == 1 || index == 3) * mid_y;
        new_TNode->UB = (index == -1) * UB + (index == 0 || index == 2) * mid_y + (index == 1 || index == 3) * UB;
    } else {
        printf("Memory allocation failed");
    }
    return new_TNode;
}

int insert(TNode* tNode, PNode* particle) {
    double mid_x = 0.5 * (tNode->LB + tNode->RB);
    double mid_y = 0.5 * (tNode->UB + tNode->DB);
    if (tNode->is_leaf) {
        if (((particle->pos_x - tNode->particle->pos_x)<1e-8 && (particle->pos_x - tNode->particle->pos_x)>-1e-8) && 
            ((particle->pos_y - tNode->particle->pos_y)<1e-8 && (particle->pos_y - tNode->particle->pos_y)>-1e-8)) {
            printf("Two particles are detected at the same location and the simulation terminates.\n");
            exit(0);
        }
        tNode->is_leaf = 0;

        int index = (tNode->particle->pos_y > mid_y) + 2 * (tNode->particle->pos_x > mid_x);
    
        tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
        tNode->child[index]->particle = tNode->particle;
        tNode->particle = malloc(sizeof(PNode));

        tNode->particle->mass = tNode->child[index]->particle->mass;
        tNode->particle->pos_x = tNode->child[index]->particle->mass * tNode->child[index]->particle->pos_x;
        tNode->particle->pos_y = tNode->child[index]->particle->mass * tNode->child[index]->particle->pos_y;
    }

    tNode->particle->pos_x += particle->mass * particle->pos_x;
    tNode->particle->pos_y += particle->mass * particle->pos_y;
    tNode->particle->mass += particle->mass;

    int index = (particle->pos_y > mid_y) + 2 * (particle->pos_x > mid_x);

    if (tNode->child[index] == NULL) {
        tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
        tNode->child[index]->particle = particle;
    } else {
        insert(tNode->child[index], particle);
    }
    return 0;
}

void barnesHut(PNode* particle, TNode* tNode, double* fx, double* fy) {
    double theta;

    if (!tNode || tNode->particle == particle) {
        return;
    }

    if (!tNode->is_leaf) {
        double mid_x = 0.5 * (tNode->LB + tNode->RB);
        double mid_y = 0.5 * (tNode->DB + tNode->UB);
        double width = (tNode->RB - tNode->LB);
        theta = width / sqrt((particle->pos_x - mid_x) * (particle->pos_x - mid_x) +
                             (particle->pos_y - mid_y) * (particle->pos_y - mid_y));
    }

    if (tNode->is_leaf || theta <= THETA_MAX) {
        double r_x = particle->pos_x - tNode->particle->pos_x;
        double r_y = particle->pos_y - tNode->particle->pos_y;
        double r_squared = r_x * r_x + r_y * r_y;
        double r_plummer = sqrt(r_squared) + EPSILON_O;
        double force_factor = -G * particle->mass * tNode->particle->mass / (r_plummer * r_plummer * r_plummer);
        *fx += force_factor * r_x;
        *fy += force_factor * r_y;
    } else if (tNode->is_leaf == 0) {
        barnesHut(particle, tNode->child[0], fx, fy);
        barnesHut(particle, tNode->child[1], fx, fy);
        barnesHut(particle, tNode->child[2], fx, fy);
        barnesHut(particle, tNode->child[3], fx, fy);
    }
}

void preOrder(TNode* tNode) {
    if (tNode == NULL || tNode->is_leaf) {
        return;
    } else {
        double inverse_mass = 1.0 / tNode->particle->mass;
        tNode->particle->pos_x *= inverse_mass;
        tNode->particle->pos_y *= inverse_mass;
        preOrder(tNode->child[0]);
        preOrder(tNode->child[1]);
        preOrder(tNode->child[2]);
        preOrder(tNode->child[3]);
    }
}

void destroy(TNode* root) {
    if (root == NULL) {
        return;
    } else if (root->is_leaf) {
        free(root);
    } else {
        destroy(root->child[0]);
        destroy(root->child[1]);
        destroy(root->child[2]);
        destroy(root->child[3]);
        free(root->particle);
        free(root);
    }
}
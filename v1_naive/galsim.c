#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    double x, y, mass, vx, vy, b;
} Particle;

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Usage: N filename nsteps\n");
        return 1;
    }

    int N = atoi(argv[1]);
    int nsteps = atoi(argv[3]);
    double dt = 0.001;
    double G = 100.0 / N;
    double eps = 1e-3;

    Particle *p = malloc(N * sizeof(Particle));
    FILE *f = fopen(argv[2], "rb");
    for (int i = 0; i < N; i++) {
        fread(&p[i].x, 8, 1, f); fread(&p[i].y, 8, 1, f);
        fread(&p[i].mass, 8, 1, f); fread(&p[i].vx, 8, 1, f);
        fread(&p[i].vy, 8, 1, f); fread(&p[i].b, 8, 1, f);
    }
    fclose(f);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    FILE *movie_file = fopen("movie.gal", "wb");
    for (int s = 0; s < nsteps; s++) {
        if (s % 1 == 0) {
            for (int i = 0; i < N; i++) {
                fwrite(&p[i].x, 8, 1, movie_file);
                fwrite(&p[i].y, 8, 1, movie_file);
                fwrite(&p[i].mass, 8, 1, movie_file);
            }
        }
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            double fx = 0, fy = 0;
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                double dx = p[j].x - p[i].x;
                double dy = p[j].y - p[i].y;
                double distSq = dx*dx + dy*dy + eps*eps;
                double invDist = 1.0 / sqrt(distSq);
                double invDist3 = invDist * invDist * invDist;
                fx += G * p[i].mass * p[j].mass * dx * invDist3;
                fy += G * p[i].mass * p[j].mass * dy * invDist3;
            }
            p[i].vx += dt * fx / p[i].mass;
            p[i].vy += dt * fy / p[i].mass;
        }
        for (int i = 0; i < N; i++) {
            p[i].x += dt * p[i].vx;
            p[i].y += dt * p[i].vy;
        }
    }
    fclose(movie_file);

    gettimeofday(&end, NULL);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Naive Simulation took %f seconds\n", elapsed);

    free(p);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dem_discr.h"

clock_t start_;

int n, d, k;
double best = 1.0;
double **points, **subset;
long long int terminal_nodes = 0;

double time_limit = 60 * 30;

double min(double a, double b) {
	if (a < b)
		return a;
	else
		return b;
}

void branch_and_bound(int cur_point, int subset_size, int p[]) {
	if ((double)(clock() - start_) / CLOCKS_PER_SEC > time_limit) {
		printf("Time limit exceeded: %.2lf seconds\n", time_limit);
		exit(0);
	}
	if (subset_size == k) {
		memset(subset, 0, sizeof subset);
		int ind = 0;
		for(int i = 0; i < n; i++) {
			if (p[i] == 1) {
				subset[ind] = points[i];
				ind++;
			}
		}
		double res = 1.0, aux;
		res = oydiscr(subset, d, k, &aux);
		best = min(res, best);
		terminal_nodes++;
		return;
	}
	if (cur_point == n || subset_size + (n - cur_point) < k) {
		return;
	}
	p[cur_point] = 0;
	branch_and_bound(cur_point + 1, subset_size, p);
	p[cur_point] = 1;
	branch_and_bound(cur_point + 1, subset_size + 1, p);
	p[cur_point] = 0;
}

int main(int argc, char** argv) {
	clock_t end;
	double coord;
	k = atoi(argv[1]);
	scanf("%d %d", &n, &d);
	points = (double**) malloc(n * sizeof(double*));
	subset = (double**) malloc(k * sizeof(double*));
	for(int i = 0; i < n; i++) {
		points[i] = (double*) malloc(d * sizeof(double));
		for(int j = 0; j < d; j++) {
			scanf("%lf", &coord);
			points[i][j] = coord;
		}
	}
	int p[n];
	memset(p, 0, n * sizeof(int));
	start_ = clock();
	branch_and_bound(0, 0, p);
	end = clock();
	printf("%lf\n", best);
	printf("Terminal nodes: %lld\n", terminal_nodes);
	printf("Time: %.2lf seconds\n", (double)(end - start_) / CLOCKS_PER_SEC);
	return 0;
}
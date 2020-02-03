#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <map>
#include <vector>

#include "dem_discr.h"

using namespace std;

clock_t start_;

int n, d, k;
double best = 1.0;
double time_limit = 60 * 30;
double **points, **subset, **points_;
int grid[520][520], grid_[520][520], grid_2[520][520];
map<double, int> map_x, map_y;

long long int n_cuts = 0;
int depth_cuts = 0;
long long int terminal_nodes = 0;
long long int bound_A = 0;
long long int bound_B = 0;
int last_bound = 0;

double min(double a, double b) {
	if (a < b)
		return a;
	else
		return b;
}

double max(double a, double b) {
	if (a > b)
		return a;
	else
		return b;
}

bool pair_comparator(pair<int, double> a, pair<int, double> b) {
	return a.second < b.second;
}

double euclidean_distance(double *a, double *b) {
	double dist = 0.0;
	for(int i = 0; i < d; i++) {
		dist += pow(a[i] - b[i], 2);
	}
	return sqrt(dist);
}

double lower_bound(int cur_point, int p[]) {
	double remaining_discrA = INT_MIN, local_discrA, cur_discrA = INT_MIN;
	double remaining_discrB = INT_MIN, local_discrB, cur_discrB = INT_MIN;
	pair<double, double> result = make_pair(0.0, 0.0);
	int n_points = 0;

	//for(int i = 0; i < cur_point; i++)
	//	n_points += p[i];

	//Current discrepancy
	for(int i = 0; i < cur_point; i++) {
		if (p[i] == 1) {
			local_discrA = INT_MIN;
			local_discrB = INT_MIN;
			local_discrA = max(local_discrA, (points[i][0] * points[i][1]) - (grid[map_y[points[i][1]]][map_x[points[i][0]]] / (double) k));
			local_discrA = max(local_discrA, (points[i][1]) - (grid[map_y[points[i][1]]][map_x[1]] / (double) k));
			local_discrA = max(local_discrA, (points[i][0]) - (grid[map_y[1]][map_x[points[i][0]]] / (double) k));
			local_discrB = max(local_discrB, (grid_2[map_y[points[i][1]]][map_x[points[i][0]]] / (double) k) - (points[i][0] * points[i][1]));
			//local_discrB = max(local_discrB, (grid_2[map_y[points[i][1]]][map_x[1]] / (double) k) - (points[i][1]));
			//local_discrB = max(local_discrB, (grid_2[map_y[1]][map_x[points[i][0]]] / (double) k) - (points[i][0]));
			for(int j = i + 1; j < cur_point; j++) {
				if (p[j] == 1) {
					local_discrA = max(local_discrA, (points[i][0] * points[j][1]) - (grid[map_y[points[j][1]]][map_x[points[i][0]]] / (double) k));
					local_discrA = max(local_discrA, (points[j][0] * points[i][1]) - (grid[map_y[points[i][1]]][map_x[points[j][0]]] / (double) k));
					local_discrB = max(local_discrB, (grid_2[map_y[points[j][1]]][map_x[points[i][0]]] / (double) k) - (points[i][0] * points[j][1]));
					local_discrB = max(local_discrB, (grid_2[map_y[points[i][1]]][map_x[points[j][0]]] / (double) k) - (points[j][0] * points[i][1]));
				}
			}
			cur_discrA = max(cur_discrA, local_discrA);
			cur_discrB = max(cur_discrB, local_discrB);
		}
	}
	//NEW
	//if (n_points == k) {
	//	return max(cur_discrA, cur_discrB);
	//}
	//Remaining discrepancy
	for (int i = cur_point; i < n; i++) {
		local_discrA = INT_MIN;
		local_discrB = INT_MIN;
		local_discrA = max(local_discrA, (points[i][0] * points[i][1]) - (grid[map_y[points[i][1]]][map_x[points[i][0]]] / (double) k));
		local_discrA = max(local_discrA, (points[i][1]) - (grid[map_y[points[i][1]]][map_x[1]] / (double) k));
		local_discrA = max(local_discrA, (points[i][0]) - (grid[map_y[1]][map_x[points[i][0]]] / (double) k));
		local_discrB = max(local_discrB, (grid_2[map_y[points[i][1]]][map_x[points[i][0]]] / (double) k) - (points[i][0] * points[i][1]));
		//local_discrB = max(local_discrB, (grid_2[map_y[points[i][1]]][map_x[1]] / (double) k) - (points[i][1]));
		//local_discrB = max(local_discrB, (grid_2[map_y[1]][map_x[points[i][0]]] / (double) k) - (points[i][0]));
		remaining_discrA = min(remaining_discrA, local_discrA);
		remaining_discrB = min(remaining_discrB, local_discrB);
	}
	if (cur_point != 0) {
		result.first = max(cur_discrA, remaining_discrA);
		result.second = max(cur_discrB, remaining_discrB);
	}
	else {
		result.first = remaining_discrA;
		result.second = remaining_discrB;
	}
	if (result.first >= result.second)
		last_bound = 1;
	else
		last_bound = 2;
	return max(result.first, result.second);
}

void bound_init() {
	vector<double> aux_x, aux_y;

	memset(grid, 0, sizeof grid);
	memset(grid_, 0, sizeof grid_);
	memset(grid_2, 0, sizeof grid_2);
	memset(subset, 0, sizeof subset);
	for(int i = 0; i < n; i++) {
		aux_x.push_back(points[i][0]);
		aux_y.push_back(points[i][1]);
	}
	aux_x.push_back(0);
	aux_x.push_back(1);
	aux_y.push_back(0);
	aux_y.push_back(1);
	sort(aux_x.begin(), aux_x.begin() + aux_x.size());
	sort(aux_y.begin(), aux_y.begin() + aux_y.size());
	
	for(int i = 0; i < n + 2; i++) {
		map_x.insert(pair<double, int>(aux_x[i], i));
		map_y.insert(pair<double, int>(aux_y[i], n + 2 - 1 - i));
	}

	for(int i = 0; i < n; i++) {
		for(int asd = map_y[points[i][1]] - 1; asd >= 0; asd--) {
			for(int j = map_x[points[i][0]] + 1; j < n + 2; j++) {
				if (grid[asd][j] < k) {
					grid[asd][j]++;
				}
				grid_[asd][j]++;
			}
		}
	}
	//SORT
	
	vector<int> layers[520];
	int n_layers = 0;
	vector<int> point_list;
	for(int i = 0; i < n; i++)
		point_list.push_back(i);
	int counter = 0;
	int points_sorted[520];
	memset(points_sorted, 0, sizeof points_sorted);
	while(counter != n) {
		for(int i = 0; i < n; i++) {
			if (points_sorted[i] == 0) {
				int is_max = 1;
				for(int j = 0; j < point_list.size(); j++) {
					if (points[j][0] > points[i][0] && points[j][1] > points[i][1] && points_sorted[j] == 0) {
						is_max = 0;
					}
				}
				if(is_max == 1) {
					layers[n_layers].push_back(i);
					points_sorted[i] = 1;
					counter++;
				}
			}
		}
		n_layers++;
	}
	counter = 0;
	for(int i = n_layers - 1; i >= 0; i--) {
		for(int j = 0; j < layers[i].size(); j++) {
			points_[counter] = points[layers[i][j]];
			counter++;
		}
	}
	points = points_;
	
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
		if(res < best) {
			printf("New best %lf at %.2lf seconds\n", res, (double)(clock() - start_) / CLOCKS_PER_SEC);
		}
		best = min(res, best);
		terminal_nodes++;
		return;
	}
	if (cur_point == n || subset_size + (n - cur_point) < k) {
		return;
	}
	if(lower_bound(cur_point, p) > best) {
		n_cuts++;
		depth_cuts += cur_point;
		if(last_bound == 1)
			bound_A++;
		else if(last_bound == 2)
			bound_B++;
		return;
	}

	p[cur_point] = 1;
	for(int asd = map_y[points[cur_point][1]]; asd >= 0; asd--) {
		for(int j = map_x[points[cur_point][0]]; j < n + 2; j++) {
			grid_2[asd][j]++;
		}
	}

	branch_and_bound(cur_point + 1, subset_size + 1, p);

	for(int asd = map_y[points[cur_point][1]]; asd >= 0; asd--) {
		for(int j = map_x[points[cur_point][0]]; j < n + 2; j++) {
			grid_2[asd][j]--;
		}
	}
	p[cur_point] = 0;

	for(int asd = map_y[points[cur_point][1]] - 1; asd >= 0; asd--) {
		for(int j = map_x[points[cur_point][0]] + 1; j < n + 2; j++) {
			grid_[asd][j]--;
			if(grid_[asd][j] < k) {
				grid[asd][j] = grid_[asd][j];
			}
		}
	}
	
	branch_and_bound(cur_point + 1, subset_size, p);

	for(int asd = map_y[points[cur_point][1]] - 1; asd >= 0; asd--) {
		for(int j = map_x[points[cur_point][0]] + 1; j < n + 2; j++) {
			grid_[asd][j]++;
			if(grid_[asd][j] >= k) {
				grid[asd][j] = k;
			}
			else
				grid[asd][j] = grid_[asd][j];
		}
	}	
}

void best_init() {
	int *chosen;
	int counter, next_point;
	double max_of_dists, point_dist;
	chosen = (int*) malloc(n * sizeof(int));
	for(int starting_point = 0; starting_point < n; starting_point++) {
		memset(chosen, 0, n * sizeof(int));
		memset(subset, 0, k * sizeof(double*));
		subset[0] = points[starting_point];
		chosen[starting_point] = 1;
		counter = 1;
		while(counter != k) {
			max_of_dists = INT_MIN;
			for(int i = 0; i < n; i++) {
				if (chosen[i] == 0) {
					point_dist = INT_MAX;
					for(int j = 0; j < n; j++) {
						if (i != j && chosen[j] == 1) {
							point_dist = min(point_dist, euclidean_distance(points[i], points[j]));
						}
					}
					if (point_dist > max_of_dists) {
						next_point = i;
						max_of_dists = point_dist;
					}
				}
			}
			subset[counter] = points[next_point];
			chosen[next_point] = 1;
			counter++;
		}
		double res = 1.0, aux;
		res = oydiscr(subset, d, k, &aux);
		best = min(res, best);
	}
	printf("Starting best: %lf\n", best);
	memset(subset, 0, k * sizeof(double*));
}

//0.176099
//Time: 201.32 seconds

//0.176099
//Time: 93.03 seconds		

int main(int argc, char** argv) {
	clock_t end;
	double coord;
	k = atoi(argv[1]);
	scanf("%d %d", &n, &d);
	points = (double**) malloc(n * sizeof(double*));
	points_ = (double**) malloc(n * sizeof(double*));
	subset = (double**) malloc(k * sizeof(double*));
	for(int i = 0; i < n; i++) {
		points[i] = (double*) malloc(d * sizeof(double));
		points_[i] = (double*) malloc(d * sizeof(double));
		for(int j = 0; j < d; j++) {
			scanf("%lf", &coord);
			points[i][j] = coord;
		}
	}
	int p[n];
	memset(p, 0, n * sizeof(int));
	
	start_ = clock();
	bound_init();

	best_init();
	branch_and_bound(0, 0, p);
	end = clock();
	printf("%lf\n", best);
	printf("Terminal nodes: %lld\n", terminal_nodes);
	printf("%lld cuts at average depth of %.2lf\n", n_cuts, depth_cuts / (double) n_cuts);
	printf("Bound A worked %lld times\n", bound_A);
	printf("Bound B worked %lld times\n", bound_B);
	printf("Time: %.2lf seconds\n", (double)(end - start_) / CLOCKS_PER_SEC);
	return 0;
}
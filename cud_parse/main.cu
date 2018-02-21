#include <iostream>
#include "Graph_modules.h"
//#include "Assembly_modules.h"
#include "Indexes.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "Parameters.h"
using namespace std;

#define N 839641 //513361830 //839641
#define MAX_SEQ_LEN 300
#define MAX_LOCATIONS 60

__global__
void process(int *device_targets, float *device_MTM_scores, int num_targets, int cnt_lines,
		float *table_A, float *table_T, float *table_C, float *table_G, int *max_path, 
		int *map_loc, int *device_sixmers, int *device_sixmer_idx, 
		int *out_test) {
	int i, j, k, h, t;
	int start, end, target_idx_start, map_loc_start, map_loc_end;
	int table_idx_start, max_idx_start;
	int diff, loc_cnt;
	float div;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	bool flag;
	for (i = index; i < cnt_lines; i += stride){
		start = device_sixmer_idx[i];
		end = device_sixmer_idx[i+1];
		out_test[i] = 0;
		for(t=0; t<num_targets; t++){
			device_MTM_scores[t] = 0.0f;
		}
		for(t=0; t<num_targets; t++){
			target_idx_start = t*MAX_SEQ_LEN;
			j=start;
			while(j<end-6){
				for(k=0; k<MAX_SEQ_LEN; k++){
					if(device_targets[target_idx_start+k+6] == -1) break;
					flag = true;
					for(h=0; h<6; h++){
						if(device_sixmers[j+h] != device_targets[target_idx_start+k]){
							flag = false;
							break;
						}
					}
					if(flag)
						device_MTM_scores[t] += 1.0f;
				}
			j = j + 6;
			}
		}
		int top_match_idx = -1;
		float top_match_score = -1.0f;
		for(t=0; t<num_targets; t++){
			if(top_match_score < device_MTM_scores[t]){
				top_match_idx = t;
				top_match_score = device_MTM_scores[t];
			}
		}

		target_idx_start = top_match_idx*MAX_SEQ_LEN;
		j=start;
		while(j<end-6){
			map_loc_start = i*MAX_LOCATIONS;
			map_loc_end = map_loc_start;
			for(t=0; t<MAX_LOCATIONS; t++){
				if(map_loc[map_loc_start+t] == -1) break;
				map_loc[map_loc_start+t] = -1;
			}
			for(k=0; k<MAX_SEQ_LEN; k++){
				if(device_targets[target_idx_start+k+6] == -1) break;
				diff = 0;
				for(h=0; h<6; h++){
					if(device_sixmers[j+h] != device_targets[target_idx_start+k])
						diff++;
					if(diff>1) break;
				}
				if(diff==0 || diff==1){
					for(h=0; h<6; h++){
						map_loc[map_loc_end+h] = device_sixmers[j+h];
					}
					map_loc_end = map_loc_end+6;
					loc_cnt++;
				}
			}
			div = 1.f/(float)loc_cnt;
			table_idx_start = i*MAX_SEQ_LEN;
			for(t=0; t<map_loc_end-map_loc_start; t++){
				if(map_loc[map_loc_start+t] == 0)
					table_A[table_idx_start+t] += div;
				else if(map_loc[map_loc_start+t] == 1)
					table_T[table_idx_start+t] += div;
				else if(map_loc[map_loc_start+t] == 2)
					table_C[table_idx_start+t] += div;
				else if(map_loc[map_loc_start+t] == 3)
					table_G[table_idx_start+t] += div;
			}
		j = j + 6;
		}

		max_idx_start = i*MAX_SEQ_LEN;
		for(t=0; t<MAX_SEQ_LEN; t++){
			max_path[max_idx_start+t] = table_A[max_idx_start+t];
			if(max_path[max_idx_start+t] < table_T[max_idx_start+t])
				max_path[max_idx_start+t] = table_T[max_idx_start+t];
			else if(max_path[max_idx_start+t] < table_C[max_idx_start+t])
				max_path[max_idx_start+t] = table_C[max_idx_start+t];
			else if(max_path[max_idx_start+t] < table_G[max_idx_start+t])
				max_path[max_idx_start+t] = table_G[max_idx_start+t];
		}

	}
}

int main() {
	Indexes *index_obj = Indexes::getInstance();
	Param* param = Param::getInstance();
    	map<char, int> ATCG_dic;
	ATCG_dic['A'] = 0;
	ATCG_dic['T'] = 1;
	ATCG_dic['C'] = 2;
	ATCG_dic['G'] = 3;
	int i, j, map_idx = 0;
	index_obj->buildIndexes();
	int num_targets = index_obj->target_sequences.size();
	int array_size = MAX_SEQ_LEN*num_targets;
	int cpu_targets[array_size];
	int outtest_targets[array_size];
	clock_t begin = clock();

	for (i = 0; i < array_size; i++) {
		cpu_targets[i] = -1;
	}
	string target_map[num_targets];
	map<string, string>::iterator it;
	for (it = index_obj->target_sequences.begin(); it != index_obj->target_sequences.end(); ++it) {
		string gene_str = it->first;
		string seq_str = it->second;
		target_map[map_idx] = gene_str;
		for(i=0; i<seq_str.size(); i++){
			cpu_targets[map_idx*MAX_SEQ_LEN+i] = ATCG_dic[seq_str[i]];
		}
		cpu_targets[map_idx*MAX_SEQ_LEN+seq_str.size()] = -1;
		map_idx++;
	}

	int *device_targets;
	cudaMalloc((void **)&device_targets, array_size*sizeof(int));
	cudaMemcpy(device_targets, cpu_targets, array_size*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(outtest_targets, device_targets, array_size*sizeof(int), cudaMemcpyDeviceToHost);

	float *device_MTM_scores;
	cudaMalloc((void **)&device_MTM_scores, num_targets*sizeof(float));
	cudaMemset(device_MTM_scores, 0.0f, num_targets*sizeof(float));



        string sequence_file = param->get_sequence_file();
	ifstream inSequnceFile(sequence_file.c_str());
	bool flag_keep_reading = true;
	string origin_gene;



	vector<string> data;
	vector<int> sixmer_list;
	vector<int> sixmer_idx;
	sixmer_idx.push_back(0);
	string line;
	int cnt_lines = 0;
        while (!inSequnceFile.eof()){
		getline(inSequnceFile, line);
		if (line[0] == 'F' || line[0] == '#' || line == "") // line.startswith('Features'):
		    continue;
		index_obj->tokenizeLine(line, data, '\t');
		if (data.size() < 5)
		    continue;
		string features = data[0];
		string fov = data[1];
		string x = data[2];
		string y = data[3];
		string xy_str = fov + "_" + x + "." + y;
		//        param->debugOutFile << xy_str << "=" << data.size() << " : ";
		for (i=4; i<data.size(); i++) {
			if (data[i] == "0" || data[i].size() != 6){
				continue;
			}
			string seq = data[i];
			if(seq.find('N')!=string::npos)
				continue;
			if(seq.size() == 6){
				for (j=0; j<6; j++)
					sixmer_list.push_back(ATCG_dic[seq[j]]);
			}
		}
		sixmer_idx.push_back(sixmer_list.size());
		cnt_lines++;
		data.clear();
	}

	int num_sixmers = sixmer_list.size();
	cout << "all_sixmer_array: " << num_sixmers << endl;
	int *cpu_sixmers = &sixmer_list[0];
	int *device_sixmers;
	cudaMalloc((void **)&device_sixmers, num_sixmers*sizeof(int));
	cudaMemcpy(device_sixmers, cpu_sixmers, num_sixmers*sizeof(int), cudaMemcpyHostToDevice);

	int num_sixmer_idx = sixmer_idx.size();
	cout << "all_sixmer_idx: " << num_sixmer_idx << endl;
	int *cpu_sixmer_idx = &sixmer_idx[0];
	int *device_sixmer_idx;
	cudaMalloc((void **)&device_sixmer_idx, num_sixmer_idx*sizeof(int));
	cudaMemcpy(device_sixmer_idx, cpu_sixmer_idx, num_sixmer_idx*sizeof(int), cudaMemcpyHostToDevice);
	
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "time elapsed for parsing: " << elapsed_secs << endl;
	clock_t begin_proc = clock();

	float *table_A;
	cudaMalloc((void **)&table_A, N*MAX_SEQ_LEN*sizeof(float));
	cudaMemset(table_A, 0.0f, array_size*sizeof(float));
	float *table_T;
	cudaMalloc((void **)&table_T, N*MAX_SEQ_LEN*sizeof(float));
	cudaMemset(table_T, 0.0f, array_size*sizeof(float));
	float *table_C;
	cudaMalloc((void **)&table_C, N*MAX_SEQ_LEN*sizeof(float));
	cudaMemset(table_C, 0.0f, array_size*sizeof(float));
	float *table_G;
	cudaMalloc((void **)&table_G, N*MAX_SEQ_LEN*sizeof(float));
	cudaMemset(table_G, 0.0f, array_size*sizeof(float));
	int *max_path;
	cudaMalloc((void **)&max_path, N*MAX_SEQ_LEN*sizeof(int));
	cudaMemset(max_path, 0, array_size*sizeof(int));

	int *map_loc;
	cudaMalloc((void **)&map_loc, N*MAX_LOCATIONS*sizeof(int));
	cudaMemset(map_loc, -1, N*MAX_LOCATIONS*sizeof(int));


	int *out_test;
	cudaMalloc((void **)&out_test, num_sixmer_idx*sizeof(int));
	int blockSize = 256;
	int numBlocks = (N + blockSize - 1) / blockSize;
	process<<<numBlocks, blockSize>>>(device_targets, device_MTM_scores, num_targets, cnt_lines,
		table_A, table_T, table_C, table_G, max_path, map_loc,
		device_sixmers, device_sixmer_idx, out_test);

	int hb[num_sixmer_idx];
	cudaMemcpy(hb, out_test, num_sixmer_idx*sizeof(int), cudaMemcpyDeviceToHost);


	//for (i = 0; i < num_sixmer_idx; i++) {
	//	cout << hb[i] << endl;
	//}
	//for (i = 0; i < sizeof(device_sixmers)/sizeof(device_sixmers[0]); i++) {
	//	cout << hb[i] << " ";
	//}

    	cudaFree(device_targets);
    	cudaFree(device_sixmers);
    	cudaFree(device_sixmer_idx);
    	cudaFree(out_test);

	clock_t end_proc = clock();
	double elapsed_secs_proc = double(end_proc - begin_proc) / CLOCKS_PER_SEC;
	cout << "time elapsed for processing: " << elapsed_secs_proc << endl;
	return 0;
}

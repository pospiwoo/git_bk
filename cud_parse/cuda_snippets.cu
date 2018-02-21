


/*
	thrust::device_vector<char> gene_d(MAX_GENE_STR * num_targets);
	thrust::copy(gene_d.begin(), gene_d.end(), gene.begin());
	thrust::device_vector<char> seq_d(MAX_SEQ_LEN * num_targets);
	thrust::copy(seq_d.begin(), seq_d.end(), seq.begin());

	//int MAX_GENE_STR = 30;
	//thrust::host_vector<char> gene(MAX_GENE_STR * num_targets);
	//thrust::host_vector<char> seq(MAX_SEQ_LEN * num_targets);

	//for(int i=0; i<100; i++){
		h_vec[0] = 'R';
		h_vec[1] = 'T';
	//}
	thrust::device_vector<char> d_vec = h_vec;
	thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());


	Indexes *instance_device;
	const unsigned int bytes = sizeof(instance_host);
	cudaMalloc((Indexes**)&instance_device,bytes);
	cudaMemcpy(instance_device,instance_host,bytes,cudaMemcpyHostToDevice);

    //instance_device->buildIndexes();
    //print_target_index(instance_device);

	map<string, string>::iterator it;
	for (it = index_obj->target_sequences.begin(); it != index_obj->target_sequences.end(); ++it) {
		string gene_str = it->first;
		string seq_str = it->second;

		strcpy(gene, gene_str.c_str());
		strcpy(seq, seq_str.c_str());

		cudaMalloc((void**)&gene,bytes);
		cudaMemcpy(instance_device,instance_host,bytes,cudaMemcpyHostToDevice);

	}

*/

/*
    clock_t begin = clock();

    Assembly *assm_obj = new Assembly();
    assm_obj->process();

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "time elapsed: " << elapsed_secs << endl;
*/

/*
__host__ __device__  void print_target_index(Indexes *index_obj){
	cout << "1111111111111111" << endl;
	map<string, map<string, vector<int> > >::iterator it;
	map<string, vector<int> >::iterator it_1;

	cout << "222222" << endl;
	for (it = index_obj->sixmer_dic_target.begin(); it != index_obj->sixmer_dic_target.end(); ++it) {
	cout << "333333333" << endl;
	    string key_seq = it->first;
	    map<string, vector<int> > gene_data = it->second;
	    cout << key_seq << " - " ;
	    for (it_1 = gene_data.begin(); it_1 != gene_data.end(); ++it_1) {
		string gene_str = it_1->first;
		vector<int> idx_list = it_1->second;
		cout << gene_str << ":" ;
		if(idx_list.size() == 1)
		    cout << idx_list[0];
		else
		    for(int i=0; i<idx_list.size(); i++)
		        cout << idx_list[i] << ",";
		cout << " ";
	    }
	    cout << endl;
	}
}
*/

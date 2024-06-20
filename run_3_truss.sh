#! /bin/bash
./SubgraphMatching.out -q ../tests/3_truss_wo_label.graph -d /home/kars1/Parallel_computation/dataset/com-amazon.ungraph.txt -output ../analysis/amazon_3_truss_si.txt
./SubgraphMatching.out -q ../tests/3_truss_wo_label.graph -d /home/kars1/Parallel_computation/dataset/com-dblp.ungraph.txt -output ../analysis/dblp_3_truss_si.txt
./SubgraphMatching.out -q ../tests/3_truss_wo_label.graph -d /home/kars1/Parallel_computation/dataset/com-youtube.ungraph.txt -output ../analysis/youtube_3_truss_si.txt
./SubgraphMatching.out -q ../tests/3_truss_wo_label.graph -d /home/kars1/Parallel_computation/dataset/com-lj.ungraph.txt -output ../analysis/lj_3_truss_si.txt
./SubgraphMatching.out -q ../tests/3_truss_wo_label.graph -d /home/kars1/Parallel_computation/dataset/com-orkut.ungraph.txt -output ../analysis/orkut_3_truss_si.txt

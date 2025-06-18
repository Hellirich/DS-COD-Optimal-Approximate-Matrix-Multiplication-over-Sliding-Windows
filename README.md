# DS-COD-Optimal-Approximate-Matrix-Multiplication-over-Sliding-Windows

You can run `demo.m` to observe the performance of the DS-COD algorithm and baselines. The `algorithms` folder contains the main implementations of hDS-COD and aDS-COD, along with their time-window variants. The core implementation of DS-COD is located in `class/FastDumpSnapshotCod.m`. Baseline methods (EH-COD and DI-COD) are available in the `eh_di_cod` folder. The `datasets` folder includes the cross-lingual datasets (APR, PAN11, EURO) used in the paper. Synthetic datasets can be generated using the `generate_synthetic_data.m` script. The `output` folder is used to store the algorithm results.


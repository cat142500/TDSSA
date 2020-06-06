# Datasets and Code of Our KDD20 Paper
Please contact the following author if you need further clarification.<br />
Author: Yue Wang<br />
Contact: cat142500@hotmail.com

## SUMMARY
This folder includes the datasets and the source code for the following paper:

Title: Truth Discovery against Strategic Sybil Attack<br />
Authors: Yue Wang, Ke Wang and Chunyan Miao<br />
Conference: KDD 2020<br />
The link to our paper will be provided later.

Please cite this work if either the datasets or the code are used by your works.

Datasets: We use two real datasets (NLP and DOG) and generate a synthetic dataset (SYN) for testing TDSSA. The details are provided in "DATASET DESCRIPTION".

Code: We implement our TDSSA framework in Java. The details are provided in "CODE DESCRIPTION", along with a demo to run TDSSA on Windows.

## DATASET DESCRIPTION

Each of the two real datasets (NLP and DOG) and the synthetic dataset (SYN) dataset is characterized by five data parameters, including task number (N), worker number (M), worker number per task (K), optional label size (L), and average worker accuracy (theta). The following table summarizes the data parameters of the NLP and DOG datasets and the tested settings of the SYN dataset.

| Data Parameters | NLP | DOG | SYN |
| ---------- | ---- | ---- | ----- |
| Task Number (N) | 1000 | 807 | 1000 to 9000 |
| Worker Number (M) | 85 | 109 | 100 to 900 |
| Worker Number per Task (K) | 20 | 10 | 5 to 20 |
| Optional Label Size (L) | 2 | 4 | 2 to 5 |
| Average Worker Accuracy (theta) | 0.8 | 0.7 | 0.7 to 1 |

### (1) RAW DATA

The original NLP dataset contained in the "NLP" folder provides four files:
1. "truth.csv" has a (task id, true label) pair in each line, representing the true label of a normal task.
2. "answer.csv" has a (task id, worker id, worker label) triple in each line, representing a worker's label on a normal task.
3. "quali_truth.csv" has a (task id, true label) pair in each line, representing the true label of a golden task. The total number of golden tasks is 20.
4. "quali.csv" has a (task id, worker id, worker label) triple in each line, representing a worker's label on a golden task.

The original DOG dataset contained in the "DOG" folder provides two files:
1. "truth.csv" has a (task id, true label) pair in each line, representing the true label of a normal task.
2. "answer.csv" has a (task id, worker id, worker label) triple in each line, representing a worker's label on a normal task.

### (2) FORMATED INPUT DATA

Based on the above raw data, we generate an "input.txt" file and a "golden.txt" file for the input of TDSSA. For each run, we also generate an "attact.txt" file and an "order.txt" file, simulating a randomized replacement of independent workers into Sybil workers and a randomized requesting order of workers. A Sybil attack is characterized by three attack parameters, including Sybil proportion (mu), deviation probability (epsilon) and attacker number (lambda). The following table summarizes the meaning and the tested setting of each attack parameter.

| Attack Parameters | Description | Setting |
| ----------------- | ----------- | ------- |
| Sybil Proportion (mu) | Percentage of Sybil workers | 0.0 to 0.8 |
| Deviation Probability (epsilon) | Probability for Sybil workers to deviate from sharing | 0.0 to 0.3 |
| Attacker Number (K) | Number of Sybil attackers | 1 to 4 |

#### i. "input.txt"
The first line contains the worker number (M), task number (N), optional label size (L) and worker number per task (K). Each of the rest lines contains a task id, the true label of the task and the (worker id, label) pairs for the task.

#### ii. "golden.txt"
The first line contains the number of initial golden tasks. The second line contains the (task id, true label) pairs of golden tasks. Each of the rest lines contains a worker id and the (task id, label) pairs of the worker on all golden tasks.

#### iii. "attack.txt"
The first line contains the Sybil proportion (mu), deviation probability (epsilon) and attacker number (lambda). For each attcker, the following two lines will be generated: one line contains the attacker id, the total number of normal and golden tasks and the (task id, label) pairs that indicate the attacker's randomized label on each task; the other line contains the attacker id, the number of Sybil workers controlled by the attacker and the id of independent workers replaced into Sybil workers of the attacker.

#### iv. "order.txt"
Each line contains a worker id, representing a request from the worker.

## CODE DESCRIPTION

### (1) DATA PREPROCESSING

The "Preprocess.java" file is for preprecessing raw data into formated input data as described above. To compile the file, user can run the following command:

javac Preprocess.java

If the information of golden tasks is provided, like the NLP dataset, the user needs to specify the dataset name, the number of runs, the Sybil proportion (mu), the deviation probability (epsilon) and the attacker number (lambda). The following command preprocess the NLP dataset into 50 runs with mu = 0.5, epsilon = 0.1 and lambda = 1.

java Preprocess NLP 50 0.5 0.1 1

If the information of golden tasks is not provided, like the DOG dataset, the user also needs to specify the number of golden tasks. The following command preprocess the DOG dataset into 50 runs with mu = 0.5, epsilon = 0.1, lambda = 1 and 20 initial golden tasks.

java Preprocess DOG 50 0.5 0.1 1 20

For the synthetic dataset, like the SYN dataset, the user needs to specify six additional parameters, including the number of golden tasks, the task number (N), the worker number (M), the optional label size (L), the worker number per task (K) and the average worker accuracy (theta). The following command preprocess the SYN dataset into 50 runs with mu = 0.5, epsilon = 0.1, lambda = 1, 20 initial golden tasks, N = 5000, M = 500, L = 4, K = 10 and theta = 0.8.

java Preprocess SYN 50 0.5 0.1 1 20 5000 500 4 10 0.8

### (2) RUNNING TDSSA

To compile the "TDSSA.java" file, user can run the following command:

javac TDSSA.java

To run the program, user needs to specify the dataset name, number of runs, batch condition (B), initial golden assignment probability (alpha), Sybil threshold (tau) and reliability threshold (delta). The following table summarizes the last four parameters used for our TDSSA framework.

| TDSSA Parameters | Description | Setting |
| ---------------- | ----------- | ------- |
| Batch Condition (B) | Percentage of Sybil workers | 5 to 20 |
| Initial Golden Assignment Probability (alpha) | Initial probability to assign a golden task to a new worker | 0.1 to 0.9 |
| Sybil Threshold (tau) | Threshold of Sybil score for banning a worker | 0.6 to 1 |
| Reliability Threshold (delta) | Threshold of reliability score for trusting a worker | 0.6 to 1 |

Taking the DOG dataset for example, the following command will run TDSSA 50 times on the DOG dataset with B = 10, alpha = 0.5, tau = 0.8 and delta = 0.8.

java TDSSA DOG 50 10 0.5 0.8 0.8

The result will be save in a "result.txt" file, showing the following four evaluation metrics:
1. A-Accuracy: the aggregation accuracy, i.e., the percentage of tasks whose aggregated label is identical to the true label.
2. E-Number: the number of golden tasks identified by the attacker.
3. T-Cost: the avergae number of golden tasks used for testing a worker.
4. Time: the running time of TDSSA.

Here is the output of the above command.

Output:<br />
Run 1 --- A-Accuracy:0.7026022304832714  E-Number:0.0  T-Cost:1.834862385321101  Time:398ms<br />
Run 2 --- A-Accuracy:0.7596034696406444  E-Number:0.0  T-Cost:1.834862385321101  Time:326ms<br />
Run 3 --- A-Accuracy:0.7050805452292441  E-Number:0.0  T-Cost:1.834862385321101  Time:304ms<br />
Run 4 --- A-Accuracy:0.6964064436183395  E-Number:0.0  T-Cost:1.834862385321101  Time:277ms<br />
Run 5 --- A-Accuracy:0.7224287484510533  E-Number:0.0  T-Cost:1.834862385321101  Time:312ms<br />
Run 6 --- A-Accuracy:0.7633209417596035  E-Number:0.0  T-Cost:1.834862385321101  Time:331ms<br />
Run 7 --- A-Accuracy:0.724907063197026  E-Number:0.0  T-Cost:1.834862385321101  Time:287ms<br />
Run 8 --- A-Accuracy:0.7633209417596035  E-Number:0.0  T-Cost:1.834862385321101  Time:308ms<br />
Run 9 --- A-Accuracy:0.7434944237918215  E-Number:0.0  T-Cost:1.834862385321101  Time:288ms<br />
Run 10 --- A-Accuracy:0.7026022304832714  E-Number:0.0  T-Cost:1.834862385321101  Time:263ms<br />
Run 11 --- A-Accuracy:0.6951672862453532  E-Number:0.0  T-Cost:1.834862385321101  Time:240ms<br />
Run 12 --- A-Accuracy:0.7137546468401487  E-Number:0.0  T-Cost:1.834862385321101  Time:295ms<br />
Run 13 --- A-Accuracy:0.7001239157372986  E-Number:0.0  T-Cost:1.834862385321101  Time:287ms<br />
Run 14 --- A-Accuracy:0.7385377942998761  E-Number:0.0  T-Cost:1.834862385321101  Time:258ms<br />
Run 15 --- A-Accuracy:0.7187112763320942  E-Number:0.0  T-Cost:1.834862385321101  Time:287ms<br />
Run 16 --- A-Accuracy:0.758364312267658  E-Number:0.0  T-Cost:1.834862385321101  Time:294ms<br />
Run 17 --- A-Accuracy:0.7001239157372986  E-Number:0.0  T-Cost:1.834862385321101  Time:219ms<br />
Run 18 --- A-Accuracy:0.6765799256505576  E-Number:0.0  T-Cost:1.834862385321101  Time:307ms<br />
Run 19 --- A-Accuracy:0.7794299876084263  E-Number:0.0  T-Cost:1.834862385321101  Time:306ms<br />
Run 20 --- A-Accuracy:0.7162329615861215  E-Number:0.0  T-Cost:1.834862385321101  Time:269ms<br />
Run 21 --- A-Accuracy:0.7038413878562577  E-Number:0.0  T-Cost:1.834862385321101  Time:268ms<br />
Run 22 --- A-Accuracy:0.7744733581164808  E-Number:0.0  T-Cost:1.9174311926605505  Time:315ms<br />
Run 23 --- A-Accuracy:0.6926889714993805  E-Number:0.0  T-Cost:1.834862385321101  Time:271ms<br />
Run 24 --- A-Accuracy:0.781908302354399  E-Number:0.0  T-Cost:1.834862385321101  Time:260ms<br />
Run 25 --- A-Accuracy:0.6902106567534076  E-Number:0.0  T-Cost:1.834862385321101  Time:269ms<br />
Run 26 --- A-Accuracy:0.7187112763320942  E-Number:0.0  T-Cost:1.834862385321101  Time:252ms<br />
Run 27 --- A-Accuracy:0.7806691449814126  E-Number:0.0  T-Cost:1.834862385321101  Time:309ms<br />
Run 28 --- A-Accuracy:0.7608426270136307  E-Number:0.0  T-Cost:1.834862385321101  Time:257ms<br />
Run 29 --- A-Accuracy:0.7125154894671624  E-Number:0.0  T-Cost:1.834862385321101  Time:276ms<br />
Run 30 --- A-Accuracy:0.7137546468401487  E-Number:0.0  T-Cost:1.834862385321101  Time:290ms<br />
Run 31 --- A-Accuracy:0.7026022304832714  E-Number:0.0  T-Cost:1.834862385321101  Time:298ms<br />
Run 32 --- A-Accuracy:0.6840148698884758  E-Number:0.0  T-Cost:1.834862385321101  Time:271ms<br />
Run 33 --- A-Accuracy:0.7087980173482032  E-Number:0.0  T-Cost:1.834862385321101  Time:253ms<br />
Run 34 --- A-Accuracy:0.7236679058240396  E-Number:0.0  T-Cost:1.834862385321101  Time:272ms<br />
Run 35 --- A-Accuracy:0.7211895910780669  E-Number:0.0  T-Cost:1.834862385321101  Time:261ms<br />
Run 36 --- A-Accuracy:0.6964064436183395  E-Number:0.0  T-Cost:1.834862385321101  Time:284ms<br />
Run 37 --- A-Accuracy:0.7286245353159851  E-Number:0.0  T-Cost:1.834862385321101  Time:255ms<br />
Run 38 --- A-Accuracy:0.7856257744733581  E-Number:0.0  T-Cost:1.834862385321101  Time:274ms<br />
Run 39 --- A-Accuracy:0.7211895910780669  E-Number:0.0  T-Cost:1.834862385321101  Time:253ms<br />
Run 40 --- A-Accuracy:0.7100371747211895  E-Number:0.0  T-Cost:1.834862385321101  Time:254ms<br />
Run 41 --- A-Accuracy:0.7063197026022305  E-Number:0.0  T-Cost:1.834862385321101  Time:280ms<br />
Run 42 --- A-Accuracy:0.7757125154894672  E-Number:0.0  T-Cost:1.834862385321101  Time:265ms<br />
Run 43 --- A-Accuracy:0.7397769516728625  E-Number:0.0  T-Cost:1.834862385321101  Time:272ms<br />
Run 44 --- A-Accuracy:0.8004956629491945  E-Number:0.0  T-Cost:1.834862385321101  Time:283ms<br />
Run 45 --- A-Accuracy:0.7434944237918215  E-Number:0.0  T-Cost:1.834862385321101  Time:271ms<br />
Run 46 --- A-Accuracy:0.7434944237918215  E-Number:0.0  T-Cost:1.834862385321101  Time:267ms<br />
Run 47 --- A-Accuracy:0.7372986369268897  E-Number:0.0  T-Cost:1.834862385321101  Time:288ms<br />
Run 48 --- A-Accuracy:0.7546468401486989  E-Number:0.0  T-Cost:1.834862385321101  Time:263ms<br />
Run 49 --- A-Accuracy:0.6344485749690211  E-Number:0.0  T-Cost:1.834862385321101  Time:331ms<br />
Run 50 --- A-Accuracy:0.724907063197026  E-Number:0.0  T-Cost:1.834862385321101  Time:262ms<br />

Average --- A-Accuracy:0.7270631970260222  E-Number:0.0  T-cost:1.8365137614678908  Time:281ms
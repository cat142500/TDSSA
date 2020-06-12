/** ExtendedTD.java provides the component of extended truth discovery for TDSSA,
 *  which iteratively infers the true label of tasks and the quality of workers
 *  with the Sybil score and reliability score of workers taken into consideration.
 */

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

public class ExtendedTD {
	private int L; // label size
	
	/* initialization */
	public ExtendedTD( int L) {
		this.L = L;
	}
	
	/* iteratively run extended truth discovery*/
	public void process(Set<Task> tasks, Set<Worker> workers) {
		// set the initial weight of workers to their accuracy on golden tasks
		for(Worker worker : workers) {
			worker.setWeight(worker.getP());
		}
		
		int iteration = 0;
		while(iteration<1000) {
			iteration++;
			int difference = 0;
			
			// extended label aggregation
			for(Task task : tasks) {
				int original_label = task.getAggregated();
				ArrayList<Worker> assigned = task.getAssigned();
				double[] votes = new double[L];
				for(int i=0; i<L; i++) {
					votes[i] = 0;
				}
				for(Worker worker : assigned) {
					int label = worker.getLabeledPairs().get(task);
					votes[label] += worker.getS()/L + (1-worker.getS())*worker.getWeight();
				}
				double max_vote = -1;
				for(int i=0; i<L; i++) {
					if(votes[i]>max_vote) {
						task.setAggregated(i);
						max_vote = votes[i];
					}
				}
				if(original_label!=task.getAggregated()) {
					difference++;
				}
			}
			
			// terminate if converge
			if(difference==0) {
				break;
			}
			
			// extended weight estimation
			for(Worker worker : workers) {
				double correct = 0;
				double count = 0;
				Map<Task, Integer> assigned = worker.getLabeledPairs();
				for(Task task : assigned.keySet()) {
					if(assigned.get(task).intValue()==task.getAggregated()) {
						correct += task.getCi();
					}
					count += task.getCi();
				}
				if(count>0) {
					worker.setWeight(correct/count);
				}
			}
		}
	}
}

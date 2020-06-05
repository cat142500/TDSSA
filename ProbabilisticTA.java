import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class ProbabilisticTA {
	private double tau; // Sybil threshold for banning workers
	private double delta; // reliability threshold for marking reliable workers
	private double alpha; // probability to assign a golden task to a new worker
	private int K; // number of workers per task
	
	
	/* initialization */
	public ProbabilisticTA(double tau, double delta, double alpha, int K) {
		this.tau = tau;
		this.delta = delta;
		this.alpha = alpha;
		this.K = K;
	}
	
	
	/* probabilically assign a task to a requesting worker based on the worker's Sybil score and reliability score */
	public Task assign(Worker worker, Set<Task> golden_tasks, Set<Task> normal_tasks) {
		if(worker.getS()<tau && worker.getR()<delta) {
			double g = alpha*(1-worker.getR())+(1-alpha)*worker.getS();
			Random rand = new Random();
			// assign a golden task with g probability
			if(rand.nextDouble()<=g) {
				Set<Task> avail_tasks = new HashSet<Task>(golden_tasks);
				avail_tasks.removeAll(worker.getLabeledPairs().keySet());
				for(Task task : avail_tasks) {
					if(worker.getAttackerID()==-1 && !worker.getPairs().keySet().contains(task)) {
						continue;
					}
					else {
						int count = task.getExpose();
						for(Worker curr_worker : task.getAssigned()) {
							if(curr_worker.getR()<delta) {
								count++;
							}
						}
						if(count>=K) {
							continue;
						}
						else {
							task.assign(worker);
							return task;
						}
					}
				}
			}
		}

		// assign a normal task
		for(Task task : normal_tasks) {
			ArrayList<Worker> assigned_workers = task.getAssigned();					
			if(!assigned_workers.contains(worker)) {
				if(worker.getAttackerID()!=-1) {
					task.assign(worker);
					return task;
				}
				else {
					if(task.getWorkers().contains(worker)) {
						task.assign(worker);
						return task;
					}
				}
			}
		}
		return null;
	}
}

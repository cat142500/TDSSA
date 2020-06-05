import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class Attacker {
	private ArrayList<Worker> sybil_workers; // Sybil workers controlled by the attacker
	private Map<Task, Integer> task_labels; // randomized label on each task
	private Map<Task, Integer> task_count; // observation times of each task
	private int attacker_id; // attacker ID
	private int K; // number of workers for each task
	private int L; // label size;
	
	/* initialization */
	public Attacker(int attacker_id, int K, int L) {
		sybil_workers = new ArrayList<Worker>();
		task_labels = new HashMap<Task, Integer>();
		task_count = new HashMap<Task, Integer>();
		this.attacker_id = attacker_id;
		this.K = K;
		this.L = L;
	}
	
	/* add a Sybil worker */
	public void addWorker(Worker worker) {
		sybil_workers.add(worker);
	}
	
	/* return Sybil workers controlled by the attacker */
	public ArrayList<Worker> getWorkers() {
		return sybil_workers;
	}
	
	/* set a randomized label on a task */
	public void setTaskLabel(Task task, Integer label) {
		task_labels.put(task, label);
	}
	
	/* return the label on a task */
	public int getTaskLabel(Task task) {
		return task_labels.get(task);
	}
	
	/* update the observation times of a task */
	public void observe(Task task) {
		if(task_count.containsKey(task)) {
			int count = task_count.get(task);
			count++;
			task_count.put(task, count);
			// label the task honestly if the task is observed for more than K times
			if(count==K+1) {
				Random rand = new Random();
				// assume the attacker has 0.8 probability to provide the true label once a golden task is identified
				if(rand.nextDouble()<=0.8) {
					setTaskLabel(task, task.getTrueLabel());
				}
				else {
					int label = rand.nextInt(L);
					while(label==task.getTrueLabel()) {
						label = rand.nextInt(L);
					}
					setTaskLabel(task, label);
				}
			}
		}
		else {
			task_count.put(task, 1);
		}
	}
	
	/* return the observation times of a task */
	public int getCount(Task task) {
		if(task_count.containsKey(task)) {
			return task_count.get(task);
		}
		return 0;
	}

	/* return the attacker ID */
	public int getAttackerId() {
		return attacker_id;
	}
}

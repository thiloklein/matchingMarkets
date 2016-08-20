import static choco.Choco.*;
import choco.cp.model.CPModel;
import choco.cp.solver.CPSolver;
import choco.kernel.model.Model;
import choco.kernel.solver.Solver;
import choco.kernel.solver.ContradictionException;
import choco.kernel.solver.constraints.integer.AbstractLargeIntSConstraint;
import choco.kernel.solver.variables.integer.IntDomainVar;
import choco.kernel.memory.IStateInt;
import choco.kernel.memory.trailing.StoredInt;
import choco.kernel.model.variables.integer.IntegerExpressionVariable;
import choco.kernel.model.variables.integer.IntegerVariable;

//
// SRN: N-ary constraint for Stable Roommates 
//

public class SRN extends AbstractLargeIntSConstraint {

    private int n;
    private int[][] rank;         // rank[i][j] = k <-> agent_i ranks agent_j as k^th choice
    private int[][] pref;         // pref[i][k] = j <-> agent_i has agent_j as k^th choice
    private int[] length;         // length of agent's preference list
    private IStateInt[] upb;      // upper bound on domain
    private IStateInt[] lwb;      // lower bound on domain
    private IntDomainVar[] agent; // the agents
    private boolean trace;

    public SRN(Solver s,IntDomainVar[] agent,int[][] pref, int[][] rank, int[] length) {
        super(agent);
        n           = agent.length;
        this.agent  = agent;
        this.pref   = pref;
        this.rank   = rank;
        this.length = length;
        upb         = new StoredInt[n];
        lwb         = new StoredInt[n];
        for (int i=0;i<n;i++){
            upb[i]  = s.getEnvironment().makeInt(length[i]);
            lwb[i]  = s.getEnvironment().makeInt(0);
        }
        trace = false;
    }

    void display(int i){
        System.out.print(i +": ");
        for (int j=0;j<n;j++) if (agent[i].getDomain().contains(j)) System.out.print(pref[i][j] +" ");
        System.out.println();
    }

    void display(){
        for (int i=0;i<n;i++) display(i);
        System.out.println();
    }

    public boolean isSatisfied(){
        if (trace){System.out.println("isSatisfied()"); display();}
        return true;
    }

    public void awake() throws ContradictionException {
        if (trace){System.out.println("awake()"); display();}
        for (int i=0;i<n;i++) awakeOnInf(i);
    }

    public void propagate() throws ContradictionException {
        if (trace){System.out.println("propagate()"); display();}
    }

    public void awakeOnInf(int i) throws ContradictionException {
        if (trace){System.out.println("awakeOnInf(" + i + ")"); display();}
        int x = agent[i].getInf(); // best (lowest) rank for agent i
        int j = pref[i][x];
        agent[j].setSup(rank[j][i]);
        for (int w=lwb[i].get();w<x;w++){
            int h = pref[i][w];
            agent[h].setSup(rank[h][i]-1);
        }
        lwb[i].set(x);
    }
    //
    // deltaMin(i)
    // lowerbound increases on agent[i]
    // agents are rejected from head of preference list and
    // all agents that agent[i] prefers to current agent
    // must find partners better than agent[i]
    //

    public void awakeOnSup(int i) throws ContradictionException {
        if (trace){System.out.println("awakeOnSup(" + i + ")"); display();}
        int x = agent[i].getSup(); // worst (largest) preference for agent[i]
        for (int y=x+1;y<=upb[i].get();y++){
            int j = pref[i][y];
            agent[j].remVal(rank[j][i]);
        }
        upb[i].set(x);
    }
    //
    // deltaMax(i)
    // upper bound of agent[i] decreases
    // agents are reject from tail of prefernce list
    //

    public void awakeOnRem(int i,int x) throws ContradictionException {
        if (trace){System.out.println("awakeOnRem(" + i + "," + x + ")"); display();}
        int j = pref[i][x];
        agent[j].remVal(rank[j][i]);
    }

    public void awakeOnInst(int i) throws ContradictionException {
        if (trace){System.out.println("awakeOnInst(" + i + "," + vars[i].getVal() + ")"); display();}
        int y = agent[i].getVal();
        for (int x = lwb[i].get();x<y;x++){
            int j = pref[i][x];
            agent[j].setSup(rank[j][i]-1);
        }
        for (int z=y+1;z<=upb[i].get();z++){
            int j = pref[i][z];
            agent[j].remVal(rank[j][i]);
        }
        lwb[i].set(y);
        upb[i].set(y);
    }
    //
    // NOTES
    //
    // x.getInf() gets lower bound
    // x.getSup() gets upper bound
    // x.getVal(i)
    // x.remVal(i)
    // x.setInf(i) updates lower bound
    // x.setSup(i) updates upper bound
}
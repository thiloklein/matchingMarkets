//
// Toolkit constraint encoding
//
import java.io.*;
import java.util.*;
import static choco.Choco.*;
import choco.cp.model.CPModel;
import choco.cp.solver.CPSolver;
import choco.kernel.model.Model;
import choco.kernel.solver.Solver;
import choco.kernel.model.variables.integer.IntegerVariable;
import choco.kernel.solver.ContradictionException;
import choco.cp.solver.search.integer.varselector.StaticVarOrder;

public class smi {

    public String sayHello(String[] args) throws IOException, ContradictionException {
    //public static void main(String[] args) throws IOException, ContradictionException {

        int n;
        int[][] rank; // rank[i][j] = k <-> agent_i ranks agent_j as k^th choice
        int[][] pref; // pref[i][k] = j <-> agent_i has agent_j as k^th choice
        int[] length; // length of agent's preference list
        Model model;
        Solver solver;
        IntegerVariable[] agent; // domain of ranks, last is unmatched
        long totalTime, modelTime, solveTime, readTime, modelSize;
        boolean search;
        int solutions, matchingSize;

        search      = true;
        totalTime   = System.currentTimeMillis();
        readTime    = System.currentTimeMillis();

        // --- read ---

        //BufferedReader fin = new BufferedReader(new FileReader(args[0]));
        String[] fin = args[0].split("n");
        //n                  = Integer.parseInt(fin.readLine());
        n                  = Integer.parseInt(fin[0]);
        pref               = new int[n][n];
        rank               = new int[n][n];
        length             = new int[n];
        for (int i=0;i<n;i++){
            //StringTokenizer st = new StringTokenizer(fin.readLine()," ");
            StringTokenizer st = new StringTokenizer(fin[i+1]," ");
            int k = 0;
            length[i] = 0;
            while (st.hasMoreTokens()){
                int j      = Integer.parseInt(st.nextToken()) - 1;
                rank[i][j] = k;
                pref[i][k] = j;
                length[i]  = length[i] + 1;
                k          = k + 1;
            }
            rank[i][i] = k;
            pref[i][k] = i;
        }
        //fin.close();

        // ---

        readTime = System.currentTimeMillis() - readTime;

        // --- build ---

        modelTime = System.currentTimeMillis();
        model     = new CPModel();
        agent     = new IntegerVariable[n];
        for (int i=0;i<n;i++) agent[i] = makeIntVar("agent_"+ i,0,length[i],"cp:enum");
        model.addVariables(agent);
        solver    = new CPSolver();
        solver.read(model);
        solver.post(new SRN(solver,solver.getVar(agent),pref,rank,length));
        modelTime  = System.currentTimeMillis() - modelTime;
        modelSize  = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024; // kilobytes

        // --- solve ---

        solutions = matchingSize = 0;
        solveTime = System.currentTimeMillis();
        solver.setVarIntSelector(new StaticVarOrder(solver,solver.getVar(agent)));

        StringBuilder sb = new StringBuilder(); // create empty StringBuilder instance

        if (solver.solve().booleanValue()) {

            // --- getMatchingSize ---
            matchingSize = 0;
            for (int i = 0; i < n; i++)
                if (solver.getVar(agent[i]).getVal() < length[i]) matchingSize++;
            matchingSize = matchingSize / 2;

            // --- displayMatching ---

            for (int i = 0; i < n; i++) {
                int j = pref[i][solver.getVar(agent[i]).getVal()];
                if (i < j) sb.append(1 + "," + (i + 1) + "," + (j + 1) + "\n");
            }
            //sb.append("\n");

            solutions = 1;
            while (solver.nextSolution().booleanValue()) {
                solutions++;

                // --- displayMatching ---

                for (int i = 0; i < n; i++) {
                    int j = pref[i][solver.getVar(agent[i]).getVal()];
                    if (i < j) sb.append(solutions + "," + (i + 1) + "," + (j + 1) + "\n");
                }
                //sb.append("\n");
            }
        }
        solveTime = System.currentTimeMillis() - solveTime;
        totalTime = System.currentTimeMillis() - totalTime;

        // --- stats ---

        sb.append("\n");
        sb.append("solutions:"+ solutions +"\n");
        if (search) sb.append("nodes:"+ solver.getNodeCount() +"\n");
        sb.append("modelTime:"+ modelTime +"\n");
        if (search) sb.append("solveTime:"+ solveTime +"\n");
        sb.append("totalTime:"+ totalTime +"\n");
        sb.append("modelSize:"+ modelSize +"\n");
        sb.append("readTime:"+ readTime +"\n");
        sb.append("matchingSize:"+ matchingSize);

        String result = sb.toString();
        //System.out.println(result);
        return(result);
    }

    public static void main(String[] args) throws IOException, ContradictionException {
    }
}
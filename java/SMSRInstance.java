//
// Allow a mixing of SM and SR using the morphing technique
// "graph morph (type B)" in AAAI-99
//
import java.util.*;

public class SMSRInstance {

    Random gen;
    int n; // number of agents, even
    int m; // number of edges
    ArrayList[] pref; // preference list of agents
    ArrayList<Integer> cliqueE; // all possible edges
    ArrayList<Integer> bicliqueE; // all possible edges in bipatrite graph
    ArrayList<Integer> mixE; // edges in mixed SM & SR, of size m
    boolean trace;

    void shuffle(ArrayList<Integer> v){
        int j,n;
        n = v.size();
        int temp;
        for (int i=n-1;i>0;i--){
            j = gen.nextInt(i+1);
            temp = v.get(i);
            v.set(i,v.get(j));
            v.set(j,temp);
        }
    }
    //
    // Knuth shuffle of integer vector
    //

    SMSRInstance(int n){
        gen       = new Random();
        this.n    = n;
        this.pref = new ArrayList[n];
        m         = (n*n)/4; // number of edges
        cliqueE   = new ArrayList<Integer>(); // all possible edges
        bicliqueE = new ArrayList<Integer>(); // all possible edges
        mixE      = new ArrayList<Integer>(); // edges in mixed SM & SR, of size m
        for (int i=0;i<n-1;i++)
            for (int j=i+1;j<n;j++)
                cliqueE.add(i*n+j);
        for (int i=0;i<n/2;i++)
            for (int j=n/2;j<n;j++)
                bicliqueE.add(i*n+j);
        for (int i=0;i<n;i++) pref[i] = new ArrayList<Integer>();
    }

    void build(String modelType,double p){
        if (modelType.equals("A")) buildTypeA(p);
        if (modelType.equals("B")) buildTypeB(p);
    }


    void buildTypeB(double p){
        int mSR                      = (int)(p*m); // number of edges to select from SR
        int mSM                      = m - mSR;  // number of edges to select from SM
        int mIntersection            = 0; // thes size of the intersection between smE and srE
        ArrayList<Integer> srE       = new ArrayList<Integer>(); // edges in SR (with same number of edges in SM)
        ArrayList<Integer> smE       = new ArrayList<Integer>(); // edges in SM
        ArrayList<Integer> deleteE   = new ArrayList<Integer>(); // edges to delete
        mixE.clear();
        //
        // create edges in complete bipartite graph smE
        //
        for (int i=0;i<n/2;i++)
            for (int j=n/2;j<n;j++)
                smE.add(i*n+j);
        //
        // create a random SRI from m shuffled edges in clique
        //
        shuffle(cliqueE);
        for (int i=0;i<m;i++) srE.add(cliqueE.get(i));
        //
        // force intersection of smE and srE in mixE
        //
        displayEdges("sm: ",smE);
        displayEdges("sr: ",srE);
        for (int e : srE)
            if (smE.contains(e)){
                mixE.add(e);
                mIntersection++;
            }
        for (Integer e : mixE){
            smE.remove((Integer)e);
            srE.remove((Integer)e);
        }
        //
        // select mSM and mSR edges and add to mix
        //
        mSR = (int)((m - mIntersection)*p);
        mSM = m - mIntersection - mSR;
        shuffle(smE);
        for (int i=0;i<mSM;i++) mixE.add(smE.get(i));
        shuffle(srE);
        for (int i=0;i<mSR;i++) mixE.add(srE.get(i));
        displayEdges("mix:",mixE);
        //
        //Create preference lists
        //
        for (int i=0;i<n;i++) pref[i].clear();
        for (Integer e : mixE){
            int v = e/n;
            int w = e%n;
            pref[v].add(w);
            pref[w].add(v);
        }
        for (int i=0;i<n;i++) shuffle(pref[i]);
    }

    void buildTypeA(double p){
        int mSR                      = (int)(p*m); // number of edges to select from SR
        int mSM                      = m - mSR;  // number of edges to select from SM
        ArrayList<Integer> additions = new ArrayList<Integer>();
        mixE.clear();
        //
        // create a random SM instance with mSM edges
        //
        shuffle(bicliqueE);
        for (int i=0;i<mSM;i++) mixE.add(bicliqueE.get(i));
        //
        // add in mSR random edges not arleady in mixE
        //
        displayEdges("sm: ",mixE);
        shuffle(cliqueE);
        for (int i=0;i<n*(n-1)/2 && mSR>0;i++){
            int edge = cliqueE.get(i);
            if (!mixE.contains(edge)){
                mixE.add(edge);
                additions.add(edge);
                mSR--;
            }
        }
        //
        //Create preference lists
        //
        displayEdges("sr: ",additions);
        displayEdges("mix:",mixE);
        for (int i=0;i<n;i++) pref[i].clear();
        for (Integer e : mixE){
            int v = e/n;
            int w = e%n;
            pref[v].add(w);
            pref[w].add(v);
        }
        for (int i=0;i<n;i++) shuffle(pref[i]);
    }

    void display(){
        System.out.println(n);
        for (int i=0;i<n;i++){
            for (int j : (ArrayList<Integer>)pref[i]) System.out.print((j+1) +" ");
            System.out.println();
        }
    }

    void displayEdges(String name,ArrayList<Integer> E){
        if (!trace) return;
        System.out.print(name +" ");
        System.out.print("|"+ E.size() +"| ");
        for (int e : E) System.out.print("("+ (e/n+1) +","+ (e%n+1) +") ");
        System.out.println();
    }

    public static void main(String args[]){
        int n    = Integer.parseInt(args[0]); // must be even
        double p = Double.parseDouble(args[1]); // (1.0-p) from SM (p=0 <-> SM), p from SR (p=1 <-> SR)
        SMSRInstance inst = new SMSRInstance(n);
        inst.trace = args.length > 3;
        if (args[2].equals("A")) inst.buildTypeA(p);
        if (args[2].equals("B")) inst.buildTypeB(p);
        inst.display();
    }
}
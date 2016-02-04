package stacey;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeUtils;
import beast.util.Randomizer;
import stacey.debugtune.Checks;
import stacey.util.Bindings;
import stacey.util.BitUnion;
import stacey.util.UnionArrays;

import java.util.ArrayList;
import java.util.List;

/**
 * Created  by Work on 14/12/2015.
 */
public class LineageRecombiner extends Operator {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> smcTreeInput =
            new Input<Tree>("smcTree",
                    "The species tree or minimal clusters tree", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<Tree>> geneTreesInput =
            new Input<List<Tree>>("geneTree",
                    "All gene trees",
                    new ArrayList<Tree>());

    public Input<RealParameter> popSFInput =
            new Input<>("popSF",
                    "The population scaling factor for the STACEY coalescent", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<PIOMSCoalescentDistribution> piomscdInput =
            new Input<PIOMSCoalescentDistribution>("piomsCoalDist",
                    "The PIOMSCoalescentDistribution",
                    Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Long> delayInput =
            new Input<>("delay",
                    "Number of times the operator is disabled.");


    private Tree sTree;
    private List<Tree> gTrees;
    private long delay = 0;
    private int callCount = 0;
    private boolean sTreeTooSmall;
    private Bindings bindings;
    private UnionArrays unionArrays;
    private ArrayList<ArrayList<ArrayList<Coalescence>>> coalescences;

    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;



    private static class Coalescence {
        int nodeNr;
        double height;

        public Coalescence(int nodeNr, double height) {
            this.nodeNr = nodeNr;
            this.height = height;
        }

        int getNodeNr() {
            return nodeNr;
        }

        double getHeight() {
            return height;
        }
    }


    private static class SegmentOfLineages {
        private final ArrayList<Integer> nodeNrs;
        private final double startHeight;
        private final double endHeight;

        public SegmentOfLineages(ArrayList<Integer> nodeNrs, double startHeight, double endHeight) {
            this.nodeNrs = nodeNrs;
            this.startHeight = startHeight;
            this.endHeight = endHeight;
        }

        public int getNLineages() {
            return nodeNrs.size();
        }

        public int getNodeNrOfLineage(int i) {
            return nodeNrs.get(i);
        }

        public double getStartHeight() {
            return startHeight;
        }

        public double getEndHeight() {
            return endHeight;
        }


    }


    @Override
    public void initAndValidate() {
        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        sTreeTooSmall = (sTree.getLeafNodeCount() < 3);
        if (delayInput != null  &&  delayInput.get() != null) {
            delay = delayInput.get().longValue();
        }
        bindings = Bindings.initialise(sTree, gTrees);
        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
        initCoalescenceArray(sTree, gTrees);
    }



    @Override
    public double proposal() {
        // TODO Input says use get(this). See TODO in CoordinatedPruneRegraft
        if (sTreeTooSmall) {
            return Double.NEGATIVE_INFINITY;
        }
        callCount++;
        if (callCount < delay) {
            return Double.NEGATIVE_INFINITY;
        }

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            sTree = smcTreeInput.get();
            gTrees = geneTreesInput.get();
            Checks.allTreesAndCompatibility(sTree, gTrees, "LineageRecombiner", "before move");
            numberofdebugchecks++;
        }

        int j = Randomizer.nextInt(gTrees.size());
        unionArrays.updateSMCTreeAndGTree(j);
        double logHR = recombineLineage(j);   //  The business
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            sTree = smcTreeInput.get();
            gTrees = geneTreesInput.get();
            Checks.allTreesAndCompatibility(sTree, gTrees, "LineageRecombiner", "after move");
            numberofdebugchecks++;
        }
        return logHR;
    }



    private void initCoalescenceArray(TreeInterface sTree, List<Tree> gTrees) {
        int numSNodes = sTree.getNodeCount();
        int numGTrees = gTrees.size();
        coalescences = new ArrayList<>(numSNodes);
        for (int n = 0; n < numSNodes; n++) {
            coalescences.add(new ArrayList<>(numGTrees));
            for (int j = 0; j < numGTrees; j++) {
                coalescences.get(n).add(new ArrayList<>(0));
            }
        }
    }


    // Returns false if gtree j doesn't fit, else updates a 3D array of heights.
    private boolean updateCoalescencesForOneGTree(int j) {
        // clear old coalescences
        for (ArrayList<ArrayList<Coalescence>> coals: coalescences) {
            coals.get(j).clear();
        }
        Node node = makeCoalescencesSubtree(j, gTrees.get(j).getRoot());
        return (node != null);
    }



    private double recombineLineage(int j) {
        double logHR = 0.0;
        Tree gTree = gTrees.get(j);
        // Choose a node which will have its parent re-attached
        int nodeNr = Randomizer.nextInt(gTree.getNodeCount());
        while (gTree.getNode(nodeNr).isRoot()) {
            nodeNr = Randomizer.nextInt(gTree.getNodeCount());
        }
        Node gNode = gTree.getNode(nodeNr);

        // find which stree branch nodeNr is in.
        int sNodeNr = unionArrays.hostNodeNrOfGNode(j, gNode);



        // TODO
        // Find all time segments going back from rp to gtree root within stree branches
        // between coalescences/stree nodes. Each segment has a time interval and a number
        // of lineages. The number of lineages does not include the (old) lineage from rp back
        // to its first coalescence x. The number of lineages decreases by one at each coalescence,
        // and increases by at least one at each stree node.
        ArrayList<SegmentOfLineages> segments = new ArrayList<>();
        updateCoalescencesForOneGTree(j);


        // Get an 'average' population size from the current popSF and the inv-gamma-mixture middlingValue().
        // In case it matters, mv and cf are fixed for the analysis.
        double mv = piomscdInput.get().getInverseGammaMixture().middlingValue();
        double cf = piomscdInput.get().coalFactor(j);
        double popSF = popSFInput.get().getValue();
        double theta = popSF * mv * cf;

        // TODO
        // Choose a new attachment (which might or might not change the topology).
        // The time is chosen from the density implied by the inhomogeneous Poisson
        // process with intensity given by the number of lineages and the population size.
        // The lineage is a randomly chosen one from those available at that time.

        // Going backwards in time from gNode, there are n(t) other lineages. (t measures time from
        // gNode, not zero.) We go backwards through stree branches to the root.
        // n(t) decreases by one at every coalescence and increases by at least one each time we
        // pass a side branch in the stree. n(t) can be zero initially, but since gNode is not the root
        // of the gene tree, n(t) must become one eventually, in the stree root, and then stays one forever.
        // The intensity of the (inhomogeneous) Poisson process is p(t) = n(t)/theta.
        // The density of the time of the new coalescence is
        //           f(t)  =  p(t) exp[ -integral_0^t p(u) du ].
        // A random branch is chosen from the n(t) available at time t.

        // TODO
        // calculate a HR based on the new and old density values (and lineage choice?)
        // The choice of gene tree node gNode is symmetric.
        // The sampling of a new attachment point is based on a density and needs to be accounted for.

        // The density along a branch at time t is f(t)/n(t)
        //           f(t)/n(t)  = (1/theta) exp[ -integral_0^t p(u) du ].
        // Since we want to take the ratio of two things like this for the HR, the (1/theta) can be ignored.
        // So the HR is (old density/new density) =
        //          exp[ -integral_0^s p(u) du ] / exp[ -integral_0^t p(u) du ].
        // where s and t are the times of the detachment and reattachment respectively.
        // This is
        //         exp[ (1/theta) integral_s^t n(u) du ]
        // and logHR is
        //             (1/theta) integral_s^t n(u) du
        // If the reattachment is older than the detachment, t>s and logHR>0.
        // The integral is easy as it is a step function.

        // TODO carry out the move

        return logHR;
    }




    //  TODO Pasted code of makeHeightsSubtree() from FitsHeights, which I thought was similar.
    // started editing...
    // I want to start at a single gtree node and time, not leaves ??
    private Node makeCoalescencesSubtree(int j, final Node hereG) {
        Node lftG = hereG.getChild(0);
        Node lftS;
        if (lftG.isLeaf()) {
            lftS = sTree.getNode(bindings.smcTipNrFromGTreeTipNr(j, lftG.getNr()));
        } else {
            lftS = makeCoalescencesSubtree(j, lftG);
            assert lftS != null;
        }
        Node rgtG = hereG.getChild(1);
        Node rgtS;
        if (rgtG.isLeaf()) {
            rgtS = sTree.getNode(bindings.smcTipNrFromGTreeTipNr(j, rgtG.getNr()));
        } else {
            rgtS = makeCoalescencesSubtree(j, rgtG);
            assert rgtS != null;
        }
        while (lftS != rgtS) {
            if (lftS.getHeight() < rgtS.getHeight()) {
                lftS = lftS.getParent();
            } else {
                rgtS = rgtS.getParent();
            }
        }
        Node hereS = lftS;
        double heightG = hereG.getHeight();
        if (hereS.getHeight() > heightG) {
            return null;
        }
        Node hostS = hereS;
        while (!hostS.isRoot()  &&  hostS.getParent().getHeight() < heightG) {
            hostS = hostS.getParent();
        }
        Coalescence coalG = new Coalescence(hereG.getNr(), heightG);
        coalescences.get(hostS.getNr()).get(j).add(coalG);
        return hereS;
    }


}

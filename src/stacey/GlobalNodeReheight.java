package stacey;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import stacey.debugtune.Checks;
import stacey.util.Bindings;
import stacey.util.BitUnion;
import stacey.util.Misc;
import stacey.util.UnionArrays;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Created by Graham Jones on 26/07/2015.
 */
public class GlobalNodeReheight extends Operator {

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<Tree> smcTreeInput =
            new Input<Tree>("smcTree",
                    "The species tree or minimal clusters tree", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<List<Tree>> geneTreesInput =
            new Input<List<Tree>>("geneTree",
                    "All gene trees",
                    new ArrayList<Tree>());

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

    private final boolean debugFlag = Boolean.valueOf(System.getProperty("stacey.debug"));
    private int numberofdebugchecks = 0;
    private final static int maxnumberofdebugchecks = 100000;




    @Override
    public void initAndValidate() throws Exception {
        /*if (smcTreeInput.get().getLeafNodeCount() < 3) {
            throw new Exception("CoordinatedPruneRegraft cannot be used if there are less than 3 minimal clusters.");
        } TODO: this upsets Beauti */
        sTree = smcTreeInput.get();
        gTrees = geneTreesInput.get();
        sTreeTooSmall = (sTree.getLeafNodeCount() < 4);
        // I want at least three internal heights, so two don't change (but could relax that)
        if (delayInput.get() != null) {
            delay = delayInput.get().longValue();
        }
        bindings = Bindings.initialise(sTree, gTrees);
        unionArrays = UnionArrays.initialise(sTree, gTrees, bindings);
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
            System.out.println("Doing debug check. numberofdebugchecks = " + numberofdebugchecks);
            Checks.allTreesAndCompatibility(sTree, gTrees, "CoordinatedNodeReheight", "before move");
            numberofdebugchecks++;
        }

        unionArrays.update();
        double logHR = doCoordinatedNodeReheight();   //  The business
        unionArrays.reset();

        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            sTree = smcTreeInput.get();
            gTrees = geneTreesInput.get();
            System.out.println("Doing debug check. numberofdebugchecks = " + numberofdebugchecks);
            Checks.allTreesAndCompatibility(sTree, gTrees, "CoordinatedNodeReheight", "after move");
            numberofdebugchecks++;
        }

        return logHR;
    }


    enum PMType {SNODES_PASS, GNODES_PASS, SNODE_GNODE_EQUAL, SNODE_ZERO, GNODE_ZERO};

    class PivotalMoment {

        PMType type;
        double displacement;
        int j;
        int mNodeO;
        int nNodeO;
        // If type is SNODES_PASS,       mNodeO, nNodeO are the sTree nodes        j is ignored
        // If type is GNODES_PASS,       j is gTree, mNodeO, nNodeO are the gTree nodes
        // If type is SNODE_GNODE_EQUAL, j is gTree, mNodeO is the sTree node, nNodeO is the gTree node
        // If type is SNODE_ZERO,        mNodeO is the sTree node                  j and nNodeO are ignored
        // If type is GNODE_ZERO,        j is gTree, mNodeO is the gTree node      nNodeO is ignored

        PivotalMoment(PMType type, int j, int mNodeO, int nNodeO) {
            this.type = type;
            this.displacement = Double.NaN;
            this.j = j;
            this.mNodeO = mNodeO;
            this.nNodeO = nNodeO;
        }

        void setDisplacement(double displacement) {
            this.displacement = displacement;
        }

    }


    /*

        BEAST node Nrs                 left-right order
    2     0     1    4     3        0     2     4    6     8
    |     |       \  |    /         |     |       \  |    /
     \   /          6    /           \   /          5    /
      \ /             \ /             \ /             \ /
       5               7               1               7
        \             /                \             /
          \         /                    \         /
            \     /                       \     /
             \  /                          \  /
              8                             3

    The following are constants during the move, and are used during it:

    The maximum height of all gTrees, maxGTreeHeight. The changes to node heights
    are a dilation/contraction with this as the fixed point.

    The heights of all sTree nodes except the one chosen to move.
    These other heights are used to set a scale.
    (The identity of the chosen node needs care, see next.)

    Once the left-right orderings have been chosen, the orderings of the tips
    are invariant. In the example above, order and invOrder are like this:
        order invOrder
    0     2      2
    1            4
    2     0      0
    3            8
    4     1      6
    5
    6     4
    7
    8     3
    order[i] maps left-right index i to tip nodeNr.
    invOrder[nodeNr] maps back to left-right index.
    The rest of order[] and invOrder[] can be calculated from these.
    After the move, the BEAST node Nrs for internal nodes may have changed.
    But the left-right position sNodeO of the chosen node does not change, so that
    is used when calculating Hastings ratio, which requires considering the reverse move.

    Note that more than one gTree node in a single gtree may change height.

    the allowed height interval
    for the reverse move, for logHR.

    */
    private double doCoordinatedNodeReheight() {
        double logHR = 0.0;
        // room for orderings, inverse orderings, weights, internal node counts,
        // and unions of nodes in trees.
        // Note only odd-numbered weight indices used.
        int [] sInvOrder = new int[sTree.getNodeCount()];
        int[][] gInvOrders = new int[gTrees.size()][];
        int [] sOrder = new int[sTree.getNodeCount()];
        int[][] gOrders = new int[gTrees.size()][];
        double [] sWeights = new double[sTree.getNodeCount()];
        double [][] gWeights = new double[gTrees.size()][];
        int sINCount = sTree.getInternalNodeCount();
        int [] gINCounts = new int[gTrees.size()];
        BitUnion [] sOUnions = new BitUnion[sTree.getNodeCount()];
        BitUnion [] sOLftUnions = new BitUnion[sTree.getNodeCount()];
        BitUnion [] sORgtUnions = new BitUnion[sTree.getNodeCount()];
        BitUnion [][] gOUnions = new BitUnion[gTrees.size()][];
        for (int j = 0; j < gInvOrders.length; j++) {
            Tree gTree = gTrees.get(j);
            gInvOrders[j] = new int[gTree.getNodeCount()];
            gOrders[j] = new int[gTree.getNodeCount()];
            gWeights[j] = new double[gTree.getNodeCount()];
            gINCounts[j] = gTree.getInternalNodeCount();
            gOUnions[j] = new BitUnion[gTree.getNodeCount()];
        }
        // Randomly orient the internal nodes to make a ordering of tips.
        makeRandomTipOrderingOfSubtree(sTree, sTree.getRoot(), 0, sOrder, sInvOrder);
        for (int j = 0; j < gInvOrders.length; j++) {
            Tree gTree = gTrees.get(j);
            makeRandomTipOrderingOfSubtree(gTree, gTree.getRoot(), 0, gOrders[j], gInvOrders[j]);
        }
        // Fill in the rest of the orderings from the tips
        fillInOrderings(sInvOrder, sOrder, gInvOrders, gOrders);

        fillInStreeUnions(sOrder, sOUnions, sOLftUnions, sORgtUnions);
        fillInGTreeUnions(gOrders, gOUnions);

        chooseWeights(sWeights, gWeights);

        boolean moveIsUp = Randomizer.nextBoolean();

        double eta;
        // now find the allowed interval which keeps move compatible
        if (moveIsUp) {
            double upperBound = findUpperBound(sOrder, gOrders, sWeights, gWeights, sINCount, gINCounts);
            if (!Double.isFinite(upperBound)  ||  !(upperBound > 0.0)) {
                System.out.println("BUG in doCoordinatedNodeReheight(), bad upperBound");
            }
            eta = Randomizer.uniform(0.0, upperBound);
        } else {
            double lowerBound = findLowerBound(sOrder, gOrders, sWeights, gWeights, sINCount, gINCounts);
            if (!Double.isFinite(lowerBound)  ||  !(lowerBound < 0.0)) {
                System.out.println("BUG in doCoordinatedNodeReheight(), bad lowerBound");
            }
            eta = Randomizer.uniform(lowerBound, 0.0);
        }

        String allTreesAsTextBefore = "";
        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            allTreesAsTextBefore = Misc.allTreesAsText(sTree, gTrees);
        }

        // edit the trees: first change heights
        for (int mO = 1; mO < sOrder.length; mO += 2) {
            double weight = sWeights[mO];
            double dh = weight * eta;
            Node node = sTree.getNode(sOrder[mO]);
            node.setHeight(node.getHeight() +  dh);
        }
        for (int j = 0; j < gTrees.size(); j++) {
            for (int mO = 1; mO < gOrders[j].length; mO += 2) {
                double weight = gWeights[j][mO];
                double dh = weight * eta;
                Node node = gTrees.get(j).getNode(gOrders[j][mO]);
                node.setHeight(node.getHeight() +  dh);
            }
        }
        // reconstruct the trees, changing their topologies as needed for the new heights
        reconstructFromHeights(sTree, sOrder);
        for (int j = 0; j < gTrees.size(); j++) {
            Tree gTree = gTrees.get(j);
            reconstructFromHeights(gTree, gOrders[j]);
        }


        if (debugFlag  &&  numberofdebugchecks < maxnumberofdebugchecks) {
            sTree = smcTreeInput.get();
            gTrees = geneTreesInput.get();
            System.out.println("Doing debug check. numberofdebugchecks = " + numberofdebugchecks);
            if (!Checks.allTreesAreOKAndCompatible(sTree, gTrees)) {
                System.out.println("Before:");
                System.out.println(allTreesAsTextBefore);

                System.out.print("sOrder[] =");
                for (int i = 0; i < sOrder.length; i++) {
                    System.out.print(" [" + i + "]=" + sOrder[i]);
                }
                System.out.println();
                for (int j = 0; j < gTrees.size(); j++) {
                    System.out.print("gOrders[" + j + "] =");
                    for (int i = 0; i < gOrders[j].length; i++) {
                        System.out.print(" [" + i + "]=" + gOrders[j][i]);
                    }
                    System.out.println();
                }
                System.out.println("eta = " + eta);

                System.out.println("After:");
                System.out.println(Misc.allTreesAsText(sTree, gTrees));
            }
        }
        return logHR;
    }






    private void chooseWeights(double [] sWeights, double [][] gWeights) {
        for (int k = 1; k < sWeights.length; k += 2) {
            sWeights[k] = Randomizer.nextDouble();
        }
        int k1 = 2 * Randomizer.nextInt(sTree.getInternalNodeCount()) + 1;
        sWeights[k1] = 1.0;
        for (int j = 0; j < gWeights.length; j++) {
            for (int k = 1; k < gWeights[j].length; k += 2) {
                gWeights[j][k] = Randomizer.nextDouble();
            }
        }
    }


    private double findLowerBound(int [] sOrder, int [][] gOrders,
                                  double [] sWeights, double [][] gWeights,
                                  int sINCount, int [] gINCounts) {
        double lower = Double.NEGATIVE_INFINITY;

        PriorityQueue<PivotalMoment> priorityQueue = makePriorityQueue(
                sOrder, gOrders, sWeights, gWeights, sINCount, gINCounts, false);

        PivotalMoment pm = priorityQueue.poll();
        while (pm != null) {
            while (pm != null) {
                switch (pm.type) {
                    case SNODES_PASS: // TODO
                        break;

                    case GNODES_PASS://  TODO update gtree BitUnions
                        break;

                    case SNODE_GNODE_EQUAL://  TODO if BitUnions olap, pm.displacement determines bound
                        break;

                    case SNODE_ZERO:// TODO
                        break;

                    case GNODE_ZERO:// TODO
                        break;

                    default:
                        System.exit(1);
                }
            }

        }
        return lower;

    }

    private double findUpperBound(int [] sOrder, int [][] gOrders,
                                 double [] sWeights, double [][] gWeights,
                                 int sINCount, int [] gINCounts) {
        double upper = Double.POSITIVE_INFINITY;

        PriorityQueue<PivotalMoment> priorityQueue = makePriorityQueue(
                sOrder, gOrders, sWeights, gWeights, sINCount, gINCounts, true);

        PivotalMoment pm = priorityQueue.poll();
        while (pm != null) {
            switch (pm.type) {
                case SNODES_PASS:// TODO
                    break;

                case GNODES_PASS:// TODO update gtree BitUnions
                    break;

                case SNODE_GNODE_EQUAL:// TODO if BitUnions olap, pm.displacement determines bound
                    break;

                default:
                    System.exit(1);
            }
        }
        return upper;
    }





    private boolean addPivotalMomentIfNeeded(PriorityQueue<PivotalMoment> priorityQueue,
                                             boolean moveIsUp,
                                             double mH, double nH, double mW, double nW,
                                             PivotalMoment pm) {
        double displacement = Double.NaN;
        if (moveIsUp) {
            if (nH > mH) {
                if (mW > nW) { displacement = (nH - mH) / (mW - nW); }
            } else {
                if (nW > mW) { displacement = (mH - nH) / (nW - mW); }
            }
        } else {
            if (nH > mH) {
                if (nW > mW) { displacement = (nH - mH) / (nW - mW); }
            } else {
                if (mW > nW) { displacement = (mH - nH) / (mW - nW); }
            }
        }
        if (!Double.isNaN(displacement)) {
            assert Double.isFinite(displacement);
            pm.setDisplacement(displacement);
            priorityQueue.add(pm);
            return true;
        } else {
            return false;
        }
    }



    private PriorityQueue<PivotalMoment> makePriorityQueue(int [] sOrder, int [][] gOrders,
                                                   double [] sWeights, double [][] gWeights,
                                                   int sINCount, int [] gINCounts,
                                                   boolean moveIsUp) {
        int nPMs = 0;
        nPMs += sINCount;
        nPMs += sINCount * (sINCount - 1) / 2;
        for (int j = 0; j < gOrders.length; j++) {
            nPMs += gINCounts[j];
            nPMs += gINCounts[j] * (gINCounts[j] - 1) / 2;
            nPMs += sINCount * gINCounts[j];
        }

        PriorityQueue<PivotalMoment> priorityQueue = new PriorityQueue<>(nPMs, PIVOTALMOMENT_TIME_ORDER);

        for (int mO = 1;  mO < sWeights.length;  mO += 2) {
            double mHeight = sTree.getNode(sOrder[mO]).getHeight();
            double mWeight = sWeights[mO];
            for (int nO = mO + 2;  nO < sWeights.length;  nO += 2) {
                double nHeight = sTree.getNode(sOrder[nO]).getHeight();
                double nWeight = sWeights[nO];
                PivotalMoment pm = new PivotalMoment(PMType.SNODES_PASS, -1, mO, nO);
                addPivotalMomentIfNeeded(priorityQueue, moveIsUp, mHeight, nHeight, mWeight, nWeight, pm);
            }
        }
        for (int j = 0; j < gOrders.length; j++) {
            for (int mO = 1;  mO < gWeights[j].length;  mO += 2) {
                double mHeight = gTrees.get(j).getNode(gOrders[j][mO]).getHeight();
                double mWeight = gWeights[j][mO];
                for (int nO = mO + 2;  nO < gWeights[j].length;  nO += 2) {
                    double nHeight = gTrees.get(j).getNode(gOrders[j][nO]).getHeight();
                    double nWeight = gWeights[j][nO];
                    PivotalMoment pm = new PivotalMoment(PMType.GNODES_PASS, j, mO, nO);
                    addPivotalMomentIfNeeded(priorityQueue, moveIsUp, mHeight, nHeight, mWeight, nWeight, pm);
                }
            }
        }
        for (int j = 0; j < gOrders.length; j++) {
            for (int mO = 1;  mO < sWeights.length;  mO += 2) {
                double mHeight = sTree.getNode(sOrder[mO]).getHeight();
                double mWeight = sWeights[mO];
                for (int nO = 1;  nO < gWeights[j].length;  nO += 2) {
                    double nHeight = gTrees.get(j).getNode(gOrders[j][nO]).getHeight();
                    double nWeight = gWeights[j][nO];
                    double displacement;
                    assert nHeight > mHeight;
                    displacement = (nHeight - mHeight) / (mWeight-nWeight);
                    PivotalMoment pm = new PivotalMoment(PMType.SNODE_GNODE_EQUAL, j, mO, nO);
                    addPivotalMomentIfNeeded(priorityQueue, moveIsUp, mHeight, nHeight, mWeight, nWeight, pm);
                }
            }
        }
        if (!moveIsUp) {
            for (int mO = 1;  mO < sWeights.length;  mO += 2) {
                double mHeight = sTree.getNode(sOrder[mO]).getHeight();
                double mWeight = sWeights[mO];
                double displacement = mHeight / mWeight;
                PivotalMoment pm = new PivotalMoment(PMType.SNODE_ZERO, -1, mO, -1);
                addPivotalMomentIfNeeded(priorityQueue, moveIsUp, mHeight, 0.0, mWeight, 0.0, pm);
            }
            for (int j = 0; j < gOrders.length; j++) {
                for (int mO = 1;  mO < gWeights[j].length;  mO += 2) {
                    double mHeight = gTrees.get(j).getNode(gOrders[j][mO]).getHeight();
                    double mWeight = gWeights[j][mO];
                    double displacement = mHeight / mWeight;
                    PivotalMoment pm = new PivotalMoment(PMType.GNODE_ZERO, j, mO, -1);
                    addPivotalMomentIfNeeded(priorityQueue, moveIsUp, mHeight, 0.0, mWeight, 0.0, pm);
                }
            }
        }

        assert priorityQueue.size() == nPMs;
        return priorityQueue;
    }




    private int makeRandomTipOrderingOfSubtree(Tree tree, Node node, int nextO, int[] order, int[] invOrder) {
        if (node.isLeaf()) {
            invOrder[node.getNr()] = nextO;
            order[nextO] = node.getNr();
            nextO += 2;
            return nextO;
        } else {
            boolean swap = Randomizer.nextBoolean();
            nextO = makeRandomTipOrderingOfSubtree(tree, node.getChild(swap ? 1 : 0), nextO, order, invOrder);
            nextO = makeRandomTipOrderingOfSubtree(tree, node.getChild(swap ? 0 : 1), nextO, order, invOrder);
            return nextO;
        }
    }


    private void fillInOrderings(int [] sInvOrder, int [] sOrder, int [][] gInvOrders, int [][] gOrders) {
        fillInInternalNodeOrder(sTree, sInvOrder, sOrder);
        for (int j = 0; j < gInvOrders.length; j++) {
            Tree gTree = gTrees.get(j);
            fillInInternalNodeOrder(gTree, gInvOrders[j], gOrders[j]);
        }
    }

    private void fillInInternalNodeOrder(Tree tree, int[] invOrder, int[] order) {
        int [][] ranges = new int[order.length][4];
        fillInInternalNodeOrderSubtree(tree.getRoot(), invOrder, order, ranges);
    }


    private void fillInInternalNodeOrderSubtree(Node node, int[] invOrder, int[] order, int[][] ranges) {
        int n = node.getNr();
        if (node.isLeaf()) {
            assert invOrder[n] % 2 == 0;
            ranges[n][0] = ranges[n][1] = ranges[n][2] = ranges[n][3] = invOrder[n];
        } else {
            fillInInternalNodeOrderSubtree(node.getChild(0), invOrder, order, ranges);
            fillInInternalNodeOrderSubtree(node.getChild(1), invOrder, order, ranges);
            int c0 = node.getChild(0).getNr();
            int c1 = node.getChild(1).getNr();
            if (ranges[c0][0] < ranges[c1][0] ) {
                ranges[n][0] = ranges[c0][0];
                ranges[n][1] = ranges[c0][3];
                ranges[n][2] = ranges[c1][0];
                ranges[n][3] = ranges[c1][3];
            } else {
                ranges[n][0] = ranges[c1][0];
                ranges[n][1] = ranges[c1][3];
                ranges[n][2] = ranges[c0][0];
                ranges[n][3] = ranges[c0][3];
            }
            assert ranges[n][1] + 2 == ranges[n][2];
            int nO = ranges[n][1] + 1;
            invOrder[node.getNr()] = nO;
            order[nO] = node.getNr();

        }
    }


    private void fillInStreeUnions(int [] sOrder, BitUnion [] sOUnions, BitUnion [] sOLftUnions, BitUnion [] sORgtUnions) {
        // tips
        for (int mO = 0; mO < sOrder.length; mO +=2) {
            Node node = sTree.getNode(sOrder[mO]);
            assert node.isLeaf();
            sOUnions[mO] = new BitUnion(sTree.getLeafNodeCount());
            sOUnions[mO].replaceWith(unionArrays.sNodeUnion(node.getNr()));
        }
        for (int mO = 1; mO < sOrder.length; mO +=2) {
            Node node = sTree.getNode(sOrder[mO]);
            sOLftUnions[mO] = new BitUnion(sTree.getLeafNodeCount());
            sORgtUnions[mO] = new BitUnion(sTree.getLeafNodeCount());
            double mOHeight = sTree.getNode(sOrder[mO]).getHeight();
            // find the first node to left of mO which is higher than mO
            // if no such node, lftIntO becomes -1
            int lftIntO = mO - 2;
            while (lftIntO > 0  &&  sTree.getNode(sOrder[lftIntO]).getHeight() < mOHeight) {
                lftIntO -= 2;
            }
            // find the first node to right of mO which is higher than mO
            // if no such node, rgtIntO becomes sOrder.length
            int rgtIntO = mO + 2;
            while (rgtIntO < sOrder.length  &&  sTree.getNode(sOrder[rgtIntO]).getHeight() < mOHeight) {
                rgtIntO += 2;
            }
            // collect species (or min clusters) from left-hand range
            // If height of sNodeO decreases, these nodes could belong to
            // the children of iO
            for (int lftTipO = mO - 1;  lftTipO >= lftIntO + 1;  lftTipO -= 2) {
                Node nodeOLft = sTree.getNode(sOrder[lftTipO]);
                assert nodeOLft.isLeaf();
                sOLftUnions[mO].insert(nodeOLft.getNr());
            }
            // collect species (or min clusters) from right-hand range
            for (int rgtTipO = mO + 1; rgtTipO <= rgtIntO - 1; rgtTipO += 2) {
                final Node nodeORgt = sTree.getNode(sOrder[rgtTipO]);
                assert node.isLeaf();
                sORgtUnions[mO].insert(node.getNr());
            }
        }
    }




    private void fillInGTreeUnions(int [][] gOrders, BitUnion [][] gOUnions) {
        for (int j = 0; j < gOrders.length; j++) {
            Tree gTree = gTrees.get(j);
            // tips
            for (int mO = 0; mO < gOrders[j].length; mO +=2) {
                Node node = sTree.getNode(gOrders[j][mO]);
                assert node.isLeaf();
                gOUnions[j][mO] = new BitUnion(sTree.getLeafNodeCount());
                gOUnions[j][mO].replaceWith(unionArrays.sNodeUnion(node.getNr()));
            }

            for (int mO = 1; mO < gOrders[j].length; mO +=2) {
                Node node = gTree.getNode(gOrders[j][mO]);
                gOUnions[j][mO] = new BitUnion(sTree.getLeafNodeCount());
                double mOHeight = sTree.getNode(gOrders[j][mO]).getHeight();
                // find the first node to left of mO which is higher than mO
                // if no such node, lftIntO becomes -1
                int lftIntO = mO - 2;
                while (lftIntO > 0  &&  sTree.getNode(gOrders[j][lftIntO]).getHeight() < mOHeight) {
                    lftIntO -= 2;
                }
                // find the first node to right of mO which is higher than mO
                // if no such node, rgtIntO becomes gOrders[j].length
                int rgtIntO = mO + 2;
                while (rgtIntO < gOrders[j].length  &&  sTree.getNode(gOrders[j][rgtIntO]).getHeight() < mOHeight) {
                    rgtIntO += 2;
                }
                // collect species (or min clusters) from left-hand range
                // If height of sNodeO decreases, these nodes could belong to
                // the children of iO
                for (int lftTipO = mO - 1;  lftTipO >= lftIntO + 1;  lftTipO -= 2) {
                    Node nodeOLft = gTree.getNode(gOrders[j][lftTipO]);
                    assert nodeOLft.isLeaf();
                    gOUnions[j][mO].insert(nodeOLft.getNr());
                }
                // collect species (or min clusters) from right-hand range
                for (int rgtTipO = mO + 1; rgtTipO <= rgtIntO - 1; rgtTipO += 2) {
                    final Node nodeORgt = sTree.getNode(gOrders[j][rgtTipO]);
                    assert node.isLeaf();
                    gOUnions[j][mO].insert(node.getNr());
                }

            }
        }
    }




    private void reconstructFromHeights(Tree tree, int [] order) {
        Node root = reconstructSubtreeFromHeights(tree, order, 0, order.length);
        root.setParent(null);
        tree.setRoot(root);
    }



    private Node reconstructSubtreeFromHeights(Tree tree, int [] order, int from, int to) {
        assert from % 2 == 0;
        assert to % 2 == 1;
        if (from + 1 == to) {
            return tree.getNode(order[from]);
        }
        int rootO = highestNode(tree, order, from, to);
        Node root = tree.getNode(order[rootO]);
        Node lft  = reconstructSubtreeFromHeights(tree, order, from, rootO);
        Node rgt  = reconstructSubtreeFromHeights(tree, order, rootO + 1, to);
        root.setChild(0, lft);
        root.getChild(0).setParent(root);
        root.setChild(1, rgt);
        root.getChild(1).setParent(root);
        return root;
    }



    private int highestNode(Tree tree, int [] order, int from, int to) {
        int rootO = -1;
        double maxHeight = Double.NEGATIVE_INFINITY;
        for (int i = from + 1;  i < to;  i+=2) {
            double h = tree.getNode(order[i]).getHeight();
            if (h > maxHeight) {
                maxHeight = h;
                rootO = i;
            }
        }
        if (rootO < 0) {
            System.out.println("BUG in highestNode()");
        }
        assert rootO >= 0;
        return rootO;
    }



    static final Comparator<Node> NODE_HEIGHT_ORDER =
            new Comparator<Node>() {
                public int compare(Node x, Node y) {
                    if (x.getHeight() == y.getHeight()) {
                        return 0;
                    } else {
                        return x.getHeight() > y.getHeight() ? 1 : -1;
                    }
                }
            };


    static final Comparator<PivotalMoment> PIVOTALMOMENT_TIME_ORDER =
            new Comparator<PivotalMoment>() {
                public int compare(PivotalMoment x, PivotalMoment y) {
                    if (x.displacement == y.displacement) {
                        return 0;
                    } else {
                        return x.displacement > y.displacement ? 1 : -1;
                    }
                }
            };

}
